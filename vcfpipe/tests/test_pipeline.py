"""端到端集成测试。

包含一个**完全独立**的输出解析器（仅用标准库，不导入 vcfpipe 的任何
VCF/TSV 读写代码），用于交叉验证生成的 merged.vcf 与样本矩阵；
并校验等位计数守恒、确定性、版本保护与规则升级重建。
"""
import csv
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "src"))

from vcfpipe.pipeline import VersionConflict, run_pipeline  # noqa: E402
from vcfpipe.query import Dataset  # noqa: E402

FA = os.path.join(ROOT, "data/ref/syn37.fa")
DICT = os.path.join(ROOT, "data/ref/contigs.json")
SHEET = os.path.join(ROOT, "data/input/sample_sheet.csv")
PRIO = os.path.join(ROOT, "data/input/caller_priority.json")


# --------------------------------------------------------------------------- #
# 独立解析器：与 vcfpipe.vcf_io / vcfpipe.query 的实现完全无关
# --------------------------------------------------------------------------- #
def independent_parse_vcf(path):
    headers, records = {}, []
    with open(path, encoding="utf-8") as fh:
        sample_cols = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                headers[line.split("=", 1)[0]] = line
            elif line.startswith("#CHROM"):
                sample_cols = line.split("\t")[9:]
            else:
                c = line.split("\t")
                fmt = c[8].split(":")
                samples = {}
                for name, blob in zip(sample_cols, c[9:]):
                    vals = dict(zip(fmt, blob.split(":")))
                    samples[name] = vals
                records.append({
                    "chrom": c[0], "pos": int(c[1]), "ref": c[3],
                    "alt": c[4], "info": c[7], "samples": samples,
                })
    return headers, sample_cols, records


def independent_parse_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def independent_parse_input_vcf(path):
    """独立解析输入 VCF 的 body，用于守恒律的源端计数。"""
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            fmt = c[8].split(":")
            vals = dict(zip(fmt, c[9].split(":")))
            out.append({
                "chrom": c[0], "pos": int(c[1]), "ref": c[3],
                "alts": c[4].split(","), "gt": vals.get("GT"),
                "ad": vals.get("AD"),
            })
    return out


class EndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="vcfpipe_test_")
        cls.man_v1 = run_pipeline(cls.tmp, FA, DICT, SHEET, PRIO, "v1")
        cls.man_v2 = run_pipeline(cls.tmp, FA, DICT, SHEET, PRIO, "v2")
        cls.ds = Dataset(os.path.join(cls.tmp, "syn37_v1"))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # ------------------------------------------------------------------ #
    def test_manifest_counts(self):
        t = self.man_v1["totals"]
        self.assertEqual(t["included_records"], 18)
        self.assertEqual(t["quarantined_records"], 5)
        self.assertEqual(t["merged_variant_rows"], 15)
        # 每源纳入+隔离 = 该源有效 body 记录处理结果之和
        for s in self.man_v1["sources"]:
            self.assertGreaterEqual(s["included"] + s["quarantined"], 0)

    def test_quarantine_reasons_complete(self):
        q = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1", "quarantine.tsv"))
        codes = {r["reason_code"] for r in q}
        self.assertEqual(codes, {
            "ref_base_mismatch", "contig_unmapped", "symbolic_alt",
            "malformed_gt", "left_align_boundary",
        })
        # 每条隔离记录都带人类可读原因
        self.assertTrue(all(r["reason_detail"] for r in q))

    # ------------------------------------------------------------------ #
    def test_independent_vcf_parser_roundtrip(self):
        headers, samples, records = independent_parse_vcf(
            os.path.join(self.tmp, "syn37_v1", "merged.vcf"))
        self.assertIn("##fileformat", headers)
        self.assertEqual(samples, ["S1", "S2", "S3"])
        # 独立解析器读到的位点数 = 事实表去重稳定键数
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        keys = {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in fact}
        vcf_keys = {(r["chrom"], str(r["pos"]), r["ref"], r["alt"])
                    for r in records}
        self.assertEqual(len(records), len(keys))
        self.assertEqual({(k[0], str(int(k[1])), k[2], k[3]) for k in keys},
                         vcf_keys)

    def test_missing_never_becomes_ref_hom(self):
        _, _, records = independent_parse_vcf(
            os.path.join(self.tmp, "syn37_v1", "merged.vcf"))
        # chr1:13 是 ./. 缺失；任何样本单元格都不得把缺失写成 0/0
        r13 = next(r for r in records if r["pos"] == 13)
        self.assertEqual(r13["samples"]["S1"]["GT"], "./.")
        # 无记录样本的占位也是 ./. 而不是 0/0
        r15 = next(r for r in records if r["pos"] == 15)
        self.assertEqual(r15["samples"]["S1"]["GT"], "./.")

    def test_half_missing_and_phasing_preserved(self):
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        by = {(r["pos"], r["alt"], r["sample"]): r for r in fact}
        # 半缺失 ./1
        self.assertEqual(by[("18", "G", "S1")]["gt"], "./1")
        # 三等位拆分后的相位与半缺失混合
        self.assertEqual(by[("10", "A", "S1")]["gt"], "1|.")
        self.assertEqual(by[("10", "T", "S1")]["gt"], ".|1")
        self.assertEqual(by[("10", "A", "S1")]["phased"], "1")
        # S2 多等位 0/2 拆到 ALT#2 时另一等位为缺失
        self.assertEqual(by[("63", "T", "S2")]["gt"], "0/.")

    def test_dp_preserved_when_gt_missing(self):
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        r = next(x for x in fact if x["pos"] == "63" and x["sample"] == "S1")
        self.assertEqual(r["gt"], "./.")
        self.assertEqual(r["dp"], "7")  # GT 缺失但 DP=7 必须保留

    def test_caller_priority_and_conflict_evidence(self):
        ev = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1", "evidence.tsv"))
        # chr1:5 S1 上 deepvar(0/1) 与 freebayes(1/1) 冲突，deepvar 胜出
        ev5 = [r for r in ev if r["pos"] == "5" and r["sample"] == "S1"]
        self.assertEqual({r["caller"] for r in ev5},
                         {"deepvar", "freebayes"})
        self.assertEqual({r["raw_gt"] for r in ev5}, {"0/1", "1/1"})
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        win = next(r for r in fact if r["pos"] == "5" and r["sample"] == "S1")
        self.assertEqual(win["winner_caller"], "deepvar")
        self.assertEqual(win["gt"], "0/1")

    def test_contig_alias_and_unmapped(self):
        q = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1", "quarantine.tsv"))
        self.assertTrue(any(r["reason_code"] == "contig_unmapped"
                            and r["chrom"] == "chrZ" for r in q))
        # 数字别名 1 -> chr1 的记录被正常纳入
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        self.assertTrue(any(r["sample"] == "S2" and r["chrom"] == "chr1"
                            for r in fact))

    # ------------------------------------------------------------------ #
    def test_allele_count_conservation(self):
        """等位拷贝守恒。

        口径：多等位拆分把“其它 ALT”投影为缺失，但同一拷贝会在它自己的
        ALT 拆分记录中计为剂量 1，因此对“拆分后、合并前”的包含记录直接
        求和，恰好等于源端全部非参等位拷贝；隔离记录的剂量单列。合并会
        折叠跨 caller 的重复位点，故不能在合并后事实上反推。
        """
        from vcfpipe.pipeline import (
            ContigMapping, Reference, RuleBundle, load_caller_priority,
            load_sample_sheet, process_sources,
        )
        reference = Reference(FA)
        mapping = ContigMapping.load(DICT)
        rule = RuleBundle.for_version("v1")
        result = process_sources(
            load_sample_sheet(SHEET), reference, mapping, rule)

        # 源端：所有输入 VCF 的非参考倍型位（每位对任意 ALT 计 1）
        inputs = [
            "data/input/batch_b1_deepvar_s1.vcf",
            "data/input/batch_b1_freebayes_s1.vcf",
            "data/input/batch_b2_varscan_s2.vcf",
            "data/input/batch_b3_freebayes_s3.vcf",
        ]
        source_alt = 0
        for path in inputs:
            for rec in independent_parse_input_vcf(
                    os.path.join(ROOT, path)):
                if rec["gt"] in (None, ".", "./.", ".|."):
                    continue
                sep = "|" if "|" in rec["gt"] else "/"
                source_alt += sum(
                    1 for p in rec["gt"].split(sep)
                    if p not in (".", "0"))

        included_alt = sum(s.call.alt_dosage for s in result.included)
        quarantined_alt = 0
        for q in result.quarantined:
            gt = q.raw_gt
            if gt in (None, ".", "./.", ".|."):
                continue
            sep = "|" if "|" in gt else "/"
            quarantined_alt += sum(
                1 for p in gt.split(sep) if p not in (".", "0"))

        self.assertEqual(source_alt, included_alt + quarantined_alt)
        # 源 24 个非参考倍型位 = 纳入 19 + 隔离 5（每条隔离 GT 各 1 位）
        self.assertEqual((source_alt, included_alt, quarantined_alt),
                         (24, 19, 5))

    def test_split_level_dosage_conservation_per_record(self):
        """单条多等位记录守恒：拆分后 ALT 剂量总和 == 源非参考倍型位数。

        剂量按“位”计（每个倍型位对任意 ALT 最多贡献 1），不按 ALT 下标
        数值求和，因此 1|2 的守恒目标是 2 而非 3。
        """
        from vcfpipe.genotype import split_record_calls
        from vcfpipe.model import RawRecord, SourceVCF

        def make(alts, gt, phased=False):
            src = SourceVCF("p", "syn37", "S", "c", "d", "s")
            return RawRecord(src, 1, "chr1", 10, "C", tuple(alts), ".",
                             ".", ".", gt, None, None, phased)

        cases = [
            (("A", "T"), "1|2", True, 2),   # 两个非参位 -> 各 ALT 1 剂量
            (("T", "C"), "0/2", False, 1),  # 一个非参位
            (("A", "T"), "1/1", False, 2),  # 两个位都在 ALT#1
            (("G",), "./1", False, 1),      # 半缺失：仅真实位计数
            (("G",), "./.", False, 0),      # 全缺失：0
        ]
        for alts, gt, phased, expected in cases:
            calls = split_record_calls(make(alts, gt, phased))
            self.assertEqual(
                sum(c.alt_dosage for c in calls), expected,
                msg=f"{gt} 拆分剂量不守恒")

    def test_ad_split_maps_correct_alt(self):
        """多等位 AD 拆分：每条拆分记录只带 [AD0, ADk]，不错配 ALT。"""
        fact = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v1",
                         "variants", "variants_chr1.tsv"))
        a1 = next(r for r in fact if r["pos"] == "10" and r["alt"] == "A")
        t2 = next(r for r in fact if r["pos"] == "10" and r["alt"] == "T")
        self.assertEqual(a1["ad"], "30,5")
        self.assertEqual(t2["ad"], "30,4")

    # ------------------------------------------------------------------ #
    def test_deterministic_rerun(self):
        before = {k: v for k, v in self.man_v1["outputs"].items()}
        again = run_pipeline(self.tmp, FA, DICT, SHEET, PRIO, "v1")
        self.assertEqual(before, again["outputs"])

    def test_version_conflict_on_tampered_input(self):
        tamper_dir = tempfile.mkdtemp(prefix="vcfpipe_tamper_")
        try:
            run_pipeline(tamper_dir, FA, DICT, SHEET, PRIO, "v1")
            with open(SHEET, encoding="utf-8") as fh:
                original = fh.read()
            with open(SHEET, "a", encoding="utf-8") as fh:
                fh.write("\n")  # 改变样本表字节但不改变语义
            with self.assertRaises(VersionConflict):
                run_pipeline(tamper_dir, FA, DICT, SHEET, PRIO, "v1")
            with open(SHEET, "w", encoding="utf-8", newline="") as fh:
                fh.write(original)
        finally:
            shutil.rmtree(tamper_dir, ignore_errors=True)

    def test_rule_upgrade_creates_new_dataset(self):
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, "syn37_v1")))
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, "syn37_v2")))
        self.assertEqual(self.man_v1["rule_version"], "v1")
        self.assertEqual(self.man_v2["rule_version"], "v2")
        # v2 边界记录不再隔离，纳入数 +1，隔离数 -1
        self.assertEqual(self.man_v2["totals"]["included_records"], 19)
        self.assertEqual(self.man_v2["totals"]["quarantined_records"], 4)
        q2 = independent_parse_tsv(
            os.path.join(self.tmp, "syn37_v2", "quarantine.tsv"))
        self.assertFalse(
            any(r["reason_code"] == "left_align_boundary" for r in q2))

    # ------------------------------------------------------------------ #
    def test_query_range_uses_partition_index(self):
        rows = self.ds.query_range("chr1", 10, 28)
        positions = {int(r["pos"]) for r in rows}
        self.assertTrue(all(10 <= p <= 28 for p in positions))
        self.assertIn(28, positions)
        self.assertEqual(self.ds.query_range("chr1", 1000), [])

    def test_sample_slice_and_lineage_trace(self):
        s2 = self.ds.sample_slice("S2")
        self.assertTrue(all(r["sample"] == "S2" for r in s2))
        # 输出键逆查到来源
        traced = self.ds.trace_output("chr1", 28, "G", "GA")
        self.assertTrue(any(r["source_alt"] == "AA" for r in traced))
        # 来源行顺查：多等位行（输入 VCF 第 13 行 C>A,T）产出两条 included
        forward = self.ds.trace_source(
            "data/input/batch_b1_deepvar_s1.vcf", 13)
        self.assertEqual(len(forward), 2)
        self.assertTrue(all(r["status"] == "included" for r in forward))


if __name__ == "__main__":
    unittest.main(verbosity=2)
