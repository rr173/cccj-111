"""家系三人组筛查的手算回归、独立穷举核对与失败例子进程回归。

覆盖任务书要求：
* 家系表全部结构性非法情形（缺列/三人相同/重复行/矛盾父母/环/缺样本/
  重复样本列），错误点名行号或人物链，不丢家庭后继续；
* 位点与调用校验（列数、POS、REF/ALT、FORMAT、GT 语法/越界、负 DP、
  非二倍体），坏记录带文件行号使命令失败且不发布输出；
* 补全枚举的三类判定（含缺失、多等位、相位无关、低深度整调用未知）；
* 独立穷举核对器（tools/trio_exhaustive_check.py，不导入产品补全/
  传递函数）对全部明细逐行复核；
* 正常数据后接坏行：后半段坏记录不留下看似完成的摘要；
* 摘要与逐行明细可对上，白名单外只计跳过数。
"""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "src")
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "trio_exhaustive_check.py")
sys.path.insert(0, SRC)

from vcfpipe.trio import (  # noqa: E402
    COMPATIBLE, INCOMPATIBLE, UNCERTAIN, UNSUPPORTED,
    TrioConfigError, TrioRecordError,
    diploid_candidates, evaluate_trio, load_families,
    read_vcf_header, run_trio, transmission_witness,
)
from vcfpipe.trio import MemberCall  # noqa: E402

DEMO_VCF = os.path.join(ROOT, "data", "trio", "trio_demo.vcf")
DEMO_FAM = os.path.join(ROOT, "data", "trio", "families_demo.tsv")
CHROMS = ("chr1", "chr2")
MIN_DP = 10


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
def make_call(role, person, gt, dp=30, min_dp=MIN_DP):
    """从 GT 文本直接构造 MemberCall（走产品解析，供单元判定使用）。"""
    if gt == ".":
        alleles, phased = None, False
    else:
        sep = "|" if "|" in gt else "/"
        parts = gt.split(sep)
        alleles = tuple(None if p == "." else int(p) for p in parts)
        if all(a is None for a in alleles):
            alleles = None
        phased = sep == "|"
    trusted = dp is not None and dp >= min_dp
    reason = "ok" if trusted else ("dp_missing" if dp is None
                                   else "dp_below_min")
    return MemberCall(role=role, person=person, raw_gt=gt, dp=dp,
                      trusted=trusted, trust_reason=reason,
                      alleles=alleles, phased=phased)


def evaluate(father_gt, mother_gt, child_gt, n_alts,
             dps=(30, 30, 30), min_dp=MIN_DP):
    calls = {
        "father": make_call("father", "F", father_gt, dps[0], min_dp),
        "mother": make_call("mother", "M", mother_gt, dps[1], min_dp),
        "child": make_call("child", "C", child_gt, dps[2], min_dp),
    }
    return evaluate_trio(calls, n_alts)


HEADER = ("##fileformat=VCFv4.2\n"
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
          "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
          "F\tM\tC\n")


def write_vcf(path, body, header=HEADER, gz=False):
    text = header + body
    if gz:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
            fh.write(text)
    else:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)


def write_families(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("family\tchild\tfather\tmother\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def locus_line(chrom, pos, ref, alt, f, m, c, dp=30, fmt="GT:DP"):
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}\t"
            f"{f}:{dp}\t{m}:{dp}\t{c}:{dp}\n")


# --------------------------------------------------------------------------- #
# 手算单元测试：补全与三类判定
# --------------------------------------------------------------------------- #
class CompletionHandCalc(unittest.TestCase):
    def test_unordered_dedup_phase_irrelevant(self):
        # 0/. 的补全 = {0/0, 0/1}；1/0 与 0/1 同一补全
        self.assertEqual(diploid_candidates((0, None), 1),
                         ((0, 0), (0, 1)))
        self.assertEqual(diploid_candidates((1, 0), 1), ((0, 1),))
        self.assertEqual(diploid_candidates((None, None), 1),
                         ((0, 0), (0, 1), (1, 1)))

    def test_multiallelic_completion_uses_full_range(self):
        # ./2 在 2 个 ALT 的位点：缺失位可取 0/1/2
        self.assertEqual(diploid_candidates((None, 2), 2),
                         ((0, 2), (1, 2), (2, 2)))
        # 编号 2 必须在本位点合法（不会被当成另一位点的 1）
        self.assertEqual(diploid_candidates((2, 2), 2), ((2, 2),))

    def test_fully_unknown_dot(self):
        self.assertEqual(diploid_candidates((None, None), 2),
                         ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)))

    def test_transmission_witness(self):
        self.assertEqual(transmission_witness((1, 1), (2, 2), (1, 2)),
                         (1, 2))
        self.assertIsNone(transmission_witness((1, 1), (0, 0), (1, 1)))


class VerdictHandCalc(unittest.TestCase):
    def test_normal_transmission_compatible(self):
        r = evaluate("0/0", "0/1", "0/1", 1)
        self.assertEqual(r.verdict, COMPATIBLE)
        self.assertEqual((r.total, r.explainable), (1, 1))
        self.assertIn("paternal=0", r.witness)

    def test_two_distinct_alts_compatible(self):
        # 父 1/1、母 2/2、孩子 1/2：多等位整体判断
        r = evaluate("1/1", "2/2", "1/2", 2)
        self.assertEqual(r.verdict, COMPATIBLE)
        self.assertEqual(r.witness.split(";")[1:],
                         ["paternal=1", "maternal=2"])

    def test_ref_ref_child_alt_incompatible(self):
        r = evaluate("0/0", "0/0", "0/1", 1)
        self.assertEqual(r.verdict, INCOMPATIBLE)
        self.assertEqual(r.explainable, 0)
        self.assertIsNone(r.witness)  # 不可解释不编造传递

    def test_partial_missing_still_incompatible(self):
        # 任务书例子：父 0/.、母 0/0、孩子 1/1 -> 确定不一致
        r = evaluate("0/.", "0/0", "1/1", 1)
        self.assertEqual(r.verdict, INCOMPATIBLE)
        self.assertEqual((r.total, r.explainable), (2, 0))

    def test_partial_missing_all_completions_compatible(self):
        # 任务书例子：二位点父母皆 0/1、孩子 0/. -> 全部补全可解释
        r = evaluate("0/1", "0/1", "0/.", 1)
        self.assertEqual(r.verdict, COMPATIBLE)
        self.assertEqual((r.total, r.explainable), (2, 2))

    def test_feasible_and_infeasible_mix_uncertain(self):
        # 父 0/.、母 0/0、孩子 0/1：补 0 可解释、补 1 不可解释
        r = evaluate("0/.", "0/0", "0/1", 1)
        self.assertEqual(r.verdict, UNCERTAIN)
        self.assertEqual((r.total, r.explainable), (2, 1))

    def test_multiallelic_missing_uncertain(self):
        # 父 ./2、母 0/0、孩子 1/1：三种父补全无一能给孩子 1/1
        r = evaluate("./2", "0/0", "1/1", 2)
        self.assertEqual(r.verdict, INCOMPATIBLE)
        self.assertEqual((r.total, r.explainable), (3, 0))
        # 同一结构孩子 0/0 时只有父补 0/2 可解释 -> uncertain
        r2 = evaluate("./2", "0/0", "0/0", 2)
        self.assertEqual(r2.verdict, UNCERTAIN)
        self.assertEqual((r2.total, r2.explainable), (3, 1))

    def test_phase_does_not_change_result(self):
        r_slash = evaluate("0/1", "0/1", "1/0", 1)
        r_pipe = evaluate("0|1", "0|1", "1|0", 1)
        self.assertEqual(r_slash.verdict, COMPATIBLE)
        self.assertEqual((r_pipe.total, r_pipe.explainable),
                         (r_slash.total, r_slash.explainable))

    def test_low_dp_whole_call_unknown_ignores_known_allele(self):
        # 父 1/1 但 DP 低于阈值：整调用未知，候选为全部 3 种二倍体，
        # 原 GT 的 1/1 不再约束补全（候选里出现了 0/0、0/1）
        r = evaluate("1/1", "0/0", "1/1", 1, dps=(5, 30, 30))
        self.assertEqual(r.candidates["father"],
                         ((0, 0), (0, 1), (1, 1)))
        # 母 0/0 只能给 0：孩子 1/1 在 3 个父补全下都凑不出 -> 不一致，
        # 但已是“3 个补全全部失败”而非单点判断
        self.assertEqual(r.verdict, INCOMPATIBLE)
        self.assertEqual((r.total, r.explainable), (3, 0))
        # 同 GT 深度足够时只有 1 个组合、同样确定不一致
        r_ok = evaluate("1/1", "0/0", "1/1", 1)
        self.assertEqual(r_ok.verdict, INCOMPATIBLE)
        self.assertEqual((r_ok.total, r_ok.explainable), (1, 0))

        # 阈值翻转：孩子 0/1 时，低深度父补 0/1、1/1 可解释（2 组）、
        # 补 0/0 不可解释（1 组）-> uncertain；深度足够（父固定 1/1）
        # 时 1/1 + 0/0 -> 0/1 本是正常传递 -> compatible
        r_flip = evaluate("1/1", "0/0", "0/1", 1, dps=(5, 30, 30))
        self.assertEqual(r_flip.verdict, UNCERTAIN)
        self.assertEqual((r_flip.total, r_flip.explainable), (3, 2))
        r_flip_ok = evaluate("1/1", "0/0", "0/1", 1)
        self.assertEqual(r_flip_ok.verdict, COMPATIBLE)

    def test_dp_missing_treated_unknown(self):
        call = make_call("father", "F", "1/1", dp=None)
        self.assertFalse(call.trusted)
        self.assertEqual(call.trust_reason, "dp_missing")

    def test_single_dot_and_dotdot_fully_unknown(self):
        for token in (".", "./.", ".|."):
            r = evaluate(token, "0/0", "0/1", 1)
            self.assertEqual(r.candidates["father"],
                             ((0, 0), (0, 1), (1, 1)), token)

    def test_unsupported_non_diploid_listed_separately(self):
        r = evaluate("1", "0/0", "0/1", 1)
        self.assertEqual(r.verdict, UNSUPPORTED)
        self.assertEqual(r.unsupported_members, ("father",))
        self.assertEqual(r.total, 0)
        r2 = evaluate("0/0", "0/0", "1/2/3", 3)
        self.assertEqual(r2.unsupported_members, ("child",))
        # 低深度的单倍体不进 unsupported（整调用未知）
        r3 = evaluate("1", "0/0", "0/1", 1, dps=(5, 30, 30))
        self.assertEqual(r3.verdict, UNCERTAIN)


# --------------------------------------------------------------------------- #
# 家系表校验
# --------------------------------------------------------------------------- #
class FamilyTableValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _fam(self, rows):
        p = os.path.join(self.tmp, "fam.tsv")
        write_families(p, rows)
        return p

    def test_valid_siblings_and_multigeneration(self):
        p = self._fam([
            ("F1", "CH1", "FA", "MO"),
            ("F1", "CH2", "FA", "MO"),
            ("F2", "FA", "GF", "GM"),
        ])
        trios = load_families(p)
        self.assertEqual(len(trios), 3)

    def test_bad_header_rejected(self):
        p = os.path.join(self.tmp, "fam.tsv")
        with open(p, "w") as fh:
            fh.write("fam\tchild\tdad\tmom\n")
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p)
        self.assertIn("表头", str(ctx.exception))

    def test_missing_column_and_not_distinct(self):
        p = self._fam([("F1", "C", "", "M")])
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p)
        self.assertIn("四列均必填", str(ctx.exception))

        p2 = self._fam([("F1", "C", "C", "M")])
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p2)
        self.assertIn("两两不同", str(ctx.exception))

    def test_duplicate_and_conflicting_parents_rejected(self):
        p = self._fam([
            ("F1", "C", "F", "M"),
            ("F1", "C", "F", "M"),
        ])
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p)
        self.assertIn("完全重复", str(ctx.exception))

        p2 = self._fam([
            ("F1", "C", "F", "M"),
            ("F9", "C", "F2", "M"),
        ])
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p2)
        self.assertIn("矛盾", str(ctx.exception))
        self.assertIn("family", str(ctx.exception))

    def test_cycle_detected_with_chain(self):
        # A 的父母含 B，B 的父母含 A -> 环
        p = self._fam([
            ("F1", "A", "B", "X"),
            ("F1", "B", "A", "Y"),
        ])
        with self.assertRaises(TrioConfigError) as ctx:
            load_families(p)
        self.assertIn("环", str(ctx.exception))
        self.assertIn("A", str(ctx.exception))

    def test_missing_person_in_vcf(self):
        p = self._fam([("F1", "GHOST", "F", "M")])
        vcf = os.path.join(self.tmp, "x.vcf")
        write_vcf(vcf, locus_line("chr1", 1, "A", "G", "0/0", "0/0", "0/0"))
        samples = read_vcf_header(vcf)
        trios = load_families(p)
        with self.assertRaises(TrioConfigError) as ctx:
            from vcfpipe.trio import check_persons
            check_persons(trios, samples, vcf)
        self.assertIn("GHOST", str(ctx.exception))

    def test_duplicate_vcf_sample_columns(self):
        vcf = os.path.join(self.tmp, "dup.vcf")
        header = HEADER.replace("F\tM\tC", "F\tM\tF")
        write_vcf(vcf, locus_line("chr1", 1, "A", "G",
                                  "0/0", "0/0", "0/0"), header=header)
        with self.assertRaises(TrioConfigError) as ctx:
            read_vcf_header(vcf)
        self.assertIn("重复", str(ctx.exception))


# --------------------------------------------------------------------------- #
# VCF 记录校验：坏行带行号失败
# --------------------------------------------------------------------------- #
class LocusValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_")
        self.fam = os.path.join(self.tmp, "fam.tsv")
        write_families(self.fam, [("F1", "C", "F", "M")])

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run_body(self, body, min_dp=MIN_DP):
        vcf = os.path.join(self.tmp, "x.vcf")
        write_vcf(vcf, body)
        out = os.path.join(self.tmp, "out")
        return run_trio(vcf, self.fam, CHROMS, min_dp, out)

    def test_bad_gt_syntax_fails_with_line(self):
        body = locus_line("chr1", 1, "A", "G", "x/1", "0/0", "0/0")
        with self.assertRaises(TrioRecordError) as ctx:
            self._run_body(body)
        self.assertIn("GT", str(ctx.exception))
        # 头 4 行 + 位点 = 第 5 行
        self.assertIn(":5", str(ctx.exception))

    def test_out_of_range_allele_fails(self):
        # 单 ALT 却出现编号 2：不能当成另一位点的 1
        body = locus_line("chr1", 1, "A", "G", "0/2", "0/0", "0/0")
        with self.assertRaises(TrioRecordError) as ctx:
            self._run_body(body)
        self.assertIn("越界", str(ctx.exception))

    def test_negative_dp_fails(self):
        vcf = os.path.join(self.tmp, "neg.vcf")
        write_vcf(vcf, "chr1\t1\t.\tA\tG\t.\t.\t.\tGT:DP\t"
                       "0/0:-3\t0/0:30\t0/0:30\n")
        out = os.path.join(self.tmp, "out")
        with self.assertRaises(TrioRecordError) as ctx:
            run_trio(vcf, self.fam, CHROMS, MIN_DP, out)
        self.assertIn("DP", str(ctx.exception))

    def test_bad_pos_ref_alt_and_column_count(self):
        bad_pos = locus_line("chr1", 0, "A", "G", "0/0", "0/0", "0/0")
        with self.assertRaises(TrioRecordError):
            self._run_body(bad_pos)

        bad_ref = locus_line("chr1", 1, "N", "G", "0/0", "0/0", "0/0")
        with self.assertRaises(TrioRecordError):
            self._run_body(bad_ref)

        bad_alt = locus_line("chr1", 1, "A", "G,G", "0/0", "0/0", "0/0")
        with self.assertRaises(TrioRecordError):
            self._run_body(bad_alt)

        five_alts = locus_line("chr1", 1, "A", "G,T,C,AA,AC",
                               "0/0", "0/0", "0/0")
        with self.assertRaises(TrioRecordError) as ctx:
            self._run_body(five_alts)
        self.assertIn("上限 4", str(ctx.exception))

        four_alts = locus_line("chr1", 1, "A", "G,T,C,AA",
                               "0/3", "1/2", "0/1")
        report = self._run_body(four_alts)
        self.assertEqual(report["selected_lines"], 1)

    def test_field_count_mismatch_fails(self):
        vcf = os.path.join(self.tmp, "fc.vcf")
        write_vcf(vcf, "chr1\t1\t.\tA\tG\t.\t.\t.\tGT:DP\t"
                       "0/0:30:99\t0/0:30\t0/0:30\n")
        out = os.path.join(self.tmp, "out")
        with self.assertRaises(TrioRecordError) as ctx:
            run_trio(vcf, self.fam, CHROMS, MIN_DP, out)
        self.assertIn("字段数", str(ctx.exception))

    def test_non_whitelist_only_skipped(self):
        body = locus_line("chrX", 1, "A", "G", "0/1", "0/1", "1/1")
        report = self._run_body(body)
        self.assertEqual(report["selected_lines"], 0)
        self.assertEqual(report["skipped_non_whitelist"], 1)
        self.assertEqual(report["detail_rows"], 0)


# --------------------------------------------------------------------------- #
# 端到端：演示数据、gzip、摘要对账与失败不发布
# --------------------------------------------------------------------------- #
class EndToEndDemo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_e2e_")
        cls.out = os.path.join(cls.tmp, "screen")
        cls.report = run_trio(DEMO_VCF, DEMO_FAM, CHROMS, MIN_DP, cls.out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _details(self):
        with open(os.path.join(self.out, "trio_details.tsv"),
                  encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_status_success_and_outputs(self):
        with open(os.path.join(self.out, "trio_status.json"),
                  encoding="utf-8") as fh:
            status = json.load(fh)
        self.assertEqual(status["status"], "success")
        # chrX 行被跳过：12 个 body 行中选中 11
        self.assertEqual(status["skipped_non_whitelist"], 1)
        self.assertEqual(status["selected_lines"], 11)

    def test_demo_scenario_verdicts(self):
        rows = {(r["line_no"], r["child"]): r for r in self._details()}

        def verdict(line_no, child):
            return rows[(str(line_no), child)]["verdict"]

        # 普通传递
        self.assertEqual(verdict(5, "CH1"), COMPATIBLE)
        # 父母全参考孩子含 ALT
        self.assertEqual(verdict(6, "CH1"), INCOMPATIBLE)
        # 父母分别贡献两个不同 ALT（多等位整体）
        self.assertEqual(verdict(7, "CH1"), COMPATIBLE)
        # 部分缺失但确定不一致（父 0/. 母 0/0 孩子 1/1）
        self.assertEqual(verdict(8, "CH1"), INCOMPATIBLE)
        # 部分缺失且全部补全一致（父母 0/1 孩子 0/.）
        self.assertEqual(verdict(9, "CH1"), COMPATIBLE)
        # 可行与不可行并存（130 位点 CH2）
        self.assertEqual(verdict(8, "CH2"), UNCERTAIN)
        # 低 DP（160 位点=第 11 行，父 1/1 DP5）：原 GT 保留、整调用未知
        low1 = rows[("11", "CH1")]
        self.assertEqual(low1["father_gt"], "1/1")
        self.assertEqual(low1["father_trust"], "dp_below_min")
        # 非二倍体单列（180 位点=第 13 行）
        self.assertEqual(verdict(13, "CH1"), UNSUPPORTED)
        # 相位无关（190 位点=第 14 行）
        self.assertEqual(verdict(14, "CH1"), COMPATIBLE)

    def test_low_dp_threshold_flip(self):
        # 同一 160 位点（VCF 第 11 行），父 GT 1/1、DP=5：
        # min-dp=10 父整调用未知 -> 候选含全部 3 种二倍体；
        #   CH2 孩子 0/0：父补 0/0、0/1 可解释、1/1 不可解释 -> uncertain
        #   CH1 孩子 1/1：母 0/0 给不出 1，3 个补全全失败 -> incompatible
        # min-dp=5  父可信固定 1/1：CH2 变为确定 compatible（父 0 给不出，
        #   实际 1/1+0/0 无法凑 0/0 -> incompatible），结论随阈值变化
        with open(os.path.join(self.out, "trio_details.tsv"),
                  newline="") as fh:
            high = {(r["line_no"], r["child"]): r
                    for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(high[("11", "CH2")]["verdict"], UNCERTAIN)
        self.assertEqual(high[("11", "CH2")]["father_trust"],
                         "dp_below_min")
        self.assertEqual(high[("11", "CH1")]["verdict"], INCOMPATIBLE)

        out_lo = os.path.join(self.tmp, "screen_lo")
        rep_lo = run_trio(DEMO_VCF, DEMO_FAM, CHROMS, 5, out_lo)
        with open(os.path.join(out_lo, "trio_details.tsv"),
                  newline="") as fh:
            low = {(r["line_no"], r["child"]): r
                   for r in csv.DictReader(fh, delimiter="\t")}
        # 阈值降到 5 后父亲可信：CH2 孩子 0/0 无法由 1/1 + 0/0 凑出
        self.assertEqual(low[("11", "CH2")]["verdict"], INCOMPATIBLE)
        self.assertEqual(low[("11", "CH2")]["father_trust"], "ok")
        # 同一行结论从 uncertain 翻成 incompatible
        self.assertNotEqual(low[("11", "CH2")]["verdict"],
                            high[("11", "CH2")]["verdict"])
        self.assertGreater(rep_lo["detail_rows"], 0)

    def test_shared_parents_siblings_differ(self):
        rows = {(r["line_no"], r["child"]): r for r in self._details()}
        # 110 位点：CH1 含 ALT 不一致，CH2 全参考一致（共享 FA/MO）
        self.assertEqual(rows[("6", "CH1")]["verdict"], INCOMPATIBLE)
        self.assertEqual(rows[("6", "CH2")]["verdict"], COMPATIBLE)

    def test_summary_matches_details(self):
        details = self._details()
        with open(os.path.join(self.out, "trio_summary.tsv"),
                  newline="") as fh:
            summary = list(csv.DictReader(fh, delimiter="\t"))
        per_trio = {}
        for r in details:
            key = (r["family"], r["child"])
            per_trio.setdefault(key, {COMPATIBLE: 0, INCOMPATIBLE: 0,
                                      UNCERTAIN: 0, UNSUPPORTED: 0})
            per_trio[key][r["verdict"]] += 1
        for row in summary:
            if row["family"] == "*TOTAL*":
                continue
            got = per_trio[(row["family"], row["child"])]
            for k in (COMPATIBLE, INCOMPATIBLE, UNCERTAIN, UNSUPPORTED):
                self.assertEqual(int(row[k]), got[k],
                                 f"{row['family']}/{row['child']} {k}")
            self.assertEqual(int(row["total"]), sum(got.values()))
        total_row = next(r for r in summary if r["family"] == "*TOTAL*")
        for k in (COMPATIBLE, INCOMPATIBLE, UNCERTAIN, UNSUPPORTED):
            self.assertEqual(int(total_row[k]),
                             sum(v[k] for v in per_trio.values()))

    def test_gzip_input_reads_identically(self):
        gz_path = os.path.join(self.tmp, "demo.vcf.gz")
        with open(DEMO_VCF, "rb") as src, gzip.open(gz_path, "wb") as dst:
            dst.write(src.read())
        out_gz = os.path.join(self.tmp, "screen_gz")
        run_trio(gz_path, DEMO_FAM, CHROMS, MIN_DP, out_gz)
        with open(os.path.join(out_gz, "trio_details.tsv"), newline="") as fh:
            gz_rows = list(csv.reader(fh, delimiter="\t"))
        with open(os.path.join(self.out, "trio_details.tsv"),
                  newline="") as fh:
            plain_rows = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(gz_rows, plain_rows)


class FailureNoPublish(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_fail_")
        self.fam = os.path.join(self.tmp, "fam.tsv")
        write_families(self.fam, [("F1", "C", "F", "M")])

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_bad_line_after_good_lines_leaves_no_output(self):
        good = locus_line("chr1", 1, "A", "G", "0/0", "0/1", "0/1")
        bad = locus_line("chr1", 2, "C", "T", "0/9", "0/0", "0/0")
        vcf = os.path.join(self.tmp, "halfbad.vcf")
        write_vcf(vcf, good + bad)
        out = os.path.join(self.tmp, "out")
        with self.assertRaises(TrioRecordError) as ctx:
            run_trio(vcf, self.fam, CHROMS, MIN_DP, out)
        self.assertIn("越界", str(ctx.exception))
        # 不发布输出目录，也不残留 staging
        self.assertFalse(os.path.exists(out))
        self.assertEqual(
            [n for n in os.listdir(self.tmp) if n.startswith(".")], [])

    def test_existing_nonempty_out_dir_rejected(self):
        vcf = os.path.join(self.tmp, "ok.vcf")
        write_vcf(vcf, locus_line("chr1", 1, "A", "G",
                                  "0/0", "0/1", "0/1"))
        out = os.path.join(self.tmp, "out")
        os.makedirs(out)
        with open(os.path.join(out, "old.txt"), "w") as fh:
            fh.write("x")
        with self.assertRaises(TrioConfigError):
            run_trio(vcf, self.fam, CHROMS, MIN_DP, out)


# --------------------------------------------------------------------------- #
# 独立穷举核对器（子进程，不导入产品补全/传递函数）
# --------------------------------------------------------------------------- #
class IndependentExhaustiveCheck(unittest.TestCase):
    def test_checker_matches_demo_details(self):
        tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_check_")
        try:
            out = os.path.join(tmp, "screen")
            run_trio(DEMO_VCF, DEMO_FAM, CHROMS, MIN_DP, out)
            proc = subprocess.run(
                [sys.executable, "-I", CHECKER,
                 "--details", os.path.join(out, "trio_details.tsv"),
                 "--vcf", DEMO_VCF, "--min-dp", str(MIN_DP)],
                capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("独立穷举核对通过", proc.stdout)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_checker_independent_enumeration_spot_check(self):
        """核对器内部枚举对若干含缺失多等位组合的手算值。"""
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "trio_exhaustive_check", CHECKER)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        # 父 ./2、母 0/0、孩子 1/1（2 ALT）：3 组合 0 可解释
        verdict, cand, total, expl = mod.evaluate(
            "./2", "0/0", "1/1", 30, 30, 30, MIN_DP, 2)
        self.assertEqual(verdict, "incompatible")
        self.assertEqual((total, expl), (3, 0))
        # 全未知三人在二等位位点：每人 3 种补全，3^3 = 27 组合；
        # 并非每组父母补全都能凑出固定孩子，故为 uncertain
        verdict, cand, total, expl = mod.evaluate(
            ".", ".", ".", 30, 30, 30, MIN_DP, 1)
        self.assertEqual(total, 27)
        self.assertEqual(verdict, "uncertain")
        self.assertGreater(expl, 0)
        self.assertLess(expl, total)
        # 低深度整调用未知：候选不受原 GT 1/1 约束，仍为全部 3 种
        self.assertEqual(len(mod.complete_unordered(None, 1)), 3)
        # brute_force 统一返回四元组；低深度时原 GT 不约束补全
        self.assertEqual(
            mod.brute_force("1/1", 5, MIN_DP, 1),
            ("ok", [(0, 0), (0, 1), (1, 1)], 0, 0))
        # 非二倍体被独立判为 unsupported
        self.assertEqual(
            mod.brute_force("1", 30, MIN_DP, 1),
            ("unsupported", [], 0, 0))


# --------------------------------------------------------------------------- #
# CLI 子进程回归（退出码与既有命令可用性）
# --------------------------------------------------------------------------- #
class CliSubprocess(unittest.TestCase):
    def test_trio_cli_exit_codes(self):
        tmp = tempfile.mkdtemp(prefix="vcfpipe_trio_cli_")
        try:
            out = os.path.join(tmp, "screen")
            ok = subprocess.run(
                [sys.executable, "-I", RUNNER, "trio",
                 "--vcf", DEMO_VCF, "--families", DEMO_FAM,
                 "--chroms", ",".join(CHROMS), "--min-dp", str(MIN_DP),
                 "--out-dir", out],
                capture_output=True, text=True)
            self.assertEqual(ok.returncode, 0, ok.stderr)
            self.assertTrue(os.path.isfile(
                os.path.join(out, "trio_status.json")))

            # 配置错误 -> 退出码 2
            bad_fam = os.path.join(tmp, "bad.tsv")
            write_families(bad_fam, [("F1", "C", "C", "M")])
            cfg = subprocess.run(
                [sys.executable, "-I", RUNNER, "trio",
                 "--vcf", DEMO_VCF, "--families", bad_fam,
                 "--chroms", ",".join(CHROMS), "--min-dp", str(MIN_DP),
                 "--out-dir", os.path.join(tmp, "x")],
                capture_output=True, text=True)
            self.assertEqual(cfg.returncode, 2)
            self.assertIn("两两不同", cfg.stderr)

            # 坏记录 -> 退出码 1
            good_fam = os.path.join(tmp, "fam.tsv")
            write_families(good_fam, [("F1", "C", "F", "M")])
            bad_vcf = os.path.join(tmp, "bad.vcf")
            write_vcf(bad_vcf, locus_line("chr1", 1, "A", "G",
                                          "0/9", "0/0", "0/0"))
            rec = subprocess.run(
                [sys.executable, "-I", RUNNER, "trio",
                 "--vcf", bad_vcf, "--families", good_fam,
                 "--chroms", ",".join(CHROMS), "--min-dp", str(MIN_DP),
                 "--out-dir", os.path.join(tmp, "y")],
                capture_output=True, text=True)
            self.assertEqual(rec.returncode, 1)
            self.assertIn("越界", rec.stderr)
            # 坏记录后不发布输出目录
            self.assertFalse(os.path.exists(os.path.join(tmp, "y")))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_existing_commands_still_available(self):
        proc = subprocess.run(
            [sys.executable, "-I", RUNNER, "--help"],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0)
        for name in ("run", "demux", "query", "assoc", "trio"):
            self.assertIn(name, proc.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
