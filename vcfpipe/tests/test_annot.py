"""annot 注释流程的手算回归与独立核对。

期望值来自两条互不复用的路径：

1. **手算/设计值**：测试里的人造小参考序列与密码子按设计直接写出，
   不调用产品的剪接或翻译函数；
2. **独立核对器**：``tools/annot_independent_check.py`` 自己解析
   FASTA/TSV、自己拼整条突变 CDS、自己翻译，逐行比对产品输出。

覆盖任务书要求的核对点：正/负链跨接缝密码子、同位点在不同异构体
中的不同后果、双碱基同时改、跨两个密码子、终止获得/丢失、CDS 端点、
外显子末端、跨内含子、无相交异构体、全非编码异构体、多 ALT、REF
错配与非法转录本定义、TSV/JSONL 两种格式、排序与流式写出。
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "src")
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "annot_independent_check.py")
DATA = os.path.join(ROOT, "data", "annot")
FA = os.path.join(ROOT, "data", "ref", "annot_syn.fa")
TXP = os.path.join(DATA, "transcripts_demo.tsv")
VAR = os.path.join(DATA, "variants_demo.tsv")
REF_MISMATCH = os.path.join(DATA, "regression_ref_mismatch.tsv")
BAD_TXP = os.path.join(DATA, "regression_bad_transcripts.tsv")

sys.path.insert(0, SRC)

from vcfpipe.annot import (  # noqa: E402
    AnnotInputError, HEADER, load_fasta, load_transcripts,
    reverse_complement, run_annot, translate_codon,
)


def cli(*args):
    return subprocess.run([sys.executable, "-I", RUNNER, *args],
                          cwd=ROOT, capture_output=True, text=True,
                          timeout=120)


def parse_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TmpCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="annot_test_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def out(self, name="out"):
        return os.path.join(self.tmp, name)

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return path


# --------------------------------------------------------------------------- #
# 1. 设计好的微型参考：只断言产品输出，不借用其剪接/翻译
# --------------------------------------------------------------------------- #
class HandCalcMicro(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="annot_micro_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return path

    """外显子 1 贡献 2 个编码碱基、外显子 2 贡献 4 个：首个密码子跨接缝。"""

    def build(self, exons, cds_start, cds_end, seq):
        fasta = self.write("r.fa", ">c1\n" + seq + "\n")
        starts = ",".join(str(s) for s, _ in exons)
        ends = ",".join(str(e) for _, e in exons)
        txp = os.path.join(self.tmp, "t.tsv")
        with open(txp, "w", encoding="utf-8", newline="") as fh:
            fh.write("transcript\tgene\tchrom\tstrand\texon_starts\t"
                     "exon_ends\tcds_start\tcds_end\n"
                     f"T\tG\tc1\t+\t{starts}\t{ends}\t"
                     f"{cds_start}\t{cds_end}\n")
        return load_fasta(fasta), txp

    PAD = "TTTT"

    def test_seam_codon_crosses_two_exons(self):
        # 外显子1 贡献 1 个编码碱基 "A"，外显子2 贡献 5 个 "TGCCC"。
        # 拼成 "ATGCCC"：首个密码子 ATG 横跨接缝，且两段自身都不整除 3。
        seq = "GG" + "A" + "GGG" + "TGCCC" + "TT"
        ref, txp = self.build([(0, 3), (6, 11)], 2, 11, seq)
        t = load_transcripts(txp, ref)[0]
        self.assertEqual(t.cds, "ATGCCC")
        self.assertEqual(t.coords, (3, 7, 8, 9, 10, 11))
        self.assertEqual(t.cds_length % 3, 0)

    def test_first_codon_spans_seam(self):
        # 外显子1 编码片段贡献 2 碱基 "AT"，外显子2 贡献 4 碱基 "GCCC"
        # 拼成 "ATGCCC"；首个密码子 ATG 横跨两段（阅读框贯穿接缝）。
        seq = "GGAT" + "AAAAAA" + "GCCC"     # [0,4) 与 [10,14)
        ref, txp = self.build([(0, 4), (10, 14)], 2, 14, seq)
        t = load_transcripts(txp, ref)[0]
        self.assertEqual(t.cds, "ATGCCC")
        self.assertEqual(t.coords, (3, 4, 11, 12, 13, 14))
        self.assertEqual(t.cds_length % 3, 0)
        # 反向互补同一条：负链拼接顺序相反，坐标按翻译方向递减
        txp_neg = self.write("tneg.tsv",
            "transcript\tgene\tchrom\tstrand\texon_starts\t"
            "exon_ends\tcds_start\tcds_end\n"
            "TN\tG\tc1\t-\t0,10\t4,14\t2\t14\n")
        tn = load_transcripts(txp_neg, ref)[0]
        self.assertEqual(tn.cds, reverse_complement(ref.sequences["c1"][10:14])
                         + reverse_complement(ref.sequences["c1"][2:4]))
        self.assertEqual(tn.coords, (14, 13, 12, 11, 4, 3))


# --------------------------------------------------------------------------- #
# 2. 演示材料：逐条手算/独立核对
# --------------------------------------------------------------------------- #
class DemoAnnotations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="annot_demo_")
        cls.summary = run_annot(FA, TXP, VAR, cls.tmp, "tsv")
        cls.rows = parse_tsv(os.path.join(cls.tmp, "annotation.tsv"))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def find(self, vid, alt_index, transcript=None):
        for r in self.rows:
            if r["id"] != vid or int(r["alt_index"]) != alt_index:
                continue
            if transcript is None or r["transcript"] == transcript:
                return r
        self.fail(f"缺少 {vid} ALT#{alt_index} {transcript}")

    # -- 跨接缝 ---------------------------------------------------------- #
    def test_positive_seam_codon(self):
        r = self.find("V_CROSS_POS", 1, "TXP_POS1")
        self.assertEqual(r["region"], "coding")
        self.assertEqual(r["codon_number"], "7")
        self.assertEqual((r["codon_pos_1"], r["codon_pos_2"],
                          r["codon_pos_3"]), ("29", "30", "91"))
        self.assertEqual((r["codon_ref"], r["codon_alt"]), ("TAC", "TGC"))
        self.assertEqual((r["aa_ref"], r["aa_alt"]), ("Y", "C"))
        self.assertEqual(r["consequence"], "missense")

    def test_negative_seam_codon_coords_in_translation_order(self):
        r = self.find("V_NEG_SEAM", 1, "TXP_NEG1")
        self.assertEqual(r["strand"], "-")
        self.assertEqual(r["codon_number"], "8")
        # 负链密码子三个坐标按翻译方向递减：61, 40, 39
        self.assertEqual((r["codon_pos_1"], r["codon_pos_2"],
                          r["codon_pos_3"]), ("61", "40", "39"))
        self.assertEqual((r["codon_ref"], r["codon_alt"]), ("TGC", "TAC"))
        self.assertEqual((r["aa_ref"], r["aa_alt"]), ("C", "Y"))
        self.assertEqual(r["consequence"], "missense")

    # -- 同一点、不同异构体 ---------------------------------------------- #
    def test_same_locus_different_isoforms(self):
        coding = self.find("V_CODING_VS_NONCODING", 1, "TXP_POS1")
        nc = self.find("V_CODING_VS_NONCODING", 1, "TXP_NONCODING")
        utr = self.find("V_CODING_VS_NONCODING", 1, "TXP_POS1_UTR")
        self.assertEqual(coding["region"], "coding")
        self.assertEqual(coding["consequence"], "missense")
        self.assertEqual(nc["region"], "noncoding_exon")
        self.assertEqual(nc["codon_ref"], "")
        self.assertEqual(nc["consequence"], "")
        # 同一基因的另一个异构体把同一位点归到 5' 非编码片段
        self.assertEqual(utr["region"], "utr")
        self.assertEqual(utr["consequence"], "")
        # 三行同属一条输入记录的一个 ALT，按 transcript 字典序排列
        self.assertLess(self.rows.index(nc), self.rows.index(coding))
        self.assertLess(self.rows.index(coding), self.rows.index(utr))

    def test_multi_alt_same_locus_independent(self):
        # V_CROSS_POS 同一位点三个 ALT：各自展开、后果互不干扰
        rows = [r for r in self.rows
                if r["id"] == "V_CROSS_POS" and r["transcript"] == "TXP_POS1"]
        self.assertEqual([r["alt"] for r in rows], ["G", "T", "C"])
        self.assertEqual([r["codon_alt"] for r in rows],
                         ["TGC", "TTC", "TCC"])
        self.assertEqual([r["aa_alt"] for r in rows], ["C", "F", "S"])
        self.assertTrue(all(r["consequence"] == "missense" for r in rows))

    # -- 双碱基同时改 ---------------------------------------------------- #
    def test_double_base_requires_joint_change(self):
        r = self.find("V_DOUBLE_BASE", 1, "TXP_POS1")
        self.assertEqual(r["codon_number"], "15")
        self.assertEqual((r["codon_ref"], r["codon_alt"]), ("CGA", "AGC"))
        self.assertEqual((r["aa_ref"], r["aa_alt"]), ("R", "S"))
        self.assertEqual(r["consequence"], "missense")
        # 两个单碱基替换各自同义：A 单独改与 G 单独改都不改变 R
        ref = load_fasta(FA)
        txp = {t.name: t for t in load_transcripts(TXP, ref)}["TXP_POS1"]
        codon = txp.cds[42:45]
        self.assertEqual(codon, "CGA")
        self.assertEqual(translate_codon("AGA"), "R")
        self.assertEqual(translate_codon("CGG"), "R")
        self.assertNotEqual(translate_codon("AGC"), "R")

    # -- 一段替换跨两个密码子 -------------------------------------------- #
    def test_edit_spanning_two_codons(self):
        rows = [r for r in self.rows
                if r["id"] == "V_TWO_CODONS" and r["transcript"] == "TXP_POS1"]
        self.assertEqual([r["codon_number"] for r in rows], ["14", "15"])
        self.assertEqual((rows[0]["codon_ref"], rows[0]["codon_alt"]),
                         ("AGA", "AGT"))
        self.assertEqual((rows[1]["codon_ref"], rows[1]["codon_alt"]),
                         ("CGA", "TAA"))
        self.assertEqual(rows[0]["aa_ref"], "R")
        self.assertEqual(rows[0]["aa_alt"], "S")
        self.assertEqual(rows[1]["consequence"], "stop_gained")

    # -- 终止获得/丢失/同义 ---------------------------------------------- #
    def test_stop_gained(self):
        r = self.find("V_STOP_GAINED", 1, "TXP_POS1")
        self.assertEqual((r["codon_ref"], r["codon_alt"]), ("TGG", "TGA"))
        self.assertEqual((r["aa_ref"], r["aa_alt"]), ("W", "*"))
        self.assertEqual(r["consequence"], "stop_gained")

    def test_stop_lost(self):
        r = self.find("V_STOP_LOST", 1, "TXP_POS1")
        self.assertEqual((r["codon_ref"], r["codon_alt"]), ("TAA", "CAA"))
        self.assertEqual((r["aa_ref"], r["aa_alt"]), ("*", "Q"))
        self.assertEqual(r["consequence"], "stop_lost")

    def test_synonymous_and_stop_to_stop(self):
        r = self.find("V_SYNONYMOUS", 1, "TXP_POS1")
        self.assertEqual(r["consequence"], "synonymous")
        # * -> * 也算 synonymous
        self.assertEqual(
            run_class("TAA", "TAG"), "synonymous")

    # -- 区域类别 -------------------------------------------------------- #
    def test_cds_endpoints(self):
        start = self.find("V_CDS_START", 1, "TXP_POS1")
        self.assertEqual((start["region"], start["cds_position"],
                          start["codon_number"]), ("coding", "1", "1"))
        self.assertEqual((start["codon_ref"], start["codon_alt"]),
                         ("ATG", "GTG"))
        end = self.find("V_CDS_END", 1, "TXP_POS1")
        # 250 是密码子 50 的第三碱基，对应 CDS 位置 150
        self.assertEqual((end["region"], end["cds_position"],
                          end["codon_number"]), ("coding", "150", "50"))
        self.assertEqual(end["consequence"], "synonymous")

    def test_exon_end_coding_last_base(self):
        r = self.find("V_EXON_END_CODING", 1, "TXP_POS1")
        self.assertEqual(r["region"], "coding")
        self.assertEqual(r["codon_number"], "37")
        self.assertEqual((r["codon_pos_1"], r["codon_pos_2"],
                          r["codon_pos_3"]), ("179", "180", "211"))

    def test_boundary_categories(self):
        for vid, expected in [
            ("V_BOUNDARY_EXON", "boundary_unsupported"),
            ("V_BOUNDARY_CODING", "boundary_unsupported"),
            ("V_BOUNDARY_SPAN", "boundary_unsupported"),
            ("V_INTRON", "intron"),
            ("V_UTR5", "utr"),
            ("V_UTR3", "utr"),
        ]:
            r = self.find(vid, 1, "TXP_POS1")
            self.assertEqual(r["region"], expected, vid)
            self.assertEqual(r["codon_number"], "", vid)

    def test_intergenic_has_empty_transcript(self):
        r = self.find("V_INTERGENIC", 1)
        self.assertEqual(r["region"], "intergenic")
        self.assertEqual(r["transcript"], "")
        self.assertEqual(r["gene"], "")
        self.assertEqual(r["strand"], "")

    def test_unsupported_edit_keeps_intersecting_transcripts(self):
        rs = [r for r in self.rows if r["id"] == "V_UNSUPPORTED"]
        self.assertEqual([r["transcript"] for r in rs],
                         ["TXP_NONCODING", "TXP_POS1", "TXP_POS1_UTR"])
        self.assertTrue(all(r["region"] == "unsupported_edit" for r in rs))
        # 混合 ALT：编码 + 长度不等
        mixed = [r for r in self.rows if r["id"] == "V_ALT_MIXED"]
        self.assertEqual(mixed[0]["region"], "noncoding_exon")  # ALT#1
        self.assertEqual(mixed[-1]["region"], "unsupported_edit")  # ALT#2

    def test_negative_intron_and_utr(self):
        self.assertEqual(self.find("V_NEG_INTRON", 1, "TXP_NEG1")["region"],
                         "intron")
        self.assertEqual(self.find("V_NEG_UTR", 1, "TXP_NEG1")["region"],
                         "utr")

    # -- “只碰到边缘最后一个碱基”也按整个区间判定 ---------------------- #
    def test_edge_last_base_still_uses_full_span(self):
        # REF 的最后一个碱基正好是外显子末碱基 / 编码首碱基：仍按
        # 完整 REF 区间判定，不看 pos 锚点。
        from vcfpipe.annot import Variant, annotate_variant
        txps = load_transcripts(TXP, load_fasta(FA))
        # 外显子 [0,30)：REF 区间 [28,30) 仍在外显子内 -> coding
        v = Variant("E_LAST_CODING", "chrP", 29, "AA", ("GG",), 3)
        rows = annotate_variant(v, txps)
        pos1 = [r for r in rows if r.transcript == "TXP_POS1"]
        self.assertEqual(pos1[0].region, "coding")
        # REF 区间 [29,31)：末碱基伸进内含子 -> boundary（哪怕只看 pos=30
        # 会误判成编码）
        v2 = Variant("E_LAST_INTRON", "chrP", 30, "AT", ("GG",), 4)
        rows2 = annotate_variant(v2, txps)
        self.assertEqual(
            [r.region for r in rows2 if r.transcript == "TXP_POS1"],
            ["boundary_unsupported"])
        # REF 区间 [8,10)：末碱基正好是非编码/编码交界前一个碱基 -> utr
        v3 = Variant("E_LAST_UTR", "chrP", 9, "TA", ("CG",), 5)
        self.assertEqual(
            [r.region for r in annotate_variant(v3, txps)
             if r.transcript == "TXP_POS1"], ["utr"])
        # REF 区间 [9,11)：跨过编码首碱基 -> boundary
        v4 = Variant("E_FIRST_CODING", "chrP", 10, "AC", ("GT",), 6)
        self.assertEqual(
            [r.region for r in annotate_variant(v4, txps)
             if r.transcript == "TXP_POS1"], ["boundary_unsupported"])

    # -- 排序与结构 ------------------------------------------------------ #
    def test_sorting_is_variant_order_alt_transcript_codon(self):
        keys = []
        for r in self.rows:
            keys.append((r["id"], int(r["alt_index"]),
                         r["transcript"] or "-",
                         int(r["codon_number"]) if r["codon_number"] else -1))
        # 输入顺序按首次出现建立
        order = []
        for k in keys:
            if k[0] not in order:
                order.append(k[0])
        self.assertEqual(order, [
            "V_CROSS_POS", "V_NEG_SEAM", "V_STOP_LOST", "V_STOP_GAINED",
            "V_SYNONYMOUS", "V_DOUBLE_BASE", "V_TWO_CODONS",
            "V_CODING_VS_NONCODING", "V_UTR5", "V_UTR3", "V_INTRON",
            "V_CDS_START", "V_CDS_END", "V_EXON_END_CODING",
            "V_BOUNDARY_EXON", "V_BOUNDARY_CODING", "V_BOUNDARY_SPAN",
            "V_INTERGENIC", "V_UNSUPPORTED", "V_ALT_MIXED",
            "V_NEG_INTRON", "V_NEG_UTR"])
        for vid in order:
            block = [k for k in keys if k[0] == vid]
            self.assertEqual(block, sorted(block), vid)

    def test_all_regions_and_consequences_present(self):
        regions = {r["region"] for r in self.rows}
        self.assertEqual(regions, {
            "coding", "utr", "noncoding_exon", "intron", "intergenic",
            "boundary_unsupported", "unsupported_edit"})
        cons = {r["consequence"] for r in self.rows if r["consequence"]}
        self.assertEqual(cons, {"synonymous", "missense", "stop_gained",
                                "stop_lost"})

    def test_summary_counts(self):
        self.assertEqual(self.summary["explanations"], len(self.rows))
        self.assertEqual(sum(self.summary["region_counts"].values()),
                         len(self.rows))
        self.assertEqual(self.summary["inputs"]["transcripts"], 4)
        self.assertEqual(self.summary["inputs"]["variants"], 22)


def run_class(ref_codon, alt_codon):
    from vcfpipe.annot import classify_consequence
    return classify_consequence(translate_codon(ref_codon),
                                translate_codon(alt_codon))


# --------------------------------------------------------------------------- #
# 3. 独立核对器 + 两种输出格式
# --------------------------------------------------------------------------- #
class IndependentCheck(TmpCase):
    def test_independent_checker_matches_tsv(self):
        out = self.out()
        r = cli("annot", "--fasta", FA, "--transcripts", TXP,
                "--variants", VAR, "--out-dir", out)
        self.assertEqual(r.returncode, 0, r.stderr)
        c = subprocess.run(
            [sys.executable, "-I", CHECKER, "--fasta", FA,
             "--transcripts", TXP, "--variants", VAR, "--out-dir", out],
            cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(c.returncode, 0, c.stderr)
        self.assertIn("核对通过", c.stdout)

    def test_jsonl_format_equals_tsv(self):
        out_t = self.out("t")
        out_j = self.out("j")
        self.assertEqual(cli("annot", "--fasta", FA, "--transcripts", TXP,
                             "--variants", VAR, "--out-dir", out_t,
                             "--format", "tsv").returncode, 0)
        self.assertEqual(cli("annot", "--fasta", FA, "--transcripts", TXP,
                             "--variants", VAR, "--out-dir", out_j,
                             "--format", "jsonl").returncode, 0)
        tsv = parse_tsv(os.path.join(out_t, "annotation.tsv"))
        with open(os.path.join(out_j, "annotation.jsonl"),
                  encoding="utf-8") as fh:
            jsonl = [json.loads(line) for line in fh if line.strip()]
        self.assertEqual(len(tsv), len(jsonl))
        for a, b in zip(tsv, jsonl):
            self.assertEqual(list(b.keys()), list(HEADER))
            for key in HEADER:
                self.assertEqual(a[key], "" if b[key] is None else str(b[key]),
                                 key)
        c = subprocess.run(
            [sys.executable, "-I", CHECKER, "--fasta", FA,
             "--transcripts", TXP, "--variants", VAR, "--out-dir", out_j],
            cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(c.returncode, 0, c.stderr)


# --------------------------------------------------------------------------- #
# 4. 失败路径（非零退出、指出文件与行号）
# --------------------------------------------------------------------------- #
class FailureModes(TmpCase):
    def run_expect_fail(self, fasta, txp, var, needle):
        out = self.out()
        r = cli("annot", "--fasta", fasta, "--transcripts", txp,
                "--variants", var, "--out-dir", out)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn(needle, r.stderr)
        self.assertFalse(os.path.exists(os.path.join(out, "annotation.tsv")))
        return r

    def test_ref_mismatch(self):
        self.run_expect_fail(FA, TXP, REF_MISMATCH, "REF 校验失败")

    def test_bad_transcript_overlap(self):
        self.run_expect_fail(FA, BAD_TXP, VAR, "外显子")

    def test_missing_contig(self):
        bad = self.write("bad.tsv",
                         "id\tchrom\tpos\tref\talt\nV\tchrZ\t1\tA\tG\n")
        self.run_expect_fail(FA, TXP, bad, "chrZ")

    def test_duplicate_fasta_name(self):
        bad = self.write("dup.fa", ">chrP\nACGT\n>chrP\nACGT\n")
        self.run_expect_fail(bad, TXP, VAR, "重复")

    def test_duplicate_variant_id(self):
        bad = self.write("v.tsv",
                         "id\tchrom\tpos\tref\talt\n"
                         "V\tchrP\t5\tA\tG\nV\tchrP\t5\tA\tC\n")
        self.run_expect_fail(FA, TXP, bad, "重复")

    def test_bad_base(self):
        bad = self.write("v.tsv",
                         "id\tchrom\tpos\tref\talt\nV\tchrP\t5\tA\tN\n")
        self.run_expect_fail(FA, TXP, bad, "A/C/G/T")

    def test_bad_cds_length_not_divisible(self):
        # 编码长度 20 + 90 + 39 = 149，不整除 3；两端都在外显子上
        txp = self.write("t.tsv",
                         "transcript\tgene\tchrom\tstrand\texon_starts\t"
                         "exon_ends\tcds_start\tcds_end\n"
                         "T\tG\tchrP\t+\t0,210\t30,300\t10\t249\n")
        self.run_expect_fail(FA, txp, VAR, "不整除 3")

    def test_noncoding_allele_but_ref_ok(self):
        # 非编码异构体合法：cds_start == cds_end
        txp = self.write("t.tsv",
                         "transcript\tgene\tchrom\tstrand\texon_starts\t"
                         "exon_ends\tcds_start\tcds_end\n"
                         "T\tG\tchrP\t+\t0\t30\t20\t20\n")
        out = self.out()
        r = cli("annot", "--fasta", FA, "--transcripts", txp,
                "--variants", VAR, "--out-dir", out)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_empty_out_dir_required(self):
        out = self.out()
        os.makedirs(out)
        with open(os.path.join(out, "old.txt"), "w", encoding="utf-8") as fh:
            fh.write("x")
        r = cli("annot", "--fasta", FA, "--transcripts", TXP,
                "--variants", VAR, "--out-dir", out)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("非空", r.stderr)

    def test_direct_run_annot_raises_input_error(self):
        with self.assertRaises(AnnotInputError):
            run_annot(FA, TXP, REF_MISMATCH, self.out())


# --------------------------------------------------------------------------- #
# 5. 流式写出：逐变异 flush，不积攒到整批结束
# --------------------------------------------------------------------------- #
class StreamingWrite(TmpCase):
    def test_rows_visible_before_batch_finishes(self):
        import vcfpipe.annot as annot_mod
        out = self.out()
        seen = []
        original = annot_mod.AnnotationWriter.emit

        def spy(self, explanations):
            original(self, explanations)
            path = os.path.join(out, "annotation.tsv")
            if os.path.isfile(path):
                with open(path, encoding="utf-8") as fh:
                    seen.append(len(fh.readlines()))

        annot_mod.AnnotationWriter.emit = spy
        try:
            run_annot(FA, TXP, VAR, out)
        finally:
            annot_mod.AnnotationWriter.emit = original
        # 第一次 emit 后文件里就应该有表头 + 数据行
        self.assertGreater(seen[0], 1)

    def test_previous_entries_still_work(self):
        # 新能力加入后，原有入口至少仍可 --help；完整回归在各自测试文件。
        for sub in ("run", "demux", "query", "assoc", "trio", "annot"):
            r = cli(sub, "--help")
            self.assertEqual(r.returncode, 0, sub)


if __name__ == "__main__":
    unittest.main(verbosity=2)
