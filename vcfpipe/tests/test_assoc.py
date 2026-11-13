"""基因级携带者关联分析的手算回归与命令行回归。

覆盖任务书要求的可手算核对点：
* 同一个人同一基因多处变异不重复计人数；
* 已确定携带者不因其它位点未知被排除；
* 相邻区间合并后完整覆盖才可补参考（只覆盖锚点/缺一段都不行）；
* 显式缺失/坏调用不会被可判读区间掩码洗成参考；
* 重叠基因与删失变异跨区间边界时的归属；
* 用**独立于产品统计实现**的小表枚举核对双侧 Fisher；
* 用含相同概率的多基因例子核对 BH 校正与（校正值, 原始概率, 基因名）排序；
* 组别冲突与一组分母归零的**实际命令结果**（子进程）；
* 关联分析对已发布目录零改动。
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "src")
RUNNER = os.path.join(ROOT, "run_pipeline.py")
sys.path.insert(0, SRC)

from vcfpipe.assoc import (  # noqa: E402
    AssocInputError, aggregate_gene_call, benjamini_hochberg,
    classify_call, decide_locus, fisher_exact_two_sided, interval_covered,
    merge_intervals, odds_ratio, run_assoc, to_jsonable,
)
from vcfpipe.pipeline import run_pipeline  # noqa: E402

FA = os.path.join(ROOT, "data/ref/syn37.fa")
DICT = os.path.join(ROOT, "data/ref/contigs.json")
SHEET = os.path.join(ROOT, "data/input/sample_sheet.csv")
PRIO = os.path.join(ROOT, "data/input/caller_priority.json")
ASSOC_DATA = os.path.join(ROOT, "data", "assoc")
CHECKED_IN_DATASET = os.path.join(ROOT, "output", "syn37_v1")


def demo_file(name):
    return os.path.join(ASSOC_DATA, name)


# --------------------------------------------------------------------------- #
# 演示队列（data/assoc/ + 管线真实输出 syn37_v1）逐条手算核对
# --------------------------------------------------------------------------- #
class DemoCohortHandCalc(unittest.TestCase):
    """演示队列的全部期望值均按 README 的语义手工推导，不引用产品实现。"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_")
        run_pipeline(cls.tmp, FA, DICT, SHEET, PRIO, "v1")
        cls.dataset = os.path.join(cls.tmp, "syn37_v1")
        cls.report = run_assoc(
            cls.dataset, demo_file("labels_demo.tsv"),
            demo_file("genes_demo.bed"), demo_file("callable_demo.bed"), 10)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # ------------------------------------------------------------------ #
    # 工具
    # ------------------------------------------------------------------ #
    def entry(self, gene):
        for e in self.report["results"] + self.report["not_tested"]:
            if e["gene"] == gene:
                return e
        self.fail(f"报告中没有 {gene}")

    def decision(self, gene, sample):
        for d in self.entry(gene)["decisions"]:
            if d["sample"] == sample:
                return d
        self.fail(f"{gene} 缺少样本 {sample} 的判定明细")

    def locus(self, gene, sample, pos, alt):
        for l in self.decision(gene, sample)["loci"]:
            if l["pos"] == pos and l["alt"] == alt:
                return l
        self.fail(f"{gene}/{sample} 缺少位点 {pos}>{alt}")

    # ------------------------------------------------------------------ #
    # 候选位点归属（含重叠基因与跨边界删失）
    # ------------------------------------------------------------------ #
    def test_candidate_loci_assignment(self):
        expected = {
            # [0,12)：pos5、pos10 两个拆分位点
            "GENE_A": [(5, "G"), (10, "A"), (10, "T")],
            # [12,28)：pos13/15/18、pos28 插入 [27,28)，以及跨右边界
            # 的删失 [27,29)（27<28 相交）
            "GENE_B": [(13, "T"), (15, "A"), (18, "G"), (28, "GA"),
                       (28, "G")],
            "GENE_C": [(46, "C"), (48, "CT")],
            "GENE_D": [(63, "C"), (63, "T")],
            # [28,31)：插入 [27,28) 恰止于左边界不相交；删失 [27,29)
            # 跨左边界 -> 仅删失归入
            "GENE_E": [(28, "G")],
            # 重叠两行合并为 [4,14)：pos5、pos10 两个、pos13，各一次
            "GENE_G": [(5, "G"), (10, "A"), (10, "T"), (13, "T")],
        }
        for gene, want in expected.items():
            got = [(l["pos"], l["alt"])
                   for l in self.entry(gene)["candidate_loci"]]
            self.assertEqual(got, want, f"{gene} 候选位点不符")
        self.assertEqual(self.entry("GENE_F")["candidate_loci"], [])

    def test_variant_belongs_to_multiple_genes(self):
        # 删失 GA>G [27,29) 同时属于 GENE_B 与 GENE_E
        self.assertIn((28, "G"), [(l["pos"], l["alt"]) for l in
                                  self.entry("GENE_B")["candidate_loci"]])
        self.assertIn((28, "G"), [(l["pos"], l["alt"]) for l in
                                  self.entry("GENE_E")["candidate_loci"]])
        # pos13 [12,13) 同时属于 GENE_B 与 GENE_G
        self.assertIn((13, "T"), [(l["pos"], l["alt"]) for l in
                                  self.entry("GENE_G")["candidate_loci"]])
        # 删失的参考区间是 [27,29)（pos-1 起、长 len(ref)）
        locus = next(l for l in self.entry("GENE_E")["candidate_loci"]
                     if l["alt"] == "G")
        self.assertEqual((locus["ref_start0"], locus["ref_end0"]), (27, 29))

    # ------------------------------------------------------------------ #
    # 两组计数与实际分母
    # ------------------------------------------------------------------ #
    def test_counts_and_denominators(self):
        # gene: (case 携带/未携带/排除/分母, control 携带/未携带/排除/分母)
        expected = {
            "GENE_A": ((2, 0, 0, 2), (1, 0, 0, 1)),
            "GENE_B": ((2, 0, 0, 2), (0, 1, 0, 1)),
            "GENE_C": ((1, 1, 0, 2), (0, 1, 0, 1)),
            "GENE_D": ((1, 0, 1, 1), (0, 1, 0, 1)),
            "GENE_E": ((1, 0, 1, 1), (0, 1, 0, 1)),
            "GENE_F": ((0, 0, 2, 0), (0, 0, 1, 0)),
            "GENE_G": ((2, 0, 0, 2), (1, 0, 0, 1)),
        }
        for gene, (case, control) in expected.items():
            e = self.entry(gene)
            for group, want in (("case", case), ("control", control)):
                got = e["counts"][group]
                self.assertEqual(
                    (got["carrier"], got["non_carrier"], got["excluded"],
                     got["denominator"]), want,
                    f"{gene}/{group} 计数不符")

    # ------------------------------------------------------------------ #
    # 关键语义
    # ------------------------------------------------------------------ #
    def test_multi_variant_person_counts_once(self):
        # S1 在 GENE_C 上有两个有效携带位点（pos46 0/1、pos48 1/1），
        # 仍只计为一个携带者
        d = self.decision("GENE_C", "S1")
        self.assertEqual(d["gene_call"], "carrier")
        carriers = [l for l in d["loci"] if l["decision"] == "carrier"]
        self.assertEqual(len(carriers), 2)
        self.assertEqual(self.entry("GENE_C")["counts"]["case"]["carrier"], 1)

    def test_carrier_not_excluded_by_other_unknown_loci(self):
        # S1 在 GENE_A：pos5 有效携带；pos10 两个位点半缺失 -> 未知，
        # 但携带结论不受未知位点影响
        d = self.decision("GENE_A", "S1")
        self.assertEqual(d["gene_call"], "carrier")
        reasons = {l["reason"] for l in d["loci"]}
        self.assertIn("unknown_half_missing", reasons)
        self.assertIn("carrier_valid_alt_call", reasons)
        self.assertEqual(self.entry("GENE_A")["counts"]["case"]["excluded"], 0)

    def test_adjacent_intervals_merge_to_full_coverage(self):
        # S2 的可判读区间 [40,46)+[46,50) 合并为 [40,50)，完整覆盖
        # GENE_C 的 [45,47) 与 [47,48) -> 两个位点都补为参考 -> 未携带
        d = self.decision("GENE_C", "S2")
        self.assertEqual(d["gene_call"], "non_carrier")
        self.assertTrue(all(l["reason"] == "reference_imputed_by_coverage"
                            for l in d["loci"]))

    def test_partial_coverage_is_not_reference(self):
        # 把 S2 的可判读区间换成只覆盖锚点的 [45,46)：[45,47) 缺一段，
        # [47,48) 完全没覆盖 -> 两个位点都未知 -> S2 被排除
        tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_part_")
        try:
            path = os.path.join(tmp, "callable.bed")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("chr1\t0\t106\tS1\nchr1\t0\t106\tS3\n"
                         "chr1\t45\t46\tS2\n")
            report = run_assoc(
                self.dataset, demo_file("labels_demo.tsv"),
                demo_file("genes_demo.bed"), path, 10)
            gene_c = next(e for e in report["results"] + report["not_tested"]
                          if e["gene"] == "GENE_C")
            s2 = next(d for d in gene_c["decisions"] if d["sample"] == "S2")
            self.assertEqual(s2["gene_call"], "excluded")
            self.assertTrue(all(l["reason"] == "unknown_not_covered"
                                for l in s2["loci"]))
            self.assertEqual(gene_c["counts"]["case"]["excluded"], 1)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_missing_call_not_masked_by_coverage(self):
        # S1 可判读区间覆盖整个 chr1，但 pos63 A>T 是显式 ./.：
        # 不能被洗成参考；A>C 无记录且被覆盖 -> 补参考；
        # 无携带且非全参考 -> S1 从 GENE_D 分母排除
        d = self.decision("GENE_D", "S1")
        self.assertEqual(d["gene_call"], "excluded")
        by_alt = {l["alt"]: l for l in d["loci"]}
        self.assertEqual(by_alt["T"]["reason"], "unknown_fully_missing")
        self.assertEqual(by_alt["C"]["reason"],
                         "reference_imputed_by_coverage")

    def test_na_in_matrix_is_not_no_variant(self):
        # S2 在 GENE_E 的删失位点无记录且无覆盖 -> 未知 -> 排除，
        # 而不是当作没有变异
        d = self.decision("GENE_E", "S2")
        self.assertEqual(d["gene_call"], "excluded")
        self.assertEqual(d["loci"][0]["reason"], "unknown_not_covered")

    def test_low_dp_call_is_unknown(self):
        # min-dp 提到 20：S2 的 pos15（DP=17）与 S1 的 pos46/48
        # （DP=13/15）失效 -> 相应人被排除
        report = run_assoc(
            self.dataset, demo_file("labels_demo.tsv"),
            demo_file("genes_demo.bed"), demo_file("callable_demo.bed"), 20)
        gene_b = next(e for e in report["results"] + report["not_tested"]
                      if e["gene"] == "GENE_B")
        s2 = next(d for d in gene_b["decisions"] if d["sample"] == "S2")
        self.assertEqual(s2["gene_call"], "excluded")
        self.assertIn("unknown_low_dp", {l["reason"] for l in s2["loci"]})
        gene_c = next(e for e in report["results"] + report["not_tested"]
                      if e["gene"] == "GENE_C")
        s1 = next(d for d in gene_c["decisions"] if d["sample"] == "S1")
        self.assertEqual(s1["gene_call"], "excluded")
        # S1 的 pos28 插入 DP=20 达到阈值仍有效 -> GENE_B 携带不受影响
        s1_b = next(d for d in gene_b["decisions"] if d["sample"] == "S1")
        self.assertEqual(s1_b["gene_call"], "carrier")

    # ------------------------------------------------------------------ #
    # 统计结果（手算值）
    # ------------------------------------------------------------------ #
    def test_fisher_and_odds_ratio_hand_values(self):
        expected = {
            # (2,0,1,0)：边际下唯一可能的表 -> p=1；ad=bc=0 -> undefined
            "GENE_A": (Fraction(1), "undefined"),
            # (2,0,0,1)：P(2)=1/3，P(1)=2/3 -> p=1/3；bc=0, ad>0 -> inf
            "GENE_B": (Fraction(1, 3), "inf"),
            # (1,1,0,1)：p_obs=2/3，全部表 <= 2/3 -> p=1；bc=0 -> inf
            "GENE_C": (Fraction(1), "inf"),
            # (1,0,0,1)：两张表各 1/2 -> p=1；bc=0 -> inf
            "GENE_D": (Fraction(1), "inf"),
            "GENE_E": (Fraction(1), "inf"),
            "GENE_G": (Fraction(1), "undefined"),
        }
        for gene, (p, or_) in expected.items():
            e = self.entry(gene)
            self.assertEqual(e["status"], "tested")
            self.assertEqual(e["fisher_p_two_sided"], p, f"{gene} p 值不符")
            self.assertEqual(e["odds_ratio"], or_, f"{gene} OR 不符")

    def test_bh_and_stable_order(self):
        # 6 个被检验基因：p 最小的是 GENE_B(1/3)，其余全为 1；
        # BH 后全部为 1（6/1 * 1/3 = 2 被截断为 1）；
        # 排序 = (校正值, 原始概率, 基因名)
        self.assertEqual(self.report["genes_tested"], 6)
        order = [e["gene"] for e in self.report["results"]]
        self.assertEqual(order, ["GENE_B", "GENE_A", "GENE_C", "GENE_D",
                                 "GENE_E", "GENE_G"])
        self.assertTrue(all(e["bh_adjusted_p"] == 1
                            for e in self.report["results"]))
        self.assertEqual([e["rank"] for e in self.report["results"]],
                         [1, 2, 3, 4, 5, 6])

    def test_no_candidate_loci_gene_marked_not_tested(self):
        e = self.entry("GENE_F")
        self.assertEqual(e["status"], "no_candidate_loci")
        self.assertIn("候选位点", e["reason"])
        # 未检验项不进入检验数
        self.assertEqual(self.report["genes_tested"], 6)
        self.assertEqual(self.report["genes_total"], 7)

    def test_unlabeled_samples_not_in_denominator(self):
        # 只标注 S1(病例) 与 S3(对照)：S2 不参与任何分母
        tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_unlab_")
        try:
            path = os.path.join(tmp, "labels.tsv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("sample\tgroup\nS1\tcase\nS3\tcontrol\n")
            report = run_assoc(
                self.dataset, path, demo_file("genes_demo.bed"),
                demo_file("callable_demo.bed"), 10)
            self.assertEqual(report["analysis_samples"],
                             {"case": ["S1"], "control": ["S3"]})
            gene_b = next(e for e in report["results"]
                          if e["gene"] == "GENE_B")
            self.assertEqual({d["sample"] for d in gene_b["decisions"]},
                             {"S1", "S3"})
            self.assertEqual(gene_b["table"],
                             {"case_carrier": 1, "case_non_carrier": 0,
                              "control_carrier": 0, "control_non_carrier": 1})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# 小构件的单元手算
# --------------------------------------------------------------------------- #
class IntervalAndCallUnits(unittest.TestCase):
    def test_merge_intervals_unions_adjacent(self):
        self.assertEqual(merge_intervals([(40, 46), (46, 50)]),
                         [(40, 50)])
        self.assertEqual(merge_intervals([(46, 50), (40, 46), (41, 45)]),
                         [(40, 50)])
        # 有缺口的区间不合并
        self.assertEqual(merge_intervals([(0, 5), (7, 10)]),
                         [(0, 5), (7, 10)])

    def test_interval_covered_requires_full_containment(self):
        merged = [(40, 50)]
        self.assertTrue(interval_covered(merged, 45, 47))
        self.assertTrue(interval_covered(merged, 40, 50))
        # 只覆盖锚点（起点）不算
        self.assertFalse(interval_covered([(40, 46)], 45, 47))
        # 缺右端一段不算
        self.assertFalse(interval_covered([(46, 50)], 45, 47))
        self.assertFalse(interval_covered([], 45, 47))

    def test_classify_call_truth_table(self):
        # 有效携带 / 有效参考
        self.assertEqual(classify_call("0/1", 20, 10),
                         ("carrier", "carrier_valid_alt_call"))
        self.assertEqual(classify_call("1|0", 17, 10),
                         ("carrier", "carrier_valid_alt_call"))
        self.assertEqual(classify_call("1/1", 20, 10),
                         ("carrier", "carrier_valid_alt_call"))
        self.assertEqual(classify_call("0/0", 20, 10),
                         ("reference", "reference_valid_call"))
        # 半缺失 / 全缺失 / 非二倍体 / 等位越界
        self.assertEqual(classify_call("./1", 11, 10)[1],
                         "unknown_half_missing")
        self.assertEqual(classify_call("0/.", 30, 10)[1],
                         "unknown_half_missing")
        self.assertEqual(classify_call("./.", None, 10)[1],
                         "unknown_fully_missing")
        self.assertEqual(classify_call("1", 30, 10)[1],
                         "unknown_non_diploid")
        self.assertEqual(classify_call("0/1/1", 30, 10)[1],
                         "unknown_non_diploid")
        self.assertEqual(classify_call("2/2", 30, 10)[1],
                         "unknown_allele_out_of_range")
        # 低深度 / 深度缺失
        self.assertEqual(classify_call("0/1", 9, 10)[1], "unknown_low_dp")
        self.assertEqual(classify_call("0/1", None, 10)[1],
                         "unknown_dp_missing")
        # 阈值边界：DP == min_dp 有效
        self.assertEqual(classify_call("0/1", 10, 10)[0], "carrier")

    def test_decide_locus_record_wins_over_coverage(self):
        # 显式坏调用即使被可判读区间覆盖也不洗成参考
        self.assertEqual(decide_locus(("./.", None), True, 10),
                         ("unknown", "unknown_fully_missing"))
        self.assertEqual(decide_locus(("./1", 11), True, 10),
                         ("unknown", "unknown_half_missing"))
        # 无记录时才看覆盖
        self.assertEqual(decide_locus(None, True, 10),
                         ("reference", "reference_imputed_by_coverage"))
        self.assertEqual(decide_locus(None, False, 10),
                         ("unknown", "unknown_not_covered"))

    def test_aggregate_gene_call(self):
        self.assertEqual(aggregate_gene_call(["unknown", "carrier"]),
                         "carrier")
        self.assertEqual(aggregate_gene_call(["reference", "reference"]),
                         "non_carrier")
        self.assertEqual(aggregate_gene_call(["reference", "unknown"]),
                         "excluded")
        self.assertEqual(aggregate_gene_call([]), "excluded")


# --------------------------------------------------------------------------- #
# 独立于产品统计实现的小表枚举，核对双侧 Fisher
# --------------------------------------------------------------------------- #
def independent_fisher_two_sided(a, b, c, d):
    """与 vcfpipe.assoc 无关的独立实现：显式枚举所有同边际表。

    用纯整数阶乘算组合数、浮点求概率；逐一枚举左上格的所有合法
    取值（固定边际下 2x2 表只有一个自由度），把概率不超过观察表
    （含浮点容差）的表的概率相加。
    """
    from math import factorial

    def choose(n, k):
        return factorial(n) // (factorial(k) * factorial(n - k))

    r1, r2 = a + b, c + d
    c1 = a + c
    n = r1 + r2

    def prob(x):
        return choose(r1, x) * choose(r2, c1 - x) / choose(n, c1)

    p_obs = prob(a)
    total = 0.0
    for x in range(0, n + 1):  # 显式枚举所有表，不合法的自动为 0
        aa, bb = x, r1 - x
        cc, dd = c1 - x, r2 - c1 + x
        if min(aa, bb, cc, dd) < 0:
            continue
        p = prob(x)
        if p <= p_obs + 1e-15:
            total += p
    return total


class FisherIndependentCheck(unittest.TestCase):
    TABLES = [
        # (a, b, c, d, 手算精确 p)
        ((2, 0, 0, 1), Fraction(1, 3)),          # P={1/3,2/3}，取 <=1/3
        ((1, 1, 0, 1), Fraction(1)),             # 两张表都 <= 2/3
        ((1, 0, 0, 1), Fraction(1)),             # 各 1/2
        ((3, 1, 1, 3), Fraction(34, 70)),        # C(4,x)^2/70，取 <=16/70
        ((0, 3, 3, 0), Fraction(2, 20)),         # 两端表各 1/20
        ((5, 0, 0, 5), Fraction(2, 252)),        # 两端表各 1/252
        ((7, 3, 2, 8), Fraction(293, 4199)),     # 11720/167960 约分
        ((0, 0, 0, 0), Fraction(1)),             # 退化：唯一表
    ]

    def test_product_matches_independent_enumeration(self):
        for (a, b, c, d), expected in self.TABLES:
            got = fisher_exact_two_sided(a, b, c, d)
            self.assertEqual(got, expected, f"表 {(a, b, c, d)} 手算值不符")
            indep = independent_fisher_two_sided(a, b, c, d)
            self.assertAlmostEqual(float(got), indep, places=12,
                                   msg=f"表 {(a, b, c, d)} 与独立枚举不符")

    def test_odds_ratio_sentinels(self):
        self.assertEqual(odds_ratio(7, 3, 2, 8), Fraction(28, 3))
        self.assertEqual(odds_ratio(2, 0, 0, 1), "inf")   # bc=0, ad>0
        self.assertEqual(odds_ratio(2, 0, 1, 0), "undefined")  # ad=bc=0
        self.assertEqual(odds_ratio(0, 0, 0, 0), "undefined")
        self.assertEqual(odds_ratio(0, 1, 1, 0), Fraction(0))  # ad=0, bc>0
        self.assertEqual(odds_ratio(1, 1, 1, 1), Fraction(1))


# --------------------------------------------------------------------------- #
# BH 校正：含相同概率的多基因例子
# --------------------------------------------------------------------------- #
class BenjaminiHochbergCheck(unittest.TestCase):
    def test_tied_p_values_get_tied_adjustment(self):
        # 3 个 1/3 + 2 个 1：m=5
        # 1/3 组：min(5/1*1/3, 5/2*1/3, 5/3*1/3)=5/9；1 组：1
        got = benjamini_hochberg([Fraction(1, 3)] * 3 + [Fraction(1)] * 2)
        self.assertEqual(got, [Fraction(5, 9)] * 3 + [Fraction(1)] * 2)

    def test_monotone_example(self):
        # p = [1/2, 1/4, 1/10]，m=3（校正值按输入顺序返回）：
        # rank3: 1/2 -> 3/3*1/2 = 1/2；
        # rank2: 1/4 -> min(3/2*1/4, 1/2) = 3/8；
        # rank1: 1/10 -> min(3/1*1/10, 3/8, 1/2) = 3/10
        got = benjamini_hochberg(
            [Fraction(1, 2), Fraction(1, 4), Fraction(1, 10)])
        self.assertEqual(got, [Fraction(1, 2), Fraction(3, 8),
                               Fraction(3, 10)])

    def test_single_test_passthrough(self):
        self.assertEqual(benjamini_hochberg([Fraction(1, 20)]),
                         [Fraction(1, 20)])


# --------------------------------------------------------------------------- #
# 手写合成数据集：多基因相同概率的端到端校正/排序 + 大表 Fisher
# --------------------------------------------------------------------------- #
def write_synthetic_dataset(root, dataset_id, samples, rows):
    """按 vcfpipe 发布格式手写一个最小数据集目录。

    rows: (chrom, pos, ref, alt, sample, gt, dp)；dp 为 None 时写 '.'。
    """
    ds_dir = os.path.join(root, dataset_id)
    os.makedirs(os.path.join(ds_dir, "matrix"))
    with open(os.path.join(ds_dir, "manifest.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"dataset_id": dataset_id}, fh)
    keys = sorted({(r[0], r[1], r[2], r[3]) for r in rows})
    with open(os.path.join(ds_dir, "matrix", "sample_gt_matrix.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["assembly", "chrom", "pos", "ref", "alt"] + list(samples))
        cell = {(r[0], r[1], r[2], r[3], r[4]): r[5] for r in rows}
        for key in keys:
            w.writerow(["synT"] + list(key) +
                       [cell.get(key + (s,), "NA") for s in samples])
    with open(os.path.join(ds_dir, "matrix", "genotypes_long.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["assembly", "chrom", "pos", "ref", "alt", "sample",
                    "gt", "ad", "dp", "phased", "winner_caller"])
        for chrom, pos, ref, alt, sample, gt, dp in sorted(rows):
            w.writerow(["synT", chrom, pos, ref, alt, sample, gt, ".",
                        "." if dp is None else dp, "0", "hand"])
    return ds_dir


def write_file(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


class SyntheticDatasets(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_syn_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def p(self, name):
        return os.path.join(self.tmp, name)

    # ------------------------------------------------------------------ #
    def test_tied_p_multi_gene_correction_and_order(self):
        # 5 个被检验基因：GENE1/2/3 的表都是 (2,0,0,1)（p=1/3），
        # GENE4/5 都是 (1,1,0,1)（p=1）；GENE6 无候选位点；
        # GENE7 病例组全排除 -> 分母为零不检验。m=5。
        rows = []
        for pos in (10, 20, 30):
            rows += [("chr1", pos, "A", "G", "P1", "0/1", 30),
                     ("chr1", pos, "A", "G", "P2", "0/1", 30),
                     ("chr1", pos, "A", "G", "N1", "0/0", 30)]
        for pos in (40, 50):
            rows += [("chr1", pos, "A", "G", "P1", "0/1", 30),
                     ("chr1", pos, "A", "G", "P2", "0/0", 30),
                     ("chr1", pos, "A", "G", "N1", "0/0", 30)]
        rows += [("chr1", 60, "A", "G", "P1", "./.", None),
                 ("chr1", 60, "A", "G", "P2", "./.", None),
                 ("chr1", 60, "A", "G", "N1", "0/0", 30)]
        ds = write_synthetic_dataset(self.tmp, "tied", ["P1", "P2", "N1"],
                                     rows)
        labels = self.p("labels.tsv")
        write_file(labels, "sample\tgroup\nP1\tcase\nP2\tcase\nN1\tcontrol\n")
        genes = self.p("genes.bed")
        write_file(genes, "".join(
            f"chr1\t{pos - 1}\t{pos + 1}\tGENE{i}\n"
            for i, pos in enumerate((10, 20, 30, 40, 50, 60), start=1))
            + "chr1\t500\t510\tGENE7\n")
        # 注意：上面 GENE6 对应 pos60 区间，GENE7 无位点
        callable_bed = self.p("callable.bed")
        write_file(callable_bed, "")

        report = run_assoc(ds, labels, genes, callable_bed, 10)
        self.assertEqual(report["genes_tested"], 5)
        order = [e["gene"] for e in report["results"]]
        # 校正值 5/9 的三个按基因名稳定排序，然后才是校正值 1 的两个
        self.assertEqual(order, ["GENE1", "GENE2", "GENE3", "GENE4",
                                 "GENE5"])
        for e in report["results"][:3]:
            self.assertEqual(e["fisher_p_two_sided"], Fraction(1, 3))
            self.assertEqual(e["bh_adjusted_p"], Fraction(5, 9))
            self.assertEqual(e["odds_ratio"], "inf")
        for e in report["results"][3:]:
            self.assertEqual(e["fisher_p_two_sided"], Fraction(1))
            self.assertEqual(e["bh_adjusted_p"], Fraction(1))
        not_tested = {e["gene"]: e for e in report["not_tested"]}
        self.assertEqual(not_tested["GENE6"]["status"], "zero_denominator")
        self.assertIn("case", not_tested["GENE6"]["reason"])
        self.assertEqual(not_tested["GENE7"]["status"], "no_candidate_loci")

    def test_larger_table_fisher_and_or(self):
        # 20 人单基因：病例 7 携带/3 未携带，对照 2/8
        # p = 11720/167960 = 293/4199；OR = 56/6 = 28/3
        samples = [f"C{i}" for i in range(1, 11)] + \
                  [f"D{i}" for i in range(1, 11)]
        rows = []
        for i in range(1, 11):
            rows.append(("chr1", 10, "A", "G", f"C{i}",
                         "0/1" if i <= 7 else "0/0", 30))
        for i in range(1, 11):
            rows.append(("chr1", 10, "A", "G", f"D{i}",
                         "0/1" if i <= 2 else "0/0", 30))
        ds = write_synthetic_dataset(self.tmp, "big", samples, rows)
        labels = self.p("labels.tsv")
        write_file(labels, "sample\tgroup\n" + "".join(
            f"{s}\t{'case' if s.startswith('C') else 'control'}\n"
            for s in samples))
        genes = self.p("genes.bed")
        write_file(genes, "chr1\t9\t11\tGENE_X\n")
        callable_bed = self.p("callable.bed")
        write_file(callable_bed, "")

        report = run_assoc(ds, labels, genes, callable_bed, 10)
        (e,) = report["results"]
        self.assertEqual(e["table"], {"case_carrier": 7,
                                      "case_non_carrier": 3,
                                      "control_carrier": 2,
                                      "control_non_carrier": 8})
        self.assertEqual(e["fisher_p_two_sided"], Fraction(293, 4199))
        self.assertAlmostEqual(float(e["fisher_p_two_sided"]),
                               independent_fisher_two_sided(7, 3, 2, 8),
                               places=12)
        self.assertEqual(e["odds_ratio"], Fraction(28, 3))
        self.assertEqual(e["bh_adjusted_p"], Fraction(293, 4199))  # m=1


# --------------------------------------------------------------------------- #
# 输入校验：所有不合法输入都必须报清楚并拒绝
# --------------------------------------------------------------------------- #
class InputValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_val_")
        run_pipeline(cls.tmp, FA, DICT, SHEET, PRIO, "v1")
        cls.dataset = os.path.join(cls.tmp, "syn37_v1")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="vcfpipe_assoc_valcase_")

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    def w(self, name, text):
        path = os.path.join(self.work, name)
        write_file(path, text)
        return path

    def assert_assoc_error(self, needle, labels=None, genes=None,
                           callable_bed=None, min_dp=10):
        labels = labels or self.w("labels.tsv",
                                  "sample\tgroup\nS1\tcase\nS3\tcontrol\n")
        genes = genes or self.w("genes.bed", "chr1\t0\t12\tG1\n")
        callable_bed = callable_bed if callable_bed is not None \
            else self.w("callable.bed", "chr1\t0\t106\tS1\n")
        with self.assertRaises(AssocInputError) as ctx:
            run_assoc(self.dataset, labels, genes, callable_bed, min_dp)
        self.assertIn(needle, str(ctx.exception))

    def test_duplicate_sample_rejected(self):
        bad = self.w("bad.tsv",
                     "sample\tgroup\nS1\tcase\nS1\tcase\nS3\tcontrol\n")
        self.assert_assoc_error("重复", labels=bad)

    def test_group_conflict_rejected(self):
        bad = self.w("bad.tsv",
                     "sample\tgroup\nS1\tcase\nS1\tcontrol\nS3\tcontrol\n")
        self.assert_assoc_error("组别冲突", labels=bad)
        self.assert_assoc_error("S1", labels=bad)

    def test_unknown_sample_rejected(self):
        bad = self.w("bad.tsv",
                     "sample\tgroup\nS1\tcase\nS9\tcontrol\n")
        self.assert_assoc_error("不存在于数据集", labels=bad)
        self.assert_assoc_error("S9", labels=bad)

    def test_invalid_group_label_rejected(self):
        bad = self.w("bad.tsv",
                     "sample\tgroup\nS1\tCASE\nS3\tcontrol\n")
        self.assert_assoc_error("不合法的组别标签", labels=bad)

    def test_labels_missing_column_rejected(self):
        bad = self.w("bad.tsv", "sample\tphenotype\nS1\tcase\n")
        self.assert_assoc_error("缺列", labels=bad)

    def test_negative_coordinate_rejected(self):
        bad = self.w("bad.bed", "chr1\t-5\t12\tG1\n")
        self.assert_assoc_error("负坐标", genes=bad)

    def test_empty_interval_rejected(self):
        bad = self.w("bad.bed", "chr1\t12\t12\tG1\n")
        self.assert_assoc_error("空区间", genes=bad)
        bad2 = self.w("bad2.bed", "chr1\t12\t5\tG1\n")
        self.assert_assoc_error("空区间", genes=bad2)

    def test_non_integer_coordinate_rejected(self):
        bad = self.w("bad.bed", "chr1\ta\t12\tG1\n")
        self.assert_assoc_error("整数", genes=bad)

    def test_bed_wrong_column_count_rejected(self):
        bad = self.w("bad.bed", "chr1\t0\t12\n")
        self.assert_assoc_error("4 列", genes=bad)

    def test_empty_gene_name_rejected(self):
        bad = self.w("bad.bed", "chr1\t0\t12\t\n")
        self.assert_assoc_error("基因名为空", genes=bad)

    def test_callable_unknown_sample_rejected(self):
        bad = self.w("bad.bed", "chr1\t0\t106\tS9\n")
        self.assert_assoc_error("不存在于数据集", callable_bed=bad)

    def test_callable_negative_coordinate_rejected(self):
        bad = self.w("bad.bed", "chr1\t-1\t50\tS1\n")
        self.assert_assoc_error("负坐标", callable_bed=bad)

    def test_negative_min_dp_rejected(self):
        self.assert_assoc_error("非负整数", min_dp=-1)

    def test_missing_dataset_file_rejected(self):
        with self.assertRaises(AssocInputError) as ctx:
            run_assoc(os.path.join(self.work, "no_such_dir"),
                      self.w("labels.tsv",
                             "sample\tgroup\nS1\tcase\nS3\tcontrol\n"),
                      self.w("genes.bed", "chr1\t0\t12\tG1\n"),
                      self.w("callable.bed", "chr1\t0\t106\tS1\n"), 10)
        self.assertIn("数据集缺少已发布文件", str(ctx.exception))


# --------------------------------------------------------------------------- #
# 实际命令结果（子进程）：组别冲突、一组分母归零、输出文件、目录零改动
# --------------------------------------------------------------------------- #
def snapshot_tree(directory):
    snap = {}
    for dirpath, _dirs, files in os.walk(directory):
        for name in files:
            full = os.path.join(dirpath, name)
            with open(full, "rb") as fh:
                snap[os.path.relpath(full, directory)] = fh.read()
    return snap


class CommandLine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="vcfpipe_assoc_cli_")
        run_pipeline(cls.tmp, FA, DICT, SHEET, PRIO, "v1")
        cls.dataset = os.path.join(cls.tmp, "syn37_v1")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def cli(self, *args):
        return subprocess.run(
            [sys.executable, "-I", RUNNER, "assoc", *args],
            cwd=ROOT, capture_output=True, text=True, timeout=120)

    def demo_args(self, labels="labels_demo.tsv"):
        return ["--dataset", self.dataset,
                "--labels", demo_file(labels),
                "--genes", demo_file("genes_demo.bed"),
                "--callable", demo_file("callable_demo.bed"),
                "--min-dp", "10"]

    def test_group_conflict_exit_2_with_clear_message(self):
        r = self.cli(*self.demo_args("labels_conflict.tsv"))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("组别冲突", r.stderr)
        self.assertIn("S1", r.stderr)
        self.assertEqual(r.stdout, "")

    def test_zero_denominator_gene_not_tested(self):
        r = self.cli(*self.demo_args("labels_zero_denominator.tsv"))
        self.assertEqual(r.returncode, 0, r.stderr)
        report = json.loads(r.stdout)
        self.assertEqual(report["genes_tested"], 5)
        tested_genes = {e["gene"] for e in report["results"]}
        self.assertNotIn("GENE_D", tested_genes)
        not_tested = {e["gene"]: e for e in report["not_tested"]}
        self.assertEqual(not_tested["GENE_D"]["status"], "zero_denominator")
        self.assertIn("case", not_tested["GENE_D"]["reason"])
        self.assertEqual(
            not_tested["GENE_D"]["counts"]["case"]["denominator"], 0)
        self.assertEqual(not_tested["GENE_F"]["status"],
                         "no_candidate_loci")

    def test_out_dir_writes_three_reports(self):
        out_dir = os.path.join(self.tmp, "assoc_out")
        r = self.cli(*self.demo_args(), "--out-dir", out_dir)
        self.assertEqual(r.returncode, 0, r.stderr)
        for name in ("assoc_report.json", "assoc_summary.tsv",
                     "assoc_decisions.tsv"):
            self.assertTrue(os.path.isfile(os.path.join(out_dir, name)),
                            f"缺少 {name}")
        with open(os.path.join(out_dir, "assoc_report.json"),
                  encoding="utf-8") as fh:
            report = json.load(fh)
        self.assertEqual(report["genes_tested"], 6)
        with open(os.path.join(out_dir, "assoc_summary.tsv"),
                  encoding="utf-8", newline="") as fh:
            rows = {row["gene"]: row for row in
                    csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(rows["GENE_B"]["fisher_p_two_sided"],
                         repr(1 / 3))
        self.assertEqual(rows["GENE_B"]["odds_ratio"], "inf")
        self.assertEqual(rows["GENE_A"]["odds_ratio"], "undefined")
        self.assertEqual(rows["GENE_F"]["status"], "no_candidate_loci")
        with open(os.path.join(out_dir, "assoc_decisions.tsv"),
                  encoding="utf-8", newline="") as fh:
            dec = list(csv.DictReader(fh, delimiter="\t"))
        hit = [d for d in dec if d["gene"] == "GENE_D"
               and d["sample"] == "S1" and d["alt"] == "T"]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit[0]["locus_decision"], "unknown")
        self.assertEqual(hit[0]["reason"], "unknown_fully_missing")

    def test_assoc_leaves_published_dataset_untouched(self):
        before = snapshot_tree(self.dataset)
        r = self.cli(*self.demo_args())
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(snapshot_tree(self.dataset), before)

    @unittest.skipUnless(os.path.isdir(CHECKED_IN_DATASET),
                         "仓库未携带预生成 output/ 目录")
    def test_checked_in_dataset_untouched_by_assoc(self):
        before = snapshot_tree(CHECKED_IN_DATASET)
        r = subprocess.run(
            [sys.executable, "-I", RUNNER, "assoc",
             "--dataset", CHECKED_IN_DATASET,
             "--labels", demo_file("labels_demo.tsv"),
             "--genes", demo_file("genes_demo.bed"),
             "--callable", demo_file("callable_demo.bed"),
             "--min-dp", "10"],
            cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(snapshot_tree(CHECKED_IN_DATASET), before)

    def test_json_report_uses_sentinels_not_pseudocounts(self):
        r = self.cli(*self.demo_args())
        self.assertEqual(r.returncode, 0, r.stderr)
        report = json.loads(r.stdout)
        by_gene = {e["gene"]: e for e in report["results"]}
        self.assertEqual(by_gene["GENE_B"]["odds_ratio"], "inf")
        self.assertEqual(by_gene["GENE_A"]["odds_ratio"], "undefined")
        self.assertIn("合成数据", report["disclaimer"])
        self.assertIn("不构成临床诊断", report["disclaimer"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
