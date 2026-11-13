"""规范化与基因型拆分的手算单元测试。

所有期望值均依据合成参考 syn37 的真实碱基手工推导，
不引用管线自身的输出，避免“自证正确”。
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "src"))

from vcfpipe.fasta import Reference
from vcfpipe.genotype import dosage_conservation, split_record_calls, split_gt
from vcfpipe.model import RawRecord, SourceVCF
from vcfpipe.normalize import RuleBundle, normalize_variant
from vcfpipe.vcf_io import MISSING

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FA = os.path.join(ROOT, "data", "ref", "syn37.fa")

V1 = RuleBundle.for_version("v1")
V2 = RuleBundle.for_version("v2")


def src(caller="deepvar", assembly="syn37", sample="S1"):
    return SourceVCF(path="x.vcf", assembly=assembly, sample=sample,
                     caller=caller, digest="d", summary="s")


def raw(chrom="chr1", pos=1, ref="A", alts=("G",), gt="0/1",
        ad=None, dp=None, phased=False, line_no=1):
    return RawRecord(source=src(), line_no=line_no, chrom=chrom, pos=pos,
                     ref=ref, alts=tuple(alts), qual=".", vfilter=".",
                     info=".", raw_gt=gt, raw_ad=ad, raw_dp=dp, phased=phased)


class NormalizeHandCalc(unittest.TestCase):
    def setUp(self):
        self.ref = Reference(FA)

    def norm(self, pos, r, a, rule=V1):
        n, reason = normalize_variant(self.ref, "chr1", pos, r, a, rule)
        return n, reason

    def test_snp_unchanged(self):
        # chr1:5 参考 A；SNP 无公共前后缀
        n, reason = self.norm(5, "A", "G")
        self.assertIsNone(reason)
        self.assertEqual((n.pos, n.ref, n.alt), (5, "A", "G"))
        self.assertEqual(n.moved, 0)

    def test_polyA_insertion_left_align_to_anchor(self):
        # pos37 A>AA：连续 A 区段 29..38。逐次同步左扩弹后缀，
        # 最终锚定到区段起点的前一个碱基 -> 28 G>GA
        n, _ = self.norm(37, "A", "AA")
        self.assertEqual((n.pos, n.ref, n.alt), (28, "G", "GA"))
        self.assertEqual(n.moved, 9)

    def test_polyA_deletion_matches_insertion_anchor(self):
        # pos37 AA>A 的删除与上面的插入共享锚点 28：28 GA>G
        n, _ = self.norm(37, "AA", "A")
        self.assertEqual((n.pos, n.ref, n.alt), (28, "GA", "G"))

    def test_repeat_T_insertion_left_alignment(self):
        # pos57 T>TT：T 区段 49..58，前一位 48=C -> 48 C>CT
        n, _ = self.norm(57, "T", "TT")
        self.assertEqual((n.pos, n.ref, n.alt), (48, "C", "CT"))

    def test_dinuc_repeat_deletion(self):
        # pos46 CG>C：末碱基 G!=C，首碱基 C==C 但 ALT 长度为 1 不可裁，
        # 因此已是最小/最左表示，保持 46 CG>C（管线事实表同值）
        n, _ = self.norm(46, "CG", "C")
        self.assertEqual((n.pos, n.ref, n.alt), (46, "CG", "C"))

    def test_deletion_with_shared_suffix_left_aligns(self):
        # pos46 CGC>CC（末 C==C）：左扩 pos45=G 得 GCGC/GCC，弹 C
        # -> GCG/GC @45；末 G!=C 停止；首 G 公共但 ALT=GC 可裁前缀
        # -> 46 CG>C。验证多步“左扩弹后缀 + 前缀裁剪”的组合行为。
        n, _ = self.norm(46, "CGC", "CC")
        self.assertEqual((n.pos, n.ref, n.alt), (46, "CG", "C"))

    def test_different_representations_converge(self):
        # 同一插入的两种书写必须规范化到同一键
        n1, _ = self.norm(37, "A", "AA")
        n2, _ = self.norm(28, "G", "GA")
        self.assertEqual((n1.pos, n1.ref, n1.alt),
                         (n2.pos, n2.ref, n2.alt))

    def test_boundary_v1_quarantine_v2_keep(self):
        # pos1 A>CA：末碱基 A 相同但无法再取 pos-1
        n1, reason = self.norm(1, "A", "CA", V1)
        self.assertIsNone(n1)
        self.assertEqual(reason, "left_align_boundary")
        n2, reason2 = self.norm(1, "A", "CA", V2)
        self.assertIsNone(reason2)
        self.assertTrue(n2.boundary_flag)
        self.assertEqual((n2.pos, n2.ref, n2.alt), (1, "A", "CA"))

    def test_minimal_representation_trims_prefix(self):
        # CGC>CTG：末碱基 C==G 不等；首碱基 C==C 且两侧均 >1，
        # 裁前缀 -> pos45 GC>TG（最小表示，首碱基 G!=T）
        n, _ = self.norm(44, "CGC", "CTG")
        self.assertEqual((n.pos, n.ref, n.alt), (45, "GC", "TG"))

    def test_ref_equals_alt_rejected(self):
        n, reason = self.norm(5, "A", "A")
        self.assertIsNone(n)
        self.assertEqual(reason, "malformed_record")


class GenotypeSplitHandCalc(unittest.TestCase):
    def test_tri_allelic_phased_dosage(self):
        # 1|2 拆到 ALT#1 -> 1|.；拆到 ALT#2 -> .|1
        rec = raw(alts=("A", "T"), gt="1|2", ad=(30, 5, 4), dp=39,
                  phased=True)
        calls = split_record_calls(rec)
        self.assertEqual(calls[0].gt_alleles, (1, None))
        self.assertEqual(calls[1].gt_alleles, (None, 1))
        self.assertTrue(all(c.phased for c in calls))
        self.assertEqual(calls[0].ad, (30, 5))
        self.assertEqual(calls[1].ad, (30, 4))
        self.assertTrue(dosage_conservation("1|2", calls))

    def test_fully_missing_not_ref_hom(self):
        rec = raw(gt="./.", ad=None, dp=None)
        calls = split_record_calls(rec)
        self.assertIsNone(calls[0].gt_alleles)  # 绝不是 (0, 0)
        self.assertTrue(calls[0].is_missing)
        self.assertEqual(calls[0].alt_dosage, 0)
        self.assertTrue(dosage_conservation("./.", calls))

    def test_half_missing_keeps_known_dose(self):
        # ./1：一个真实 ALT 剂量 + 一个缺失
        rec = raw(gt="./1", ad=(6, 5), dp=11)
        calls = split_record_calls(rec)
        self.assertEqual(calls[0].gt_alleles, (None, 1))
        self.assertEqual(calls[0].alt_dosage, 1)
        self.assertEqual(calls[0].ref_dosage, 0)
        self.assertTrue(dosage_conservation("./1", calls))

    def test_other_alt_becomes_missing_not_zero(self):
        # 三等位，GT=0/2：投影到 ALT#1 时等位 2 必须变缺失而非 0
        rec = raw(alts=("T", "C"), gt="0/2", ad=(20, 3, 7), dp=30)
        calls = split_record_calls(rec)
        self.assertEqual(calls[0].gt_alleles, (0, None))
        self.assertEqual(calls[1].gt_alleles, (0, 1))
        # 剂量守恒：源里 1 个非参等位 -> 拆分后总和 1
        self.assertEqual(sum(c.alt_dosage for c in calls), 1)
        self.assertTrue(dosage_conservation("0/2", calls))

    def test_ref_hom_dosage_zero(self):
        rec = raw(gt="0/0", ad=(20, 0), dp=20)
        calls = split_record_calls(rec)
        self.assertEqual(calls[0].gt_alleles, (0, 0))
        self.assertEqual(calls[0].alt_dosage, 0)
        self.assertEqual(calls[0].ref_dosage, 2)

    def test_dp_preserved_even_when_gt_missing(self):
        rec = raw(gt="./.", ad=None, dp=7)
        calls = split_record_calls(rec)
        self.assertIsNone(calls[0].gt_alleles)
        self.assertEqual(calls[0].dp, 7)

    def test_phased_separator_roundtrip(self):
        rec = raw(gt="1|0", ad=(8, 9), dp=17, phased=True)
        calls = split_record_calls(rec)
        self.assertEqual(calls[0].gt_alleles, (1, 0))
        self.assertTrue(calls[0].phased)

    def test_no_gt_field_stays_no_call(self):
        rec = RawRecord(source=src(), line_no=1, chrom="chr1", pos=5,
                        ref="A", alts=("G",), qual=".", vfilter=".",
                        info=".", raw_gt=None, raw_ad=None, raw_dp=None,
                        phased=False)
        calls = split_record_calls(rec)
        self.assertIsNone(calls[0].gt_alleles)


if __name__ == "__main__":
    unittest.main(verbosity=2)
