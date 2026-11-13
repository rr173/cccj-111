#!/usr/bin/env python3
"""生成注释演示的全部合成材料（确定性、可重跑）。

产物
----
* ``data/ref/annot_syn.fa``：两条合成 contig。
    - ``chrP`` 按“正链密码子 + 填充内含子”设计；
    - ``chrM`` 由同一设计序列的反向互补片段拼成，专用于负链异构体。
* ``data/annot/transcripts_demo.tsv``：三个异构体
    - ``TXP_POS1`` 正链、三个外显子、含跨接缝密码子与内部终止密码子；
    - ``TXP_NEG1`` 负链、两个外显子、含跨接缝密码子与内部终止密码子；
    - ``TXP_POS1_UTR`` 正链、同外显子但 CDS 更短（``[118,250)``），
      chrP:95 在它里面是 5' 非编码片段；
   - ``TXP_NONCODING`` 正链、与 TXP_POS1 相同外显子、
      ``cds_start == cds_end``（非编码异构体）。
* ``data/annot/variants_demo.tsv``：覆盖编码/非编码/边界/内含子/基因间
  以及多 ALT、双碱基、跨两个密码子、终止获得/丢失等演示点。
* ``data/annot/regression_ref_mismatch.tsv``、
  ``data/annot/regression_bad_transcripts.tsv``：失败路径回归输入。

输入表的 REF 由脚本从组装好的 contig 直接截取，避免手抄坐标漂移；
ALT 显式给出并校验“异于 REF、彼此互异”。脚本只写输入，不做注释
计算；全部期望值由 ``tools/annot_independent_check.py`` 独立推导。
"""
from __future__ import annotations

import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
REF_DIR = os.path.join(ROOT, "data", "ref")

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


# --------------------------------------------------------------------------- #
# 设计序列
# --------------------------------------------------------------------------- #
# 正链 CDS：50 个密码子。第 1 位 ATG，第 15 位 CGA（R），第 16 位 TAA（*），
# 第 17 位 TGG（W），第 50 位 TCG（S）。第 7、37 密码子横跨外显子接缝。
PLUS_CODONS = (
    "ATG GCT TTC AAA GGA CCC TAC GAT TGG CTA "
    "AGC GGT CTC AGA CGA TAA TGG GCT GAC TTT "
    "GGC ATC CTG TCC AAC GTT CCT GAA CGT ACT "
    "GCA TAT TTG ACG TGC CAT GTA TCT AGG CCA "
    "GAG TGG ATA GGG TTC AAC GTT CCT GAA TCG"
).split()
# 负链 CDS：14 个密码子（已是翻译方向）。第 8 位 TGC（C）横跨接缝，
# 第 13 位 TGA（*）。
MINUS_CODONS = "ATG GCT TCG TGG TAC GAC GGT TGC GAA TTC AAA GGT TGA GCT".split()

assert len(PLUS_CODONS) == 50 and all(len(c) == 3 for c in PLUS_CODONS)
assert len(MINUS_CODONS) == 14 and all(len(c) == 3 for c in MINUS_CODONS)

PLUS_CDS = "".join(PLUS_CODONS)
MINUS_CDS = "".join(MINUS_CODONS)

# chrP 外显子（零起点半开）：[0,30) [90,180) [210,300)。
# 编码范围 [10,250)：三段编码贡献 20 + 90 + 40 = 150 碱基（整除 3）。
P_PREFIX = "ACGTACGTAC"                     # [0,10)
P_EXON1_CODING = PLUS_CDS[:20]              # [10,30)
P_INTRON1 = "TGCATGCA" * 7 + "TGCA"         # [30,90)  60 碱基
P_EXON2_CODING = PLUS_CDS[20:110]           # [90,180) 90 碱基
P_INTRON2 = "ACGT" * 7 + "AC"               # [180,210) 30 碱基
P_EXON3_CODING = PLUS_CDS[110:150]          # [210,250) 40 碱基
P_EXON3_TAIL = "ACGT" * 12 + "AC"           # [250,300) 50 碱基 3' 非编码
P_TRAILING = "ACGT" * 5                     # [300,320) contig 尾部

# chrM 外显子（零起点半开）：[0,40) [60,104)。
# 负链编码片段 [20,40) 与 [60,82)；按坐标递减做反向互补拼接 = MINUS_CDS
# （22 + 20 = 42 碱基，整除 3）。第 8 密码子跨 [60,82)/[20,40) 接缝。
M_PREFIX = "ACGTACGTACGTACGTACGT"           # [0,20)
M_G1 = reverse_complement(MINUS_CDS[22:42])  # [20,40)
M_INTRON = "TGCATGCATGCATGCATGCA"           # [40,60)
M_G2 = reverse_complement(MINUS_CDS[0:22])  # [60,82)
M_TAIL = "ACGT" * 5 + "AC"                  # [82,104) 22 碱基非编码
M_TRAILING = "ACGT" * 24                    # [104,200)


def build_plus() -> str:
    return "".join([
        P_PREFIX, P_EXON1_CODING, P_INTRON1,
        P_EXON2_CODING, P_INTRON2,
        P_EXON3_CODING, P_EXON3_TAIL, P_TRAILING,
    ])


def build_minus() -> str:
    return "".join([M_PREFIX, M_G1, M_INTRON, M_G2, M_TAIL, M_TRAILING])


PLUS = build_plus()
MINUS = build_minus()

# 设计断言：外显子边界、编码贡献、负链拼接方向必须与转录本表一致。
assert len(PLUS) == 320, len(PLUS)
assert len(MINUS) == 200, len(MINUS)
assert PLUS[10:30] == PLUS_CDS[:20]
assert PLUS[90:180] == PLUS_CDS[20:110]
assert PLUS[210:250] == PLUS_CDS[110:150]
assert reverse_complement(MINUS[60:82]) + reverse_complement(MINUS[20:40]) == MINUS_CDS

TRANSCRIPTS = """transcript\tgene\tchrom\tstrand\texon_starts\texon_ends\tcds_start\tcds_end
TXP_POS1\tGENEPOS\tchrP\t+\t0,90,210\t30,180,300\t10\t250
TXP_POS1_UTR\tGENEPOS\tchrP\t+\t0,90,210\t30,180,300\t118\t250
TXP_NEG1\tGENENEG\tchrM\t-\t0,60\t40,104\t20\t82
TXP_NONCODING\tGENENC\tchrP\t+\t0,90,210\t30,180,300\t0\t0
"""

# --------------------------------------------------------------------------- #
# 变异点
# --------------------------------------------------------------------------- #
CONFIGS = {"chrP": PLUS, "chrM": MINUS}

# (id, chrom, pos, ref_len, alts, 备注)；REF 从 CONFIGS 实截，不手抄。
VARIANT_SPECS = [
    # 正链跨接缝密码子 7（TAC，基因组 29/30/91），同一点三个 ALT
    ("V_CROSS_POS", "chrP", 30, 1, "G,T,C", "接缝密码子第 2 碱基"),
    # 负链跨接缝密码子 8（TGC，基因组 61/40/39）
    ("V_NEG_SEAM", "chrM", 40, 1, "T", "负链接缝密码子第 2 碱基"),
    # 终止丢失：密码子 16 TAA -> CAA（T>C）
    ("V_STOP_LOST", "chrP", 116, 1, "C", "TAA->CAA"),
    # 终止获得：密码子 17 TGG -> TGA（G>A）
    ("V_STOP_GAINED", "chrP", 121, 1, "A", "TGG->TGA"),
    # 同义：密码子 2 GCT -> GCC（第 3 碱基 T>C）
    ("V_SYNONYMOUS", "chrP", 16, 1, "C", "GCT->GCC"),
    # 双碱基同时改：密码子 15 CGA -> AGC（单改各自同义，合改错义）
    ("V_DOUBLE_BASE", "chrP", 113, 3, "AGC", "CGA->AGC"),
    # 一段替换跨密码子 14/15
    ("V_TWO_CODONS", "chrP", 112, 3, "TTA", "AGA|CGA -> ?"),
    # 同一位置在编码异构体与全非编码异构体中的不同区域类别
    ("V_CODING_VS_NONCODING", "chrP", 95, 1, "A", "TGG->AGG"),
    # 非编码片段：5' 与 3'
    ("V_UTR5", "chrP", 5, 1, "G", "5' 非编码"),
    ("V_UTR3", "chrP", 260, 1, "T", "3' 非编码"),
    # 单一内含子
    ("V_INTRON", "chrP", 50, 1, "C", "内含子"),
    # CDS 两端点：首碱基（密码子 1）与末碱基（密码子 50）
    ("V_CDS_START", "chrP", 11, 1, "G", "ATG->GTG"),
    ("V_CDS_END", "chrP", 250, 1, "A", "TCG->TCA"),
    # 外显子末端、跨接缝密码子 37 的最后一个外显子 2 碱基
    ("V_EXON_END_CODING", "chrP", 180, 1, "C", "GTA->GCA"),
    # REF 从外显子末端伸入内含子 -> boundary_unsupported
    ("V_BOUNDARY_EXON", "chrP", 30, 2, "CG", "外显子/内含子边界"),
    # REF 跨非编码/编码交界 -> boundary_unsupported
    ("V_BOUNDARY_CODING", "chrP", 10, 2, "GT", "非编码/编码边界"),
    # REF 从最后一个外显子碱基伸到转录本跨度之外 -> boundary_unsupported
    ("V_BOUNDARY_SPAN", "chrP", 300, 2, "AA", "跨出转录本跨度"),
    # 不与任何异构体相交 -> intergenic
    ("V_INTERGENIC", "chrP", 310, 1, "A", "基因间"),
    # 长度不等 -> unsupported_edit（仍保留相交 transcript）
    ("V_UNSUPPORTED", "chrP", 95, 1, "TA", "长度不等"),
    # 同一记录多 ALT：一个编码后果、一个长度不等
    ("V_ALT_MIXED", "chrP", 96, 1, "A,GA", "混合 ALT"),
    # 负链内含子与非编码片段
    ("V_NEG_INTRON", "chrM", 50, 1, "A", "负链内含子"),
    ("V_NEG_UTR", "chrM", 90, 1, "G", "负链 3' 非编码"),
]


def _check_variant(chrom: str, pos: int, ref_len: int, alts: str) -> str:
    seq = CONFIGS[chrom]
    if not 1 <= pos <= len(seq) or pos - 1 + ref_len > len(seq):
        raise AssertionError(f"{chrom}:{pos} 越界（len={len(seq)}）")
    ref = seq[pos - 1:pos - 1 + ref_len]
    if set(ref) - set("ACGT"):
        raise AssertionError(f"{chrom}:{pos} REF {ref!r} 非 ACGT")
    choices = alts.split(",")
    if len(set(choices)) != len(choices):
        raise AssertionError(f"{chrom}:{pos} ALT 重复: {alts}")
    if ref in choices:
        raise AssertionError(f"{chrom}:{pos} ALT 不得等于 REF {ref}")
    if any(set(a) - set("ACGT") or not a for a in choices):
        raise AssertionError(f"{chrom}:{pos} ALT 含非法碱基: {alts}")
    return ref


def main() -> int:
    os.makedirs(REF_DIR, exist_ok=True)

    def wrap(seq: str, width: int = 60) -> str:
        return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

    fasta = (">chrP 合成正链设计\n" + wrap(PLUS) + "\n\n"
             + ">chrM 合成负链设计\n" + wrap(MINUS) + "\n")
    path = os.path.join(REF_DIR, "annot_syn.fa")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(fasta)
    print(f"wrote {path}")

    with open(os.path.join(HERE, "transcripts_demo.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write(TRANSCRIPTS)

    rows = []
    for vid, chrom, pos, ref_len, alts, _note in VARIANT_SPECS:
        ref = _check_variant(chrom, pos, ref_len, alts)
        rows.append("\t".join([vid, chrom, str(pos), ref, alts]))
    with open(os.path.join(HERE, "variants_demo.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("id\tchrom\tpos\tref\talt\n")
        fh.write("\n".join(rows) + "\n")

    with open(os.path.join(HERE, "regression_ref_mismatch.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("id\tchrom\tpos\tref\talt\nV_REF_BAD\tchrP\t30\tC\tG\n")

    # 外显子重叠（[0,30) 与 [20,60)）且 CDS 端点不在外显子上 -> 必须拒绝
    with open(os.path.join(HERE, "regression_bad_transcripts.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write(
            "transcript\tgene\tchrom\tstrand\texon_starts\texon_ends\t"
            "cds_start\tcds_end\n"
            "TXP_BAD\tGENEBAD\tchrP\t+\t0,20\t30,60\t10\t40\n")
    print(f"wrote variants ({len(rows)} 行) 与回归输入")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
