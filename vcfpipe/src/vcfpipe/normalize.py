"""变异规范化：最左对齐（left alignment）与最小表示（minimal representation）。

采用与 Tan et al. 2015 (vt) / bcftools norm 等价的“先同步左扩、再弹公共
后缀、最后裁公共前缀”算法，正确处理 REF 长度为 1 的锚定插入：

    while 任一等位长度 > 1:
      if 两端最后碱基相同:
          # 长度为 1 的等位无法再弹后缀：先把 pos-1 的参考碱基同步补到
          # REF 与 ALT 左端（pos 左移 1），再各弹一个末碱基
          if pos == 1: 边界策略（见 RuleBundle）
          ref = base(pos-1) + ref ; alt = base(pos-1) + alt ; pos -= 1
          ref, alt = ref[:-1], alt[:-1]
      elif 两端首碱基相同且两者长度都 > 1:
          裁掉各自首碱基 ; pos += 1
      else:
          break

规则版本差异（由 RuleBundle 控制）：
  * v1：左移需要 pos-1 参考碱基但已到 contig 起点（pos == 1）时，
    记录无法完成规范化 -> 隔离（left_align_boundary）。
  * v2：同样情况下保留“已部分左移”的最小表示并加 BOUNDARY 标记，
    不隔离，以减少边界重复序列造成的丢失。
"""
from __future__ import annotations

from dataclasses import dataclass

from .fasta import Reference


@dataclass(frozen=True)
class RuleBundle:
    rule_version: str
    boundary_policy: str  # "quarantine" (v1) | "keep_flagged" (v2)

    @classmethod
    def for_version(cls, version: str) -> "RuleBundle":
        table = {
            "v1": cls("v1", "quarantine"),
            "v2": cls("v2", "keep_flagged"),
        }
        if version not in table:
            raise ValueError(f"未知规范化规则版本: {version}")
        return table[version]


@dataclass(frozen=True)
class Normalized:
    pos: int
    ref: str
    alt: str
    boundary_flag: bool
    moved: int  # 相对原始位置的左移量（调试/追踪用）


def normalize_variant(
    ref: Reference, chrom: str, pos: int, ref_allele: str, alt_allele: str,
    rule: RuleBundle,
) -> tuple[Normalized | None, str | None]:
    """返回 (Normalized, None) 或 (None, 原因码)。"""
    ref_allele = ref_allele.upper()
    alt_allele = alt_allele.upper()
    start_pos = pos
    boundary_hit = False

    if ref_allele == alt_allele or not ref_allele or not alt_allele:
        return None, "malformed_record"

    r, a = ref_allele, alt_allele

    while len(r) > 1 or len(a) > 1:
        if r[-1] == a[-1]:
            if pos <= 1:
                boundary_hit = True
                if rule.boundary_policy == "quarantine":
                    return None, "left_align_boundary"
                break  # keep_flagged：停止左移，随后做前缀裁剪
            prev_base = ref.base(chrom, pos - 1)
            # 标准 vt 左对齐：两等位同步左扩一个锚定碱基，再无条件
            # 各弹一个公共后缀。同步左扩保证 REF/ALT 等长且 REF 永不空。
            r = prev_base + r
            a = prev_base + a
            r, a = r[:-1], a[:-1]
            pos -= 1
        elif len(r) > 1 and len(a) > 1 and r[0] == a[0]:
            # 只有两侧都能保留至少一个碱基时才可裁公共前缀
            r, a = r[1:], a[1:]
            pos += 1
        else:
            break

    # 主循环的前缀分支已带“两侧长度均 >1”守卫，终止时不存在可再裁的
    # 公共前缀；此处只做防御性校验，不再进行任何裁剪。
    if r == a or not r or not a:
        return None, "malformed_record"

    return (
        Normalized(
            pos=pos, ref=r, alt=a, boundary_flag=boundary_hit,
            moved=start_pos - pos,
        ),
        None,
    )
