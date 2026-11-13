"""管线内使用的不可变数据记录。

所有记录均为 frozen dataclass，保证合并与重跑过程中不会被意外修改，
也便于做确定性序列化（排序 + 固定字段顺序）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# --------------------------------------------------------------------------- #
# 输入侧
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SourceVCF:
    """输入 VCF 的头摘要（来自显式 META 头行，不做文件名猜测）。"""

    path: str
    assembly: str
    sample: str
    caller: str
    digest: str
    summary: str
    contig_lengths: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class RawRecord:
    """VCF 中解析出的一条（可能是多等位的）原始记录。"""

    source: SourceVCF
    line_no: int
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    qual: str
    vfilter: str
    info: str
    raw_gt: Optional[str]
    raw_ad: Optional[tuple[Optional[int], ...]]
    raw_dp: Optional[int]
    phased: bool


@dataclass(frozen=True)
class GenotypeCall:
    """单等位拆分后保留的样本证据。

    注意语义：
      * ``gt_alleles`` 为 None 表示整条基因型缺失（./. 或 .|.），
        不得转成 (0, 0) 参考纯合。
      * 半缺失（如 ./1）保留 None 与真实等位混合，剂量只计真实等位。
      * ``phased`` 保留相位，单倍体/缺失同样记录。
    """

    gt_alleles: Optional[tuple[Optional[int], ...]]
    ad: Optional[tuple[Optional[int], Optional[int]]]
    dp: Optional[int]
    phased: bool

    @property
    def is_missing(self) -> bool:
        return self.gt_alleles is None

    @property
    def alt_dosage(self) -> int:
        """该 ALT 上的真实等位剂量；缺失等位不计。"""
        if self.gt_alleles is None:
            return 0
        return sum(1 for a in self.gt_alleles if a == 1)

    @property
    def ref_dosage(self) -> int:
        if self.gt_alleles is None:
            return 0
        return sum(1 for a in self.gt_alleles if a == 0)


@dataclass(frozen=True)
class SplitRecord:
    """多等位拆分 + 规范化后的单 ALT 变异（尚未合并）。"""

    source: SourceVCF
    line_no: int
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: str
    vfilter: str
    info: str
    call: GenotypeCall

    def key(self) -> tuple[str, str, int, str, str]:
        return (self.assembly, self.chrom, self.pos, self.ref, self.alt)


@dataclass(frozen=True)
class QuarantineRecord:
    """无法纳入事实表的记录及机器可读原因。"""

    source_path: str
    digest: str
    line_no: int
    assembly: str
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    raw_gt: Optional[str]
    reason_code: str
    reason_detail: str


# --------------------------------------------------------------------------- #
# 输出侧
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Evidence:
    """同一规范化位点上一个 caller 来源的字段证据（冲突时全部保留）。"""

    source_path: str
    digest: str
    line_no: int
    caller: str
    qual: str
    vfilter: str
    info: str
    raw_gt: Optional[str]
    raw_ad: Optional[tuple[Optional[int], ...]]
    raw_dp: Optional[int]
    phased: bool


@dataclass(frozen=True)
class MergedVariant:
    """事实表中的一行：稳定键 + 胜出样本基因型 + 全部来源证据。"""

    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str
    sample: str
    winner_caller: str
    call: GenotypeCall
    evidence: tuple[Evidence, ...]

    def key(self) -> tuple[str, str, int, str, str]:
        return (self.assembly, self.chrom, self.pos, self.ref, self.alt)


# 隔离原因码（稳定标识符，规则版本间不随意改名）
REASON_REF_MISMATCH = "ref_base_mismatch"
REASON_CONTIG_UNMAPPED = "contig_unmapped"
REASON_SYMBOLIC_ALT = "symbolic_alt"
REASON_BAD_GT = "malformed_gt"
REASON_BOUNDARY = "left_align_boundary"
REASON_BAD_RECORD = "malformed_record"
