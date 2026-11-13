"""多等位记录拆分时的 GT / AD / DP 处理。

拆分规则（与 bcftools norm -m- 的等位视角一致）：

对每条 ALT 等位 k（原始 ALT 下标从 1 开始）生成一条单 ALT 记录，
其二值化 GT 为：
    原始 0      -> 0（参考）
    原始 k      -> 1（本 ALT）
    其它非 0/k  -> '.'（对该 ALT 而言是“其它等位”，半缺失语义）
    原始 '.'    -> '.'（缺失，绝不能改写成 0）

AD 只取 [AD[0], AD[k]]；DP 保持原值。
"""
from __future__ import annotations

from typing import Optional

from .model import GenotypeCall, RawRecord
from .vcf_io import MISSING


def split_gt(
    raw_gt: Optional[str],
    gt_alleles: Optional[tuple[Optional[int], ...]],
    phased: bool,
    alt_index: int,
) -> GenotypeCall:
    """把多等位 GT 投影到单个 ALT 上。

    ``alt_index`` 为该 ALT 在原始 ALT 列表中的 1-based 下标。
    """
    if raw_gt is None:
        # 源记录根本没有 GT：保留“无调用”，不伪造参考纯合
        return GenotypeCall(gt_alleles=None, ad=None, dp=None, phased=phased)
    if gt_alleles is None:
        # 整体缺失 ./. 或 .|.
        return GenotypeCall(gt_alleles=None, ad=None, dp=None, phased=phased)
    if all(a is None for a in gt_alleles):
        # 所有倍型位均缺失（./.、.|.、.）：整体缺失，
        # 必须与半缺失（./1）区分，绝不能折叠成参考纯合
        return GenotypeCall(gt_alleles=None, ad=None, dp=None, phased=phased)

    remapped: list[Optional[int]] = []
    for allele in gt_alleles:
        if allele is None:
            remapped.append(None)
        elif allele == 0:
            remapped.append(0)
        elif allele == alt_index:
            remapped.append(1)
        else:
            # 属于其它 ALT：对本 ALT 不可判定，记为缺失而非 0
            remapped.append(None)
    return GenotypeCall(tuple(remapped), None, None, phased)


def split_record_calls(record: RawRecord) -> list[GenotypeCall]:
    """为一条原始记录的每个 ALT 生成拆分后的 GenotypeCall（含 AD/DP）。"""
    n_alts = len(record.alts)
    phased = record.phased

    # 整体缺失（无 GT、'.'、'./.'、'.|.'）：每个 ALT 都得到缺失调用，
    # 保留 DP 但不携带 AD，绝不伪造为参考纯合
    token = record.raw_gt
    fully_missing = (
        token is None or token == MISSING
        or all(p == MISSING for p in token.replace("|", "/").split("/"))
    )
    if fully_missing:
        return [
            GenotypeCall(gt_alleles=None, ad=None, dp=record.raw_dp,
                         phased=phased)
            for _ in range(n_alts)
        ]

    # 与 vcf_io.iter_records 保持一致地还原 GT 结构
    sep = "|" if phased else "/"
    gt_alleles = tuple(
        None if part == MISSING else int(part)
        for part in token.split(sep)
    )

    calls: list[GenotypeCall] = []
    for idx in range(1, n_alts + 1):
        call = split_gt(record.raw_gt, gt_alleles, phased, idx)
        ad: Optional[tuple[Optional[int], Optional[int]]] = None
        if record.raw_ad is not None and idx < len(record.raw_ad):
            ad = (record.raw_ad[0], record.raw_ad[idx])
        calls.append(
            GenotypeCall(
                gt_alleles=call.gt_alleles,
                ad=ad,
                dp=record.raw_dp,
                phased=phased,
            )
        )
    return calls


def dosage_conservation(
    raw_gt: Optional[str], split_calls: list[GenotypeCall]
) -> bool:
    """守恒检查：源 GT 中每个非参考、非缺失等位，恰好对应一个拆分 ALT 剂量。

    用于测试断言；半缺失 GT 的缺失部分不参与计数。
    """
    if raw_gt is None or raw_gt == MISSING:
        return all(c.alt_dosage == 0 for c in split_calls)
    sep = "|" if "|" in raw_gt else "/"
    source_alt_copies = [
        int(p) for p in raw_gt.split(sep) if p not in (MISSING, "0")
    ]
    source_total = len(source_alt_copies)
    split_total = sum(c.alt_dosage for c in split_calls)
    return source_total == split_total
