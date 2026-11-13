#!/usr/bin/env python3
"""独立的小范围穷举核对器（independent brute-force cross-check）。

本脚本**完全不导入** vcfpipe.trio 的补全/传递/判定函数，只用 Python
标准库重新实现一遍枚举逻辑，对含缺失的多等位三人组逐一穷举并核对：

1. 从 GT 文本独立解析（``.``=完全未知、``.`` 等位位=该位缺失、``|``
   与 ``/`` 等价），深度不足的整调用按完全未知处理（原 GT 与原因另存）；
2. 把每个人缺失的等位位用该位点全部合法编号 0..n_alts 补全，收集
   **不同的无序二倍体基因型**（``0/1`` 与 ``1/0`` 视为同一个）；
3. 笛卡尔积枚举父×母×孩子的全部组合，对每组直接检查孩子两个等位能否
   分别在父亲、母亲的基因型中各取一个（多重集合比较），独立计数；
4. 按 0 可解释 / 全部可解释 / 部分可解释归 incompatible / compatible /
   uncertain，并与生产命令写出的 trio_details.tsv 逐行核对总数、可解释
   数、候选基因型集合与判定类别。

用法::

    python3 -I tools/trio_exhaustive_check.py \
        --details output/trio_screen/trio_details.tsv \
        --vcf data/trio/trio_demo.vcf \
        --min-dp 10

退出码 0 表示逐行全部一致；任一不符即非零并打印差异。仅用于测试核对。
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from itertools import product


def open_text(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def parse_gt_independent(token: str):
    """独立 GT 解析，返回 (等位列表或 None 表示完全未知, 倍型长度)。"""
    if token == ".":
        return None, 2
    sep = "|" if "|" in token else "/"
    parts = token.split(sep)
    alleles = [None if p == "." else int(p) for p in parts]
    if all(a is None for a in alleles):
        return None, len(alleles)
    return alleles, len(alleles)


def complete_unordered(alleles, n_alts: int):
    """独立补全：缺失位枚举全部合法编号，无序去重，排序后返回。"""
    if alleles is None:
        alleles = [None, None]
    pools = [tuple(range(n_alts + 1)) if a is None else (a,)
             for a in alleles]
    seen = set()
    for combo in product(*pools):
        seen.add(tuple(sorted(combo)))
    return sorted(seen)


def can_make_child(father_gt, mother_gt, child_gt):
    """独立传递检查：孩子多重集合能否由父一个 + 母一个等位组成。"""
    target = tuple(sorted(child_gt))
    for paternal in father_gt:
        for maternal in mother_gt:
            if tuple(sorted((paternal, maternal))) == target:
                return True
    return False


def brute_force(raw_gt: str, dp, min_dp: int, n_alts: int):
    """返回统一四元组 (verdict 标记, 候选列表, 总组合, 可解释组合)。

    非二倍体可信调用标记为 ``"unsupported"``；其余情形标记为 ``"ok"``，
    总/可解释组合留给 :func:`evaluate` 枚举，这里只产出单人候选集合。
    """
    alleles, ploidy = parse_gt_independent(raw_gt)
    low_depth = dp is None or dp < min_dp
    if not low_depth and alleles is not None and ploidy != 2:
        return "unsupported", [], 0, 0
    if low_depth:
        # 低深度/缺 DP：整调用未知，已知等位不再约束补全
        candidates = complete_unordered(None, n_alts)
    else:
        candidates = complete_unordered(alleles, n_alts)
    return "ok", candidates, 0, 0


def evaluate(father_gt, mother_gt, child_gt,
             father_dp, mother_dp, child_dp, min_dp, n_alts):
    f = brute_force(father_gt, father_dp, min_dp, n_alts)
    m = brute_force(mother_gt, mother_dp, min_dp, n_alts)
    c = brute_force(child_gt, child_dp, min_dp, n_alts)
    if f[0] == "unsupported" or m[0] == "unsupported" or c[0] == "unsupported":
        return "unsupported", (f[1], m[1], c[1]), 0, 0
    fc, mc, cc = f[1], m[1], c[1]
    total = 0
    explainable = 0
    for fg, mg, cg in product(fc, mc, cc):
        total += 1
        if can_make_child(fg, mg, cg):
            explainable += 1
    if explainable == 0:
        verdict = "incompatible"
    elif explainable == total:
        verdict = "compatible"
    else:
        verdict = "uncertain"
    return verdict, (fc, mc, cc), total, explainable


def read_vcf_loci(path: str):
    """独立读取 VCF 选中位点：line_no -> (n_alts, {sample: (gt, dp)})。"""
    loci = {}
    with open_text(path) as fh:
        line_no = 0
        samples = []
        for raw in fh:
            line_no += 1
            line = raw.rstrip("\n\r")
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            n_alts = 0 if cols[4] == "." else len(cols[4].split(","))
            fmt = cols[8].split(":")
            calls = {}
            for name, blob in zip(samples, cols[9:]):
                vals = dict(zip(fmt, blob.split(":")))
                dp_raw = vals.get("DP", ".")
                dp = None if dp_raw in (".", "") else int(dp_raw)
                calls[name] = (vals["GT"], dp)
            loci[line_no] = (n_alts, calls)
    return loci


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--details", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--min-dp", required=True, type=int)
    args = ap.parse_args()

    if not os.path.isfile(args.details):
        print(f"明细文件不存在: {args.details}", file=sys.stderr)
        return 2
    if not os.path.isfile(args.vcf):
        print(f"VCF 不存在: {args.vcf}", file=sys.stderr)
        return 2

    loci = read_vcf_loci(args.vcf)
    failures = []
    checked = 0
    with open(args.details, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            checked += 1
            line_no = int(row["line_no"])
            n_alts, calls = loci[line_no]
            f_gt, f_dp = calls[row["father"]]
            m_gt, m_dp = calls[row["mother"]]
            c_gt, c_dp = calls[row["child"]]
            verdict, cand, total, explainable = evaluate(
                f_gt, m_gt, c_gt, f_dp, m_dp, c_dp,
                args.min_dp, n_alts)

            if verdict != row["verdict"]:
                failures.append(
                    f"行 {line_no} {row['family']}/{row['child']}: 判定 "
                    f"{row['verdict']} != 独立 {verdict}")
                continue
            if verdict == "unsupported":
                continue
            if total != int(row["total_combinations"]) or \
                    explainable != int(row["explainable_combinations"]):
                failures.append(
                    f"行 {line_no} {row['family']}/{row['child']}: 计数 "
                    f"({row['total_combinations']},"
                    f"{row['explainable_combinations']}) != 独立 "
                    f"({total},{explainable})")
            for role, got, col in (
                ("father", cand[0], "father_candidates"),
                ("mother", cand[1], "mother_candidates"),
                ("child", cand[2], "child_candidates"),
            ):
                got_text = ",".join("/".join(map(str, g)) for g in got)
                if got_text != row[col]:
                    failures.append(
                        f"行 {line_no} {row['family']}/{row['child']} "
                        f"{role} 候选 {row[col]!r} != 独立 {got_text!r}")

    if failures:
        print("独立穷举核对发现差异：", file=sys.stderr)
        for msg in failures:
            print(" -", msg, file=sys.stderr)
        print(f"共核对 {checked} 行，{len(failures)} 处不符", file=sys.stderr)
        return 1
    print(f"独立穷举核对通过：{checked} 行明细的候选集合、总组合数、"
          f"可解释组合数与判定类别全部一致。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
