#!/usr/bin/env python3
"""dedup 规模演示：三万篇合成语料上的候选过滤与精确比较计数。

生成（固定种子，确定性；阈值 1/2）：

* 28,000 篇稀疏文档：私有词表互不交叉 -> 应几乎不产生候选；
* 300 个近重复簇 x 5：基准 25 token，成员替换 2 个 token
  （与基准 17/29 >= 1/2）-> 每簇 1 保留 + 4 冗余；
* 200 篇模板同文（同一文本的大小写/全角/标点变体）-> 1 保留 + 199 冗余；
* 300 篇“高频模板但未达阈值”：共享 18 个模板片段 + 16 个私有片段
  （两两 18/50 < 1/2）-> 全部保留；模板片段进入每篇前缀，两两互为
  候选，制造 44,850 对真实候选——展示“候选多时可慢但不漏”；
* 200 篇滑动窗口：相邻两篇 22/34 >= 1/2、隔篇 16/40 < 1/2
  -> 100 保留 + 100 冗余；
* 100 篇短文本 -> too_short。

随后用真实命令（run_pipeline.py dedup，子进程）处理，报告墙钟、
子进程峰值内存（ru_maxrss）、候选对数与精确比较数，并对照全体配
对规模 N(N-1)/2，证明稀疏部分确实经过索引过滤。本脚本只证明这份
输入上的行为；小语料的完整正确性由 tests/test_dedup.py 与
tools/dedup_independent_check.py 负责。

用法::

    python3 benchmarks/run_dedup_benchmark.py
    python3 benchmarks/run_dedup_benchmark.py --sparse 40000 --keep
"""
from __future__ import annotations

import argparse
import json
import os
import random
import resource
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNNER = os.path.join(ROOT, "run_pipeline.py")

SEED = 20261113
THRESHOLD = "1/2"

# 各家族规模
N_CLUSTERS = 300          # 近重复簇数（每簇 1 基准 + 4 成员）
CLUSTER_MEMBERS = 4
N_TPL_SAME = 200          # 模板同文组
N_TPL_FAMILY = 300        # 高频模板未达阈值家族
N_SLIDING = 200           # 滑动窗口
N_TOO_SHORT = 100         # 短文本

TPL_SAME_TOKENS = ("notice user 用户 order 订单 shipped 发货 tracking "
                   "单号 contact 联系 support 客服 thanks").split()
SEPARATORS = [" ", "  ", "，", "、", "。", " , ", "；"]


def render_tpl_same(rng: random.Random) -> str:
    """同一文本的大小写/全角/标点变体（规范化后 token 序列相同）。"""
    parts = []
    for tok in TPL_SAME_TOKENS:
        if tok.isascii():
            roll = rng.random()
            if roll < 0.3:
                tok = tok.upper()
            elif roll < 0.5:
                tok = tok.capitalize()
            if rng.random() < 0.4:
                tok = "".join(chr(ord(c) + 0xFEE0) for c in tok)  # 全角
        parts.append(tok)
    out = []
    for i, tok in enumerate(parts):
        out.append(tok)
        if i + 1 < len(parts):
            out.append(rng.choice(SEPARATORS))
    return "".join(out)


def generate(work: str, n_sparse: int) -> str:
    """流式生成合成语料；返回语料路径。"""
    rng = random.Random(SEED)
    path = os.path.join(work, "corpus.jsonl")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        def emit(rid, quality, text):
            fh.write(json.dumps({"id": rid, "quality": quality,
                                 "text": text}, ensure_ascii=False) + "\n")

        # 1) 稀疏文档：私有词表，互不交叉
        for i in range(n_sparse):
            toks = [f"w{i:05d}x{j:02d}" for j in range(28)]
            emit(f"sparse-{i:05d}", 100 + (i % 50), " ".join(toks))

        # 2) 近重复簇：基准 25 token，成员替换 2 个 token
        for c in range(N_CLUSTERS):
            base = [f"c{c:04d}w{j:02d}" for j in range(25)]
            emit(f"clu{c:04d}-base", 2000, " ".join(base))
            for m in range(1, CLUSTER_MEMBERS + 1):
                toks = base[:]
                for pos in (2 * m + 1, 2 * m + 11):
                    toks[pos] = f"c{c:04d}m{m}x{pos:02d}"
                emit(f"clu{c:04d}-m{m}", 1900 - m, " ".join(toks))

        # 3) 模板同文（大小写/全角/标点变体）
        for i in range(N_TPL_SAME):
            emit(f"tplsame-{i:04d}", 3000 - i, render_tpl_same(rng))

        # 4) 高频模板家族：共享 18 个模板片段 + 16 个私有片段
        tpl = [f"tpl{j:02d}" for j in range(20)]  # 18 个内部片段
        for i in range(N_TPL_FAMILY):
            left = [f"tf{i:04d}a{j:02d}" for j in range(8)]
            right = [f"tf{i:04d}b{j:02d}" for j in range(8)]
            emit(f"tplfam-{i:04d}", 500, " ".join(left + tpl + right))

        # 5) 滑动窗口：步长 6、窗宽 30
        for k in range(N_SLIDING):
            toks = [f"g{j:04d}" for j in range(6 * k, 6 * k + 30)]
            emit(f"slide-{k:04d}", 4000 - k, " ".join(toks))

        # 6) 短文本
        shorts = ["", "只", "两个 token", "！！！", "x"]
        for i in range(N_TOO_SHORT):
            emit(f"short-{i:04d}", 1, shorts[i % len(shorts)])
    return path


def verify(out_dir: str, n_sparse: int) -> list[str]:
    """按生成器的先验预期核对处置表与汇总；返回差异列表。"""
    failures = []
    with open(os.path.join(out_dir, "dispositions.jsonl"),
              encoding="utf-8") as fh:
        disps = [json.loads(line) for line in fh]
    expect = {"documents": n_sparse + N_CLUSTERS * 5 + N_TPL_SAME
              + N_TPL_FAMILY + N_SLIDING + N_TOO_SHORT,
              "kept": n_sparse + N_CLUSTERS + 1 + N_TPL_FAMILY
              + N_SLIDING // 2,
              "too_short": N_TOO_SHORT,
              "redundant": N_CLUSTERS * CLUSTER_MEMBERS + N_TPL_SAME - 1
              + N_SLIDING // 2}
    counts = {"kept": 0, "too_short": 0, "redundant": 0}
    for d in disps:
        counts[d["disposition"]] += 1
        rid = d["id"]
        if rid.startswith("clu") and "-m" in rid:
            # 簇成员必须冗余且见证为簇基准
            c = rid.split("-")[0]
            if d["disposition"] != "redundant" or \
                    d["witness_id"] != f"{c}-base":
                failures.append(f"{rid}: 应冗余于 {c}-base，实际 {d}")
        elif rid.startswith("tplsame-"):
            if rid == "tplsame-0000":
                if d["disposition"] != "kept":
                    failures.append("tplsame-0000 应保留")
            elif d["disposition"] != "redundant" or \
                    d["witness_id"] != "tplsame-0000":
                failures.append(f"{rid}: 应冗余于 tplsame-0000，实际 {d}")
        elif rid.startswith("tplfam-"):
            if d["disposition"] != "kept":
                failures.append(f"{rid}: 模板家族成员应保留，实际 {d}")
        elif rid.startswith("slide-"):
            k = int(rid.split("-")[1])
            if k % 2 == 0:
                if d["disposition"] != "kept":
                    failures.append(f"{rid}: 滑动窗口偶数篇应保留")
            elif d["disposition"] != "redundant" or \
                    d["witness_id"] != f"slide-{k - 1:04d}":
                failures.append(f"{rid}: 应冗余于 slide-{k - 1:04d}")
        elif rid.startswith("short-"):
            if d["disposition"] != "too_short":
                failures.append(f"{rid}: 应为 too_short")
        elif rid.startswith("sparse-") or rid.endswith("-base"):
            if d["disposition"] != "kept":
                failures.append(f"{rid}: 应保留")
    if len(disps) != expect["documents"]:
        failures.append(f"处置表篇数 {len(disps)} != {expect['documents']}")
    for key in ("kept", "too_short", "redundant"):
        if counts[key] != expect[key]:
            failures.append(f"{key}: {counts[key]} != 预期 {expect[key]}")

    with open(os.path.join(out_dir, "summary.json"), encoding="utf-8") as fh:
        summary = json.load(fh)
    for key in ("documents", "kept", "too_short", "redundant"):
        if summary[key] != expect[key]:
            failures.append(f"summary.{key}: {summary[key]} != {expect[key]}")
    # 模板家族 + 簇成员 + 模板同文贡献的候选对数有解析下界
    guaranteed = (N_TPL_FAMILY * (N_TPL_FAMILY - 1) // 2
                  + N_CLUSTERS * CLUSTER_MEMBERS + N_TPL_SAME - 1)
    if summary["candidate_pairs"] < guaranteed:
        failures.append(f"candidate_pairs {summary['candidate_pairs']} "
                        f"< 解析下界 {guaranteed}")
    if summary["exact_comparisons"] > summary["candidate_pairs"]:
        failures.append("exact_comparisons 超过 candidate_pairs")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sparse", type=int, default=28000,
                    help="稀疏文档篇数（默认 28000，总篇数 30300 起）")
    ap.add_argument("--keep", action="store_true",
                    help="保留临时目录（默认运行完删除）")
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix="dedup_bench_")
    try:
        t0 = time.perf_counter()
        corpus = generate(work, args.sparse)
        gen_s = time.perf_counter() - t0
        in_mb = os.path.getsize(corpus) / 1e6
        n_docs = (args.sparse + N_CLUSTERS * 5 + N_TPL_SAME
                  + N_TPL_FAMILY + N_SLIDING + N_TOO_SHORT)

        out_dir = os.path.join(work, "out")
        t0 = time.perf_counter()
        proc = subprocess.run(
            [sys.executable, "-I", RUNNER, "dedup",
             "--corpus", corpus, "--threshold", THRESHOLD,
             "--out-dir", out_dir])
        dedup_s = time.perf_counter() - t0
        if proc.returncode != 0:
            print("生产命令失败", file=sys.stderr)
            return 1
        child_rss_mb = resource.getrusage(
            resource.RUSAGE_CHILDREN).ru_maxrss / 1024

        failures = verify(out_dir, args.sparse)
        with open(os.path.join(out_dir, "summary.json"),
                  encoding="utf-8") as fh:
            summary = json.load(fh)

        report = {
            "documents": n_docs,
            "input_mb": round(in_mb, 1),
            "generate_seconds": round(gen_s, 2),
            "dedup_seconds": round(dedup_s, 2),
            "docs_per_second": int(n_docs / dedup_s),
            "child_max_rss_mb": round(child_rss_mb, 1),
            "kept": summary["kept"],
            "too_short": summary["too_short"],
            "redundant": summary["redundant"],
            "distinct_shingles": summary["distinct_shingles"],
            "candidate_pairs": summary["candidate_pairs"],
            "exact_comparisons": summary["exact_comparisons"],
            "pairwise_total": summary["pairwise_total"],
            "candidate_vs_all_pairs": (
                f"{summary['candidate_pairs'] / summary['pairwise_total']:.6%}"),
            "expectation_check": "通过" if not failures else "失败",
            "work_dir": work if args.keep else "(已清理)",
        }
        print(json.dumps(report, indent=2, ensure_ascii=False))
        if failures:
            print("先验预期核对发现差异：", file=sys.stderr)
            for msg in failures[:20]:
                print(" -", msg, file=sys.stderr)
            return 1
        return 0
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print(f"临时目录保留于 {work}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
