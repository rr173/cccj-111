#!/usr/bin/env python3
"""dedup 手写小语料演练：先列三词片段与预期处置，再运行命令读回证据。

语料：``data/corpus/demo_small.jsonl``（21 篇，阈值 1/2）。脚本分三段：

1. 手工推导：逐篇打印 token 序列与三词片段集合，给出预期处置表
   （保留 / too_short / redundant→见证 与相似度），全部写死在脚本里，
   不来自运行结果；
2. 实际运行：调用 ``vcfpipe.dedup.run_dedup`` 处理同一份语料，逐项
   核对处置表、保留清单与原始行（含附带字段）是否逐字节读回；
3. 拒绝路径：重复 id 与非法阈值必须报错并拒绝作业，不产出貌似成功
   的结果（不创建输出目录）。

用法::

    python3 -I tools/run_dedup_demo.py
    python3 -I tools/run_dedup_demo.py --out-dir output/dedup_demo

退出码 0 表示实际结果与手工预期完全一致。
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from vcfpipe.dedup import (  # noqa: E402
    DISPOSITION_KEPT, DISPOSITION_REDUNDANT, DISPOSITION_TOO_SHORT,
    DedupInputError, run_dedup, shingle_triples, tokenize,
)

CORPUS = os.path.join(ROOT, "data", "corpus", "demo_small.jsonl")
DEFAULT_OUT = os.path.join(ROOT, "output", "dedup_demo")
THRESHOLD = "1/2"

# --------------------------------------------------------------------------- #
# 手工推导的预期处置（阈值 1/2；先于此运行结果存在）
# ---------------------------------------------------------------------------
# (id, quality, 预期处置, 见证, 相似度)
EXPECTED = [
    # 非传递链：a..e / b..f / c..g；中间归第一篇(2/4=1/2，恰好等于也算)，
    # 第三篇与第一篇只有 1/5，必须留下（不能沿链误删）
    ("chain-a", 40, DISPOSITION_KEPT, None, None),
    ("chain-b", 30, DISPOSITION_REDUNDANT, "chain-a", "2/4"),
    # 一个冗余同时命中两个保留项：multi-r 对两者都是 2/4，同分取最早
    ("multi-k1", 29, DISPOSITION_KEPT, None, None),
    ("multi-k2", 28, DISPOSITION_KEPT, None, None),
    # 全角 + 大小写：NFKC+casefold 后与 fw-1 完全相同；同质量按 id 决胜
    ("fw-1", 27, DISPOSITION_KEPT, None, None),
    ("fw-2", 27, DISPOSITION_REDUNDANT, "fw-1", "2/2"),
    ("multi-r", 26, DISPOSITION_REDUNDANT, "multi-k1", "2/4"),
    # 刚低于阈值：交集 1、并集 3 -> 1/3 < 1/2，两篇都保留
    ("edge-a", 25, DISPOSITION_KEPT, None, None),
    ("edge-b", 24, DISPOSITION_KEPT, None, None),
    # 标点（含全角标点）全作分隔，连续分隔不产生空 token
    ("punct-1", 23, DISPOSITION_KEPT, None, None),
    ("punct-2", 22, DISPOSITION_REDUNDANT, "punct-1", "2/2"),
    # 链的第三篇：与 chain-a 仅 1/5，chain-b 已冗余不参与 -> 保留
    ("chain-c", 21, DISPOSITION_KEPT, None, None),
    # 重复片段只算一个：na×5+batman 与 na×3+batman 的片段集合相同
    ("dup-1", 20, DISPOSITION_KEPT, None, None),
    ("dup-2", 19, DISPOSITION_REDUNDANT, "dup-1", "2/2"),
    # 短文本：保留原文但不入索引；两篇空文不因此互相删除
    ("short-empty-1", 18, DISPOSITION_TOO_SHORT, None, None),
    ("short-empty-2", 17, DISPOSITION_TOO_SHORT, None, None),
    ("short-two", 16, DISPOSITION_TOO_SHORT, None, None),
    # 高频模板：8 个公共 token -> 6 个公共片段，6/8 >= 1/2
    ("tpl-1", 15, DISPOSITION_KEPT, None, None),
    ("tpl-2", 14, DISPOSITION_REDUNDANT, "tpl-1", "6/8"),
    # 互不交叉的词表：无候选
    ("uniq-1", 13, DISPOSITION_KEPT, None, None),
    ("uniq-2", 12, DISPOSITION_KEPT, None, None),
]

EXPECTED_SUMMARY = {
    "documents": 21, "kept": 12, "too_short": 3, "redundant": 6,
}


def read_corpus(path):
    records = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            records.append((line.rstrip("\n"), json.loads(line)))
    return records


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=DEFAULT_OUT,
                    help="输出目录（须不存在或为空）")
    args = ap.parse_args()

    records = read_corpus(CORPUS)
    expected_by_id = {e[0]: e for e in EXPECTED}

    print("=" * 72)
    print(f"第一段：手工推导（语料 {os.path.relpath(CORPUS, ROOT)}，"
          f"阈值 {THRESHOLD}）")
    print("=" * 72)
    for _raw, obj in records:
        tokens = tokenize(obj["text"])
        triples = shingle_triples(tokens)
        exp = expected_by_id[obj["id"]]
        print(f"\n[{obj['id']}] quality={obj['quality']}")
        print(f"  原文: {obj['text']!r}")
        print(f"  token({len(tokens)}): {' '.join(tokens) if tokens else '(无)'}")
        if triples:
            for t in sorted(triples):
                print(f"    片段: {' / '.join(t)}")
        else:
            print("    片段: (不足三个 token，无片段)")
        if exp[2] == DISPOSITION_REDUNDANT:
            print(f"  预期: redundant -> {exp[3]} (相似度 {exp[4]})")
        else:
            print(f"  预期: {exp[2]}")

    print("\n" + "=" * 72)
    print("第二段：实际运行并读回证据")
    print("=" * 72)
    summary = run_dedup(CORPUS, THRESHOLD, args.out_dir)
    out = args.out_dir

    failures = []

    # --- 处置表逐篇核对 -------------------------------------------------
    with open(os.path.join(out, "dispositions.jsonl"),
              encoding="utf-8") as fh:
        disps = [json.loads(line) for line in fh]
    order = [d["id"] for d in disps]
    expected_order = [e[0] for e in sorted(
        EXPECTED, key=lambda e: (-e[1], e[0]))]
    if order != expected_order:
        failures.append(f"清理顺序不符: {order} != {expected_order}")
    for d in disps:
        exp = expected_by_id[d["id"]]
        got = (d["disposition"], d.get("witness_id"), d.get("similarity"))
        want = (exp[2], exp[3], exp[4])
        if got != want:
            failures.append(f"{d['id']}: 处置 {got} != 预期 {want}")
    print(f"处置表: {len(disps)} 篇逐篇核对"
          f"（{'一致' if not failures else '有差异'}）")
    for d in disps:
        line = json.dumps(d, ensure_ascii=False)
        print(f"  {line}")

    # --- 保留清单与原始行逐字节读回 --------------------------------------
    with open(os.path.join(out, "kept.jsonl"), "rb") as fh:
        kept_lines = fh.read().decode("utf-8").splitlines()
    raw_by_id = {obj["id"]: raw for raw, obj in records}
    expected_kept_ids = [e[0] for e in sorted(
        (e for e in EXPECTED if e[2] != DISPOSITION_REDUNDANT),
        key=lambda e: (-e[1], e[0]))]
    got_kept_ids = [json.loads(line)["id"] for line in kept_lines]
    if got_kept_ids != expected_kept_ids:
        failures.append(f"保留清单不符: {got_kept_ids} != {expected_kept_ids}")
    for line in kept_lines:
        rid = json.loads(line)["id"]
        if line != raw_by_id[rid]:
            failures.append(f"{rid}: kept.jsonl 行与原始行不一致")
    print(f"\nkept.jsonl: {len(kept_lines)} 行（正常保留 + too_short），"
          f"原始行含附带字段逐字节读回"
          f"（{'一致' if not failures else '有差异'}）")
    for line in kept_lines:
        print(f"  {line}")

    # --- 汇总计数 ---------------------------------------------------------
    with open(os.path.join(out, "summary.json"), encoding="utf-8") as fh:
        summary_on_disk = json.load(fh)
    for key, want in EXPECTED_SUMMARY.items():
        if summary_on_disk[key] != want:
            failures.append(f"summary.{key}: {summary_on_disk[key]} != {want}")
    print(f"\nsummary.json: {json.dumps(summary_on_disk, ensure_ascii=False, sort_keys=True)}")

    # --- 第三段：拒绝路径（重复 id / 非法阈值不产出貌似成功） ------------
    print("\n" + "=" * 72)
    print("第三段：拒绝路径（不能跳过坏输入当作清理成功）")
    print("=" * 72)
    work = tempfile.mkdtemp(prefix="dedup_demo_reject_")
    try:
        dup_corpus = os.path.join(work, "dup.jsonl")
        with open(dup_corpus, "w", encoding="utf-8") as fh:
            fh.write('{"id":"a","quality":1,"text":"x y z"}\n')
            fh.write('{"id":"b","quality":1,"text":"p q r"}\n')
            fh.write('{"id":"a","quality":2,"text":"m n o"}\n')
        dup_out = os.path.join(work, "dup_out")
        try:
            run_dedup(dup_corpus, THRESHOLD, dup_out)
            failures.append("重复 id 未被拒绝")
            print("重复 id: 未拒绝（错误！）")
        except DedupInputError as exc:
            print(f"重复 id -> 拒绝（指出位置）: {exc}")
            if os.path.exists(dup_out):
                failures.append("重复 id 却创建了输出目录")
        bad_out = os.path.join(work, "bad_out")
        try:
            run_dedup(CORPUS, "0.5", bad_out)
            failures.append("非法阈值未被拒绝")
            print("非法阈值: 未拒绝（错误！）")
        except DedupInputError as exc:
            print(f"非法阈值 0.5 -> 拒绝: {exc}")
            if os.path.exists(bad_out):
                failures.append("非法阈值却创建了输出目录")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if failures:
        print("\n演练发现差异：", file=sys.stderr)
        for msg in failures:
            print(" -", msg, file=sys.stderr)
        return 1
    print("\n演练通过：手工预期与实际处置、保留清单、原始行、汇总计数"
          "全部一致；重复 id 与非法阈值均被拒绝且未产出结果。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
