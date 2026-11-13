#!/usr/bin/env python3
"""dedup 独立慢速核对器（independent brute-force cross-check）。

本脚本**完全不导入** vcfpipe.dedup 的分词、候选检索与交集函数，只按
规格用标准库重新实现：

* 分词：NFKC 后 casefold；Unicode 类别以 L/N 开头的连续段为 token；
* 片段：连续三 token 的集合（Python tuple 集合，原值比较）；
* 清理：按 quality 降序、id 字典序，对每篇**扫描全部在先保留项**
  （无任何前缀/长度过滤），精确 Jaccard 用整数交叉乘法判定，同分
  按清理顺序取最早——与生产命令的磁盘前缀倒排 + 长度过滤互相对照，
  专门检查过滤层没有漏掉边界对。

自检模式（默认）对多个有理阈值分别合成含“恰好等于 / 刚低于 / 刚高
于”阈值边界对的语料，叠加非传递链、高频模板、全角/大小写、重复
token、短文本等情形，并对每个语料做多种打乱行序，逐张处置表（含
见证、交集/并集、相似度）、kept.jsonl 原始行与 summary 计数全面
比对::

    python3 -I tools/dedup_independent_check.py
    python3 -I tools/dedup_independent_check.py --shuffles 5

单点核对（检查一份已产出的结果）::

    python3 -I tools/dedup_independent_check.py \
        --corpus data/corpus/demo_small.jsonl \
        --out-dir output/dedup_demo --threshold 1/2

退出码 0 表示全部一致；任一不符即非零并打印差异。仅用于测试核对。
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
import unicodedata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNNER = os.path.join(ROOT, "run_pipeline.py")

THRESHOLDS = ["1/2", "1/3", "2/3", "3/5", "2/5", "1/1"]


# --------------------------------------------------------------------------- #
# 独立实现（从规格重写，不调用生产函数）
# --------------------------------------------------------------------------- #
def ind_tokens(text: str) -> list[str]:
    """独立分词：NFKC -> casefold -> L*/N* 连续段。"""
    norm = unicodedata.normalize("NFKC", text).casefold()
    tokens, current = [], []
    for ch in norm:
        cat = unicodedata.category(ch)
        if cat.startswith("L") or cat.startswith("N"):
            current.append(ch)
        elif current:
            tokens.append("".join(current))
            current = []
    if current:
        tokens.append("".join(current))
    return tokens


def ind_shingles(tokens: list[str]) -> set[tuple[str, str, str]]:
    """独立片段集合：连续三 token，重复只算一个。"""
    if len(tokens) < 3:
        return set()
    return {tuple(tokens[i:i + 3]) for i in range(len(tokens) - 2)}


def brute_force(records: list[dict], p: int, q: int):
    """慢速全量清理：每篇扫描所有在先保留项（无任何过滤）。

    返回 {id: 处置字典}。处置字典含片段集合，供逐字段核对。
    """
    order = sorted(records, key=lambda r: (-r["quality"], r["id"]))
    kept: list[tuple[str, set]] = []          # 仅正常保留项参与比较
    result: dict[str, dict] = {}
    for rec in order:
        rid = rec["id"]
        tokens = ind_tokens(rec["text"])
        shingles = ind_shingles(tokens)
        if not shingles:
            result[rid] = {"disposition": "too_short",
                           "token_count": len(tokens),
                           "shingles": shingles}
            continue
        best = None  # (inter, union, witness_id, witness_size)
        for k_id, k_set in kept:
            inter = len(shingles & k_set)
            if inter == 0:
                continue
            union = len(shingles) + len(k_set) - inter
            # 严格更大才替换 -> 同分时保留清理顺序最早者
            if best is None or inter * best[1] > best[0] * union:
                best = (inter, union, k_id, len(k_set))
        if best is not None and best[0] * q >= p * best[1]:
            result[rid] = {"disposition": "redundant",
                           "shingles": shingles,
                           "witness_id": best[2],
                           "witness_shingle_count": best[3],
                           "intersection": best[0],
                           "union": best[1]}
        else:
            kept.append((rid, shingles))
            result[rid] = {"disposition": "kept", "shingles": shingles}
    return result, [r["id"] for r in order]


# --------------------------------------------------------------------------- #
# 合成语料：边界对 + 各情形
# --------------------------------------------------------------------------- #
class Synth:
    """全局唯一 token 发生器：含私有 token 的片段在全库唯一。"""

    def __init__(self) -> None:
        self.c = 0

    def uniq(self, n: int) -> list[str]:
        toks = [f"u{self.c + i}" for i in range(n)]
        self.c += n
        return toks


def find_sizes(p: int, q: int, relation: str):
    """找最小 (n, m, s)：|A|=n、|B|=m、|A∩B|=s，使 s*q 与 p*u 满足关系。

    u = n+m-s；relation: eq 恰好等于 / below 刚低于(差1) / above 刚高于(差1)。
    """
    for s in range(1, 80):
        for n in range(s, 80):
            for m in range(s, 80):
                u = n + m - s
                lhs, rhs = s * q, p * u
                if relation == "eq" and lhs == rhs:
                    return n, m, s
                if relation == "below" and lhs == rhs - 1:
                    return n, m, s
                if relation == "above" and lhs == rhs + 1:
                    return n, m, s
    return None


def make_pair(syn: Synth, n: int, m: int, s: int):
    """两篇文档的 token 序列：片段集合大小 n、m，交集恰好 s。

    共享段 s+2 个 token -> s 个内部片段；其余 token 全局唯一，含私有
    token 的片段不会跨文档重合，故交集恰为共享段内部片段。
    """
    shared = syn.uniq(s + 2)
    a = shared + syn.uniq(n - s)
    b = shared + syn.uniq(m - s)
    return a, b


def build_corpus(p: int, q: int) -> list[dict]:
    """为阈值 p/q 合成一份覆盖边界与各情形的语料。"""
    syn = Synth()
    recs: list[dict] = []

    def add(tokens=None, quality=0, rid=None, text=None):
        recs.append({
            "id": rid or f"doc{len(recs):04d}",
            "quality": quality,
            "text": " ".join(tokens) if text is None else text,
        })

    # 1) 边界三对：恰好等于 / 刚低于 / 刚高于（专门打前缀+长度过滤的边界）
    for relation, tag in (("eq", "eq"), ("below", "below"),
                          ("above", "above")):
        sizes = find_sizes(p, q, relation)
        if sizes is None:  # 如 1/1 不存在“刚高于”
            continue
        n, m, s = sizes
        ta, tb = make_pair(syn, n, m, s)
        add(ta, 100, f"{tag}-a")
        add(tb, 90, f"{tag}-b")

    # 2) 非传递三篇链：J(a,b)=1/2、J(b,c)=1/2、J(a,c)=1/5
    chain = syn.uniq(7)
    add(chain[0:5], 30, "chain-a")
    add(chain[1:6], 20, "chain-b")
    add(chain[2:7], 10, "chain-c")

    # 3) 高频模板簇：公共 8 token（6 片段）+ 不同长度的私有尾
    tpl = syn.uniq(8)
    for i, tail in enumerate((0, 0, 2, 6, 10)):
        add(tpl + syn.uniq(tail), 50 - i, f"tpl-{i}")

    # 4) 全角/大小写/标点变体（独立分词同样要规范化）
    add(quality=40, rid="norm-1", text="Ｈｅｌｌｏ Ｗｏｒｌｄ Ｆｏｏ Ｂａｒ")
    add(quality=39, rid="norm-2", text="HELLO,, world!!foo...BAR")
    add(quality=38, rid="norm-3", text="Straße Ｘｙｚ １２３ 末端")
    add(quality=37, rid="norm-4", text="STRASSE XYZ 123 末端")

    # 5) 重复 token：片段集合去重
    add(quality=36, rid="rep-1", text="na na na na na batman")
    add(quality=35, rid="rep-2", text="na na na batman")

    # 6) 同质量同文：按 id 字典序决胜
    add(quality=34, rid="tie-b", text="same text here today")
    add(quality=34, rid="tie-a", text="same text here today")

    # 7) 短文本：空文两篇、单 token、两 token、纯标点
    add(quality=33, rid="short-empty-1", text="")
    add(quality=32, rid="short-empty-2", text="")
    add(quality=31, rid="short-one", text="solo")
    add(quality=30, rid="short-two", text="just two")
    add(quality=29, rid="short-punct", text="！！！。。。")

    # 8) 互不交叉的独立文档
    for i in range(5):
        add(syn.uniq(8), 10 - i, f"solo-{i}")

    return recs


# --------------------------------------------------------------------------- #
# 生产命令运行与比对
# --------------------------------------------------------------------------- #
def run_production(corpus_path: str, out_dir: str, threshold: str) -> None:
    proc = subprocess.run(
        [sys.executable, "-I", RUNNER, "dedup",
         "--corpus", corpus_path, "--threshold", threshold,
         "--out-dir", out_dir],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"生产命令失败（退出码 {proc.returncode}）:\n{proc.stderr}")


def compare(corpus_path: str, out_dir: str, p: int, q: int) -> list[str]:
    """比对一份生产输出与慢速全量结果；返回差异描述列表。"""
    failures: list[str] = []
    records = []
    raw_by_id = {}
    with open(corpus_path, "rb") as fh:
        for raw in fh:
            line = raw.decode("utf-8").rstrip("\n")
            obj = json.loads(line)
            records.append(obj)
            raw_by_id[obj["id"]] = line

    brute, order = brute_force(records, p, q)

    # --- 处置表 ---------------------------------------------------------
    with open(os.path.join(out_dir, "dispositions.jsonl"),
              encoding="utf-8") as fh:
        disps = [json.loads(line) for line in fh]
    if [d["id"] for d in disps] != order:
        failures.append("处置表顺序与清理顺序不符")
    if [d["seq"] for d in disps] != list(range(1, len(order) + 1)):
        failures.append("处置表 seq 不是 1..N 的清理顺序")
    for d in disps:
        rid = d["id"]
        b = brute[rid]
        if d["disposition"] != b["disposition"]:
            failures.append(f"{rid}: 处置 {d['disposition']} != "
                            f"独立 {b['disposition']}")
            continue
        if b["disposition"] == "too_short":
            if d["token_count"] != b["token_count"] or d["shingle_count"] != 0:
                failures.append(f"{rid}: too_short 计数不符")
        elif b["disposition"] == "redundant":
            want = {
                "witness_id": b["witness_id"],
                "witness_shingle_count": b["witness_shingle_count"],
                "intersection": b["intersection"],
                "union": b["union"],
                "similarity": f"{b['intersection']}/{b['union']}",
                "shingle_count": len(b["shingles"]),
            }
            for key, val in want.items():
                if d.get(key) != val:
                    failures.append(f"{rid}: {key}={d.get(key)} != 独立 {val}")
        else:  # kept
            if d["shingle_count"] != len(b["shingles"]):
                failures.append(f"{rid}: 片段数 {d['shingle_count']} != "
                                f"独立 {len(b['shingles'])}")
            best = d.get("best_below_threshold")
            if best is not None:
                # 所述在先保留项的相似度必须属实且确实低于阈值
                w = brute[best["id"]]
                inter = len(brute[rid]["shingles"] & w["shingles"])
                union = len(brute[rid]["shingles"]) + len(w["shingles"]) - inter
                if (best["intersection"], best["union"]) != (inter, union):
                    failures.append(f"{rid}: best_below_threshold 数值不属实")
                if not inter * q < p * union:
                    failures.append(f"{rid}: best_below_threshold 达到阈值")

    # --- kept.jsonl：保留清单 + 原始行逐字节 -----------------------------
    with open(os.path.join(out_dir, "kept.jsonl"), "rb") as fh:
        kept_lines = fh.read().decode("utf-8").splitlines()
    expected_kept = [rid for rid in order
                     if brute[rid]["disposition"] != "redundant"]
    got_kept = [json.loads(line)["id"] for line in kept_lines]
    if got_kept != expected_kept:
        failures.append(f"kept.jsonl 清单 {got_kept} != 独立 {expected_kept}")
    for line in kept_lines:
        rid = json.loads(line)["id"]
        if line != raw_by_id[rid]:
            failures.append(f"{rid}: kept.jsonl 行与原始行不一致")

    # --- summary 计数 -----------------------------------------------------
    with open(os.path.join(out_dir, "summary.json"), encoding="utf-8") as fh:
        summary = json.load(fh)
    counts = {"kept": 0, "too_short": 0, "redundant": 0}
    for b in brute.values():
        counts[b["disposition"]] += 1
    want = {
        "documents": len(records),
        "kept": counts["kept"],
        "too_short": counts["too_short"],
        "redundant": counts["redundant"],
        "kept_jsonl_lines": counts["kept"] + counts["too_short"],
        "pairwise_total": len(records) * (len(records) - 1) // 2,
        "threshold": f"{p}/{q}",
    }
    for key, val in want.items():
        if summary.get(key) != val:
            failures.append(f"summary.{key}: {summary.get(key)} != 独立 {val}")
    if summary["exact_comparisons"] > summary["candidate_pairs"]:
        failures.append("summary: 精确比较数超过候选对数")
    return failures


def parse_threshold(text: str) -> tuple[int, int]:
    a, b = text.split("/")
    return int(a), int(b)


def self_check(shuffles: int) -> int:
    failures: list[str] = []
    checked = 0
    with tempfile.TemporaryDirectory(prefix="dedup_indcheck_") as work:
        for threshold in THRESHOLDS:
            p, q = parse_threshold(threshold)
            records = build_corpus(p, q)
            for seed in range(shuffles):
                rng = random.Random(20261113 + seed)
                shuffled = records[:]
                rng.shuffle(shuffled)
                corpus_path = os.path.join(work, "corpus.jsonl")
                with open(corpus_path, "w", encoding="utf-8",
                          newline="") as fh:
                    for rec in shuffled:
                        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                out_dir = os.path.join(work, "out")
                run_production(corpus_path, out_dir, threshold)
                for msg in compare(corpus_path, out_dir, p, q):
                    failures.append(f"[阈值 {threshold} 打乱 {seed}] {msg}")
                checked += 1
                # 清理本轮输出，避免目录非空
                for name in os.listdir(out_dir):
                    os.unlink(os.path.join(out_dir, name))
                os.rmdir(out_dir)
    if failures:
        print("独立慢速核对发现差异：", file=sys.stderr)
        for msg in failures:
            print(" -", msg, file=sys.stderr)
        print(f"共核对 {checked} 张处置表，{len(failures)} 处不符",
              file=sys.stderr)
        return 1
    print(f"独立慢速核对通过：{len(THRESHOLDS)} 个阈值 × {shuffles} 种打乱"
          f"共 {checked} 张处置表，逐篇处置、见证、交并集、相似度、"
          f"保留清单与原始行、汇总计数全部一致。")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", help="单点核对：语料 JSONL")
    ap.add_argument("--out-dir", help="单点核对：生产输出目录")
    ap.add_argument("--threshold", help="单点核对：阈值 p/q")
    ap.add_argument("--shuffles", type=int, default=3,
                    help="自检模式每个阈值的打乱行序次数（默认 3）")
    args = ap.parse_args()

    if args.corpus or args.out_dir or args.threshold:
        if not (args.corpus and args.out_dir and args.threshold):
            print("单点核对需要同时给出 --corpus/--out-dir/--threshold",
                  file=sys.stderr)
            return 2
        p, q = parse_threshold(args.threshold)
        failures = compare(args.corpus, args.out_dir, p, q)
        if failures:
            print("单点核对发现差异：", file=sys.stderr)
            for msg in failures:
                print(" -", msg, file=sys.stderr)
            return 1
        print("单点核对通过：处置表、保留清单、原始行与汇总计数"
              "与独立慢速结果一致。")
        return 0

    return self_check(args.shuffles)


if __name__ == "__main__":
    raise SystemExit(main())
