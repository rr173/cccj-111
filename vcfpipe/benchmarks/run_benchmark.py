"""demux 流式处理的规模/内存演示：生成合成双端 FASTQ 并计时、测内存。

用法（默认 12 万对、48 个样本、K=8）：

    python3 benchmarks/run_benchmark.py
    python3 benchmarks/run_benchmark.py --pairs 500000 --samples 96

输入由固定种子的生成器合成（仅合成序列），流式写出 gzip；随后直接调用
vcfpipe.demux.run_demux 处理，报告墙钟、吞吐与进程峰值内存
（resource.ru_maxrss）。流式实现的特征是：内存占用不随对数增长。
"""
import argparse
import gzip
import json
import os
import random
import resource
import shutil
import sys
import tempfile
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

from vcfpipe.demux import run_demux  # noqa: E402

SEED = 20261113
BC_LEN = 8
UMI_LEN = 8
INSERT_LEN = 92   # R1 插入段
R2_LEN = 100


def generate(pairs: int, samples: int, work: str):
    """流式生成合成输入；返回 (barcodes_tsv, r1_path, r2_path)。"""
    rng = random.Random(SEED)
    barcodes = {}
    while len(barcodes) < samples:
        bc = "".join(rng.choices("ACGT", k=BC_LEN))
        barcodes.setdefault(bc, f"S{len(barcodes):03d}")
    tsv = os.path.join(work, "barcodes.tsv")
    with open(tsv, "w", encoding="utf-8", newline="") as fh:
        fh.write("barcode\tsample\n")
        for bc, sample in barcodes.items():
            fh.write(f"{bc}\t{sample}\n")
    bc_list = sorted(barcodes)

    r1_path = os.path.join(work, "R1.fastq.gz")
    r2_path = os.path.join(work, "R2.fastq.gz")
    qual1 = "I" * (BC_LEN + UMI_LEN + INSERT_LEN)
    qual2 = "I" * R2_LEN
    with gzip.GzipFile(r1_path, "wb", mtime=0) as f1, \
            gzip.GzipFile(r2_path, "wb", mtime=0) as f2:
        for i in range(pairs):
            roll = rng.random()
            if roll < 0.85:          # 登记条码（其中 10% 带一个错配）
                bc = bc_list[i % samples]
                if rng.random() < 0.10:
                    pos = rng.randrange(BC_LEN)
                    alt = rng.choice(
                        [b for b in "ACGT" if b != bc[pos]])
                    bc = bc[:pos] + alt + bc[pos + 1:]
            elif roll < 0.90:        # 观测条码含 N
                bc = bc_list[i % samples][:4] + "N" + \
                    "".join(rng.choices("ACGT", k=3))
            else:                    # 随机条码：绝大多数无法分配
                bc = "".join(rng.choices("ACGT", k=BC_LEN))
            umi = "".join(rng.choices("ACGT", k=UMI_LEN))
            ins = "".join(rng.choices("ACGT", k=INSERT_LEN))
            r2 = "".join(rng.choices("ACGT", k=R2_LEN))
            name = f"bench{i:07d}"
            f1.write(f"@{name}/1\n{bc}{umi}{ins}\n+\n{qual1}\n"
                     .encode("ascii"))
            f2.write(f"@{name}/2\n{r2}\n+\n{qual2}\n".encode("ascii"))
    return tsv, r1_path, r2_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", type=int, default=120000)
    ap.add_argument("--samples", type=int, default=48)
    ap.add_argument("--max-open", type=int, default=8)
    ap.add_argument("--keep", action="store_true",
                    help="保留临时目录（默认运行完删除）")
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix="demux_bench_")
    try:
        t0 = time.perf_counter()
        tsv, r1_path, r2_path = generate(args.pairs, args.samples, work)
        gen_s = time.perf_counter() - t0
        in_mb = (os.path.getsize(r1_path) + os.path.getsize(r2_path)) / 1e6

        t0 = time.perf_counter()
        summary = run_demux(r1_path, r2_path, tsv,
                            os.path.join(work, "out"),
                            umi_length=UMI_LEN, max_mismatches=1,
                            max_open=args.max_open)
        demux_s = time.perf_counter() - t0
        assert summary["status"] == "success", summary["error"]
        assert summary["pairs_completed"] == args.pairs

        # 抽查：所有去向的对数之和 == 输入对数；最大去向文件可完整读回
        total = sum(summary["destinations"].values())
        assert total == args.pairs, (total, args.pairs)
        top = max(summary["destinations"], key=summary["destinations"].get)
        n_back = 0
        with gzip.open(os.path.join(
                work, "out", "fastq", top, "R1.fastq.gz"), "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    n_back += 1
        assert n_back == summary["destinations"][top]

        rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        report = {
            "pairs": args.pairs,
            "samples": args.samples,
            "max_open": args.max_open,
            "input_gzip_mb": round(in_mb, 1),
            "generate_seconds": round(gen_s, 2),
            "demux_seconds": round(demux_s, 2),
            "pairs_per_second": int(args.pairs / demux_s),
            "max_rss_mb": round(rss_mb, 1),
            "assigned": total - summary["destinations"].get(
                "undetermined", 0),
            "undetermined": summary["destinations"].get("undetermined", 0),
            "unassigned_reasons": summary["unassigned_reasons"],
            "destinations_count": len(summary["destinations"]),
            "work_dir": work if args.keep else "(已清理)",
        }
        print(json.dumps(report, indent=2, ensure_ascii=False))
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print(f"临时目录保留于 {work}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
