"""命令行入口。

用法
----
python -m vcfpipe.cli run \
    --fasta data/ref/syn37.fa --contigs data/ref/contigs.json \
    --samples data/input/sample_sheet.csv --priority data/input/caller_priority.json \
    --out-root output --rule-version v1
"""
from __future__ import annotations

import argparse
import json

from .pipeline import VersionConflict, run_pipeline
from .query import Dataset


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="vcfpipe")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="运行管线")
    run.add_argument("--fasta", required=True)
    run.add_argument("--contigs", required=True)
    run.add_argument("--samples", required=True)
    run.add_argument("--priority", required=True)
    run.add_argument("--out-root", required=True)
    run.add_argument("--rule-version", default="v1", choices=["v1", "v2"])

    q = sub.add_parser("query", help="对已生成数据集做只读查询")
    q.add_argument("--dataset", required=True)
    q.add_argument("--mode", required=True,
                   choices=["range", "sample", "trace-output",
                            "trace-source", "quarantine", "matrix"])
    q.add_argument("--chrom")
    q.add_argument("--start", type=int)
    q.add_argument("--end", type=int)
    q.add_argument("--sample")
    q.add_argument("--ref")
    q.add_argument("--alt")
    q.add_argument("--source-path")
    q.add_argument("--line-no", type=int)
    q.add_argument("--reason")

    args = parser.parse_args(argv)

    if args.cmd == "run":
        try:
            manifest = run_pipeline(
                args.out_root, args.fasta, args.contigs, args.samples,
                args.priority, args.rule_version)
        except VersionConflict as exc:
            parser.exit(2, f"版本冲突: {exc}\n")
        print(json.dumps({
            "dataset_id": manifest["dataset_id"],
            "totals": manifest["totals"],
        }, indent=2, ensure_ascii=False))
        return 0

    ds = Dataset(args.dataset)
    if args.mode == "range":
        if not args.chrom or args.start is None:
            parser.error("range 模式需要 --chrom 与 --start")
        result = ds.query_range(args.chrom, args.start, args.end)
    elif args.mode == "sample":
        if not args.sample:
            parser.error("sample 模式需要 --sample")
        result = ds.sample_slice(args.sample)
    elif args.mode == "trace-output":
        result = ds.trace_output(args.chrom, args.start, args.ref, args.alt)
    elif args.mode == "trace-source":
        result = ds.trace_source(args.source_path, args.line_no)
    elif args.mode == "quarantine":
        result = ds.quarantine(args.reason)
    else:
        result = ds.matrix_wide()
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
