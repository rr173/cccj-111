"""命令行入口。

用法
----
python -m vcfpipe.cli run \
    --fasta data/ref/syn37.fa --contigs data/ref/contigs.json \
    --samples data/input/sample_sheet.csv --priority data/input/caller_priority.json \
    --out-root output --rule-version v1

python -m vcfpipe.cli annot \
    --fasta data/ref/annot_syn.fa \
    --transcripts data/annot/transcripts_demo.tsv \
    --variants data/annot/variants_demo.tsv \
    --out-dir output/annot_demo --format tsv

python -m vcfpipe.cli demux \
    --r1 data/demux/R1.fastq.gz --r2 data/demux/R2.fastq.gz \
    --barcodes data/demux/barcodes.tsv --out-dir output/demux_demo \
    --umi-length 4 --max-mismatches 1 --max-open 4
"""
from __future__ import annotations

import argparse
import json
import sys

from .annot import AnnotInputError, run_annot
from .assoc import AssocInputError, run_assoc, to_jsonable, write_report_files
from .demux import (
    DemuxInputError, MAX_MAX_OPEN, MAX_UMI_LENGTH, MIN_MAX_OPEN,
    MIN_UMI_LENGTH, run_demux,
)
from .pipeline import VersionConflict, run_pipeline
from .query import Dataset
from .trio import TrioConfigError, TrioRecordError, run_trio


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

    d = sub.add_parser("demux", help="测序数据接收：校验双端 FASTQ、"
                                     "按条码拆分混合文库、提取 UMI")
    d.add_argument("--r1", required=True, help="R1 FASTQ（普通或 gzip）")
    d.add_argument("--r2", required=True, help="R2 FASTQ（普通或 gzip）")
    d.add_argument("--barcodes", required=True,
                   help="条码登记表（TSV，表头 barcode/sample）")
    d.add_argument("--out-dir", required=True,
                   help="接收输出目录（必须不存在或为空）")
    d.add_argument("--umi-length", required=True, type=int, metavar="U",
                   help=f"UMI 长度，合法范围 [{MIN_UMI_LENGTH}, "
                        f"{MAX_UMI_LENGTH}]")
    d.add_argument("--max-mismatches", required=True, type=int, metavar="M",
                   help="允许的条码错配数，合法范围 [0, 登记条码长度]")
    d.add_argument("--max-open", required=True, type=int, metavar="K",
                   help=f"同时打开的去向输出对上限（undetermined 也占一个），"
                        f"合法范围 [{MIN_MAX_OPEN}, {MAX_MAX_OPEN}]")

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

    a = sub.add_parser("assoc", help="基因级携带者关联分析"
                                     "（只读已发布数据集，不重新规范化）")
    a.add_argument("--dataset", required=True, help="已发布数据集目录")
    a.add_argument("--labels", required=True,
                   help="病例/对照标签表（TSV 表头: sample, group）")
    a.add_argument("--genes", required=True,
                   help="基因区间 BED（chrom/start/end/gene，零起点半开）")
    a.add_argument("--callable", required=True, dest="callable_bed",
                   help="每样本可判读区间 BED（chrom/start/end/sample）")
    a.add_argument("--min-dp", required=True, type=int,
                   help="有效调用所需的最小 DP（达到该值才有效）")
    a.add_argument("--out-dir",
                   help="可选：写出 assoc_report.json 与汇总/明细 TSV")

    an = sub.add_parser(
        "annot",
        help="转录本编码后果注释（本地合成材料；从变异走到氨基酸）")
    an.add_argument("--fasta", required=True)
    an.add_argument("--transcripts", required=True,
                    help="转录本 TSV（8 列表头，见 docs/README.md）")
    an.add_argument("--variants", required=True,
                    help="变异 TSV（id/chrom/pos/ref/alt 五列表头）")
    an.add_argument("--out-dir", required=True,
                    help="输出目录（必须不存在或为空）")
    an.add_argument("--format", default="tsv", choices=["tsv", "jsonl"],
                    help="解释表格式：逐行 TSV（默认）或 JSONL")

    t = sub.add_parser(
        "trio", help="家系三人组遗传一致性筛查（多样本 VCF 逐位点补全"
                     "枚举；只读不改写基因型，不做亲缘/临床结论）")
    t.add_argument("--vcf", required=True,
                   help="联合分型的多样本 VCF（普通或 gzip，按内容识别）")
    t.add_argument("--families", required=True,
                   help="家系表 TSV，表头恰为 family/child/father/mother")
    t.add_argument("--chroms", required=True,
                   help="常染色体名称白名单，逗号分隔；白名单外的行只"
                        "计跳过数（不把性染色体套二倍体规则）")
    t.add_argument("--min-dp", required=True, type=int,
                   help="每人每位点可信调用所需的最小 DP（正整数）；"
                        "DP 缺失或低于该值时整调用视为未知")
    t.add_argument("--out-dir", required=True,
                   help="输出目录（必须不存在或为空；staging 成功后原子发布）")

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

    if args.cmd == "annot":
        try:
            summary = run_annot(args.fasta, args.transcripts, args.variants,
                                args.out_dir, args.format)
        except AnnotInputError as exc:
            parser.exit(2, f"注释输入错误: {exc}\n")
        except OSError as exc:
            parser.exit(1, f"注释写出失败: {exc}\n")
        print(json.dumps(summary, indent=2, ensure_ascii=False,
                         sort_keys=True))
        return 0

    if args.cmd == "demux":
        try:
            summary = run_demux(
                args.r1, args.r2, args.barcodes, args.out_dir,
                umi_length=args.umi_length,
                max_mismatches=args.max_mismatches,
                max_open=args.max_open)
        except DemuxInputError as exc:
            parser.exit(2, f"接收配置错误: {exc}\n")
        except OSError as exc:
            parser.exit(1, f"接收失败（输出错误）: {exc}\n")
        print(json.dumps(summary, indent=2, ensure_ascii=False,
                         sort_keys=True))
        if summary["status"] != "success":
            print(f"接收失败: {summary['error']['message']}", file=sys.stderr)
            return 1
        return 0

    if args.cmd == "assoc":
        try:
            report = run_assoc(args.dataset, args.labels, args.genes,
                               args.callable_bed, args.min_dp)
        except AssocInputError as exc:
            parser.exit(2, f"关联分析输入错误: {exc}\n")
        if args.out_dir:
            write_report_files(report, args.out_dir)
        print(json.dumps(to_jsonable(report), indent=2, ensure_ascii=False))
        return 0

    if args.cmd == "trio":
        chroms = tuple(c.strip() for c in args.chroms.split(",") if c.strip())
        try:
            report = run_trio(
                args.vcf, args.families, chroms, args.min_dp, args.out_dir)
        except TrioConfigError as exc:
            parser.exit(2, f"家系筛查配置错误: {exc}\n")
        except TrioRecordError as exc:
            parser.exit(1, f"家系筛查记录错误: {exc}\n")
        print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
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
