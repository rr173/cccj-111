"""端到端管线编排与确定性输出。

阶段
----
1. 载入参考 FASTA、contig 字典、样本表、caller 优先级；
2. 解析每份 VCF 的显式头声明并校验样本列；
3. 逐条记录：contig 别名解析 -> 符号 ALT 检查 -> 参考碱基校验
   -> 多等位拆分 -> 最左对齐/最小表示；失败进入隔离清单；
4. 按 (assembly, chrom, pos, REF, ALT, sample) 合并，
   caller 优先级显式可配，冲突保留全部来源证据；
5. 写分区事实表、样本矩阵、证据表、追踪表、合并 VCF、清单。

版本保护
--------
数据集目录名为 ``{target_assembly}_{rule_version}``。同目录重跑必须
逐字节得到相同清单，否则抛出 :class:`VersionConflict`；换 assembly 或
规则版本会得到新的 dataset_id，绝不原地覆盖旧版本。
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass

from .contigs import ContigMapping
from .fasta import Reference
from .genotype import split_record_calls
from .model import (
    Evidence, MergedVariant, QuarantineRecord, SplitRecord,
    REASON_BAD_GT, REASON_BAD_RECORD, REASON_BOUNDARY,
    REASON_CONTIG_UNMAPPED, REASON_REF_MISMATCH, REASON_SYMBOLIC_ALT,
)
from .normalize import RuleBundle, normalize_variant
from .vcf_io import MISSING, iter_records, parse_source, write_vcf

SCHEMA_VERSION = "1.0"


class VersionConflict(RuntimeError):
    """同 dataset_id 下已有清单与本次结果不一致。"""


# --------------------------------------------------------------------------- #
# 小工具
# --------------------------------------------------------------------------- #
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_sample_sheet(path: str) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    required = {"sample_id", "batch", "vcf_path"}
    for i, row in enumerate(rows, 2):
        missing = required - set(row)
        if missing:
            raise ValueError(f"样本表第 {i} 行缺列: {sorted(missing)}")
    return rows


def load_caller_priority(path: str) -> dict[str, int]:
    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    order = raw["caller_priority"]  # 索引 0 = 最高优先级
    return {caller: i for i, caller in enumerate(order)}


# --------------------------------------------------------------------------- #
# 核心处理
# --------------------------------------------------------------------------- #
@dataclass
class ProcessResult:
    included: list[SplitRecord]
    quarantined: list[QuarantineRecord]
    lineage: list[dict[str, str]]
    sources: list


def _quarantine(rec, reason_code, detail, lineage, alts=None):
    q = QuarantineRecord(
        source_path=rec.source.path, digest=rec.source.digest,
        line_no=rec.line_no, assembly=rec.source.assembly,
        chrom=rec.chrom, pos=rec.pos, ref=rec.ref,
        alts=alts if alts is not None else rec.alts,
        raw_gt=rec.raw_gt, reason_code=reason_code, reason_detail=detail,
    )
    lineage.append({
        "source_path": rec.source.path, "source_line": str(rec.line_no),
        "source_chrom": rec.chrom, "source_pos": str(rec.pos),
        "source_ref": rec.ref,
        "source_alt": ",".join(q.alts) if q.alts else MISSING,
        "output_chrom": "", "output_pos": "", "output_ref": "",
        "output_alt": "", "sample": rec.source.sample,
        "status": "quarantined", "reason": reason_code,
    })
    return q


def process_sources(
    sample_rows: list[dict[str, str]],
    reference: Reference,
    mapping: ContigMapping,
    rule: RuleBundle,
) -> ProcessResult:
    included: list[SplitRecord] = []
    quarantined: list[QuarantineRecord] = []
    lineage: list[dict[str, str]] = []
    sources = []

    for row in sample_rows:
        vcf_path = row["vcf_path"]
        digest = sha256_file(vcf_path)
        header = parse_source(vcf_path, digest)
        src = header.source
        sources.append(src)

        if header.sample_columns != (row["sample_id"],):
            raise ValueError(
                f"{vcf_path}: 头声明/样本列 {header.sample_columns} "
                f"与样本表 {row['sample_id']} 不一致"
            )

        for rec in iter_records(header):
            # 1) contig 别名 -> 目标 assembly
            target_chrom = mapping.resolve(src.assembly, rec.chrom)
            if target_chrom is None or not reference.has_contig(target_chrom):
                quarantined.append(_quarantine(
                    rec, REASON_CONTIG_UNMAPPED,
                    f"contig {src.assembly}:{rec.chrom} 无法唯一映射到 "
                    f"{mapping.target_assembly}", lineage))
                continue

            # 2) 符号/非碱基 ALT（<DEL>、* 等）本管线不支持
            if any(not a or a[0] in "<*" or not all(
                    c in "ACGTN" for c in a) for a in rec.alts):
                quarantined.append(_quarantine(
                    rec, REASON_SYMBOLIC_ALT,
                    f"不支持的符号/非碱基 ALT: {','.join(rec.alts)}", lineage))
                continue
            if not all(c in "ACGTN" for c in rec.ref):
                quarantined.append(_quarantine(
                    rec, REASON_BAD_RECORD, f"REF 含非法碱基: {rec.ref}",
                    lineage))
                continue

            # 3) 参考碱基校验（在原始坐标/原始 REF 上）
            ok, actual = reference.verify_ref(target_chrom, rec.pos, rec.ref)
            if not ok:
                quarantined.append(_quarantine(
                    rec, REASON_REF_MISMATCH,
                    f"REF={rec.ref} 但参考为 {actual}", lineage))
                continue

            # 4) GT 结构校验：等位下标必须落在 [0, ALT 数]，缺失除外
            if rec.raw_gt is not None and rec.raw_gt != MISSING:
                sep = "|" if rec.phased else "/"
                try:
                    gt_idx = [None if p == MISSING else int(p)
                              for p in rec.raw_gt.split(sep)]
                except ValueError:
                    quarantined.append(_quarantine(
                        rec, REASON_BAD_GT,
                        f"无法解析 GT={rec.raw_gt}", lineage))
                    continue
                bad = [i for i in gt_idx if i is not None
                       and not 0 <= i <= len(rec.alts)]
                if bad:
                    quarantined.append(_quarantine(
                        rec, REASON_BAD_GT,
                        f"GT={rec.raw_gt} 含超出 ALT 数 {len(rec.alts)} 的下标",
                        lineage))
                    continue

            # 5) 多等位拆分（GT/AD/DP/相位/缺失语义保留）
            calls = split_record_calls(rec)

            for alt_index, (alt, call) in enumerate(
                    zip(rec.alts, calls), start=1):
                norm, reason = normalize_variant(
                    reference, target_chrom, rec.pos, rec.ref, alt, rule)
                if norm is None:
                    code = REASON_BOUNDARY if reason == REASON_BOUNDARY \
                        else REASON_BAD_RECORD
                    quarantined.append(_quarantine(
                        rec, code,
                        f"ALT#{alt_index}={alt} 规范化失败: {reason}",
                        lineage, alts=(alt,)))
                    continue

                split = SplitRecord(
                    source=src, line_no=rec.line_no,
                    assembly=mapping.target_assembly, chrom=target_chrom,
                    pos=norm.pos, ref=norm.ref, alt=norm.alt,
                    qual=rec.qual, vfilter=rec.vfilter,
                    info=rec.info + (";BOUNDARY" if norm.boundary_flag else ""),
                    call=call,
                )
                included.append(split)
                lineage.append({
                    "source_path": src.path,
                    "source_line": str(rec.line_no),
                    "source_chrom": rec.chrom, "source_pos": str(rec.pos),
                    "source_ref": rec.ref, "source_alt": alt,
                    "output_chrom": target_chrom, "output_pos": str(norm.pos), "output_ref": norm.ref,
                    "output_alt": norm.alt, "sample": src.sample,
                    "status": "included", "reason": "",
                })

    return ProcessResult(included, quarantined, lineage, sources)


# --------------------------------------------------------------------------- #
# 合并
# --------------------------------------------------------------------------- #
def merge_records(
    included: list[SplitRecord], priority: dict[str, int]
) -> list[MergedVariant]:
    groups: dict[tuple, list[SplitRecord]] = {}
    for rec in included:
        groups.setdefault(rec.key() + (rec.source.sample,), []).append(rec)

    merged: list[MergedVariant] = []
    for group_key, recs in groups.items():
        key, sample = group_key[:-1], group_key[-1]
        ordered = sorted(
            recs,
            key=lambda r: (
                priority.get(r.source.caller, len(priority)),
                r.source.digest, r.line_no,
            ),
        )
        winner = ordered[0]
        evidence = tuple(
            Evidence(
                source_path=r.source.path, digest=r.source.digest,
                line_no=r.line_no, caller=r.source.caller, qual=r.qual,
                vfilter=r.vfilter, info=r.info,
                raw_gt=_alleles_to_raw(r.call),
                raw_ad=r.call.ad, raw_dp=r.call.dp, phased=r.call.phased,
            )
            for r in ordered
        )
        merged.append(MergedVariant(
            assembly=key[0], chrom=key[1], pos=key[2], ref=key[3], alt=key[4],
            sample=sample, winner_caller=winner.source.caller,
            call=winner.call, evidence=evidence,
        ))
    merged.sort(key=lambda m: (m.chrom, m.pos, m.ref, m.alt, m.sample))
    return merged


def _alleles_to_raw(call) -> str:
    from .vcf_io import gt_to_string
    return gt_to_string(call)


# --------------------------------------------------------------------------- #
# 写出
# --------------------------------------------------------------------------- #
def _write_tsv(path: str, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _gt_cell(call) -> str:
    from .vcf_io import gt_to_string
    return gt_to_string(call)


def write_outputs(
    out_dir: str, merged: list[MergedVariant],
    quarantined: list[QuarantineRecord], lineage: list[dict[str, str]],
    sources: list, reference: Reference, mapping: ContigMapping,
    rule: RuleBundle, priority: dict[str, int],
) -> dict[str, str]:
    os.makedirs(os.path.join(out_dir, "variants"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "matrix"), exist_ok=True)
    written: dict[str, str] = {}

    # 1) 分区事实表（按染色体分区，区内按 pos 排序）+ 范围索引
    by_chrom: dict[str, list[MergedVariant]] = {}
    for m in merged:
        by_chrom.setdefault(m.chrom, []).append(m)

    index_rows = []
    for chrom in sorted(by_chrom):
        part = sorted(by_chrom[chrom], key=lambda m: (m.pos, m.alt, m.sample))
        part_path = os.path.join("variants", f"variants_{chrom}.tsv")
        rows = [[
            m.assembly, m.chrom, str(m.pos), m.ref, m.alt, m.sample,
            m.winner_caller, _gt_cell(m.call),
            "." if m.call.ad is None else
            ",".join("." if x is None else str(x) for x in m.call.ad),
            "." if m.call.dp is None else str(m.call.dp),
            "1" if m.call.phased else "0",
        ] for m in part]
        _write_tsv(
            os.path.join(out_dir, part_path),
            ["assembly", "chrom", "pos", "ref", "alt", "sample",
             "winner_caller", "gt", "ad", "dp", "phased"],
            rows,
        )
        index_rows.append([
            chrom, part_path, str(part[0].pos), str(part[-1].pos), str(len(part)),
        ])
        written[part_path] = sha256_file(os.path.join(out_dir, part_path))
    _write_tsv(
        os.path.join(out_dir, "variants", "_index.tsv"),
        ["chrom", "path", "min_pos", "max_pos", "record_count"],
        sorted(index_rows),
    )
    written["variants/_index.tsv"] = sha256_file(
        os.path.join(out_dir, "variants", "_index.tsv"))

    # 2) 样本矩阵：长表 + 宽表（GT；NA = 该样本在此位点无记录）
    long_rows = [[
        m.assembly, m.chrom, str(m.pos), m.ref, m.alt, m.sample,
        _gt_cell(m.call),
        "." if m.call.ad is None else
        ",".join("." if x is None else str(x) for x in m.call.ad),
        "." if m.call.dp is None else str(m.call.dp),
        "1" if m.call.phased else "0", m.winner_caller,
    ] for m in merged]
    _write_tsv(
        os.path.join(out_dir, "matrix", "genotypes_long.tsv"),
        ["assembly", "chrom", "pos", "ref", "alt", "sample", "gt",
         "ad", "dp", "phased", "winner_caller"], long_rows,
    )
    written["matrix/genotypes_long.tsv"] = sha256_file(
        os.path.join(out_dir, "matrix", "genotypes_long.tsv"))

    sample_names = sorted({m.sample for m in merged})
    keys = sorted({m.key() for m in merged})
    cell: dict[tuple, str] = {}
    for m in merged:
        cell[(m.key(), m.sample)] = _gt_cell(m.call)
    wide_rows = []
    for key in keys:
        wide_rows.append(list(key[:1]) + list(key[1:]) + [
            cell.get((key, s), "NA") for s in sample_names])
    _write_tsv(
        os.path.join(out_dir, "matrix", "sample_gt_matrix.tsv"),
        ["assembly", "chrom", "pos", "ref", "alt"] + sample_names, wide_rows,
    )
    written["matrix/sample_gt_matrix.tsv"] = sha256_file(
        os.path.join(out_dir, "matrix", "sample_gt_matrix.tsv"))

    # 3) 冲突证据表
    ev_rows = []
    for m in merged:
        for e in m.evidence:
            ev_rows.append([
                m.assembly, m.chrom, str(m.pos), m.ref, m.alt, m.sample,
                e.caller, e.source_path, str(e.line_no), e.qual, e.vfilter,
                e.info, e.raw_gt or MISSING,
                MISSING if e.raw_ad is None else
                ",".join(MISSING if x is None else str(x) for x in e.raw_ad),
                MISSING if e.raw_dp is None else str(e.raw_dp),
                "1" if e.phased else "0",
            ])
    _write_tsv(
        os.path.join(out_dir, "evidence.tsv"),
        ["assembly", "chrom", "pos", "ref", "alt", "sample", "caller",
         "source_path", "source_line", "qual", "filter", "info", "raw_gt",
         "raw_ad", "raw_dp", "phased"], ev_rows,
    )
    written["evidence.tsv"] = sha256_file(os.path.join(out_dir, "evidence.tsv"))

    # 4) source-to-output 追踪
    lin_header = ["source_path", "source_line", "source_chrom", "source_pos",
                  "source_ref", "source_alt", "output_chrom", "output_pos",
                  "output_ref", "output_alt", "sample", "status", "reason"]
    _write_tsv(
        os.path.join(out_dir, "lineage.tsv"), lin_header,
        [[r[h] for h in lin_header] for r in sorted(
            lineage, key=lambda r: (r["source_path"], int(r["source_line"]),
                                    r["source_alt"]))],
    )
    written["lineage.tsv"] = sha256_file(os.path.join(out_dir, "lineage.tsv"))

    # 5) 隔离清单
    q_rows = [[
        q.source_path, q.digest, str(q.line_no), q.assembly, q.chrom,
        str(q.pos), q.ref, ",".join(q.alts), q.raw_gt or MISSING,
        q.reason_code, q.reason_detail,
    ] for q in sorted(
        quarantined, key=lambda q: (q.source_path, q.line_no, q.alts))]
    _write_tsv(
        os.path.join(out_dir, "quarantine.tsv"),
        ["source_path", "source_digest", "source_line", "assembly", "chrom",
         "pos", "ref", "alts", "raw_gt", "reason_code", "reason_detail"],
        q_rows,
    )
    written["quarantine.tsv"] = sha256_file(
        os.path.join(out_dir, "quarantine.tsv"))

    # 6) 合并 VCF
    vcf_path = os.path.join(out_dir, "merged.vcf")
    write_vcf(vcf_path, mapping.target_assembly, merged, sources)
    written["merged.vcf"] = sha256_file(vcf_path)
    return written


# --------------------------------------------------------------------------- #
# 清单与入口
# --------------------------------------------------------------------------- #
def build_manifest(
    dataset_id: str, mapping: ContigMapping, rule: RuleBundle,
    priority: dict[str, int], reference: Reference, sources: list,
    result: ProcessResult, written: dict[str, str],
    sample_sheet_path: str, contig_dict_path: str,
) -> dict:
    per_source: dict[str, dict[str, int]] = {}
    for rec in result.included:
        d = per_source.setdefault(rec.source.path,
                                  {"included": 0, "quarantined": 0})
        d["included"] += 1
    for q in result.quarantined:
        d = per_source.setdefault(q.source_path,
                                  {"included": 0, "quarantined": 0})
        d["quarantined"] += 1

    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": dataset_id,
        "target_assembly": mapping.target_assembly,
        "rule_version": rule.rule_version,
        "boundary_policy": rule.boundary_policy,
        "caller_priority": [c for c, _ in sorted(priority.items(),
                                                 key=lambda kv: kv[1])],
        "reference": {
            "path": reference.path,
            "file_sha256": reference.digest.file_sha256,
            "contig_lengths": reference.digest.contig_lengths,
            "contig_md5": reference.digest.contig_md5,
        },
        "contig_dict": {
            "path": contig_dict_path,
            "file_sha256": sha256_file(contig_dict_path),
        },
        "sample_sheet": {
            "path": sample_sheet_path,
            "file_sha256": sha256_file(sample_sheet_path),
        },
        "sources": sorted([
            {
                "path": s.path,
                "file_sha256": s.digest,
                "assembly": s.assembly,
                "sample": s.sample,
                "caller": s.caller,
                "summary": s.summary,
                "included": per_source.get(s.path, {}).get("included", 0),
                "quarantined": per_source.get(s.path, {}).get(
                    "quarantined", 0),
            } for s in sources
        ], key=lambda s: s["path"]),
        "totals": {
            "included_records": len(result.included),
            "quarantined_records": len(result.quarantined),
            # merged_variant_rows 由调用方在合并后回填
            "merged_variant_rows": 0,
        },
        "outputs": dict(sorted(written.items())),
    }


def run_pipeline(
    out_root: str, fasta_path: str, contig_dict_path: str,
    sample_sheet_path: str, priority_path: str, rule_version: str = "v1",
) -> dict:
    reference = Reference(fasta_path)
    mapping = ContigMapping.load(contig_dict_path)
    rule = RuleBundle.for_version(rule_version)
    priority = load_caller_priority(priority_path)
    sample_rows = load_sample_sheet(sample_sheet_path)

    dataset_id = f"{mapping.target_assembly}_{rule.rule_version}"
    out_dir = os.path.join(out_root, dataset_id)
    os.makedirs(out_dir, exist_ok=True)

    result = process_sources(sample_rows, reference, mapping, rule)
    merged = merge_records(result.included, priority)
    written = write_outputs(
        out_dir, merged, result.quarantined, result.lineage, result.sources,
        reference, mapping, rule, priority)

    manifest = build_manifest(
        dataset_id, mapping, rule, priority, reference, result.sources,
        result, written, sample_sheet_path, contig_dict_path)
    # merged_variant_rows 需要在合并后修正
    manifest["totals"]["merged_variant_rows"] = len(merged)

    manifest_path = os.path.join(out_dir, "manifest.json")
    # 统一以“渲染文本 + 单个换行”作为落盘形态，比较两侧都规范化，
    # 避免历史文件的尾换行造成假性版本冲突。
    payload = (json.dumps(manifest, indent=2, sort_keys=True,
                          ensure_ascii=False) + "\n").encode("utf-8")
    if os.path.exists(manifest_path):
        with open(manifest_path, "rb") as fh:
            old = fh.read()
        if old != payload:
            raise VersionConflict(
                f"{manifest_path} 已存在且内容不同；拒绝原地覆盖。"
                f"更换 assembly/规则版本会生成新 dataset_id。")
    with open(manifest_path, "wb") as fh:
        fh.write(payload)
    return manifest
