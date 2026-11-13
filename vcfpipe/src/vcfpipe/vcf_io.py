"""VCF 解析与写出（标准库实现，无第三方依赖）。

输入 VCF 必须显式声明（缺失即报错，不做隐式推断）：
    ##ASSEMBLY=<...>     ##SAMPLE=<ID=...>     ##CALLER=<Name=...>
    ##FILE_SUMMARY=<...>

解析范围覆盖本管线需要的 FORMAT 字段：GT / AD / DP；
其余 FORMAT 字段被忽略但不会报错。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional

from .model import RawRecord, SourceVCF

MISSING = "."


@dataclass(frozen=True)
class VcfHeader:
    source: SourceVCF
    sample_columns: tuple[str, ...]
    format_keys: tuple[str, ...]
    raw_headers: tuple[str, ...]


# --------------------------------------------------------------------------- #
# 头解析
# --------------------------------------------------------------------------- #
def _meta_value(line: str, key: str) -> Optional[str]:
    """从 ``##KEY=<...>`` 头行抽取裸值或 ID/Name 属性。"""
    body = line[2:].rstrip(">")
    if "=" not in body:
        return None
    head, _, rest = body.partition("=")
    if head != key:
        return None
    rest = rest.strip()
    if rest.startswith("<"):
        rest = rest[1:]
        # 简单属性扫描（合成材料不含嵌套引号）
        for part in rest.rstrip(">").split(","):
            if "=" in part:
                k, _, v = part.partition("=")
                if k in ("ID", "Name", "Value"):
                    return v.strip('"')
        return None
    return rest


def parse_source(path: str, digest: str) -> VcfHeader:
    assembly = sample = caller = summary = None
    sample_columns: list[str] = []
    format_keys: list[str] = []
    contig_lengths: dict[str, int] = {}
    raw_headers: list[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("##"):
                raw_headers.append(line)
                if line.startswith("##ASSEMBLY=") or line.startswith("##ASSEMBLY<"):
                    assembly = _meta_value(line, "ASSEMBLY")
                elif line.startswith("##SAMPLE"):
                    sample = _meta_value(line, "SAMPLE")
                elif line.startswith("##CALLER"):
                    caller = _meta_value(line, "CALLER")
                elif line.startswith("##FILE_SUMMARY"):
                    summary = _meta_value(line, "FILE_SUMMARY")
                elif line.startswith("##FORMAT=<"):
                    fid = _meta_value(line, "FORMAT")
                    if fid:
                        format_keys.append(fid)
                elif line.startswith("##contig=<"):
                    body = line[len("##contig=<"):].rstrip(">")
                    attrs = dict(
                        p.split("=", 1) for p in body.split(",") if "=" in p
                    )
                    if "ID" in attrs:
                        contig_lengths[attrs["ID"]] = int(
                            attrs.get("length", "0")
                        )
                continue
            if line.startswith("#CHROM"):
                cols = line[1:].split("\t")
                if len(cols) < 8:
                    raise ValueError(f"{path}: #CHROM 头列数不足")
                sample_columns = cols[9:] if len(cols) > 9 else []
                break
            raise ValueError(f"{path}: 遇到非头行且缺少 #CHROM: {line!r}")

    missing = [
        n for n, v in (
            ("ASSEMBLY", assembly), ("SAMPLE", sample),
            ("CALLER", caller), ("FILE_SUMMARY", summary),
        ) if v is None
    ]
    if missing:
        raise ValueError(f"{path}: VCF 缺少显式声明头: {', '.join(missing)}")
    if not sample_columns:
        raise ValueError(f"{path}: VCF 必须包含至少一个样本列")

    source = SourceVCF(
        path=path, assembly=assembly, sample=sample, caller=caller,
        digest=digest, summary=summary, contig_lengths=contig_lengths,
    )
    return VcfHeader(source, tuple(sample_columns), tuple(format_keys),
                     tuple(raw_headers))


# --------------------------------------------------------------------------- #
# FORMAT 字段
# --------------------------------------------------------------------------- #
def _parse_gt(token: str) -> tuple[Optional[tuple[Optional[int], ...]], bool]:
    """解析 GT。返回 (等位元组或 None 表示整体缺失, 是否相位化)。"""
    if token == MISSING:
        return None, False
    phased = "|" in token
    sep = "|" if phased else "/"
    alleles: list[Optional[int]] = []
    for part in token.split(sep):
        if part == MISSING:
            alleles.append(None)
        else:
            alleles.append(int(part))
    return tuple(alleles), phased


def _parse_ad(token: str) -> Optional[tuple[Optional[int], ...]]:
    if token in (MISSING, ""):
        return None
    out: list[Optional[int]] = []
    for part in token.split(","):
        out.append(None if part == MISSING else int(part))
    return tuple(out)


def _parse_dp(token: str) -> Optional[int]:
    if token in (MISSING, ""):
        return None
    return int(token)


# --------------------------------------------------------------------------- #
# 记录解析
# --------------------------------------------------------------------------- #
def iter_records(header: VcfHeader) -> Iterator[RawRecord]:
    src = header.source
    with open(src.path, "r", encoding="utf-8") as fh:
        line_no = 0
        in_body = False
        for line in fh:
            line_no += 1
            line = line.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                in_body = True
                continue
            if not in_body:
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                raise ValueError(f"{src.path}:{line_no} 列数不足 8")
            chrom, pos_s, _id, ref, alt_s, qual, vfilter, info = cols[:8]
            alts = tuple() if alt_s == MISSING else tuple(alt_s.split(","))
            fmt = cols[8].split(":") if len(cols) > 8 and cols[8] != MISSING else []
            sample_vals = cols[9].split(":") if len(cols) > 9 else []
            values = dict(zip(fmt, sample_vals))

            raw_gt = values.get("GT")
            gt_alleles = None
            phased = False
            if raw_gt is not None:
                gt_alleles, phased = _parse_gt(raw_gt)
            raw_ad = _parse_ad(values["AD"]) if "AD" in values else None
            raw_dp = _parse_dp(values["DP"]) if "DP" in values else None

            yield RawRecord(
                source=src, line_no=line_no, chrom=chrom, pos=int(pos_s),
                ref=ref.upper(), alts=tuple(a.upper() for a in alts),
                qual=qual, vfilter=vfilter, info=info,
                raw_gt=raw_gt, raw_ad=raw_ad, raw_dp=raw_dp, phased=phased,
            )


# --------------------------------------------------------------------------- #
# 写出
# --------------------------------------------------------------------------- #
def gt_to_string(call) -> str:
    """GenotypeCall -> VCF GT 字符串；缺失语义完整保留。"""
    if call.gt_alleles is None:
        return "./." if not call.phased else ".|."
    sep = "|" if call.phased else "/"
    return sep.join(MISSING if a is None else str(a) for a in call.gt_alleles)


def ad_to_string(ad) -> str:
    if ad is None:
        return MISSING
    return ",".join(MISSING if x is None else str(x) for x in ad)


def write_vcf(path: str, assembly: str, merged: list, sources: list[SourceVCF]) -> None:
    """写出合并后的多样本 VCF（样本与位点均排序，逐字节确定性输出）。

    注意：不能使用集合或依赖字典插入序，否则不同运行间头/行顺序会漂移。
    """
    sample_names = tuple(sorted({m.sample for m in merged}))
    # (key) -> {sample: MergedVariant}，键按稳定元组排序
    by_key: dict[tuple, dict[str, object]] = {}
    for m in merged:
        by_key.setdefault(m.key(), {})[m.sample] = m

    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f'##ASSEMBLY=<Value={assembly}>\n')
        fh.write('##INFO=<ID=WINNER,Number=1,Type=String,Description="Winning caller">\n')
        fh.write('##INFO=<ID=EVIDENCE,Number=.,Type=String,Description="Conflicting caller evidence: caller|rawGT|AD|DP">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        # 来源头按 (caller, path, digest) 去重排序
        seen = set()
        for s in sorted(sources, key=lambda x: (x.caller, x.path, x.digest)):
            stamp = (s.caller, s.path, s.digest)
            if stamp in seen:
                continue
            seen.add(stamp)
            fh.write(f'##source=<ID={s.caller},Digest={s.digest},File={s.path}>\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(sample_names) + "\n")

        def _ev_string(e) -> str:
            ad_part = (
                "+".join(MISSING if x is None else str(x) for x in e.raw_ad)
                if e.raw_ad else MISSING
            )
            return (f"{e.caller}|{e.raw_gt or MISSING}|{ad_part}|"
                    f"{e.raw_dp if e.raw_dp is not None else MISSING}")

        for key in sorted(by_key):
            per_sample = by_key[key]
            _assembly, chrom, pos, ref, alt = key
            rows: list[str] = []
            for sample in sample_names:
                m = per_sample.get(sample)
                if m is None:
                    rows.append("./.")
                    continue
                gt = gt_to_string(m.call)
                ad = ad_to_string(m.call.ad)
                dp = MISSING if m.call.dp is None else str(m.call.dp)
                rows.append(f"{gt}:{ad}:{dp}")
            # 以样本名排序后的第一个有记录样本决定 WINNER
            first = next(per_sample[s] for s in sample_names if s in per_sample)
            info = (
                f"WINNER={first.winner_caller};EVIDENCE="
                + ",".join(_ev_string(e) for e in first.evidence)
            )
            fh.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT:AD:DP\t"
                + "\t".join(rows) + "\n"
            )
