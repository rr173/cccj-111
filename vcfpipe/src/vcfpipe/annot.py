"""从变异走到氨基酸：本地、确定性的转录本编码后果注释。

定位
----
研究者交来三份**合成/本地**材料——变异表（TSV）、转录本定义（TSV）
与参考 FASTA——本模块逐条给出“输入 id × ALT × 转录本”的解释：变异落
在哪个区域、若在编码区则改动哪个密码子的哪几个碱基、替换前后密码子
与氨基酸，以及后果类别。全部计算在本机完成，只用 Python 标准库，不
查询任何在线注释库，也不依赖第三方后果注释器。

输入契约（字段固定，便于其他程序调用）
-------------------------------------
1. **转录本 TSV**，表头恰为::

       transcript  gene  chrom  strand  exon_starts  exon_ends  cds_start  cds_end

   * 每行一个异构体；``strand`` ∈ ``+`` / ``-``；
   * ``exon_starts`` / ``exon_ends`` 为按坐标**递增**的逗号整数列表，
     成对形成**零起点、右开**的非空外显子区间，两个列表长度必须相等，
     外显子之间不得重叠（相邻允许）；
   * CDS 取“外显子 ∩ ``[cds_start, cds_end)``”；``cds_start == cds_end``
     表示非编码异构体；非空编码范围的两端（``cds_start`` 与
     ``cds_end - 1``）都必须落在某个外显子上；
   * 拼接后的 CDS 总长须整除 3（单个片段不必整除 3）。

2. **变异 TSV**，表头恰为::

       id  chrom  pos  ref  alt

   * ``pos`` 从 1 计数；``ref`` / ``alt`` 写**基因组正链**字母，限
     A/C/G/T 且非空；``alt`` 可用逗号给出多个选择，ALT 之间互异且
     均不得等于 REF；``id`` 在本表内唯一。

3. **FASTA**：支持折行，``>`` 头行第一个空白字段为序列名；碱基限
   A/C/G/T；序列名不得重复。

任何结构错误、缺失 contig、越界、外显子乱序/重叠、两个 exon 列表
长度不等、非法 CDS、非法碱基或 REF 与 FASTA 不符，都在**开写输出前**
以带“文件 + 行号”的 :class:`AnnotInputError` 结束，CLI 返回非零。

区域与后果
----------
对每条变异、每个 ALT、每个相交异构体（相交按“外显子最小起点至最大
终点的跨度”判断，且只用整个 REF 区间，不看锚点）：

* REF 完整落在同一外显子的**编码片段**内 -> ``coding``；
* REF 完整落在同一外显子的**非编码片段**内 -> ``utr``；
* 异构体本身非编码（``cds_start == cds_end``）且 REF 完整落在其外显子
  内 -> ``noncoding_exon``；
* REF 完整落在单一内含子内 -> ``intron``；
* REF 横跨外显子/内含子、编码/非编码，或越出转录本跨度 ->
  ``boundary_unsupported``；
* 没有任何相交异构体 -> ``intergenic``（transcript 为空）。

长度不等的 REF/ALT 统一优先给 ``unsupported_edit``：保留相交异构体，
一个都没有时 transcript 为空。编码后果只支持**等长碱基替换**：把整
段替换一次性施加到该转录本的 CDS，再翻译触及的密码子——同一密码子
改两处时答案由两处共同决定，绝不由两个单碱基结果拼接。

编码到达
--------
外显子按基因组正链递增给出。先取与编码范围相交的片段；正链按坐标
递增拼接，负链按坐标递减拼接各片段的**反向互补**。拼接完成后**从
头**按三联体切分，阅读框贯穿接缝。负链 ALT 反向互补后施加，输出
的三个基因组坐标按翻译方向排列。

后果类别（标准遗传密码表，``*`` 表终止）：
``synonymous``（含 ``*``→``*``）、``stop_gained``、``stop_lost``、
``missense``。不推测起始丢失、单倍型或致病性。

输出与排序
----------
``<out-dir>/annotation.tsv``（或 ``annotation.jsonl``）：**逐条写出
并逐变异 flush**，不在整批结束后才落盘；``summary.json`` 在成功完成
后写出。TSV 列固定，见 :data:`HEADER`。排序键:
变异在输入表中的原顺序 → ALT 列表顺序 → transcript 字典序 → 密码子
编号。非编码/暂不支持项同样占一行，密码子相关列为空并给出原因。
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass

# --------------------------------------------------------------------------- #
# 稳定标识符（输出的一部分，不随意改名）
# --------------------------------------------------------------------------- #
REGION_CODING = "coding"
REGION_UTR = "utr"
REGION_NONCODING_EXON = "noncoding_exon"
REGION_INTRON = "intron"
REGION_INTERGENIC = "intergenic"
REGION_BOUNDARY = "boundary_unsupported"
REGION_UNSUPPORTED_EDIT = "unsupported_edit"

REGIONS = (
    REGION_CODING, REGION_UTR, REGION_NONCODING_EXON, REGION_INTRON,
    REGION_INTERGENIC, REGION_BOUNDARY, REGION_UNSUPPORTED_EDIT,
)

CONSEQUENCE_SYNONYMOUS = "synonymous"
CONSEQUENCE_STOP_GAINED = "stop_gained"
CONSEQUENCE_STOP_LOST = "stop_lost"
CONSEQUENCE_MISSENSE = "missense"

CONSEQUENCES = (
    CONSEQUENCE_SYNONYMOUS, CONSEQUENCE_STOP_GAINED,
    CONSEQUENCE_STOP_LOST, CONSEQUENCE_MISSENSE,
)

TRANSCRIPT_HEADER = (
    "transcript", "gene", "chrom", "strand", "exon_starts", "exon_ends",
    "cds_start", "cds_end",
)
VARIANT_HEADER = ("id", "chrom", "pos", "ref", "alt")

HEADER = (
    "id", "chrom", "pos", "ref", "alt", "alt_index",
    "transcript", "gene", "strand", "region",
    "cds_position", "codon_number", "codon_pos_1", "codon_pos_2",
    "codon_pos_3", "codon_ref", "codon_alt", "aa_ref", "aa_alt",
    "consequence", "reason",
)

OUTPUT_TSV = "annotation.tsv"
OUTPUT_JSONL = "annotation.jsonl"
SUMMARY_JSON = "summary.json"

_BASES = set("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# 标准遗传密码表（NCBI table 1）。
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
assert len(CODON_TABLE) == 64


class AnnotInputError(ValueError):
    """输入结构/内容错误；消息自带来源文件与行号。"""


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_COMPLEMENT)


def translate_codon(codon: str) -> str:
    if len(codon) != 3 or codon not in CODON_TABLE:
        raise AnnotInputError(f"内部错误: 非法密码子 {codon!r}")
    return CODON_TABLE[codon]


# --------------------------------------------------------------------------- #
# 参考 FASTA
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Reference:
    path: str
    sequences: dict[str, str]
    file_sha256: str

    def has(self, chrom: str) -> bool:
        return chrom in self.sequences

    def length(self, chrom: str) -> int:
        return len(self.sequences[chrom])

    def slice1(self, chrom: str, start: int, end: int) -> str:
        """1-based 闭区间 [start, end]。"""
        return self.sequences[chrom][start - 1:end]


def load_fasta(path: str) -> Reference:
    if not os.path.isfile(path):
        raise AnnotInputError(f"FASTA 不存在: {path}")
    sequences: dict[str, str] = {}
    order: list[str] = []
    current: str | None = None
    header_line = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].split() else ""
                if not name:
                    raise AnnotInputError(
                        f"FASTA {path} 第 {line_no} 行: 头行缺少序列名")
                if name in sequences:
                    raise AnnotInputError(
                        f"FASTA {path} 第 {line_no} 行: 序列名 {name!r} "
                        f"重复（首次出现在第 {header_line} 行）")
                sequences[name] = ""
                order.append(name)
                current = name
                header_line = line_no
                continue
            if current is None:
                raise AnnotInputError(
                    f"FASTA {path} 第 {line_no} 行: 第一个 > 头之前出现"
                    f"序列行")
            seq = line.strip().upper()
            bad = sorted(set(seq) - _BASES)
            if bad:
                raise AnnotInputError(
                    f"FASTA {path} 第 {line_no} 行: contig {current} 含"
                    f"非法碱基 {''.join(bad)!r}（只接受 A/C/G/T）")
            sequences[current] += seq
    if not sequences:
        raise AnnotInputError(f"FASTA {path} 无任何 contig")
    with open(path, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    return Reference(path=path, sequences=sequences, file_sha256=digest)


# --------------------------------------------------------------------------- #
# 转录本
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Transcript:
    name: str
    gene: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]   # 零起点半开
    cds_start: int
    cds_end: int
    cds: str                             # 剪接后的 CDS（翻译方向）
    coords: tuple[int, ...]              # CDS 每个偏移对应的 1-based 基因组坐标
    pieces: tuple[tuple[int, int, int], ...]  # (CDS 偏移, 片段起0, 片段止0)
    span: tuple[int, int]                # 最小起点、最大终点（零起点半开）
    line_no: int

    @property
    def noncoding(self) -> bool:
        return self.cds_start == self.cds_end

    @property
    def cds_length(self) -> int:
        return len(self.cds)


def _parse_int_list(text: str, path: str, line_no: int, field: str) -> list[int]:
    parts = text.split(",")
    if not text or any(p.strip() == "" for p in parts):
        raise AnnotInputError(
            f"转录本表 {path} 第 {line_no} 行: {field} 为空或含空项: "
            f"{text!r}")
    values = []
    for p in parts:
        if not p.lstrip("-").isdigit() or not p.isascii():
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: {field} 含非整数项 {p!r}")
        values.append(int(p))
    return values


def _matching_exon(exons: tuple[tuple[int, int], ...],
                   base0: int) -> tuple[int, int] | None:
    """返回包含零起点碱基 base0 的外显子（按外显子区间判断）。"""
    for s, e in exons:
        if s <= base0 < e:
            return s, e
    return None


def _build_cds(seq: str, exons: tuple[tuple[int, int], ...], strand: str,
               cds_start: int, cds_end: int
               ) -> tuple[str, tuple[int, ...], tuple[tuple[int, int, int], ...]]:
    pieces: list[tuple[int, int]] = []
    for s, e in exons:
        cs = max(s, cds_start)
        ce = min(e, cds_end)
        if cs < ce:
            pieces.append((cs, ce))
    if not pieces:
        return "", (), ()
    ordered = pieces if strand == "+" else list(reversed(pieces))
    blocks: list[str] = []
    coords: list[int] = []
    piece_map: list[tuple[int, int, int]] = []
    offset = 0
    for s, e in ordered:
        block = seq[s:e]
        if strand == "-":
            block = reverse_complement(block)
        blocks.append(block)
        piece_map.append((offset, s, e))
        if strand == "+":
            coords.extend(range(s + 1, e + 1))
        else:
            coords.extend(range(e, s, -1))
        offset += e - s
    return "".join(blocks), tuple(coords), tuple(piece_map)


def load_transcripts(path: str, reference: Reference) -> tuple[Transcript, ...]:
    if not os.path.isfile(path):
        raise AnnotInputError(f"转录本表不存在: {path}")
    try:
        with open(path, "r", encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except OSError as exc:
        raise AnnotInputError(f"转录本表无法读取 {path}: {exc}") from exc
    if not rows:
        raise AnnotInputError(f"转录本表 {path} 为空")
    header = [c.strip() for c in rows[0]]
    if tuple(header) != TRANSCRIPT_HEADER:
        raise AnnotInputError(
            f"转录本表 {path} 第 1 行: 表头必须恰为 "
            f"{' / '.join(TRANSCRIPT_HEADER)}，实际为 {rows[0]!r}")

    out: list[Transcript] = []
    seen: dict[str, int] = {}
    for line_no, raw in enumerate(rows[1:], start=2):
        if not raw or all(c.strip() == "" for c in raw):
            continue
        if len(raw) != len(TRANSCRIPT_HEADER):
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: 必须恰为 "
                f"{len(TRANSCRIPT_HEADER)} 列，实际 {len(raw)} 列")
        (name, gene, chrom, strand, starts_s, ends_s,
         cds_start_s, cds_end_s) = (c.strip() for c in raw)

        if not name:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: transcript 不能为空")
        if name in seen:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: transcript {name!r} "
                f"重复（首次出现在第 {seen[name]} 行）")
        seen[name] = line_no
        if not gene:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: gene 不能为空")
        if not chrom:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: chrom 不能为空")
        if not reference.has(chrom):
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: 转录本 {name} 的 contig "
                f"{chrom!r} 不在 FASTA {reference.path} 中")
        if strand not in ("+", "-"):
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: strand 必须为 + 或 -，"
                f"实际 {strand!r}")

        starts = _parse_int_list(starts_s, path, line_no, "exon_starts")
        ends = _parse_int_list(ends_s, path, line_no, "exon_ends")
        if len(starts) != len(ends):
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: exon_starts 长度 "
                f"{len(starts)} 与 exon_ends 长度 {len(ends)} 不等")
        if not starts:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: 至少需要一个外显子")
        exons: list[tuple[int, int]] = []
        for i, (s, e) in enumerate(zip(starts, ends)):
            if s < 0 or e < 0:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 外显子坐标不得"
                    f"为负（第 {i + 1} 对 {s},{e}）")
            if s >= e:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 外显子 "
                    f"[{s},{e}) 为空（要求 start < end）")
            if e > reference.length(chrom):
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 外显子终点 {e} "
                    f"越出 contig {chrom} 长度 {reference.length(chrom)}")
            if exons and s < exons[-1][1]:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 外显子起点乱序或"
                    f"与前一外显子重叠（{s} < {exons[-1][1]}；外显子须"
                    f"按坐标递增且互不重叠）")
            if exons and s < exons[-1][0]:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 外显子起点未按"
                    f"坐标递增（{s} < {exons[-1][0]}）")
            exons.append((s, e))
        exons_t = tuple(exons)

        if not cds_start_s.lstrip("-").isdigit() or not cds_start_s.isascii() \
                or not cds_end_s.lstrip("-").isdigit() \
                or not cds_end_s.isascii():
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: cds_start/cds_end 必须"
                f"为非负整数，实际 {cds_start_s!r},{cds_end_s!r}")
        cds_start, cds_end = int(cds_start_s), int(cds_end_s)
        if cds_start < 0 or cds_end < 0:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: CDS 坐标不得为负")
        if cds_start > cds_end:
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: cds_start {cds_start} "
                f"大于 cds_end {cds_end}")
        if cds_end > reference.length(chrom):
            raise AnnotInputError(
                f"转录本表 {path} 第 {line_no} 行: cds_end {cds_end} "
                f"越出 contig {chrom} 长度 {reference.length(chrom)}")
        if cds_start != cds_end:
            # 非空编码范围：两端都必须落在外显子上（第一个与最后一个
            # 编码碱基，零起点分别为 cds_start 与 cds_end - 1）。
            if _matching_exon(exons_t, cds_start) is None:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 非法 CDS——首个"
                    f"编码碱基 {cds_start}（零起点）不在任何外显子上")
            if _matching_exon(exons_t, cds_end - 1) is None:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 非法 CDS——末个"
                    f"编码碱基 {cds_end - 1}（零起点）不在任何外显子上")

        seq = reference.sequences[chrom]
        cds, coords, piece_map = _build_cds(
            seq, exons_t, strand, cds_start, cds_end)
        if cds_start != cds_end:
            if not cds:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 非法 CDS——非空编码"
                    f"范围与外显子交集为空")
            if len(cds) % 3 != 0:
                raise AnnotInputError(
                    f"转录本表 {path} 第 {line_no} 行: 非法 CDS——拼接"
                    f"长度 {len(cds)} 不整除 3")
        span = (exons_t[0][0], exons_t[-1][1])
        out.append(Transcript(
            name=name, gene=gene, chrom=chrom, strand=strand,
            exons=exons_t, cds_start=cds_start, cds_end=cds_end,
            cds=cds, coords=coords, pieces=piece_map, span=span,
            line_no=line_no))
    if not out:
        raise AnnotInputError(f"转录本表 {path} 没有任何数据行")
    return tuple(out)


# --------------------------------------------------------------------------- #
# 变异
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Variant:
    id: str
    chrom: str
    pos: int                 # 1-based，REF 首碱基
    ref: str                 # 基因组正链
    alts: tuple[str, ...]
    line_no: int

    @property
    def start0(self) -> int:
        return self.pos - 1

    @property
    def end0(self) -> int:
        return self.pos - 1 + len(self.ref)


def load_variants(path: str, reference: Reference) -> tuple[Variant, ...]:
    if not os.path.isfile(path):
        raise AnnotInputError(f"变异表不存在: {path}")
    try:
        with open(path, "r", encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except OSError as exc:
        raise AnnotInputError(f"变异表无法读取 {path}: {exc}") from exc
    if not rows:
        raise AnnotInputError(f"变异表 {path} 为空")
    header = [c.strip() for c in rows[0]]
    if tuple(header) != VARIANT_HEADER:
        raise AnnotInputError(
            f"变异表 {path} 第 1 行: 表头必须恰为 "
            f"{' / '.join(VARIANT_HEADER)}，实际为 {rows[0]!r}")

    out: list[Variant] = []
    seen: dict[str, int] = {}
    for line_no, raw in enumerate(rows[1:], start=2):
        if not raw or all(c.strip() == "" for c in raw):
            continue
        if len(raw) != len(VARIANT_HEADER):
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: 必须恰为 "
                f"{len(VARIANT_HEADER)} 列，实际 {len(raw)} 列")
        (vid, chrom, pos_s, ref, alt_s) = (c.strip() for c in raw)
        ref = ref.upper()
        alt_s = alt_s.upper()
        if not vid:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: id 不能为空")
        if vid in seen:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: id {vid!r} 重复"
                f"（首次出现在第 {seen[vid]} 行）")
        seen[vid] = line_no
        if not chrom:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: chrom 不能为空")
        if not reference.has(chrom):
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: contig {chrom!r} "
                f"不在 FASTA {reference.path} 中")
        if not pos_s.isdigit() or not pos_s.isascii() or int(pos_s) < 1:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: pos {pos_s!r} 不是"
                f"从 1 开始的正整数")
        pos = int(pos_s)
        if not ref or set(ref) - _BASES:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: REF {ref!r} 必须是非空"
                f" A/C/G/T 序列（大写）")
        if not alt_s:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: ALT 不能为空")
        alts = alt_s.split(",")
        if any(not a or set(a) - _BASES for a in alts):
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: ALT {alt_s!r} 含空项或"
                f"非 A/C/G/T 字符")
        if len(set(alts)) != len(alts):
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: ALT {alt_s!r} 存在重复")
        if ref in alts:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: ALT 不得等于 REF {ref}")
        if pos - 1 + len(ref) > reference.length(chrom):
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: REF 区间 {chrom}:{pos}-"
                f"{pos + len(ref) - 1} 越出 contig 长度 "
                f"{reference.length(chrom)}")
        actual = reference.slice1(chrom, pos, pos + len(ref) - 1)
        if actual != ref:
            raise AnnotInputError(
                f"变异表 {path} 第 {line_no} 行: REF 校验失败——"
                f"{chrom}:{pos}-{pos + len(ref) - 1} 参考序列为 {actual!r}，"
                f"输入 REF 为 {ref!r}")
        out.append(Variant(id=vid, chrom=chrom, pos=pos, ref=ref,
                           alts=tuple(alts), line_no=line_no))
    if not out:
        raise AnnotInputError(f"变异表 {path} 没有任何数据行")
    return tuple(out)


# --------------------------------------------------------------------------- #
# 区域分类与编码后果
# --------------------------------------------------------------------------- #
def _intersects(txp: Transcript, chrom: str, v0: int, v1: int) -> bool:
    return (txp.chrom == chrom and v0 < txp.span[1]
            and txp.span[0] < v1)


def _in_single_intron(txp: Transcript, v0: int, v1: int) -> bool:
    for left, right in zip(txp.exons, txp.exons[1:]):
        gap_start, gap_end = left[1], right[0]
        if gap_start <= v0 and v1 <= gap_end:
            return True
    return False


def classify_region(txp: Transcript, v0: int, v1: int
                    ) -> tuple[str, str]:
    """返回 (区域类别, 原因)。v0/v1 为零起点半开 REF 区间。"""
    if not (txp.span[0] <= v0 and v1 <= txp.span[1]):
        return (REGION_BOUNDARY,
                "REF 区间越出该转录本的外显子跨度")
    exon = None
    for s, e in txp.exons:
        if s <= v0 and v1 <= e:
            exon = (s, e)
            break
    if exon is not None:
        if txp.noncoding:
            return (REGION_NONCODING_EXON,
                    "该异构体 CDS 为空（cds_start == cds_end），"
                    "外显子全部非编码")
        cs = max(exon[0], txp.cds_start)
        ce = min(exon[1], txp.cds_end)
        if ce <= cs:
            return (REGION_UTR,
                    "REF 完整落在该外显子的非编码片段（UTR），不翻译")
        if v1 <= cs:
            return (REGION_UTR,
                    "REF 完整落在编码区上游的非编码片段（UTR），不翻译")
        if v0 >= ce:
            return (REGION_UTR,
                    "REF 完整落在编码区下游的非编码片段（UTR），不翻译")
        if cs <= v0 and v1 <= ce:
            return REGION_CODING, ""
        return (REGION_BOUNDARY,
                "REF 横跨同一外显子的编码片段与非编码片段")
    if _in_single_intron(txp, v0, v1):
        return REGION_INTRON, "REF 完整落在单一内含子内"
    return (REGION_BOUNDARY,
            "REF 横跨外显子与内含子，或不在单一外显子/内含子内")


def _locate_piece(txp: Transcript, v0: int, v1: int
                  ) -> tuple[int, int, int] | None:
    for off, s, e in txp.pieces:
        if s <= v0 and v1 <= e:
            return off, s, e
    return None


def classify_consequence(aa_ref: str, aa_alt: str) -> str:
    if aa_ref == aa_alt:
        return CONSEQUENCE_SYNONYMOUS
    if aa_ref != "*" and aa_alt == "*":
        return CONSEQUENCE_STOP_GAINED
    if aa_ref == "*" and aa_alt != "*":
        return CONSEQUENCE_STOP_LOST
    return CONSEQUENCE_MISSENSE


@dataclass(frozen=True)
class Explanation:
    id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    alt_index: int
    transcript: str
    gene: str
    strand: str
    region: str
    cds_position: int | None
    codon_number: int | None
    codon_positions: tuple[int, int, int] | None
    codon_ref: str | None
    codon_alt: str | None
    aa_ref: str | None
    aa_alt: str | None
    consequence: str | None
    reason: str

    def as_dict(self) -> dict:
        return {
            "id": self.id, "chrom": self.chrom, "pos": self.pos,
            "ref": self.ref, "alt": self.alt, "alt_index": self.alt_index,
            "transcript": self.transcript, "gene": self.gene,
            "strand": self.strand, "region": self.region,
            "cds_position": self.cds_position,
            "codon_number": self.codon_number,
            "codon_pos_1": (None if self.codon_positions is None
                            else self.codon_positions[0]),
            "codon_pos_2": (None if self.codon_positions is None
                            else self.codon_positions[1]),
            "codon_pos_3": (None if self.codon_positions is None
                            else self.codon_positions[2]),
            "codon_ref": self.codon_ref, "codon_alt": self.codon_alt,
            "aa_ref": self.aa_ref, "aa_alt": self.aa_alt,
            "consequence": self.consequence, "reason": self.reason,
        }

    def as_row(self) -> list[str]:
        def s(v):
            return "" if v is None else str(v)
        return [
            self.id, self.chrom, str(self.pos), self.ref, self.alt,
            str(self.alt_index), self.transcript, self.gene, self.strand,
            self.region, s(self.cds_position), s(self.codon_number),
            s(None if self.codon_positions is None
              else self.codon_positions[0]),
            s(None if self.codon_positions is None
              else self.codon_positions[1]),
            s(None if self.codon_positions is None
              else self.codon_positions[2]),
            s(self.codon_ref), s(self.codon_alt), s(self.aa_ref),
            s(self.aa_alt), s(self.consequence), self.reason,
        ]


def _noncoding_explanations(variant: Variant, alt_index: int, alt: str,
                            txps: tuple[Transcript, ...], region: str,
                            reason: str) -> list[Explanation]:
    out = []
    for txp in txps:
        out.append(Explanation(
            id=variant.id, chrom=variant.chrom, pos=variant.pos,
            ref=variant.ref, alt=alt, alt_index=alt_index,
            transcript=txp.name, gene=txp.gene, strand=txp.strand,
            region=region, cds_position=None, codon_number=None,
            codon_positions=None, codon_ref=None, codon_alt=None,
            aa_ref=None, aa_alt=None, consequence=None, reason=reason))
    return out


def explain_coding(variant: Variant, alt_index: int, alt: str,
                   txp: Transcript) -> list[Explanation]:
    """等长替换的编码后果：整段替换一次性施加到 CDS 再按密码子解释。"""
    v0, v1 = variant.start0, variant.end0
    piece = _locate_piece(txp, v0, v1)
    if piece is None:
        raise AnnotInputError(
            f"内部错误: 变异 {variant.id} 的编码区间未映射到 CDS 片段")
    off, s, e = piece
    if txp.strand == "+":
        cds_off = off + (v0 - s)
        ref_cds = variant.ref
    else:
        # 负链：REF 末端对应片段反向互补后的起点；REF 需反向互补。
        cds_off = off + (e - v1)
        ref_cds = reverse_complement(variant.ref)
    alt_cds = alt if txp.strand == "+" else reverse_complement(alt)
    if len(ref_cds) != len(alt_cds):
        raise AnnotInputError(
            f"内部错误: {variant.id} 的等长替换在 CDS 方向长度不一致")
    if cds_off + len(ref_cds) > len(txp.cds):
        raise AnnotInputError(
            f"内部错误: {variant.id} 在转录本 {txp.name} 越过 CDS 末端")

    mutated = txp.cds[:cds_off] + alt_cds + txp.cds[cds_off + len(ref_cds):]
    first_codon = cds_off // 3 + 1
    last_codon = (cds_off + len(ref_cds) - 1) // 3 + 1

    out: list[Explanation] = []
    for codon_no in range(first_codon, last_codon + 1):
        start = (codon_no - 1) * 3
        codon_ref = txp.cds[start:start + 3]
        codon_alt = mutated[start:start + 3]
        aa_ref = translate_codon(codon_ref)
        aa_alt = translate_codon(codon_alt)
        positions = (txp.coords[start], txp.coords[start + 1],
                     txp.coords[start + 2])
        out.append(Explanation(
            id=variant.id, chrom=variant.chrom, pos=variant.pos,
            ref=variant.ref, alt=alt, alt_index=alt_index,
            transcript=txp.name, gene=txp.gene, strand=txp.strand,
            region=REGION_CODING,
            cds_position=cds_off + 1, codon_number=codon_no,
            codon_positions=positions, codon_ref=codon_ref,
            codon_alt=codon_alt, aa_ref=aa_ref, aa_alt=aa_alt,
            consequence=classify_consequence(aa_ref, aa_alt), reason=""))
    return out


def annotate_variant(variant: Variant, txps: tuple[Transcript, ...]
                     ) -> list[Explanation]:
    """一个输入记录的全部分项（ALT 列表顺序 → transcript 字典序）。"""
    # transcript 字典序：先按名称排序，输出顺序稳定。
    ordered = tuple(sorted(txps, key=lambda t: t.name))
    v0, v1 = variant.start0, variant.end0
    intersecting = tuple(
        t for t in ordered if _intersects(t, variant.chrom, v0, v1))

    # 输出顺序：ALT 列表顺序优先于 transcript 字典序。
    rows: list[Explanation] = []
    for alt_index, alt in enumerate(variant.alts, start=1):
        if len(alt) != len(variant.ref):
            # 长度不等的编辑优先给 unsupported_edit；保留相交 transcript。
            if intersecting:
                rows.extend(_noncoding_explanations(
                    variant, alt_index, alt, intersecting,
                    REGION_UNSUPPORTED_EDIT,
                    f"REF 长度 {len(variant.ref)} 与 ALT 长度 {len(alt)} "
                    f"不等；只支持等长碱基替换"))
            else:
                rows.append(Explanation(
                    id=variant.id, chrom=variant.chrom, pos=variant.pos,
                    ref=variant.ref, alt=alt, alt_index=alt_index,
                    transcript="", gene="", strand="",
                    region=REGION_UNSUPPORTED_EDIT, cds_position=None,
                    codon_number=None, codon_positions=None,
                    codon_ref=None, codon_alt=None, aa_ref=None,
                    aa_alt=None, consequence=None,
                    reason=f"REF 长度 {len(variant.ref)} 与 ALT 长度 "
                           f"{len(alt)} 不等；只支持等长碱基替换，且"
                           f"没有任何相交异构体"))
            continue

        if not intersecting:
            rows.append(Explanation(
                id=variant.id, chrom=variant.chrom, pos=variant.pos,
                ref=variant.ref, alt=alt, alt_index=alt_index,
                transcript="", gene="", strand="",
                region=REGION_INTERGENIC, cds_position=None,
                codon_number=None, codon_positions=None, codon_ref=None,
                codon_alt=None, aa_ref=None, aa_alt=None, consequence=None,
                reason="REF 区间不与任何异构体的外显子跨度相交"))
            continue

        for txp in intersecting:
            region, why = classify_region(txp, v0, v1)
            if region == REGION_CODING:
                rows.extend(explain_coding(variant, alt_index, alt, txp))
            else:
                rows.extend(_noncoding_explanations(
                    variant, alt_index, alt, (txp,), region, why))
    return rows


# --------------------------------------------------------------------------- #
# 输出（流式）
# --------------------------------------------------------------------------- #
def _prepare_out_dir(out_dir: str) -> None:
    if os.path.lexists(out_dir):
        if not os.path.isdir(out_dir):
            raise AnnotInputError(f"输出路径 {out_dir} 已存在且不是目录")
        if os.listdir(out_dir):
            raise AnnotInputError(
                f"输出目录 {out_dir} 非空：为避免混入不同批次，"
                f"请换一个新目录或先清理")
    else:
        os.makedirs(out_dir)


class AnnotationWriter:
    """逐变异写出并 flush；支持 TSV 与 JSONL。"""

    def __init__(self, out_dir: str, fmt: str):
        self.fmt = fmt
        self.path = os.path.join(
            out_dir, OUTPUT_TSV if fmt == "tsv" else OUTPUT_JSONL)
        self.fh = open(self.path, "w", encoding="utf-8", newline="")
        self.count = 0
        if fmt == "tsv":
            self._tsv = csv.writer(self.fh, delimiter="\t",
                                   lineterminator="\n")
            self._tsv.writerow(HEADER)
            self.fh.flush()

    def emit(self, explanations: list[Explanation]) -> None:
        for e in explanations:
            if self.fmt == "tsv":
                self._tsv.writerow(e.as_row())
            else:
                self.fh.write(json.dumps(e.as_dict(), ensure_ascii=False,
                                         sort_keys=False) + "\n")
            self.count += 1
        if explanations:
            self.fh.flush()

    def close(self) -> None:
        self.fh.close()


def _summary(reference: Reference, transcripts: tuple[Transcript, ...],
             variants: tuple[Variant, ...], counts: dict[str, int],
             consequence_counts: dict[str, int], fmt: str,
             out_path: str) -> dict:
    return {
        "status": "success",
        "format": fmt,
        "output": os.path.basename(out_path),
        "inputs": {
            "fasta": reference.path,
            "fasta_sha256": reference.file_sha256,
            "transcripts": len(transcripts),
            "variants": len(variants),
        },
        "transcripts": [
            {"transcript": t.name, "gene": t.gene, "chrom": t.chrom,
             "strand": t.strand, "exons": len(t.exons),
             "cds_length": t.cds_length,
             "noncoding": t.noncoding}
            for t in transcripts
        ],
        "explanations": sum(counts.values()),
        "region_counts": counts,
        "consequence_counts": consequence_counts,
    }


def run_annot(fasta: str, transcripts_path: str, variants_path: str,
              out_dir: str, fmt: str = "tsv") -> dict:
    """校验全部输入 → 逐条解释并流式写出 → 写 summary.json。

    输入错误在打开输出前即抛出 :class:`AnnotInputError`；输出 I/O 错误
    以 ``OSError`` 冒泡（CLI 映射为非零退出码），已写出的行保留供排查。
    """
    if fmt not in ("tsv", "jsonl"):
        raise AnnotInputError(f"输出格式必须为 tsv 或 jsonl，实际 {fmt!r}")
    reference = load_fasta(fasta)
    txps = load_transcripts(transcripts_path, reference)
    variants = load_variants(variants_path, reference)
    # 先完整校验再开写：合法输入不会半途因数据错误中止。
    _prepare_out_dir(out_dir)

    counts = {r: 0 for r in REGIONS}
    consequence_counts = {c: 0 for c in CONSEQUENCES}
    writer = AnnotationWriter(out_dir, fmt)
    try:
        for variant in variants:
            rows = annotate_variant(variant, txps)
            for row in rows:
                counts[row.region] = counts.get(row.region, 0) + 1
                if row.consequence:
                    consequence_counts[row.consequence] += 1
            writer.emit(rows)
    finally:
        writer.close()

    summary = _summary(reference, txps, variants, counts,
                       consequence_counts, fmt, writer.path)
    with open(os.path.join(out_dir, SUMMARY_JSON), "w",
              encoding="utf-8") as fh:
        fh.write(json.dumps(summary, indent=2, ensure_ascii=False,
                            sort_keys=True) + "\n")
    return summary
