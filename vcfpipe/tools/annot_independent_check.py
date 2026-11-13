#!/usr/bin/env python3
"""独立核对器：直接从 FASTA 拼接整条突变 CDS 再翻译，与产品输出逐行比对。

**不导入** ``vcfpipe.annot`` 的任何代码：FASTA 解析、转录本/变异 TSV
解析、CDS 剪接、坐标换算、遗传密码表与分类规则都在这里独立重写一遍。
对每个输入变异、每个 ALT、每个命名转录本，独立推导：

* 区域类别（含边界/内含子/基因间/长度不等）；
* 编码项的 CDS 位置、密码子编号、三个基因组坐标、替换前后密码子与
  氨基酸；
* 后果类别。

并把结果与 ``<out-dir>/annotation.tsv``（或 ``.jsonl``）逐行、按同一
排序键比对；任何字段不符即打印差异并返回非零。

用法::

    python3 -I tools/annot_independent_check.py \
        --fasta data/ref/annot_syn.fa \
        --transcripts data/annot/transcripts_demo.tsv \
        --variants data/annot/variants_demo.tsv \
        --out-dir output/annot_demo

仅用于测试核对，不参与生产输出。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASES = set("ACGT")


def revcomp(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


# 独立遗传密码表：由 NCBI table 1 的碱基顺序展开，不引用产品字典。
_AAS = ("FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG")


def _codon_table() -> dict[str, str]:
    table, k = {}, 0
    for first in "TCAG":
        for second in "TCAG":
            for third in "TCAG":
                table[first + second + third] = _AAS[k]
                k += 1
    return table


TABLE = _codon_table()
assert len(TABLE) == 64


def translate(seq: str) -> str:
    if len(seq) % 3:
        raise ValueError("CDS 长度不整除 3")
    return "".join(TABLE[seq[i:i + 3]] for i in range(0, len(seq), 3))


# --------------------------------------------------------------------------- #
# 独立解析
# --------------------------------------------------------------------------- #
def read_fasta(path: str) -> dict[str, str]:
    seqs, name = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            else:
                seqs[name] += line
    return seqs


def read_tsv(path: str) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Txp:
    def __init__(self, row: dict, seqs: dict[str, str]):
        self.name = row["transcript"]
        self.gene = row["gene"]
        self.chrom = row["chrom"]
        self.strand = row["strand"]
        starts = [int(x) for x in row["exon_starts"].split(",")]
        ends = [int(x) for x in row["exon_ends"].split(",")]
        self.exons = list(zip(starts, ends))
        self.cds_start = int(row["cds_start"])
        self.cds_end = int(row["cds_end"])
        self.seq = seqs[self.chrom]
        self._build()

    def _build(self):
        pieces = []
        for s, e in self.exons:
            a, b = max(s, self.cds_start), min(e, self.cds_end)
            if a < b:
                pieces.append((a, b))
        order = pieces if self.strand == "+" else list(reversed(pieces))
        blocks = []
        self.coord_of_offset = []
        for a, b in order:
            block = self.seq[a:b]
            if self.strand == "-":
                block = revcomp(block)
            blocks.append(block)
            if self.strand == "+":
                self.coord_of_offset.extend(range(a + 1, b + 1))
            else:
                self.coord_of_offset.extend(range(b, a, -1))
        self.cds = "".join(blocks)
        self.pieces = []
        off = 0
        for a, b in order:
            self.pieces.append((off, a, b))
            off += b - a
        self.span = (self.exons[0][0], self.exons[-1][1])

    @property
    def noncoding(self) -> bool:
        return self.cds_start == self.cds_end


def region_of(t: Txp, start0: int, end0: int) -> tuple[str, str]:
    if not (t.span[0] <= start0 and end0 <= t.span[1]):
        return "boundary_unsupported", "越出跨度"
    for s, e in t.exons:
        if s <= start0 and end0 <= e:
            if t.noncoding:
                return "noncoding_exon", "全非编码"
            cs, ce = max(s, t.cds_start), min(e, t.cds_end)
            if ce <= cs or end0 <= cs or start0 >= ce:
                return "utr", "非编码片段"
            if cs <= start0 and end0 <= ce:
                return "coding", ""
            return "boundary_unsupported", "编码/非编码交界"
    for (_, gap0), (gap1, _) in zip(t.exons, t.exons[1:]):
        if gap0 <= start0 and end0 <= gap1:
            return "intron", "单一内含子"
    return "boundary_unsupported", "跨外显子/内含子"


def class_of(ref_aa: str, alt_aa: str) -> str:
    if ref_aa == alt_aa:
        return "synonymous"
    if ref_aa != "*" and alt_aa == "*":
        return "stop_gained"
    if ref_aa == "*" and alt_aa != "*":
        return "stop_lost"
    return "missense"


def coding_rows(variant: dict, alt: str, alt_index: int,
                t: Txp) -> list[dict]:
    start0 = variant["pos"] - 1
    end0 = start0 + len(variant["ref"])
    off = None
    for p_off, a, b in t.pieces:
        if a <= start0 and end0 <= b:
            off = p_off
            break
    if off is None:
        raise AssertionError("编码变异未落在单一编码片段")
    if t.strand == "+":
        cds_off = off + (start0 - a)
        ref_cds, alt_cds = variant["ref"], alt
    else:
        cds_off = off + (b - end0)
        ref_cds, alt_cds = revcomp(variant["ref"]), revcomp(alt)
    assert cds_off >= 0 and len(ref_cds) == len(alt_cds)
    mutated = t.cds[:cds_off] + alt_cds + t.cds[cds_off + len(ref_cds):]
    out = []
    cds_position = cds_off + 1          # 首个被改动碱基的 CDS 位置（1 起）
    for codon_no in range(cds_off // 3 + 1,
                          (cds_off + len(ref_cds) - 1) // 3 + 2):
        s = (codon_no - 1) * 3
        c_ref, c_alt = t.cds[s:s + 3], mutated[s:s + 3]
        out.append({
            "id": variant["id"], "chrom": variant["chrom"],
            "pos": variant["pos"], "ref": variant["ref"], "alt": alt,
            "alt_index": alt_index, "transcript": t.name, "gene": t.gene,
            "strand": t.strand, "region": "coding",
            "cds_position": cds_position, "codon_number": codon_no,
            "codon_pos_1": t.coord_of_offset[s],
            "codon_pos_2": t.coord_of_offset[s + 1],
            "codon_pos_3": t.coord_of_offset[s + 2],
            "codon_ref": c_ref, "codon_alt": c_alt,
            "aa_ref": TABLE[c_ref], "aa_alt": TABLE[c_alt],
            "consequence": class_of(TABLE[c_ref], TABLE[c_alt]),
            "reason": "",
        })
    return out


def expected_rows(fasta: str, transcripts: str, variants: str) -> list[dict]:
    seqs = read_fasta(fasta)
    txps = sorted((Txp(r, seqs) for r in read_tsv(transcripts)),
                  key=lambda t: t.name)
    rows: list[dict] = []
    for vr in read_tsv(variants):
        variant = {
            "id": vr["id"], "chrom": vr["chrom"], "pos": int(vr["pos"]),
            "ref": vr["ref"],
        }
        start0 = variant["pos"] - 1
        end0 = start0 + len(variant["ref"])
        hits = [t for t in txps
                if t.chrom == variant["chrom"]
                and start0 < t.span[1] and t.span[0] < end0]
        for alt_index, alt in enumerate(vr["alt"].split(","), start=1):
            if len(alt) != len(variant["ref"]):
                if not hits:
                    rows.append({
                        "id": variant["id"], "chrom": variant["chrom"],
                        "pos": variant["pos"], "ref": variant["ref"],
                        "alt": alt, "alt_index": alt_index,
                        "transcript": "", "gene": "", "strand": "",
                        "region": "unsupported_edit",
                        "cds_position": None, "codon_number": None,
                        "codon_pos_1": None, "codon_pos_2": None,
                        "codon_pos_3": None, "codon_ref": None,
                        "codon_alt": None, "aa_ref": None, "aa_alt": None,
                        "consequence": None, "reason": "长度不等且无相交",
                    })
                for t in hits:
                    rows.append({
                        "id": variant["id"], "chrom": variant["chrom"],
                        "pos": variant["pos"], "ref": variant["ref"],
                        "alt": alt, "alt_index": alt_index,
                        "transcript": t.name, "gene": t.gene,
                        "strand": t.strand, "region": "unsupported_edit",
                        "cds_position": None, "codon_number": None,
                        "codon_pos_1": None, "codon_pos_2": None,
                        "codon_pos_3": None, "codon_ref": None,
                        "codon_alt": None, "aa_ref": None, "aa_alt": None,
                        "consequence": None, "reason": "长度不等",
                    })
                continue
            if not hits:
                rows.append({
                    "id": variant["id"], "chrom": variant["chrom"],
                    "pos": variant["pos"], "ref": variant["ref"],
                    "alt": alt, "alt_index": alt_index,
                    "transcript": "", "gene": "", "strand": "",
                    "region": "intergenic", "cds_position": None,
                    "codon_number": None, "codon_pos_1": None,
                    "codon_pos_2": None, "codon_pos_3": None,
                    "codon_ref": None, "codon_alt": None, "aa_ref": None,
                    "aa_alt": None, "consequence": None,
                    "reason": "无相交异构体",
                })
                continue
            for t in hits:
                region, why = region_of(t, start0, end0)
                if region == "coding":
                    rows.extend(coding_rows(variant, alt, alt_index, t))
                else:
                    rows.append({
                        "id": variant["id"], "chrom": variant["chrom"],
                        "pos": variant["pos"], "ref": variant["ref"],
                        "alt": alt, "alt_index": alt_index,
                        "transcript": t.name, "gene": t.gene,
                        "strand": t.strand, "region": region,
                        "cds_position": None, "codon_number": None,
                        "codon_pos_1": None, "codon_pos_2": None,
                        "codon_pos_3": None, "codon_ref": None,
                        "codon_alt": None, "aa_ref": None, "aa_alt": None,
                        "consequence": None, "reason": why,
                    })
    return rows


def _norm(v):
    if isinstance(v, int):
        return str(v)
    if v is None:
        return ""
    return v


def _key(row: dict) -> tuple:
    return (row["id"], str(row["alt_index"]), row["transcript"] or "-",
            row["codon_number"])


def load_actual(out_dir: str) -> list[dict]:
    tsv = os.path.join(out_dir, "annotation.tsv")
    jsonl = os.path.join(out_dir, "annotation.jsonl")
    if os.path.isfile(tsv):
        return read_tsv(tsv)
    rows = []
    with open(jsonl, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--variants", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args(argv)

    expected = expected_rows(args.fasta, args.transcripts, args.variants)
    actual = load_actual(args.out_dir)
    problems = []
    if len(expected) != len(actual):
        problems.append(f"行数不符：期望 {len(expected)}，产品 {len(actual)}")

    keys = ["id", "chrom", "pos", "ref", "alt", "alt_index", "transcript",
            "gene", "strand", "region", "cds_position", "codon_number",
            "codon_pos_1", "codon_pos_2", "codon_pos_3", "codon_ref",
            "codon_alt", "aa_ref", "aa_alt", "consequence"]
    n = min(len(expected), len(actual))
    for i in range(n):
        e, a = expected[i], actual[i]
        diffs = []
        for k in keys:
            if _norm(e.get(k)) != _norm(a.get(k)):
                diffs.append(f"{k}: 期望 {e.get(k)!r} 产品 {a.get(k)!r}")
        if diffs:
            problems.append(f"第 {i + 1} 行（{_key(e)}）: " + "; ".join(diffs))
    if problems:
        print("不一致：", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print(f"核对通过：{len(expected)} 行逐字段一致（独立拼接整条突变 CDS）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
