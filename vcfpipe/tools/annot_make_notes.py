#!/usr/bin/env python3
"""生成“核对笔记”：先摆出外显子序列、剪接后的 CDS 和预期翻译，再运行
新入口并读回解释。

本脚本**不导入产品注释函数**：FASTA 解析、密码子表、CDS 剪接与翻译
在这里独立实现（只复用标准库），用于生成笔记里的“预期”部分；随后用
子进程调用产品入口，逐条核对输出与预期。不写第三份拼接实现——预期
部分的实现与 ``tools/annot_independent_check.py`` 同源，方便研究者
对照。

用法::

    python3 -I tools/annot_make_notes.py [--out output/annot_notes.md]
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(ROOT, "run_pipeline.py")
FA = os.path.join(ROOT, "data", "ref", "annot_syn.fa")
A = os.path.join(ROOT, "data", "annot")
TXP = os.path.join(A, "transcripts_demo.tsv")
VAR = os.path.join(A, "variants_demo.tsv")

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


_AAS = ("FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG")
TABLE = {}
for i, f in enumerate("TCAG"):
    for j, s in enumerate("TCAG"):
        for k, t in enumerate("TCAG"):
            TABLE[f + s + t] = _AAS[i * 16 + j * 4 + k]


def translate(seq: str) -> str:
    return "".join(TABLE[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3))


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line
    return seqs


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Tx:
    def __init__(self, row, seqs):
        self.name = row["transcript"]
        self.gene = row["gene"]
        self.chrom = row["chrom"]
        self.strand = row["strand"]
        self.starts = [int(x) for x in row["exon_starts"].split(",")]
        self.ends = [int(x) for x in row["exon_ends"].split(",")]
        self.cds_start = int(row["cds_start"])
        self.cds_end = int(row["cds_end"])
        self.seq = seqs[self.chrom]
        pieces = []
        self.genomic_pieces = []
        for a, b in zip(self.starts, self.ends):
            s, e = max(a, self.cds_start), min(b, self.cds_end)
            if s < e:
                pieces.append((s, e))
                self.genomic_pieces.append((s, e))
        order = pieces if self.strand == "+" else list(reversed(pieces))
        self.blocks = []
        self.coords = []
        for a, b in order:
            block = self.seq[a:b]
            if self.strand == "-":
                block = rc(block)
            self.blocks.append((a, b, block))
            if self.strand == "+":
                self.coords += list(range(a + 1, b + 1))
            else:
                self.coords += list(range(b, a, -1))
        self.cds = "".join(b for _, _, b in self.blocks)
        self.span = (self.starts[0], self.ends[-1])

    @property
    def protein(self):
        return translate(self.cds)



def run_product(out_dir):
    r = subprocess.run(
        [sys.executable, "-I", RUNNER, "annot", "--fasta", FA,
         "--transcripts", TXP, "--variants", VAR, "--out-dir", out_dir],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise SystemExit(r.stderr)
    with open(os.path.join(out_dir, "annotation.tsv"), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def seq_block(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "output",
                                                  "annot_notes.md"))
    ap.add_argument("--run-out-dir", default=os.path.join(
        ROOT, "output", "annot_notes_run"))
    args = ap.parse_args(argv)

    if os.path.isdir(args.run_out_dir):
        shutil.rmtree(args.run_out_dir)

    seqs = read_fasta(FA)
    txps = {t.name: t for t in (Tx(r, seqs) for r in read_tsv(TXP))}
    variants = read_tsv(VAR)
    rows = run_product(args.run_out_dir)

    lines = []
    lines.append("# 注释核对笔记（合成材料，本地计算）\n")
    lines.append("全部序列为合成材料；本笔记先用独立代码摆出外显子片段、"
                 "剪接后的 CDS 与预期翻译，再运行产品入口并读回解释表。\n")

    for name in ("TXP_POS1", "TXP_POS1_UTR", "TXP_NEG1"):
        t = txps[name]
        lines.append(f"## {name}（{t.gene}，{t.chrom}，{t.strand} 链）\n")
        lines.append(f"外显子（零起点半开）：`{list(zip(t.starts, t.ends))}`；"
                     f"编码范围 `[{t.cds_start},{t.cds_end})`。\n")
        lines.append("参与剪接的编码片段（基因组 1-based 闭区间，按拼接"
                     "顺序）：\n")
        for a, b, block in t.blocks:
            lines.append(f"- `{t.chrom}:{a + 1}-{b}`（{b - a} 碱基）"
                         f"→ `{block}`")
        lines.append("\n剪接后的 CDS：\n")
        lines.append("```\n" + seq_block(t.cds) + "\n```\n")
        prot = t.protein
        lines.append(f"预期翻译（{len(prot)} 个氨基酸）：\n")
        lines.append("```\n" + prot + "\n```\n")
        lines.append("密码子与基因组坐标（按翻译方向）：\n")
        lines.append("| # | 密码子 | 氨基酸 | 基因组坐标 |")
        lines.append("|---|--------|--------|------------|")
        for i in range(0, len(t.cds) - 2, 3):
            codon = t.cds[i:i + 3]
            cs = t.coords[i:i + 3]
            lines.append(f"| {i // 3 + 1} | {codon} | {TABLE[codon]} | "
                         f"{cs[0]}, {cs[1]}, {cs[2]} |")
        lines.append("")

    # 跨接缝密码子的专门说明
    t = txps["TXP_POS1"]
    lines.append("### 正链跨接缝：密码子 7\n")
    lines.append("编码片段 1 只有 20 碱基（10..30 零起点），最后的 `TAC` "
                 "首碱基 = chrP:29、第二碱基 = chrP:30（外显子 1 末碱基），"
                 "第三碱基 = chrP:91（外显子 2 首碱基）。拼接后从头切三"
                 "联体，阅读框贯穿接缝。\n")
    t = txps["TXP_NEG1"]
    lines.append("### 负链跨接缝：密码子 8\n")
    lines.append("负链片段按坐标递减拼接：先 `chrM:61-82` 的反向互补，"
                 "再接 `chrM:21-40` 的反向互补。密码子 8 三个坐标按翻译"
                 "方向为 61、40、39：首碱基来自高坐标片段，后两碱基来自"
                 "低坐标片段。\n")

    # 解释表
    lines.append("## 运行入口并读回解释\n")
    lines.append("```bash\npython3 -I run_pipeline.py annot \\\n"
                 "  --fasta data/ref/annot_syn.fa \\\n"
                 "  --transcripts data/annot/transcripts_demo.tsv \\\n"
                 "  --variants data/annot/variants_demo.tsv \\\n"
                 "  --out-dir output/annot_notes_run --format tsv\n```\n")
    lines.append("解释表（省略同名空列；空值记 `-`）：\n")
    lines.append("| id | ALT# | transcript | region | 密码子# | 三坐标 | "
                 "REF>ALT 密码子 | AA | 后果 |")
    lines.append("|----|------|-----------|--------|---------|--------|"
                 "----------------|----|------|")
    for r in rows:
        parts = [r[k] for k in ("codon_pos_1", "codon_pos_2", "codon_pos_3")]
        coords = ",".join(parts) if all(parts) else "-"
        codon_pair = ("-"
                      if not r["codon_ref"]
                      else f"{r['codon_ref']}>{r['codon_alt']}")
        aa_pair = ("-" if not r["aa_ref"]
                   else f"{r['aa_ref']}>{r['aa_alt']}")
        lines.append(
            f"| {r['id']} | {r['alt_index']} | {r['transcript'] or '-'} | "
            f"{r['region']} | {r['codon_number'] or '-'} | {coords} | "
            f"{codon_pair} | {aa_pair} | "
            f"{r['consequence'] or '-'} |")
    lines.append("")

    # 重点核对项
    lines.append("## 重点核对项\n")
    by = {(r["id"], r["alt_index"], r["transcript"]): r for r in rows}
    def get(vid, ai=1, tx="TXP_POS1"):
        return by[(vid, str(ai), tx)]

    items = [
        ("正链跨接缝密码子",
         get("V_CROSS_POS", 1), "TAC→TGC，Y→C，密码子 7 横跨 chrP:29/30/91"),
        ("负链跨接缝密码子、坐标方向",
         get("V_NEG_SEAM", 1, "TXP_NEG1"),
         "TGC→TAC，C→Y，坐标按翻译方向 61,40,39"),
        ("同位点在两种异构体中的不同后果",
         None, "chrP:95 T>A：TXP_POS1 为 coding（TGG→AGG，W→R）；"
               "同基因短 CDS 异构体 TXP_POS1_UTR 为 utr；"
               "TXP_NONCODING 为 noncoding_exon"),
        ("双碱基同时改变才得到氨基酸",
         get("V_DOUBLE_BASE", 1),
         "CGA→AGC：单独改第 1 碱基（AGA）仍为 R，单独改第 3 碱基（CGG）"
         "仍为 R；两处一起改才是 S"),
        ("跨两个密码子的替换",
         None, "chrP:112 ACG>TTA 同时触及密码子 14（AGA→AGT，R→S）与"
               "密码子 15（CGA→TAA，R→*）"),
        ("终止获得",
         get("V_STOP_GAINED", 1), "TGG→TGA，W→*"),
        ("终止丢失",
         get("V_STOP_LOST", 1), "TAA→CAA，*→Q"),
        ("同义（含 *→*）",
         get("V_SYNONYMOUS", 1), "GCT→GCC，A→A；规则上 TAA→TAG 也为 synonymous"),
        ("CDS 首碱基",
         get("V_CDS_START", 1), "ATG→GTG，密码子 1"),
        ("CDS 末碱基",
         get("V_CDS_END", 1), "TCG→TCA，密码子 50，同义"),
        ("外显子末端仍属编码",
         get("V_EXON_END_CODING", 1), "chrP:180 是外显子 2 末碱基，"
         "属密码子 37（179,180,211）"),
    ]
    for label, row, expect in items:
        if row is not None:
            got = (f"{row['region']} / {row['codon_number'] or '-'} / "
                   f"{row['consequence'] or '-'}")
        else:
            got = "见解释表"
        lines.append(f"- **{label}**：{expect}。输出：`{got}`")
    lines.append("")

    lines.append("## 区域类别预期\n")
    lines.append("| 变异 | transcript | 预期类别 | 说明 |")
    lines.append("|------|-----------|---------|------|")
    cats = [
        ("V_CDS_START", "TXP_POS1", "coding", "CDS 首碱基，两端都在外显子上"),
        ("V_CDS_END", "TXP_POS1", "coding", "CDS 末碱基"),
        ("V_EXON_END_CODING", "TXP_POS1", "coding", "外显子末碱基仍编码"),
        ("V_BOUNDARY_EXON", "TXP_POS1", "boundary_unsupported",
         "REF 从外显子末碱基伸入内含子"),
        ("V_BOUNDARY_CODING", "TXP_POS1", "boundary_unsupported",
         "REF 横跨非编码/编码交界"),
        ("V_BOUNDARY_SPAN", "TXP_POS1", "boundary_unsupported",
         "REF 伸到转录本跨度之外"),
        ("V_INTRON", "TXP_POS1", "intron", "单一内含子"),
        ("V_UTR5", "TXP_POS1", "utr", "5' 非编码片段"),
        ("V_UTR3", "TXP_POS1", "utr", "3' 非编码片段"),
        ("V_INTERGENIC", "（空）", "intergenic", "无相交异构体"),
        ("V_CODING_VS_NONCODING", "TXP_NONCODING", "noncoding_exon",
         "全非编码异构体"),
        ("V_UNSUPPORTED", "TXP_POS1", "unsupported_edit",
         "REF/ALT 长度不等，保留相交 transcript"),
    ]
    for vid, tx, region, why in cats:
        got = "—"
        match = [r for r in rows if r["id"] == vid
                 and (r["transcript"] == tx or (tx == "（空）" and not r["transcript"]))]
        if match:
            got = match[0]["region"]
        ok = "✓" if got == region else f"✗（实际 {got}）"
        lines.append(f"| {vid} | {tx} | {region} | {why} {ok} |")
    lines.append("")

    lines.append("## 独立核对\n")
    lines.append("```bash\npython3 -I tools/annot_independent_check.py \\\n"
                 "  --fasta data/ref/annot_syn.fa \\\n"
                 "  --transcripts data/annot/transcripts_demo.tsv \\\n"
                 "  --variants data/annot/variants_demo.tsv \\\n"
                 "  --out-dir output/annot_notes_run\n```\n")
    lines.append("核对器不导入产品函数：自己解析 FASTA/TSV、自己拼整条突变 "
                 "CDS、自己翻译，按同一排序逐字段比对。\n")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {args.out}")

    c = subprocess.run(
        [sys.executable, "-I",
         os.path.join(ROOT, "tools", "annot_independent_check.py"),
         "--fasta", FA, "--transcripts", TXP, "--variants", VAR,
         "--out-dir", args.run_out_dir],
        cwd=ROOT, capture_output=True, text=True, timeout=120)
    print(c.stdout.strip())
    return 0 if c.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
