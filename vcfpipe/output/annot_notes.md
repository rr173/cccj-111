# 注释核对笔记（合成材料，本地计算）

全部序列为合成材料；本笔记先用独立代码摆出外显子片段、剪接后的 CDS 与预期翻译，再运行产品入口并读回解释表。

## TXP_POS1（GENEPOS，chrP，+ 链）

外显子（零起点半开）：`[(0, 30), (90, 180), (210, 300)]`；编码范围 `[10,250)`。

参与剪接的编码片段（基因组 1-based 闭区间，按拼接顺序）：

- `chrP:11-30`（20 碱基）→ `ATGGCTTTCAAAGGACCCTA`
- `chrP:91-180`（90 碱基）→ `CGATTGGCTAAGCGGTCTCAGACGATAATGGGCTGACTTTGGCATCCTGTCCAACGTTCCTGAACGTACTGCATATTTGACGTGCCATGT`
- `chrP:211-250`（40 碱基）→ `ATCTAGGCCAGAGTGGATAGGGTTCAACGTTCCTGAATCG`

剪接后的 CDS：

```
ATGGCTTTCAAAGGACCCTACGATTGGCTAAGCGGTCTCAGACGATAATGGGCTGACTTT
GGCATCCTGTCCAACGTTCCTGAACGTACTGCATATTTGACGTGCCATGTATCTAGGCCA
GAGTGGATAGGGTTCAACGTTCCTGAATCG
```

预期翻译（50 个氨基酸）：

```
MAFKGPYDWLSGLRR*WADFGILSNVPERTAYLTCHVSRPEWIGFNVPES
```

密码子与基因组坐标（按翻译方向）：

| # | 密码子 | 氨基酸 | 基因组坐标 |
|---|--------|--------|------------|
| 1 | ATG | M | 11, 12, 13 |
| 2 | GCT | A | 14, 15, 16 |
| 3 | TTC | F | 17, 18, 19 |
| 4 | AAA | K | 20, 21, 22 |
| 5 | GGA | G | 23, 24, 25 |
| 6 | CCC | P | 26, 27, 28 |
| 7 | TAC | Y | 29, 30, 91 |
| 8 | GAT | D | 92, 93, 94 |
| 9 | TGG | W | 95, 96, 97 |
| 10 | CTA | L | 98, 99, 100 |
| 11 | AGC | S | 101, 102, 103 |
| 12 | GGT | G | 104, 105, 106 |
| 13 | CTC | L | 107, 108, 109 |
| 14 | AGA | R | 110, 111, 112 |
| 15 | CGA | R | 113, 114, 115 |
| 16 | TAA | * | 116, 117, 118 |
| 17 | TGG | W | 119, 120, 121 |
| 18 | GCT | A | 122, 123, 124 |
| 19 | GAC | D | 125, 126, 127 |
| 20 | TTT | F | 128, 129, 130 |
| 21 | GGC | G | 131, 132, 133 |
| 22 | ATC | I | 134, 135, 136 |
| 23 | CTG | L | 137, 138, 139 |
| 24 | TCC | S | 140, 141, 142 |
| 25 | AAC | N | 143, 144, 145 |
| 26 | GTT | V | 146, 147, 148 |
| 27 | CCT | P | 149, 150, 151 |
| 28 | GAA | E | 152, 153, 154 |
| 29 | CGT | R | 155, 156, 157 |
| 30 | ACT | T | 158, 159, 160 |
| 31 | GCA | A | 161, 162, 163 |
| 32 | TAT | Y | 164, 165, 166 |
| 33 | TTG | L | 167, 168, 169 |
| 34 | ACG | T | 170, 171, 172 |
| 35 | TGC | C | 173, 174, 175 |
| 36 | CAT | H | 176, 177, 178 |
| 37 | GTA | V | 179, 180, 211 |
| 38 | TCT | S | 212, 213, 214 |
| 39 | AGG | R | 215, 216, 217 |
| 40 | CCA | P | 218, 219, 220 |
| 41 | GAG | E | 221, 222, 223 |
| 42 | TGG | W | 224, 225, 226 |
| 43 | ATA | I | 227, 228, 229 |
| 44 | GGG | G | 230, 231, 232 |
| 45 | TTC | F | 233, 234, 235 |
| 46 | AAC | N | 236, 237, 238 |
| 47 | GTT | V | 239, 240, 241 |
| 48 | CCT | P | 242, 243, 244 |
| 49 | GAA | E | 245, 246, 247 |
| 50 | TCG | S | 248, 249, 250 |

## TXP_POS1_UTR（GENEPOS，chrP，+ 链）

外显子（零起点半开）：`[(0, 30), (90, 180), (210, 300)]`；编码范围 `[118,250)`。

参与剪接的编码片段（基因组 1-based 闭区间，按拼接顺序）：

- `chrP:119-180`（62 碱基）→ `TGGGCTGACTTTGGCATCCTGTCCAACGTTCCTGAACGTACTGCATATTTGACGTGCCATGT`
- `chrP:211-250`（40 碱基）→ `ATCTAGGCCAGAGTGGATAGGGTTCAACGTTCCTGAATCG`

剪接后的 CDS：

```
TGGGCTGACTTTGGCATCCTGTCCAACGTTCCTGAACGTACTGCATATTTGACGTGCCAT
GTATCTAGGCCAGAGTGGATAGGGTTCAACGTTCCTGAATCG
```

预期翻译（34 个氨基酸）：

```
WADFGILSNVPERTAYLTCHVSRPEWIGFNVPES
```

密码子与基因组坐标（按翻译方向）：

| # | 密码子 | 氨基酸 | 基因组坐标 |
|---|--------|--------|------------|
| 1 | TGG | W | 119, 120, 121 |
| 2 | GCT | A | 122, 123, 124 |
| 3 | GAC | D | 125, 126, 127 |
| 4 | TTT | F | 128, 129, 130 |
| 5 | GGC | G | 131, 132, 133 |
| 6 | ATC | I | 134, 135, 136 |
| 7 | CTG | L | 137, 138, 139 |
| 8 | TCC | S | 140, 141, 142 |
| 9 | AAC | N | 143, 144, 145 |
| 10 | GTT | V | 146, 147, 148 |
| 11 | CCT | P | 149, 150, 151 |
| 12 | GAA | E | 152, 153, 154 |
| 13 | CGT | R | 155, 156, 157 |
| 14 | ACT | T | 158, 159, 160 |
| 15 | GCA | A | 161, 162, 163 |
| 16 | TAT | Y | 164, 165, 166 |
| 17 | TTG | L | 167, 168, 169 |
| 18 | ACG | T | 170, 171, 172 |
| 19 | TGC | C | 173, 174, 175 |
| 20 | CAT | H | 176, 177, 178 |
| 21 | GTA | V | 179, 180, 211 |
| 22 | TCT | S | 212, 213, 214 |
| 23 | AGG | R | 215, 216, 217 |
| 24 | CCA | P | 218, 219, 220 |
| 25 | GAG | E | 221, 222, 223 |
| 26 | TGG | W | 224, 225, 226 |
| 27 | ATA | I | 227, 228, 229 |
| 28 | GGG | G | 230, 231, 232 |
| 29 | TTC | F | 233, 234, 235 |
| 30 | AAC | N | 236, 237, 238 |
| 31 | GTT | V | 239, 240, 241 |
| 32 | CCT | P | 242, 243, 244 |
| 33 | GAA | E | 245, 246, 247 |
| 34 | TCG | S | 248, 249, 250 |

## TXP_NEG1（GENENEG，chrM，- 链）

外显子（零起点半开）：`[(0, 40), (60, 104)]`；编码范围 `[20,82)`。

参与剪接的编码片段（基因组 1-based 闭区间，按拼接顺序）：

- `chrM:61-82`（22 碱基）→ `ATGGCTTCGTGGTACGACGGTT`
- `chrM:21-40`（20 碱基）→ `GCGAATTCAAAGGTTGAGCT`

剪接后的 CDS：

```
ATGGCTTCGTGGTACGACGGTTGCGAATTCAAAGGTTGAGCT
```

预期翻译（14 个氨基酸）：

```
MASWYDGCEFKG*A
```

密码子与基因组坐标（按翻译方向）：

| # | 密码子 | 氨基酸 | 基因组坐标 |
|---|--------|--------|------------|
| 1 | ATG | M | 82, 81, 80 |
| 2 | GCT | A | 79, 78, 77 |
| 3 | TCG | S | 76, 75, 74 |
| 4 | TGG | W | 73, 72, 71 |
| 5 | TAC | Y | 70, 69, 68 |
| 6 | GAC | D | 67, 66, 65 |
| 7 | GGT | G | 64, 63, 62 |
| 8 | TGC | C | 61, 40, 39 |
| 9 | GAA | E | 38, 37, 36 |
| 10 | TTC | F | 35, 34, 33 |
| 11 | AAA | K | 32, 31, 30 |
| 12 | GGT | G | 29, 28, 27 |
| 13 | TGA | * | 26, 25, 24 |
| 14 | GCT | A | 23, 22, 21 |

### 正链跨接缝：密码子 7

编码片段 1 只有 20 碱基（10..30 零起点），最后的 `TAC` 首碱基 = chrP:29、第二碱基 = chrP:30（外显子 1 末碱基），第三碱基 = chrP:91（外显子 2 首碱基）。拼接后从头切三联体，阅读框贯穿接缝。

### 负链跨接缝：密码子 8

负链片段按坐标递减拼接：先 `chrM:61-82` 的反向互补，再接 `chrM:21-40` 的反向互补。密码子 8 三个坐标按翻译方向为 61、40、39：首碱基来自高坐标片段，后两碱基来自低坐标片段。

## 运行入口并读回解释

```bash
python3 -I run_pipeline.py annot \
  --fasta data/ref/annot_syn.fa \
  --transcripts data/annot/transcripts_demo.tsv \
  --variants data/annot/variants_demo.tsv \
  --out-dir output/annot_notes_run --format tsv
```

解释表（省略同名空列；空值记 `-`）：

| id | ALT# | transcript | region | 密码子# | 三坐标 | REF>ALT 密码子 | AA | 后果 |
|----|------|-----------|--------|---------|--------|----------------|----|------|
| V_CROSS_POS | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CROSS_POS | 1 | TXP_POS1 | coding | 7 | 29,30,91 | TAC>TGC | Y>C | missense |
| V_CROSS_POS | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_CROSS_POS | 2 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CROSS_POS | 2 | TXP_POS1 | coding | 7 | 29,30,91 | TAC>TTC | Y>F | missense |
| V_CROSS_POS | 2 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_CROSS_POS | 3 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CROSS_POS | 3 | TXP_POS1 | coding | 7 | 29,30,91 | TAC>TCC | Y>S | missense |
| V_CROSS_POS | 3 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_NEG_SEAM | 1 | TXP_NEG1 | coding | 8 | 61,40,39 | TGC>TAC | C>Y | missense |
| V_STOP_LOST | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_STOP_LOST | 1 | TXP_POS1 | coding | 16 | 116,117,118 | TAA>CAA | *>Q | stop_lost |
| V_STOP_LOST | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_STOP_GAINED | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_STOP_GAINED | 1 | TXP_POS1 | coding | 17 | 119,120,121 | TGG>TGA | W>* | stop_gained |
| V_STOP_GAINED | 1 | TXP_POS1_UTR | coding | 1 | 119,120,121 | TGG>TGA | W>* | stop_gained |
| V_SYNONYMOUS | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_SYNONYMOUS | 1 | TXP_POS1 | coding | 2 | 14,15,16 | GCT>GCC | A>A | synonymous |
| V_SYNONYMOUS | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_DOUBLE_BASE | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_DOUBLE_BASE | 1 | TXP_POS1 | coding | 15 | 113,114,115 | CGA>AGC | R>S | missense |
| V_DOUBLE_BASE | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_TWO_CODONS | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_TWO_CODONS | 1 | TXP_POS1 | coding | 14 | 110,111,112 | AGA>AGT | R>S | missense |
| V_TWO_CODONS | 1 | TXP_POS1 | coding | 15 | 113,114,115 | CGA>TAA | R>* | stop_gained |
| V_TWO_CODONS | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_CODING_VS_NONCODING | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CODING_VS_NONCODING | 1 | TXP_POS1 | coding | 9 | 95,96,97 | TGG>AGG | W>R | missense |
| V_CODING_VS_NONCODING | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_UTR5 | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_UTR5 | 1 | TXP_POS1 | utr | - | - | - | - | - |
| V_UTR5 | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_UTR3 | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_UTR3 | 1 | TXP_POS1 | utr | - | - | - | - | - |
| V_UTR3 | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_INTRON | 1 | TXP_NONCODING | intron | - | - | - | - | - |
| V_INTRON | 1 | TXP_POS1 | intron | - | - | - | - | - |
| V_INTRON | 1 | TXP_POS1_UTR | intron | - | - | - | - | - |
| V_CDS_START | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CDS_START | 1 | TXP_POS1 | coding | 1 | 11,12,13 | ATG>GTG | M>V | missense |
| V_CDS_START | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_CDS_END | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_CDS_END | 1 | TXP_POS1 | coding | 50 | 248,249,250 | TCG>TCA | S>S | synonymous |
| V_CDS_END | 1 | TXP_POS1_UTR | coding | 34 | 248,249,250 | TCG>TCA | S>S | synonymous |
| V_EXON_END_CODING | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_EXON_END_CODING | 1 | TXP_POS1 | coding | 37 | 179,180,211 | GTA>GCA | V>A | missense |
| V_EXON_END_CODING | 1 | TXP_POS1_UTR | coding | 21 | 179,180,211 | GTA>GCA | V>A | missense |
| V_BOUNDARY_EXON | 1 | TXP_NONCODING | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_EXON | 1 | TXP_POS1 | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_EXON | 1 | TXP_POS1_UTR | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_CODING | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_BOUNDARY_CODING | 1 | TXP_POS1 | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_CODING | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_BOUNDARY_SPAN | 1 | TXP_NONCODING | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_SPAN | 1 | TXP_POS1 | boundary_unsupported | - | - | - | - | - |
| V_BOUNDARY_SPAN | 1 | TXP_POS1_UTR | boundary_unsupported | - | - | - | - | - |
| V_INTERGENIC | 1 | - | intergenic | - | - | - | - | - |
| V_UNSUPPORTED | 1 | TXP_NONCODING | unsupported_edit | - | - | - | - | - |
| V_UNSUPPORTED | 1 | TXP_POS1 | unsupported_edit | - | - | - | - | - |
| V_UNSUPPORTED | 1 | TXP_POS1_UTR | unsupported_edit | - | - | - | - | - |
| V_ALT_MIXED | 1 | TXP_NONCODING | noncoding_exon | - | - | - | - | - |
| V_ALT_MIXED | 1 | TXP_POS1 | coding | 9 | 95,96,97 | TGG>TAG | W>* | stop_gained |
| V_ALT_MIXED | 1 | TXP_POS1_UTR | utr | - | - | - | - | - |
| V_ALT_MIXED | 2 | TXP_NONCODING | unsupported_edit | - | - | - | - | - |
| V_ALT_MIXED | 2 | TXP_POS1 | unsupported_edit | - | - | - | - | - |
| V_ALT_MIXED | 2 | TXP_POS1_UTR | unsupported_edit | - | - | - | - | - |
| V_NEG_INTRON | 1 | TXP_NEG1 | intron | - | - | - | - | - |
| V_NEG_UTR | 1 | TXP_NEG1 | utr | - | - | - | - | - |

## 重点核对项

- **正链跨接缝密码子**：TAC→TGC，Y→C，密码子 7 横跨 chrP:29/30/91。输出：`coding / 7 / missense`
- **负链跨接缝密码子、坐标方向**：TGC→TAC，C→Y，坐标按翻译方向 61,40,39。输出：`coding / 8 / missense`
- **同位点在两种异构体中的不同后果**：chrP:95 T>A：TXP_POS1 为 coding（TGG→AGG，W→R）；同基因短 CDS 异构体 TXP_POS1_UTR 为 utr；TXP_NONCODING 为 noncoding_exon。输出：`见解释表`
- **双碱基同时改变才得到氨基酸**：CGA→AGC：单独改第 1 碱基（AGA）仍为 R，单独改第 3 碱基（CGG）仍为 R；两处一起改才是 S。输出：`coding / 15 / missense`
- **跨两个密码子的替换**：chrP:112 ACG>TTA 同时触及密码子 14（AGA→AGT，R→S）与密码子 15（CGA→TAA，R→*）。输出：`见解释表`
- **终止获得**：TGG→TGA，W→*。输出：`coding / 17 / stop_gained`
- **终止丢失**：TAA→CAA，*→Q。输出：`coding / 16 / stop_lost`
- **同义（含 *→*）**：GCT→GCC，A→A；规则上 TAA→TAG 也为 synonymous。输出：`coding / 2 / synonymous`
- **CDS 首碱基**：ATG→GTG，密码子 1。输出：`coding / 1 / missense`
- **CDS 末碱基**：TCG→TCA，密码子 50，同义。输出：`coding / 50 / synonymous`
- **外显子末端仍属编码**：chrP:180 是外显子 2 末碱基，属密码子 37（179,180,211）。输出：`coding / 37 / missense`

## 区域类别预期

| 变异 | transcript | 预期类别 | 说明 |
|------|-----------|---------|------|
| V_CDS_START | TXP_POS1 | coding | CDS 首碱基，两端都在外显子上 ✓ |
| V_CDS_END | TXP_POS1 | coding | CDS 末碱基 ✓ |
| V_EXON_END_CODING | TXP_POS1 | coding | 外显子末碱基仍编码 ✓ |
| V_BOUNDARY_EXON | TXP_POS1 | boundary_unsupported | REF 从外显子末碱基伸入内含子 ✓ |
| V_BOUNDARY_CODING | TXP_POS1 | boundary_unsupported | REF 横跨非编码/编码交界 ✓ |
| V_BOUNDARY_SPAN | TXP_POS1 | boundary_unsupported | REF 伸到转录本跨度之外 ✓ |
| V_INTRON | TXP_POS1 | intron | 单一内含子 ✓ |
| V_UTR5 | TXP_POS1 | utr | 5' 非编码片段 ✓ |
| V_UTR3 | TXP_POS1 | utr | 3' 非编码片段 ✓ |
| V_INTERGENIC | （空） | intergenic | 无相交异构体 ✓ |
| V_CODING_VS_NONCODING | TXP_NONCODING | noncoding_exon | 全非编码异构体 ✓ |
| V_UNSUPPORTED | TXP_POS1 | unsupported_edit | REF/ALT 长度不等，保留相交 transcript ✓ |

## 独立核对

```bash
python3 -I tools/annot_independent_check.py \
  --fasta data/ref/annot_syn.fa \
  --transcripts data/annot/transcripts_demo.tsv \
  --variants data/annot/variants_demo.tsv \
  --out-dir output/annot_notes_run
```

核对器不导入产品函数：自己解析 FASTA/TSV、自己拼整条突变 CDS、自己翻译，按同一排序逐字段比对。

