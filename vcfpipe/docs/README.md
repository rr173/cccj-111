# vcfpipe：确定性 VCF 变异规范化与队列整合管线

把来自不同实验批次、不同 caller、不同 assembly 的小型 VCF，整合为统一
参考版本的**变异事实表**与**样本基因型矩阵**。仅依赖 Python 3.11 标准库。

## 1. 目录结构

```
vcfpipe/
├── src/vcfpipe/
│   ├── model.py       不可变数据记录与隔离原因码
│   ├── fasta.py       参考 FASTA 读取、碱基校验、文件/contig 摘要
│   ├── contigs.py     contig 字典与跨 assembly 别名解析
│   ├── vcf_io.py      VCF 头/记录解析与确定性 VCF 写出
│   ├── genotype.py    多等位 GT/AD/DP 拆分、剂量与缺失/相位语义
│   ├── normalize.py   最左对齐 + 最小表示（RuleBundle v1/v2）
│   ├── pipeline.py    端到端编排、合并、分区、清单与版本保护
│   ├── query.py       范围查询 / 样本切片 / source-to-output 追踪
│   └── cli.py         命令行入口
├── data/
│   ├── ref/           合成参考 syn37/syn38 与 contig 字典
│   └── input/         4 份合成 VCF、样本表、caller 优先级
├── tests/             手算单元测试 + 端到端集成测试（含独立解析器）
├── run_pipeline.py    本地 runner（-I 隔离模式下定位 src 包）
└── output/            运行后生成，按 {assembly}_{rule} 分数据集版本
```

## 2. 处理阶段

1. **头声明校验**：每份 VCF 必须显式给出 `##ASSEMBLY`、`##SAMPLE`、
   `##CALLER`、`##FILE_SUMMARY`，缺失即报错，不做文件名猜测。
2. **contig 别名**：经 contig 字典映射到目标 assembly；无法唯一映射
   （显式 `null` 或未声明）→ `contig_unmapped` 隔离。
3. **符号 ALT / 非法碱基**：`<DEL>`、`*`、非 `ACGTN` →
   `symbolic_alt` / `malformed_record` 隔离。
4. **参考碱基校验**：原始坐标上的 REF 必须与 FASTA 一致，否则
   `ref_base_mismatch` 隔离。
5. **GT 结构校验**：等位下标越界（如单 ALT 却出现 `0/2`）→
   `malformed_gt` 隔离。
6. **多等位拆分**：每个 ALT 生成一条单 ALT 记录，GT/AD/DP/相位按下述
   规则保留。
7. **规范化**：最左对齐 + 最小表示。
8. **合并**：稳定键 `(assembly, chrom, pos, REF, ALT, sample)` 分组，
   按显式 caller 优先级选胜出基因型，全部来源证据保留在 `evidence.tsv`。
9. **确定性写出**：分区事实表、长/宽样本矩阵、证据表、追踪表、隔离表、
   合并 VCF、清单。

## 3. 关键语义

### 3.1 多等位拆分（genotype.py）

对 ALT 下标 k（1-based）二值化 GT：

| 源等位 | 拆分后 |
|---|---|
| `0` 参考 | `0` |
| `k` 本 ALT | `1` |
| 其它非零 ALT | `.`（对本 ALT 不可判定，**不是** 0） |
| `.` 缺失 | `.` |

- 整体缺失（`./.`、`.|.`、无 GT）→ `gt_alleles=None`，**绝不**改写为
  参考纯合 `(0,0)`；DP 保留、AD 置空。
- 半缺失（`./1`）→ 保留 `(None, 1)`，剂量只计真实位。
- AD 只取 `[AD0, ADk]`，DP 整条保留，相位（`|` vs `/`）保留。
- 剂量守恒：源 GT 的每个非参考倍型位，恰好对应某条拆分 ALT 的 1 个剂量。

### 3.2 最左对齐（normalize.py，vt/bcftools 等价）

```
while 任一等位长度 > 1:
  if 两端末碱基相同:
      if pos==1: 边界策略（v1 隔离 / v2 保留并打 BOUNDARY）
      两等位同步左扩 base(pos-1)，各弹一个公共后缀；pos -= 1
  elif 两端首碱基相同且两侧长度都 >1:
      各裁一个公共前缀；pos += 1
  else: break
```

### 3.3 规则版本

| 版本 | 左移到 contig 起点仍需继续时 |
|---|---|
| `v1` | `left_align_boundary` 隔离 |
| `v2` | 保留部分左移结果，INFO 加 `BOUNDARY`，不隔离 |

数据集目录 = `{target_assembly}_{rule_version}`。换 assembly 或规则版本
得到新 dataset_id；同 dataset_id 重跑若清单字节不同则抛
`VersionConflict`，绝不原地覆盖。

## 4. 手算示例（syn37 chr1）

- 多等位 `chr1:10 C>A,T`，GT `1|2`，AD `30,5,4`：
  - ALT A → `1|.`，AD `30,5`；ALT T → `.|1`，AD `30,4`（相位保留）。
- poly-A 插入 `37 A>AA`：A 区段 29–38，逐次左移 9 次锚定到
  `28 G>GA`；同一插入的另一书写 `28 G>GA` 规范化到完全相同的键。
- poly-A 删除 `37 AA>A` → `28 GA>G`。
- poly-T 插入 `57 T>TT`：T 区段 49–58，前一位 48=C → `48 C>CT`。
- 边界 `pos1 A>CA`：末碱基相同但无 pos-1；v1 隔离，v2 保留并标记。
- 参考不符 `pos9 REF=T`（实际为 A）→ `ref_base_mismatch`。
- 别名 `1`→`chr1` 正常纳入；syn38 的 `chrZ` 目标缺失 → 隔离。
- 跨 caller 冲突 `chr1:5 A>G`：deepvar `0/1` vs freebayes `1/1`，
  优先级 deepvar > freebayes，胜出 `0/1`，两条证据均落 `evidence.tsv`。

## 5. 运行

```bash
# 运行（v1 / v2）
python3 -I run_pipeline.py run \
  --fasta data/ref/syn37.fa --contigs data/ref/contigs.json \
  --samples data/input/sample_sheet.csv \
  --priority data/input/caller_priority.json \
  --out-root output --rule-version v1

# 只读查询
python3 -I run_pipeline.py query --dataset output/syn37_v1 \
  --mode range --chrom chr1 --start 28 --end 48
python3 -I run_pipeline.py query --dataset output/syn37_v1 \
  --mode sample --sample S2
python3 -I run_pipeline.py query --dataset output/syn37_v1 \
  --mode trace-output --chrom chr1 --start 28 --ref G --alt GA
```

## 6. 测试

```bash
python3 -I -m unittest discover -s tests -v
```

34 个测试覆盖：手算多等位拆分、重复序列左对齐、参考碱基不符、contig
别名、半缺失 GT、相位与剂量映射、跨 caller 冲突、规则升级重建、守恒律，
以及一个**不导入本项目 VCF/TSV 读写代码**的独立解析器对生成 VCF 与
矩阵的交叉验证。
