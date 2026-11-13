# vcfpipe：确定性 VCF 变异规范化与队列整合管线

把来自不同实验批次、不同 caller、不同 assembly 的小型 VCF，整合为统一
参考版本的**变异事实表**与**样本基因型矩阵**；上游带一个**测序数据接收
环节**（`demux`）：校验双端 FASTQ、按条码把混合文库拆回各样本并提取
UMI，拆出的每样本 read 对即可交给外部比对/变异调用。仅依赖
Python 3.11 标准库。

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
│   ├── assoc.py       基因级携带者关联分析（Fisher 精确检验 + BH 校正）
│   ├── demux.py       测序数据接收：FASTQ 校验、按条码拆分、UMI 提取
│   └── cli.py         命令行入口
├── data/
│   ├── ref/           合成参考 syn37/syn38 与 contig 字典
│   ├── input/         4 份合成 VCF、样本表、caller 优先级
│   ├── assoc/         关联分析演示输入（标签表、基因/可判读区间 BED）
│   └── demux/         接收演示输入（合成序列）：10 对手算样例 +
│                      alternating/ 去向交替样例 + 确定性生成脚本
├── benchmarks/        run_benchmark.py：十万对级生成输入的规模/内存演示
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
得到新 dataset_id。一个已发布的 dataset_id 永久对应同一份完整内容：

- **原子发布**：全部输出（含 manifest.json）先写入 `output/` 下的隐藏
  staging 目录 `.{dataset_id}.staging-*`，再由**单个 `os.rename`** 提交
  为最终目录。查询侧任意时刻只能看到旧完整版本或新完整版本，读到
  “部分新分区 + 旧清单”的中间态不存在。
- **重跑语义**：相同输入重跑得到逐字节相同的清单 → 幂等成功，已发布
  目录不被触碰；任何实质差异 → 抛 `VersionConflict`（退出码 2），
  拒绝前后旧目录的每个文件与查询结果完全一致，绝不原地覆盖。
- **中断处理**：写到一半失败或进程被终止时，垃圾只留在隐藏 staging
  目录里，查询链不可见；首次构建失败不会留下看似可查询的半成品。
  下一次构建在持锁状态下清理本数据集的 staging 残留（此时它们只
  可能属于已退出的进程）。
- **并发互斥**：同一 dataset_id 的构建经 `flock`（锁文件
  `output/.{dataset_id}.lock`）串行化。两个进程同时写同一版本时，
  只有一个完整结果成为发布版本：内容相同的另一个进程确认已存在而
  成功，内容不同的另一个进程失败且不影响胜出结果。锁随 fd 关闭
  释放，进程退出/被 SIGKILL 时由内核回收，不会留下永久阻塞后续
  运行的占用；锁文件本身可安全留存。

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

## 6. 基因级携带者关联分析（assoc）

`assoc` 子命令回答“携带某个基因变异的人，在病例组是否比对照组更
多”。它**只读已发布数据集**（`manifest.json` + 样本宽表表头 + 基因型
长表），不重新跑 VCF 规范化，对已发布目录零改动。分析单位是
人 × 基因。这是**合成数据上的关联筛查演示，不输出临床诊断结论**。

### 6.1 三份分析输入（格式定义）

1. **病例/对照标签表**（TSV，带表头）：

   ```
   sample	group
   S1	case
   S2	control
   ```

   `group` 只能为 `case`/`control`。重复样本、同一样本组别冲突、
   数据集中不存在的样本都是错误，全部列出并以退出码 2 拒绝——
   绝不悄悄扩充研究人数。只有标签表中出现的人进入分析；数据集中
   未标注的样本不参与任何分母。

2. **基因区间 BED**（无表头，`#` 起注释行）：
   `chrom<TAB>start<TAB>end<TAB>gene`。零起点、右端不包含；同一
   基因可多行（区间先取并集）。拒绝负坐标、空区间（`end <= start`）
   与空基因名。

3. **每样本可判读区间 BED**（无表头）：
   `chrom<TAB>start<TAB>end<TAB>sample`。坐标语义同上；第 4 列必须
   是数据集中存在的样本。含义：该样本在该区域有足够测序覆盖。
   同一样本同一染色体的重叠/相邻区间先合并成并集再判定。

### 6.2 判定语义

- **候选位点** = 数据集中出现过的规范化变异键 `(chrom, pos, ref,
  alt)`；位点占用的参考区间为 `[pos-1, pos-1+len(ref))`（零起点
  半开），与某基因任一区间相交即归入该基因。同一基因的重叠区间
  不造成重复计数；一个位点可以属于不同基因。
- **位点有事实记录**：只有完整二倍体、每个等位 ∈ {0,1} 且 DP ≥
  `--min-dp` 的调用才有效；有效调用含 ALT → 携带，否则参考。
  半缺失、全缺失、非二倍体、低深度 → 无法判断；即使可判读区间
  覆盖它，也**不**把显式坏调用改写为参考。
- **位点无记录**（样本矩阵中的 NA）：仅当可判读区间并集**完整**
  覆盖位点参考区间才补为参考；只覆盖锚点、区间缺一段都是无法
  判断——NA 绝不直接当作没有变异。
- **人 × 基因聚合**：任一候选位点携带 → 携带（不受其它未知位点
  影响）；无携带且全部候选位点为参考 → 未携带；其余从该基因的
  分母排除。一个人可以在不同基因使用不同的分母资格；一人同一
  基因有多个变异仍只计为一个携带者。

### 6.3 统计

- 每个可检验基因一张 2×2 表，**双侧 Fisher 精确检验**：固定行列
  边际下，概率不大于观察表概率的所有表的概率和（`fractions.
  Fraction` 精确计算，纯标准库，干净环境可运行）。
- **交叉乘积优势比** `ad/bc`：零分母时区分 `inf`（分子非零）与
  `undefined`（分子也为零），不用任意伪计数替换。
- 任一组有效分母为零 → 不检验并标明；没有候选位点的基因标明
  原因。未检验项不进入检验数 m。
- 只对实际做过检验的基因统一执行 **Benjamini–Hochberg 校正**，
  按（校正值, 原始概率, 基因名）稳定排序。

### 6.4 演示（合成队列 syn37_v1，样本 S1/S2/S3）

```bash
python3 -I run_pipeline.py assoc \
  --dataset output/syn37_v1 \
  --labels data/assoc/labels_demo.tsv \
  --genes data/assoc/genes_demo.bed \
  --callable data/assoc/callable_demo.bed \
  --min-dp 10
```

演示输入的设计（均可手算核对，`tests/test_assoc.py` 逐条回归）：

| 基因 | 区间（0 起点半开） | 候选位点 | 病例（S1,S2） | 对照（S3） |
|---|---|---|---|---|
| GENE_A | chr1:0-12 | pos5、pos10×2 | 2 携带 | 1 携带 |
| GENE_B | chr1:12-28 | pos13/15/18、pos28 插入+删失 | 2 携带 | 1 未携带 |
| GENE_C | chr1:44-49 | pos46、pos48 | 1 携带 1 未携带 | 1 未携带 |
| GENE_D | chr1:60-65 | pos63×2 | 1 携带 1 排除 | 1 未携带 |
| GENE_E | chr1:28-31 | pos28 删失 | 1 携带 1 排除 | 1 未携带 |
| GENE_F | chr2:0-50 | 无候选位点 | — | — |
| GENE_G | chr1:4-10 ∪ 8-14 | pos5、pos10×2、pos13 | 2 携带 | 1 携带 |

要点：S1 在 GENE_C 有两个携带位点只计一人；S1 在 GENE_A 的携带
结论不受 pos10 半缺失位点影响；S2 的相邻可判读区间
`[40,46)+[46,50)` 合并后完整覆盖 GENE_C 两个位点才补参考；S1 在
GENE_D 的显式 `./.` 不被全覆盖的可判读区间洗成参考（排除）；删失
`GA>G` 占参考区间 `[27,29)`，跨过 GENE_B 右端与 GENE_E 左端，
同时归入两基因，而插入 `G>GA` 占 `[27,28)` 恰止于 GENE_E 左边界
不属于 GENE_E；GENE_G 两行重叠区间合并后 pos5 不重复计数。

输出（节选；完整 JSON 含每个样本 × 位点的判定明细，`--out-dir`
可另写 `assoc_summary.tsv` / `assoc_decisions.tsv`）：

```
genes_tested: 6
1 GENE_B  table=(2,0,0,1)  p=0.3333333333333333  OR=inf        adj=1.0
2 GENE_A  table=(2,0,1,0)  p=1.0                 OR=undefined  adj=1.0
3 GENE_C  table=(1,1,0,1)  p=1.0                 OR=inf        adj=1.0
4 GENE_D  table=(1,0,0,1)  p=1.0                 OR=inf        adj=1.0
5 GENE_E  table=(1,0,0,1)  p=1.0                 OR=inf        adj=1.0
6 GENE_G  table=(2,0,1,0)  p=1.0                 OR=undefined  adj=1.0
GENE_F    no_candidate_loci（数据集中没有与该基因区间相交的候选位点）
```

错误输入的实际命令结果：

```bash
# 组别冲突：退出码 2，报清样本与行号
$ python3 -I run_pipeline.py assoc --dataset output/syn37_v1 \
    --labels data/assoc/labels_conflict.tsv --genes data/assoc/genes_demo.bed \
    --callable data/assoc/callable_demo.bed --min-dp 10
关联分析输入错误: 标签表 data/assoc/labels_conflict.tsv 校验失败:
第 4 行: 样本 S1 组别冲突（先前为 case，本行为 control）
（退出码 2）

# 一组分母归零：该基因不检验，不进入 BH 的检验数
$ python3 -I run_pipeline.py assoc --dataset output/syn37_v1 \
    --labels data/assoc/labels_zero_denominator.tsv \
    --genes data/assoc/genes_demo.bed --callable data/assoc/callable_demo.bed \
    --min-dp 10
genes_tested: 5
GENE_D  zero_denominator（有效分母为零的组: case；该基因不参与检验）
GENE_F  no_candidate_loci
```

## 7. 测序数据接收（demux）：混合文库拆分与 UMI 提取

`demux` 子命令是变异管线的**上游接收环节**。实验室交来成对的
R1/R2 FASTQ，R1 的 5' 端依次是 B 碱基的样本条码和 U 碱基的 UMI
（标记原始分子的短序列）。本环节把混合文库拆回各个样本、剥离
条码/UMI 并把 UMI 记入两端头行；拆出的每样本 R1/R2 对即可交给
外部比对器与变异调用，产出的 VCF 再由 `run` 子命令规范化整合。
本环节不调用比对器，只保证接收到的数据结构完整、配对严格、
去向可审计。演示与基准只用**合成序列**。

**环境与启动**：无需安装任何第三方包——与管线其余部分一样只用
Python 3.11 标准库（`gzip`/`zlib`/`csv`/`json` 等）。在干净
Python 3.11 环境中克隆本仓库后，于项目根目录直接执行下文命令
即可（`run_pipeline.py` 会在 `-I` 隔离模式下自行定位 `src` 包）。

### 7.1 输入与参数

```bash
python3 -I run_pipeline.py demux \
  --r1 R1.fastq.gz --r2 R2.fastq.gz \
  --barcodes barcodes.tsv \
  --umi-length 4 --max-mismatches 1 --max-open 4 \
  --out-dir output/demux_demo
```

* `--r1/--r2`：普通或 gzip 压缩的 FASTQ，按内容魔数识别，与文件
  名无关；支持多 member gzip（例如追加写出的文件）。
* `--barcodes`：TSV 登记表，表头恰为 `barcode`、`sample` 两列。
  条码只能由 A/C/G/T 组成（小写归一为大写）且全部等长，其长度即
  B；一个样本可登记多个条码，同一条码不得重复登记；样本名必须能
  安全用作文件名——白名单 `[A-Za-z0-9._-]`、须以字母或数字开头，
  显式拒绝路径分隔符与 `..`，`undetermined` 为保留名。
* `--umi-length U`：UMI 长度，合法范围 **[0, 64]**；U=0 时不追加
  RX 标签。
* `--max-mismatches M`：允许的条码错配数，合法范围 **[0, B]**。
* `--max-open K`：同时打开的去向输出对上限（undetermined 也占一
  个），合法范围 **[1, 256]**。
* `--out-dir`：接收输出目录，必须不存在或为空——**非空目录在开始
  前拒绝**（退出码 2），避免把两批接收混进同一目录。

### 7.2 校验与配对语义

* FASTQ 按四行一条记录解析：头行 `@`、序列、`+`、质量。检查记录
  完整（截断即失败）、序列与质量等长、碱基 ∈ ACGTN（小写归一）、
  质量字符 ∈ Phred 33..126。
* 双端**逐条锁步对应**，不分别过滤再凑对：read 名 = 头行第一个
  空白分隔字段去掉末尾 `/1`、`/2`，两端必须相等；头行其余信息
  原样保留。
* 两端**同时**耗尽才算正常结束。任一端多出记录、损坏的 gzip 尾
  部、截断的记录、名字不一致，都使本次接收失败，错误带**文件、
  记录序号与原因**（gzip 流逐 member 手工解析并自验 CRC32/ISIZE，
  损坏被定位到正在读取的记录）。
* 全程流式：任何时刻只持有当前记录，不整载输入。

### 7.3 拆分判定（每对）

取 R1 前 B 碱基为观测条码、随后 U 碱基为 UMI，按序判定：

| # | 条件 | 去向（原因码） |
|---|---|---|
| 1 | `len(R1) < B+U` | undetermined（`read_too_short`） |
| 2 | `len(R1) == B+U`，截取后没有剩余序列 | undetermined（`no_insert`） |
| 3 | UMI 含 N | undetermined（`umi_contains_n`） |
| 4 | 最近登记条码的 Hamming 距离 > M | undetermined（`barcode_unmatched`） |
| 5 | 距离合格但最近条码属于多个不同样本 | undetermined（`ambiguous_sample`） |
| 6 | 最近条码（可多个并列）同属一个样本 | 分配该样本 |

观测值里的 N 与任何登记碱基都算一次错配。只有**恰好**达到最小距
离的登记条码算“最近”——不把阈值内所有条码都当最近，也不按登记
顺序取第一个；同一样本的多条码并列可以分配，跨样本并列不能猜。

### 7.4 输出

```
<out-dir>/
├── fastq/<去向>/R1.fastq.gz 与 R2.fastq.gz   每个去向一对
├── demux_table.tsv        逐对去向表（随处理流式写出）
└── summary.json           批次摘要
```

* 分配成功的对：R1 序列/质量同步去掉前 B+U 个字符，R2 不裁剪；
  提取到的 UMI 以 `RX:Z:<值>` 附在**两端**头行末尾（第一个空白
  字段不变，两端仍能认出同一 read）。无法分配的对两端**原样**
  写入 `undetermined`。同一去向内维持输入顺序。
* 同时打开的去向输出对不超过 K 个（LRU 回收）；关闭后重开以
  **追加**（新 gzip member）继续写——标准 gzip 读取器能一次读回
  全部内容，不会截掉之前的 read，也不重压旧文件。
* `demux_table.tsv` 列：`pair_index / read_name / observed_barcode
  / umi / min_distance / nearest_samples / destination / reason`；
  短 read 缺少的值留空。
* `summary.json`：已完整处理的对数、各去向数量、各无法分配原因
  的数量、运行状态。**只有双端均正常耗尽且全部输出成功关闭后**
  才写 `success`；输入损坏或输出错误时，已生成内容保留供排查，
  摘要明确 `failed` 并带实际进度与错误位置。
* 退出码：0 成功；1 接收失败（数据损坏/输出错误，摘要为
  `failed`）；2 接收配置错误（登记表、参数、目标目录）。

### 7.5 手算演示（data/demux，B=4，U=4，M=1）

登记表：`AAAA→ALPHA`、`AAAT→ALPHA`（同一样本多条码）、
`AACC→BETA`。10 对合成 read 的预期（`tests/test_demux.py`
逐字段回归）：

| read | R1 5' 端 | 判定 |
|---|---|---|
| read1 | `AAAA`+UMI | 精确匹配 → ALPHA，d=0 |
| read2 | `AATA` | 单错配纠正 → ALPHA，d=1 |
| read3 | `AAAG` | 距 AAAA、AAAT 同为 1，同属 ALPHA → 同样本并列，分配 |
| read4 | `AACA` | 距 AAAA(ALPHA)、AACC(BETA) 同为 1 → 跨样本并列，undetermined |
| read5 | `TTTT` | 最近 d=3 > 1 → undetermined |
| read6 | `AANA` | N 算一次错配，d=1 → ALPHA |
| read7 | `AAAA`+`ACGN` | UMI 含 N → undetermined |
| read8 | `AAA`（3 碱基） | 不足 B+U → undetermined |
| read9 | `AAAA`+`ACGT`（恰 8） | 截取后无剩余 → undetermined |
| read10 | `AACC` | 精确匹配 → BETA |

运行与直接解压核对：

```bash
python3 -I run_pipeline.py demux \
  --r1 data/demux/R1.fastq.gz --r2 data/demux/R2.fastq.gz \
  --barcodes data/demux/barcodes.tsv \
  --umi-length 4 --max-mismatches 1 --max-open 4 \
  --out-dir output/demux_demo

# summary：ALPHA=4、BETA=1、undetermined=5，五种无法分配原因各 1
gzip -dc output/demux_demo/fastq/ALPHA/R1.fastq.gz
# @read1/1 1:N:0:FCX RX:Z:ACGT      <- 原头行信息保留，UMI 附在末尾
# TTTTGGGG                          <- 前 8 个碱基（条码+UMI）已剥离
# +
# IJKLMNOP                          <- 质量同步裁掉前 8 个字符
gzip -dc output/demux_demo/fastq/ALPHA/R2.fastq.gz | head -4
# @read1/2 2:N:0:FCX RX:Z:ACGT      <- R2 头行同样带 RX，不裁剪
```

### 7.6 句柄回收演示（K=1，data/demux/alternating）

12 对在 P1、P2、undetermined 三个去向间反复交替，`--max-open 1`
迫使每个去向被反复关闭重开：

```bash
python3 -I run_pipeline.py demux \
  --r1 data/demux/alternating/R1.fastq.gz \
  --r2 data/demux/alternating/R2.fastq.gz \
  --barcodes data/demux/alternating/barcodes.tsv \
  --umi-length 2 --max-mismatches 0 --max-open 1 \
  --out-dir output/demux_alt

python3 - <<'EOF'
import gzip
for d in ("P1", "P2", "undetermined"):
    for mate in ("R1", "R2"):
        p = f"output/demux_alt/fastq/{d}/{mate}.fastq.gz"
        with gzip.open(p, "rt") as fh:      # 标准读取器读回全部 member
            names = [l.split()[0][1:] for l in fh if l.startswith("@")]
        print(d, mate, names)
EOF
# P1 R1 ['alt01/1', 'alt04/1', 'alt07/1', 'alt10/1']   （保序）
# P1 R2 ['alt01/2', 'alt04/2', 'alt07/2', 'alt10/2']   （成对）
# ... 每个文件含 4 个 gzip member = 被回收重开 4 次，内容完整
```

### 7.7 规模与内存（合成基准）

`benchmarks/run_benchmark.py` 生成固定种子的合成双端 FASTQ（48 个
样本、8bp 条码、8bp UMI、R1 插入 92bp、R2 100bp；约 85% 可分配、
10% 含 N、5% 随机条码），随后流式拆分并报告墙钟与峰值内存
（`resource.ru_maxrss`）：

```bash
python3 -I benchmarks/run_benchmark.py --pairs 120000
python3 -I benchmarks/run_benchmark.py --pairs 480000
```

本机实测（K=8，去向 49 个 > K，句柄池持续回收）：

| 对数 | 输入 gzip | 拆分耗时 | 吞吐 | 峰值内存 |
|---|---|---|---|---|
| 120,000 | 8.9 MB | ~10 s | ~1.2 万对/s | **18.5 MB** |
| 480,000 | 35.5 MB | ~79 s | ~0.6 万对/s | **18.5 MB** |

输入放大 4 倍，峰值内存不变——处理是流式的，内存占用与输入规模
无关。输出守恒（各去向对数之和 == 输入对数）与最大去向文件的可
完整读回由脚本自动断言。

## 8. 测试

```bash
python3 -I -m unittest discover -s tests -v
```

34 个手算/集成测试覆盖：手算多等位拆分、重复序列左对齐、参考碱基不符、
contig 别名、半缺失 GT、相位与剂量映射、跨 caller 冲突、规则升级重建、
守恒律，以及一个**不导入本项目 VCF/TSV 读写代码**的独立解析器对生成
VCF 与矩阵的交叉验证；另有 9 个发布完整性回归（子进程驱动）：真实
GT/AD 变更后的冲突重跑逐文件比对、写入中途 SIGKILL、相同/不同内容
并发竞争同一版本、已有分区目录与另一规则版本的查询可用性。

关联分析另有 46 个回归（`tests/test_assoc.py`）：演示队列逐条手算
（一人多位点不重复计、携带不被未知拖累、相邻区间合并覆盖、显式坏
调用不被掩码、重叠基因与跨边界删失归属）、独立于产品统计实现的小表
枚举核对双侧 Fisher、含相同概率的多基因 BH 校正与排序、全部非法输入
的报错，以及组别冲突/一组分母归零/已发布目录零改动的子进程命令回归。

接收环节另有 68 个回归（`tests/test_demux.py`）：登记表全部非法情形
（非 ACGT/不等长/重复登记/路径穿越样本名）、参数合法范围、FASTQ 逐条
解析（截断、序列质量不等长、非法碱基/质量字符、CRLF、gzip 尾部截断/
CRC 损坏/尾部垃圾的精确记录序号、多 member 读回）、双端错名与一端先
结束的失败摘要、10 对手算样例逐字段核对（含解压核对 read 名/序列/
质量/RX 标签）、K=1 句柄回收后的保序成对与 member 计数、句柄池单元
行为、U=0 与错配边界、目标目录规则，以及 CLI 退出码（0/1/2）的子进程
回归。
