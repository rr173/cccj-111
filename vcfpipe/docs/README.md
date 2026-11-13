# vcfpipe：确定性 VCF 变异规范化与队列整合管线

把来自不同实验批次、不同 caller、不同 assembly 的小型 VCF，整合为统一
参考版本的**变异事实表**与**样本基因型矩阵**；上游带一个**测序数据接收
环节**（`demux`）：校验双端 FASTQ、按条码把混合文库拆回各样本并提取
UMI，拆出的每样本 read 对即可交给外部比对/变异调用。另有**转录本编码
后果注释**（`annot`）：从变异表、转录本定义与 FASTA 出发，在本机给出
每个转录本上改动到的密码子与氨基酸，不查询在线注释库。仅依赖
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
│   ├── trio.py        家系三人组遗传一致性筛查（补全枚举，流式）
│   ├── annot.py       转录本编码后果注释（剪接、坐标、密码子后果）
│   └── cli.py         命令行入口
├── data/
│   ├── ref/           合成参考 syn37/syn38、annot_syn（注释演示用）
│   ├── input/         4 份合成 VCF、样本表、caller 优先级
│   ├── assoc/         关联分析演示输入（标签表、基因/可判读区间 BED）
│   ├── trio/          家系筛查演示：联合多样本 VCF + 四列家系表
│   ├── annot/         注释演示：转录本表、变异表、失败路径输入与
│   │                  确定性生成脚本 make_demo_inputs.py
│   └── demux/         接收演示输入（合成序列）：10 对手算样例 +
│                      alternating/ 去向交替样例 + 确定性生成脚本
├── tools/
│   ├── trio_exhaustive_check.py  独立小范围穷举核对器（不导入产品补全/
│   │                              传递函数，逐行复核明细）
│   ├── annot_independent_check.py 独立拼接整条突变 CDS 再翻译的核对器
│   ├── annot_make_notes.py        生成“核对笔记”（外显子→CDS→翻译）
│   └── run_annot_demo.py          运行注释演示并保存失败路径回归
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

## 8. 转录本编码后果注释（annot）：从变异走到氨基酸

`annot` 子命令读取**变异表**、**转录本定义**与**参考 FASTA**，逐条给出
“输入 id × ALT × 转录本”的解释：变异落在哪个区域；若在编码区，改动
了哪个密码子、哪三个基因组碱基、替换前后的密码子与氨基酸，以及后果
类别。全部材料是本地合成序列，**在本机计算**；不查询在线注释库，也
不依赖第三方后果注释器。仅依赖 Python 3.11 标准库。

跨接缝密码子、负链坐标、同一点在编码/非编码异构体中的不同后果，都
以**整段 REF 区间**判定（只看 `pos` 锚点会把这些差别丢掉）。

### 8.1 安装与运行

与管线其余部分相同：无需第三方包，在干净 Python 3.11 环境中克隆仓库，
于项目根目录执行：

```bash
python3 -I run_pipeline.py annot \
  --fasta data/ref/annot_syn.fa \
  --transcripts data/annot/transcripts_demo.tsv \
  --variants data/annot/variants_demo.tsv \
  --out-dir output/annot_demo --format tsv
```

* `--format tsv`（默认）或 `jsonl`；两种格式字段一一对应。
* `--out-dir` 必须不存在或为空；解释表**逐条写出并 flush**，不会攒到
  整批结束才落盘。`summary.json` 只在全部输入通过校验并成功写出后生成。
* 缓存参考 FASTA / 转录本表是允许的；逐条产出不缓存。
* 退出码：0 成功；2 输入结构/内容错误（消息带文件与行号）；1 输出
  写出失败（I/O）。

### 8.2 输入格式（字段固定）

**转录本 TSV**，表头恰为八列：

```
transcript	 gene	 chrom	 strand	 exon_starts	 exon_ends	 cds_start	 cds_end
```

| 字段 | 取值规则 |
|---|---|
| `transcript` | 本表内唯一异构体名，非空 |
| `gene` | 非空基因名 |
| `chrom` | 必须在 FASTA 中存在 |
| `strand` | `+` 或 `-` |
| `exon_starts` / `exon_ends` | 逗号分隔整数，**按坐标递增**；成对形成零起点、右开的非空外显子 `[start,end)`；两列表长度相等；外显子互不重叠（相邻允许） |
| `cds_start` / `cds_end` | 零起点半开；`cds_start == cds_end` 表示**非编码异构体**；非空时两端碱基（`cds_start` 与 `cds_end-1`）都必须在外显子上，且拼接总长整除 3（单片段不必整除） |

**变异 TSV**，表头恰为五列：

```
id	 chrom	 pos	 ref	 alt
```

| 字段 | 取值规则 |
|---|---|
| `id` | 本表内唯一，非空 |
| `chrom` | 必须在 FASTA 中存在 |
| `pos` | 从 1 计数 |
| `ref` / `alt` | 基因组**正链**字母，A/C/G/T（大小写归一为大写），至少 1 碱基 |
| `alt` | 可含逗号分隔的多个选择；ALT 之间互异且均异于 REF |

**FASTA**：支持折行，`>` 头行第一个空白字段为序列名；序列名不得重复；
碱基限 A/C/G/T。

### 8.3 校验与失败语义

开始写输出**之前**整体校验；以下情形一律结束运行、返回非零，并指出
来源文件与行号：

* FASTA 序列名重复、contig 缺失、序列含非法碱基；
* 转录本外显子乱序/重叠、两个 exon 列表长度不等、负坐标或空外显子、
  非法 CDS（端点不在外显子、拼接长度不整除 3、`cds_start > cds_end`）、
  contig 不存在、strand 非法、重复 transcript；
* 变异 id 重复、`pos` 非正整数、非法碱基、ALT 与 REF 重复/相同、
  contig 不存在、REF 区间越界；
* **REF 校验失败**：输入 REF 与 FASTA 对应区段不相等，报出实际参考串。

### 8.4 区域分类

按“相交异构体”逐一判断：异构体外显子跨度 `[最小起点, 最大终点)` 与
完整 REF 区间相交才纳入；同一 `chrom` 才可能相交（不同 contig 不混）。
只看 `pos` 锚点会丢掉以下差别，本实现按**整段区间**判定：

| 类别 | 条件 |
|---|---|
| `coding` | REF 完整落在同一外显子的编码片段（`[cds_start,cds_end)` 的交集）内 |
| `utr` | REF 完整落在同一外显子的非编码片段内 |
| `noncoding_exon` | 异构体本身非编码（`cds_start == cds_end`），REF 完整落在其外显子内 |
| `intron` | REF 完整落在单一内含子内 |
| `intergenic` | 没有任何相交异构体（`transcript` 为空） |
| `boundary_unsupported` | REF 横跨外显子/内含子、编码/非编码，或越出转录本跨度 |
| `unsupported_edit` | REF 与 ALT 长度不等；**优先**于其他分类，保留相交 transcript；找不到相交者 transcript 为空 |

上述规则对整个 REF 区间生效，包括只碰到边缘最后一个碱基的情形。

### 8.5 编码后果

* 只支持**等长碱基替换**。把整段替换**一次性**施加到该转录本的 CDS，
  再翻译触及的密码子。同一密码子改两处时，答案由两处共同决定，**不
  能用两个单碱基注释拼成答案**（`CGA→AGC` 是两个碱基合改才得到的
  错义，两个单改各自同义）。
* 不同 ALT、不同输入记录、不同异构体**独立处理**；一条记录多个 ALT
  时按 ALT 列表顺序各自展开；重叠异构体各自给出解释。
* 一段替换触及两个密码子时**各给一项解释**。

CDS 的重建：外显子在基因组正链给出，先取与编码范围相交的片段；**正链**
按坐标递增拼接，**负链**按坐标递减拼接各片段的**反向互补**；整条拼好
后**从头**分成三联体，阅读框贯穿接缝（例如两段分别贡献 2 和 4 个编码
碱基时，首个密码子横跨两段，第二段不从新密码子开头算）。负链 ALT 反向
互补后施加；输出的三个基因组坐标按**翻译方向**排列（负链为递减坐标）。

后果类别用标准遗传密码表（NCBI table 1），`*` 表示终止：

| 前后氨基酸 | 类别 |
|---|---|
| 未变（含 `*`→`*`） | `synonymous` |
| 非终止 → `*` | `stop_gained` |
| `*` → 非终止 | `stop_lost` |
| 其余 | `missense` |

不推测起始丢失、单倍型或致病性。

### 8.6 输出与排序

```
<out-dir>/annotation.tsv      逐行解释（或 annotation.jsonl）
<out-dir>/summary.json        区域/后果计数、输入摘要与 transcripts 列表
```

TSV 列固定为：

```
id chrom pos ref alt alt_index transcript gene strand region
cds_position codon_number codon_pos_1 codon_pos_2 codon_pos_3
codon_ref codon_alt aa_ref aa_alt consequence reason
```

* `alt_index` 从 1 起，对应该记录 ALT 列表序号；
* `cds_position` 是**首个被改动碱基**在剪接 CDS 中的 1 起位置；
  `codon_number` 从 1 起；`codon_pos_1..3` 是该密码子三个基因组坐标，
  正链递增、负链递减（按翻译方向）；
* `codon_ref` / `codon_alt` 是替换前后的密码子，`aa_ref` / `aa_alt` 为
  对应氨基酸；
* 非编码/暂不支持项同样占一行，密码子相关列留空，`reason` 说明原因；
  `intergenic` 的 `transcript`/`gene`/`strand` 为空。

排序键：**变异在输入表中的原顺序 → ALT 列表顺序 → transcript 字典序 →
密码子编号**。JSONL 每行一个对象，字段与 TSV 列一一对应；空值写 `null`。

### 8.7 演示与核对笔记

```bash
python3 -I data/annot/make_demo_inputs.py          # 生成合成材料（确定性）
python3 -I tools/run_annot_demo.py                 # 运行演示 + 保存失败回归
python3 -I tools/annot_make_notes.py               # 生成 output/annot_notes.md
python3 -I tools/annot_independent_check.py \
  --fasta data/ref/annot_syn.fa \
  --transcripts data/annot/transcripts_demo.tsv \
  --variants data/annot/variants_demo.tsv \
  --out-dir output/annot_demo
```

`output/annot_notes.md` 是“核对笔记”：先摆出外显子片段、剪接后的 CDS
与预期翻译，再运行入口并读回解释，覆盖：

* 正链跨接缝密码子 7（`chrP:29/30/91`，`TAC`）；
* 负链跨接缝密码子 8（`chrM:61/40/39`，坐标按翻译方向递减）；
* 同位点在不同异构体中的不同后果（`chrP:95 T>A`：`TXP_POS1` 为
  `coding`（`TGG→AGG`，W→R），同基因短 CDS 异构体 `TXP_POS1_UTR` 为
  `utr`，全非编码异构体 `TXP_NONCODING` 为 `noncoding_exon`）；
* 双碱基同时改才得到的氨基酸（`CGA→AGC`，R→S，两个单改各自同义）；
* 一段替换跨两个密码子（`chrP:112 ACG>TTA` → 密码子 14 与 15）；
* 终止获得（`TGG→TGA`）与丢失（`TAA→CAA`），以及同义（`GCT→GCC`，
  含 `*`→`*` 也按规则记 synonymous）；
* CDS 端点、外显子末端、跨内含子、无相交异构体、全非编码异构体。

`tools/annot_independent_check.py` **完全不导入** `vcfpipe.annot` 的代码：
自己解析 FASTA/TSV、自己拼整条突变 CDS、自己翻译，逐行逐字段比对产品
输出。`output/annot_regression/` 保存两类失败回归的实际结果：

```bash
# REF 错配：退出码 2，报出参考序列与输入 REF
$ python3 -I run_pipeline.py annot --fasta data/ref/annot_syn.fa \
    --transcripts data/annot/transcripts_demo.tsv \
    --variants data/annot/regression_ref_mismatch.tsv --out-dir /tmp/x
注释输入错误: 变异表 .../regression_ref_mismatch.tsv 第 2 行: REF 校验失败——
chrP:30-30 参考序列为 'A'，输入 REF 为 'C'

# 错误转录本定义（外显子重叠）：退出码 2，报出外显子与行号
$ python3 -I run_pipeline.py annot --fasta data/ref/annot_syn.fa \
    --transcripts data/annot/regression_bad_transcripts.tsv \
    --variants data/annot/variants_demo.tsv --out-dir /tmp/y
注释输入错误: 转录本表 .../regression_bad_transcripts.tsv 第 2 行: 外显子
起点乱序或与前一外显子重叠（20 < 30；外显子须按坐标递增且互不重叠）
```

### 8.8 回归测试

`tests/test_annot.py` 的 35 个测试覆盖：微型参考上“首密码子跨接缝”的
剪接与两段各自不整除 3 的情形、正/负链接缝与坐标方向、同位点不同
异构体、多 ALT、双碱基合改、跨两密码子、终止获得/丢失、同义（含
`*`→`*`）、CDS 两端点、外显子末端、boundary/intron/utr/intergenic/
noncoding_exon/unsupported_edit 全部类别、排序键、TSV 与 JSONL 等值、
独立核对器逐行一致、流式写出（首条 emit 后文件已可见），以及 REF 错配、
外显子重叠、缺失 contig、重复 FASTA 名、重复 id、非法碱基、CDS 长度
不整除、非空输出目录的子进程退出码与报错行号。原有入口（`run`、
`demux`、`query`、`assoc`、`trio`）的测试继续全部通过。

## 9. 测试

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

注释另有 35 个回归（`tests/test_annot.py`）：微型参考上首密码子跨接缝的
剪接、正/负链接缝与坐标方向、同位点不同异构体、多 ALT、双碱基合改、
跨两密码子、终止获得/丢失、同义（含 `*`→`*`）、CDS 端点、外显子末端、
全部区域类别、排序键、TSV 与 JSONL 等值、独立核对器逐行一致、流式
写出，以及 REF 错配/外显子重叠/重复 FASTA 名/非法碱基等失败路径的
退出码与行号（详见第 8.8 节）。

家系筛查另有 41 个回归（`tests/test_trio.py`）：补全去重与相位无关、
含缺失多等位的三类判定、低深度整调用未知与阈值翻转、非二倍体单列、
家系表全部结构性非法（缺列/三人相同/重复/矛盾父母/环/缺样本/重复样本
列）、坏 GT 语法/越界/负 DP/列数与 REF/ALT 校验的精确行号、正常后接
坏行不发布、gzip 与普通输入逐字节一致、摘要与明细对账，以及 CLI 退出
码（0/1/2）子进程回归；另由 `tools/trio_exhaustive_check.py` 做不导入
产品补全/传递函数的独立穷举逐行核对。

## 10. 家系三人组遗传一致性筛查（trio）

`trio` 子命令读取**联合分型的多样本 VCF**（普通或 gzip，按魔数识别，
与文件名无关）与一份四列家系表，逐位点、逐三人组解释二倍体遗传是否
一致。它**不依赖**规范化管线、**不拆分多等位位点**（整条多等位记录
作为一个位点判断）、**不改写输入基因型**。这是合成研究数据的质量
筛查：不据此声称亲缘关系，也不输出临床诊断。仅依赖 Python 3.11
标准库。

### 10.1 输入与参数

```bash
python3 -I run_pipeline.py trio \
    --vcf data/trio/trio_demo.vcf \
    --families data/trio/families_demo.tsv \
    --chroms chr1,chr2 --min-dp 10 \
    --out-dir output/trio_demo
```

* `--vcf`：联合分型多样本 VCF，普通或 gzip（内容识别）。
* `--families`：TSV，表头恰为 `family child father mother`。
  - 四列均必填，每行 child/father/mother 两两不同；
  - 同一子女只能登记一行父母：完全重复行、跨 family 的矛盾父母都
    拒绝（family 只是分组名，不能掩盖同一人的矛盾父母）；
  - 沿父母→子女的边不允许成环，报错给出人物链；
  - 涉及的每个人都必须在 VCF 样本列；重复样本列、找不到的人点名报出。
  - 任何关系非法都整体失败（退出码 2），不丢弃家庭后继续给成功报告。
* `--chroms`：常染色体名称白名单（逗号分隔）。白名单外的 VCF 行只
  计入跳过数，不把性染色体套二倍体规则。
* `--min-dp`：每人每位点可信调用所需最小 DP（正整数，达到即有效）。
* `--out-dir`：必须不存在或为空；先写隐藏 staging，全部输出关闭后
  单次 `os.replace` 发布，再写成功状态——后半段坏记录不会留下看似
  完成的摘要。

### 10.2 位点与调用校验（选中行）

* 列数与样本数吻合；POS 为正整数；REF 与各 ALT 是互不重复的非空
  A/C/G/T 序列；ALT 至多 4 个；
* FORMAT 必须含 GT、键不重复，各样本字段数与 FORMAT 键数一致
  （FORMAT 顺序不固定，如 `DP:GT` 亦可）；
* GT 等位编号必须落在当前 REF/ALT 的合法区间 `0..ALT数` 内——编号 2
  不会被当成另一条位点的编号 1；坏 GT 语法、混用 `|` 与 `/`、越界
  编号、负 DP、损坏选中行均使命令失败（退出码 1）并给出**文件与行号**；
* 有效的非二倍体调用（深度足够但不是两个等位位，如单倍体 `1`、
  三倍体 `1/2/3`）单列为 `unsupported` 并指明成员，**不**截成两个
  等位再算，也不混入三个判定类别。

### 10.3 缺失与深度语义

* DP 缺失或 `< --min-dp` 时，该人该位点的**整个调用视为未知**：原 GT
  原样保留，同时记录原因 `dp_missing`/`dp_below_min`；其已知等位不再
  约束补全（低深度者的已知等位不能偷偷约束结果）。
* 深度足够时 `0/.`、`./1` 保留已知那一半；单个 `.` 与 `./.`、`.|.`
  都表示完全未知。

### 10.4 判定定义（精确问题）

把父亲、母亲、孩子各自缺失的等位值，用该位点**全部合法编号**补全，
得到所有不同的**无序**二倍体基因型（相位 `|` 不影响计算，`0/1` 与
`1/0` 是同一个补全）。枚举父×母×孩子的全部完整组合，逐一检查孩子
能否由父亲贡献一个等位、母亲贡献一个等位组成（多重集合比较）：

| 条件 | 判定 |
|---|---|
| 所有补全组合都能解释 | `compatible` |
| 一个组合都解释不了 | `incompatible` |
| 部分能、部分不能 | `uncertain` |
| 可信的非二倍体调用 | `unsupported`（单列，不计入上三类） |

可解释时附**一组**父方/母方传递见证；不可解释时不编造传递。不会只
找到一种可行组合就断言一致，也不会遇到一个点就放弃整个家庭。

### 10.5 输出（写到 `--out-dir`）

* `trio_details.tsv`：按输入顺序逐位点、逐三人组的明细，列含
  `line_no/chrom/pos/ref/alt`、`family/child/father/mother`、三人的
  原 GT（`*_gt`）、原 DP（`*_dp`）、可信状态（`*_trust`）、**实际
  使用的候选基因型集合**（`*_candidates`）、`total_combinations`、
  `explainable_combinations`、`verdict`、`unsupported_members`、
  `witness`。
* `trio_summary.tsv`：每个三人组一行的四类计数与 `total`，末行
  `*TOTAL*` 为全量合计；逐行明细可与摘要逐格对上。
* `trio_status.json`：**完整读完输入并关闭输出后**才写的成功状态，
  含选中/跳过位点数、各三人组与全量计数、输出路径与免责声明。

### 10.6 演示（data/trio，min-dp=10，白名单 chr1,chr2）

家系表含共享父母的两个孩子 CH1/CH2（父母 FA/MO）与多代关系
（FA 的父母 GF/GM）。下表为 CH1 三人组在各 chr1 位点的手算结论：

| VCF 行 | 位点 | 父 GT | 母 GT | 孩子 GT | 总/可解释 | 判定 | 说明 |
|---|---|---|---|---|---|---|---|
| 5 | 100 A>G | 0/0 | 0/1 | 0/1 | 1/1 | compatible | 普通传递（母给 1） |
| 6 | 110 C>T | 0/0 | 0/0 | 0/1 | 1/0 | incompatible | 父母全参考而孩子含 ALT |
| 7 | 120 A>G,T | 1/1 | 2/2 | 1/2 | 1/1 | compatible | 父母分别贡献两个不同 ALT（多等位整体） |
| 8 | 130 G>A | 0/. | 0/0 | 1/1 | 2/0 | incompatible | 部分缺失仍确定不一致（母给不出 1） |
| 9 | 140 T>C | 0/1 | 0/1 | 0/. | 2/2 | compatible | 部分缺失但全部补全可解释 |
| 8（CH2） | 130 | 0/. | 0/0 | 0/1 | 2/1 | uncertain | 可行与不可行并存 |
| 10 | 150 C>A,T | ./2 | 0/0 | 1/1 | 3/0 | incompatible | 多等位缺失，3 个父补全全失败 |
| 11 | 160 G>A | 1/1（DP5） | 0/0 | 1/1 | 3/0 | incompatible | 父低深度整调用未知；CH2（孩子 0/0）同行则 3/2 uncertain |
| 12 | 170 A>G | 1/1 | 0/0 | 1/1 | 1/0 | incompatible | 深度足够的确定不一致 |
| 13 | 180 T>G | 1（单倍体） | 0/0 | 0/1 | — | unsupported | 非二倍体单列，成员 father |
| 14 | 190 C>G | 0\|1 | 0/1 | 1\|0 | 1/1 | compatible | 相位不影响，0/1 与 1/0 同一补全 |

`chrX:200` 不在白名单，只计跳过；`chr2:300` 在白名单内正常判断。
**共享父母的两个孩子结论不同**：110 位点 CH1 含 ALT 为
incompatible、CH2 全参考为 compatible。

**低 DP 导致的结论变化（阈值翻转）**：160 位点父 `1/1`、DP=5。
`--min-dp 10` 时父整调用未知，候选为 `0/0,0/1,1/1`：CH2 孩子 `0/0`
得到 uncertain（3 组合中 2 可解释）；把阈值降到 `--min-dp 5` 后父亲
可信、固定为 `1/1`，同一行翻转为 incompatible（`1/1 + 0/0` 凑不出
`0/0`）。原 GT `1/1` 在两种阈值下都原样保留在明细里。

### 10.7 独立穷举核对器

`tools/trio_exhaustive_check.py` **不导入** `vcfpipe.trio` 的任何补全/
传递/判定函数，只用标准库重新实现解析、无序补全与笛卡尔积枚举，对
明细逐行复核候选集合、总组合数、可解释数与判定类别：

```bash
python3 -I tools/trio_exhaustive_check.py \
    --details output/trio_demo/trio_details.tsv \
    --vcf data/trio/trio_demo.vcf --min-dp 10
# 独立穷举核对通过：33 行明细的候选集合、总组合数、可解释组合数与
# 判定类别全部一致。
```

### 10.8 失败例子与退出码

退出码：0 成功；1 记录/流损坏（坏 GT、越界、负 DP、损坏行、IO）；
2 配置/关系错误（家系表、白名单、阈值、缺样本、循环、重复样本列）。

```bash
# 循环家系：报人物链，退出码 2
#   家系筛查配置错误: ... 父母→子女关系存在环: A -> B -> A
# 找不到的人：点名样本与 VCF 样本列，退出码 2
# 越界 GT：报文件:行号与合法区间，退出码 1
#   家系筛查记录错误: x.vcf:4 样本 F 的 GT 等位编号 9 越界：... 0..1
# 正常数据后接坏行：退出码 1，out-dir 不发布、无 staging 残留
```

### 10.9 流式与内存

处理全程流式：任何时刻只持有当前 VCF 行与每三人组的四个计数器，
内存只与家庭数有关，不随位点数增长。5 万个位点的合成 VCF 实测约
3.1 秒、峰值 RSS 约 11 MB（机器相关），明细 50,000 行、摘要合计与
明细逐格一致。

