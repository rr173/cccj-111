"""vcfpipe - 确定性 VCF 变异规范化与队列整合管线。

模块划分
--------
model       数据记录（源记录、基因型、规范化记录、隔离记录）
fasta       参考 FASTA 读取与参考碱基校验
contigs     contig 字典与跨 assembly 别名解析
vcf_io      VCF 头/记录解析与写出
genotype    多等位 GT/AD 拆分与剂量映射
normalize   最左对齐 + 最小表示
pipeline    端到端编排、合并、分区与清单
query       变异查询 / 样本切片 / source-to-output 追踪
assoc       基因级携带者关联分析（Fisher 精确检验 + BH 校正）
demux       测序数据接收：双端 FASTQ 校验、按条码拆分、UMI 提取
trio        家系三人组遗传一致性筛查（补全枚举，流式）
annot       转录本编码后果注释（剪接、坐标换算、密码子后果）
replan      分区重算规划：内容指纹账本、依赖闭包规划与完成确认
"""

__version__ = "1.2.0"
