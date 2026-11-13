"""重新生成 data/demux 下的演示输入（确定性，可重复运行核对）。

演示一（本目录）：B=4, U=4, m=1，10 对，覆盖精确匹配、单错配纠正、
同一样本多条码并列、不同样本并列、超出错配限值、观测条码含 N、
UMI 含 N、短 read、无剩余插入——全部可手算核对（docs/README.md
第 7 节给出每一对的预期去向）。

演示二（alternating/）：3 个去向（P1、P2、undetermined）反复交替的
12 对，配合 --max-open 1 演示句柄回收后各分流文件仍保序且成对。

运行：python3 data/demux/make_demo_inputs.py
"""
import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------- #
# 演示一：10 对手算样例（B=4, U=4）
# --------------------------------------------------------------------------- #
BARCODES_TSV = "barcode\tsample\nAAAA\tALPHA\nAAAT\tALPHA\nAACC\tBETA\n"

# (name, r1_seq, r2_seq, r1_qual, r2_qual)；R1 = 条码(4) + UMI(4) + 插入
DEMO_READS = [
    # 精确匹配 AAAA -> ALPHA，d=0；质量串特意用递增值便于手算裁剪位置
    ("read1", "AAAA" + "ACGT" + "TTTTGGGG", "AACCGGTT",
     "ABCDEFGHIJKLMNOP", "QRSTUVWX"),
    # AATA 距 AAAA 为 1（其它更远）-> 单错配纠正到 ALPHA
    ("read2", "AATA" + "GGGG" + "CCCCAAAA", "TTGGCCAA", "I" * 16, "I" * 8),
    # AAAG 距 AAAA、AAAT 同为 1，两条码同属 ALPHA -> 同样本并列，仍分配
    ("read3", "AAAG" + "TTTT" + "ACACACAC", "GGCCAATT", "I" * 16, "I" * 8),
    # AACA 距 AAAA(ALPHA)、AACC(BETA) 同为 1 -> 跨样本并列，不能猜
    ("read4", "AACA" + "CCCC" + "GGGGTTTT", "CCAATTGG", "I" * 16, "I" * 8),
    # TTTT 最近为 AAAT 距离 3 > 1 -> 超出错配限值
    ("read5", "TTTT" + "AAAA" + "CCCCCCCC", "AATTGGCC", "I" * 16, "I" * 8),
    # AANA：观测中的 N 对任何登记碱基算一次错配，距 AAAA 为 1 -> ALPHA
    ("read6", "AANA" + "CGTA" + "GGCCGGCC", "TTAACCGG", "I" * 16, "I" * 8),
    # 条码精确匹配但 UMI 含 N -> umi_contains_n
    ("read7", "AAAA" + "ACGN" + "TTTTAAAA", "GGAATTCC", "I" * 16, "I" * 8),
    # R1 仅 3 碱基 < B+U=8 -> read_too_short（观测条码都取不齐）
    ("read8", "AAA", "CCGGTTAA", "III", "I" * 8),
    # R1 恰为 B+U=8，截取后没有剩余序列 -> no_insert
    ("read9", "AAAA" + "ACGT", "TTCCAAGG", "I" * 8, "I" * 8),
    # AACC 精确匹配 -> BETA
    ("read10", "AACC" + "GGCC" + "TTTTCCCC", "AAAATTTT", "I" * 16, "I" * 8),
]

# --------------------------------------------------------------------------- #
# 演示二：去向反复交替（B=4, U=2, m=0），用于 --max-open 1
# --------------------------------------------------------------------------- #
ALT_BARCODES_TSV = "barcode\tsample\nAAAA\tP1\nCCCC\tP2\n"

ALT_BARCODES = ["AAAA", "CCCC", "TTTT"]  # TTTT 在 m=0 下无法分配
ALT_UMIS = ["AC", "CA", "GT", "TG", "AA", "CC", "GG", "TT", "AG", "GA",
            "CT", "TC"]
ALT_INSERTS = ["GATTAC", "ACGTAC", "TTGCAA", "CCATGG", "GGTAAC", "ATATGC",
               "CGCGAT", "TACGTA", "AAGCTT", "GGCCTA", "TTAAGC", "CATGCA"]
ALT_R2 = ["AACCGG", "TTGGCC", "GGCCAA", "CCAATT", "AATTGG", "TTAACC",
          "GGAATT", "CCGGTT", "TTCCAA", "AAAACC", "GGGACC", "TTTACC"]


def _records(demo_reads):
    """[(r1_text, r2_text)]，头行带 /1、/2 后缀与额外仪器字段。"""
    out = []
    for name, s1, s2, q1, q2 in demo_reads:
        r1 = f"@{name}/1 1:N:0:FCX\n{s1}\n+\n{q1}\n"
        r2 = f"@{name}/2 2:N:0:FCX\n{s2}\n+\n{q2}\n"
        out.append((r1, r2))
    return out


def _write(path, text, gz):
    if gz:
        # mtime=0：输出字节确定，演示可重复生成逐字节一致
        with gzip.GzipFile(path, "wb", mtime=0) as fh:
            fh.write(text.encode("utf-8"))
    else:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)


def _write_pair(directory, records, gz):
    r1 = "".join(a for a, _ in records)
    r2 = "".join(b for _, b in records)
    suffix = ".fastq.gz" if gz else ".fastq"
    _write(os.path.join(directory, "R1" + suffix), r1, gz)
    _write(os.path.join(directory, "R2" + suffix), r2, gz)


def main():
    # 演示一：plain（便于查看）+ gzip（演示运行使用）
    with open(os.path.join(HERE, "barcodes.tsv"), "w", encoding="utf-8",
              newline="") as fh:
        fh.write(BARCODES_TSV)
    records = _records(DEMO_READS)
    _write_pair(HERE, records, gz=False)
    _write_pair(HERE, records, gz=True)

    # 演示二：去向交替，仅 gzip
    alt_dir = os.path.join(HERE, "alternating")
    os.makedirs(alt_dir, exist_ok=True)
    with open(os.path.join(alt_dir, "barcodes.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        fh.write(ALT_BARCODES_TSV)
    alt_reads = []
    for i in range(12):
        name = f"alt{i + 1:02d}"
        r1_seq = ALT_BARCODES[i % 3] + ALT_UMIS[i] + ALT_INSERTS[i]
        alt_reads.append((name, r1_seq, ALT_R2[i],
                          "I" * len(r1_seq), "I" * len(ALT_R2[i])))
    _write_pair(alt_dir, _records(alt_reads), gz=True)

    print(f"演示输入已写入 {HERE}（demo 10 对 + alternating 12 对）")


if __name__ == "__main__":
    main()
