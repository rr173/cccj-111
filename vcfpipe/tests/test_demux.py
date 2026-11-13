"""demux 接收流程的手算/集成回归。

覆盖：登记表与参数校验、FASTQ 逐条解析（含 gzip 损坏精确定位）、双端
配对失败语义、10 对手算样例的逐字段核对、K=1 句柄回收后的保序成对、
失败摘要与目标目录规则、CLI 退出码。
"""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RUNNER = os.path.join(ROOT, "run_pipeline.py")
DEMO = os.path.join(ROOT, "data", "demux")

sys.path.insert(0, os.path.join(ROOT, "src"))

from vcfpipe.demux import (  # noqa: E402
    DemuxInputError, FastqFormatError, _DestinationPool, iter_fastq,
    load_barcode_table, run_demux,
)

BARCODES_4 = "barcode\tsample\nAAAA\tALPHA\nAAAT\tALPHA\nAACC\tBETA\n"


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
def fq_record(name, seq, qual=None, suffix="", extra=""):
    qual = qual if qual is not None else "I" * len(seq)
    header = f"@{name}{suffix}" + (f" {extra}" if extra else "")
    return f"{header}\n{seq}\n+\n{qual}\n"


def fq_pair(name, s1, s2, q1=None, q2=None):
    return (fq_record(name, s1, q1, "/1", "1:N:0:FCX"),
            fq_record(name, s2, q2, "/2", "2:N:0:FCX"))


def write_text(path, text):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)


def write_gz(path, text):
    with gzip.GzipFile(path, "wb", mtime=0) as fh:
        fh.write(text.encode("utf-8"))


def read_gz_text(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def read_gz_records(path):
    """标准 gzip 读取器读回全部记录：[(header, seq, plus, qual)]。"""
    lines = read_gz_text(path).splitlines()
    assert len(lines) % 4 == 0, f"{path}: 记录数不是 4 的倍数"
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def count_gzip_members(path):
    """数 gzip member 个数（追加写一次 = 多一个 member）。"""
    with open(path, "rb") as fh:
        raw = fh.read()
    n, buf = 0, raw
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    while buf:
        dec.decompress(buf)
        assert dec.eof, f"{path}: member 不完整"
        n += 1
        buf = dec.unused_data
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    return n


class TmpDirCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="demux_test_")
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def path(self, *parts):
        return os.path.join(self.tmp, *parts)

    def write_barcodes(self, text=BARCODES_4):
        p = self.path("barcodes.tsv")
        write_text(p, text)
        return p

    def write_pair(self, r1_text, r2_text, gz=False):
        p1, p2 = self.path("R1.fastq"), self.path("R2.fastq")
        (write_gz if gz else write_text)(p1, r1_text)
        (write_gz if gz else write_text)(p2, r2_text)
        return p1, p2


# --------------------------------------------------------------------------- #
# 条码登记表校验
# --------------------------------------------------------------------------- #
class BarcodeTableTest(TmpDirCase):
    def load(self, text):
        return load_barcode_table(self.write_barcodes(text))

    def test_valid_multi_barcode_per_sample(self):
        reg = self.load(BARCODES_4)
        self.assertEqual(reg.barcode_length, 4)
        self.assertEqual(reg.barcode_to_sample["AAAA"], "ALPHA")
        self.assertEqual(reg.barcode_to_sample["AAAT"], "ALPHA")
        self.assertEqual(reg.barcode_to_sample["AACC"], "BETA")

    def test_header_column_order_free(self):
        reg = self.load("sample\tbarcode\nALPHA\tAAAA\n")
        self.assertEqual(reg.barcode_to_sample, {"AAAA": "ALPHA"})

    def test_lowercase_barcode_normalized(self):
        reg = self.load("barcode\tsample\nacgt\tALPHA\n")
        self.assertEqual(reg.barcode_to_sample, {"ACGT": "ALPHA"})

    def test_reject_non_acgt_barcode(self):
        with self.assertRaisesRegex(DemuxInputError, "第 2 行.*A/C/G/T"):
            self.load("barcode\tsample\nACGX\tALPHA\n")

    def test_reject_empty_barcode(self):
        with self.assertRaisesRegex(DemuxInputError, "第 2 行"):
            self.load("barcode\tsample\n\tALPHA\n")

    def test_reject_unequal_lengths(self):
        with self.assertRaisesRegex(DemuxInputError, "第 3 行.*长度 5.*4"):
            self.load("barcode\tsample\nAAAA\tALPHA\nACGTA\tBETA\n")

    def test_reject_duplicate_barcode_same_sample(self):
        with self.assertRaisesRegex(DemuxInputError, "第 3 行.*重复登记"):
            self.load("barcode\tsample\nAAAA\tALPHA\nAAAA\tALPHA\n")

    def test_reject_duplicate_barcode_across_samples(self):
        with self.assertRaisesRegex(DemuxInputError, "重复登记.*ALPHA"):
            self.load("barcode\tsample\nAAAA\tALPHA\nAAAA\tBETA\n")

    def test_reject_unsafe_sample_names(self):
        for bad in ("../evil", "a/b", "a\\b", ".", "..", "undetermined",
                    "Undetermined", "has space", "-leading", ".hidden",
                    "中文名", "semi;colon"):
            with self.assertRaises(DemuxInputError, msg=bad):
                self.load(f"barcode\tsample\nAAAA\t{bad}\n")

    def test_accept_safe_sample_names(self):
        reg = self.load("barcode\tsample\nAAAA\tS1.v2_3-4\n")
        self.assertEqual(reg.barcode_to_sample["AAAA"], "S1.v2_3-4")

    def test_reject_bad_header(self):
        for header in ("barcode\n", "barcode\tsample\textra\n",
                       "code\tsample\n"):
            with self.assertRaisesRegex(DemuxInputError, "表头"):
                self.load(header + "AAAA\tALPHA\n")

    def test_reject_empty_table(self):
        with self.assertRaisesRegex(DemuxInputError, "没有任何条码"):
            self.load("barcode\tsample\n")

    def test_reject_wrong_column_count(self):
        with self.assertRaisesRegex(DemuxInputError, "第 2 行.*2 列"):
            self.load("barcode\tsample\nAAAA\n")

    def test_missing_file(self):
        with self.assertRaisesRegex(DemuxInputError, "无法读取"):
            load_barcode_table(self.path("nope.tsv"))


# --------------------------------------------------------------------------- #
# 参数合法范围
# --------------------------------------------------------------------------- #
class ParamRangeTest(TmpDirCase):
    def setUp(self):
        super().setUp()
        self.bc = self.write_barcodes()
        self.r1, self.r2 = self.write_pair(
            fq_pair("r1", "AAAACCCCGGGG", "TTTT")[0]
            + fq_pair("r2", "AAAACCCCGGGG", "TTTT")[0],
            fq_pair("r1", "AAAACCCCGGGG", "TTTT")[1]
            + fq_pair("r2", "AAAACCCCGGGG", "TTTT")[1],
        )

    def invoke(self, **kw):
        args = dict(umi_length=4, max_mismatches=1, max_open=2)
        args.update(kw)
        return run_demux(self.r1, self.r2, self.bc,
                         self.path("out"), **args)

    def test_umi_length_range(self):
        for bad in (-1, 65, 1000):
            with self.assertRaisesRegex(DemuxInputError, "UMI 长度"):
                self.invoke(umi_length=bad)

    def test_max_mismatches_range(self):
        for bad in (-1, 5):  # B=4
            with self.assertRaisesRegex(DemuxInputError, "错配上限"):
                self.invoke(max_mismatches=bad)

    def test_max_open_range(self):
        for bad in (0, 257, -3):
            with self.assertRaisesRegex(DemuxInputError, "输出对上限"):
                self.invoke(max_open=bad)

    def test_config_error_creates_no_output(self):
        with self.assertRaises(DemuxInputError):
            self.invoke(umi_length=65)
        self.assertFalse(os.path.exists(self.path("out")))

    def test_boundary_values_accepted(self):
        s = self.invoke(umi_length=0, max_mismatches=4, max_open=1)
        self.assertEqual(s["status"], "success")
        self.assertEqual(s["pairs_completed"], 2)


# --------------------------------------------------------------------------- #
# FASTQ 解析
# --------------------------------------------------------------------------- #
class FastqParserTest(TmpDirCase):
    def parse(self, text, gz=False, name="R1.fastq"):
        p = self.path(name)
        (write_gz if gz else write_text)(p, text)
        return list(iter_fastq(p))

    def test_plain_and_gzip_roundtrip(self):
        text = fq_pair("r1", "ACGT", "TTGG")[0] + fq_record("r2/1", "CC")
        for gz in (False, True):
            recs = self.parse(text, gz=gz)
            self.assertEqual(len(recs), 2)
            self.assertEqual(recs[0].read_name, "r1")
            self.assertEqual(recs[1].read_name, "r2")

    def test_gzip_detected_by_content_not_name(self):
        # gzip 内容但文件名不带 .gz
        recs = self.parse(fq_record("r/1", "ACGT"), gz=True, name="x.txt")
        self.assertEqual(recs[0].seq, "ACGT")
        # 普通文本但文件名带 .gz
        recs = self.parse(fq_record("r/1", "ACGT"), gz=False, name="x.gz")
        self.assertEqual(recs[0].seq, "ACGT")

    def test_crlf_and_missing_final_newline(self):
        text = "@r/1\r\nACGT\r\n+\r\nIIII\r\n@r2/1\nACGT\n+\nIIII"
        recs = self.parse(text)
        self.assertEqual([r.read_name for r in recs], ["r", "r2"])
        self.assertEqual(recs[0].seq, "ACGT")

    def test_lowercase_bases_normalized(self):
        recs = self.parse(fq_record("r/1", "acgtn"))
        self.assertEqual(recs[0].seq, "ACGTN")

    def test_read_name_strips_suffix_keeps_extra(self):
        recs = self.parse("@INST:1:FC:1:1101:1000:1000/1 1:N:0:FCX extra\n"
                          "ACGT\n+\nIIII\n")
        self.assertEqual(recs[0].read_name,
                         "INST:1:FC:1:1101:1000:1000")
        self.assertEqual(recs[0].header,
                         "@INST:1:FC:1:1101:1000:1000/1 1:N:0:FCX extra")

    def test_empty_file_yields_nothing(self):
        self.assertEqual(self.parse(""), [])

    def test_truncated_record(self):
        with self.assertRaisesRegex(FastqFormatError,
                                    r"第 1 条记录.*不足 4 行"):
            self.parse("@r/1\nACGT\n+\n")

    def test_truncated_second_record_reports_number(self):
        text = fq_record("ok/1", "ACGT") + "@bad/1\nACGT\n"
        with self.assertRaisesRegex(FastqFormatError,
                                    r"第 2 条记录.*不足 4 行"):
            self.parse(text)

    def test_seq_qual_length_mismatch(self):
        with self.assertRaisesRegex(FastqFormatError,
                                    r"第 1 条记录.*长度 4.*长度 3"):
            self.parse("@r/1\nACGT\n+\nIII\n")

    def test_illegal_base(self):
        with self.assertRaisesRegex(FastqFormatError, r"非法碱基字符 'X'"):
            self.parse("@r/1\nACXT\n+\nIIII\n")

    def test_illegal_quality_chars(self):
        for bad in (" ", "\x1f", "\x7f"):
            with self.assertRaisesRegex(FastqFormatError, "非法字符"):
                self.parse(f"@r/1\nACGT\n+\nIII{bad}\n")

    def test_missing_at_and_plus(self):
        with self.assertRaisesRegex(FastqFormatError, "缺少 @"):
            self.parse("r/1\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(FastqFormatError, r"缺少 \+"):
            self.parse("@r/1\nACGT\n-\nIIII\n")

    def test_empty_read_name(self):
        with self.assertRaisesRegex(FastqFormatError, "缺少 read 名"):
            self.parse("@\nACGT\n+\nIIII\n")

    def test_error_reports_second_record(self):
        text = fq_record("ok/1", "ACGT") + "@bad/1\nACXT\n+\nIIII\n"
        with self.assertRaisesRegex(FastqFormatError, "第 2 条记录"):
            self.parse(text)

    def test_truncated_gzip_tail_reports_record_after_last(self):
        p = self.path("R1.fastq.gz")
        write_gz(p, "".join(fq_record(f"r{i}/1", "ACGT")
                            for i in range(3)))
        with open(p, "rb") as fh:
            data = fh.read()
        with open(p, "wb") as fh:
            fh.write(data[:-8])  # 截掉 CRC/ISIZE 尾部
        with self.assertRaisesRegex(FastqFormatError, "第 4 条记录"):
            list(iter_fastq(p))

    def test_corrupt_gzip_crc(self):
        p = self.path("R1.fastq.gz")
        write_gz(p, fq_record("r/1", "ACGT"))
        with open(p, "rb") as fh:
            raw = bytearray(fh.read())
        raw[-5] ^= 0xFF  # 破坏 CRC
        with open(p, "wb") as fh:
            fh.write(bytes(raw))
        with self.assertRaisesRegex(FastqFormatError, "第 2 条记录.*gzip"):
            list(iter_fastq(p))

    def test_garbage_after_gzip_member(self):
        p = self.path("R1.fastq.gz")
        write_gz(p, fq_record("r/1", "ACGT"))
        with open(p, "ab") as fh:
            fh.write(b"garbage!")
        with self.assertRaisesRegex(FastqFormatError, "第 2 条记录.*gzip"):
            list(iter_fastq(p))

    def test_multi_member_gzip_reads_all(self):
        p = self.path("R1.fastq.gz")
        with gzip.GzipFile(p, "wb", mtime=0) as fh:
            fh.write(fq_record("a/1", "AAAA").encode())
        with gzip.GzipFile(p, "ab", mtime=0) as fh:  # 追加第二个 member
            fh.write(fq_record("b/1", "CCCC").encode())
        recs = list(iter_fastq(p))
        self.assertEqual([r.read_name for r in recs], ["a", "b"])


# --------------------------------------------------------------------------- #
# 双端配对失败语义
# --------------------------------------------------------------------------- #
class PairingFailureTest(TmpDirCase):
    def setUp(self):
        super().setUp()
        self.bc = self.write_barcodes()

    def invoke(self, r1_text, r2_text, gz=False):
        r1, r2 = self.write_pair(r1_text, r2_text, gz=gz)
        return run_demux(r1, r2, self.bc, self.path("out"),
                         umi_length=4, max_mismatches=1, max_open=2)

    def test_name_mismatch_fails_with_record_number(self):
        r1 = (fq_pair("a", "AAAACCCCGGGG", "TT")[0]
              + fq_pair("b", "AAAACCCCGGGG", "TT")[0]
              + fq_pair("c", "AAAACCCCGGGG", "TT")[0])
        r2 = (fq_pair("a", "AAAACCCCGGGG", "TT")[1]
              + fq_pair("X", "AAAACCCCGGGG", "TT")[1]
              + fq_pair("c", "AAAACCCCGGGG", "TT")[1])
        s = self.invoke(r1, r2)
        self.assertEqual(s["status"], "failed")
        err = s["error"]
        self.assertEqual(err["type"], "FastqFormatError")
        self.assertEqual(err["record_no"], 2)
        self.assertIn("'b'", err["message"])
        self.assertIn("'X'", err["message"])
        self.assertEqual(s["pairs_completed"], 1)  # 只有第 1 对完整处理
        # 失败摘要落盘，已生成的部分输出保留供排查
        with open(self.path("out", "summary.json")) as fh:
            disk = json.load(fh)
        self.assertEqual(disk["status"], "failed")
        self.assertEqual(disk["pairs_completed"], 1)
        recs = read_gz_records(self.path(
            "out", "fastq", "ALPHA", "R1.fastq.gz"))
        self.assertEqual(len(recs), 1)

    def test_r1_longer_fails(self):
        r1 = (fq_pair("a", "AAAACCCCGGGG", "TT")[0]
              + fq_pair("b", "AAAACCCCGGGG", "TT")[0])
        r2 = fq_pair("a", "AAAACCCCGGGG", "TT")[1]
        s = self.invoke(r1, r2)
        self.assertEqual(s["status"], "failed")
        self.assertEqual(s["error"]["record_no"], 2)
        self.assertIn("R1 存在多余记录", s["error"]["message"])
        self.assertTrue(s["error"]["file"].endswith("R1.fastq"))

    def test_r2_longer_fails(self):
        r1 = fq_pair("a", "AAAACCCCGGGG", "TT")[0]
        r2 = (fq_pair("a", "AAAACCCCGGGG", "TT")[1]
              + fq_pair("b", "AAAACCCCGGGG", "TT")[1])
        s = self.invoke(r1, r2)
        self.assertEqual(s["status"], "failed")
        self.assertEqual(s["error"]["record_no"], 2)
        self.assertIn("R2 存在多余记录", s["error"]["message"])

    def test_truncated_gzip_tail_fails_reception(self):
        r1p, r2p = self.write_pair(
            "".join(fq_pair(f"r{i}", "AAAACCCCGGGG", "TT")[0]
                    for i in range(3)),
            "".join(fq_pair(f"r{i}", "AAAACCCCGGGG", "TT")[1]
                    for i in range(3)), gz=True)
        # 截断 R2 的 gzip 尾部
        with open(r2p, "rb") as fh:
            data = fh.read()
        with open(r2p, "wb") as fh:
            fh.write(data[:-8])
        s = run_demux(r1p, r2p, self.bc, self.path("out"),
                      umi_length=4, max_mismatches=1, max_open=2)
        self.assertEqual(s["status"], "failed")
        self.assertEqual(s["error"]["record_no"], 4)
        self.assertTrue(s["error"]["file"].endswith("R2.fastq"))
        self.assertEqual(s["pairs_completed"], 3)


# --------------------------------------------------------------------------- #
# 演示一：10 对手算样例（直接核对仓库里的演示数据）
# --------------------------------------------------------------------------- #
class DemoHandCalcTest(TmpDirCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="demux_demo_test_")
        cls.out = os.path.join(cls.tmp, "out")
        cls.summary = run_demux(
            os.path.join(DEMO, "R1.fastq.gz"),
            os.path.join(DEMO, "R2.fastq.gz"),
            os.path.join(DEMO, "barcodes.tsv"),
            cls.out, umi_length=4, max_mismatches=1, max_open=4)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_summary_counts(self):
        s = self.summary
        self.assertEqual(s["status"], "success")
        self.assertIsNone(s["error"])
        self.assertEqual(s["pairs_completed"], 10)
        self.assertEqual(s["destinations"],
                         {"ALPHA": 4, "BETA": 1, "undetermined": 5})
        self.assertEqual(s["unassigned_reasons"], {
            "ambiguous_sample": 1, "barcode_unmatched": 1,
            "no_insert": 1, "read_too_short": 1, "umi_contains_n": 1,
        })
        self.assertEqual(s["parameters"], {
            "barcode_length": 4, "umi_length": 4,
            "max_mismatches": 1, "max_open": 4,
        })

    def test_table_rows_hand_computed(self):
        with open(os.path.join(self.out, "demux_table.tsv"),
                  newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(rows[0], [
            "pair_index", "read_name", "observed_barcode", "umi",
            "min_distance", "nearest_samples", "destination", "reason"])
        self.assertEqual(rows[1:], [
            ["1", "read1", "AAAA", "ACGT", "0", "ALPHA", "ALPHA", ""],
            ["2", "read2", "AATA", "GGGG", "1", "ALPHA", "ALPHA", ""],
            ["3", "read3", "AAAG", "TTTT", "1", "ALPHA", "ALPHA", ""],
            ["4", "read4", "AACA", "CCCC", "1", "ALPHA,BETA",
             "undetermined", "ambiguous_sample"],
            ["5", "read5", "TTTT", "AAAA", "3", "ALPHA",
             "undetermined", "barcode_unmatched"],
            ["6", "read6", "AANA", "CGTA", "1", "ALPHA", "ALPHA", ""],
            ["7", "read7", "AAAA", "ACGN", "", "",
             "undetermined", "umi_contains_n"],
            ["8", "read8", "", "", "", "",
             "undetermined", "read_too_short"],
            ["9", "read9", "AAAA", "ACGT", "", "",
             "undetermined", "no_insert"],
            ["10", "read10", "AACC", "GGCC", "0", "BETA", "BETA", ""],
        ])

    def test_alpha_output_records(self):
        recs = read_gz_records(os.path.join(
            self.out, "fastq", "ALPHA", "R1.fastq.gz"))
        self.assertEqual([r[0] for r in recs], [
            "@read1/1 1:N:0:FCX RX:Z:ACGT",
            "@read2/1 1:N:0:FCX RX:Z:GGGG",
            "@read3/1 1:N:0:FCX RX:Z:TTTT",
            "@read6/1 1:N:0:FCX RX:Z:CGTA",
        ])
        # R1 序列/质量同步裁掉前 B+U=8 个字符
        self.assertEqual([r[1] for r in recs],
                         ["TTTTGGGG", "CCCCAAAA", "ACACACAC", "GGCCGGCC"])
        self.assertEqual([r[3] for r in recs],
                         ["IJKLMNOP", "IIIIIIII", "IIIIIIII", "IIIIIIII"])

    def test_alpha_r2_untrimmed_but_tagged(self):
        recs = read_gz_records(os.path.join(
            self.out, "fastq", "ALPHA", "R2.fastq.gz"))
        self.assertEqual([r[0] for r in recs], [
            "@read1/2 2:N:0:FCX RX:Z:ACGT",
            "@read2/2 2:N:0:FCX RX:Z:GGGG",
            "@read3/2 2:N:0:FCX RX:Z:TTTT",
            "@read6/2 2:N:0:FCX RX:Z:CGTA",
        ])
        self.assertEqual([r[1] for r in recs],
                         ["AACCGGTT", "TTGGCCAA", "GGCCAATT", "TTAACCGG"])
        self.assertEqual(recs[0][3], "QRSTUVWX")

    def test_beta_output(self):
        recs = read_gz_records(os.path.join(
            self.out, "fastq", "BETA", "R1.fastq.gz"))
        self.assertEqual([r[0] for r in recs],
                         ["@read10/1 1:N:0:FCX RX:Z:GGCC"])
        self.assertEqual(recs[0][1], "TTTTCCCC")

    def test_undetermined_records_byte_identical_to_input(self):
        # undetermined 两端内容保持原样（无 RX、不裁剪）
        with gzip.open(os.path.join(DEMO, "R1.fastq.gz"), "rt") as fh:
            src = fh.read().splitlines()
        idx = {l[1:].split("/")[0]: src[i:i + 4]
               for i, l in enumerate(src) if l.startswith("@")}
        recs = read_gz_records(os.path.join(
            self.out, "fastq", "undetermined", "R1.fastq.gz"))
        names = [r[0].split("/")[0][1:] for r in recs]
        self.assertEqual(names, ["read4", "read5", "read7", "read8",
                                 "read9"])
        for rec in recs:
            name = rec[0].split("/")[0][1:]
            self.assertEqual(list(rec), idx[name])
            self.assertNotIn("RX:Z:", rec[0])

    def test_outputs_paired_per_destination(self):
        for dest in ("ALPHA", "BETA", "undetermined"):
            r1 = read_gz_records(os.path.join(
                self.out, "fastq", dest, "R1.fastq.gz"))
            r2 = read_gz_records(os.path.join(
                self.out, "fastq", dest, "R2.fastq.gz"))
            self.assertEqual(len(r1), len(r2))
            for a, b in zip(r1, r2):
                self.assertEqual(a[0].split()[0].replace("/1", ""),
                                 b[0].split()[0].replace("/2", ""))


# --------------------------------------------------------------------------- #
# 演示二：K=1 句柄回收后仍保序成对
# --------------------------------------------------------------------------- #
class AlternatingK1Test(TmpDirCase):
    def run_alt(self):
        return run_demux(
            os.path.join(DEMO, "alternating", "R1.fastq.gz"),
            os.path.join(DEMO, "alternating", "R2.fastq.gz"),
            os.path.join(DEMO, "alternating", "barcodes.tsv"),
            self.path("out"), umi_length=2, max_mismatches=0, max_open=1)

    def test_k1_recycling_preserves_order_and_pairing(self):
        s = self.run_alt()
        self.assertEqual(s["status"], "success")
        self.assertEqual(s["destinations"],
                         {"P1": 4, "P2": 4, "undetermined": 4})
        expected = {
            "P1": ["alt01", "alt04", "alt07", "alt10"],
            "P2": ["alt02", "alt05", "alt08", "alt11"],
            "undetermined": ["alt03", "alt06", "alt09", "alt12"],
        }
        for dest, names in expected.items():
            for mate, suffix in (("R1", "/1"), ("R2", "/2")):
                p = self.path("out", "fastq", dest, f"{mate}.fastq.gz")
                recs = read_gz_records(p)  # 标准读取器一次读回全部 member
                # 头行第一个空白字段 = read 名（分配成功的行尾还有 RX 标签）
                self.assertEqual(
                    [r[0].split()[0][1:-2] for r in recs], names)
                self.assertTrue(all(r[0].split()[0].endswith(suffix)
                                    for r in recs))
                # 4 次轮换 × K=1 -> 该去向被关闭重开 4 次 -> 4 个 member
                self.assertEqual(count_gzip_members(p), 4)

    def test_undetermined_also_counts_against_k(self):
        # K=1 时 undetermined 同样参与回收：成员数证明它被反复重开
        self.run_alt()
        p = self.path("out", "fastq", "undetermined", "R1.fastq.gz")
        self.assertEqual(count_gzip_members(p), 4)


# --------------------------------------------------------------------------- #
# 句柄池单元行为
# --------------------------------------------------------------------------- #
class PoolUnitTest(TmpDirCase):
    def test_lru_never_exceeds_limit_and_appends(self):
        pool = _DestinationPool(self.tmp, max_open=2)
        writes = []
        for dest in ("a", "b", "c", "a", "b", "c", "a"):
            pool.write(dest, f"@{dest}/1\nA\n+\nI\n", f"@{dest}/2\nT\n+\nI\n")
            writes.append(dest)
            self.assertLessEqual(pool.open_count, 2)
        pool.close_all()
        self.assertEqual(pool.open_count, 0)
        for dest in ("a", "b", "c"):
            recs = read_gz_records(self.path(
                "fastq", dest, "R1.fastq.gz"))
            self.assertEqual([r[0] for r in recs],
                             [f"@{dest}/1"] * writes.count(dest))

    def test_k1_reopen_appends_not_rewrites(self):
        pool = _DestinationPool(self.tmp, max_open=1)
        pool.write("x", "@r1/1\nA\n+\nI\n", "@r1/2\nT\n+\nI\n")
        pool.write("y", "@r2/1\nC\n+\nI\n", "@r2/2\nG\n+\nI\n")  # 挤掉 x
        pool.write("x", "@r3/1\nG\n+\nI\n", "@r3/2\nC\n+\nI\n")  # 重开 x
        pool.close_all()
        recs = read_gz_records(self.path("fastq", "x", "R1.fastq.gz"))
        self.assertEqual([r[0] for r in recs], ["@r1/1", "@r3/1"])
        self.assertEqual(count_gzip_members(
            self.path("fastq", "x", "R1.fastq.gz")), 2)


# --------------------------------------------------------------------------- #
# 判定边界
# --------------------------------------------------------------------------- #
class DecisionEdgeTest(TmpDirCase):
    def setUp(self):
        super().setUp()
        self.bc = self.write_barcodes()

    def run_one(self, r1_seq, r2_seq="TTTT", **kw):
        args = dict(umi_length=4, max_mismatches=1, max_open=2)
        args.update(kw)
        r1, r2 = self.write_pair(fq_pair("x", r1_seq, r2_seq)[0],
                                 fq_pair("x", r1_seq, r2_seq)[1])
        out = self.path("out")
        s = run_demux(r1, r2, self.bc, out, **args)
        self.assertEqual(s["status"], "success")
        with open(os.path.join(out, "demux_table.tsv"), newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        return s, rows[1], out

    def test_umi_length_zero_no_rx_tag(self):
        s, row, out = self.run_one("AAAACCCC", umi_length=0)
        self.assertEqual(row[7], "")  # assigned
        self.assertEqual(row[6], "ALPHA")
        self.assertEqual(row[3], "")  # umi 列空
        recs = read_gz_records(os.path.join(
            out, "fastq", "ALPHA", "R1.fastq.gz"))
        self.assertEqual(recs[0][0], "@x/1 1:N:0:FCX")  # 无 RX 标签
        self.assertEqual(recs[0][1], "CCCC")  # 仍裁掉条码 B=4

    def test_max_mismatches_equal_distance_assigns(self):
        # AATA 距 AAAA 恰为 1，m=1 -> 分配（“超过限值”才无法分配）
        _, row, _ = self.run_one("AATA" + "CCCC" + "GGGG", max_mismatches=1)
        self.assertEqual(row[6], "ALPHA")

    def test_max_mismatches_equals_barcode_length_is_legal(self):
        _, row, _ = self.run_one("TTTT" + "CCCC" + "GGGG", max_mismatches=4)
        self.assertEqual(row[6], "ALPHA")  # d=4 <= 4，最近为 AAAA/AAAT 同属 ALPHA

    def test_n_in_observed_counts_one_mismatch(self):
        # AANA 距 AAAA 为 1（N 算一次错配）-> 分配
        _, row, _ = self.run_one("AANA" + "CCCC" + "GGGG")
        self.assertEqual(row[6], "ALPHA")
        self.assertEqual(row[4], "1")

    def test_short_read_partial_barcode_recorded(self):
        # B=4 <= len=6 < B+U=8：观测条码可取齐，UMI 取不齐
        _, row, _ = self.run_one("AAAACC")
        self.assertEqual(row[7], "read_too_short")
        self.assertEqual(row[2], "AAAA")
        self.assertEqual(row[3], "")

    def test_short_read_below_barcode_length_blanks(self):
        _, row, _ = self.run_one("AA")
        self.assertEqual(row[7], "read_too_short")
        self.assertEqual(row[2], "")

    def test_assigned_pair_names_still_match_after_tagging(self):
        _, _, out = self.run_one("AAAACCCCGGGG")
        r1 = read_gz_records(os.path.join(out, "fastq", "ALPHA",
                                          "R1.fastq.gz"))
        r2 = read_gz_records(os.path.join(out, "fastq", "ALPHA",
                                          "R2.fastq.gz"))
        self.assertEqual(r1[0][0].split()[0], "@x/1")
        self.assertEqual(r2[0][0].split()[0], "@x/2")
        self.assertTrue(r1[0][0].endswith("RX:Z:CCCC"))
        self.assertTrue(r2[0][0].endswith("RX:Z:CCCC"))


# --------------------------------------------------------------------------- #
# 目标目录规则
# --------------------------------------------------------------------------- #
class OutDirTest(TmpDirCase):
    def setUp(self):
        super().setUp()
        self.bc = self.write_barcodes()
        self.r1, self.r2 = self.write_pair(
            fq_pair("x", "AAAACCCCGGGG", "TTTT")[0],
            fq_pair("x", "AAAACCCCGGGG", "TTTT")[1])

    def invoke(self, out):
        return run_demux(self.r1, self.r2, self.bc, out,
                         umi_length=4, max_mismatches=1, max_open=2)

    def test_nonempty_dir_rejected_before_start(self):
        out = self.path("out")
        os.makedirs(out)
        write_text(os.path.join(out, "keep.txt"), "old")
        with self.assertRaisesRegex(DemuxInputError, "非空"):
            self.invoke(out)
        self.assertEqual(os.listdir(out), ["keep.txt"])  # 原内容未被触碰

    def test_empty_existing_dir_accepted(self):
        out = self.path("out")
        os.makedirs(out)
        self.assertEqual(self.invoke(out)["status"], "success")

    def test_out_path_is_file_rejected(self):
        out = self.path("out")
        write_text(out, "x")
        with self.assertRaisesRegex(DemuxInputError, "不是目录"):
            self.invoke(out)

    def test_missing_input_fastq_rejected(self):
        with self.assertRaisesRegex(DemuxInputError, "不存在"):
            run_demux(self.r1, self.path("nope.fastq"), self.bc,
                      self.path("out"), umi_length=4, max_mismatches=1,
                      max_open=2)
        self.assertFalse(os.path.exists(self.path("out")))


# --------------------------------------------------------------------------- #
# CLI 端到端（子进程，核对退出码）
# --------------------------------------------------------------------------- #
class CliTest(TmpDirCase):
    def run_cli(self, *args):
        return subprocess.run(
            [sys.executable, "-I", RUNNER, "demux", *args],
            cwd=ROOT, capture_output=True, text=True)

    def test_demo_success_exit0(self):
        out = self.path("out")
        r = self.run_cli(
            "--r1", "data/demux/R1.fastq.gz",
            "--r2", "data/demux/R2.fastq.gz",
            "--barcodes", "data/demux/barcodes.tsv",
            "--umi-length", "4", "--max-mismatches", "1",
            "--max-open", "4", "--out-dir", out)
        self.assertEqual(r.returncode, 0, r.stderr)
        s = json.loads(r.stdout)
        self.assertEqual(s["status"], "success")
        self.assertEqual(s["pairs_completed"], 10)
        with open(os.path.join(out, "summary.json")) as fh:
            disk = json.load(fh)
        self.assertEqual(disk["status"], "success")

    def test_bad_barcodes_exit2(self):
        bc = self.path("bad.tsv")
        write_text(bc, "barcode\tsample\nACGX\tS1\n")
        r = self.run_cli(
            "--r1", "data/demux/R1.fastq.gz",
            "--r2", "data/demux/R2.fastq.gz",
            "--barcodes", bc,
            "--umi-length", "4", "--max-mismatches", "1",
            "--max-open", "4", "--out-dir", self.path("out"))
        self.assertEqual(r.returncode, 2)
        self.assertIn("接收配置错误", r.stderr)

    def test_param_out_of_range_exit2(self):
        r = self.run_cli(
            "--r1", "data/demux/R1.fastq.gz",
            "--r2", "data/demux/R2.fastq.gz",
            "--barcodes", "data/demux/barcodes.tsv",
            "--umi-length", "4", "--max-mismatches", "99",
            "--max-open", "4", "--out-dir", self.path("out"))
        self.assertEqual(r.returncode, 2)
        self.assertIn("错配上限", r.stderr)

    def test_corrupt_input_exit1_with_failed_summary(self):
        # 截断 R2 的 gzip 尾部
        bad = self.path("R2.fastq.gz")
        with open(os.path.join(DEMO, "R2.fastq.gz"), "rb") as fh:
            data = fh.read()
        with open(bad, "wb") as fh:
            fh.write(data[:-8])
        out = self.path("out")
        r = self.run_cli(
            "--r1", "data/demux/R1.fastq.gz", "--r2", bad,
            "--barcodes", "data/demux/barcodes.tsv",
            "--umi-length", "4", "--max-mismatches", "1",
            "--max-open", "4", "--out-dir", out)
        self.assertEqual(r.returncode, 1)
        s = json.loads(r.stdout)
        self.assertEqual(s["status"], "failed")
        self.assertEqual(s["error"]["type"], "FastqFormatError")
        self.assertEqual(s["error"]["record_no"], 11)
        with open(os.path.join(out, "summary.json")) as fh:
            disk = json.load(fh)
        self.assertEqual(disk["status"], "failed")
        self.assertEqual(disk["pairs_completed"], 10)

    def test_nonempty_out_dir_exit2(self):
        out = self.path("out")
        os.makedirs(out)
        write_text(os.path.join(out, "x"), "x")
        r = self.run_cli(
            "--r1", "data/demux/R1.fastq.gz",
            "--r2", "data/demux/R2.fastq.gz",
            "--barcodes", "data/demux/barcodes.tsv",
            "--umi-length", "4", "--max-mismatches", "1",
            "--max-open", "4", "--out-dir", out)
        self.assertEqual(r.returncode, 2)
        self.assertIn("非空", r.stderr)


if __name__ == "__main__":
    unittest.main()
