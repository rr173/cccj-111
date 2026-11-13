"""测序数据接收：双端 FASTQ 严格校验、按条码拆分混合文库、UMI 提取。

定位
----
实验室交付的成对 R1/R2 FASTQ（R1 的 5' 端为 B 碱基样本条码 + U 碱基
UMI）在本环节完成接收：逐条结构校验、双端按名锁步配对、按登记条码把
混合文库拆回各样本、剥离条码/UMI 并把 UMI 以 ``RX:Z:`` 附到两端头行。
拆出的每样本 R1/R2 对即可交给外部比对与变异调用，产出的 VCF 再由本
项目的 ``run`` 子命令做规范化整合（本环节不调用比对器）。

输入契约
--------
* 两份普通或 gzip 压缩的 FASTQ（按魔数识别，与文件名无关；支持多
  member gzip，例如追加写出的文件）；
* 条码登记表：TSV，表头恰为 ``barcode``、``sample`` 两列。条码只含
  A/C/G/T（小写归一为大写）且全部等长；一个样本可登记多个条码，同一
  条码不得重复登记；样本名必须是安全的文件名（白名单
  ``[A-Za-z0-9._-]``、须以字母或数字开头、``undetermined`` 为保留名）；
* 参数合法范围：``umi_length`` ∈ [0, 64]；``max_mismatches`` ∈
  [0, B]（B 为登记条码长度）；``max_open`` ∈ [1, 256]（同时打开的
  去向输出对数上限，undetermined 也占一个）。

FASTQ 解析与配对
----------------
四行一条记录：头行 ``@``、序列、``+``、质量。检查记录完整、序列与
质量等长、碱基 ∈ ACGTN（小写归一为大写）、质量字符 ∈ [!,~]（Phred
33..126）。双端逐条锁步对应——不分别过滤再凑对：read 名 = 头行第一
个空白分隔字段去掉末尾 ``/1``、``/2``，两端必须相等；头行其余信息
原样保留。两端**同时**耗尽才算正常结束：任一端多出记录、gzip 尾部
损坏、记录截断、名字不一致都使本次接收失败，错误带文件、记录序号
与原因。

拆分判定（每对）
----------------
取 R1 前 B 碱基为观测条码、随后 U 碱基为 UMI，按序判定：

1. ``len(R1) < B+U``                -> ``read_too_short``
2. ``len(R1) == B+U``（截取后无剩余）-> ``no_insert``
3. UMI 含 N                          -> ``umi_contains_n``
4. 最近登记条码的 Hamming 距离 > 上限 -> ``barcode_unmatched``
   （观测值中的 N 与任何登记碱基都算一次错配）
5. 距离合格但最近条码属于多个不同样本 -> ``ambiguous_sample``
6. 最近条码（可多个并列）同属一个样本 -> 分配该样本

只有**恰好**达到最小距离的登记条码算“最近”（不把阈值内所有条码都当
最近，也不按登记顺序取第一个）。

输出
----
* ``fastq/<去向>/R1.fastq.gz`` 与 ``R2.fastq.gz``：分配的对 R1 序列/
  质量同步去掉前 B+U 个字符、两端头行追加 ``RX:Z:<UMI>``（U=0 时不
  追加）；未分配的对两端原样写入 ``undetermined``。同一去向内保持
  输入顺序。同时打开的去向不超过 K 个（LRU 回收），重开以**追加**
  （新 gzip member）继续写——标准 gzip 读取器可一次读回全部内容，
  不会截掉之前的 read，也不重压旧文件；
* ``demux_table.tsv``：逐对去向表（序号、read 名、观测条码、UMI、
  最小距离、最近样本集合、去向、原因），随处理流式写出；短 read
  缺少的值留空；
* ``summary.json``：批次摘要（完整处理对数、各去向数量、各无法分配
  原因数量、运行状态）。只有双端均正常耗尽且全部输出成功关闭后才
  写 ``success``；输入损坏或输出错误时保留已生成内容供排查，摘要
  明确 ``failed`` 与实际进度。目标目录非空会在开始前拒绝。
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import struct
import zlib
from dataclasses import dataclass, field

# --------------------------------------------------------------------------- #
# 参数合法范围（CLI 帮助、校验与文档共用同一组常量）
# --------------------------------------------------------------------------- #
MIN_UMI_LENGTH = 0
MAX_UMI_LENGTH = 64
MIN_MAX_OPEN = 1
MAX_MAX_OPEN = 256

UNDETERMINED = "undetermined"

# 无法分配原因码（稳定标识符）
REASON_READ_TOO_SHORT = "read_too_short"
REASON_NO_INSERT = "no_insert"
REASON_UMI_CONTAINS_N = "umi_contains_n"
REASON_BARCODE_UNMATCHED = "barcode_unmatched"
REASON_AMBIGUOUS_SAMPLE = "ambiguous_sample"

SUMMARY_SCHEMA_VERSION = "1.0"

TABLE_COLUMNS = ("pair_index", "read_name", "observed_barcode", "umi",
                 "min_distance", "nearest_samples", "destination", "reason")

_SAMPLE_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}\Z")
_BARCODE_RE = re.compile(r"[ACGT]+\Z")
_BAD_BASE_RE = re.compile(r"[^ACGTNacgtn]")
_BAD_QUAL_RE = re.compile(r"[^\x21-\x7e]")  # Phred 33..126，即可打印 ASCII

_GZIP_MAGIC = b"\x1f\x8b"
_CHUNK = 65536


class DemuxInputError(ValueError):
    """接收配置错误：登记表、参数、样本名、目标目录（开始处理前抛出）。"""


class FastqFormatError(RuntimeError):
    """FASTQ 数据损坏：携带文件、记录序号与原因。"""

    def __init__(self, path: str, record_no: int, reason: str):
        self.path = path
        self.record_no = record_no
        self.reason = reason
        super().__init__(f"{path}: 第 {record_no} 条记录: {reason}")


class _StreamBroken(Exception):
    """gzip 流损坏/截断（内部使用，iter_fastq 转成 FastqFormatError）。"""


# --------------------------------------------------------------------------- #
# 行泵：把（可能 gzip 压缩的）文件按行产出，损坏精确定位
# --------------------------------------------------------------------------- #
class _LinePump:
    """按行产出普通或 gzip 压缩文件的内容。

    gzip 逐 member 手工解析（头部、raw deflate 数据、CRC32/ISIZE 尾部），
    流损坏/截断的错误**延迟**到已解压内容被消费完之后才抛出，使调用方
    能给出准确的记录序号。支持多 member gzip（本模块追加写出的输出
    文件即为此形态）。
    """

    def __init__(self, path: str):
        self.path = path
        self._fh = open(path, "rb")
        self._pending = self._fh.read(2)  # 已读出、未处理的压缩/原始字节
        self._is_gz = self._pending == _GZIP_MAGIC
        self._buf = bytearray()           # 已解压、未按行切出的字节
        self._file_eof = False
        self._broken: _StreamBroken | None = None
        # gzip member 解析状态机：header -> deflate -> trailer -> header ...
        self._gz_state = "header"
        self._header_buf = bytearray()
        self._trailer_buf = bytearray()
        self._dec = None                  # 当前 member 的 raw-deflate 解压器
        self._crc = 0
        self._isize = 0

    def close(self) -> None:
        self._fh.close()

    # -- gzip member 解析 ------------------------------------------------- #
    def _parse_header(self) -> bool:
        """尝试从 _header_buf 解析一个 member 头；完整返回 True。

        解析成功后 _header_buf 中剩余的是 deflate 数据，由调用方回收。
        """
        buf = self._header_buf
        if len(buf) < 10:
            return False
        if buf[:2] != _GZIP_MAGIC or buf[2] != 8:
            raise _StreamBroken("gzip member 头非法（魔数/压缩方法错误）")
        flg = buf[3]
        if flg & 0xE0:
            raise _StreamBroken("gzip member 头含保留标志位")
        pos = 10
        if flg & 0x04:                    # FEXTRA
            if len(buf) < pos + 2:
                return False
            xlen = buf[pos] | (buf[pos + 1] << 8)
            pos += 2 + xlen
            if len(buf) < pos:
                return False
        if flg & 0x08:                    # FNAME（零结尾）
            end = buf.find(0, pos)
            if end < 0:
                return False
            pos = end + 1
        if flg & 0x10:                    # FCOMMENT（零结尾）
            end = buf.find(0, pos)
            if end < 0:
                return False
            pos = end + 1
        if flg & 0x02:                    # FHCRC：头 CRC32 的低 16 位
            if len(buf) < pos + 2:
                return False
            want = buf[pos] | (buf[pos + 1] << 8)
            got = zlib.crc32(bytes(buf[:pos])) & 0xFFFF
            if want != got:
                raise _StreamBroken("gzip member 头 CRC 校验失败")
            pos += 2
        del buf[:pos]
        return True

    def _fill_gz(self) -> None:
        while self._pending and self._broken is None:
            if self._gz_state == "header":
                # 头可能跨块：先攒进 _header_buf；解析成功后剩余字节
                # 属于 deflate 数据，放回 _pending 继续消化
                self._header_buf += self._pending
                self._pending = b""
                if self._parse_header():
                    self._pending = bytes(self._header_buf)
                    self._header_buf.clear()
                    self._dec = zlib.decompressobj(-15)  # raw deflate
                    self._crc = 0
                    self._isize = 0
                    self._gz_state = "deflate"
            elif self._gz_state == "deflate":
                data, self._pending = self._pending, b""
                out = self._dec.decompress(data)
                self._crc = zlib.crc32(out, self._crc)
                self._isize = (self._isize + len(out)) & 0xFFFFFFFF
                self._buf += out
                if self._dec.eof:
                    self._pending = self._dec.unused_data
                    self._dec = None
                    self._gz_state = "trailer"
            else:  # trailer：8 字节 CRC32 + ISIZE（小端）
                take = min(8 - len(self._trailer_buf), len(self._pending))
                self._trailer_buf += self._pending[:take]
                self._pending = self._pending[take:]
                if len(self._trailer_buf) == 8:
                    crc, isize = struct.unpack("<II", self._trailer_buf)
                    self._trailer_buf.clear()
                    if crc != self._crc or isize != self._isize:
                        self._broken = _StreamBroken(
                            "gzip 尾部校验失败"
                            "（CRC32/长度不符，文件损坏）")
                        return
                    self._gz_state = "header"

    def _fill(self) -> None:
        """向 _buf 补充字节；流错误记入 _broken 延迟抛出。"""
        if not self._pending and not self._file_eof:
            chunk = self._fh.read(_CHUNK)
            if chunk:
                self._pending = chunk
            else:
                self._file_eof = True
        if not self._is_gz:
            self._buf += self._pending
            self._pending = b""
            return
        try:
            self._fill_gz()
        except (zlib.error, _StreamBroken) as exc:
            self._pending = b""
            self._broken = (exc if isinstance(exc, _StreamBroken)
                            else _StreamBroken(f"gzip 流损坏: {exc}"))
            return
        if self._file_eof and not self._pending and self._broken is None:
            # 文件读完但 member 未完整结束 -> 截断
            if self._gz_state == "deflate":
                self._broken = _StreamBroken(
                    "gzip 流在 member 结束前终止（文件截断或尾部损坏）")
            elif self._gz_state == "trailer":
                # 8 字节尾部未凑齐（包括一字节都未到）即截断
                self._broken = _StreamBroken("gzip 尾部不完整（文件截断）")
            elif self._gz_state == "header" and self._header_buf:
                self._broken = _StreamBroken(
                    "gzip member 头不完整（文件截断）")

    def readline(self) -> bytes:
        """返回含换行符的一行；末尾无换行的残余也作为一行；EOF 返回 b""。"""
        while True:
            idx = self._buf.find(b"\n")
            if idx >= 0:
                line = bytes(self._buf[:idx + 1])
                del self._buf[:idx + 1]
                return line
            if self._broken is not None:
                # 已解压内容已消费完：此时才暴露流错误，记录序号才准确
                err, self._broken = self._broken, None
                raise err
            if self._file_eof and not self._pending:
                line = bytes(self._buf)
                self._buf.clear()
                return line
            self._fill()


# --------------------------------------------------------------------------- #
# FASTQ 记录与解析
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class FastqRecord:
    """一条 FASTQ 记录（头行/+/质量行原样保留，序列归一为大写）。"""

    header: str  # 完整头行，含 @
    seq: str
    plus: str    # 完整 + 行
    qual: str

    @property
    def read_name(self) -> str:
        """头行第一个空白分隔字段去掉末尾 /1、/2（双端配对键）。"""
        token = self.header[1:].split(None, 1)[0]
        if token.endswith(("/1", "/2")):
            token = token[:-2]
        return token

    def render(self) -> str:
        return f"{self.header}\n{self.seq}\n{self.plus}\n{self.qual}\n"


def _chomp(line: str) -> str:
    if line.endswith("\n"):
        line = line[:-1]
    if line.endswith("\r"):
        line = line[:-1]
    return line


def _parse_record(path: str, rec_no: int,
                  raw_lines: list[bytes]) -> FastqRecord:
    try:
        header, seq, plus, qual = (
            _chomp(raw.decode("utf-8")) for raw in raw_lines)
    except UnicodeDecodeError as exc:
        raise FastqFormatError(
            path, rec_no, f"含非 UTF-8 字节: {exc}") from exc
    if not header.startswith("@"):
        raise FastqFormatError(path, rec_no, "头行缺少 @ 起始符")
    if not header[1:].split(None, 1):
        raise FastqFormatError(path, rec_no, "头行缺少 read 名")
    if not plus.startswith("+"):
        raise FastqFormatError(path, rec_no, "第三行缺少 + 起始符")
    if len(seq) != len(qual):
        raise FastqFormatError(
            path, rec_no,
            f"序列长度 {len(seq)} 与质量长度 {len(qual)} 不一致")
    bad = _BAD_BASE_RE.search(seq)
    if bad:
        raise FastqFormatError(
            path, rec_no, f"序列含非法碱基字符 {bad.group(0)!r}")
    bad = _BAD_QUAL_RE.search(qual)
    if bad:
        raise FastqFormatError(
            path, rec_no,
            f"质量含非法字符 U+{ord(bad.group(0)):04X}"
            f"（合法范围 Phred 33..126）")
    return FastqRecord(header=header, seq=seq.upper(), plus=plus, qual=qual)


def iter_fastq(path: str):
    """逐条产出 FastqRecord；一切损坏都转成 FastqFormatError（带记录序号）。"""
    pump = _LinePump(path)
    rec_no = 0
    try:
        while True:
            try:
                raw = pump.readline()
            except (_StreamBroken, OSError) as exc:
                raise FastqFormatError(path, rec_no + 1, str(exc)) from exc
            if raw == b"":
                return
            rec_no += 1
            lines = [raw]
            try:
                for _ in range(3):
                    line = pump.readline()
                    if line == b"":
                        raise FastqFormatError(
                            path, rec_no, "记录截断：记录不足 4 行")
                    lines.append(line)
            except (_StreamBroken, OSError) as exc:
                raise FastqFormatError(path, rec_no, str(exc)) from exc
            yield _parse_record(path, rec_no, lines)
    finally:
        pump.close()


# --------------------------------------------------------------------------- #
# 条码登记表
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BarcodeRegistry:
    """登记条码 -> 样本 的索引；全部条码统一长度 B。"""

    barcode_length: int
    barcode_to_sample: dict[str, str]

    def nearest(self, observed: str) -> tuple[int, tuple[str, ...]]:
        """返回 (最小 Hamming 距离, 达到该距离的样本集合（排序去重）)。

        观测值中的 N 与任何登记碱基都算一次错配（登记侧只含 ACGT，
        直接 != 比较即满足）；只有恰好达到最小距离的条码算“最近”。
        """
        best: int | None = None
        best_samples: set[str] = set()
        for registered, sample in self.barcode_to_sample.items():
            dist = 0
            for obs_base, reg_base in zip(observed, registered):
                if obs_base != reg_base:
                    dist += 1
                    if best is not None and dist > best:
                        break  # 已不可能进入最近集合
            else:
                if best is None or dist < best:
                    best = dist
                    best_samples = {sample}
                elif dist == best:
                    best_samples.add(sample)
        assert best is not None  # 登记表非空（载入时保证）
        return best, tuple(sorted(best_samples))


def validate_sample_name(name: str, where: str) -> None:
    """样本名必须能安全用作文件名（白名单 + 显式穿越/保留名检查）。"""
    if not name:
        raise DemuxInputError(f"{where}: 样本名为空")
    if name in (".", "..") or "/" in name or "\\" in name or "\x00" in name:
        raise DemuxInputError(
            f"{where}: 样本名 {name!r} 含路径分隔符/穿越序列，非法")
    if name.lower() == UNDETERMINED:
        raise DemuxInputError(
            f"{where}: 样本名 {name!r} 是接收流程的保留去向名")
    if not _SAMPLE_NAME_RE.match(name):
        raise DemuxInputError(
            f"{where}: 样本名 {name!r} 非法：仅允许字母/数字/._-，"
            f"须以字母或数字开头，最长 128 字符")


def load_barcode_table(path: str) -> BarcodeRegistry:
    """载入并校验条码登记表（TSV，表头恰为 barcode、sample 两列）。"""
    try:
        with open(path, "r", encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except OSError as exc:
        raise DemuxInputError(f"无法读取条码登记表 {path}: {exc}") from exc
    if not rows:
        raise DemuxInputError(f"{path}: 登记表为空")
    header = rows[0]
    if sorted(header) != ["barcode", "sample"]:
        raise DemuxInputError(
            f"{path}: 表头必须恰为 barcode、sample 两列，实际为 {header!r}")
    col_bc, col_sample = header.index("barcode"), header.index("sample")

    barcode_to_sample: dict[str, str] = {}
    barcode_length: int | None = None
    for lineno, row in enumerate(rows[1:], start=2):
        if not row:
            continue  # 容忍空行（如文件末尾的空白行）
        if len(row) != 2:
            raise DemuxInputError(
                f"{path}: 第 {lineno} 行应有 2 列，实际 {len(row)} 列")
        barcode = row[col_bc].strip().upper()
        sample = row[col_sample].strip()
        if not _BARCODE_RE.match(barcode):
            raise DemuxInputError(
                f"{path}: 第 {lineno} 行条码 {row[col_bc]!r} 非法："
                f"只能由 A/C/G/T 组成且非空")
        if barcode_length is None:
            barcode_length = len(barcode)
        elif len(barcode) != barcode_length:
            raise DemuxInputError(
                f"{path}: 第 {lineno} 行条码长度 {len(barcode)} "
                f"与先前登记的 {barcode_length} 不一致")
        validate_sample_name(sample, f"{path}: 第 {lineno} 行")
        if barcode in barcode_to_sample:
            raise DemuxInputError(
                f"{path}: 第 {lineno} 行条码 {barcode} 重复登记"
                f"（先前已登记给样本 {barcode_to_sample[barcode]}）")
        barcode_to_sample[barcode] = sample
    if not barcode_to_sample:
        raise DemuxInputError(f"{path}: 登记表没有任何条码")
    return BarcodeRegistry(barcode_length=barcode_length,
                           barcode_to_sample=barcode_to_sample)


# --------------------------------------------------------------------------- #
# 拆分判定
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Decision:
    """一对 read 的拆分结论。"""

    destination: str                 # 样本名或 undetermined
    reason: str                      # "" 表示成功分配
    observed_barcode: str            # 短 read 不足 B 时为 ""
    umi: str                         # 未提取到时为 ""
    min_distance: int | None         # 未尝试分配时为 None
    nearest_samples: tuple[str, ...]  # 未尝试分配时为 ()

    @property
    def assigned(self) -> bool:
        return self.reason == ""


def decide(registry: BarcodeRegistry, r1_seq: str, umi_length: int,
           max_mismatches: int) -> Decision:
    """按模块 docstring 的判定顺序对一对 read 的 R1 序列作出结论。"""
    bc_len = registry.barcode_length
    observed = r1_seq[:bc_len] if len(r1_seq) >= bc_len else ""
    if len(r1_seq) < bc_len + umi_length:
        return Decision(UNDETERMINED, REASON_READ_TOO_SHORT, observed,
                        "", None, ())
    umi = r1_seq[bc_len:bc_len + umi_length]
    if len(r1_seq) == bc_len + umi_length:
        return Decision(UNDETERMINED, REASON_NO_INSERT, observed,
                        umi, None, ())
    if "N" in umi:
        return Decision(UNDETERMINED, REASON_UMI_CONTAINS_N, observed,
                        umi, None, ())
    min_dist, samples = registry.nearest(observed)
    if min_dist > max_mismatches:
        return Decision(UNDETERMINED, REASON_BARCODE_UNMATCHED, observed,
                        umi, min_dist, samples)
    if len(samples) > 1:
        return Decision(UNDETERMINED, REASON_AMBIGUOUS_SAMPLE, observed,
                        umi, min_dist, samples)
    return Decision(samples[0], "", observed, umi, min_dist, samples)


# --------------------------------------------------------------------------- #
# 去向输出对：LRU 句柄池
# --------------------------------------------------------------------------- #
class _DestinationPool:
    """最多 K 个去向同时打开的 LRU 句柄池（每个去向 = R1/R2 一对）。

    重开去向时用追加模式写出新的 gzip member（不重压旧内容），标准
    gzip 读取器能把多 member 文件一次读回全部记录。
    """

    def __init__(self, root: str, max_open: int):
        self._root = root
        self._max_open = max_open
        self._open: dict[str, tuple] = {}  # 插入序即 LRU 顺序

    @property
    def open_count(self) -> int:
        return len(self._open)

    def write(self, destination: str, r1_text: str, r2_text: str) -> None:
        handles = self._open.get(destination)
        if handles is None:
            if len(self._open) >= self._max_open:
                oldest = next(iter(self._open))
                self._close_pair(self._open.pop(oldest))
            handles = self._open_pair(destination)
            self._open[destination] = handles
        else:
            # 命中即提级为最近使用（弹出重插到末尾）
            self._open[destination] = self._open.pop(destination)
        f_r1, f_r2 = handles
        f_r1.write(r1_text.encode("utf-8"))
        f_r2.write(r2_text.encode("utf-8"))

    def _open_pair(self, destination: str):
        dest_dir = os.path.join(self._root, "fastq", destination)
        os.makedirs(dest_dir, exist_ok=True)
        # 追加模式 + mtime=0：重开继续写且输出字节确定
        f_r1 = gzip.GzipFile(
            os.path.join(dest_dir, "R1.fastq.gz"), "ab", mtime=0)
        f_r2 = gzip.GzipFile(
            os.path.join(dest_dir, "R2.fastq.gz"), "ab", mtime=0)
        return f_r1, f_r2

    @staticmethod
    def _close_pair(handles) -> None:
        f_r1, f_r2 = handles
        f_r1.close()
        f_r2.close()

    def close_all(self) -> None:
        """关闭全部句柄；首个关闭错误照常抛出（接收按失败处理）。"""
        first_error: OSError | None = None
        while self._open:
            _, handles = self._open.popitem()
            try:
                self._close_pair(handles)
            except OSError as exc:
                if first_error is None:
                    first_error = exc
        if first_error is not None:
            raise first_error

    def close_quietly(self) -> None:
        try:
            self.close_all()
        except OSError:
            pass


# --------------------------------------------------------------------------- #
# 统计与摘要
# --------------------------------------------------------------------------- #
@dataclass
class _Stats:
    pairs_completed: int = 0
    destinations: dict[str, int] = field(default_factory=dict)
    unassigned_reasons: dict[str, int] = field(default_factory=dict)

    def record(self, decision: Decision) -> None:
        self.pairs_completed += 1
        self.destinations[decision.destination] = (
            self.destinations.get(decision.destination, 0) + 1)
        if not decision.assigned:
            self.unassigned_reasons[decision.reason] = (
                self.unassigned_reasons.get(decision.reason, 0) + 1)


def _build_summary(*, status, error, stats, r1_path, r2_path,
                   barcodes_path, out_dir, registry, umi_length,
                   max_mismatches, max_open) -> dict:
    if error is None:
        error_obj = None
    elif isinstance(error, FastqFormatError):
        error_obj = {
            "type": "FastqFormatError",
            "file": error.path,
            "record_no": error.record_no,
            "message": str(error),
        }
    else:
        error_obj = {"type": type(error).__name__, "message": str(error)}
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "status": status,
        "error": error_obj,
        "pairs_completed": stats.pairs_completed,
        "destinations": dict(sorted(stats.destinations.items())),
        "unassigned_reasons": dict(sorted(stats.unassigned_reasons.items())),
        "parameters": {
            "barcode_length": registry.barcode_length,
            "umi_length": umi_length,
            "max_mismatches": max_mismatches,
            "max_open": max_open,
        },
        "inputs": {"r1": r1_path, "r2": r2_path, "barcodes": barcodes_path},
        "outputs": {
            "directory": out_dir,
            "table": "demux_table.tsv",
            "fastq_layout": "fastq/<去向>/R1.fastq.gz 与 R2.fastq.gz",
        },
    }


def _write_summary(out_dir: str, summary: dict) -> None:
    payload = (json.dumps(summary, indent=2, sort_keys=True,
                          ensure_ascii=False) + "\n").encode("utf-8")
    with open(os.path.join(out_dir, "summary.json"), "wb") as fh:
        fh.write(payload)


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def _validate_params(registry: BarcodeRegistry, umi_length: int,
                     max_mismatches: int, max_open: int) -> None:
    if not MIN_UMI_LENGTH <= umi_length <= MAX_UMI_LENGTH:
        raise DemuxInputError(
            f"UMI 长度 {umi_length} 超出合法范围 "
            f"[{MIN_UMI_LENGTH}, {MAX_UMI_LENGTH}]")
    if not 0 <= max_mismatches <= registry.barcode_length:
        raise DemuxInputError(
            f"条码错配上限 {max_mismatches} 超出合法范围 "
            f"[0, {registry.barcode_length}]（登记条码长度）")
    if not MIN_MAX_OPEN <= max_open <= MAX_MAX_OPEN:
        raise DemuxInputError(
            f"同时打开输出对上限 {max_open} 超出合法范围 "
            f"[{MIN_MAX_OPEN}, {MAX_MAX_OPEN}]")


def _prepare_out_dir(out_dir: str) -> None:
    if os.path.lexists(out_dir):
        if not os.path.isdir(out_dir):
            raise DemuxInputError(f"目标 {out_dir} 已存在且不是目录")
        if os.listdir(out_dir):
            raise DemuxInputError(
                f"目标目录 {out_dir} 非空：为避免混淆批次，"
                f"接收只写入全新空目录")
    else:
        os.makedirs(out_dir)


def _process_pairs(r1_path, r2_path, registry, umi_length, max_mismatches,
                   pool, writer, stats) -> None:
    """双端锁步逐对处理（流式，不整载任何一端）。"""
    it_r1 = iter_fastq(r1_path)
    it_r2 = iter_fastq(r2_path)
    pair_index = 0
    while True:
        rec_r1 = next(it_r1, None)
        rec_r2 = next(it_r2, None)
        if rec_r1 is None and rec_r2 is None:
            return  # 双端同时耗尽：唯一正常结束
        if rec_r1 is None:
            raise FastqFormatError(
                r2_path, pair_index + 1,
                "R2 存在多余记录：R1 已先耗尽（双端记录数不一致）")
        if rec_r2 is None:
            raise FastqFormatError(
                r1_path, pair_index + 1,
                "R1 存在多余记录：R2 已先耗尽（双端记录数不一致）")
        pair_index += 1
        if rec_r1.read_name != rec_r2.read_name:
            raise FastqFormatError(
                r1_path, pair_index,
                f"双端 read 名不一致: R1={rec_r1.read_name!r} "
                f"R2={rec_r2.read_name!r}")

        decision = decide(registry, rec_r1.seq, umi_length, max_mismatches)
        if decision.assigned:
            if umi_length > 0:
                tag = f" RX:Z:{decision.umi}"
                head_r1, head_r2 = rec_r1.header + tag, rec_r2.header + tag
            else:
                head_r1, head_r2 = rec_r1.header, rec_r2.header
            trim = registry.barcode_length + umi_length
            out_r1 = (f"{head_r1}\n{rec_r1.seq[trim:]}\n{rec_r1.plus}\n"
                      f"{rec_r1.qual[trim:]}\n")
            out_r2 = (f"{head_r2}\n{rec_r2.seq}\n{rec_r2.plus}\n"
                      f"{rec_r2.qual}\n")
        else:
            out_r1, out_r2 = rec_r1.render(), rec_r2.render()

        pool.write(decision.destination, out_r1, out_r2)
        writer.writerow([
            pair_index, rec_r1.read_name, decision.observed_barcode,
            decision.umi,
            "" if decision.min_distance is None else decision.min_distance,
            ",".join(decision.nearest_samples),
            decision.destination, decision.reason,
        ])
        stats.record(decision)


def run_demux(r1_path: str, r2_path: str, barcodes_path: str, out_dir: str,
              *, umi_length: int, max_mismatches: int, max_open: int) -> dict:
    """执行接收，返回批次摘要（同时写入 ``<out_dir>/summary.json``）。

    配置/登记问题抛 :class:`DemuxInputError`（此时未创建任何输出）；
    数据损坏或输出错误不外抛——已生成内容保留供排查，返回并落盘的
    摘要标记 ``failed`` 且带实际进度。只有双端均正常耗尽且全部输出
    成功关闭后，摘要才为 ``success``。
    """
    registry = load_barcode_table(barcodes_path)
    _validate_params(registry, umi_length, max_mismatches, max_open)
    for path in (r1_path, r2_path):
        if not os.path.isfile(path):
            raise DemuxInputError(f"输入 FASTQ 不存在: {path}")
    _prepare_out_dir(out_dir)

    stats = _Stats()
    pool = _DestinationPool(out_dir, max_open)
    table_path = os.path.join(out_dir, "demux_table.tsv")
    try:
        with open(table_path, "w", encoding="utf-8", newline="") as table_fh:
            writer = csv.writer(table_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(TABLE_COLUMNS)
            _process_pairs(r1_path, r2_path, registry, umi_length,
                           max_mismatches, pool, writer, stats)
        pool.close_all()
    except (FastqFormatError, OSError) as exc:
        pool.close_quietly()
        summary = _build_summary(
            status="failed", error=exc, stats=stats, r1_path=r1_path,
            r2_path=r2_path, barcodes_path=barcodes_path, out_dir=out_dir,
            registry=registry, umi_length=umi_length,
            max_mismatches=max_mismatches, max_open=max_open)
        _write_summary(out_dir, summary)
        return summary

    summary = _build_summary(
        status="success", error=None, stats=stats, r1_path=r1_path,
        r2_path=r2_path, barcodes_path=barcodes_path, out_dir=out_dir,
        registry=registry, umi_length=umi_length,
        max_mismatches=max_mismatches, max_open=max_open)
    _write_summary(out_dir, summary)
    return summary
