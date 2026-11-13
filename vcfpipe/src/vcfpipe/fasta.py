"""参考 FASTA 读取与参考碱基校验。

只实现教学/管线所需的最小子集：
  * 解析 ``>name`` 头（取第一个空白前的 token 作为 contig 名）；
  * 序列大小写保留，但比较时统一大写（合成材料本身全大写）；
  * 1-based 闭区间取子串；
  * 计算文件级 SHA-256 与每条 contig 的长度/摘要，供清单记录。
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass


@dataclass(frozen=True)
class FastaDigest:
    path: str
    file_sha256: str
    contig_lengths: dict[str, int]
    contig_md5: dict[str, str]


class Reference:
    """内存型参考序列访问器（合成材料很小，直接全量载入）。"""

    def __init__(self, path: str):
        self.path = path
        self._seq: dict[str, list[str]] = {}
        order: list[str] = []
        current: str | None = None
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n\r")
                if not line:
                    continue
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    if current in self._seq:
                        raise ValueError(f"FASTA 中 contig 重复: {current}")
                    self._seq[current] = []
                    order.append(current)
                else:
                    if current is None:
                        raise ValueError("FASTA 在第一个 > 头之前出现序列行")
                    self._seq[current].append(line.strip())
        if not self._seq:
            raise ValueError(f"FASTA 无任何 contig: {path}")
        self._order = tuple(order)
        self._joined = {c: "".join(parts).upper() for c, parts in self._seq.items()}
        self.digest = self._digest()

    # ------------------------------------------------------------------ #
    @property
    def contigs(self) -> tuple[str, ...]:
        return self._order

    def length(self, chrom: str) -> int:
        return len(self._joined[chrom])

    def base(self, chrom: str, pos: int) -> str:
        """取 1-based 单碱基。越界抛出 IndexError。"""
        seq = self._joined[chrom]
        if pos < 1 or pos > len(seq):
            raise IndexError(f"参考位置越界: {chrom}:{pos} (len={len(seq)})")
        return seq[pos - 1]

    def slice(self, chrom: str, start: int, end: int) -> str:
        """1-based 闭区间 [start, end]。"""
        seq = self._joined[chrom]
        if start < 1 or end > len(seq) or start > end:
            raise IndexError(
                f"参考切片越界: {chrom}:{start}-{end} (len={len(seq)})"
            )
        return seq[start - 1:end]

    def has_contig(self, chrom: str) -> bool:
        return chrom in self._joined

    # ------------------------------------------------------------------ #
    def verify_ref(self, chrom: str, pos: int, ref: str) -> tuple[bool, str]:
        """校验 VCF REF 是否与参考序列一致，返回 (是否一致, 实际参考串)。"""
        end = pos + len(ref) - 1
        try:
            actual = self.slice(chrom, pos, end)
        except IndexError:
            return False, "<out_of_range>"
        return actual == ref.upper(), actual

    # ------------------------------------------------------------------ #
    def _digest(self) -> FastaDigest:
        with open(self.path, "rb") as fh:
            file_hash = hashlib.sha256(fh.read()).hexdigest()
        lengths = {c: len(s) for c, s in self._joined.items()}
        md5s = {c: hashlib.md5(s.encode("ascii")).hexdigest() for c, s in
                self._joined.items()}
        return FastaDigest(self.path, file_hash, lengths, md5s)
