"""只读查询 API：变异范围查询、样本切片、source-to-output 追踪。

直接消费管线写出的 TSV/索引/清单，不依赖运行时内存状态，
可对已落盘的数据集版本独立使用。
"""
from __future__ import annotations

import csv
import json
import os
from typing import Optional


class Dataset:
    def __init__(self, dataset_dir: str):
        self.dir = dataset_dir
        with open(os.path.join(dataset_dir, "manifest.json"),
                  "r", encoding="utf-8") as fh:
            self.manifest = json.load(fh)
        self.dataset_id = self.manifest["dataset_id"]

    # ------------------------------------------------------------------ #
    def _read_tsv(self, rel_path: str) -> list[dict[str, str]]:
        with open(os.path.join(self.dir, rel_path),
                  "r", encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    # ------------------------------------------------------------------ #
    def query_range(
        self, chrom: str, start: int, end: Optional[int] = None,
    ) -> list[dict[str, str]]:
        """按位点范围读取分区事实表；end 缺省时取单点。

        先查 ``_index.tsv`` 决定需要打开哪些分区，避免全表扫描。
        """
        end = start if end is None else end
        index = {r["chrom"]: r for r in
                 self._read_tsv(os.path.join("variants", "_index.tsv"))}
        if chrom not in index:
            return []
        entry = index[chrom]
        if end < int(entry["min_pos"]) or start > int(entry["max_pos"]):
            return []
        rows = self._read_tsv(entry["path"])
        return [r for r in rows if start <= int(r["pos"]) <= end]

    # ------------------------------------------------------------------ #
    def sample_slice(self, sample: str) -> list[dict[str, str]]:
        """样本切片：该样本在整个数据集中的全部基因型行。"""
        rows = self._read_tsv(os.path.join("matrix", "genotypes_long.tsv"))
        return [r for r in rows if r["sample"] == sample]

    def matrix_wide(self) -> list[dict[str, str]]:
        return self._read_tsv(os.path.join("matrix", "sample_gt_matrix.tsv"))

    # ------------------------------------------------------------------ #
    def trace_output(self, chrom: str, pos: int, ref: str,
                     alt: str) -> list[dict[str, str]]:
        """稳定键 -> 来源记录（lineage 正向追踪的逆查）。"""
        rows = self._read_tsv("lineage.tsv")
        return [r for r in rows if r["status"] == "included"
                and r["output_chrom"] == chrom and int(r["output_pos"]) == pos
                and r["output_ref"] == ref and r["output_alt"] == alt]

    def trace_source(self, source_path: str,
                     line_no: Optional[int] = None) -> list[dict[str, str]]:
        """源文件（可指定行号）-> 全部输出/隔离去向。"""
        rows = self._read_tsv("lineage.tsv")
        out = [r for r in rows if r["source_path"] == source_path]
        if line_no is not None:
            out = [r for r in out if int(r["source_line"]) == line_no]
        return out

    def quarantine(self, reason_code: Optional[str] = None
                   ) -> list[dict[str, str]]:
        rows = self._read_tsv("quarantine.tsv")
        if reason_code is None:
            return rows
        return [r for r in rows if r["reason_code"] == reason_code]
