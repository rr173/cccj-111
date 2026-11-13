"""contig 字典与跨 assembly 别名解析。

字典文件为 JSON，结构：

{
  "target_assembly": "syn37",
  "aliases": {
     "syn37": {"chr1": "chr1", "1": "chr1", ...},
     "syn38": {"chr1": "chr1", "chrZ": null}
  }
}

* 别名值为字符串 => 唯一映射到目标 assembly 的 contig；
* 别名值为 null  => 该 contig 在目标 assembly 上不存在/无法唯一映射，
  引用它的记录进入隔离清单（contig_unmapped）；
* 完全未声明的 contig 同样按 contig_unmapped 隔离。
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ContigMapping:
    target_assembly: str
    aliases: dict[str, dict[str, Optional[str]]]

    @classmethod
    def load(cls, path: str) -> "ContigMapping":
        with open(path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)
        target = raw["target_assembly"]
        aliases = raw["aliases"]
        # 规范化：只允许 dict[str, dict[str, str|null]]
        clean: dict[str, dict[str, Optional[str]]] = {}
        for asm, table in aliases.items():
            clean[asm] = {k: (v if v is None else str(v)) for k, v in table.items()}
        return cls(target, clean)

    def resolve(self, assembly: str, chrom: str) -> Optional[str]:
        """返回目标 assembly contig 名；无法唯一映射时返回 None。"""
        table = self.aliases.get(assembly)
        if table is None:
            return None
        if chrom not in table:
            return None
        return table[chrom]
