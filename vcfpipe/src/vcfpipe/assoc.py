"""基因级携带者关联分析（case/control carrier association）。

直接读取**已发布**的 vcfpipe 数据集（manifest + 样本宽表表头 + 基因型
长表），不重新运行 VCF 规范化。分析单位是 人 × 基因。仅依赖标准库；
统计量用 ``fractions.Fraction`` 精确计算，不用浮点误差冒险，也不用
任意伪计数。

三份分析输入（格式由本模块定义）
--------------------------------
1. **病例/对照标签表**（TSV，带表头）::

       sample<TAB>group
       S1<TAB>case
       S2<TAB>control

   * ``group`` 只能为 ``case`` / ``control``；
   * 重复样本、同一样本组别冲突、数据集中不存在的样本都是错误，
     全部列出并拒绝运行——绝不悄悄扩充研究人数；
   * 只有标签表中出现的人进入分析；数据集中未标注的样本不参与
     任何分母。

2. **基因区间**（BED，无表头，``#`` 起注释行）::

       chrom<TAB>start<TAB>end<TAB>gene

   零起点、右端不包含；同一基因可多行（区间先取并集再用）；拒绝
   负坐标、空区间（``end <= start``）与空基因名。

3. **每样本可判读区间**（BED，无表头）::

       chrom<TAB>start<TAB>end<TAB>sample

   坐标语义同上；第 4 列必须是数据集中存在的样本。含义：该样本在
   该区域有足够测序覆盖。同一样本同一染色体的重叠/相邻区间先合并
   成并集再判定——相邻区间拼合出的完整覆盖才算完整覆盖。

位点归属与人 × 基因判定
-----------------------
* 候选位点 = 数据集中出现过的规范化变异键 ``(chrom, pos, ref, alt)``；
  位点占用的参考区间为 ``[pos-1, pos-1+len(ref))``（零起点半开），
  与某基因任一区间相交即归入该基因；同一基因的重叠区间不造成重复
  计数；一个位点可以属于不同基因。
* 位点有该样本的事实记录时：只有完整二倍体且每个等位 ∈ {0,1} 且
  DP ≥ 命令阈值的调用才有效；有效调用含 ALT → 携带，否则参考。
  半缺失、全缺失、非二倍体、等位越界或低深度 → 无法判断；即使可
  判读区间覆盖它，也**不**把显式坏调用改写为参考。
* 位点没有该样本记录时：仅当可判读区间并集**完整**覆盖位点参考
  区间才补为参考；只覆盖锚点、区间缺一段都是无法判断（样本矩阵
  里的 NA 绝不直接当作没有变异）。
* 人 × 基因聚合：任一候选位点为携带 → 携带（不受其它未知位点
  影响）；无携带且全部候选位点为参考 → 未携带；其余从该基因的
  分母排除。一个人可以在不同基因使用不同的分母资格；一人同一
  基因有多个变异仍只计为一个携带者。

统计
----
* 每个可检验基因一张 2×2 表，双侧 Fisher 精确检验：固定行列边际
  下，概率不大于观察表概率的所有表的概率和；
* 交叉乘积优势比 ``ad/bc``：零分母时区分 ``inf``（分子非零）与
  ``undefined``（分子也为零），不用伪计数替换；
* 任一组有效分母为零 → 不检验；没有候选位点的基因标明原因；
* 只对实际做过检验的基因统一执行 Benjamini–Hochberg 校正，按
  （校正值, 原始概率, 基因名）稳定排序；未检验项不进入检验数 m。

本命令是**合成数据上的关联筛查演示**，不输出临床诊断结论。
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Optional

GROUP_CASE = "case"
GROUP_CONTROL = "control"
GROUPS = (GROUP_CASE, GROUP_CONTROL)

DISCLAIMER = "合成数据上的关联筛查演示，不构成临床诊断结论。"

# 位点级判定结果
CARRIER = "carrier"
REFERENCE = "reference"
UNKNOWN = "unknown"

# 人 × 基因聚合结果
NON_CARRIER = "non_carrier"
EXCLUDED = "excluded"


class AssocInputError(ValueError):
    """关联分析输入（标签表/区间文件/数据集读取）不合法。"""


# --------------------------------------------------------------------------- #
# 数据集读取（只读已发布文件，不重新规范化）
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, order=True)
class Locus:
    """数据集中的一个规范化变异键。"""

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def ref_interval(self) -> tuple[int, int]:
        """零起点半开参考占用区间 [pos-1, pos-1+len(ref))。"""
        return (self.pos - 1, self.pos - 1 + len(self.ref))


@dataclass(frozen=True)
class DatasetCalls:
    dataset_id: str
    samples: tuple[str, ...]                # 数据集样本名册（宽表表头）
    loci: tuple[Locus, ...]                 # 候选位点（排序去重）
    calls: dict                             # (sample, Locus) -> (gt, dp)


def load_dataset_calls(dataset_dir: str) -> DatasetCalls:
    """从已发布数据集读取样本名册、候选位点与事实基因型记录。"""

    def _need(rel: str) -> str:
        path = os.path.join(dataset_dir, rel)
        if not os.path.isfile(path):
            raise AssocInputError(f"数据集缺少已发布文件: {path}")
        return path

    with open(_need("manifest.json"), encoding="utf-8") as fh:
        dataset_id = json.load(fh).get("dataset_id", "<unknown>")

    with open(_need(os.path.join("matrix", "sample_gt_matrix.tsv")),
              encoding="utf-8", newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"), None)
    if not header or header[:5] != ["assembly", "chrom", "pos", "ref", "alt"] \
            or len(header) < 6:
        raise AssocInputError(
            f"{dataset_dir}: 样本宽表表头不符合 vcfpipe 发布格式")
    samples = tuple(header[5:])

    loci: set[Locus] = set()
    calls: dict = {}
    with open(_need(os.path.join("matrix", "genotypes_long.tsv")),
              encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            locus = Locus(row["chrom"], int(row["pos"]), row["ref"],
                          row["alt"])
            loci.add(locus)
            dp = None if row["dp"] == "." else int(row["dp"])
            calls[(row["sample"], locus)] = (row["gt"], dp)
    return DatasetCalls(dataset_id, samples, tuple(sorted(loci)), calls)


# --------------------------------------------------------------------------- #
# 输入一：病例/对照标签表
# --------------------------------------------------------------------------- #
def load_labels(path: str, dataset_samples) -> dict[str, str]:
    """读取并严格校验标签表；任何重复/冲突/未知样本都报清楚并拒绝。"""
    known = set(dataset_samples)
    try:
        with open(path, encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fieldnames = reader.fieldnames or []
            rows = list(reader)
    except OSError as exc:
        raise AssocInputError(f"无法读取标签表 {path}: {exc}") from exc

    missing_cols = {"sample", "group"} - set(fieldnames)
    if missing_cols:
        raise AssocInputError(
            f"标签表 {path} 缺列: {sorted(missing_cols)}（需要 sample/group "
            f"两列表头）")

    errors: list[str] = []
    labels: dict[str, str] = {}
    for i, row in enumerate(rows, start=2):  # 第 1 行是表头
        sample = (row.get("sample") or "").strip()
        group = (row.get("group") or "").strip()
        if not sample:
            errors.append(f"第 {i} 行: 样本名为空")
            continue
        if group not in GROUPS:
            errors.append(
                f"第 {i} 行: 不合法的组别标签 {group!r}"
                f"（只能是 {'/'.join(GROUPS)}）")
            continue
        if sample in labels:
            prev = labels[sample]
            if prev != group:
                errors.append(
                    f"第 {i} 行: 样本 {sample} 组别冲突"
                    f"（先前为 {prev}，本行为 {group}）")
            else:
                errors.append(f"第 {i} 行: 样本 {sample} 重复出现")
            continue
        if sample not in known:
            errors.append(
                f"第 {i} 行: 样本 {sample} 不存在于数据集"
                f"（数据集样本: {', '.join(sorted(known))}）")
            continue
        labels[sample] = group
    if errors:
        raise AssocInputError(f"标签表 {path} 校验失败:\n"
                              + "\n".join(errors))
    if not labels:
        raise AssocInputError(f"标签表 {path} 没有任何有效样本行")
    return labels


# --------------------------------------------------------------------------- #
# 输入二/三：BED 区间文件（零起点、右端不包含）
# --------------------------------------------------------------------------- #
def _read_bed_lines(path: str) -> list[tuple[int, list[str]]]:
    rows: list[tuple[int, list[str]]] = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line_no, line in enumerate(fh, start=1):
                line = line.rstrip("\n\r")
                if not line.strip() or line.startswith("#"):
                    continue
                rows.append((line_no, [c.strip() for c in
                                       line.split("\t")]))
    except OSError as exc:
        raise AssocInputError(f"无法读取区间文件 {path}: {exc}") from exc
    return rows


def _check_interval(path: str, line_no: int, chrom: str, start_s: str,
                    end_s: str, errors: list[str]
                    ) -> Optional[tuple[str, int, int]]:
    """校验一条 BED 区间；不合法时记录错误并返回 None。"""
    if not chrom:
        errors.append(f"第 {line_no} 行: chrom 为空")
        return None
    try:
        start, end = int(start_s), int(end_s)
    except ValueError:
        errors.append(f"第 {line_no} 行: 坐标必须是整数"
                      f"（得到 {start_s!r}, {end_s!r}）")
        return None
    if start < 0 or end < 0:
        errors.append(f"第 {line_no} 行: 拒绝负坐标 [{start}, {end})")
        return None
    if end <= start:
        errors.append(f"第 {line_no} 行: 拒绝空区间 [{start}, {end})"
                      f"（零起点半开语义要求 end > start）")
        return None
    return chrom, start, end


def merge_intervals(intervals) -> list[tuple[int, int]]:
    """零起点半开区间的并集：重叠与首尾相接（end == start）都合并。"""
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def interval_covered(merged, start: int, end: int) -> bool:
    """[start, end) 是否被并集的某个连通段**完整**覆盖。"""
    return any(s <= start and end <= e for s, e in merged)


def load_genes_bed(path: str) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """基因区间：gene -> chrom -> 已并集化的区间列表。"""
    genes: dict[str, dict[str, list[tuple[int, int]]]] = {}
    errors: list[str] = []
    for line_no, cols in _read_bed_lines(path):
        if len(cols) != 4:
            errors.append(f"第 {line_no} 行: 需要 4 列"
                          f" chrom/start/end/gene，实际 {len(cols)} 列")
            continue
        chrom, start_s, end_s, gene = cols
        iv = _check_interval(path, line_no, chrom, start_s, end_s, errors)
        if iv is None:
            continue
        if not gene:
            errors.append(f"第 {line_no} 行: 基因名为空")
            continue
        _, start, end = iv
        genes.setdefault(gene, {}).setdefault(chrom, []).append((start, end))
    if errors:
        raise AssocInputError(f"基因区间文件 {path} 校验失败:\n"
                              + "\n".join(errors))
    if not genes:
        raise AssocInputError(f"基因区间文件 {path} 没有任何有效区间")
    # 同一基因的重叠/相邻区间先合并，避免同一候选位点被重复计数
    return {g: {c: merge_intervals(ivs) for c, ivs in by_chrom.items()}
            for g, by_chrom in genes.items()}


def load_callable_bed(path: str, dataset_samples
                      ) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """每样本可判读区间：sample -> chrom -> 已并集化的区间列表。"""
    known = set(dataset_samples)
    per_sample: dict[str, dict[str, list[tuple[int, int]]]] = {}
    errors: list[str] = []
    for line_no, cols in _read_bed_lines(path):
        if len(cols) != 4:
            errors.append(f"第 {line_no} 行: 需要 4 列"
                          f" chrom/start/end/sample，实际 {len(cols)} 列")
            continue
        chrom, start_s, end_s, sample = cols
        iv = _check_interval(path, line_no, chrom, start_s, end_s, errors)
        if iv is None:
            continue
        if not sample:
            errors.append(f"第 {line_no} 行: 样本名为空")
            continue
        if sample not in known:
            errors.append(
                f"第 {line_no} 行: 可判读区间的样本 {sample} 不存在于数据集"
                f"（数据集样本: {', '.join(sorted(known))}）")
            continue
        _, start, end = iv
        per_sample.setdefault(sample, {}).setdefault(chrom, []).append(
            (start, end))
    if errors:
        raise AssocInputError(f"可判读区间文件 {path} 校验失败:\n"
                              + "\n".join(errors))
    # 空的可判读文件合法：语义为“任何样本都没有可判读区域”，
    # 此时所有无记录位点均为无法判断，绝不补参考。
    return {s: {c: merge_intervals(ivs) for c, ivs in by_chrom.items()}
            for s, by_chrom in per_sample.items()}


# --------------------------------------------------------------------------- #
# 位点级与人 × 基因级判定
# --------------------------------------------------------------------------- #
def parse_gt_alleles(gt: str) -> tuple[Optional[int], ...]:
    """把长表 GT 字符串解析为等位元组（None = 缺失等位）。"""
    alleles: list[Optional[int]] = []
    for part in gt.replace("|", "/").split("/"):
        if part == ".":
            alleles.append(None)
        else:
            try:
                alleles.append(int(part))
            except ValueError as exc:
                raise AssocInputError(
                    f"数据集中无法解析的 GT: {gt!r}") from exc
    return tuple(alleles)


def classify_call(gt: str, dp: Optional[int], min_dp: int
                  ) -> tuple[str, str]:
    """有位点事实记录时的判定：返回 (decision, reason)。

    只有完整二倍体、每个等位 ∈ {0,1} 且 DP 达标的调用才有效；
    显式坏调用一律“无法判断”，绝不因可判读区间覆盖而改写为参考。
    """
    alleles = parse_gt_alleles(gt)
    if all(a is None for a in alleles):
        return UNKNOWN, "unknown_fully_missing"
    if len(alleles) != 2:
        return UNKNOWN, "unknown_non_diploid"
    if any(a is None for a in alleles):
        return UNKNOWN, "unknown_half_missing"
    if any(a not in (0, 1) for a in alleles):
        return UNKNOWN, "unknown_allele_out_of_range"
    if dp is None:
        return UNKNOWN, "unknown_dp_missing"
    if dp < min_dp:
        return UNKNOWN, "unknown_low_dp"
    if 1 in alleles:
        return CARRIER, "carrier_valid_alt_call"
    return REFERENCE, "reference_valid_call"


def decide_locus(record, covered: bool, min_dp: int) -> tuple[str, str]:
    """单个 (样本, 位点) 判定。``record`` 为 (gt, dp) 或 None（无记录）。"""
    if record is not None:
        gt, dp = record
        return classify_call(gt, dp, min_dp)
    if covered:
        return REFERENCE, "reference_imputed_by_coverage"
    return UNKNOWN, "unknown_not_covered"


def aggregate_gene_call(decisions) -> str:
    """人 × 基因聚合：任一携带即携带；无携带且全参考才未携带。"""
    decisions = list(decisions)
    if any(d == CARRIER for d in decisions):
        return CARRIER
    if decisions and all(d == REFERENCE for d in decisions):
        return NON_CARRIER
    return EXCLUDED


# --------------------------------------------------------------------------- #
# 统计：双侧 Fisher 精确检验 / 交叉乘积优势比 / BH 校正（全部精确分数）
# --------------------------------------------------------------------------- #
def _hypergeom_pmf(x: int, r1: int, r2: int, c1: int) -> Fraction:
    """固定行边际 r1/r2、列边际 c1 时，左上格为 x 的超几何概率。"""
    return Fraction(comb(r1, x) * comb(r2, c1 - x), comb(r1 + r2, c1))


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> Fraction:
    """双侧 Fisher 精确检验。

    表布局：``[[a, b], [c, d]]``（行 = 病例/对照，列 = 携带/未携带）。
    双侧概率 = 固定行列边际下，概率不大于观察表概率的所有表的概率和。
    """
    r1, r2 = a + b, c + d
    c1 = a + c
    lo = max(0, c1 - r2)
    hi = min(r1, c1)
    p_obs = _hypergeom_pmf(a, r1, r2, c1)
    total = Fraction(0)
    for x in range(lo, hi + 1):
        p = _hypergeom_pmf(x, r1, r2, c1)
        if p <= p_obs:
            total += p
    return min(total, Fraction(1))  # 精确分数本不超 1，min 仅作防御


def odds_ratio(a: int, b: int, c: int, d: int):
    """交叉乘积优势比 ad/bc。

    零分母时：分子非零 → ``"inf"``；分子也为零 → ``"undefined"``。
    不使用任何伪计数。
    """
    num, den = a * d, b * c
    if den == 0:
        return "inf" if num > 0 else "undefined"
    return Fraction(num, den)


def benjamini_hochberg(p_values) -> list[Fraction]:
    """Benjamini–Hochberg 校正；返回与输入等长、顺序对应的校正值。

    只对实际做过检验的基因调用（调用方保证 m = len(p_values) 就是
    检验总数，未检验项不得混入）。相同原始概率得到相同校正值。
    """
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    adjusted = [Fraction(1)] * m
    best = Fraction(1)
    for rank in range(m, 0, -1):  # 从最大 p 向最小回推，取后缀最小值
        i = order[rank - 1]
        val = p_values[i] * m / rank
        if val < best:
            best = val
        adjusted[i] = min(best, Fraction(1))
    return adjusted


# --------------------------------------------------------------------------- #
# 编排
# --------------------------------------------------------------------------- #
def run_assoc(dataset_dir: str, labels_path: str, genes_path: str,
              callable_path: str, min_dp: int) -> dict:
    """执行基因级携带者关联分析，返回结构化报告（p 值为精确分数）。"""
    if min_dp < 0:
        raise AssocInputError(f"--min-dp 必须为非负整数，得到 {min_dp}")
    ds = load_dataset_calls(dataset_dir)
    labels = load_labels(labels_path, ds.samples)
    genes = load_genes_bed(genes_path)
    callable_map = load_callable_bed(callable_path, ds.samples)

    case_samples = sorted(s for s, g in labels.items() if g == GROUP_CASE)
    control_samples = sorted(
        s for s, g in labels.items() if g == GROUP_CONTROL)
    analysis_samples = ([(s, GROUP_CASE) for s in case_samples]
                        + [(s, GROUP_CONTROL) for s in control_samples])

    loci_by_chrom: dict[str, list[Locus]] = {}
    for locus in ds.loci:
        loci_by_chrom.setdefault(locus.chrom, []).append(locus)

    tested: list[dict] = []
    not_tested: list[dict] = []

    for gene in sorted(genes):
        # 候选位点：参考占用区间与基因任一（已并集化）区间相交
        candidates: list[Locus] = []
        for chrom, intervals in sorted(genes[gene].items()):
            for locus in loci_by_chrom.get(chrom, ()):
                vs, ve = locus.ref_interval
                if any(vs < e and s < ve for s, e in intervals):
                    candidates.append(locus)
        candidates.sort()

        counts = {g: {CARRIER: 0, NON_CARRIER: 0, EXCLUDED: 0}
                  for g in GROUPS}
        decisions: list[dict] = []
        for sample, group in analysis_samples:
            if not candidates:
                counts[group][EXCLUDED] += 1
                continue
            sample_callable = callable_map.get(sample, {})
            locus_rows = []
            for locus in candidates:
                record = ds.calls.get((sample, locus))
                vs, ve = locus.ref_interval
                covered = interval_covered(
                    sample_callable.get(locus.chrom, ()), vs, ve)
                decision, reason = decide_locus(record, covered, min_dp)
                row = {"chrom": locus.chrom, "pos": locus.pos,
                       "ref": locus.ref, "alt": locus.alt,
                       "ref_start0": vs, "ref_end0": ve,
                       "decision": decision, "reason": reason}
                if record is not None:
                    row["gt"] = record[0]
                    row["dp"] = record[1]
                locus_rows.append(row)
            gene_call = aggregate_gene_call(
                r["decision"] for r in locus_rows)
            counts[group][gene_call] += 1
            decisions.append({"sample": sample, "group": group,
                              "gene_call": gene_call, "loci": locus_rows})

        def _denom(g: str) -> int:
            return counts[g][CARRIER] + counts[g][NON_CARRIER]

        entry = {
            "gene": gene,
            "candidate_loci": [
                {"chrom": l.chrom, "pos": l.pos, "ref": l.ref, "alt": l.alt,
                 "ref_start0": l.ref_interval[0],
                 "ref_end0": l.ref_interval[1]}
                for l in candidates],
            "counts": {g: {**counts[g], "denominator": _denom(g)}
                       for g in GROUPS},
            "decisions": decisions,
        }
        if not candidates:
            entry["status"] = "no_candidate_loci"
            entry["reason"] = "数据集中没有与该基因区间相交的候选位点"
            not_tested.append(entry)
            continue

        a = counts[GROUP_CASE][CARRIER]
        b = counts[GROUP_CASE][NON_CARRIER]
        c = counts[GROUP_CONTROL][CARRIER]
        d = counts[GROUP_CONTROL][NON_CARRIER]
        entry["table"] = {"case_carrier": a, "case_non_carrier": b,
                          "control_carrier": c, "control_non_carrier": d}
        if _denom(GROUP_CASE) == 0 or _denom(GROUP_CONTROL) == 0:
            zero = [g for g in GROUPS if _denom(g) == 0]
            entry["status"] = "zero_denominator"
            entry["reason"] = ("有效分母为零的组: " + ", ".join(zero)
                               + "；该基因不参与检验")
            not_tested.append(entry)
            continue

        entry["status"] = "tested"
        entry["fisher_p_two_sided"] = fisher_exact_two_sided(a, b, c, d)
        entry["odds_ratio"] = odds_ratio(a, b, c, d)
        tested.append(entry)

    # BH 校正只覆盖实际检验的基因；未检验项不进入检验数 m
    adjusted = benjamini_hochberg([e["fisher_p_two_sided"] for e in tested])
    for entry, q in zip(tested, adjusted):
        entry["bh_adjusted_p"] = q
    tested.sort(key=lambda e: (e["bh_adjusted_p"],
                               e["fisher_p_two_sided"], e["gene"]))
    for rank, entry in enumerate(tested, start=1):
        entry["rank"] = rank
    not_tested.sort(key=lambda e: e["gene"])

    return {
        "dataset_id": ds.dataset_id,
        "dataset_dir": dataset_dir,
        "min_dp": min_dp,
        "analysis_samples": {GROUP_CASE: case_samples,
                             GROUP_CONTROL: control_samples},
        "genes_total": len(genes),
        "genes_tested": len(tested),
        "results": tested,
        "not_tested": not_tested,
        "disclaimer": DISCLAIMER,
    }


# --------------------------------------------------------------------------- #
# 输出
# --------------------------------------------------------------------------- #
def to_jsonable(obj):
    """把报告中的 Fraction 递归转成 float，其余原样保留。"""
    if isinstance(obj, Fraction):
        return float(obj)
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    return obj


def _fmt_value(x) -> str:
    if x is None:
        return ""
    if isinstance(x, Fraction):
        return repr(float(x))
    return str(x)


def write_report_files(report: dict, out_dir: str) -> list[str]:
    """写出 assoc_report.json + 汇总/明细 TSV；返回写出的文件路径。"""
    os.makedirs(out_dir, exist_ok=True)
    written = []

    report_path = os.path.join(out_dir, "assoc_report.json")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(to_jsonable(report), indent=2,
                            ensure_ascii=False) + "\n")
    written.append(report_path)

    summary_path = os.path.join(out_dir, "assoc_summary.tsv")
    with open(summary_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene", "status",
            "case_carrier", "case_non_carrier", "case_excluded",
            "case_denominator",
            "control_carrier", "control_non_carrier", "control_excluded",
            "control_denominator",
            "fisher_p_two_sided", "odds_ratio", "bh_adjusted_p", "reason",
        ])
        for e in report["results"] + report["not_tested"]:
            cc, tc = e["counts"][GROUP_CASE], e["counts"][GROUP_CONTROL]
            writer.writerow([
                e["gene"], e["status"],
                cc[CARRIER], cc[NON_CARRIER], cc[EXCLUDED], cc["denominator"],
                tc[CARRIER], tc[NON_CARRIER], tc[EXCLUDED], tc["denominator"],
                _fmt_value(e.get("fisher_p_two_sided")),
                _fmt_value(e.get("odds_ratio")),
                _fmt_value(e.get("bh_adjusted_p")),
                e.get("reason", ""),
            ])
    written.append(summary_path)

    decisions_path = os.path.join(out_dir, "assoc_decisions.tsv")
    with open(decisions_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene", "sample", "group", "gene_call",
            "chrom", "pos", "ref", "alt", "ref_start0", "ref_end0",
            "locus_decision", "reason", "gt", "dp",
        ])
        for e in report["results"] + report["not_tested"]:
            for dec in e["decisions"]:
                for locus in dec["loci"]:
                    writer.writerow([
                        e["gene"], dec["sample"], dec["group"],
                        dec["gene_call"], locus["chrom"], locus["pos"],
                        locus["ref"], locus["alt"], locus["ref_start0"],
                        locus["ref_end0"], locus["decision"],
                        locus["reason"], locus.get("gt", ""),
                        "" if locus.get("dp") is None else locus["dp"],
                    ])
    written.append(decisions_path)
    return written
