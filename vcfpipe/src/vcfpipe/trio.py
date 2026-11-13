"""家系三人组遗传一致性筛查（family trio Mendelian consistency screen）。

定位
----
读取联合分型的**多样本 VCF**（普通或 gzip，按魔数识别，与文件名无关）
与一份四列 TSV 家系表（``family/child/father/mother``），逐位点、逐
三人组解释二倍体遗传是否一致。本命令**不依赖**规范化管线，不拆分多
等位位点（整条多等位记录作为一个位点判断），也**不改写输入基因型**。
仅用于合成研究数据的质量筛查：不据此声称亲缘关系，也不输出临床诊断。

家系表校验（任何一项不合法都拒绝运行，不丢弃家庭后继续）
--------------------------------------------------------
* 表头恰为 ``family child father mother``，四列均必填；
* 每行三人（child/father/mother）两两不同；
* 同一子女只能登记一行父母：完全重复行、跨 family 的矛盾父母都拒绝，
  family 只是分组名，不能用来掩盖同一人的矛盾；
* 沿父母→子女的边不允许成环，报错给出人物链；
* 涉及的每个人都必须出现在 VCF 样本列；VCF 重复样本列、找不到的人
  都点名报出。

位点与调用校验
--------------
* 常染色体名称由 ``--chroms`` 白名单显式给出；白名单外的 VCF 行只
  计入跳过数，不套二倍体规则（性染色体不参与）；
* 选中行要求列数与样本数吻合、POS 为正整数、REF 与各 ALT 是互不
  重复的非空 A/C/G/T 序列、ALT 至多 4 个；FORMAT 必须含 GT，键不得
  重复，各样本字段数与 FORMAT 一致；
* GT 等位编号必须落在当前 REF/ALT 的合法区间内（编号 2 不会被当成
  另一条位点的编号 1）；坏 GT 语法、越界编号、负 DP、损坏选中行均
  使本次命令失败并给出文件与行号；
* 有效的非二倍体调用（深度足够但不是两个等位位）单列为
  ``unsupported``，绝不截成两个等位再算。

缺失与深度语义
--------------
* DP 缺失或低于 ``--min-dp`` 时，该人该位点的整个调用视为**未知**：
  原 GT 原样保留，同时记录不可信原因（``dp_missing``/``dp_below_min``），
  其已知等位不再偷偷约束补全；
* 深度足够时，``0/.``、``./1`` 保留已知那一半，``.`` 与 ``./.``
  表示完全未知。

判定定义（精确问题）
--------------------
把父亲、母亲、孩子各自缺失的等位值，用该位点**全部合法编号**补全，
得到所有不同的**无序**二倍体基因型（相位 ``|`` 不影响计算，``0/1``
与 ``1/0`` 是同一个补全）。枚举父×母×孩子的全部完整组合，检查孩子
能否由父亲贡献一个等位、母亲贡献一个等位组成：

* 所有组合都能解释       -> ``compatible``
* 一个组合都解释不了     -> ``incompatible``
* 部分能、部分不能       -> ``uncertain``

可解释时附一组父方/母方传递见证；不可解释时不编造传递。处理全程
流式：按输入顺序逐位点写明细，内存不随位点数增长；完整读完输入并
关闭明细/摘要输出后才写成功状态，后半段的坏记录不会留下看似完成的
摘要（先写隐藏 staging，成功后单次 rename 发布）。
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from itertools import product

# 判定类别（稳定标识符，不随意改名）
COMPATIBLE = "compatible"
INCOMPATIBLE = "incompatible"
UNCERTAIN = "uncertain"
UNSUPPORTED = "unsupported"
VERDICTS = (COMPATIBLE, INCOMPATIBLE, UNCERTAIN, UNSUPPORTED)

# 调用可信状态
TRUST_OK = "ok"
TRUST_DP_MISSING = "dp_missing"
TRUST_DP_LOW = "dp_below_min"

_BASE_RE = re.compile(r"[ACGT]+")
_INT_RE = re.compile(r"\d+")
_GZIP_MAGIC = b"\x1f\x8b"

DETAIL_FILE = "trio_details.tsv"
SUMMARY_FILE = "trio_summary.tsv"
STATUS_FILE = "trio_status.json"
TOTAL_ROW = "*TOTAL*"

DETAIL_HEADER = [
    "line_no", "chrom", "pos", "ref", "alt",
    "family", "child", "father", "mother",
    "father_gt", "mother_gt", "child_gt",
    "father_dp", "mother_dp", "child_dp",
    "father_trust", "mother_trust", "child_trust",
    "father_candidates", "mother_candidates", "child_candidates",
    "total_combinations", "explainable_combinations",
    "verdict", "unsupported_members", "witness",
]
SUMMARY_HEADER = [
    "family", "child", "father", "mother",
    COMPATIBLE, INCOMPATIBLE, UNCERTAIN, UNSUPPORTED, "total",
]


class TrioError(Exception):
    """家系筛查错误基类。"""


class TrioConfigError(TrioError):
    """配置/关系/表头类错误（CLI 退出码 2）。"""


class TrioRecordError(TrioError):
    """VCF 记录或流损坏类错误（CLI 退出码 1）。"""


# --------------------------------------------------------------------------- #
# 家系表
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Trio:
    """一行登记的三人组。"""

    family: str
    child: str
    father: str
    mother: str
    line_no: int

    def persons(self) -> tuple[str, str, str]:
        return self.child, self.father, self.mother


def load_families(path: str) -> tuple[Trio, ...]:
    """解析并做不依赖 VCF 的结构性校验，错误一次性列全后拒绝。"""
    if not os.path.isfile(path):
        raise TrioConfigError(f"家系表不存在: {path}")

    errors: list[str] = []
    trios: list[Trio] = []
    try:
        with open(path, encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except OSError as exc:
        raise TrioConfigError(f"家系表无法读取 {path}: {exc}") from exc

    if not rows:
        raise TrioConfigError(f"家系表 {path} 为空：缺少表头与登记行")
    header = [c.strip() for c in rows[0]]
    if header != ["family", "child", "father", "mother"]:
        raise TrioConfigError(
            f"家系表 {path} 第 1 行: 表头必须恰为 "
            f"family<TAB>child<TAB>father<TAB>mother，实际为 {rows[0]!r}")

    # child -> (father, mother, 首次登记行号, family)
    registered: dict[str, tuple[str, str, int, str]] = {}
    seen_rows: dict[tuple[str, str, str, str], int] = {}

    for idx, raw in enumerate(rows[1:], start=2):
        if not raw or all(c.strip() == "" for c in raw):
            continue  # 允许表尾空行
        if len(raw) != 4:
            errors.append(f"家系表 {path} 第 {idx} 行: 必须恰为 4 列，"
                          f"实际 {len(raw)} 列")
            continue
        family, child, father, mother = (c.strip() for c in raw)
        if not all((family, child, father, mother)):
            errors.append(f"家系表 {path} 第 {idx} 行: 四列均必填，"
                          f"不得有空值")
            continue
        if len({child, father, mother}) != 3:
            errors.append(f"家系表 {path} 第 {idx} 行: child/father/mother "
                          f"必须两两不同（{child}/{father}/{mother}）")
            continue
        key = (family, child, father, mother)
        if key in seen_rows:
            errors.append(f"家系表 {path} 第 {idx} 行: 与第 "
                          f"{seen_rows[key]} 行完全重复")
            continue
        seen_rows[key] = idx
        if child in registered:
            pf, pm, pline, pfamily = registered[child]
            if (pf, pm) == (father, mother):
                errors.append(
                    f"家系表 {path} 第 {idx} 行: 子女 {child} 已在第 "
                    f"{pline} 行登记过同一对父母（family {pfamily}），"
                    f"同一子女只能登记一行")
            else:
                errors.append(
                    f"家系表 {path} 第 {idx} 行: 子女 {child} 的父母与第 "
                    f"{pline} 行矛盾（第 {pline} 行: {pf}/{pm}，本行: "
                    f"{father}/{mother}）；family 分组名不能掩盖同一人的"
                    f"矛盾父母")
            continue
        registered[child] = (father, mother, idx, family)
        trios.append(Trio(family, child, father, mother, idx))

    if errors:
        raise TrioConfigError("\n".join(errors))

    cycle = _detect_cycle(trios)
    if cycle is not None:
        raise TrioConfigError(
            f"家系表 {path}: 父母→子女关系存在环: {' -> '.join(cycle)}")
    if not trios:
        raise TrioConfigError(f"家系表 {path} 没有任何有效登记行")
    return tuple(trios)


def _detect_cycle(trios: tuple[Trio, ...]) -> list[str] | None:
    """DFS 检测 parent->child 有向图中的环，返回环上的人物链。"""
    adjacency: dict[str, set[str]] = {}
    nodes: set[str] = set()
    for t in trios:
        nodes.update(t.persons())
        adjacency.setdefault(t.father, set()).add(t.child)
        adjacency.setdefault(t.mother, set()).add(t.child)
        for p in t.persons():
            adjacency.setdefault(p, set())

    WHITE, GRAY, BLACK = 0, 1, 2
    color = {n: WHITE for n in nodes}
    stack: list[str] = []

    def dfs(node: str) -> list[str] | None:
        color[node] = GRAY
        stack.append(node)
        for nxt in sorted(adjacency[node]):
            if color[nxt] == GRAY:
                start = stack.index(nxt)
                return stack[start:] + [nxt]
            if color[nxt] == WHITE:
                found = dfs(nxt)
                if found is not None:
                    return found
        stack.pop()
        color[node] = BLACK
        return None

    for n in sorted(nodes):
        if color[n] == WHITE:
            found = dfs(n)
            if found is not None:
                return found
    return None


# --------------------------------------------------------------------------- #
# VCF 流式读取
# --------------------------------------------------------------------------- #
def _open_text(path: str):
    """按魔数打开普通或 gzip 文本文件（与文件名后缀无关）。"""
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
    except OSError as exc:
        raise TrioConfigError(f"VCF 无法打开 {path}: {exc}") from exc
    if magic == _GZIP_MAGIC:
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def read_vcf_header(path: str) -> tuple[str, ...]:
    """读取 #CHROM 头，返回样本列；重复样本列与缺样本由调用方校验。"""
    with _open_text(path) as fh:
        line_no = 0
        for raw in fh:
            line_no += 1
            line = raw.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) < 8:
                    raise TrioConfigError(
                        f"{path}:{line_no} #CHROM 头列数不足 8")
                if len(cols) == 8:
                    raise TrioConfigError(
                        f"{path}:{line_no} VCF 没有样本列，无法做家系筛查")
                if len(cols) == 9:
                    raise TrioConfigError(
                        f"{path}:{line_no} 有 FORMAT 列但没有任何样本列")
                samples = tuple(cols[9:])
                dupes = sorted({s for s in samples if samples.count(s) > 1})
                if dupes:
                    raise TrioConfigError(
                        f"{path}:{line_no} VCF 样本列存在重复: "
                        f"{', '.join(dupes)}")
                return samples
            raise TrioRecordError(
                f"{path}:{line_no} #CHROM 之前出现非法头行: {line!r}")
    raise TrioConfigError(f"{path}: 文件中找不到 #CHROM 头行")


def check_persons(trios: tuple[Trio, ...], samples: tuple[str, ...],
                  vcf_path: str) -> None:
    """家系中的每个人都必须是 VCF 样本列之一，错误一次性列全。"""
    sample_set = set(samples)
    missing: list[str] = []
    seen: set[str] = set()
    for t in trios:
        for person in t.persons():
            if person not in sample_set and person not in seen:
                seen.add(person)
                missing.append(person)
    if missing:
        raise TrioConfigError(
            f"家系表中的人物不在 VCF {vcf_path} 样本列中: "
            f"{', '.join(missing)}；VCF 样本列为 {', '.join(samples)}")


# --------------------------------------------------------------------------- #
# 选中行与调用解析
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MemberCall:
    """三人组中一个成员在一个位点的解析结果。"""

    role: str               # father / mother / child
    person: str
    raw_gt: str             # 原样保留的 GT 文本
    dp: int | None          # 原样保留的 DP（None = 缺失）
    trusted: bool           # DP 存在且 >= 阈值
    trust_reason: str | None
    alleles: tuple[int | None, ...] | None  # None = 完全未知（含 "."）
    phased: bool


def _parse_gt(token: str, n_alts: int, path: str, line_no: int,
              person: str) -> tuple[tuple[int | None, ...] | None, bool]:
    """解析 GT 并做语法/越界校验。返回 (等位元组或 None, 是否相位化)。"""
    if token == "":
        raise TrioRecordError(
            f"{path}:{line_no} 样本 {person} 的 GT 为空")
    if token == ".":
        return None, False
    if "|" in token and "/" in token:
        raise TrioRecordError(
            f"{path}:{line_no} 样本 {person} 的 GT {token!r} 混用了 "
            f"相位分隔符 | 与 /")
    sep = "|" if "|" in token else "/"
    parts = token.split(sep)
    alleles: list[int | None] = []
    for part in parts:
        if part == "":
            raise TrioRecordError(
                f"{path}:{line_no} 样本 {person} 的 GT {token!r} 含空位")
        if part == ".":
            alleles.append(None)
            continue
        if not _INT_RE.fullmatch(part):
            raise TrioRecordError(
                f"{path}:{line_no} 样本 {person} 的 GT {token!r} 含非法"
                f"等位编号 {part!r}")
        value = int(part)
        if value > n_alts:
            raise TrioRecordError(
                f"{path}:{line_no} 样本 {person} 的 GT 等位编号 {value} "
                f"越界：该位点 REF+ALT 合法编号为 0..{n_alts}")
        alleles.append(value)
    if all(a is None for a in alleles):
        return None, (sep == "|")
    return tuple(alleles), (sep == "|")


def _parse_dp(token: str, path: str, line_no: int,
              person: str) -> int | None:
    if token == "":
        raise TrioRecordError(
            f"{path}:{line_no} 样本 {person} 的 DP 字段为空")
    if token == ".":
        return None
    if not _INT_RE.fullmatch(token):
        raise TrioRecordError(
            f"{path}:{line_no} 样本 {person} 的 DP {token!r} 非法："
            f"必须是非负整数或 .（负深度/非数字均不接受）")
    return int(token)


def _validate_locus(cols: list[str], samples: tuple[str, ...],
                    path: str, line_no: int) -> tuple[int, str,
                                                      tuple[str, ...],
                                                      dict[str, dict[str, str]],
                                                      list[str]]:
    """校验一个白名单内位点，返回 (pos, ref, alts, 每样本字段表, FORMAT)。"""
    n_samples = len(samples)
    if len(cols) != 9 + n_samples:
        raise TrioRecordError(
            f"{path}:{line_no} 列数 {len(cols)} 与 8 固定列 + FORMAT + "
            f"{n_samples} 个样本列（应共 {9 + n_samples} 列）不符")
    chrom, pos_s, _vid, ref, alt_s = cols[0], cols[1], cols[2], cols[3], cols[4]
    if not _INT_RE.fullmatch(pos_s) or int(pos_s) <= 0:
        raise TrioRecordError(
            f"{path}:{line_no} POS {pos_s!r} 不是正整数")
    pos = int(pos_s)
    if not _BASE_RE.fullmatch(ref):
        raise TrioRecordError(
            f"{path}:{line_no} REF {ref!r} 不是非空 A/C/G/T 序列")
    if alt_s == ".":
        alts: tuple[str, ...] = ()
    else:
        alts = tuple(alt_s.split(","))
        if len(alts) > 4:
            raise TrioRecordError(
                f"{path}:{line_no} ALT 数量 {len(alts)} 超过上限 4（多等位"
                f"位点整体判断，但最多接受 4 个 ALT）")
        if any(a == "" for a in alts):
            raise TrioRecordError(
                f"{path}:{line_no} ALT {alt_s!r} 含空等位")
        if any(not _BASE_RE.fullmatch(a) for a in alts):
            raise TrioRecordError(
                f"{path}:{line_no} ALT {alt_s!r} 含非 A/C/G/T 序列")
        if len(set(alts)) != len(alts) or ref in alts:
            raise TrioRecordError(
                f"{path}:{line_no} REF {ref!r} 与各 ALT {alts} 必须互不"
                f"重复")
    format_keys = cols[8].split(":")
    if "" in format_keys:
        raise TrioRecordError(
            f"{path}:{line_no} FORMAT {cols[8]!r} 含空键")
    if len(set(format_keys)) != len(format_keys):
        raise TrioRecordError(
            f"{path}:{line_no} FORMAT {cols[8]!r} 含重复键")
    if "GT" not in format_keys:
        raise TrioRecordError(
            f"{path}:{line_no} FORMAT {cols[8]!r} 缺少必需的 GT 键")
    per_sample: dict[str, dict[str, str]] = {}
    for person, blob in zip(samples, cols[9:]):
        fields = blob.split(":")
        if len(fields) != len(format_keys):
            raise TrioRecordError(
                f"{path}:{line_no} 样本 {person} 字段数 {len(fields)} 与 "
                f"FORMAT 键数 {len(format_keys)} 不符")
        per_sample[person] = dict(zip(format_keys, fields))
    return pos, ref, alts, per_sample, format_keys


# --------------------------------------------------------------------------- #
# 补全与传递判定（产品实现；独立核对器不得调用本节函数）
# --------------------------------------------------------------------------- #
def diploid_candidates(alleles: tuple[int | None, ...],
                       n_alts: int) -> tuple[tuple[int, ...], ...]:
    """把缺失位用全部合法编号补全，返回排序去重后的无序二倍体基因型。"""
    options = [range(n_alts + 1) if a is None else (a,) for a in alleles]
    out = {tuple(sorted((x, y))) for x, y in product(*options)}
    return tuple(sorted(out))


def full_unknown_candidates(n_alts: int) -> tuple[tuple[int, ...], ...]:
    """完全未知（含低深度整调用未知）的全量无序二倍体补全。"""
    return diploid_candidates((None, None), n_alts)


def _format_gt(alleles: tuple[int, ...]) -> str:
    return "/".join(str(a) for a in alleles)


def transmission_witness(father_gt: tuple[int, ...],
                         mother_gt: tuple[int, ...],
                         child_gt: tuple[int, ...]
                         ) -> tuple[int, int] | None:
    """孩子能否由父亲一个等位 + 母亲一个等位组成；返回 (父等位, 母等位)。"""
    for paternal in father_gt:
        for maternal in mother_gt:
            if tuple(sorted((paternal, maternal))) == child_gt:
                return paternal, maternal
    return None


@dataclass(frozen=True)
class TrioVerdict:
    verdict: str
    candidates: dict[str, tuple[tuple[int, ...], ...]]
    total: int
    explainable: int
    unsupported_members: tuple[str, ...]
    witness: str | None


def evaluate_trio(calls: dict[str, MemberCall], n_alts: int) -> TrioVerdict:
    """对一个三人组在一个位点做补全枚举与三类判定。"""
    roles = ("father", "mother", "child")
    unsupported = tuple(
        c.role for c in (calls[r] for r in roles)
        if c.trusted and c.alleles is not None and len(c.alleles) != 2
    )
    if unsupported:
        return TrioVerdict(
            verdict=UNSUPPORTED,
            candidates={r: () for r in roles},
            total=0, explainable=0,
            unsupported_members=unsupported, witness=None)

    candidates: dict[str, tuple[tuple[int, ...], ...]] = {}
    for role in roles:
        call = calls[role]
        if not call.trusted or call.alleles is None:
            # 整调用未知：低深度者的已知等位不参与，全量补全
            candidates[role] = full_unknown_candidates(n_alts)
        else:
            candidates[role] = diploid_candidates(call.alleles, n_alts)

    total = 0
    explainable = 0
    witness: str | None = None
    for fg, mg, cg in product(candidates["father"],
                              candidates["mother"],
                              candidates["child"]):
        total += 1
        got = transmission_witness(fg, mg, cg)
        if got is not None:
            explainable += 1
            if witness is None:
                paternal, maternal = got
                witness = (
                    f"{_format_gt(fg)}|{_format_gt(mg)}|"
                    f"{_format_gt(cg)};paternal={paternal};maternal={maternal}"
                )

    if explainable == 0:
        verdict = INCOMPATIBLE
    elif explainable == total:
        verdict = COMPATIBLE
    else:
        verdict = UNCERTAIN
    return TrioVerdict(
        verdict=verdict, candidates=candidates, total=total,
        explainable=explainable, unsupported_members=(), witness=witness)


# --------------------------------------------------------------------------- #
# 端到端流式筛查
# --------------------------------------------------------------------------- #
def _prepare_out_dir(out_dir: str) -> str:
    if os.path.lexists(out_dir):
        if not os.path.isdir(out_dir):
            raise TrioConfigError(f"目标 {out_dir} 已存在且不是目录")
        if os.listdir(out_dir):
            raise TrioConfigError(
                f"目标目录 {out_dir} 非空：筛查只写入全新空目录")
    parent = os.path.dirname(os.path.abspath(out_dir))
    os.makedirs(parent, exist_ok=True)
    return tempfile.mkdtemp(
        prefix=f".{os.path.basename(out_dir)}.staging-", dir=parent)


def _build_member(role: str, person: str, fields: dict[str, str],
                  min_dp: int, n_alts: int, path: str,
                  line_no: int) -> MemberCall:
    raw_gt = fields["GT"]
    alleles, phased = _parse_gt(raw_gt, n_alts, path, line_no, person)
    dp = _parse_dp(fields["DP"], path, line_no, person) \
        if "DP" in fields else None
    if dp is None:
        trusted, reason = False, TRUST_DP_MISSING
    elif dp < min_dp:
        trusted, reason = False, TRUST_DP_LOW
    else:
        trusted, reason = True, TRUST_OK
    return MemberCall(
        role=role, person=person, raw_gt=raw_gt, dp=dp,
        trusted=trusted, trust_reason=reason,
        alleles=alleles, phased=phased)


def run_trio(vcf_path: str, family_path: str, chroms: tuple[str, ...],
             min_dp: int, out_dir: str) -> dict:
    """执行筛查并原子发布输出，返回可 JSON 化的运行报告。"""
    if not chroms:
        raise TrioConfigError("--chroms 至少需要一个常染色体名称")
    if len(set(chroms)) != len(chroms):
        raise TrioConfigError(f"--chroms 存在重复名称: {chroms}")
    if min_dp < 1:
        raise TrioConfigError(f"--min-dp 必须是正整数，实际为 {min_dp}")
    if not os.path.isfile(vcf_path):
        raise TrioConfigError(f"VCF 不存在: {vcf_path}")

    trios = load_families(family_path)
    samples = read_vcf_header(vcf_path)
    check_persons(trios, samples, vcf_path)

    chrom_set = set(chroms)
    staging = _prepare_out_dir(out_dir)
    detail_path = os.path.join(staging, DETAIL_FILE)
    summary_path = os.path.join(staging, SUMMARY_FILE)

    # 每个三人组一个计数器；内存只与家庭数有关，不随位点数增长
    counters = [dict.fromkeys(VERDICTS, 0) for _ in trios]
    totals = dict.fromkeys(VERDICTS, 0)
    body_lines = 0
    selected_lines = 0
    skipped_lines = 0

    detail_fh = open(detail_path, "w", encoding="utf-8", newline="")
    summary_fh = open(summary_path, "w", encoding="utf-8", newline="")
    published = False
    try:
        detail_writer = csv.writer(detail_fh, delimiter="\t",
                                   lineterminator="\n")
        detail_writer.writerow(DETAIL_HEADER)
        summary_writer = csv.writer(summary_fh, delimiter="\t",
                                    lineterminator="\n")
        summary_writer.writerow(SUMMARY_HEADER)

        with _open_text(vcf_path) as fh:
            line_no = 0
            in_body = False
            for raw in fh:
                line_no += 1
                line = raw.rstrip("\n\r")
                if line == "":
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    in_body = True
                    continue
                if not in_body:
                    raise TrioRecordError(
                        f"{vcf_path}:{line_no} #CHROM 之前出现非法行: "
                        f"{line!r}")
                body_lines += 1
                cols = line.split("\t")
                if cols[0] not in chrom_set:
                    skipped_lines += 1
                    continue

                selected_lines += 1
                pos, ref, alts, per_sample, _fmt = _validate_locus(
                    cols, samples, vcf_path, line_no)
                n_alts = len(alts)
                alt_text = ",".join(alts) if alts else "."

                for idx, trio in enumerate(trios):
                    calls = {
                        "father": _build_member(
                            "father", trio.father,
                            per_sample[trio.father], min_dp, n_alts,
                            vcf_path, line_no),
                        "mother": _build_member(
                            "mother", trio.mother,
                            per_sample[trio.mother], min_dp, n_alts,
                            vcf_path, line_no),
                        "child": _build_member(
                            "child", trio.child,
                            per_sample[trio.child], min_dp, n_alts,
                            vcf_path, line_no),
                    }
                    result = evaluate_trio(calls, n_alts)
                    counters[idx][result.verdict] += 1
                    totals[result.verdict] += 1
                    row = [
                        line_no, cols[0], pos, ref, alt_text,
                        trio.family, trio.child, trio.father, trio.mother,
                        calls["father"].raw_gt,
                        calls["mother"].raw_gt,
                        calls["child"].raw_gt,
                        "." if calls["father"].dp is None
                        else calls["father"].dp,
                        "." if calls["mother"].dp is None
                        else calls["mother"].dp,
                        "." if calls["child"].dp is None
                        else calls["child"].dp,
                        calls["father"].trust_reason,
                        calls["mother"].trust_reason,
                        calls["child"].trust_reason,
                        ",".join(_format_gt(g)
                                 for g in result.candidates["father"]),
                        ",".join(_format_gt(g)
                                 for g in result.candidates["mother"]),
                        ",".join(_format_gt(g)
                                 for g in result.candidates["child"]),
                        result.total, result.explainable,
                        result.verdict,
                        ",".join(result.unsupported_members),
                        result.witness or "",
                    ]
                    detail_writer.writerow(row)

        # 输入完整读完、明细落盘后才写摘要
        detail_fh.close()
        for trio, counts in zip(trios, counters):
            row_total = sum(counts.values())
            summary_writer.writerow([
                trio.family, trio.child, trio.father, trio.mother,
                counts[COMPATIBLE], counts[INCOMPATIBLE],
                counts[UNCERTAIN], counts[UNSUPPORTED], row_total,
            ])
        summary_writer.writerow([
            TOTAL_ROW, "", "", "",
            totals[COMPATIBLE], totals[INCOMPATIBLE],
            totals[UNCERTAIN], totals[UNSUPPORTED],
            sum(totals.values()),
        ])
        summary_fh.close()

        # 输入完整读完、明细与摘要均已关闭后，才构造成功状态并写入
        # staging；后半段任何坏记录都会在到达此处前抛出，绝不留下看似
        # 完成的摘要。三个文件齐备后单次 rename 原子发布。
        report = {
            "status": "success",
            "vcf_path": vcf_path,
            "family_path": family_path,
            "chrom_whitelist": list(chroms),
            "min_dp": min_dp,
            "vcf_body_lines": body_lines,
            "selected_lines": selected_lines,
            "skipped_non_whitelist": skipped_lines,
            "trio_count": len(trios),
            "detail_rows": sum(totals.values()),
            "totals": totals,
            "summary_by_trio": [
                {
                    "family": t.family, "child": t.child,
                    "father": t.father, "mother": t.mother,
                    **counts,
                    "total": sum(counts.values()),
                }
                for t, counts in zip(trios, counters)
            ],
            "outputs": {
                "details": os.path.join(out_dir, DETAIL_FILE),
                "summary": os.path.join(out_dir, SUMMARY_FILE),
                "status": os.path.join(out_dir, STATUS_FILE),
            },
            "disclaimer": "合成研究数据的质量筛查，不改写输入基因型，"
                          "不据此声称亲缘关系或临床诊断。",
        }
        with open(os.path.join(staging, STATUS_FILE), "w",
                  encoding="utf-8") as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2,
                      sort_keys=True)
            fh.write("\n")
        os.replace(staging, out_dir)
        published = True
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise TrioRecordError(f"读取/写出失败 {vcf_path}: {exc}") from exc
    finally:
        # 失败时关闭句柄并清理 staging，绝不留下看似完成的输出；
        # 成功发布后 staging 已不存在，rmtree 为空操作
        for fh_ in (detail_fh, summary_fh):
            if not fh_.closed:
                fh_.close()
        if not published:
            shutil.rmtree(staging, ignore_errors=True)
    return report
