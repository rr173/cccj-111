"""分区重算规划（replan）：内容指纹账本、依赖闭包规划与完成确认。

定位
----
实验室每天接入一批数据：清洗表按天生成，周滚动表读前几天，月表要等整月
的清洗分区。补交一天数据后，需要回答三个问题——指定目标分区**还需算
什么**、**眼下能开工什么**、**缺哪份数据**。本模块只做这件规划的事：
实际清洗/分析程序由外部执行，本模块不启动作业，也不推断测序结果。

概念
----
* **节点**：``source``（外部数据，账本只登记其内容指纹）或 ``derived``
  （外部计算产物，账本保存最近一次成功计算的回执）。每个节点有唯一名
  称与粒度 ``day``（YYYY-MM-DD）或 ``month``（YYYY-MM），分区标识严格
  按公历计算。
* **依赖映射**（只有三种）：
  - ``identity``：两端粒度相同，输入分区 = 输出分区；
  - ``day_window``：两端都是 day，给定整数 ``lo<=hi<=0``，取输出日期
    加 lo 到加 hi 的全部日分区（两端包含）；
  - ``month_days``：输出 month、输入 day，取该月每一天。
  展开后相同输入分区只计一次；配置中节点次序、依赖次序不影响含义。
* **账本**（SQLite，只管理内容指纹，不搬运数据文件）：source 分区登记
  64 位十六进制 SHA256 或被明确撤回（撤回 = 缺数据，绝不等价于一份空
  表；摘要相同 = 内容相同，不按登记次数制造新版本）；derived 保存最近
  一次成功计算的输出摘要、当时的 recipe/依赖定义、实际读到的每一个直
  接输入分区与摘要。账本退出重开后仍可继续使用，也允许装载修改过的配
  置——回执逐节点记录自己的规则，无关节点的配置调整不会让全库过期。

状态定义
--------
* source 分区：``present``（已登记）/ ``missing``（未登记或已撤回）。
* derived 分区：
  - ``current``：回执存在、本节点规则未变、直接输入全部有效且摘要与
    成功回执逐项相等；
  - ``ready``：无回执或回执过期，但直接输入已经全部 current——眼下
    可以开工，规划会为它产出输入票据；
  - ``waiting``：还得等待上游派生分区重算；
  - ``blocked``：沿任何依赖最终缺失了 source 分区，缺失链点名那一天。

确认的事务语义
--------------
计算者拿票据与新输出摘要来确认时，程序在**一次数据库事务**内重新检查
票据对应的规则与所有直接输入：它们须仍然有效且与票据相符，才替换该分
区回执。拿到票据后源数据改了、上游失效了或本节点 recipe 改了，都拒收
旧结果、保留原回执（退出码 1），重新规划即可。票据只是"申领时的状态
快照"，不能绕过节点、粒度与依赖完整性检查。

传播以内容为准
--------------
下游回执记录的是上游分区的**输出摘要**。上游重算后交回相同摘要时，下
游输入实际上没变，直接重新成为 current，无须重跑；摘要不同才进入
ready。检查下游时不只看它记录的摘要仍在数据库里，还要求生产该摘要的
上游回执如今依然有效。
"""
from __future__ import annotations

import calendar
import json
import os
import re
import sqlite3
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta, timezone
from urllib.parse import quote

SCHEMA_VERSION = "1"
TICKET_VERSION = 1

DAY_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
MONTH_RE = re.compile(r"^[0-9]{4}-[0-9]{2}$")
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

MAP_IDENTITY = "identity"
MAP_DAY_WINDOW = "day_window"
MAP_MONTH_DAYS = "month_days"
MAP_KINDS = (MAP_IDENTITY, MAP_DAY_WINDOW, MAP_MONTH_DAYS)

KIND_SOURCE = "source"
KIND_DERIVED = "derived"
GRAN_DAY = "day"
GRAN_MONTH = "month"


# --------------------------------------------------------------------------- #
# 错误类型（exit_code 与 CLI 退出码对应）
# --------------------------------------------------------------------------- #
class ReplanError(Exception):
    """重算规划操作被拒绝。"""

    exit_code = 2


class ConfigError(ReplanError):
    """任务依赖 JSON 非法（引用不明、成环、映射与粒度不符等）。"""


class UsageError(ReplanError):
    """用法/状态非法（分区标识、摘要格式、账本缺失、日期超界、票据结构等）。"""


class ConfirmRejected(ReplanError):
    """完成确认被拒绝：票据与当前规则/输入不符，原回执保留。"""

    exit_code = 1

    def __init__(self, message, detail=None):
        super().__init__(message)
        self.detail = detail or {}


# --------------------------------------------------------------------------- #
# 分区标识与摘要
# --------------------------------------------------------------------------- #
def validate_partition(granularity, part):
    """校验分区标识与粒度相符且为真实公历日期；非法时抛 UsageError。"""
    if granularity == GRAN_DAY:
        if not isinstance(part, str) or not DAY_RE.match(part):
            raise UsageError(
                f"分区标识 {part!r} 与 day 粒度不符（需要严格的 YYYY-MM-DD）")
        try:
            date(int(part[0:4]), int(part[5:7]), int(part[8:10]))
        except ValueError:
            raise UsageError(
                f"分区标识 {part!r} 不是公历中真实存在的日期") from None
    else:
        if not isinstance(part, str) or not MONTH_RE.match(part):
            raise UsageError(
                f"分区标识 {part!r} 与 month 粒度不符（需要严格的 YYYY-MM）")
        if not 1 <= int(part[5:7]) <= 12:
            raise UsageError(f"分区标识 {part!r} 的月份超出 01..12")


def normalize_sha256(text, what="SHA256 摘要"):
    """校验 64 位十六进制摘要并归一为小写。"""
    if not isinstance(text, str) or not SHA256_RE.match(text):
        raise UsageError(f"{what}必须是 64 位十六进制字符串，实际为 {text!r}")
    return text.lower()


# --------------------------------------------------------------------------- #
# 任务依赖配置
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Dependency:
    node: str
    map: str
    lo: int | None = None
    hi: int | None = None


@dataclass
class NodeCfg:
    name: str
    kind: str                # KIND_SOURCE | KIND_DERIVED
    granularity: str         # GRAN_DAY | GRAN_MONTH
    recipe: str | None = None
    dependencies: tuple = ()  # tuple[Dependency]，已按规范序排列并去重


@dataclass
class Config:
    nodes: dict              # name -> NodeCfg


_NODE_KEYS = {"name", "type", "granularity", "recipe", "dependencies"}
_DEP_KEYS = {"node", "map", "lo", "hi"}
_TOP_KEYS = {"nodes", "description"}


def load_config_file(path):
    """读取并整体校验任务依赖 JSON；任何非法都解释清楚并拒绝。"""
    try:
        with open(path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)
    except OSError as exc:
        raise UsageError(f"无法读取配置 {path}: {exc}") from None
    except json.JSONDecodeError as exc:
        raise ConfigError(f"配置 {path} 不是合法 JSON: {exc}") from None
    return parse_config(raw, path)


def parse_config(raw, path="<config>"):
    if not isinstance(raw, dict) or not isinstance(raw.get("nodes"), list):
        raise ConfigError(f"配置 {path} 必须是含 nodes 列表的 JSON 对象")
    extra = set(raw) - _TOP_KEYS
    if extra:
        raise ConfigError(f"配置 {path} 含未知顶层字段: {sorted(extra)}")
    nodes = {}
    for i, entry in enumerate(raw["nodes"]):
        node = _parse_node(entry, f"{path} nodes[{i}]")
        if node.name in nodes:
            raise ConfigError(f"配置 {path}: 节点名 {node.name!r} 重复定义")
        nodes[node.name] = node
    if not nodes:
        raise ConfigError(f"配置 {path}: nodes 为空")
    for node in nodes.values():
        for dep in node.dependencies:
            target = nodes.get(dep.node)
            if target is None:
                raise ConfigError(
                    f"节点 {node.name}: 依赖引用了未定义的节点 {dep.node!r}")
            _check_mapping(node, target, dep)
    _check_cycles(nodes)
    return Config(nodes)


def _parse_node(entry, where):
    if not isinstance(entry, dict):
        raise ConfigError(f"{where}: 节点定义必须是对象")
    extra = set(entry) - _NODE_KEYS
    if extra:
        raise ConfigError(f"{where}: 未知字段 {sorted(extra)}")
    name = entry.get("name")
    if not isinstance(name, str) or not NAME_RE.match(name):
        raise ConfigError(
            f"{where}: 节点名 {name!r} 非法"
            "（需匹配 [A-Za-z0-9][A-Za-z0-9._-]*）")
    kind = entry.get("type")
    if kind not in (KIND_SOURCE, KIND_DERIVED):
        raise ConfigError(
            f"节点 {name}: type 必须是 source 或 derived，实际为 {kind!r}")
    gran = entry.get("granularity")
    if gran not in (GRAN_DAY, GRAN_MONTH):
        raise ConfigError(
            f"节点 {name}: granularity 必须是 day 或 month，实际为 {gran!r}")
    if kind == KIND_SOURCE:
        if "recipe" in entry or "dependencies" in entry:
            raise ConfigError(
                f"节点 {name}: source 节点不带 recipe 或 dependencies")
        return NodeCfg(name, kind, gran)
    recipe = entry.get("recipe")
    if not isinstance(recipe, str) or not recipe.strip():
        raise ConfigError(
            f"节点 {name}: derived 节点必须给出非空 recipe 字符串")
    deps_raw = entry.get("dependencies")
    if not isinstance(deps_raw, list) or not deps_raw:
        raise ConfigError(
            f"节点 {name}: derived 节点至少依赖一项"
            "（dependencies 须为非空列表）")
    deps, seen = [], set()
    for j, d in enumerate(deps_raw):
        dep = _parse_dependency(d, f"节点 {name} dependencies[{j}]")
        key = (dep.node, dep.map, dep.lo, dep.hi)
        if key not in seen:  # 完全重复的依赖不影响含义，去重即可
            seen.add(key)
            deps.append(dep)
    deps.sort(key=_dep_sort_key)
    return NodeCfg(name, kind, gran, recipe=recipe, dependencies=tuple(deps))


def _parse_dependency(d, where):
    if not isinstance(d, dict):
        raise ConfigError(f"{where}: 依赖必须是对象")
    extra = set(d) - _DEP_KEYS
    if extra:
        raise ConfigError(f"{where}: 未知字段 {sorted(extra)}")
    node = d.get("node")
    if not isinstance(node, str) or not node:
        raise ConfigError(f"{where}: 依赖缺少节点名 node")
    map_kind = d.get("map")
    if map_kind not in MAP_KINDS:
        raise ConfigError(
            f"{where}: map 必须是 {list(MAP_KINDS)} 之一，实际为 {map_kind!r}")
    if map_kind == MAP_DAY_WINDOW:
        if "lo" not in d or "hi" not in d:
            raise ConfigError(f"{where}: day_window 需要整数 lo 与 hi")
        lo, hi = d["lo"], d["hi"]
        if (isinstance(lo, bool) or isinstance(hi, bool)
                or not isinstance(lo, int) or not isinstance(hi, int)):
            raise ConfigError(f"{where}: day_window 的 lo/hi 必须是整数")
        if not lo <= hi <= 0:
            raise ConfigError(
                f"{where}: day_window 需要 lo<=hi<=0，实际 lo={lo} hi={hi}")
        return Dependency(node, map_kind, lo, hi)
    if "lo" in d or "hi" in d:
        raise ConfigError(f"{where}: {map_kind} 映射不接受 lo/hi")
    return Dependency(node, map_kind)


def _check_mapping(node, target, dep):
    pair = f"节点 {node.name}({node.granularity}) 经 {dep.map} 依赖 " \
           f"{target.name}({target.granularity})"
    if dep.map == MAP_IDENTITY:
        if node.granularity != target.granularity:
            raise ConfigError(f"{pair}: identity 要求两端粒度相同")
    elif dep.map == MAP_DAY_WINDOW:
        if node.granularity != GRAN_DAY or target.granularity != GRAN_DAY:
            raise ConfigError(f"{pair}: day_window 要求两端都是 day")
    else:  # month_days
        if node.granularity != GRAN_MONTH or target.granularity != GRAN_DAY:
            raise ConfigError(f"{pair}: month_days 要求输出 month、输入 day")


def _check_cycles(nodes):
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {n: WHITE for n in nodes}
    stack = []

    def dfs(name):
        color[name] = GRAY
        stack.append(name)
        for dep in nodes[name].dependencies:
            t = dep.node
            if color[t] == GRAY:
                chain = stack[stack.index(t):] + [t]
                raise ConfigError("节点依赖存在环: " + " -> ".join(chain))
            if color[t] == WHITE:
                dfs(t)
        stack.pop()
        color[name] = BLACK

    for n in nodes:
        if color[n] == WHITE:
            dfs(n)


def _dep_sort_key(d):
    return (d.node, d.map, d.lo if d.lo is not None else 0,
            d.hi if d.hi is not None else 0)


def canonical_rule(cfg, node_name):
    """本节点规则的规范形式：recipe + 粒度 + 依赖定义（含输入端粒度）。

    回执里保存的就是它；规划与确认时按字符串/结构逐项比对，配置中节点
    次序、依赖次序不影响结果。
    """
    node = cfg.nodes[node_name]
    deps = [{
        "node": d.node,
        "map": d.map,
        "lo": d.lo,
        "hi": d.hi,
        "input_granularity": cfg.nodes[d.node].granularity,
    } for d in node.dependencies]
    deps.sort(key=lambda x: (x["node"], x["map"],
                             x["lo"] if x["lo"] is not None else 0,
                             x["hi"] if x["hi"] is not None else 0))
    return {
        "granularity": node.granularity,
        "recipe": node.recipe,
        "dependencies": deps,
    }


def canonical_json(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False)


# --------------------------------------------------------------------------- #
# 依赖展开
# --------------------------------------------------------------------------- #
def expand_inputs(cfg, node_name, part):
    """展开 derived 节点某输出分区的全部直接输入。

    返回按 (节点, 分区) 排序、去重后的 [(node, partition)]。日期运算超界
    （如 0001-01-01 再往前取窗口）解释清楚并拒绝。
    """
    node = cfg.nodes[node_name]
    out = set()
    for dep in node.dependencies:
        if dep.map == MAP_IDENTITY:
            out.add((dep.node, part))
        elif dep.map == MAP_DAY_WINDOW:
            base = date(int(part[0:4]), int(part[5:7]), int(part[8:10]))
            for k in range(dep.lo, dep.hi + 1):
                try:
                    d = base + timedelta(days=k)
                except OverflowError:
                    raise UsageError(
                        f"节点 {node_name} 分区 {part} 的 day_window 依赖 "
                        f"{dep.node}: 日期运算超界（{part} 加 {k} 天超出"
                        "公历可表示范围）") from None
                out.add((dep.node, d.isoformat()))
        else:  # month_days：该月每一天（闰年 2 月 29 天，普通年 28 天）
            y, m = int(part[0:4]), int(part[5:7])
            for day in range(1, calendar.monthrange(y, m)[1] + 1):
                out.add((dep.node, f"{y:04d}-{m:02d}-{day:02d}"))
    return sorted(out)


# --------------------------------------------------------------------------- #
# 账本（SQLite）
# --------------------------------------------------------------------------- #
_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS source_partitions (
  node TEXT NOT NULL,
  part TEXT NOT NULL,
  sha256 TEXT NOT NULL,
  registered_at TEXT NOT NULL,
  PRIMARY KEY (node, part)
);
CREATE TABLE IF NOT EXISTS derived_receipts (
  node TEXT NOT NULL,
  part TEXT NOT NULL,
  output_sha256 TEXT NOT NULL,
  rule_json TEXT NOT NULL,
  inputs_json TEXT NOT NULL,
  confirmed_at TEXT NOT NULL,
  PRIMARY KEY (node, part)
);
"""


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _check_schema(conn, path):
    try:
        row = conn.execute(
            "SELECT value FROM meta WHERE key='schema_version'").fetchone()
    except sqlite3.Error:
        raise UsageError(f"{path} 不是有效的重算规划账本（缺 meta 表）") from None
    if row is None or row[0] != SCHEMA_VERSION:
        raise UsageError(
            f"账本 {path} 的 schema 版本与本工具不符（需要 {SCHEMA_VERSION}）")


def open_ledger(path, create=False, readonly=False):
    """打开账本。readonly 用于规划（规划绝不能写库）；create 用于登记/撤回。"""
    exists = os.path.exists(path)
    if readonly:
        if not exists:
            raise UsageError(f"账本不存在: {path}（请先 register 登记来源分区）")
        uri = "file:" + quote(os.path.abspath(path)) + "?mode=ro"
        conn = sqlite3.connect(uri, uri=True, isolation_level=None)
        _check_schema(conn, path)
        return conn
    if not create and not exists:
        raise UsageError(f"账本不存在: {path}（请先 register 登记来源分区）")
    conn = sqlite3.connect(path, isolation_level=None)
    conn.executescript(_SCHEMA)
    row = conn.execute(
        "SELECT value FROM meta WHERE key='schema_version'").fetchone()
    if row is None:
        conn.execute("INSERT INTO meta(key, value) VALUES('schema_version', ?)",
                     (SCHEMA_VERSION,))
    elif row[0] != SCHEMA_VERSION:
        conn.close()
        raise UsageError(
            f"账本 {path} 的 schema 版本 {row[0]!r} 与本工具不符"
            f"（需要 {SCHEMA_VERSION}）")
    return conn


def get_source_digest(conn, node, part):
    row = conn.execute(
        "SELECT sha256 FROM source_partitions WHERE node=? AND part=?",
        (node, part)).fetchone()
    return row[0] if row else None


def register_source(conn, node, part, digest):
    """登记 source 分区摘要。相同摘要 = 相同内容，不产生新版本。"""
    conn.execute("BEGIN IMMEDIATE")
    try:
        old = get_source_digest(conn, node, part)
        changed = old != digest
        if changed:
            conn.execute(
                "INSERT INTO source_partitions(node, part, sha256,"
                " registered_at) VALUES(?,?,?,?) "
                "ON CONFLICT(node, part) DO UPDATE SET"
                " sha256=excluded.sha256,"
                " registered_at=excluded.registered_at",
                (node, part, digest, _now()))
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return {"changed": changed, "previous_sha256": old}


def retract_source(conn, node, part):
    """撤回 source 分区（表示缺数据；绝不等价于登记一份空表）。"""
    conn.execute("BEGIN IMMEDIATE")
    try:
        cur = conn.execute(
            "DELETE FROM source_partitions WHERE node=? AND part=?",
            (node, part))
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return cur.rowcount > 0


def get_receipt(conn, node, part):
    row = conn.execute(
        "SELECT output_sha256, rule_json, inputs_json, confirmed_at"
        " FROM derived_receipts WHERE node=? AND part=?",
        (node, part)).fetchone()
    if row is None:
        return None
    return {
        "output_sha256": row[0],
        "rule": json.loads(row[1]),
        "inputs": json.loads(row[2]),
        "confirmed_at": row[3],
    }


def store_receipt(conn, node, part, output_sha256, rule, inputs):
    """替换该分区回执（须在调用方的事务内）。"""
    conn.execute(
        "INSERT INTO derived_receipts(node, part, output_sha256, rule_json,"
        " inputs_json, confirmed_at) VALUES(?,?,?,?,?,?) "
        "ON CONFLICT(node, part) DO UPDATE SET"
        " output_sha256=excluded.output_sha256,"
        " rule_json=excluded.rule_json,"
        " inputs_json=excluded.inputs_json,"
        " confirmed_at=excluded.confirmed_at",
        (node, part, output_sha256, canonical_json(rule),
         canonical_json(inputs), _now()))


# --------------------------------------------------------------------------- #
# 状态求值
# --------------------------------------------------------------------------- #
@dataclass
class Status:
    state: str                 # present|missing|current|ready|waiting|blocked
    digest: str | None         # 有效内容摘要（present/current 时非空）
    reason: dict = field(default_factory=dict)


class Evaluator:
    """对闭包内每个 (节点, 分区) 求状态；只读账本，带备忘。"""

    def __init__(self, cfg, conn):
        self.cfg = cfg
        self.conn = conn
        self.memo = {}

    def evaluate(self, node, part):
        key = (node, part)
        st = self.memo.get(key)
        if st is None:
            st = self._evaluate(node, part)
            self.memo[key] = st
        return st

    def _evaluate(self, node, part):
        nc = self.cfg.nodes[node]
        if nc.kind == KIND_SOURCE:
            digest = get_source_digest(self.conn, node, part)
            if digest is None:
                return Status("missing", None, {
                    "code": "missing",
                    "detail": f"source 分区 {node}:{part} 未登记（或已撤回）",
                    "missing_chains": [[node, part]],
                })
            return Status("present", digest, {"code": "present"})

        inputs = expand_inputs(self.cfg, node, part)
        dep_status = {ip: self.evaluate(*ip) for ip in inputs}
        unavailable = {ip: s for ip, s in dep_status.items()
                       if s.digest is None}

        # 缺失链：直接缺 source，或某个直接输入自身已 blocked
        chains = []
        for ip, s in dep_status.items():
            if s.state == "missing":
                chains.append([[node, part], [ip[0], ip[1]]])
            elif s.state == "blocked":
                for ch in s.reason["missing_chains"]:
                    chains.append([[node, part]] + ch)

        if unavailable:
            if chains:
                return Status("blocked", None,
                              {"code": "blocked", "missing_chains": chains})
            pend = [{"node": ip[0], "partition": ip[1], "status": s.state}
                    for ip, s in sorted(dep_status.items())
                    if s.digest is None]
            return Status("waiting", None,
                          {"code": "waiting", "pending_inputs": pend})

        # 直接输入全部有效：回执规则未变且摘要逐项相等才是 current
        receipt = get_receipt(self.conn, node, part)
        rule = canonical_rule(self.cfg, node)
        if (receipt is not None and receipt["rule"] == rule
                and _receipt_inputs_match(receipt, inputs, dep_status)):
            return Status("current", receipt["output_sha256"], {
                "code": "current",
                "output_sha256": receipt["output_sha256"],
                "confirmed_at": receipt["confirmed_at"],
            })
        return Status("ready", None,
                      _ready_reason(receipt, rule, inputs, dep_status))


def _receipt_inputs_match(receipt, inputs, dep_status):
    recorded = {(i["node"], i["partition"]): i["sha256"]
                for i in receipt["inputs"]}
    if len(recorded) != len(inputs) or len(receipt["inputs"]) != len(inputs):
        return False
    for ip in inputs:
        if recorded.get(ip) != dep_status[ip].digest:
            return False
    return True


def _ready_reason(receipt, rule, inputs, dep_status):
    """ready 的具体差异：从未计算 / 规则变了 / 输入摘要变了（可叠加）。"""
    if receipt is None:
        return {"code": "no_receipt", "detail": "从未成功计算过该分区"}
    old = receipt["rule"]
    recipe_changed = old.get("recipe") != rule["recipe"]
    deps_changed = (old.get("dependencies") != rule["dependencies"]
                    or old.get("granularity") != rule["granularity"])
    diffs = []
    recorded = {(i["node"], i["partition"]): i["sha256"]
                for i in receipt["inputs"]}
    for ip in inputs:
        rec = recorded.get(ip)
        cur = dep_status[ip].digest
        if rec != cur:
            diffs.append({"node": ip[0], "partition": ip[1],
                          "recorded_sha256": rec, "current_sha256": cur})
    for ip in sorted(set(recorded) - set(inputs)):
        diffs.append({"node": ip[0], "partition": ip[1],
                      "recorded_sha256": recorded[ip], "current_sha256": None,
                      "note": "回执中的输入不在当前规则展开内"})
    code = "rule_changed" if (recipe_changed or deps_changed) \
        else "inputs_changed"
    return {
        "code": code,
        "recipe_changed": recipe_changed,
        "dependencies_changed": deps_changed,
        "input_diffs": diffs,
    }


# --------------------------------------------------------------------------- #
# 规划
# --------------------------------------------------------------------------- #
def parse_target(cfg, text):
    if ":" not in text:
        raise UsageError(f"目标 {text!r} 需要 节点:分区标识 格式")
    node, part = text.split(":", 1)
    if node not in cfg.nodes:
        raise UsageError(f"目标 {text!r}: 节点 {node!r} 未在配置中定义")
    validate_partition(cfg.nodes[node].granularity, part)
    return node, part


def make_ticket(cfg, evaluator, node, part):
    """为 ready 项生成输入票据：本节点规则 + 稳定顺序的直接输入指纹。"""
    inputs = expand_inputs(cfg, node, part)
    fps = []
    for ip in inputs:
        st = evaluator.evaluate(*ip)
        fps.append({"node": ip[0], "partition": ip[1], "sha256": st.digest})
    return {
        "ticket_version": TICKET_VERSION,
        "node": node,
        "partition": part,
        "rule": canonical_rule(cfg, node),
        "inputs": fps,
    }


def build_plan(cfg, conn, targets):
    """沿声明递归展开目标真正需要的输入，结果只覆盖这个闭包。

    规划只读账本，不把任何任务标成已算。返回 (plan文档, tickets字典)。
    """
    ev = Evaluator(cfg, conn)
    parsed = list(dict.fromkeys(parse_target(cfg, t) for t in targets))
    for node, part in parsed:
        ev.evaluate(node, part)

    items = []
    tickets = {}
    for (node, part) in sorted(ev.memo):
        st = ev.memo[(node, part)]
        nc = cfg.nodes[node]
        item = {"node": node, "partition": part, "kind": nc.kind,
                "status": st.state}
        if nc.kind == KIND_SOURCE:
            if st.state == "present":
                item["sha256"] = st.digest
            else:
                item["reason"] = {"code": "missing",
                                  "detail": st.reason["detail"]}
        else:
            if st.state == "current":
                item["output_sha256"] = st.digest
            else:
                item["reason"] = _format_reason(st)
            if st.state == "ready":
                ticket = make_ticket(cfg, ev, node, part)
                tickets[(node, part)] = ticket
                item["ticket"] = ticket
        items.append(item)

    summary = {}
    for it in items:
        summary[it["status"]] = summary.get(it["status"], 0) + 1
    ready = [f"{n}:{p}" for (n, p) in sorted(ev.memo)
             if ev.memo[(n, p)].state == "ready"]
    plan = {
        "targets": [{"node": n, "partition": p} for n, p in parsed],
        "summary": summary,
        "ready": ready,
        "items": items,
    }
    return plan, tickets


def _format_reason(st):
    r = st.reason
    if r["code"] == "blocked":
        return {
            "code": "blocked",
            "missing": [
                {"chain": [{"node": n, "partition": p} for n, p in chain],
                 "missing_source": {"node": chain[-1][0],
                                    "partition": chain[-1][1]}}
                for chain in r["missing_chains"]
            ],
        }
    return r


def _write_tickets(tickets, directory, plan):
    os.makedirs(directory, exist_ok=True)
    for (node, part), ticket in sorted(tickets.items()):
        path = os.path.join(directory, f"ticket_{node}_{part}.json")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(ticket, indent=2, sort_keys=True,
                                ensure_ascii=False) + "\n")
        for item in plan["items"]:
            if item["node"] == node and item["partition"] == part:
                item["ticket_path"] = path


# --------------------------------------------------------------------------- #
# 完成确认（单事务重验）
# --------------------------------------------------------------------------- #
def confirm(cfg, conn, ticket, output_sha256):
    """凭票据与新输出摘要替换回执；规则/输入任一不符即拒收，原回执保留。"""
    node, part, trule, tinputs = _validate_ticket(cfg, ticket)
    out_digest = normalize_sha256(output_sha256, "输出摘要")

    conn.execute("BEGIN IMMEDIATE")
    try:
        rule = canonical_rule(cfg, node)
        if trule != rule:
            raise ConfirmRejected(
                f"票据规则与当前配置不符（节点 {node} 的 recipe/依赖定义"
                "已变化），拒收旧结果，原回执保留",
                detail={"node": node, "partition": part,
                        "cause": "rule_changed",
                        "ticket_rule": trule, "current_rule": rule})
        inputs = expand_inputs(cfg, node, part)
        tmap = {(i["node"], i["partition"]): i["sha256"] for i in tinputs}
        if len(tmap) != len(tinputs) or set(tmap) != set(inputs):
            raise ConfirmRejected(
                f"票据输入清单与节点 {node} 当前规则展开不一致，拒收旧结果",
                detail={"node": node, "partition": part,
                        "cause": "inputs_mismatch",
                        "ticket_inputs": [list(ip) for ip in sorted(tmap)],
                        "expected_inputs": [list(ip) for ip in inputs]})
        ev = Evaluator(cfg, conn)
        for ip in inputs:
            st = ev.evaluate(*ip)
            want = tmap[ip]
            if st.digest is None:
                raise ConfirmRejected(
                    f"直接输入 {ip[0]}:{ip[1]} 已失效（当前状态 {st.state}），"
                    "拒收旧结果，原回执保留",
                    detail={"node": node, "partition": part,
                            "cause": "input_invalid",
                            "input": {"node": ip[0], "partition": ip[1],
                                      "status": st.state}})
            if st.digest != want:
                raise ConfirmRejected(
                    f"直接输入 {ip[0]}:{ip[1]} 的内容摘要已变化"
                    f"（票据 {want}，当前 {st.digest}），拒收旧结果，"
                    "原回执保留",
                    detail={"node": node, "partition": part,
                            "cause": "input_digest_changed",
                            "input": {"node": ip[0], "partition": ip[1],
                                      "ticket_sha256": want,
                                      "current_sha256": st.digest}})
        store_receipt(conn, node, part, out_digest, rule, tinputs)
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return {"node": node, "partition": part, "output_sha256": out_digest,
            "inputs": tinputs, "recorded": True}


def _validate_ticket(cfg, ticket):
    """票据结构完整性：节点、粒度、依赖字段一个都不能绕过。"""
    if not isinstance(ticket, dict):
        raise UsageError("票据必须是 JSON 对象")
    if ticket.get("ticket_version") != TICKET_VERSION:
        raise UsageError(f"票据版本不符（需要 ticket_version={TICKET_VERSION}）")
    node = ticket.get("node")
    if not isinstance(node, str) or node not in cfg.nodes:
        raise UsageError(f"票据引用了未定义的节点 {node!r}")
    nc = cfg.nodes[node]
    if nc.kind != KIND_DERIVED:
        raise UsageError(
            f"节点 {node} 是 source：源分区用 register 登记，不存在完成确认")
    part = ticket.get("partition")
    if not isinstance(part, str):
        raise UsageError("票据缺 partition 字段")
    validate_partition(nc.granularity, part)
    rule = _normalize_ticket_rule(ticket.get("rule"))
    tinputs = ticket.get("inputs")
    if not isinstance(tinputs, list):
        raise UsageError("票据 inputs 必须是列表")
    norm_inputs = []
    for i, e in enumerate(tinputs):
        if not isinstance(e, dict):
            raise UsageError(f"票据 inputs[{i}] 必须是对象")
        en, ep = e.get("node"), e.get("partition")
        if not isinstance(en, str) or en not in cfg.nodes:
            raise UsageError(f"票据 inputs[{i}] 引用了未定义的节点 {en!r}")
        if not isinstance(ep, str):
            raise UsageError(f"票据 inputs[{i}] 缺 partition 字段")
        validate_partition(cfg.nodes[en].granularity, ep)
        norm_inputs.append({
            "node": en, "partition": ep,
            "sha256": normalize_sha256(e.get("sha256"),
                                       f"票据 inputs[{i}] 的摘要")})
    norm_inputs.sort(key=lambda x: (x["node"], x["partition"]))
    return node, part, rule, norm_inputs


def _normalize_ticket_rule(rule):
    if not isinstance(rule, dict):
        raise UsageError("票据缺 rule 对象")
    recipe = rule.get("recipe")
    gran = rule.get("granularity")
    deps = rule.get("dependencies")
    if not isinstance(recipe, str):
        raise UsageError("票据 rule.recipe 必须是字符串")
    if gran not in (GRAN_DAY, GRAN_MONTH):
        raise UsageError("票据 rule.granularity 必须是 day 或 month")
    if not isinstance(deps, list):
        raise UsageError("票据 rule.dependencies 必须是列表")
    norm_deps = []
    for i, d in enumerate(deps):
        if not isinstance(d, dict):
            raise UsageError(f"票据 rule.dependencies[{i}] 必须是对象")
        dn, dm = d.get("node"), d.get("map")
        if not isinstance(dn, str) or not dn:
            raise UsageError(f"票据 rule.dependencies[{i}] 缺节点名")
        if dm not in MAP_KINDS:
            raise UsageError(
                f"票据 rule.dependencies[{i}] 的 map 非法: {dm!r}")
        lo, hi = d.get("lo"), d.get("hi")
        if dm == MAP_DAY_WINDOW:
            if (isinstance(lo, bool) or isinstance(hi, bool)
                    or not isinstance(lo, int) or not isinstance(hi, int)
                    or not lo <= hi <= 0):
                raise UsageError(
                    f"票据 rule.dependencies[{i}] 的 lo/hi 非法"
                    "（day_window 需要整数 lo<=hi<=0）")
        else:
            lo = hi = None
        ig = d.get("input_granularity")
        if ig not in (GRAN_DAY, GRAN_MONTH):
            raise UsageError(
                f"票据 rule.dependencies[{i}] 的 input_granularity 非法")
        norm_deps.append({"node": dn, "map": dm, "lo": lo, "hi": hi,
                          "input_granularity": ig})
    norm_deps.sort(key=lambda x: (x["node"], x["map"],
                                  x["lo"] if x["lo"] is not None else 0,
                                  x["hi"] if x["hi"] is not None else 0))
    return {"granularity": gran, "recipe": recipe, "dependencies": norm_deps}


# --------------------------------------------------------------------------- #
# 命令分发（由 cli 调用）
# --------------------------------------------------------------------------- #
def dispatch(args):
    cfg = load_config_file(args.config)
    sub = args.replan_cmd

    if sub == "load":
        return _config_summary(cfg)

    if sub == "register":
        digest = normalize_sha256(args.sha256)
        node = _require_node(cfg, args.node, KIND_SOURCE)
        validate_partition(node.granularity, args.partition)
        conn = open_ledger(args.ledger, create=True)
        try:
            r = register_source(conn, node.name, args.partition, digest)
        finally:
            conn.close()
        return {"node": node.name, "partition": args.partition,
                "sha256": digest, **r}

    if sub == "retract":
        node = _require_node(cfg, args.node, KIND_SOURCE)
        validate_partition(node.granularity, args.partition)
        conn = open_ledger(args.ledger, create=True)
        try:
            existed = retract_source(conn, node.name, args.partition)
        finally:
            conn.close()
        return {"node": node.name, "partition": args.partition,
                "retracted": True, "was_registered": existed}

    if sub == "plan":
        conn = open_ledger(args.ledger, readonly=True)
        try:
            plan, tickets = build_plan(cfg, conn, args.target)
        finally:
            conn.close()
        if args.tickets_dir:
            _write_tickets(tickets, args.tickets_dir, plan)
        return plan

    if sub == "confirm":
        ticket = _read_ticket(args.ticket)
        conn = open_ledger(args.ledger)
        try:
            return confirm(cfg, conn, ticket, args.output_sha256)
        finally:
            conn.close()

    raise UsageError(f"未知 replan 子命令 {sub!r}")


def _require_node(cfg, name, kind):
    node = cfg.nodes.get(name)
    if node is None:
        raise UsageError(f"节点 {name!r} 未在配置中定义")
    if node.kind != kind:
        raise UsageError(
            f"节点 {name} 是 {node.kind}，该命令只接受 {kind} 节点")
    return node


def _read_ticket(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except OSError as exc:
        raise UsageError(f"无法读取票据 {path}: {exc}") from None
    except json.JSONDecodeError as exc:
        raise UsageError(f"票据 {path} 不是合法 JSON: {exc}") from None


def _config_summary(cfg):
    nodes = []
    edges = []
    for name in sorted(cfg.nodes):
        nc = cfg.nodes[name]
        entry = {"name": nc.name, "type": nc.kind,
                 "granularity": nc.granularity}
        if nc.kind == KIND_DERIVED:
            entry["recipe"] = nc.recipe
            entry["dependencies"] = []
            for d in nc.dependencies:
                dep = {"node": d.node, "map": d.map}
                if d.map == MAP_DAY_WINDOW:
                    dep["lo"], dep["hi"] = d.lo, d.hi
                entry["dependencies"].append(dep)
                edges.append({"from": d.node, "to": nc.name, "map": d.map})
        nodes.append(entry)
    edges.sort(key=lambda e: (e["from"], e["to"], e["map"]))
    return {"nodes": nodes, "edges": edges}
