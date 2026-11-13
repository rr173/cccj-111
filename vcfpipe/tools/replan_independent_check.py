#!/usr/bin/env python3
"""分区重算规划的独立朴素核对器。

**不导入** `vcfpipe.replan` 的任何代码：只用标准库重新实现依赖展开与
回执比较，直接读取配置 JSON 与账本 SQLite 的行，对产品 `plan` 的输出
逐条复核：

* 闭包一致性：核对器自己从目标出发朴素展开到的 (节点, 分区) 集合，必须
  与 plan 覆盖的集合完全相同；
* 状态一致性：每个分区的 present/missing/current/ready/waiting/blocked
  按定义独立重算（回执规则逐项比对 recipe/依赖定义、直接输入摘要逐项
  与当前有效内容比对），必须与 plan 一致；
* ready 列表与票据：ready 集合一致，且每个 ready 项票据里的规则与直接
  输入指纹（按稳定顺序）与独立重算结果一致。

用法：
    python3 -I tools/replan_independent_check.py \
        --config pipeline.json --ledger ledger.db --plan plan.json

退出码：0 一致；1 不一致（逐条打印差异）；2 用法错误。
"""
import argparse
import calendar
import json
import sqlite3
import sys
from datetime import date, timedelta


def load_config(path):
    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    nodes = {}
    for n in raw["nodes"]:
        nodes[n["name"]] = n
    return nodes


def load_ledger(path):
    conn = sqlite3.connect(path)
    sources = {}
    for node, part, sha in conn.execute(
            "SELECT node, part, sha256 FROM source_partitions"):
        sources[(node, part)] = sha
    receipts = {}
    for node, part, out, rule_json, inputs_json in conn.execute(
            "SELECT node, part, output_sha256, rule_json, inputs_json"
            " FROM derived_receipts"):
        receipts[(node, part)] = {
            "output_sha256": out,
            "rule": json.loads(rule_json),
            "inputs": json.loads(inputs_json),
        }
    conn.close()
    return sources, receipts


def expand(nodes, name, part):
    """朴素依赖展开：逐依赖枚举，集合去重，排序返回。"""
    out = set()
    for dep in nodes[name].get("dependencies", []):
        kind = dep["map"]
        if kind == "identity":
            out.add((dep["node"], part))
        elif kind == "day_window":
            y, m, dd = int(part[0:4]), int(part[5:7]), int(part[8:10])
            base = date(y, m, dd)
            for k in range(dep["lo"], dep["hi"] + 1):
                out.add((dep["node"], (base + timedelta(days=k)).isoformat()))
        elif kind == "month_days":
            y, m = int(part[0:4]), int(part[5:7])
            for dd in range(1, calendar.monthrange(y, m)[1] + 1):
                out.add((dep["node"], f"{y:04d}-{m:02d}-{dd:02d}"))
        else:
            raise ValueError(f"未知映射 {kind!r}")
    return sorted(out)


def config_rule(nodes, name):
    """从配置独立整理出节点规则的规范形式（与回执里的 rule 比对用）。"""
    nd = nodes[name]
    deps = []
    for d in nd.get("dependencies", []):
        deps.append({
            "node": d["node"],
            "map": d["map"],
            "lo": d.get("lo"),
            "hi": d.get("hi"),
            "input_granularity": nodes[d["node"]]["granularity"],
        })
    deps.sort(key=lambda x: (x["node"], x["map"],
                             x["lo"] if x["lo"] is not None else 0,
                             x["hi"] if x["hi"] is not None else 0))
    return {"granularity": nd["granularity"], "recipe": nd.get("recipe"),
            "dependencies": deps}


class NaiveEvaluator:
    """按状态定义独立重算；返回 (state, effective_digest)。"""

    def __init__(self, nodes, sources, receipts):
        self.nodes = nodes
        self.sources = sources
        self.receipts = receipts
        self.memo = {}

    def evaluate(self, node, part):
        key = (node, part)
        if key in self.memo:
            return self.memo[key]
        nd = self.nodes[node]
        if nd["type"] == "source":
            st = (("present", self.sources[key]) if key in self.sources
                  else ("missing", None))
            self.memo[key] = st
            return st
        inputs = expand(self.nodes, node, part)
        ist = {ip: self.evaluate(*ip) for ip in inputs}
        if all(s[1] is not None for s in ist.values()):
            rec = self.receipts.get(key)
            if rec is not None and self._rule_ok(node, rec) \
                    and self._inputs_ok(rec, inputs, ist):
                st = ("current", rec["output_sha256"])
            else:
                st = ("ready", None)
        elif any(s[0] in ("missing", "blocked") for s in ist.values()):
            st = ("blocked", None)
        else:
            st = ("waiting", None)
        self.memo[key] = st
        return st

    def _rule_ok(self, node, rec):
        return rec["rule"] == config_rule(self.nodes, node)

    @staticmethod
    def _inputs_ok(rec, inputs, ist):
        recorded = {(i["node"], i["partition"]): i["sha256"]
                    for i in rec["inputs"]}
        if len(recorded) != len(inputs) or len(rec["inputs"]) != len(inputs):
            return False
        return all(recorded.get(ip) == ist[ip][1] for ip in inputs)


def main(argv=None):
    ap = argparse.ArgumentParser(description="replan 独立朴素核对器")
    ap.add_argument("--config", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--plan", required=True, help="产品 plan 的 JSON 输出")
    args = ap.parse_args(argv)

    nodes = load_config(args.config)
    sources, receipts = load_ledger(args.ledger)
    with open(args.plan, "r", encoding="utf-8") as fh:
        plan = json.load(fh)

    ev = NaiveEvaluator(nodes, sources, receipts)
    for t in plan["targets"]:
        ev.evaluate(t["node"], t["partition"])

    problems = []

    expected = {f"{n}:{p}": st for (n, p), (st, _) in ev.memo.items()}
    actual = {f"{it['node']}:{it['partition']}": it["status"]
              for it in plan["items"]}
    only_mine = sorted(set(expected) - set(actual))
    only_plan = sorted(set(actual) - set(expected))
    if only_mine:
        problems.append(f"闭包不一致：plan 缺少 {only_mine}")
    if only_plan:
        problems.append(f"闭包不一致：plan 多出 {only_plan}")
    for key in sorted(set(expected) & set(actual)):
        if expected[key] != actual[key]:
            problems.append(
                f"{key}: 独立重算为 {expected[key]}，plan 为 {actual[key]}")

    my_ready = sorted(k for k, v in expected.items() if v == "ready")
    if my_ready != plan.get("ready"):
        problems.append(
            f"ready 列表不一致：独立重算 {my_ready}，plan {plan.get('ready')}")

    for it in plan["items"]:
        if it["status"] != "ready":
            continue
        node, part = it["node"], it["partition"]
        ticket = it.get("ticket") or {}
        want_inputs = sorted(
            (ip[0], ip[1], ev.memo[ip][1]) for ip in expand(nodes, node, part))
        got_inputs = sorted(
            (e.get("node"), e.get("partition"), e.get("sha256"))
            for e in ticket.get("inputs", []))
        if got_inputs != want_inputs:
            problems.append(f"{node}:{part} 票据输入指纹与独立重算不一致")
        if ticket.get("rule") != config_rule(nodes, node):
            problems.append(f"{node}:{part} 票据规则与配置不一致")

    if problems:
        print("独立朴素核对发现不一致：")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"独立朴素核对通过：{len(expected)} 个闭包分区的状态、ready 列表"
          "与全部票据均与产品规划一致。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
