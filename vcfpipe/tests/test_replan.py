"""分区重算规划（replan）的回归测试。

覆盖：
* 配置校验：未知引用、环、映射与粒度不符、source/derived 结构、
  lo/hi 规则、未知字段、节点与依赖次序不影响含义；
* 分区标识：严格 YYYY-MM-DD / YYYY-MM、真实公历日期、粒度匹配；
* 依赖展开：identity、day_window 两端包含、month_days 闰年/普通年、
  跨依赖去重、日期运算超界；
* 账本：登记幂等/替换、撤回、重开持久、规划只读；
* 状态机：current/ready/waiting/blocked 及具体差异与缺失链；
* 内容传播：相同摘要恢复 current、不同摘要进入 ready、撤回/补回；
* 确认：接受、改源/改 recipe/上游失效/篡改票据均拒收且旧回执不变；
* CLI 退出码（0/1/2）与票据文件；
* 独立朴素核对器（不调用产品规划函数）对小账本的一致性核对。
"""
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "src")
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "replan_independent_check.py")
sys.path.insert(0, SRC)

from vcfpipe import replan  # noqa: E402
from vcfpipe.replan import (  # noqa: E402
    ConfigError, ConfirmRejected, UsageError,
)


def dig(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def demo_cfg(rolling_recipe="rolling-v1", qc_recipe="qc-v1",
             clean_recipe="clean-v1"):
    return {
        "nodes": [
            {"name": "raw", "type": "source", "granularity": "day"},
            {"name": "clean", "type": "derived", "granularity": "day",
             "recipe": clean_recipe,
             "dependencies": [{"node": "raw", "map": "identity"}]},
            {"name": "rolling", "type": "derived", "granularity": "day",
             "recipe": rolling_recipe,
             "dependencies": [{"node": "clean", "map": "day_window",
                               "lo": -2, "hi": 0}]},
            {"name": "monthly", "type": "derived", "granularity": "month",
             "recipe": "monthly-v1",
             "dependencies": [{"node": "clean", "map": "month_days"}]},
            {"name": "qc", "type": "derived", "granularity": "day",
             "recipe": qc_recipe,
             "dependencies": [{"node": "raw", "map": "identity"}]},
        ]
    }


class LedgerCase(unittest.TestCase):
    """每个用例一个临时目录 + 内存级配置。"""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="vcfpipe_replan_test_")
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.db = os.path.join(self.dir, "ledger.db")
        self.cfg = replan.parse_config(demo_cfg())
        self.conn = replan.open_ledger(self.db, create=True)
        self.addCleanup(self.conn.close)

    def register(self, day, content):
        return replan.register_source(self.conn, "raw", day, dig(content))

    def plan(self, *targets):
        return replan.build_plan(self.cfg, self.conn, list(targets))

    def status(self, plan, node, part):
        doc = plan[0] if isinstance(plan, tuple) else plan
        for it in doc["items"]:
            if it["node"] == node and it["partition"] == part:
                return it["status"]
        return None

    def item(self, plan, node, part):
        doc = plan[0] if isinstance(plan, tuple) else plan
        for it in doc["items"]:
            if it["node"] == node and it["partition"] == part:
                return it
        return None

    def confirm_ready(self, plan, out_of):
        doc, tickets = plan
        for key in doc["ready"]:
            node, part = key.split(":", 1)
            replan.confirm(self.cfg, self.conn, tickets[(node, part)],
                           out_of(node, part))


# --------------------------------------------------------------------------- #
# 配置校验
# --------------------------------------------------------------------------- #
class TestConfig(unittest.TestCase):
    def test_valid_config_and_order_invariance(self):
        a = replan.parse_config(demo_cfg())
        shuffled = demo_cfg()
        shuffled["nodes"] = list(reversed(shuffled["nodes"]))
        b = replan.parse_config(shuffled)
        for name in a.nodes:
            self.assertEqual(replan.canonical_rule(a, name),
                             replan.canonical_rule(b, name))

    def test_duplicate_node_rejected(self):
        cfg = demo_cfg()
        cfg["nodes"].append(cfg["nodes"][0].copy())
        with self.assertRaisesRegex(ConfigError, "重复定义"):
            replan.parse_config(cfg)

    def test_unknown_dependency_rejected(self):
        cfg = {"nodes": [
            {"name": "a", "type": "derived", "granularity": "day",
             "recipe": "x",
             "dependencies": [{"node": "ghost", "map": "identity"}]}]}
        with self.assertRaisesRegex(ConfigError, "未定义的节点 'ghost'"):
            replan.parse_config(cfg)

    def test_cycle_rejected_with_chain(self):
        cfg = {"nodes": [
            {"name": "a", "type": "derived", "granularity": "day",
             "recipe": "a",
             "dependencies": [{"node": "b", "map": "identity"}]},
            {"name": "b", "type": "derived", "granularity": "day",
             "recipe": "b",
             "dependencies": [{"node": "a", "map": "identity"}]}]}
        with self.assertRaisesRegex(ConfigError, "环: a -> b -> a"):
            replan.parse_config(cfg)

    def test_self_loop_rejected(self):
        cfg = {"nodes": [
            {"name": "a", "type": "derived", "granularity": "day",
             "recipe": "a",
             "dependencies": [{"node": "a", "map": "identity"}]}]}
        with self.assertRaisesRegex(ConfigError, "环: a -> a"):
            replan.parse_config(cfg)

    def test_identity_granularity_mismatch(self):
        cfg = {"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "m", "type": "derived", "granularity": "month",
             "recipe": "m",
             "dependencies": [{"node": "s", "map": "identity"}]}]}
        with self.assertRaisesRegex(ConfigError, "identity 要求两端粒度相同"):
            replan.parse_config(cfg)

    def test_day_window_requires_day_day(self):
        cfg = {"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "m", "type": "derived", "granularity": "month",
             "recipe": "m",
             "dependencies": [{"node": "s", "map": "day_window",
                               "lo": -1, "hi": 0}]}]}
        with self.assertRaisesRegex(ConfigError, "day_window 要求两端都是 day"):
            replan.parse_config(cfg)

    def test_month_days_requires_month_out_day_in(self):
        cfg = {"nodes": [
            {"name": "s", "type": "source", "granularity": "month"},
            {"name": "d", "type": "derived", "granularity": "day",
             "recipe": "d",
             "dependencies": [{"node": "s", "map": "month_days"}]}]}
        with self.assertRaisesRegex(ConfigError,
                                    "month_days 要求输出 month、输入 day"):
            replan.parse_config(cfg)

    def test_source_must_not_have_recipe_or_deps(self):
        for extra in ({"recipe": "x"}, {"dependencies": []}):
            cfg = {"nodes": [
                dict({"name": "s", "type": "source", "granularity": "day"},
                     **extra),
                {"name": "d", "type": "derived", "granularity": "day",
                 "recipe": "d",
                 "dependencies": [{"node": "s", "map": "identity"}]}]}
            with self.assertRaisesRegex(ConfigError, "source 节点不带"):
                replan.parse_config(cfg)

    def test_derived_needs_recipe_and_deps(self):
        with self.assertRaisesRegex(ConfigError, "非空 recipe"):
            replan.parse_config({"nodes": [
                {"name": "d", "type": "derived", "granularity": "day",
                 "dependencies": [{"node": "s", "map": "identity"}]}]})
        with self.assertRaisesRegex(ConfigError, "至少依赖一项"):
            replan.parse_config({"nodes": [
                {"name": "d", "type": "derived", "granularity": "day",
                 "recipe": "x", "dependencies": []}]})

    def test_day_window_lo_hi_rules(self):
        base = {"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "d", "type": "derived", "granularity": "day",
             "recipe": "d", "dependencies": []}]}
        for lo, hi in [(0, 1), (1, 0), (-2, -3), (True, 0), ("-1", 0)]:
            cfg = json.loads(json.dumps(base))
            cfg["nodes"][1]["dependencies"] = [
                {"node": "s", "map": "day_window", "lo": lo, "hi": hi}]
            with self.assertRaises(ConfigError, msg=f"lo={lo} hi={hi}"):
                replan.parse_config(cfg)
        # identity 不接受 lo/hi
        cfg = json.loads(json.dumps(base))
        cfg["nodes"][1]["dependencies"] = [
            {"node": "s", "map": "identity", "lo": 0, "hi": 0}]
        with self.assertRaisesRegex(ConfigError, "不接受 lo/hi"):
            replan.parse_config(cfg)

    def test_unknown_fields_rejected(self):
        cfg = demo_cfg()
        cfg["nodes"][0]["comment"] = "x"
        with self.assertRaisesRegex(ConfigError, "未知字段"):
            replan.parse_config(cfg)
        cfg = demo_cfg()
        cfg["nodes"][1]["dependencies"][0]["note"] = 1
        with self.assertRaisesRegex(ConfigError, "未知字段"):
            replan.parse_config(cfg)

    def test_bad_node_name_rejected(self):
        cfg = {"nodes": [
            {"name": "a/b", "type": "source", "granularity": "day"}]}
        with self.assertRaisesRegex(ConfigError, "节点名"):
            replan.parse_config(cfg)


# --------------------------------------------------------------------------- #
# 分区标识与展开
# --------------------------------------------------------------------------- #
class TestPartitionsAndExpansion(unittest.TestCase):
    def test_partition_validation(self):
        replan.validate_partition("day", "2024-02-29")   # 闰日合法
        replan.validate_partition("month", "2024-02")
        for gran, part in [("day", "2023-02-29"), ("day", "2024-2-8"),
                           ("day", "2024-02-30"), ("day", "2024-02"),
                           ("month", "2024-02-28"), ("month", "2024-13"),
                           ("month", "2024-2")]:
            with self.assertRaises(UsageError, msg=f"{gran} {part}"):
                replan.validate_partition(gran, part)

    def test_identity_and_window_expansion(self):
        cfg = replan.parse_config(demo_cfg())
        self.assertEqual(replan.expand_inputs(cfg, "clean", "2024-02-28"),
                         [("raw", "2024-02-28")])
        self.assertEqual(
            replan.expand_inputs(cfg, "rolling", "2024-03-01"),
            [("clean", "2024-02-28"), ("clean", "2024-02-29"),
             ("clean", "2024-03-01")])

    def test_month_days_leap_and_common(self):
        cfg = replan.parse_config(demo_cfg())
        leap = replan.expand_inputs(cfg, "monthly", "2024-02")
        common = replan.expand_inputs(cfg, "monthly", "2023-02")
        self.assertEqual(len(leap), 29)
        self.assertIn(("clean", "2024-02-29"), leap)
        self.assertEqual(len(common), 28)
        self.assertNotIn(("clean", "2023-02-29"), common)

    def test_expansion_dedup_across_dependencies(self):
        cfg = replan.parse_config({"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "d", "type": "derived", "granularity": "day",
             "recipe": "d",
             "dependencies": [
                 {"node": "s", "map": "day_window", "lo": -2, "hi": 0},
                 {"node": "s", "map": "day_window", "lo": -1, "hi": 0},
                 {"node": "s", "map": "identity"}]}]})
        self.assertEqual(
            replan.expand_inputs(cfg, "d", "2024-03-10"),
            [("s", "2024-03-08"), ("s", "2024-03-09"), ("s", "2024-03-10")])

    def test_date_overflow_explained(self):
        cfg = replan.parse_config(demo_cfg())
        with self.assertRaisesRegex(UsageError, "超界"):
            replan.expand_inputs(cfg, "rolling", "0001-01-01")


# --------------------------------------------------------------------------- #
# 账本行为
# --------------------------------------------------------------------------- #
class TestLedger(LedgerCase):
    def test_register_idempotent_and_replace(self):
        r1 = self.register("2024-02-28", "v1")
        self.assertTrue(r1["changed"])
        r2 = self.register("2024-02-28", "v1")
        self.assertFalse(r2["changed"])          # 相同摘要不制造新版本
        r3 = self.register("2024-02-28", "v2")
        self.assertTrue(r3["changed"])
        self.assertEqual(r3["previous_sha256"], dig("v1"))

    def test_retract(self):
        self.register("2024-02-28", "v1")
        self.assertTrue(replan.retract_source(self.conn, "raw", "2024-02-28"))
        self.assertFalse(replan.retract_source(self.conn, "raw", "2024-02-28"))
        self.assertIsNone(
            replan.get_source_digest(self.conn, "raw", "2024-02-28"))

    def test_ledger_survives_reopen(self):
        self.register("2024-02-28", "v1")
        self.conn.close()
        conn = replan.open_ledger(self.db, readonly=True)
        self.assertEqual(replan.get_source_digest(conn, "raw", "2024-02-28"),
                         dig("v1"))
        conn.close()
        self.conn = replan.open_ledger(self.db)

    def test_plan_is_read_only(self):
        self.register("2024-02-28", "v1")
        before = self._dump()
        self.plan("clean:2024-02-28", "monthly:2024-02")
        self.assertEqual(before, self._dump())   # 规划不把任何任务标成已算

    def _dump(self):
        s = list(self.conn.execute(
            "SELECT * FROM source_partitions ORDER BY 1, 2"))
        d = list(self.conn.execute(
            "SELECT * FROM derived_receipts ORDER BY 1, 2"))
        return s, d

    def test_missing_ledger_rejected(self):
        with self.assertRaisesRegex(UsageError, "账本不存在"):
            replan.open_ledger(os.path.join(self.dir, "nope.db"),
                               readonly=True)


# --------------------------------------------------------------------------- #
# 状态机与内容传播
# --------------------------------------------------------------------------- #
class TestStatusFlow(LedgerCase):
    FEB = [f"2024-02-{i:02d}" for i in range(1, 30)]
    DAYS = FEB + ["2024-03-01", "2024-03-02"]
    TARGETS = ["rolling:2024-02-28", "rolling:2024-03-02", "monthly:2024-02"]

    def setUp(self):
        super().setUp()
        for day in self.DAYS:
            self.register(day, f"raw-v1:{day}")

    def compute_all(self):
        for _ in range(6):
            plan = self.plan(*self.TARGETS)
            if not plan[0]["ready"]:
                return plan
            self.confirm_ready(
                plan, lambda n, p: dig(f"{n}-out:{p}"))
        raise AssertionError("规划轮次超限")

    def test_initial_statuses(self):
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(self.status(plan, "raw", "2024-02-28"), "present")
        self.assertEqual(self.status(plan, "clean", "2024-02-28"), "ready")
        self.assertEqual(self.item(plan, "clean", "2024-02-28")
                         ["reason"]["code"], "no_receipt")
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"),
                         "waiting")
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "waiting")
        pend = self.item(plan, "rolling", "2024-02-28")["reason"]
        self.assertEqual(pend["code"], "waiting")
        self.assertTrue(any(p["node"] == "clean"
                            for p in pend["pending_inputs"]))

    def test_missing_source_blocks_with_chain(self):
        replan.retract_source(self.conn, "raw", "2024-02-15")
        plan, _ = self.plan("monthly:2024-02")
        self.assertEqual(self.status(plan, "raw", "2024-02-15"), "missing")
        self.assertEqual(self.status(plan, "clean", "2024-02-15"), "blocked")
        it = self.item(plan, "monthly", "2024-02")
        self.assertEqual(it["status"], "blocked")
        chains = [[f"{x['node']}:{x['partition']}" for x in m["chain"]]
                  for m in it["reason"]["missing"]]
        self.assertIn(["monthly:2024-02", "clean:2024-02-15",
                       "raw:2024-02-15"], chains)

    def test_all_current_after_compute(self):
        plan, _ = self.compute_all()
        self.assertEqual(plan["ready"], [])
        self.assertLessEqual(set(plan["summary"]), {"present", "current"})

    def test_backfill_affects_window_and_month_only(self):
        self.compute_all()
        self.register("2024-02-28", "raw-v2:2024-02-28")
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(self.status(plan, "clean", "2024-02-28"), "ready")
        reason = self.item(plan, "clean", "2024-02-28")["reason"]
        self.assertEqual(reason["code"], "inputs_changed")
        diff = reason["input_diffs"][0]
        self.assertEqual((diff["node"], diff["partition"]),
                         ("raw", "2024-02-28"))
        self.assertEqual(diff["recorded_sha256"], dig("raw-v1:2024-02-28"))
        self.assertEqual(diff["current_sha256"], dig("raw-v2:2024-02-28"))
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"),
                         "waiting")
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "waiting")
        # 窗口 [02-29..03-02] 不含 02-28：不波及
        self.assertEqual(self.status(plan, "rolling", "2024-03-02"),
                         "current")

    def test_same_output_digest_restores_current(self):
        self.compute_all()
        self.register("2024-02-28", "raw-v2:2024-02-28")
        plan, tickets = self.plan(*self.TARGETS)
        # clean 重算，交回与以前相同的输出摘要
        replan.confirm(self.cfg, self.conn,
                       tickets[("clean", "2024-02-28")],
                       dig("clean-out:2024-02-28"))
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(plan["ready"], [])
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"),
                         "current")
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "current")

    def test_different_output_digest_makes_downstream_ready(self):
        self.compute_all()
        self.register("2024-02-28", "raw-v2:2024-02-28")
        plan, tickets = self.plan(*self.TARGETS)
        replan.confirm(self.cfg, self.conn,
                       tickets[("clean", "2024-02-28")],
                       dig("clean-out-v2:2024-02-28"))
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"), "ready")
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "ready")
        self.assertEqual(self.status(plan, "rolling", "2024-03-02"),
                         "current")

    def test_retract_then_restore_same_content(self):
        self.compute_all()
        replan.retract_source(self.conn, "raw", "2024-02-29")
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "blocked")
        self.assertEqual(self.status(plan, "rolling", "2024-03-02"),
                         "blocked")
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"),
                         "current")  # 窗口 [02-26..02-28] 不含 02-29
        self.register("2024-02-29", "raw-v1:2024-02-29")  # 补回同样内容
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(plan["ready"], [])
        self.assertLessEqual(set(plan["summary"]), {"present", "current"})

    def test_recipe_change_only_affects_its_node(self):
        self.compute_all()
        # 先算好 qc，再只改它的 recipe：主链闭包不受影响
        plan, tickets = self.plan("qc:2024-02-28")
        replan.confirm(self.cfg, self.conn, tickets[("qc", "2024-02-28")],
                       dig("qc-out:2024-02-28"))
        self.cfg = replan.parse_config(demo_cfg(qc_recipe="qc-v2"))
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(plan["ready"], [])
        plan_qc, _ = self.plan("qc:2024-02-28")
        it = self.item(plan_qc, "qc", "2024-02-28")
        self.assertEqual(it["status"], "ready")
        self.assertEqual(it["reason"]["code"], "rule_changed")
        self.assertTrue(it["reason"]["recipe_changed"])
        # 改 clean 的 recipe：clean 与其下游都按规则失效
        self.cfg = replan.parse_config(demo_cfg(clean_recipe="clean-v2"))
        plan, _ = self.plan(*self.TARGETS)
        self.assertEqual(self.status(plan, "clean", "2024-02-28"), "ready")
        self.assertEqual(self.status(plan, "rolling", "2024-02-28"),
                         "waiting")
        self.assertEqual(self.status(plan, "monthly", "2024-02"), "waiting")

    def test_multi_target_merges_shared_partitions(self):
        plan, _ = self.plan("rolling:2024-02-29", "monthly:2024-02")
        keys = [(it["node"], it["partition"]) for it in plan["items"]]
        self.assertEqual(len(keys), len(set(keys)))   # 合并展示，无重复
        self.assertEqual(keys.count(("clean", "2024-02-28")), 1)

    def test_ready_sorted_and_ticket_stable(self):
        plan, tickets = self.plan("monthly:2024-02")
        self.assertEqual(plan["ready"], sorted(plan["ready"]))
        tk = tickets[("clean", "2024-02-28")]
        keys = [(i["node"], i["partition"]) for i in tk["inputs"]]
        self.assertEqual(keys, sorted(keys))
        self.assertEqual(tk["rule"]["recipe"], "clean-v1")
        self.assertEqual(tk["inputs"][0]["sha256"], dig("raw-v1:2024-02-28"))


# --------------------------------------------------------------------------- #
# 完成确认
# --------------------------------------------------------------------------- #
class TestConfirm(LedgerCase):
    def setUp(self):
        super().setUp()
        self.register("2024-02-26", "raw-v1:2024-02-26")
        self.register("2024-02-27", "raw-v1:2024-02-27")
        self.register("2024-02-28", "raw-v1:2024-02-28")
        # 先算好 clean
        for day in ("2024-02-26", "2024-02-27", "2024-02-28"):
            plan, tickets = self.plan(f"clean:{day}")
            replan.confirm(self.cfg, self.conn, tickets[("clean", day)],
                           dig(f"clean-out:{day}"))

    def plan(self, *targets):
        return replan.build_plan(self.cfg, self.conn, list(targets))

    def ticket(self, node, part):
        plan, tickets = self.plan(f"{node}:{part}")
        self.assertIn((node, part), tickets)
        return tickets[(node, part)]

    def test_confirm_accept_writes_receipt(self):
        tk = self.ticket("rolling", "2024-02-28")
        replan.confirm(self.cfg, self.conn, tk, dig("rolling-out:2024-02-28"))
        rec = replan.get_receipt(self.conn, "rolling", "2024-02-28")
        self.assertEqual(rec["output_sha256"], dig("rolling-out:2024-02-28"))
        self.assertEqual(rec["rule"]["recipe"], "rolling-v1")
        self.assertEqual(len(rec["inputs"]), 3)

    def test_source_change_after_ticket_rejected(self):
        tk = self.ticket("rolling", "2024-02-28")
        self.register("2024-02-28", "raw-v2:2024-02-28")  # 拿票后改源
        with self.assertRaisesRegex(ConfirmRejected, "已失效|已变化"):
            replan.confirm(self.cfg, self.conn, tk,
                           dig("rolling-out:2024-02-28"))
        self.assertIsNone(
            replan.get_receipt(self.conn, "rolling", "2024-02-28"))

    def test_recipe_change_after_ticket_rejected_and_preserved(self):
        tk = self.ticket("rolling", "2024-02-28")
        replan.confirm(self.cfg, self.conn, tk, dig("rolling-out:2024-02-28"))
        tk2 = json.loads(json.dumps(tk))  # 旧票据（规则 rolling-v1）
        self.cfg = replan.parse_config(demo_cfg(rolling_recipe="rolling-v2"))
        with self.assertRaisesRegex(ConfirmRejected, "规则"):
            replan.confirm(self.cfg, self.conn, tk2, dig("rolling-out2"))
        rec = replan.get_receipt(self.conn, "rolling", "2024-02-28")
        self.assertEqual(rec["output_sha256"], dig("rolling-out:2024-02-28"))
        self.assertEqual(rec["rule"]["recipe"], "rolling-v1")  # 旧回执不变

    def test_upstream_invalidated_rejected(self):
        tk = self.ticket("rolling", "2024-02-28")
        # 上游 clean 的输入被撤回 → clean 失效 → rolling 确认必须拒收
        replan.retract_source(self.conn, "raw", "2024-02-27")
        with self.assertRaisesRegex(ConfirmRejected, "已失效"):
            replan.confirm(self.cfg, self.conn, tk, dig("rolling-out"))

    def test_tampered_ticket_rejected(self):
        tk = self.ticket("rolling", "2024-02-28")
        # 篡改输入摘要
        bad = json.loads(json.dumps(tk))
        bad["inputs"][0]["sha256"] = dig("forged")
        with self.assertRaises(ConfirmRejected):
            replan.confirm(self.cfg, self.conn, bad, dig("out"))
        # 额外塞进一个输入
        bad2 = json.loads(json.dumps(tk))
        bad2["inputs"].append({"node": "raw", "partition": "2024-02-25",
                               "sha256": dig("extra")})
        with self.assertRaises(ConfirmRejected):
            replan.confirm(self.cfg, self.conn, bad2, dig("out"))
        # 改掉节点名（结构完整性）
        bad3 = json.loads(json.dumps(tk))
        bad3["node"] = "ghost"
        with self.assertRaises(UsageError):
            replan.confirm(self.cfg, self.conn, bad3, dig("out"))
        # source 节点不存在完成确认
        bad4 = json.loads(json.dumps(tk))
        bad4["node"] = "raw"
        with self.assertRaisesRegex(UsageError, "source"):
            replan.confirm(self.cfg, self.conn, bad4, dig("out"))

    def test_confirm_rejects_bad_output_digest(self):
        tk = self.ticket("rolling", "2024-02-28")
        with self.assertRaises(UsageError):
            replan.confirm(self.cfg, self.conn, tk, "not-a-sha")

    def test_confirm_does_not_freshen_stale_inputs(self):
        """确认时看到的新摘要不能写进一份用旧数据算出的回执。"""
        tk = self.ticket("rolling", "2024-02-28")
        self.register("2024-02-28", "raw-v2:2024-02-28")
        # clean[2024-02-28] 已过期；确认 rolling 必须拒收而不是采纳新摘要
        with self.assertRaises(ConfirmRejected):
            replan.confirm(self.cfg, self.conn, tk, dig("rolling-out"))
        rec = replan.get_receipt(self.conn, "clean", "2024-02-28")
        self.assertEqual(rec["inputs"][0]["sha256"], dig("raw-v1:2024-02-28"))


# --------------------------------------------------------------------------- #
# CLI 子进程回归
# --------------------------------------------------------------------------- #
class TestCli(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="vcfpipe_replan_cli_")
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.config = os.path.join(self.dir, "cfg.json")
        with open(self.config, "w", encoding="utf-8") as fh:
            json.dump(demo_cfg(), fh)
        self.db = os.path.join(self.dir, "ledger.db")

    def run_cli(self, args, expect=0):
        proc = subprocess.run(
            [sys.executable, "-I", RUNNER, "replan"] + args,
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, expect,
                         f"{args}\nstdout={proc.stdout}\nstderr={proc.stderr}")
        return proc

    def test_full_flow_and_exit_codes(self):
        self.run_cli(["load", "--config", self.config])
        self.run_cli(["register", "--config", self.config,
                      "--ledger", self.db, "--node", "raw",
                      "--partition", "2024-02-28", "--sha256", dig("v1")])
        tickets = os.path.join(self.dir, "tickets")
        proc = self.run_cli(["plan", "--config", self.config,
                             "--ledger", self.db,
                             "--target", "clean:2024-02-28",
                             "--tickets-dir", tickets])
        plan = json.loads(proc.stdout)
        self.assertEqual(plan["ready"], ["clean:2024-02-28"])
        ticket_path = plan["items"][0]["ticket_path"]
        self.assertTrue(os.path.exists(ticket_path))
        self.run_cli(["confirm", "--config", self.config, "--ledger", self.db,
                      "--ticket", ticket_path,
                      "--output-sha256", dig("clean-out")])
        proc = self.run_cli(["plan", "--config", self.config,
                             "--ledger", self.db,
                             "--target", "clean:2024-02-28"])
        self.assertEqual(json.loads(proc.stdout)["summary"],
                         {"current": 1, "present": 1})
        # 拿票后改源 → 确认退出码 1
        self.run_cli(["register", "--config", self.config,
                      "--ledger", self.db, "--node", "raw",
                      "--partition", "2024-02-28", "--sha256", dig("v2")])
        self.run_cli(["confirm", "--config", self.config, "--ledger", self.db,
                      "--ticket", ticket_path,
                      "--output-sha256", dig("clean-out")], expect=1)
        # 非法日期 / 未知节点 / 缺账本 → 退出码 2
        self.run_cli(["register", "--config", self.config,
                      "--ledger", self.db, "--node", "raw",
                      "--partition", "2023-02-29", "--sha256", dig("x")],
                     expect=2)
        self.run_cli(["plan", "--config", self.config,
                      "--ledger", self.db, "--target", "ghost:2024-02-28"],
                     expect=2)
        self.run_cli(["plan", "--config", self.config,
                      "--ledger", os.path.join(self.dir, "none.db"),
                      "--target", "clean:2024-02-28"], expect=2)
        # 撤回 → 0；blocked 出现在规划中
        self.run_cli(["retract", "--config", self.config,
                      "--ledger", self.db, "--node", "raw",
                      "--partition", "2024-02-28"])
        proc = self.run_cli(["plan", "--config", self.config,
                             "--ledger", self.db,
                             "--target", "clean:2024-02-28"])
        self.assertEqual(json.loads(proc.stdout)["summary"],
                         {"blocked": 1, "missing": 1})

    def test_bad_config_exit_2(self):
        bad = os.path.join(self.dir, "bad.json")
        with open(bad, "w", encoding="utf-8") as fh:
            json.dump({"nodes": [
                {"name": "a", "type": "derived", "granularity": "day",
                 "recipe": "a",
                 "dependencies": [{"node": "a", "map": "identity"}]}]}, fh)
        proc = self.run_cli(["load", "--config", bad], expect=2)
        self.assertIn("环", proc.stderr)


# --------------------------------------------------------------------------- #
# 独立朴素核对器
# --------------------------------------------------------------------------- #
class TestIndependentCheck(LedgerCase):
    def _scenario(self):
        """构造一个小账本：部分已算、部分待算、一条缺失链。"""
        for day in ("2024-02-25", "2024-02-26", "2024-02-27", "2024-02-28"):
            self.register(day, f"raw-v1:{day}")
        for day in ("2024-02-25", "2024-02-26", "2024-02-27"):
            plan, tickets = self.plan(f"clean:{day}")
            replan.confirm(self.cfg, self.conn, tickets[("clean", day)],
                           dig(f"clean-out:{day}"))
        plan, tickets = self.plan("rolling:2024-02-27")
        replan.confirm(self.cfg, self.conn,
                       tickets[("rolling", "2024-02-27")],
                       dig("rolling-out:2024-02-27"))
        # 留下：clean[2024-02-28] ready、rolling[2024-02-28] waiting、
        # monthly[2024-02] blocked（2 月其余天未登记）
        return self.plan("rolling:2024-02-28", "monthly:2024-02",
                         "rolling:2024-02-27")

    def test_checker_agrees_with_product(self):
        plan, _ = self._scenario()
        plan_path = os.path.join(self.dir, "plan.json")
        with open(plan_path, "w", encoding="utf-8") as fh:
            json.dump(plan, fh)
        cfg_path = os.path.join(self.dir, "cfg.json")
        with open(cfg_path, "w", encoding="utf-8") as fh:
            json.dump(demo_cfg(), fh)
        proc = subprocess.run(
            [sys.executable, "-I", CHECKER, "--config", cfg_path,
             "--ledger", self.db, "--plan", plan_path],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("核对通过", proc.stdout)

    def test_checker_detects_tampered_plan(self):
        plan, _ = self._scenario()
        for it in plan["items"]:
            if it["status"] == "ready":
                it["status"] = "current"     # 篡改一处状态
                break
        plan_path = os.path.join(self.dir, "plan.json")
        with open(plan_path, "w", encoding="utf-8") as fh:
            json.dump(plan, fh)
        cfg_path = os.path.join(self.dir, "cfg.json")
        with open(cfg_path, "w", encoding="utf-8") as fh:
            json.dump(demo_cfg(), fh)
        proc = subprocess.run(
            [sys.executable, "-I", CHECKER, "--config", cfg_path,
             "--ledger", self.db, "--plan", plan_path],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("不一致", proc.stdout)


if __name__ == "__main__":
    unittest.main()
