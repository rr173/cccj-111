#!/usr/bin/env python3
"""分区重算规划（replan）合成账本操作稿。

在临时目录里驱动 `run_pipeline.py replan` 的完整演示，逐步展示：

  1. 全部算完后无任务待做；
  2. 补交 2 月 28 日只影响窗口内恰当的 rolling 日期与月表，不波及无关日；
  3. 闰年 2 月需要 29 天，普通年只需 28 天；
  4. 撤回月末一天后阻塞链点名那一天，补回后可以继续；
  5. clean 重算交回相同摘要时下游直接恢复 current，不同摘要才进入 ready；
  6. 票据发出后改源或改 recipe，旧确认被拒绝且旧回执不变；
  7. 只调整无关支线（qc 的 recipe）不引起重算；
  8. 循环定义与不匹配的粒度被实际拒绝。

依赖图（两条支线共用 clean）：

                 identity                day_window [-2,0]
  raw[日] ─────────────────> clean[日] ─────────────────────> rolling[日]
  (source)                    │                             输出 d 读 clean[d-2..d]
                              │ month_days
                              v
                            monthly[月]    输出 YYYY-MM 读当月全部 clean 日分区

  另有无关支线：qc[日] --identity--> raw[日]（nothing depends on qc）。

运行：python3 -I tools/run_replan_demo.py
"""
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "replan_independent_check.py")

FAILURES = []


def dig(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def step(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def check(cond, message):
    if cond:
        print(f"  ✓ {message}")
    else:
        FAILURES.append(message)
        print(f"  ✗ 不符合预期: {message}")


class Demo:
    def __init__(self, workdir):
        self.dir = workdir
        self.ledger = os.path.join(workdir, "ledger.db")
        self.tickets = os.path.join(workdir, "tickets")
        self.config = os.path.join(workdir, "pipeline.json")
        self.write_config(rolling_recipe="rolling-v1", qc_recipe="qc-v1")

    def write_config(self, rolling_recipe, qc_recipe):
        cfg = {
            "nodes": [
                {"name": "raw", "type": "source", "granularity": "day"},
                {"name": "clean", "type": "derived", "granularity": "day",
                 "recipe": "clean-v1",
                 "dependencies": [{"node": "raw", "map": "identity"}]},
                {"name": "rolling", "type": "derived", "granularity": "day",
                 "recipe": rolling_recipe,
                 "dependencies": [{"node": "clean", "map": "day_window",
                                   "lo": -2, "hi": 0}]},
                {"name": "monthly", "type": "derived",
                 "granularity": "month", "recipe": "monthly-v1",
                 "dependencies": [{"node": "clean", "map": "month_days"}]},
                {"name": "qc", "type": "derived", "granularity": "day",
                 "recipe": qc_recipe,
                 "dependencies": [{"node": "raw", "map": "identity"}]},
            ]
        }
        with open(self.config, "w", encoding="utf-8") as fh:
            json.dump(cfg, fh, ensure_ascii=False, indent=2)

    # ---------------- CLI 包装 ---------------- #
    def run(self, args, expect=0, keep_json=None):
        cmd = [sys.executable, "-I", RUNNER, "replan"] + args
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != expect:
            raise AssertionError(
                f"命令 {' '.join(args)} 退出码 {proc.returncode} != {expect}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
        out = proc.stdout.strip()
        data = json.loads(out) if out else {}
        if keep_json:
            with open(keep_json, "w", encoding="utf-8") as fh:
                fh.write(out)
        return data

    def register(self, day, content):
        r = self.run(["register", "--config", self.config,
                      "--ledger", self.ledger, "--node", "raw",
                      "--partition", day, "--sha256", dig(content)])
        return r

    def retract(self, day):
        return self.run(["retract", "--config", self.config,
                         "--ledger", self.ledger, "--node", "raw",
                         "--partition", day])

    def plan(self, targets, keep_json=None):
        args = ["plan", "--config", self.config, "--ledger", self.ledger]
        for t in targets:
            args += ["--target", t]
        args += ["--tickets-dir", self.tickets]
        return self.run(args, keep_json=keep_json)

    def confirm(self, node, part, out_digest, ticket_path=None, expect=0):
        if ticket_path is None:
            ticket_path = os.path.join(self.tickets,
                                       f"ticket_{node}_{part}.json")
        return self.run(["confirm", "--config", self.config,
                         "--ledger", self.ledger, "--ticket", ticket_path,
                         "--output-sha256", out_digest], expect=expect)

    # ---------------- 展示辅助 ---------------- #
    def receipt_digest(self, node, part):
        conn = sqlite3.connect(self.ledger)
        row = conn.execute(
            "SELECT output_sha256 FROM derived_receipts"
            " WHERE node=? AND part=?", (node, part)).fetchone()
        conn.close()
        return row[0] if row else None

    @staticmethod
    def brief(plan):
        """按状态分组打印闭包内各项（摘要形式）。"""
        by = {}
        for it in plan["items"]:
            by.setdefault(it["status"], []).append(
                f"{it['node']}[{it['partition']}]")
        for state in ("present", "missing", "current", "ready", "waiting",
                      "blocked"):
            if state in by:
                print(f"    {state:8s} ({len(by[state])}): "
                      + ", ".join(by[state]))

    @staticmethod
    def item(plan, node, part):
        for it in plan["items"]:
            if it["node"] == node and it["partition"] == part:
                return it
        return None

    def compute_ready(self, plan, out_digest_of):
        """模拟外部计算者：对当前 plan 里每个 ready 项拿票据确认一轮。"""
        for key in plan["ready"]:
            node, part = key.split(":", 1)
            self.confirm(node, part, out_digest_of(node, part))


def main():
    workdir = tempfile.mkdtemp(prefix="vcfpipe_replan_demo_")
    print(f"工作目录（合成账本，结束后可自行删除）: {workdir}")
    d = Demo(workdir)

    feb = [f"2024-02-{i:02d}" for i in range(1, 30)]          # 闰年 29 天
    extra_days = ["2024-03-01", "2024-03-02"]
    all_days = feb + extra_days
    targets = ([f"rolling:{x}" for x in
                ["2024-02-28", "2024-02-29", "2024-03-01", "2024-03-02"]]
               + ["monthly:2024-02", "qc:2024-02-28"])
    chain_targets = targets[:-1]  # 不含 qc 的主链目标

    clean_out = {day: dig(f"clean-out:{day}") for day in all_days}
    rolling_out = {day: dig(f"rolling-out:{day}") for day in all_days}
    monthly_out = dig("monthly-out:2024-02")
    qc_out = dig("qc-out:2024-02-28")

    # ------------------------------------------------------------------ #
    step("第 1 步：登记来源 → 逐轮规划/确认 → 全部算完后无任务待做")
    for day in all_days:
        d.register(day, f"raw-v1:{day}")
    print("  已登记 raw 分区 2024-02-01..2024-02-29、2024-03-01、2024-03-02")
    for rnd in range(1, 6):
        plan = d.plan(targets)
        print(f"  第 {rnd} 轮规划: ready={len(plan['ready'])} "
              f"summary={plan['summary']}")
        if not plan["ready"]:
            break
        d.compute_ready(plan, lambda n, p: (
            clean_out[p] if n == "clean" else
            rolling_out[p] if n == "rolling" else
            monthly_out if n == "monthly" else qc_out))
    plan = d.plan(targets)
    Demo.brief(plan)
    check(plan["ready"] == [] and
          set(plan["summary"]) <= {"present", "current"},
          "全部算完后无任务待做（闭包内只有 present/current）")

    # ------------------------------------------------------------------ #
    step("第 2 步：补交 2 月 28 日 → 只影响窗口内恰当的 rolling 与月表")
    d.register("2024-02-28", "raw-v2:2024-02-28")
    plan = d.plan(targets)
    Demo.brief(plan)
    ready_set = set(plan["ready"])
    waiting = {f"{it['node']}:{it['partition']}"
               for it in plan["items"] if it["status"] == "waiting"}
    check(ready_set == {"clean:2024-02-28", "qc:2024-02-28"},
          f"只有 clean[2024-02-28] 与依赖 raw 的 qc 可开工: {sorted(ready_set)}")
    check(waiting == {"rolling:2024-02-28", "rolling:2024-02-29",
                      "rolling:2024-03-01", "monthly:2024-02"},
          "窗口 [d-2,d] 含 02-28 的 rolling(02-28/02-29/03-01) 与月表转为等待")
    it = Demo.item(plan, "rolling", "2024-03-02")
    check(it and it["status"] == "current",
          "rolling[2024-03-02]（窗口 02-29..03-02）不受影响，仍 current")
    it = Demo.item(plan, "clean", "2024-02-28")
    diff = it["reason"]["input_diffs"][0]
    check(it["reason"]["code"] == "inputs_changed"
          and diff["node"] == "raw" and diff["partition"] == "2024-02-28",
          "ready 原因给出具体差异：raw[2024-02-28] 摘要 旧→新")

    # ------------------------------------------------------------------ #
    step("第 3 步：clean 重算交回相同摘要 → 下游直接恢复 current")
    d.confirm("clean", "2024-02-28", clean_out["2024-02-28"])  # 与以前相同
    d.confirm("qc", "2024-02-28", qc_out)
    plan = d.plan(targets)
    Demo.brief(plan)
    check(plan["ready"] == [] and
          set(plan["summary"]) <= {"present", "current"},
          "clean 输出摘要未变 → rolling/monthly 输入实际没变，无须重跑")

    # ------------------------------------------------------------------ #
    step("第 4 步：clean 重算交回不同摘要 → 下游进入 ready")
    d.register("2024-02-28", "raw-v3:2024-02-28")
    plan = d.plan(targets)
    clean_out_v2 = dig("clean-out-v2:2024-02-28")
    d.confirm("clean", "2024-02-28", clean_out_v2)             # 与以前不同
    d.confirm("qc", "2024-02-28", qc_out)
    plan = d.plan(targets)
    Demo.brief(plan)
    ready_set = set(plan["ready"])
    check(ready_set == {"rolling:2024-02-28", "rolling:2024-02-29",
                        "rolling:2024-03-01", "monthly:2024-02"},
          f"clean 摘要变化后进入 ready 的恰为窗口内 rolling 与月表: {sorted(ready_set)}")
    it = Demo.item(plan, "rolling", "2024-03-02")
    check(it and it["status"] == "current", "rolling[2024-03-02] 仍不受影响")

    # ------------------------------------------------------------------ #
    step("第 5 步：票据发出后改源 → 旧确认被拒绝，旧回执不变")
    old_receipt = d.receipt_digest("rolling", "2024-02-28")
    kept = os.path.join(workdir, "kept_ticket_rolling_2024-02-28.json")
    shutil.copyfile(os.path.join(
        d.tickets, "ticket_rolling_2024-02-28.json"), kept)
    d.register("2024-02-28", "raw-v4:2024-02-28")   # 拿到票据后源数据又改了
    rej = d.confirm("rolling", "2024-02-28", rolling_out["2024-02-28"],
                    ticket_path=kept, expect=1)
    print(f"  拒绝原因: {rej['message']}")
    check(rej.get("cause") == "input_invalid"
          and rej["input"]["node"] == "clean",
          "拒收原因点名失效的直接输入 clean[2024-02-28]")
    check(d.receipt_digest("rolling", "2024-02-28") == old_receipt,
          "旧回执未被改动")

    # ------------------------------------------------------------------ #
    step("第 6 步：票据发出后改 recipe → 旧确认被拒绝，旧回执不变")
    d.register("2024-02-28", "raw-v3:2024-02-28")   # 源恢复，clean 重新有效
    plan = d.plan(targets)
    kept2 = os.path.join(workdir, "kept_ticket2_rolling_2024-02-28.json")
    shutil.copyfile(os.path.join(
        d.tickets, "ticket_rolling_2024-02-28.json"), kept2)
    d.write_config(rolling_recipe="rolling-v2", qc_recipe="qc-v1")  # 改 recipe
    rej = d.confirm("rolling", "2024-02-28", rolling_out["2024-02-28"],
                    ticket_path=kept2, expect=1)
    print(f"  拒绝原因: {rej['message']}")
    check(rej.get("cause") == "rule_changed", "拒收原因为规则已变化")
    check(d.receipt_digest("rolling", "2024-02-28") == old_receipt,
          "旧回执仍未被改动")
    plan = d.plan(targets)
    it = Demo.item(plan, "rolling", "2024-02-28")
    check(it["reason"]["code"] == "rule_changed"
          and it["reason"]["recipe_changed"] is True,
          "重新规划指出 recipe_changed=true")
    rolling_out_v2 = {day: dig(f"rolling-out-v2:{day}") for day in all_days}
    monthly_out_v2 = dig("monthly-out-v2:2024-02")
    # recipe 变了，rolling 全部分区（含未受源影响的 03-02）都需按新规则重算
    for day in ["2024-02-28", "2024-02-29", "2024-03-01", "2024-03-02"]:
        d.confirm("rolling", day, rolling_out_v2[day])
    d.confirm("monthly", "2024-02", monthly_out_v2)
    plan = d.plan(targets)
    check(plan["ready"] == [] and
          set(plan["summary"]) <= {"present", "current"},
          "按新规则确认后全部恢复 current")

    # ------------------------------------------------------------------ #
    step("第 7 步：只调整无关支线（qc 的 recipe）不引起重算")
    d.write_config(rolling_recipe="rolling-v2", qc_recipe="qc-v2")
    plan = d.plan(chain_targets)
    check(plan["ready"] == [] and
          set(plan["summary"]) <= {"present", "current"},
          "主链目标闭包内全部 current，qc 的调整没有波及")
    plan_qc = d.plan(["qc:2024-02-28"])
    it = Demo.item(plan_qc, "qc", "2024-02-28")
    check(it["status"] == "ready" and it["reason"]["recipe_changed"] is True,
          "失效是精确的：只有 qc 自己因 recipe 变化需要重算")
    d.confirm("qc", "2024-02-28", dig("qc-out-v2:2024-02-28"))

    # ------------------------------------------------------------------ #
    step("第 8 步：撤回月末一天 → 阻塞链点名那一天；补回后可以继续")
    d.retract("2024-02-29")
    plan = d.plan(chain_targets)
    Demo.brief(plan)
    it = Demo.item(plan, "monthly", "2024-02")
    chains = [[f"{x['node']}[{x['partition']}]" for x in m["chain"]]
              for m in it["reason"]["missing"]]
    print(f"  monthly[2024-02] 缺失链: {chains}")
    check(it["status"] == "blocked" and
          any(c[-1] == "raw[2024-02-29]" for c in chains),
          "阻塞链点名缺失的正是 raw[2024-02-29]")
    for day, expect in [("2024-02-28", "current"), ("2024-02-29", "blocked"),
                        ("2024-03-01", "blocked"), ("2024-03-02", "blocked")]:
        it = Demo.item(plan, "rolling", day)
        check(it and it["status"] == expect,
              f"rolling[{day}] 状态为 {expect}")
    d.register("2024-02-29", "raw-v1:2024-02-29")   # 补回（内容相同）
    plan = d.plan(chain_targets)
    check(plan["ready"] == [] and
          set(plan["summary"]) <= {"present", "current"},
          "补回内容相同的分区后，整条链直接恢复 current，可以继续")

    # ------------------------------------------------------------------ #
    step("第 9 步：闰年 2 月需要 29 天，普通年只需 28 天")
    feb2023 = [f"2023-02-{i:02d}" for i in range(1, 29)]
    for day in feb2023:
        d.register(day, f"raw-v1:{day}")
    p2024 = d.plan(["monthly:2024-02"])
    p2023 = d.plan(["monthly:2023-02"])
    n2024 = sum(1 for it in p2024["items"] if it["node"] == "clean")
    n2023 = sum(1 for it in p2023["items"] if it["node"] == "clean")
    print(f"  monthly[2024-02] 闭包内 clean 日分区数 = {n2024}")
    print(f"  monthly[2023-02] 闭包内 clean 日分区数 = {n2023}")
    check(n2024 == 29, "闰年 2024-02 展开为 29 个日分区（含 2024-02-29）")
    check(n2023 == 28, "普通年 2023-02 展开为 28 个日分区")
    bad = d.run(["register", "--config", d.config, "--ledger", d.ledger,
                 "--node", "raw", "--partition", "2023-02-29",
                 "--sha256", dig("x")], expect=2)
    print("  登记 2023-02-29 被拒绝（普通年没有这一天）✓")

    # ------------------------------------------------------------------ #
    step("第 10 步：循环定义与不匹配的粒度被实际拒绝")
    bad_cfgs = {
        "cycle": {"nodes": [
            {"name": "a", "type": "derived", "granularity": "day",
             "recipe": "a", "dependencies": [{"node": "b", "map": "identity"}]},
            {"name": "b", "type": "derived", "granularity": "day",
             "recipe": "b", "dependencies": [{"node": "a", "map": "identity"}]},
        ]},
        "gran_identity": {"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "m", "type": "derived", "granularity": "month",
             "recipe": "m", "dependencies": [{"node": "s", "map": "identity"}]},
        ]},
        "gran_window": {"nodes": [
            {"name": "s", "type": "source", "granularity": "day"},
            {"name": "m", "type": "derived", "granularity": "month",
             "recipe": "m",
             "dependencies": [{"node": "s", "map": "day_window",
                               "lo": -1, "hi": 0}]},
        ]},
        "gran_month_days": {"nodes": [
            {"name": "s", "type": "source", "granularity": "month"},
            {"name": "d", "type": "derived", "granularity": "day",
             "recipe": "d",
             "dependencies": [{"node": "s", "map": "month_days"}]},
        ]},
    }
    for name, cfg in bad_cfgs.items():
        path = os.path.join(workdir, f"bad_{name}.json")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(cfg, fh)
        proc = subprocess.run(
            [sys.executable, "-I", RUNNER, "replan", "load",
             "--config", path], capture_output=True, text=True)
        check(proc.returncode == 2, f"配置 {name} 被拒绝（退出码 2）")
        print(f"    {name}: {proc.stderr.strip()}")
    proc = subprocess.run(
        [sys.executable, "-I", RUNNER, "replan", "plan",
         "--config", d.config, "--ledger", d.ledger,
         "--target", "rolling:0001-01-01"], capture_output=True, text=True)
    check(proc.returncode == 2 and "超界" in proc.stderr,
          "日期运算超界被解释并拒绝（rolling[0001-01-01] 窗口越出公历）")
    print(f"    overflow: {proc.stderr.strip()}")

    # ------------------------------------------------------------------ #
    step("第 11 步：独立朴素核对器复核最终账本")
    plan_json = os.path.join(workdir, "final_plan.json")
    d.plan(chain_targets, keep_json=plan_json)
    proc = subprocess.run(
        [sys.executable, "-I", CHECKER, "--config", d.config,
         "--ledger", d.ledger, "--plan", plan_json],
        capture_output=True, text=True)
    print("  " + proc.stdout.strip().replace("\n", "\n  "))
    check(proc.returncode == 0, "独立核对器与产品规划结果一致")

    # ------------------------------------------------------------------ #
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"演示发现 {len(FAILURES)} 处不符合预期：")
        for m in FAILURES:
            print(f"  - {m}")
        return 1
    print("演示完成：全部检查通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
