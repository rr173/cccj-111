#!/usr/bin/env python3
"""运行注释演示并保存回归记录。

依次执行：
1. 演示材料 -> ``output/annot_demo``（annotation.tsv + summary.json）；
2. 独立核对器逐行比对；
3. 失败路径（同位点多 ALT 已在演示内；REF 错配、错误转录本定义）
   分别运行并保存退出码与 stderr 到 ``output/annot_regression/``。

退出码 0 表示演示成功、核对一致，且两类错误输入都以非零状态被拒绝。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "annot_independent_check.py")
A = os.path.join(ROOT, "data", "annot")
FA = os.path.join(ROOT, "data", "ref", "annot_syn.fa")
OUT = os.path.join(ROOT, "output", "annot_demo")
REG = os.path.join(ROOT, "output", "annot_regression")


def run(args, **kw):
    return subprocess.run([sys.executable, "-I", *args], cwd=ROOT,
                          capture_output=True, text=True, timeout=120, **kw)


def main() -> int:
    shutil.rmtree(OUT, ignore_errors=True)
    shutil.rmtree(REG, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    demo = run([RUNNER, "annot",
                "--fasta", FA,
                "--transcripts", os.path.join(A, "transcripts_demo.tsv"),
                "--variants", os.path.join(A, "variants_demo.tsv"),
                "--out-dir", OUT])
    if demo.returncode != 0:
        print(demo.stderr, file=sys.stderr)
        return 1

    check = run([CHECKER,
                 "--fasta", FA,
                 "--transcripts", os.path.join(A, "transcripts_demo.tsv"),
                 "--variants", os.path.join(A, "variants_demo.tsv"),
                 "--out-dir", OUT])
    print(check.stdout.strip())
    if check.returncode != 0:
        print(check.stderr, file=sys.stderr)
        return 1

    os.makedirs(REG, exist_ok=True)
    cases = [
        ("ref_mismatch",
         os.path.join(A, "regression_ref_mismatch.tsv"),
         os.path.join(A, "transcripts_demo.tsv")),
        ("bad_transcript",
         os.path.join(A, "variants_demo.tsv"),
         os.path.join(A, "regression_bad_transcripts.tsv")),
    ]
    ok = True
    for name, vfile, tfile in cases:
        path = os.path.join(REG, f"{name}.stderr.txt")
        r = run([RUNNER, "annot",
                 "--fasta", FA,
                 "--transcripts", tfile,
                 "--variants", vfile,
                 "--out-dir", os.path.join(REG, f"{name}_out")])
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(f"exit_code={r.returncode}\n{r.stderr}")
        if r.returncode == 0:
            print(f"错误输入 {name} 未被拒绝", file=sys.stderr)
            ok = False
        else:
            print(f"{name}: exit={r.returncode} -> {r.stderr.splitlines()[0]}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
