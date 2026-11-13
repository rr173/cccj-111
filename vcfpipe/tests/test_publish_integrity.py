"""发布完整性回归：版本冲突、中断与并发下的整批结果可见性。

这些测试全部通过**子进程**驱动 CLI（``run_pipeline.py``），模拟真实
使用场景，而不是只断言异常类型：

1. 发布后改动真实 GT/AD 再重跑 -> 拒绝，且拒绝前后旧目录逐文件字节
   一致、查询结果不变；
2. 输出写入中途被 SIGKILL -> 不留半成品/残留锁，重启构建可成功；
3. 相同/不同来源内容并发竞争同一目的版本 -> 只有一个完整结果成为
   发布版本，最终清单、矩阵与 VCF 同属一批结果；
4. 已有分区目录与另一规则版本的查询始终可用。

为什么旧回归没抓住这个缺陷
--------------------------
旧测试 ``test_version_conflict_on_tampered_input`` 有两个盲区：

* 它只断言“抛出版本冲突”，从不在异常之后重新读取数据集。旧实现
  先 ``write_outputs`` 覆盖已发布目录、后比对清单，因此“异常照抛、
  数据已毁”的组合根本无法被该断言观测到；
* 它使用的扰动是往样本表末尾追加一个空行——只改变清单里的样本表
  摘要，不改变任何输出文件内容。即便旧目录被整体覆盖，覆盖后的
  字节也与旧内容完全相同，破坏本身不可见。

本文件的回归改用真实 GT/AD 变更（输出内容必然不同），在拒绝后逐
文件比对字节并重新查询；同时用 SIGKILL 与并发子进程覆盖提交边界
的中断与竞争语义。
"""
import csv
import fcntl
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SRC = os.path.join(ROOT, "src")
RUNNER = os.path.join(ROOT, "run_pipeline.py")
DATA = os.path.join(ROOT, "data")

sys.path.insert(0, SRC)
from vcfpipe.query import Dataset  # noqa: E402

# 与 docs/README.md 一致的相对路径调用方式（cwd = 工作目录）
RUN_ARGS = ["run", "--fasta", "data/ref/syn37.fa",
            "--contigs", "data/ref/contigs.json",
            "--samples", "data/input/sample_sheet.csv",
            "--priority", "data/input/caller_priority.json"]

# 写入中途自杀的故障注入包装器：第 __KILL_AT__ 次 _write_tsv 时 SIGKILL
# 自己，模拟“输出写到一半进程被终止”。占位符用 replace 注入，避免
# 与代码中的花括号冲突。
KILL_WRAPPER = r'''
import os, signal, sys
sys.path.insert(0, "__SRC__")
from vcfpipe import pipeline as p

_orig_write_tsv = p._write_tsv
_calls = [0]

def _killing_write_tsv(path, header, rows):
    _calls[0] += 1
    if _calls[0] == __KILL_AT__:
        os.kill(os.getpid(), signal.SIGKILL)
    return _orig_write_tsv(path, header, rows)

p._write_tsv = _killing_write_tsv
p.run_pipeline("__OUT_ROOT__", "data/ref/syn37.fa", "data/ref/contigs.json",
               "data/input/sample_sheet.csv",
               "data/input/caller_priority.json", "v1")
'''


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    return sha256_bytes(read_bytes(path))


def read_bytes(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def snapshot_tree(directory: str) -> dict:
    """目录下全部文件的 {相对路径: 字节} 快照，用于逐文件比对。"""
    snap = {}
    for dirpath, _dirs, files in os.walk(directory):
        for name in files:
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, directory)
            with open(full, "rb") as fh:
                snap[rel] = fh.read()
    return snap


class PublishIntegrity(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="vcfpipe_pub_")
        shutil.copytree(DATA, os.path.join(self.work, "data"))
        self.out_root = os.path.join(self.work, "out")
        self.dataset_dir = os.path.join(self.out_root, "syn37_v1")

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    # ------------------------------------------------------------------ #
    # 子进程助手
    # ------------------------------------------------------------------ #
    def cli(self, cwd, *args):
        return subprocess.run(
            [sys.executable, "-I", RUNNER, *args],
            cwd=cwd, capture_output=True, text=True, timeout=120)

    def build_v1(self, cwd=None, out_root=None, rule="v1"):
        cwd = cwd or self.work
        out_root = out_root or "out"
        return self.cli(cwd, *(RUN_ARGS + ["--out-root", out_root,
                                           "--rule-version", rule]))

    def query_range(self, dataset, chrom, start, end, cwd=None):
        r = self.cli(cwd or self.work, "query", "--dataset", dataset,
                     "--mode", "range", "--chrom", chrom,
                     "--start", str(start), "--end", str(end))
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def tamper_deepvar_chr1_5(self, vcf_path):
        """真实改动来源记录：chr1:5 的 GT 0/1 -> 1/1，AD 12,8 -> 12,9。"""
        with open(vcf_path, encoding="utf-8") as fh:
            text = fh.read()
        old = "chr1\t5\t.\tA\tG\t38.2\tPASS\t.\tGT:AD:DP\t0/1:12,8:20"
        new = "chr1\t5\t.\tA\tG\t38.2\tPASS\t.\tGT:AD:DP\t1/1:12,9:20"
        self.assertEqual(text.count(old), 1, "篡改锚点必须唯一命中")
        with open(vcf_path, "w", encoding="utf-8", newline="") as fh:
            fh.write(text.replace(old, new))

    def assert_manifest_consistent(self, dataset_dir):
        """清单摘要与磁盘文件一一对应：整批结果同源，无撕裂、无多余文件。"""
        with open(os.path.join(dataset_dir, "manifest.json"),
                  encoding="utf-8") as fh:
            manifest = json.load(fh)
        outputs = manifest["outputs"]
        for rel, digest in outputs.items():
            self.assertEqual(
                sha256_file(os.path.join(dataset_dir, rel)), digest,
                f"{rel} 与清单摘要不符（可能混入了不同批次的内容）")
        on_disk = set()
        for dirpath, _dirs, files in os.walk(dataset_dir):
            for name in files:
                on_disk.add(os.path.relpath(
                    os.path.join(dirpath, name), dataset_dir))
        self.assertEqual(on_disk, set(outputs) | {"manifest.json"},
                         "数据集目录含清单未声明的文件（半成品残留）")
        return manifest

    def assert_no_staging_left(self):
        for name in os.listdir(self.out_root):
            self.assertNotIn("staging", name,
                             f"staging 残留未被清理: {name}")

    def assert_lock_free(self, dataset_id="syn37_v1"):
        """进程退出后不得留下永久占用：测试进程应能立即独占该锁。"""
        lock_path = os.path.join(self.out_root, f".{dataset_id}.lock")
        self.assertTrue(os.path.exists(lock_path))
        fd = os.open(lock_path, os.O_RDWR)
        try:  # 若锁仍被占用，LOCK_NB 会抛 BlockingIOError
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        finally:
            os.close(fd)

    def write_kill_wrapper(self, kill_at, out_root):
        path = os.path.join(self.work, f"kill_build_{kill_at}.py")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(KILL_WRAPPER
                     .replace("__SRC__", SRC)
                     .replace("__OUT_ROOT__", out_root)
                     .replace("__KILL_AT__", str(kill_at)))
        return path

    def run_kill_wrapper(self, kill_at, out_root=None):
        wrapper = self.write_kill_wrapper(kill_at, out_root or self.out_root)
        return subprocess.run([sys.executable, "-I", wrapper],
                              cwd=self.work, capture_output=True, text=True,
                              timeout=120)

    # ------------------------------------------------------------------ #
    # 1. 发布后改动真实 GT/AD 再重跑：拒绝且旧内容分毫不动
    # ------------------------------------------------------------------ #
    def test_conflict_rerun_preserves_files_and_queries(self):
        r = self.build_v1()
        self.assertEqual(r.returncode, 0, r.stderr)

        before = snapshot_tree(self.dataset_dir)
        q_before = self.query_range("out/syn37_v1", "chr1", 5, 5)
        gt_before = {row["sample"]: row["gt"] for row in q_before}
        self.assertEqual(gt_before["S1"], "0/1")  # deepvar 胜出

        # 改动真实 GT/AD 后向原输出位置重跑 -> 版本冲突（退出码 2）
        self.tamper_deepvar_chr1_5(
            os.path.join(self.work, "data", "input",
                         "batch_b1_deepvar_s1.vcf"))
        r = self.build_v1()
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("版本冲突", r.stderr)

        # 拒绝后旧目录逐文件字节一致（含 manifest.json 与全部分区/矩阵/VCF）
        self.assertEqual(snapshot_tree(self.dataset_dir), before)
        # 查询链返回拒绝前的结果：仍读到旧的 0/1，而不是篡改后的 1/1
        self.assertEqual(self.query_range("out/syn37_v1", "chr1", 5, 5),
                         q_before)
        # 清单摘要与磁盘文件仍然同源
        self.assert_manifest_consistent(self.dataset_dir)
        self.assert_no_staging_left()
        self.assert_lock_free()

    def test_identical_rerun_succeeds_and_is_noop(self):
        self.assertEqual(self.build_v1().returncode, 0)
        before = snapshot_tree(self.dataset_dir)
        q_before = self.query_range("out/syn37_v1", "chr1", 1, 106)

        # 相同输入再次执行：成功且已发布目录不被触碰
        r = self.build_v1()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(snapshot_tree(self.dataset_dir), before)
        self.assertEqual(self.query_range("out/syn37_v1", "chr1", 1, 106),
                         q_before)
        self.assert_no_staging_left()

    # ------------------------------------------------------------------ #
    # 2. 写入中途被 SIGKILL：无半成品、无残留占用，重启可成功
    # ------------------------------------------------------------------ #
    def test_kill_mid_write_first_build_leaves_no_partial_dataset(self):
        r = self.run_kill_wrapper(kill_at=3)  # 第 3 个 TSV 写入时自杀
        self.assertEqual(r.returncode, -signal.SIGKILL, r.stderr)

        # 首次构建被终止：不能留下看似可查询的半成品目录
        self.assertFalse(os.path.exists(self.dataset_dir),
                         "被终止的首次构建留下了部分发布目录")
        # 进程退出后锁被内核回收，不阻塞后续运行
        self.assert_lock_free()

        # 重新启动构建：成功、完整、可查询；staging 残留被清理
        r = self.build_v1()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assert_manifest_consistent(self.dataset_dir)
        rows = self.query_range("out/syn37_v1", "chr1", 5, 5)
        self.assertEqual({x["sample"]: x["gt"] for x in rows}["S1"], "0/1")
        self.assert_no_staging_left()

    def test_kill_mid_write_rebuild_preserves_published_dataset(self):
        self.assertEqual(self.build_v1().returncode, 0)
        before = snapshot_tree(self.dataset_dir)

        # 重跑（相同输入）在写入中途被终止：已发布目录必须分毫不动
        r = self.run_kill_wrapper(kill_at=4)
        self.assertEqual(r.returncode, -signal.SIGKILL, r.stderr)
        self.assertEqual(snapshot_tree(self.dataset_dir), before)
        self.assert_lock_free()

        # 再次正常重跑：幂等成功，内容仍与最初发布一致
        r = self.build_v1()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(snapshot_tree(self.dataset_dir), before)
        self.assert_manifest_consistent(self.dataset_dir)
        self.assert_no_staging_left()

    def test_failed_first_build_leaves_no_dataset_dir(self):
        # 输入非法（删掉强制头声明）-> 构建失败，不得留下半成品
        vcf = os.path.join(self.work, "data", "input",
                           "batch_b1_deepvar_s1.vcf")
        with open(vcf, encoding="utf-8") as fh:
            lines = fh.readlines()
        with open(vcf, "w", encoding="utf-8", newline="") as fh:
            fh.writelines(l for l in lines
                          if not l.startswith("##ASSEMBLY"))
        r = self.build_v1()
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(os.path.exists(self.dataset_dir))
        self.assert_no_staging_left()
        self.assert_lock_free()

    # ------------------------------------------------------------------ #
    # 3. 并发竞争同一目的版本
    # ------------------------------------------------------------------ #
    def _make_input_variant(self, name, tamper=False):
        variant_dir = os.path.join(self.work, name)
        shutil.copytree(DATA, os.path.join(variant_dir, "data"))
        if tamper:
            self.tamper_deepvar_chr1_5(os.path.join(
                variant_dir, "data", "input", "batch_b1_deepvar_s1.vcf"))
        return variant_dir

    def _race(self, cwd_a, cwd_b, race_out):
        """两个进程同时向 race_out 构建同一 dataset_id。"""
        args = RUN_ARGS + ["--out-root", race_out, "--rule-version", "v1"]
        pa = subprocess.Popen([sys.executable, "-I", RUNNER, *args],
                              cwd=cwd_a, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        pb = subprocess.Popen([sys.executable, "-I", RUNNER, *args],
                              cwd=cwd_b, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
        out_a, err_a = pa.communicate(timeout=120)
        out_b, err_b = pb.communicate(timeout=120)
        return (pa.returncode, err_a), (pb.returncode, err_b)

    def test_concurrent_identical_content_both_confirm(self):
        """相同内容并发：两个进程都应确认成功，结果与单独构建逐字节一致。"""
        for i in range(3):  # 多轮提高交错概率；锁保证结果确定
            race_out = os.path.join(self.work, f"race_same_{i}")
            (rc_a, err_a), (rc_b, err_b) = self._race(
                self.work, self.work, race_out)
            self.assertEqual(rc_a, 0, err_a)
            self.assertEqual(rc_b, 0, err_b)
            dataset = os.path.join(race_out, "syn37_v1")
            self.assert_manifest_consistent(dataset)
            # 与单独构建的内容逐字节一致
            solo = os.path.join(self.work, f"solo_same_{i}")
            self.assertEqual(self.build_v1(out_root=solo).returncode, 0)
            self.assertEqual(
                snapshot_tree(dataset),
                snapshot_tree(os.path.join(solo, "syn37_v1")))

    def test_concurrent_different_content_single_winner(self):
        """不同内容并发：只有一个完整结果成为发布版本。"""
        dir_a = self._make_input_variant("a", tamper=False)   # 0/1
        dir_b = self._make_input_variant("b", tamper=True)    # 1/1

        # 各自单独构建，取得两份“合法批次”的清单字节
        solo_a = os.path.join(self.work, "solo_a")
        solo_b = os.path.join(self.work, "solo_b")
        self.assertEqual(self.build_v1(cwd=dir_a, out_root=solo_a)
                         .returncode, 0)
        self.assertEqual(self.build_v1(cwd=dir_b, out_root=solo_b)
                         .returncode, 0)
        man_a = read_bytes(os.path.join(solo_a, "syn37_v1", "manifest.json"))
        man_b = read_bytes(os.path.join(solo_b, "syn37_v1", "manifest.json"))
        self.assertNotEqual(man_a, man_b)

        # 竞争同一目的版本：恰一个胜出（rc 0），另一个被拒绝（rc 2）
        race_out = os.path.join(self.work, "race_diff")
        (rc_a, err_a), (rc_b, err_b) = self._race(dir_a, dir_b, race_out)
        self.assertEqual(sorted([rc_a, rc_b]), [0, 2],
                         f"rc_a={rc_a} rc_b={rc_b}\n{err_a}\n{err_b}")

        # 最终清单必须逐字节等于某一个合法批次，不能是两者的混合
        final = read_bytes(os.path.join(race_out, "syn37_v1",
                                        "manifest.json"))
        self.assertIn(final, (man_a, man_b),
                      "发布清单不属于任何一个竞争批次（撕裂写入）")
        a_won = final == man_a
        self.assertEqual(rc_a if a_won else rc_b, 0)
        loser_err = err_b if a_won else err_a
        self.assertIn("版本冲突", loser_err)

        # 清单摘要、矩阵与 VCF 同属胜出批次；查询与之一致
        dataset = os.path.join(race_out, "syn37_v1")
        self.assert_manifest_consistent(dataset)
        rows = self.query_range(dataset, "chr1", 5, 5)
        self.assertEqual({x["sample"]: x["gt"] for x in rows}["S1"],
                         "0/1" if a_won else "1/1")
        # 胜出内容可被相同输入再次确认（幂等）
        winner_dir = dir_a if a_won else dir_b
        r = self.build_v1(cwd=winner_dir, out_root=race_out)
        self.assertEqual(r.returncode, 0, r.stderr)

    # ------------------------------------------------------------------ #
    # 4. 已有分区目录与另一规则版本的查询继续可用
    # ------------------------------------------------------------------ #
    def test_sibling_rule_version_and_old_partitions_remain_queryable(self):
        self.assertEqual(self.build_v1().returncode, 0)
        q_v1_before = self.query_range("out/syn37_v1", "chr1", 1, 106)

        # v1 冲突重跑被拒后，v2（另一规则版本 = 另一 dataset_id）照常构建
        self.tamper_deepvar_chr1_5(os.path.join(
            self.work, "data", "input", "batch_b1_deepvar_s1.vcf"))
        self.assertEqual(self.build_v1().returncode, 2)
        r = self.build_v1(rule="v2")
        self.assertEqual(r.returncode, 0, r.stderr)

        # v1 仍是拒绝前的旧内容；v2 独立可查（边界记录在 v2 被纳入）
        self.assertEqual(self.query_range("out/syn37_v1", "chr1", 1, 106),
                         q_v1_before)
        v2_rows = self.query_range("out/syn37_v2", "chr1", 1, 1)
        self.assertTrue(any(row["pos"] == "1" for row in v2_rows),
                        "v2 应包含 v1 被隔离的边界记录")
        self.assert_manifest_consistent(
            os.path.join(self.out_root, "syn37_v2"))

    @unittest.skipUnless(
        os.path.isdir(os.path.join(ROOT, "output", "syn37_v1")),
        "仓库未携带预生成 output/ 目录")
    def test_checked_in_datasets_remain_queryable(self):
        """仓库中已生成的分区目录（修复前写出）读取方式不变。"""
        for dataset_id in ("syn37_v1", "syn37_v2"):
            ds = Dataset(os.path.join(ROOT, "output", dataset_id))
            self.assertEqual(ds.dataset_id, dataset_id)
            rows = ds.query_range("chr1", 28, 48)
            self.assertTrue(rows, f"{dataset_id} 已有分区应仍可查询")
            self.assertTrue(ds.sample_slice("S2"))
            self.assertTrue(ds.matrix_wide())
            # 清单摘要与磁盘文件一致（历史目录本身未被破坏）
            for rel, digest in ds.manifest["outputs"].items():
                self.assertEqual(
                    sha256_file(os.path.join(ds.dir, rel)), digest)


if __name__ == "__main__":
    unittest.main(verbosity=2)
