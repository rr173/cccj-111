"""dedup 正文语料清理的手算/集成回归。

覆盖：阈值恰好相等与刚低于、重复片段只算一个、全角与大小写、标点
分词、质量相同的 id 决胜、一个冗余文档命中两个保留项时的相似度与
同分决胜、非传递三篇链、短文本不互删、重复 id 与非法输入拒绝、
附带字段逐字节保留、输入行序不影响结果、CLI 退出码。
"""
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RUNNER = os.path.join(ROOT, "run_pipeline.py")
CHECKER = os.path.join(ROOT, "tools", "dedup_independent_check.py")

sys.path.insert(0, os.path.join(ROOT, "src"))

from vcfpipe.dedup import (  # noqa: E402
    DedupInputError, encode_shingle, length_bounds, parse_threshold,
    prefix_length, run_dedup, shingle_triples, shingles_of, tokenize,
)


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
class Workspace:
    """每个用例一个临时目录，跑完即删。"""

    def __init__(self, test):
        self.dir = tempfile.mkdtemp(prefix="dedup_test_")
        test.addCleanup(shutil.rmtree, self.dir, True)
        self.corpus = os.path.join(self.dir, "corpus.jsonl")
        self.out = os.path.join(self.dir, "out")

    def write_corpus(self, records):
        """records: dict 列表或已序列化的字符串行列表。"""
        lines = []
        for rec in records:
            if isinstance(rec, str):
                lines.append(rec)
            else:
                lines.append(json.dumps(rec, ensure_ascii=False))
        with open(self.corpus, "w", encoding="utf-8", newline="") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))
        return lines

    def run(self, threshold="1/2"):
        return run_dedup(self.corpus, threshold, self.out)

    def dispositions(self):
        with open(os.path.join(self.out, "dispositions.jsonl"),
                  encoding="utf-8") as fh:
            return [json.loads(line) for line in fh]

    def kept_lines(self):
        with open(os.path.join(self.out, "kept.jsonl"), "rb") as fh:
            return fh.read().decode("utf-8").splitlines()

    def by_id(self):
        return {d["id"]: d for d in self.dispositions()}


def rec(rid, quality, text, **extra):
    return {"id": rid, "quality": quality, "text": text, **extra}


# --------------------------------------------------------------------------- #
# 比较口径：分词 / 片段 / 编码
# --------------------------------------------------------------------------- #
class TestTokenize(unittest.TestCase):
    def test_punctuation_and_separators(self):
        # 标点（含全角）全作分隔；连续分隔不产生空 token
        self.assertEqual(tokenize("a,,b..c  d!!!"), ["a", "b", "c", "d"])
        self.assertEqual(tokenize("a，b。c！d"), ["a", "b", "c", "d"])
        self.assertEqual(tokenize("！！！。。。"), [])
        self.assertEqual(tokenize(""), [])
        self.assertEqual(tokenize("  x  "), ["x"])

    def test_fullwidth_and_casefold(self):
        # 全角拉丁 -> ASCII；大小写折叠
        self.assertEqual(tokenize("Ｈｅｌｌｏ Ｗｏｒｌｄ"),
                         ["hello", "world"])
        self.assertEqual(tokenize("Straße Ｘ"), ["strasse", "x"])
        self.assertEqual(tokenize("STRASSE"), ["strasse"])

    def test_letter_and_number_categories(self):
        # L* 与 N* 都是 token 字符；①(NFKC->1)、²(NFKC->2) 先规范化
        self.assertEqual(tokenize("第①版 v2 ²"), ["第1版", "v2", "2"])
        # 下划线是 Pc（连接符），不是 L/N -> 分隔
        self.assertEqual(tokenize("a_b"), ["a", "b"])

    def test_shingles_dedup(self):
        triples = shingle_triples(["na", "na", "na", "na", "na", "bat"])
        self.assertEqual(triples, {("na", "na", "na"),
                                   ("na", "na", "bat")})
        self.assertEqual(shingle_triples(["a", "b"]), set())
        self.assertEqual(shingle_triples([]), set())

    def test_encode_order_matches_tuple_order(self):
        # 编码文本的字典序必须与三元组逐元素字典序一致
        rng = random.Random(20261113)
        alphabet = ["a", "ab", "abc", "b", "z", "0", "09", "é", "条",
                    "模板", "²", "ﬀ", "ｘ"]  # 含全角与预组合字符
        triples = set()
        for _ in range(3000):
            triples.add(tuple(rng.choice(alphabet) for _ in range(3)))
        by_tuple = sorted(triples)
        by_encoded = sorted(encode_shingle(t) for t in triples)
        self.assertEqual([encode_shingle(t) for t in by_tuple], by_encoded)

    def test_prefix_and_bounds(self):
        # 前缀长度 n - ceil(pn/q) + 1
        self.assertEqual(prefix_length(3, 1, 2), 2)
        self.assertEqual(prefix_length(4, 1, 2), 3)
        self.assertEqual(prefix_length(5, 1, 1), 1)
        self.assertEqual(prefix_length(1, 1, 2), 1)
        # 长度界 [ceil(pn/q), floor(qn/p)]
        self.assertEqual(length_bounds(3, 1, 2), (2, 6))
        self.assertEqual(length_bounds(4, 1, 1), (4, 4))
        self.assertEqual(length_bounds(5, 2, 3), (4, 7))


class TestThresholdParsing(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(parse_threshold("1/2"), (1, 2))
        self.assertEqual(parse_threshold("1/1"), (1, 1))
        self.assertEqual(parse_threshold("01/02"), (1, 2))
        self.assertEqual(parse_threshold("3/5"), (3, 5))

    def test_invalid(self):
        for bad in ("0/1", "1/0", "2/1", "0.5", "1/2/3", "abc", "/2",
                    "1/", "1 /2", "1/ 2", "-1/2", "1/-2", "", "1", "½"):
            with self.assertRaises(DedupInputError, msg=bad):
                parse_threshold(bad)


# --------------------------------------------------------------------------- #
# 清理语义（手算样例）
# --------------------------------------------------------------------------- #
class TestCleaningSemantics(unittest.TestCase):
    def test_nontransitive_chain(self):
        # 规格样例：阈值 1/2，质量依次降低的 a..e / b..f / c..g
        ws = Workspace(self)
        ws.write_corpus([
            rec("doc1", 3, "a b c d e"),
            rec("doc2", 2, "b c d e f"),
            rec("doc3", 1, "c d e f g"),
        ])
        summary = ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["doc1"]["disposition"], "kept")
        # 中间一篇归第一篇：交集 2、并集 4，恰好等于阈值也算重复
        self.assertEqual(d["doc2"]["disposition"], "redundant")
        self.assertEqual(d["doc2"]["witness_id"], "doc1")
        self.assertEqual((d["doc2"]["intersection"], d["doc2"]["union"]),
                         (2, 4))
        self.assertEqual(d["doc2"]["similarity"], "2/4")
        # 第三篇与第一篇只有 1/5，必须留下（不能沿链误删）
        self.assertEqual(d["doc3"]["disposition"], "kept")
        self.assertEqual(summary["kept"], 2)
        self.assertEqual(summary["redundant"], 1)

    def test_threshold_exactly_equal_is_duplicate(self):
        ws = Workspace(self)
        # |A|=|B|=3，交集 2 -> 2/4 = 1/2，恰好等于阈值
        ws.write_corpus([
            rec("a", 2, "x1 x2 x3 x4 ra"),
            rec("b", 1, "x1 x2 x3 x4 rb"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["b"]["disposition"], "redundant")
        self.assertEqual(d["b"]["witness_id"], "a")
        self.assertEqual(d["b"]["similarity"], "2/4")

    def test_just_below_threshold_is_kept(self):
        ws = Workspace(self)
        # 交集 1、并集 3 -> 1/3 < 1/2
        ws.write_corpus([
            rec("a", 2, "s1 s2 s3 ua"),
            rec("b", 1, "s1 s2 s3 ub"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["a"]["disposition"], "kept")
        self.assertEqual(d["b"]["disposition"], "kept")
        # 保留行说明在先保留项未达阈值：候选中最高相似度 1/3
        best = d["b"]["best_below_threshold"]
        self.assertIsNotNone(best)
        self.assertEqual(best["id"], "a")
        self.assertEqual(best["similarity"], "1/3")

    def test_duplicate_shingles_counted_once(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("a", 2, "na na na na na batman"),   # 片段集合 2 项
            rec("b", 1, "na na na batman"),          # 同样的 2 项
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["a"]["shingle_count"], 2)
        self.assertEqual(d["b"]["shingle_count"], 2)
        self.assertEqual(d["b"]["disposition"], "redundant")
        self.assertEqual(d["b"]["similarity"], "2/2")

    def test_fullwidth_and_casefold_dedup(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("fw", 2, "Ｈｅｌｌｏ Ｗｏｒｌｄ Ｆｏｏ Ｂａｒ"),
            rec("plain", 1, "HELLO WORLD FOO BAR"),
            rec("esz1", 2, "Straße x y"),
            rec("esz2", 1, "STRASSE x y"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["plain"]["disposition"], "redundant")
        self.assertEqual(d["plain"]["witness_id"], "fw")
        self.assertEqual(d["esz2"]["disposition"], "redundant")
        self.assertEqual(d["esz2"]["witness_id"], "esz1")

    def test_punctuation_variant_dedup(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("a", 2, "red, green, blue; yellow"),
            rec("b", 1, "red。。green！！blue……yellow"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["b"]["disposition"], "redundant")
        self.assertEqual(d["b"]["witness_id"], "a")
        self.assertEqual(d["b"]["similarity"], "2/2")

    def test_quality_desc_then_id_order(self):
        ws = Workspace(self)
        # 同质量同文：id 字典序在先者保留
        ws.write_corpus([
            rec("tie-b", 5, "same text here today"),
            rec("tie-a", 5, "same text here today"),
            rec("low", 1, "same text here today"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["tie-a"]["disposition"], "kept")
        self.assertEqual(d["tie-b"]["disposition"], "redundant")
        self.assertEqual(d["tie-b"]["witness_id"], "tie-a")
        self.assertEqual(d["low"]["witness_id"], "tie-a")
        # 处置表顺序：tie-a, tie-b, low
        self.assertEqual([x["id"] for x in ws.dispositions()],
                         ["tie-a", "tie-b", "low"])

    def test_witness_is_highest_similarity(self):
        ws = Workspace(self)
        # R 对 K1 为 2/4，对 K2 为 2/3 -> 见证 K2
        ws.write_corpus([
            rec("K1", 5, "x1 x2 x3 x4 k1"),
            rec("K2", 4, "x2 x3 x4 r"),
            rec("R", 3, "x1 x2 x3 x4 r"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["K1"]["disposition"], "kept")
        self.assertEqual(d["K2"]["disposition"], "kept")
        self.assertEqual(d["R"]["disposition"], "redundant")
        self.assertEqual(d["R"]["witness_id"], "K2")
        self.assertEqual(d["R"]["similarity"], "2/3")

    def test_witness_tie_picks_earliest(self):
        ws = Workspace(self)
        # R 对 K1、K2 都是 2/4 -> 同分取清理顺序最早的 K1
        ws.write_corpus([
            rec("K1", 5, "x1 x2 x3 x4 k1"),
            rec("K2", 4, "x2 x3 x4 r k2"),
            rec("R", 3, "x1 x2 x3 x4 r"),
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["K2"]["disposition"], "kept")
        self.assertEqual(d["R"]["disposition"], "redundant")
        self.assertEqual(d["R"]["witness_id"], "K1")
        self.assertEqual(d["R"]["similarity"], "2/4")

    def test_redundant_not_in_index(self):
        ws = Workspace(self)
        # B 冗余于 A；C 只与 B 达到阈值。B 不入索引 -> C 必须保留
        ws.write_corpus([
            rec("A", 3, "a b c d"),          # 片段 {abc,bcd}
            rec("B", 2, "a b c d x"),        # 与 A: 2/3 >= 1/2 -> 冗余
            rec("C", 1, "b c d x y"),        # 与 B: 2/4；与 A: 1/4 -> 保留
        ])
        ws.run("1/2")
        d = ws.by_id()
        self.assertEqual(d["B"]["disposition"], "redundant")
        self.assertEqual(d["B"]["witness_id"], "A")
        self.assertEqual(d["C"]["disposition"], "kept")

    def test_too_short_never_merged(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("empty-1", 4, ""),
            rec("empty-2", 3, ""),
            rec("two", 2, "just two"),
            rec("punct", 1, "！！！"),
        ])
        summary = ws.run("1/2")
        d = ws.by_id()
        for rid in ("empty-1", "empty-2", "two", "punct"):
            self.assertEqual(d[rid]["disposition"], "too_short", rid)
        self.assertEqual(summary["too_short"], 4)
        self.assertEqual(summary["redundant"], 0)
        # 短文本不进索引：不产生任何候选
        self.assertEqual(summary["candidate_pairs"], 0)
        self.assertEqual(summary["exact_comparisons"], 0)
        # 但原文全部保留
        self.assertEqual(len(ws.kept_lines()), 4)

    def test_threshold_one_requires_identical(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("a", 3, "x1 x2 x3 x4"),
            rec("b", 2, "x1 x2 x3 x4"),      # 完全相同 -> 1/1 冗余
            rec("c", 1, "x1 x2 x3 y9"),      # 2/3 < 1 -> 保留
        ])
        ws.run("1/1")
        d = ws.by_id()
        self.assertEqual(d["b"]["disposition"], "redundant")
        self.assertEqual(d["c"]["disposition"], "kept")


# --------------------------------------------------------------------------- #
# 输出契约
# --------------------------------------------------------------------------- #
class TestOutputContract(unittest.TestCase):
    def test_extra_fields_preserved_byte_identical(self):
        ws = Workspace(self)
        raw = ('{"id":"x","quality":1,"text":"a b c d",'
               '"source":"wβ","tags":["α",2],"meta":{"n":1e2},"u":"①"}')
        ws.write_corpus([raw, rec("y", 2, "e f g h")])
        ws.run("1/2")
        kept = ws.kept_lines()
        self.assertEqual(len(kept), 2)
        self.assertIn(raw, kept)  # 原始行（含附带字段）逐字节读回

    def test_crlf_lines_roundtrip(self):
        # CRLF 行尾：\r 属于原始行内容，kept.jsonl 逐字节保留
        ws = Workspace(self)
        with open(ws.corpus, "wb") as fh:
            fh.write(b'{"id":"a","quality":2,"text":"a b c d"}\r\n')
            fh.write(b'{"id":"b","quality":1,"text":"e f g h"}\r\n')
        ws.run("1/2")
        with open(os.path.join(ws.out, "kept.jsonl"), "rb") as fh:
            data = fh.read()
        self.assertEqual(
            data,
            b'{"id":"a","quality":2,"text":"a b c d"}\r\n'
            b'{"id":"b","quality":1,"text":"e f g h"}\r\n')

    def test_unicode_id_order(self):
        # 同质量按 id 的 Unicode 码点字典序（与 UTF-8 字节序一致）
        ws = Workspace(self)
        ws.write_corpus([
            rec("z", 5, "same text here"),
            rec("é", 5, "same text here"),
            rec("e", 5, "same text here"),
            rec("条", 5, "same text here"),
        ])
        ws.run("1/2")
        order = [d["id"] for d in ws.dispositions()]
        self.assertEqual(order, sorted(["z", "é", "e", "条"]))
        self.assertEqual(ws.dispositions()[0]["disposition"], "kept")
        for d in ws.dispositions()[1:]:
            self.assertEqual(d["disposition"], "redundant")
            self.assertEqual(d["witness_id"], "e")

    def test_shuffle_invariance(self):
        records = [
            rec("chain-a", 30, "alpha beta gamma delta epsilon",
                source="hand"),
            rec("chain-b", 20, "beta gamma delta epsilon zeta"),
            rec("chain-c", 10, "gamma delta epsilon zeta eta"),
            rec("fw-1", 15, "Ｈｅｌｌｏ Ｗｏｒｌｄ Ｆｏｏ Ｂａｒ"),
            rec("fw-2", 15, "HELLO WORLD FOO BAR"),
            rec("tpl-1", 4, "通用 模板 条款 第一 第二 第三 第四 第五 甲"),
            rec("tpl-2", 3, "通用 模板 条款 第一 第二 第三 第四 第五 乙"),
            rec("short-1", 7, ""),
            rec("short-2", 6, "just two"),
            rec("uniq", 2, "orchid maple river quartz lantern"),
        ]
        outputs = []
        for seed in range(4):
            rng = random.Random(seed)
            shuffled = records[:]
            rng.shuffle(shuffled)
            ws = Workspace(self)
            ws.write_corpus(shuffled)
            ws.run("1/2")
            with open(os.path.join(ws.out, "summary.json"),
                      encoding="utf-8") as fh:
                summary = json.load(fh)
            outputs.append((
                [json.dumps(d, sort_keys=True) for d in ws.dispositions()],
                ws.kept_lines(),
                summary,
            ))
        for i in range(1, len(outputs)):
            self.assertEqual(outputs[0], outputs[i],
                             "重排输入后处置表/保留清单/计数必须一致")

    def test_dispositions_cover_every_doc_in_clean_order(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("b", 5, "a b c d"),
            rec("a", 5, "a b c d"),
            rec("s", 9, ""),
        ])
        ws.run("1/2")
        disps = ws.dispositions()
        self.assertEqual([d["id"] for d in disps], ["s", "a", "b"])
        self.assertEqual([d["seq"] for d in disps], [1, 2, 3])
        self.assertEqual([d["disposition"] for d in disps],
                         ["too_short", "kept", "redundant"])

    def test_empty_corpus(self):
        ws = Workspace(self)
        ws.write_corpus([])
        summary = ws.run("1/2")
        self.assertEqual(summary["documents"], 0)
        self.assertEqual(summary["kept"], 0)
        self.assertEqual(summary["candidate_pairs"], 0)
        self.assertEqual(ws.kept_lines(), [])
        self.assertEqual(ws.dispositions(), [])


# --------------------------------------------------------------------------- #
# 输入拒绝：指出位置并拒绝这次作业
# --------------------------------------------------------------------------- #
class TestInputRejection(unittest.TestCase):
    def assert_rejected(self, lines, needle, threshold="1/2"):
        ws = Workspace(self)
        ws.write_corpus(lines)
        with self.assertRaises(DedupInputError) as ctx:
            ws.run(threshold)
        self.assertIn(needle, str(ctx.exception))
        # 拒绝时不创建输出目录（不能貌似成功）
        self.assertFalse(os.path.exists(ws.out))

    def test_duplicate_id(self):
        self.assert_rejected(
            [rec("a", 1, "x y z"), rec("b", 1, "x y z"), rec("a", 2, "p q r")],
            "第 3 行")
        # 报错指出首次出现位置
        ws = Workspace(self)
        ws.write_corpus([rec("a", 1, "x y z"), rec("a", 2, "p q r")])
        with self.assertRaises(DedupInputError) as ctx:
            ws.run()
        self.assertIn("第 2 行", str(ctx.exception))
        self.assertIn("第 1 行", str(ctx.exception))

    def test_bad_json_line(self):
        self.assert_rejected(['{"id":"a","quality":1,"text":"x"}',
                              '{broken json'],
                             "第 2 行")

    def test_blank_line_rejected(self):
        self.assert_rejected(['{"id":"a","quality":1,"text":"x"}', ""],
                             "第 2 行")

    def test_non_object_line(self):
        self.assert_rejected(['["id","a"]'], "第 1 行")

    def test_bad_field_types(self):
        cases = [
            ({"id": "a", "quality": True, "text": "x y z"}, "quality"),
            ({"id": "a", "quality": 1.5, "text": "x y z"}, "quality"),
            ({"id": "a", "quality": 2 ** 63, "text": "x y z"}, "quality"),
            ({"id": "a", "quality": -(2 ** 63) - 1, "text": "x"}, "quality"),
            ({"id": "", "quality": 1, "text": "x y z"}, "id"),
            ({"id": 7, "quality": 1, "text": "x y z"}, "id"),
            ({"quality": 1, "text": "x y z"}, "id"),
            ({"id": "a", "text": "x y z"}, "quality"),
            ({"id": "a", "quality": 1}, "text"),
            ({"id": "a", "quality": 1, "text": None}, "text"),
            ({"id": "a", "quality": 1, "text": 3}, "text"),
        ]
        for obj, needle in cases:
            with self.subTest(obj=obj):
                self.assert_rejected([obj], needle)

    def test_int64_boundaries_accepted(self):
        ws = Workspace(self)
        ws.write_corpus([
            rec("lo", -(2 ** 63), "a b c d"),
            rec("hi", 2 ** 63 - 1, "e f g h"),
        ])
        summary = ws.run("1/2")
        self.assertEqual(summary["documents"], 2)

    def test_multiple_errors_reported_with_positions(self):
        ws = Workspace(self)
        ws.write_corpus([
            "{broken",
            rec("ok", 1, "a b c"),
            "{also broken",
        ])
        with self.assertRaises(DedupInputError) as ctx:
            ws.run()
        msg = str(ctx.exception)
        self.assertIn("第 1 行", msg)
        self.assertIn("第 3 行", msg)
        self.assertIn("共 2 处", msg)

    def test_invalid_threshold_rejected_before_work(self):
        ws = Workspace(self)
        ws.write_corpus([rec("a", 1, "x y z")])
        for bad in ("0/1", "2/1", "1/0", "0.5", "x"):
            with self.assertRaises(DedupInputError):
                ws.run(bad)
        self.assertFalse(os.path.exists(ws.out))

    def test_out_dir_must_be_empty(self):
        ws = Workspace(self)
        ws.write_corpus([rec("a", 1, "x y z")])
        os.makedirs(ws.out)
        with open(os.path.join(ws.out, "junk"), "w") as fh:
            fh.write("x")
        with self.assertRaises(DedupInputError):
            ws.run()


# --------------------------------------------------------------------------- #
# CLI 端到端
# --------------------------------------------------------------------------- #
class TestCli(unittest.TestCase):
    def run_cli(self, *cli_args):
        return subprocess.run(
            [sys.executable, "-I", RUNNER, "dedup", *cli_args],
            capture_output=True, text=True)

    def test_cli_success_and_exit_codes(self):
        work = tempfile.mkdtemp(prefix="dedup_cli_")
        self.addCleanup(shutil.rmtree, work, True)
        corpus = os.path.join(work, "c.jsonl")
        with open(corpus, "w", encoding="utf-8") as fh:
            fh.write('{"id":"a","quality":2,"text":"a b c d"}\n')
            fh.write('{"id":"b","quality":1,"text":"a b c d"}\n')
        out = os.path.join(work, "out")
        proc = self.run_cli("--corpus", corpus, "--threshold", "1/2",
                            "--out-dir", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(proc.stdout)
        self.assertEqual(summary["redundant"], 1)
        self.assertTrue(os.path.isfile(
            os.path.join(out, "dispositions.jsonl")))

        # 重复 id -> 退出码 2，且指出位置
        bad = os.path.join(work, "bad.jsonl")
        with open(bad, "w", encoding="utf-8") as fh:
            fh.write('{"id":"a","quality":1,"text":"x"}\n')
            fh.write('{"id":"a","quality":2,"text":"y"}\n')
        proc = self.run_cli("--corpus", bad, "--threshold", "1/2",
                            "--out-dir", os.path.join(work, "out2"))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("第 2 行", proc.stderr)
        self.assertFalse(os.path.exists(os.path.join(work, "out2")))

        # 非法阈值 -> 退出码 2
        proc = self.run_cli("--corpus", corpus, "--threshold", "0.5",
                            "--out-dir", os.path.join(work, "out3"))
        self.assertEqual(proc.returncode, 2)

    def test_cli_demo_corpus(self):
        work = tempfile.mkdtemp(prefix="dedup_cli_demo_")
        self.addCleanup(shutil.rmtree, work, True)
        proc = self.run_cli(
            "--corpus", os.path.join(ROOT, "data", "corpus",
                                     "demo_small.jsonl"),
            "--threshold", "1/2", "--out-dir", os.path.join(work, "out"))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads(proc.stdout)
        self.assertEqual(
            (summary["documents"], summary["kept"],
             summary["too_short"], summary["redundant"]),
            (21, 12, 3, 6))


class TestIndependentChecker(unittest.TestCase):
    """独立慢速核对器：不调用生产分词/候选/交集函数的全量对照。"""

    def test_checker_single_point_on_demo(self):
        work = tempfile.mkdtemp(prefix="dedup_check_demo_")
        self.addCleanup(shutil.rmtree, work, True)
        corpus = os.path.join(ROOT, "data", "corpus", "demo_small.jsonl")
        out = os.path.join(work, "out")
        proc = subprocess.run(
            [sys.executable, "-I", RUNNER, "dedup", "--corpus", corpus,
             "--threshold", "1/2", "--out-dir", out],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        proc = subprocess.run(
            [sys.executable, "-I", CHECKER, "--corpus", corpus,
             "--out-dir", out, "--threshold", "1/2"],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_checker_self_check_thresholds_and_shuffles(self):
        # 多阈值 × 打乱行序的合成语料整张处置表对照
        proc = subprocess.run(
            [sys.executable, "-I", CHECKER, "--shuffles", "2"],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("核对通过", proc.stdout)


if __name__ == "__main__":
    unittest.main()
