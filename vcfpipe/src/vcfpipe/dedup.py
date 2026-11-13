"""正文语料清理：三词片段 Jaccard 去重 + 磁盘前缀候选检索。

定位
----
采集队列里的转载（加减几个词）与模板正文（多来源重复发送）在交给
后续训练前经本环节去重：逐篇与“已保留”正文比较，达到相似度阈值
的记为 redundant 并指向见证；短文本只留档不参与比较。本环节只做
正文清理，不涉及变异计算；输入为 JSONL 语料，每行一个对象。仅用
Python 标准库 + SQLite（临时库落盘），不依赖外部搜索服务。

输入契约
--------
* 每行一个 JSON 对象（UTF-8）：``id`` 非空字符串且全库唯一；
  ``quality`` 为 64 位有符号整数（不接受 bool）；``text`` 为字符
  串；其余字段任意，随保留记录**原样**交付（kept.jsonl 逐字节回写
  原始行，规范化只用于比较，不改写正文与附带字段）。
* 无法解析的行、字段类型错误、重复 id：逐条指出行号并**拒绝本次
  作业**（退出码 2），不创建输出目录，不把跳过当作清理成功。
* 阈值参数形如 ``p/q``（p、q 正整数且 p<=q）；相似度判断全部使用
  整数交叉乘法（i*q >= p*u），不经过浮点，边界不靠碰巧。

比较口径
--------
1. ``text`` 先做 Unicode NFKC，再 casefold（仅用于比较）；
2. token = Unicode 类别以 L 或 N 开头的连续字符段，其余字符一律作
   分隔，连续分隔不产生空 token；
3. 文档表示为连续三个 token 组成的片段（三词 shingle）的**集合**：
   同一片段重复出现只算一个；比较的是片段**原值**（JSON 编码文本
   落库、等值连接），不使用任何摘要/指纹，不存在碰撞误判；
4. 少于三个 token 的记录标 ``too_short``：保留原文但不加入去重索
   引，也不与任何记录互相比较（两篇空文不会因此被合并）。

清理顺序与判定
--------------
按 ``quality`` 从大到小、同分按 ``id`` 的 Unicode 字典序逐篇考察
（码点序；SQLite BINARY 对 UTF-8 的字节比较与码点序一致）。输入
文件的行序不影响任何决定。对每篇正常文档 A：若没有任何**已保留**
正常文档 B 使 Jaccard(A,B) = |A∩B|/|A∪B| >= p/q（恰好等于也算重
复），则保留 A；否则记 redundant，见证为已保留文档中 Jaccard 最
高者，同分时按清理顺序取最早。冗余文档从此不进入保留索引、不参
与后续判定——不是先求连通分量再每组留一篇。

候选检索（前缀过滤 + 长度过滤：完备过滤，不漏真候选）
------------------------------------------------------
先统计每个片段的文档频次 df（出现在几篇文档中；too_short 没有片
段，自然不参与），按 (df 升序, 片段三元组字典序) 建立整个作业共
用的严格序，片段的位次记 rank。文档 X 的**前缀** P(X) = 按 rank
排序后的前 |X| - ceil(p*|X|/q) + 1 个片段。已保留文档的前缀片段
进磁盘倒排 prefix_index；新文档 A 只拿 P(A) 的片段去倒排里取候选
（磁盘临时表去重），再施加长度过滤
``ceil(p*|A|/q) <= |B| <= floor(q*|A|/p)``，最后对幸存者计算精确
交并集。

* 前缀为什么不漏（达到门槛必共享前缀片段）：若
  Jaccard(A,B) >= p/q，则 |A∩B| >= p*|A∪B|/q；又 |A∪B| >= |A|、
  |B|，故 |A∩B| >= ceil(p*|A|/q) =: tA 且 |A∩B| >= ceil(p*|B|/q)
  =: tB。设 x 为 P(A) 中 rank 最大的元素（A 中 rank 比 x 大的元素
  恰有 tA-1 个），y 为 P(B) 中 rank 最大的元素，不妨 x 在共享序上
  不大于 y。A∩B 中 rank 大于 x 的元素至多 tA-1 个，故至少有
  |A∩B| - (tA-1) >= 1 个元素 z 的 rank 不超过 x。这个 z ∈ A 且
  rank(z) <= rank(x)，所以 z ∈ P(A)；z ∈ B 且
  rank(z) <= rank(x) <= rank(y)，所以 z ∈ P(B)。于是
  P(A) ∩ P(B) 非空：A 的前缀探针必命中 B 已索引的前缀片段。
* 长度界为什么必要（取整）：Jaccard <= min(|A|,|B|)/max(|A|,|B|)，
  要它 >= p/q，当 |B| <= |A| 时需 q*|B| >= p*|A| 即
  |B| >= ceil(p*|A|/q)；当 |B| >= |A| 时需 p*|B| <= q*|A| 即
  |B| <= floor(q*|A|/p)。两式合并即上述闭区间。

不用抽样、MinHash 估计或“候选截前几个”代替完备检索；高频模板造
成大量真实候选时只会变慢，不会漏结果。冗余项与 too_short 项不进
入保留索引。

输出（--out-dir，须不存在或为空）
---------------------------------
* ``kept.jsonl``：保留记录（正常保留 + too_short）的**原始行**，
  按清理顺序写出；id/quality/text 与附带字段逐字节保留；
* ``dispositions.jsonl``：逐篇处置表，按清理顺序。kept 行给出候
  选对数、精确比较对数与候选中相似度最高者（未达阈值；未被列为
  候选的在先保留项经前缀+长度过滤证明不可能达到阈值）；redundant
  行给出见证 id、双方片段数、交集数、并集数与精确相似度 i/u；
  too_short 行给出 token 数；
* ``summary.json``：阈值、文档数、保留数、短文本数、冗余数、前
  缀命中去重后的候选对数、真正执行精确交集的对数、全体配对规模
  等计数。

临时文件
--------
工作库（语料原文、片段、频次/位次表、前缀倒排、候选表）默认建在
系统临时目录，作业结束（无论成败）即删除；``--keep-work`` 可保留
用于排查。工作库是纯粹的派生状态，删除不影响任何输出。
"""
from __future__ import annotations

import json
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import unicodedata
from functools import lru_cache

SUMMARY_SCHEMA_VERSION = "1.0"

DISPOSITION_KEPT = "kept"
DISPOSITION_TOO_SHORT = "too_short"
DISPOSITION_REDUNDANT = "redundant"

INT64_MIN = -(2 ** 63)
INT64_MAX = 2 ** 63 - 1

# 校验失败时最多逐条列出的错误数（错误总数另计，避免坏文件撑爆内存）
MAX_VALIDATION_ERRORS = 20

_THRESHOLD_RE = re.compile(r"([0-9]+)/([0-9]+)\Z")

# SQLite 单条语句的变量上限留足余量后的分批大小
_BATCH = 500


class DedupInputError(ValueError):
    """语料清理输入/配置错误：阈值、语料校验、输出目录（开始前抛出）。"""


# --------------------------------------------------------------------------- #
# 阈值与比较口径
# --------------------------------------------------------------------------- #
def parse_threshold(text: str) -> tuple[int, int]:
    """解析 ``p/q`` 阈值：p、q 正整数且 p<=q；返回 (p, q)。"""
    m = _THRESHOLD_RE.fullmatch(text) if isinstance(text, str) else None
    if not m:
        raise DedupInputError(
            f"阈值必须形如 p/q（p、q 为正整数且 p<=q），得到: {text!r}")
    p, q = int(m.group(1)), int(m.group(2))
    if p < 1 or q < 1:
        raise DedupInputError(f"阈值分子分母必须为正整数: {text!r}")
    if p > q:
        raise DedupInputError(f"阈值要求 p<=q: {text!r}")
    return p, q


def normalize_text(text: str) -> str:
    """比较用规范化：Unicode NFKC 后 casefold（不改写原文）。"""
    return unicodedata.normalize("NFKC", text).casefold()


@lru_cache(maxsize=None)
def _category(ch: str) -> str:
    return unicodedata.category(ch)


def tokenize(text: str) -> list[str]:
    """L*/N* 类别的连续字符段为 token，其余字符一律作分隔。

    连续分隔不产生空 token。先 NFKC 再 casefold 后逐字符判定类别。
    """
    tokens = []
    current = []
    for ch in normalize_text(text):
        if _category(ch)[0] in ("L", "N"):
            current.append(ch)
        elif current:
            tokens.append("".join(current))
            current = []
    if current:
        tokens.append("".join(current))
    return tokens


def shingle_triples(tokens: list[str]) -> set[tuple[str, str, str]]:
    """连续三 token 片段的集合（三元组形式，重复出现只算一个）。"""
    if len(tokens) < 3:
        return set()
    return {tuple(tokens[i:i + 3]) for i in range(len(tokens) - 2)}


def encode_shingle(triple: tuple[str, str, str]) -> str:
    """片段三元组 -> 规范文本（紧凑 JSON 数组编码）。

    token 只含 L*/N* 字符（码点均 >= 0x30），大于 JSON 结构字符
    ``"`` ``,`` ``[`` ``]``，且 token 内不会出现需转义的字符，因此
    编码文本的字典序与三元组逐元素字典序一致（UTF-8 字节序 =
    Unicode 码点序），SQLite BINARY 排序即可实现“片段三元组字典
    序”。等值比较基于编码原文，不存在摘要碰撞。
    """
    return json.dumps(list(triple), ensure_ascii=False, separators=(",", ":"))


def shingles_of(tokens: list[str]) -> set[str]:
    """文档的片段集合（编码文本形式）。"""
    return {encode_shingle(t) for t in shingle_triples(tokens)}


def prefix_length(n: int, p: int, q: int) -> int:
    """前缀长度：n - ceil(p*n/q) + 1（n>=1、p<=q 时恒 >= 1）。"""
    return n - (p * n + q - 1) // q + 1


def length_bounds(n: int, p: int, q: int) -> tuple[int, int]:
    """与大小为 n 的集合可能达到 p/q 的候选大小范围 [lo, hi]。

    lo = ceil(p*n/q)，hi = floor(q*n/p)；推导见模块 docstring。
    """
    return (p * n + q - 1) // q, (q * n) // p


# --------------------------------------------------------------------------- #
# 输入校验
# --------------------------------------------------------------------------- #
class _ErrorLog:
    """逐条收集校验错误（封顶 MAX_VALIDATION_ERRORS 条），最后统一拒绝。"""

    def __init__(self) -> None:
        self.count = 0
        self.messages: list[str] = []

    def add(self, line_no: int, reason: str) -> None:
        self.count += 1
        if len(self.messages) < MAX_VALIDATION_ERRORS:
            self.messages.append(f"第 {line_no} 行: {reason}")

    def raise_if_any(self) -> None:
        if not self.count:
            return
        more = ""
        if self.count > len(self.messages):
            more = f"（其余 {self.count - len(self.messages)} 处从略）"
        raise DedupInputError(
            f"语料校验失败，共 {self.count} 处错误{more}:\n"
            + "\n".join(self.messages))


def _validate_record(obj) -> str | None:
    """校验一条已解析的 JSON 值；返回错误描述或 None。"""
    if not isinstance(obj, dict):
        return "不是 JSON 对象"
    if "id" not in obj:
        return "缺少字段 id"
    if not isinstance(obj["id"], str) or obj["id"] == "":
        return "id 必须是非空字符串"
    if "quality" not in obj:
        return "缺少字段 quality"
    quality = obj["quality"]
    if isinstance(quality, bool) or not isinstance(quality, int):
        return "quality 必须是 64 位有符号整数（不接受 bool/浮点）"
    if not (INT64_MIN <= quality <= INT64_MAX):
        return "quality 超出 64 位有符号整数范围"
    if "text" not in obj:
        return "缺少字段 text"
    if not isinstance(obj["text"], str):
        return "text 必须是字符串"
    return None


# --------------------------------------------------------------------------- #
# 工作库结构
# --------------------------------------------------------------------------- #
def _create_schema(db: sqlite3.Connection) -> None:
    db.executescript(
        """
        CREATE TABLE docs(
            doc_pk         INTEGER PRIMARY KEY,
            line_no        INTEGER NOT NULL,
            id             TEXT NOT NULL UNIQUE,
            quality        INTEGER NOT NULL,
            line           TEXT NOT NULL,
            token_count    INTEGER NOT NULL,
            shingle_count  INTEGER NOT NULL,
            status         TEXT,
            seq            INTEGER
        );
        -- 每篇文档的片段集合（集合语义：装载前已在内存去重）
        CREATE TABLE doc_shingles(
            doc_pk  INTEGER NOT NULL,
            shingle TEXT NOT NULL,
            PRIMARY KEY (doc_pk, shingle)
        );
        -- 已保留文档的前缀倒排：片段 -> 保留文档
        CREATE TABLE prefix_index(
            shingle TEXT NOT NULL,
            doc_pk  INTEGER NOT NULL,
            PRIMARY KEY (shingle, doc_pk)
        );
        -- 当前文档的前缀探针
        CREATE TEMP TABLE cur_prefix(shingle TEXT PRIMARY KEY);
        -- 前缀命中去重后的候选（磁盘去重，不累积成 Python 容器）
        CREATE TEMP TABLE cand(doc_pk INTEGER PRIMARY KEY);
        -- 长度过滤后的候选（逐对计算精确交集）
        CREATE TEMP TABLE cand_len(doc_pk INTEGER PRIMARY KEY);
        """
    )


def _load_corpus(db: sqlite3.Connection, corpus_path: str,
                 errors: _ErrorLog) -> int:
    """第一遍：流式读取语料，校验并落库（原文 + 片段集合）。

    任何坏行只记录位置与原因，最后由 :meth:`_ErrorLog.raise_if_any`
    统一拒绝本次作业。返回合法记录数。
    """
    n_docs = 0
    with open(corpus_path, "rb") as fh:
        for line_no, raw in enumerate(fh, 1):
            try:
                line = raw.decode("utf-8")
            except UnicodeDecodeError as exc:
                errors.add(line_no, f"不是合法 UTF-8: {exc}")
                continue
            if line.endswith("\n"):
                line = line[:-1]
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                errors.add(line_no, f"JSON 解析失败: {exc}")
                continue
            problem = _validate_record(obj)
            if problem is not None:
                errors.add(line_no, problem)
                continue
            tokens = tokenize(obj["text"])
            shingle_set = shingles_of(tokens)
            try:
                cur = db.execute(
                    "INSERT INTO docs(line_no, id, quality, line,"
                    "                  token_count, shingle_count)"
                    " VALUES (?,?,?,?,?,?)",
                    (line_no, obj["id"], obj["quality"], line,
                     len(tokens), len(shingle_set)))
            except sqlite3.IntegrityError:
                first = db.execute(
                    "SELECT line_no FROM docs WHERE id = ?",
                    (obj["id"],)).fetchone()
                errors.add(line_no, f"id 重复: {obj['id']!r}"
                                    f"（首次出现于第 {first[0]} 行）")
                continue
            doc_pk = cur.lastrowid
            n_docs += 1
            if shingle_set:
                db.executemany(
                    "INSERT INTO doc_shingles(doc_pk, shingle)"
                    " VALUES (?,?)",
                    [(doc_pk, s) for s in shingle_set])
    return n_docs


def _build_rank(db: sqlite3.Connection) -> None:
    """全局片段序：(df 升序, 片段三元组字典序) -> 行号 rank。"""
    db.executescript(
        """
        CREATE INDEX idx_ds_shingle ON doc_shingles(shingle, doc_pk);
        CREATE TABLE shingle_rank AS
            SELECT shingle, df,
                   ROW_NUMBER() OVER (ORDER BY df, shingle) AS rank
            FROM (SELECT shingle, COUNT(*) AS df
                  FROM doc_shingles GROUP BY shingle);
        CREATE UNIQUE INDEX idx_rank_shingle ON shingle_rank(shingle);
        """
    )


def _build_clean_order(db: sqlite3.Connection) -> None:
    """清理顺序物化：quality 降序、id 字典序（rowid 即清理序号 seq）。

    排序在 SQLite 内完成（必要时落临时文件），不把全文件读成
    Python 列表；物化后逐行扫描不受后续 UPDATE docs 影响。
    """
    db.executescript(
        """
        CREATE TEMP TABLE clean_order AS
            SELECT doc_pk, id, token_count, shingle_count
            FROM docs ORDER BY quality DESC, id ASC;
        """
    )


# --------------------------------------------------------------------------- #
# 第二遍：优先级遍历 + 候选检索 + 判定
# --------------------------------------------------------------------------- #
def _write_disp(fh, disp: dict) -> None:
    fh.write(json.dumps(disp, ensure_ascii=False) + "\n")


def _clean(db: sqlite3.Connection, p: int, q: int, disp_fh) -> dict:
    """按清理顺序逐篇判定，流式写出处置记录；返回检索/处置计数。"""
    stats = {"kept": 0, "too_short": 0, "redundant": 0,
             "candidate_pairs": 0, "exact_comparisons": 0}
    rows = db.execute(
        "SELECT rowid, doc_pk, id, token_count, shingle_count"
        " FROM clean_order ORDER BY rowid")
    for seq, doc_pk, doc_id, token_count, n in rows:
        if n == 0:
            # 少于三个 token：留档但不入索引、不参与比较
            db.execute("UPDATE docs SET status = ?, seq = ?"
                       " WHERE doc_pk = ?",
                       (DISPOSITION_TOO_SHORT, seq, doc_pk))
            stats["too_short"] += 1
            _write_disp(disp_fh, {
                "seq": seq, "id": doc_id,
                "disposition": DISPOSITION_TOO_SHORT,
                "token_count": token_count, "shingle_count": 0,
                "note": "少于三个 token：保留原文，不加入去重索引，"
                        "不参与任何比较"})
            continue

        # --- 前缀探针（全局共享序的前 n-ceil(pn/q)+1 项） ---------------
        # CROSS JOIN 强制以本文档的片段为外层（temp 表无统计信息，
        # 不能让规划器把 shingle_rank 全表当外层扫）
        plen = prefix_length(n, p, q)
        prefix = [r[0] for r in db.execute(
            "SELECT s.shingle FROM doc_shingles s"
            " CROSS JOIN shingle_rank r ON r.shingle = s.shingle"
            " WHERE s.doc_pk = ? ORDER BY r.rank LIMIT ?",
            (doc_pk, plen))]
        db.execute("DELETE FROM cur_prefix")
        db.executemany("INSERT INTO cur_prefix(shingle) VALUES (?)",
                       [(s,) for s in prefix])

        # --- 候选：共享前缀片段对应的已保留文档（磁盘去重） -------------
        # CROSS JOIN 强制以前缀探针（十余项）为外层、按 (shingle) 主键
        # 探测倒排；否则规划器会把整个 prefix_index 当外层逐篇全扫
        db.execute("DELETE FROM cand")
        db.execute(
            "INSERT OR IGNORE INTO cand(doc_pk)"
            " SELECT pi.doc_pk FROM cur_prefix cp"
            " CROSS JOIN prefix_index pi ON pi.shingle = cp.shingle")
        n_cand = db.execute("SELECT COUNT(*) FROM cand").fetchone()[0]
        stats["candidate_pairs"] += n_cand

        # --- 长度过滤：ceil(pn/q) <= |B| <= floor(qn/p) -----------------
        # CROSS JOIN 强制以候选小表为外层、按 rowid 主键探测 docs
        lo, hi = length_bounds(n, p, q)
        db.execute("DELETE FROM cand_len")
        db.execute(
            "INSERT OR IGNORE INTO cand_len(doc_pk)"
            " SELECT c.doc_pk FROM cand c"
            " CROSS JOIN docs d ON d.doc_pk = c.doc_pk"
            " WHERE d.shingle_count >= ? AND d.shingle_count <= ?",
            (lo, hi))
        n_exact = db.execute("SELECT COUNT(*) FROM cand_len").fetchone()[0]
        stats["exact_comparisons"] += n_exact

        # --- 精确交集（原值等值连接，结果分批流回，不累积成容器） ------
        # CROSS JOIN 强制以本文档片段为外层、按 (shingle, doc_pk) 索引
        # 探测候选侧
        best = None  # (inter, union, seq, doc_pk, id, shingle_count)
        if n_exact:
            for b_pk, b_id, b_m, b_seq, inter in db.execute(
                    "SELECT b.doc_pk, d.id, d.shingle_count, d.seq,"
                    "       COUNT(*)"
                    " FROM doc_shingles a"
                    " CROSS JOIN doc_shingles b ON b.shingle = a.shingle"
                    " JOIN docs d ON d.doc_pk = b.doc_pk"
                    " WHERE a.doc_pk = ?"
                    "   AND b.doc_pk IN (SELECT doc_pk FROM cand_len)"
                    " GROUP BY b.doc_pk", (doc_pk,)):
                union = n + b_m - inter
                if best is None or inter * best[1] > best[0] * union or \
                        (inter * best[1] == best[0] * union
                         and b_seq < best[2]):
                    best = (inter, union, b_seq, b_pk, b_id, b_m)

        if best is not None and best[0] * q >= p * best[1]:
            # 达到阈值（整数交叉乘法，恰好等于也算重复）-> redundant
            inter, union, _, _b_pk, b_id, b_m = best
            db.execute("UPDATE docs SET status = ?, seq = ?"
                       " WHERE doc_pk = ?",
                       (DISPOSITION_REDUNDANT, seq, doc_pk))
            stats["redundant"] += 1
            _write_disp(disp_fh, {
                "seq": seq, "id": doc_id,
                "disposition": DISPOSITION_REDUNDANT,
                "shingle_count": n,
                "witness_id": b_id,
                "witness_shingle_count": b_m,
                "intersection": inter,
                "union": union,
                "similarity": f"{inter}/{union}",
                "candidates": n_cand,
                "examined": n_exact})
        else:
            # 保留：前缀片段入倒排，参与后续判定
            db.execute("UPDATE docs SET status = ?, seq = ?"
                       " WHERE doc_pk = ?",
                       (DISPOSITION_KEPT, seq, doc_pk))
            db.execute(
                "INSERT INTO prefix_index(shingle, doc_pk)"
                " SELECT shingle, ? FROM cur_prefix", (doc_pk,))
            stats["kept"] += 1
            _write_disp(disp_fh, {
                "seq": seq, "id": doc_id,
                "disposition": DISPOSITION_KEPT,
                "shingle_count": n,
                "candidates": n_cand,
                "examined": n_exact,
                "best_below_threshold": (
                    {"id": best[4],
                     "intersection": best[0],
                     "union": best[1],
                     "similarity": f"{best[0]}/{best[1]}"}
                    if best is not None else None),
                "note": "在先保留项均未达到阈值（未列为候选者经前缀"
                        "+长度过滤证明不可能达到）"})
    return stats


def _write_kept(db: sqlite3.Connection, path: str) -> int:
    """保留记录（正常保留 + too_short）的原始行，按清理顺序写出。"""
    n = 0
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for (line,) in db.execute(
                "SELECT line FROM docs"
                " WHERE status IN (?, ?) ORDER BY seq",
                (DISPOSITION_KEPT, DISPOSITION_TOO_SHORT)):
            fh.write(line + "\n")
            n += 1
    return n


# --------------------------------------------------------------------------- #
# 作业入口
# --------------------------------------------------------------------------- #
def run_dedup(corpus_path: str, threshold_text: str, out_dir: str,
              keep_work: bool = False) -> dict:
    """正文语料清理作业：校验 -> 全局频次/位次 -> 优先级清理 -> 写出。

    返回 summary 字典；输入/配置错误抛 :class:`DedupInputError`。
    """
    p, q = parse_threshold(threshold_text)
    if not os.path.isfile(corpus_path):
        raise DedupInputError(f"语料文件不存在: {corpus_path}")
    if os.path.exists(out_dir):
        if not os.path.isdir(out_dir):
            raise DedupInputError(f"输出路径已存在且不是目录: {out_dir}")
        if os.listdir(out_dir):
            raise DedupInputError(f"输出目录已存在且非空: {out_dir}")

    work = tempfile.mkdtemp(prefix="vcfpipe_dedup_")
    try:
        db = sqlite3.connect(os.path.join(work, "work.db"),
                             isolation_level=None)
        try:
            db.execute("PRAGMA journal_mode = OFF")
            db.execute("PRAGMA synchronous = OFF")
            db.execute("PRAGMA cache_size = -65536")  # 64MB，上限有界
            _create_schema(db)

            errors = _ErrorLog()
            db.execute("BEGIN")
            n_docs = _load_corpus(db, corpus_path, errors)
            db.execute("COMMIT")
            # 有任何坏行：指出位置并拒绝这次作业（不创建输出目录）
            errors.raise_if_any()

            # 注意：executescript 会先隐式提交挂起事务，这两步不包显式
            # 事务（工作库是派生状态，journal 已关，无需原子性）
            _build_rank(db)
            _build_clean_order(db)

            os.makedirs(out_dir, exist_ok=True)
            disp_path = os.path.join(out_dir, "dispositions.jsonl")
            with open(disp_path, "w", encoding="utf-8",
                      newline="") as disp_fh:
                db.execute("BEGIN")
                stats = _clean(db, p, q, disp_fh)
                db.execute("COMMIT")
            kept_lines = _write_kept(db, os.path.join(out_dir,
                                                      "kept.jsonl"))
            distinct = db.execute(
                "SELECT COUNT(*) FROM shingle_rank").fetchone()[0]

            summary = {
                "schema": SUMMARY_SCHEMA_VERSION,
                "threshold": f"{p}/{q}",
                "threshold_p": p,
                "threshold_q": q,
                "documents": n_docs,
                "kept": stats["kept"],
                "too_short": stats["too_short"],
                "redundant": stats["redundant"],
                "candidate_pairs": stats["candidate_pairs"],
                "exact_comparisons": stats["exact_comparisons"],
                "pairwise_total": n_docs * (n_docs - 1) // 2,
                "distinct_shingles": distinct,
                "kept_jsonl_lines": kept_lines,
                "outputs": {
                    "kept": "kept.jsonl",
                    "dispositions": "dispositions.jsonl",
                    "summary": "summary.json",
                },
            }
            with open(os.path.join(out_dir, "summary.json"), "w",
                      encoding="utf-8", newline="") as fh:
                json.dump(summary, fh, indent=2, ensure_ascii=False,
                          sort_keys=True)
                fh.write("\n")
            return summary
        finally:
            db.close()
    finally:
        if keep_work:
            print(f"临时工作库保留于 {work}", file=sys.stderr)
        else:
            shutil.rmtree(work, ignore_errors=True)
