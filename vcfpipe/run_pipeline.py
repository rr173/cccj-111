"""项目本地运行入口：保证 -I 隔离模式下也能找到 src 包。"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from vcfpipe.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
