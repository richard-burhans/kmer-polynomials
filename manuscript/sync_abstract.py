#!/usr/bin/env python3
"""Keep ONE source for the abstract: copy it from paper.md into ../README.md.

    python3 sync_abstract.py           write the abstract into ../README.md
    python3 sync_abstract.py --check   exit 1 if ../README.md's abstract differs from paper.md's

Both files carry the abstract between the same two marker comments. Only text between the markers is touched.

Inline math: paper.md uses pandoc's $...$. GitHub renders $...$ unreliably (it fails next to punctuation such as an
em-dash), so the README copy rewrites each $...$ as GitHub's robust $`...`$ form. --check compares after that rewrite.
"""
import pathlib
import re
import sys

#: a pandoc inline-math span: $ not followed by space, content without $, $ not preceded by space
MATH = re.compile(r"(?<![\\$`])\$(?=[^\s`])([^$\n]+?)(?<=[^\s`])\$(?![`$])")  # skips spans already in $`...`$ form


def github_math(text: str) -> str:
    """Rewrite pandoc $...$ inline math as GitHub's $`...`$ form."""
    return MATH.sub(lambda m: "$`" + m.group(1) + "`$", text)

BEGIN, END = "<!-- ABSTRACT:BEGIN -->", "<!-- ABSTRACT:END -->"
HERE = pathlib.Path(__file__).resolve().parent
PAPER, README = HERE / "paper.md", HERE.parent / "README.md"


def between(text: str, where: pathlib.Path) -> tuple[int, int]:
    """Start and end offsets of the text between the markers. Exits if the markers are missing or out of order."""
    a, b = text.find(BEGIN), text.find(END)
    if a < 0 or b < 0 or b < a or text.count(BEGIN) != 1 or text.count(END) != 1:
        sys.exit(f"{where}: needs exactly one {BEGIN} followed by one {END}")
    return a + len(BEGIN), b


def main() -> int:
    paper, readme = PAPER.read_text(encoding="utf-8"), README.read_text(encoding="utf-8")
    pa, pb = between(paper, PAPER)
    ra, rb = between(readme, README)
    abstract = github_math(paper[pa:pb])
    if "--check" in sys.argv[1:]:
        if readme[ra:rb] != abstract:
            print(f"{README.name}: abstract differs from {PAPER.name}; run `make readme`")
            return 1
        print("abstract in sync")
        return 0
    README.write_text(readme[:ra] + abstract + readme[rb:], encoding="utf-8")
    print(f"wrote the abstract into {README}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
