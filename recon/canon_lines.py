"""Canonicalize DNA lines: each line becomes min(line, reverse_complement(line)). stdin -> stdout.

Used to ask how much of km-peeler's (strand-specific) extended-syncmer difference between two assemblies is only
orientation: canonicalize both extended-syncmer sets, sort -u, and count the difference again.
"""
import sys

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def main() -> int:
    out = sys.stdout.buffer
    for line in sys.stdin.buffer:
        s = line.rstrip(b"\n")
        r = s.translate(COMP)[::-1]
        out.write((s if s <= r else r) + b"\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
