"""Shared helpers for the KP-4 / KP-13 spectral experiments (IDEAS.md KP-4, KP-13).

Encoding (IDEAS.md section 0): A = 1, C = i, G = -i, T = -1, so complement is negation. N (and anything else) is 0.
`set_encoding("eisenstein")` switches to A = 1, C = w, G = -w, T = -1 with w = exp(2 pi i / 3) (IDEAS.md section 7,
results/kp7_eisenstein_2026-11-10.md): complement is still negation, and squared distances become 1 for transitions
(A-G, C-T), 3 for the A-C / G-T transversions and 4 for complements.
Sequences are kept as uint8 codes (A0 C1 G2 T3, N4) and mapped to complex only per segment, to bound memory.
"""
from __future__ import annotations

import pathlib

import numpy as np

CODE = np.full(256, 4, dtype=np.uint8)
for _ch, _v in zip("ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    CODE[ord(_ch)] = _v
#: complex value per code; index 4 (N) is 0
VAL = np.array([1, 1j, -1j, -1, 0], dtype=np.complex64)
_W = complex(-0.5, 3 ** 0.5 / 2)
#: encodings by name: value per code A, C, G, T, N
ENCODINGS = {"quarter": np.array([1, 1j, -1j, -1, 0], dtype=np.complex64),
             "eisenstein": np.array([1, _W, -_W, -1, 0], dtype=np.complex64)}
ENCODING = "quarter"


def set_encoding(name: str) -> None:
    """Select the base -> complex map used by to_complex (and read by callers through `values()`)."""
    global VAL, ENCODING
    VAL = ENCODINGS[name]
    ENCODING = name


def values(dtype=np.complex128) -> np.ndarray:
    """The current encoding's values for codes A, C, G, T (no N)."""
    return VAL[:4].astype(dtype)


#: complement code (A<->T, C<->G, N stays N)
COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)


def read_fasta_one(path: pathlib.Path) -> tuple[str, np.ndarray, np.ndarray]:
    """The first record: (name, codes uint8, soft-masked bool)."""
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    break
                name = line[1:].split()[0].decode()
                continue
            chunks.append(line.strip())
    raw = np.frombuffer(b"".join(chunks), dtype=np.uint8)
    lower = (raw >= ord("a")) & (raw <= ord("z"))
    return name, CODE[raw], lower


def revcomp(codes: np.ndarray, lower: np.ndarray | None = None):
    rc = COMP[codes[::-1]]
    return (rc, lower[::-1].copy()) if lower is not None else rc


def to_complex(codes: np.ndarray, lower: np.ndarray | None = None, mask_soft: bool = False) -> np.ndarray:
    v = VAL[codes]
    if mask_soft and lower is not None:
        v = np.where(lower, 0, v).astype(np.complex64)
    return v


def read_chains(path: pathlib.Path, tname: str, qname: str):
    """Chains between tname and qname: list of dicts with strand, score and gapless blocks (t, q, size).

    For a '-' chain, q is in reverse-complement coordinates of the query (UCSC chain format), so the block aligns
    target[t:t+size] with revcomp(query)[q:q+size].
    """
    out, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("chain"):
                f = line.split()
                if cur is not None:
                    out.append(cur)
                cur = None
                if f[2] == tname and f[7] == qname and f[4] == "+":
                    cur = {"score": int(f[1]), "strand": f[9], "t": int(f[5]), "q": int(f[10]), "blocks": []}
                continue
            if cur is None:
                continue
            f = line.split()
            if not f:
                out.append(cur)
                cur = None
                continue
            size = int(f[0])
            cur["blocks"].append((cur["t"], cur["q"], size))
            if len(f) == 3:
                cur["t"] += size + int(f[1])
                cur["q"] += size + int(f[2])
    if cur is not None:
        out.append(cur)
    return out
