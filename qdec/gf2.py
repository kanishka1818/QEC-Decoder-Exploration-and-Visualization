"""
gf2.py - linear algebra over the two-element field GF(2) = {0, 1}.

What this file is FOR
---------------------
Every question we ask about a stabilizer code ("do the checks commute?",
"how many logical qubits?", "what are the logical operators?") is a question
about 0/1 matrices where 1 + 1 = 0.  numpy does ordinary arithmetic, so we
need our own small toolkit: row reduction, rank, null space, and "is this
vector inside that row space?".  Later the OSD decoder will reuse `solve`.

Everything here is plain numpy, no QEC libraries.

Reference: any coding-theory text, e.g. MacWilliams & Sloane,
"The Theory of Error-Correcting Codes" (1977), ch. 1; the same routines
appear in Roffe's `ldpc` package (mod2 module), which we only use to
cross-check in the tests.
"""

from __future__ import annotations

import numpy as np


def as_gf2(M) -> np.ndarray:
    """Return a copy of M as a uint8 array with every entry reduced mod 2."""
    return (np.array(M, dtype=np.int64) % 2).astype(np.uint8)


def row_reduce(M) -> tuple[np.ndarray, list[int]]:
    """
    Gaussian elimination mod 2 (reduced row-echelon form).

    Returns
    -------
    R : the reduced matrix (same shape as M)
    pivots : list of pivot column indices, one per non-zero row of R
    """
    R = as_gf2(M)
    n_rows, n_cols = R.shape
    pivots: list[int] = []
    row = 0                                   # next row that needs a pivot
    for col in range(n_cols):
        if row >= n_rows:
            break
        # find a row at or below `row` with a 1 in this column
        candidates = np.nonzero(R[row:, col])[0]
        if candidates.size == 0:
            continue                          # no pivot in this column
        pivot_row = row + candidates[0]
        if pivot_row != row:                  # swap it up into place
            R[[row, pivot_row]] = R[[pivot_row, row]]
        # clear every OTHER 1 in this column (adding rows = XOR mod 2)
        ones = np.nonzero(R[:, col])[0]
        for r in ones:
            if r != row:
                R[r] ^= R[row]
        pivots.append(col)
        row += 1
    return R, pivots


def rank(M) -> int:
    """Rank of M over GF(2) = number of pivots after row reduction."""
    return len(row_reduce(M)[1])


def nullspace(M) -> np.ndarray:
    """
    Basis of {v : M v = 0 (mod 2)}, returned as the ROWS of a matrix.

    Method: in reduced row-echelon form every non-pivot ("free") column gives
    one basis vector: set that free variable to 1, the other free variables
    to 0, and read the pivot variables off the reduced rows.
    """
    R, pivots = row_reduce(M)
    n_cols = R.shape[1]
    free = [c for c in range(n_cols) if c not in set(pivots)]
    basis = np.zeros((len(free), n_cols), dtype=np.uint8)
    for k, f in enumerate(free):
        basis[k, f] = 1
        for r, p in enumerate(pivots):        # pivot var = value in column f
            basis[k, p] = R[r, f]
    return basis


def in_row_space(M, v) -> bool:
    """True if vector v is a mod-2 sum of rows of M (adding v keeps rank)."""
    return rank(np.vstack([as_gf2(M), as_gf2(v)[None, :]])) == rank(M)


def solve(M, b) -> np.ndarray | None:
    """
    One solution x of M x = b (mod 2), or None if there is none.
    Free variables are set to 0.  (Used later by the OSD decoder.)
    """
    M = as_gf2(M)
    b = as_gf2(b).reshape(-1, 1)
    R, pivots = row_reduce(np.hstack([M, b]))
    n = M.shape[1]
    if n in pivots:                           # pivot in the b column => 0 = 1
        return None
    x = np.zeros(n, dtype=np.uint8)
    for r, p in enumerate(pivots):
        x[p] = R[r, n]
    return x
