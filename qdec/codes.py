"""
codes.py - build bivariate-bicycle (BB) codes from two polynomials.

What this file is FOR
---------------------
It turns a short recipe - a torus size (l, m) and two polynomials A, B in
the variables x and y - into everything a decoder needs:

    * the check matrices H_X and H_Z,
    * the number of physical qubits n and logical qubits k,
    * a set of logical operators (X and Z, paired up),
    * a picture-friendly position for every qubit and check on the torus.

One recipe covers the toric code, the [[72,12,6]] code and IBM's gross code
[[144,12,12]], so all three can be studied with the same tools.

The recipe (Bravyi et al., Nature 627, 778 (2024), arXiv:2308.07915)
-----------------------------------------------------------------------
    x = (cyclic shift on l)  (x)  (identity on m)        an lm x lm matrix
    y = (identity on l)      (x)  (cyclic shift on m)    an lm x lm matrix
    A, B = sums of monomials x^a y^b                      (mod 2)
    H_X = [ A | B ]          H_Z = [ B^T | A^T ]

Because x and y commute, A B = B A, so H_X H_Z^T = A B + B A = 0 (mod 2):
every X check commutes with every Z check, for ANY choice of A and B.

Qubit layout
------------
The torus has l x m "cells", cell (i, j) holds four things:
    one L data qubit, one R data qubit, one X check, one Z check.
Data qubit index = i*m + j for L, and l*m + i*m + j for R.
On the picture grid (period 2l rows x 2m columns) we place:
    X check (2i,   2j)      R qubit (2i,   2j+1)
    L qubit (2i+1, 2j)      Z check (2i+1, 2j+1)
With our shift convention (x moves i -> i-1) the toric code A = 1 + x,
B = 1 + y comes out with every check touching its 4 grid neighbours, i.e.
Kitaev's toric code drawn the usual way.  (Bravyi et al. use the opposite
shift direction; that just mirrors the torus and gives the same code.)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .gf2 import as_gf2, in_row_space, nullspace, rank, row_reduce

Monomial = tuple[int, int]          # (a, b) stands for x^a y^b


# --------------------------------------------------------------------------
# polynomial -> matrix
# --------------------------------------------------------------------------
def cyclic_shift(size: int) -> np.ndarray:
    """size x size matrix S with S[a, (a-1) mod size] = 1 (moves index down by one)."""
    S = np.zeros((size, size), dtype=np.uint8)
    for a in range(size):
        S[a, (a - 1) % size] = 1
    return S


def shift_matrices(l: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """The two commuting lm x lm matrices x and y of the BB construction."""
    x = np.kron(cyclic_shift(l), np.eye(m, dtype=np.uint8)).astype(np.uint8)
    y = np.kron(np.eye(l, dtype=np.uint8), cyclic_shift(m)).astype(np.uint8)
    return x, y


def poly_matrix(terms: list[Monomial], l: int, m: int) -> np.ndarray:
    """Matrix of the polynomial sum_k x^{a_k} y^{b_k} (mod 2)."""
    x, y = shift_matrices(l, m)
    M = np.zeros((l * m, l * m), dtype=np.int64)
    for a, b in terms:
        # x^a y^b as a matrix product; matrix_power works on ints
        M += np.linalg.matrix_power(x.astype(np.int64), a) @ \
             np.linalg.matrix_power(y.astype(np.int64), b)
    return as_gf2(M)


def poly_str(terms: list[Monomial]) -> str:
    """Human-readable polynomial, e.g. [(3,0),(0,1)] -> 'x^3 + y'."""
    def mono(a, b):
        parts = []
        if a: parts.append("x" if a == 1 else f"x^{a}")
        if b: parts.append("y" if b == 1 else f"y^{b}")
        return "".join(parts) or "1"
    return " + ".join(mono(a, b) for a, b in terms)


# --------------------------------------------------------------------------
# GF(2) helper used for pairing logical operators
# --------------------------------------------------------------------------
def _gf2_inverse(G: np.ndarray) -> np.ndarray:
    """Inverse of a square invertible 0/1 matrix over GF(2)."""
    k = G.shape[0]
    R, pivots = row_reduce(np.hstack([as_gf2(G), np.eye(k, dtype=np.uint8)]))
    if pivots[:k] != list(range(k)):
        raise ValueError("matrix is not invertible over GF(2)")
    return R[:, k:]


def _independent_mod(span: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    """
    From `candidates` (rows), keep the ones that are NOT already explained by
    `span` (rows) or by candidates kept earlier.  These are representatives of
    the quotient space  span(candidates) / span(span).
    """
    kept = []
    current = as_gf2(span)
    r = rank(current)
    for v in candidates:
        trial = np.vstack([current, v[None, :]])
        r_new = rank(trial)
        if r_new > r:                      # v adds something new
            kept.append(v)
            current, r = trial, r_new
    return np.array(kept, dtype=np.uint8).reshape(-1, span.shape[1])


# --------------------------------------------------------------------------
# the code object
# --------------------------------------------------------------------------
@dataclass
class BBCode:
    """A bivariate-bicycle code built from (l, m, A, B)."""

    name: str
    l: int
    m: int
    A: list[Monomial]
    B: list[Monomial]
    HX: np.ndarray = field(init=False, repr=False)
    HZ: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        Am = poly_matrix(self.A, self.l, self.m)
        Bm = poly_matrix(self.B, self.l, self.m)
        self.HX = np.hstack([Am, Bm])                    # X checks: [A | B]
        self.HZ = np.hstack([Bm.T, Am.T])                # Z checks: [B^T | A^T]
        self._logicals = None                            # computed lazily

    # ---- sizes ----------------------------------------------------------
    @property
    def n(self) -> int:
        """Number of data qubits (L block + R block)."""
        return 2 * self.l * self.m

    @property
    def k(self) -> int:
        """Logical qubits: k = n - rank(H_X) - rank(H_Z) over GF(2)."""
        return self.n - rank(self.HX) - rank(self.HZ)

    def commutes(self) -> bool:
        """True if every X check commutes with every Z check (H_X H_Z^T = 0)."""
        return not np.any((self.HX.astype(np.int64) @ self.HZ.T.astype(np.int64)) % 2)

    def summary(self) -> str:
        return (f"{self.name}: l={self.l}, m={self.m}, A = {poly_str(self.A)}, "
                f"B = {poly_str(self.B)}  ->  n = {self.n}, k = {self.k}, "
                f"check weight = {int(self.HX.sum(1).max())}")

    # ---- logical operators ---------------------------------------------
    def logicals(self) -> tuple[np.ndarray, np.ndarray]:
        """
        (LX, LZ): k X-type and k Z-type logical operators as 0/1 rows, paired
        so that LX[i] and LZ[j] overlap an odd number of times only if i == j.

        X logicals = vectors that pass every Z check (null space of H_Z) but
        are not products of X checks (not in the row space of H_X).  Same for
        Z logicals with the roles swapped.  Then we fix up the pairing with
        the inverse of the overlap matrix.
        """
        if self._logicals is None:
            LX = _independent_mod(self.HX, nullspace(self.HZ))
            LZ = _independent_mod(self.HZ, nullspace(self.HX))
            G = (LX.astype(np.int64) @ LZ.T.astype(np.int64)) % 2  # overlap matrix
            LZ = (_gf2_inverse(G).T.astype(np.int64) @ LZ.astype(np.int64)) % 2
            self._logicals = (LX.astype(np.uint8), LZ.astype(np.uint8))
        return self._logicals

    # ---- positions on the torus picture ---------------------------------
    #
    # Cells are elements g = (i, j) of the group Z_l x Z_m.  The monomial
    # x^a y^b moves a cell by v(x^a y^b) = (-a, -b) (our shift convention):
    #     X check c touches L qubits c + v(A_k) and R qubits c + v(B_k)
    #     Z check c touches L qubits c - v(B_k) and R qubits c - v(A_k)
    #
    # "toric layout" (Bravyi et al. 2024, Sec. "Tanner graph and layout"):
    # pick two terms A_p, A_q and two terms B_p, B_q and use
    #     s = v(A_p) - v(A_q)   as the "one row down" step
    #     t = v(B_p) - v(B_q)   as the "one column right" step.
    # If s and t together reach every cell exactly once (ord(s)*ord(t) = lm)
    # we can draw the code on an ord(s) x ord(t) grid where each check has
    # FOUR nearest-neighbour edges (from A_p, A_q, B_p, B_q) and every other
    # monomial becomes a long-range edge.  The toric code then looks like the
    # textbook picture, and the gross code looks like "toric code + 2 long
    # edges per check".  If no such choice exists we fall back to rows = i,
    # columns = j ("plain" layout).
    def _v(self, mono: Monomial) -> np.ndarray:
        """Cell displacement of a monomial x^a y^b."""
        return np.array([-mono[0] % self.l, -mono[1] % self.m])

    def _order(self, step: np.ndarray) -> int:
        """How many times `step` must be added to come back to (0, 0)."""
        g, k = step.copy(), 1
        while np.any(g % [self.l, self.m]):
            g, k = g + step, k + 1
        return k

    def layout(self) -> dict:
        """
        Work out (and cache) the drawing layout.  Returns a dict with
        'kind' ('toric' or 'plain'), grid size 'rows' x 'cols', the cell
        coordinate of every group element, and the chosen terms.
        """
        if getattr(self, "_layout", None) is not None:
            return self._layout
        mod = np.array([self.l, self.m])
        for p, q in [(p, q) for p in range(len(self.A)) for q in range(len(self.A)) if p != q]:
            s = (self._v(self.A[p]) - self._v(self.A[q])) % mod
            for bp, bq in [(a, b) for a in range(len(self.B)) for b in range(len(self.B)) if a != b]:
                t = (self._v(self.B[bp]) - self._v(self.B[bq])) % mod
                rows, cols = self._order(s), self._order(t)
                if rows * cols != self.l * self.m:
                    continue
                coord = {}
                for a in range(rows):
                    for b in range(cols):
                        coord[tuple((a * s + b * t) % mod)] = (a, b)
                if len(coord) == self.l * self.m:       # every cell reached once
                    self._layout = dict(kind="toric", rows=rows, cols=cols, coord=coord,
                                        Ap=self.A[p], Bp=self.B[bp])
                    return self._layout
        coord = {(i, j): (i, j) for i in range(self.l) for j in range(self.m)}
        self._layout = dict(kind="plain", rows=self.l, cols=self.m, coord=coord,
                            Ap=(0, 0), Bp=(0, 0))
        return self._layout

    def grid_shape(self) -> tuple[int, int]:
        """(rows, cols) of the picture grid = 2 x (layout rows, layout cols)."""
        lay = self.layout()
        return 2 * lay["rows"], 2 * lay["cols"]

    def _cells(self):
        """Group element (i, j) of every cell in index order i*m + j."""
        return [np.array([i, j]) for i in range(self.l) for j in range(self.m)]

    def data_positions(self) -> np.ndarray:
        """(n, 2) array of (row, col) grid positions for data qubits L then R."""
        lay, mod = self.layout(), np.array([self.l, self.m])
        pos = []
        for g in self._cells():                       # L qubit g sits in cell g - v(A_p), below its check
            a, b = lay["coord"][tuple((g - self._v(lay["Ap"])) % mod)]
            pos.append((2 * a + 1, 2 * b))
        for g in self._cells():                       # R qubit g sits in cell g - v(B_p), right of its check
            a, b = lay["coord"][tuple((g - self._v(lay["Bp"])) % mod)]
            pos.append((2 * a, 2 * b + 1))
        return np.array(pos)

    def check_positions(self, kind: str) -> np.ndarray:
        """(l*m, 2) grid positions of the X checks (kind='X') or Z checks ('Z')."""
        lay, mod = self.layout(), np.array([self.l, self.m])
        pos = []
        for g in self._cells():
            if kind == "X":                           # X check c sits in cell c
                a, b = lay["coord"][tuple(g % mod)]
                pos.append((2 * a, 2 * b))
            else:                                     # Z check c sits in cell c - v(B_p) - v(A_p)
                a, b = lay["coord"][tuple((g - self._v(lay["Bp"]) - self._v(lay["Ap"])) % mod)]
                pos.append((2 * a + 1, 2 * b + 1))
        return np.array(pos)

    def H(self, kind: str) -> np.ndarray:
        """H_X for kind='X', H_Z for kind='Z'."""
        return self.HX if kind == "X" else self.HZ


# --------------------------------------------------------------------------
# the three codes we study
# --------------------------------------------------------------------------
def toric_code(L: int) -> BBCode:
    """Kitaev's toric code on an L x L torus: A = 1 + x, B = 1 + y.  [[2L^2, 2, L]]"""
    return BBCode(f"toric {L}x{L}", L, L, [(0, 0), (1, 0)], [(0, 0), (0, 1)])


def bb72() -> BBCode:
    """[[72,12,6]] BB code (Bravyi et al. 2024, Table 3): l=m=6, A=x^3+y+y^2, B=y^3+x+x^2."""
    return BBCode("BB [[72,12,6]]", 6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)])


def gross_code() -> BBCode:
    """IBM's gross code [[144,12,12]]: l=12, m=6, same A and B as bb72."""
    return BBCode("gross [[144,12,12]]", 12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)])
