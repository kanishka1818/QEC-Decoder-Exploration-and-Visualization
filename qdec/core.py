"""
core.py - one object that every quantum code in this project becomes.

What this file is FOR
---------------------
Repetition, Shor, Steane, the 5-qubit code, surface codes, toric, hypergraph
products and bivariate-bicycle codes look very different on paper.  Here
they all become the SAME kind of object, so every picture, circuit and
decoder works on every code:

    StabilizerCode
      .Sx, .Sz      the stabilizer generators in "symplectic" form:
                    generator i = product over qubits q of X^Sx[i,q] Z^Sz[i,q]
                    (X and Z on the same qubit = Y)
      .kinds        'X', 'Z' or 'mixed' per generator
      .geometry     where every qubit and every generator is drawn
      .info         family, plain description, references

From these we compute (all over GF(2), with our own gf2.py):
    n, k, whether it is CSS, H_X / H_Z, logical operators (paired),
    the syndrome of any Pauli error, and - for small codes - the distance
    by brute force.

Why symplectic?  A Pauli on n qubits is a pair of 0/1 vectors (x | z).
Two Paulis commute exactly when x1.z2 + z1.x2 = 0 (mod 2).  That single
formula covers CSS and non-CSS codes (the 5-qubit code) alike.
Reference: Gottesman, PhD thesis (1997), arXiv:quant-ph/9705052, ch. 3;
Nielsen & Chuang, sec. 10.5.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from .gf2 import as_gf2, in_row_space, nullspace, rank, row_reduce


# --------------------------------------------------------------------------
# geometry: where things are drawn
# --------------------------------------------------------------------------
@dataclass
class Geometry:
    """
    Drawing information.  Coordinates are (x, y) with y growing DOWN the
    page (screen convention), the same in matplotlib views and the browser.

    qubit_pos : (n, 2)  position of every data qubit
    check_pos : (m, 2)  position of every stabilizer generator's "label"
                        (a plaquette centre, a vertex, or an ancilla spot)
    shape     : per generator, how to draw it:
                'polygon' = filled face through its qubits (local codes)
                'star'    = node at check_pos with a spoke to each qubit
                            (codes with long-range checks, e.g. gross)
    period    : (Px, Py) if the picture wraps around (torus), else None
    """

    qubit_pos: np.ndarray
    check_pos: np.ndarray
    shape: list[str]
    period: tuple[float, float] | None = None

    def nearest_image(self, origin, point) -> np.ndarray:
        """Copy of `point` closest to `origin` when the picture wraps (torus)."""
        point = np.asarray(point, dtype=float)
        if self.period is None:
            return point
        P = np.asarray(self.period, dtype=float)
        d = (point - np.asarray(origin, dtype=float) + P / 2) % P - P / 2
        return np.asarray(origin, dtype=float) + d


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def symplectic_product(a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    """
    Commutation table between Pauli rows a (r1 x 2n) and b (r2 x 2n):
    entry 1 = anticommute, 0 = commute.   ax.bz + az.bx  (mod 2)
    """
    a = np.atleast_2d(a).astype(np.int64)
    b = np.atleast_2d(b).astype(np.int64)
    return (a[:, :n] @ b[:, n:].T + a[:, n:] @ b[:, :n].T) % 2


def pauli_string(x: np.ndarray, z: np.ndarray) -> str:
    """0/1 vectors -> 'XIZY...' string."""
    return "".join("IXZY"[int(a) + 2 * int(b)] for a, b in zip(x, z))


def _independent_mod(span: np.ndarray, candidates: np.ndarray) -> np.ndarray:
    """Keep candidates that add something new on top of `span` (quotient representatives)."""
    kept, current = [], as_gf2(span)
    r = rank(current) if current.size else 0
    for v in candidates:
        trial = np.vstack([current, v[None, :]]) if current.size else v[None, :]
        r_new = rank(trial)
        if r_new > r:
            kept.append(v)
            current, r = trial, r_new
    return np.array(kept, dtype=np.uint8).reshape(-1, candidates.shape[1])


def _gf2_inverse(G: np.ndarray) -> np.ndarray:
    k = G.shape[0]
    R, piv = row_reduce(np.hstack([as_gf2(G), np.eye(k, dtype=np.uint8)]))
    if piv[:k] != list(range(k)):
        raise ValueError("not invertible over GF(2)")
    return R[:, k:]


# --------------------------------------------------------------------------
# the code object
# --------------------------------------------------------------------------
@dataclass
class StabilizerCode:
    name: str
    Sx: np.ndarray                       # (m, n) X part of each generator
    Sz: np.ndarray                       # (m, n) Z part of each generator
    geometry: Geometry
    info: dict = field(default_factory=dict)   # family, description, references, d_known ...

    def __post_init__(self):
        self.Sx = as_gf2(self.Sx)
        self.Sz = as_gf2(self.Sz)
        self._logicals = None

    # ---- construction shortcut for CSS codes ------------------------------
    @classmethod
    def css(cls, name, HX, HZ, geometry, info=None):
        """Build a CSS code: X generators = rows of HX, Z generators = rows of HZ."""
        HX, HZ = as_gf2(HX), as_gf2(HZ)
        n = HX.shape[1]
        Sx = np.vstack([HX, np.zeros((HZ.shape[0], n), np.uint8)])
        Sz = np.vstack([np.zeros((HX.shape[0], n), np.uint8), HZ])
        return cls(name, Sx, Sz, geometry, info or {})

    # ---- basic numbers ----------------------------------------------------
    @property
    def n(self) -> int:
        return self.Sx.shape[1]

    @property
    def m(self) -> int:
        """Number of listed generators (may include dependent ones, e.g. toric code)."""
        return self.Sx.shape[0]

    @property
    def S(self) -> np.ndarray:
        """(m, 2n) symplectic matrix [Sx | Sz]."""
        return np.hstack([self.Sx, self.Sz])

    @property
    def kinds(self) -> list[str]:
        """'X', 'Z' or 'mixed' for every generator."""
        out = []
        for x, z in zip(self.Sx, self.Sz):
            out.append("X" if not z.any() else "Z" if not x.any() else "mixed")
        return out

    @property
    def is_css(self) -> bool:
        return "mixed" not in self.kinds

    @property
    def HX(self) -> np.ndarray:
        """Rows of the X-type generators (CSS codes)."""
        idx = [i for i, t in enumerate(self.kinds) if t == "X"]
        return self.Sx[idx]

    @property
    def HZ(self) -> np.ndarray:
        idx = [i for i, t in enumerate(self.kinds) if t == "Z"]
        return self.Sz[idx]

    def index_of(self, kind: str) -> list[int]:
        """Generator indices of a given kind ('X' or 'Z'); HX row r = generator index_of('X')[r]."""
        return [i for i, t in enumerate(self.kinds) if t == kind]

    @property
    def k(self) -> int:
        """k = n - rank(S): each independent generator removes one qubit's worth of freedom."""
        return self.n - rank(self.S)

    def commutes(self) -> bool:
        """All generators must commute with each other."""
        return not symplectic_product(self.S, self.S, self.n).any()

    def generator_string(self, i: int) -> str:
        return pauli_string(self.Sx[i], self.Sz[i])

    def support(self, i: int) -> np.ndarray:
        """Qubits that generator i acts on."""
        return np.nonzero(self.Sx[i] | self.Sz[i])[0]

    # ---- syndromes --------------------------------------------------------
    def syndrome(self, ex=None, ez=None) -> np.ndarray:
        """
        Which generators an error X^ex Z^ez flips (1 = anticommutes).
        For CSS codes: X errors (ex) are seen by Z generators, Z errors by X generators.
        """
        ex = np.zeros(self.n, np.uint8) if ex is None else as_gf2(ex)
        ez = np.zeros(self.n, np.uint8) if ez is None else as_gf2(ez)
        return symplectic_product(self.S, np.concatenate([ex, ez])[None, :], self.n)[:, 0].astype(np.uint8)

    # ---- logical operators ------------------------------------------------
    def logicals(self) -> tuple[np.ndarray, np.ndarray]:
        """
        (LX, LZ): k pairs of logical operators as (k, 2n) symplectic rows.
        Pair i anticommutes with itself' partner only: LX_i * LZ_j anticommute iff i == j,
        and LX's commute among themselves, LZ's among themselves.

        CSS codes: LX are X-type (null space of H_Z, not in row space of H_X),
        LZ are Z-type - the familiar picture.  Other codes: general symplectic
        Gram-Schmidt on the normaliser (things that commute with every generator).
        """
        if self._logicals is not None:
            return self._logicals
        n = self.n
        if self.is_css:
            HX, HZ = self.HX, self.HZ
            ker = lambda H: nullspace(H) if H.size else np.eye(n, dtype=np.uint8)
            span = lambda H: H if H.size else np.zeros((0, n), np.uint8)
            lx = _independent_mod(span(HX), ker(HZ))     # X logicals
            lz = _independent_mod(span(HZ), ker(HX))     # Z logicals
            G = (lx.astype(np.int64) @ lz.T.astype(np.int64)) % 2
            lz = (_gf2_inverse(G).T.astype(np.int64) @ lz.astype(np.int64)) % 2
            z = np.zeros_like(lx)
            LX = np.hstack([lx, z]).astype(np.uint8)
            LZ = np.hstack([np.zeros_like(lz), lz]).astype(np.uint8)
        else:
            # normaliser = {v : v commutes with all generators} = nullspace of [Sz | Sx]
            N = nullspace(np.hstack([self.Sz, self.Sx]))
            reps = list(_independent_mod(self.S, N))
            LXl, LZl = [], []
            sp = lambda a, b: int(symplectic_product(a, b, n)[0, 0])
            while reps:                                   # symplectic Gram-Schmidt
                a = reps.pop(0)
                j = next(i for i, b in enumerate(reps) if sp(a, b))
                b = reps.pop(j)
                LXl.append(a); LZl.append(b)
                reps = [(c + sp(c, b) * a + sp(c, a) * b) % 2 for c in reps]
            LX, LZ = np.array(LXl, np.uint8), np.array(LZl, np.uint8)
        self._logicals = (LX, LZ)
        return self._logicals

    # ---- distance (brute force, small codes only) ---------------------------
    def distance(self, max_weight: int = 12, max_candidates: int = 3_000_000) -> int | None:
        """
        Smallest weight of a logical operator (a Pauli that commutes with every
        generator but is not a stabilizer).  Brute force over supports of
        growing size; returns None if nothing found up to max_weight or the
        search would exceed max_candidates.

        Test used: a Pauli is a NON-TRIVIAL logical iff it commutes with all
        generators and anticommutes with at least one of the 2k logicals.
        """
        n = self.n
        LX, LZ = self.logicals()
        allL = np.vstack([LX, LZ])
        if self.is_css:
            # X-type and Z-type parts can be searched separately: d = min(dX, dZ)
            best = None
            for H_other, Lopp, kind in [(self.HZ, LZ[:, n:], "X"), (self.HX, LX[:, :n], "Z")]:
                dk = self._css_min_weight(H_other, Lopp, max_weight, max_candidates)
                if dk is None:
                    return None
                best = dk if best is None else min(best, dk)
            return best
        total = 0
        for w in range(1, max_weight + 1):
            for supp in itertools.combinations(range(n), w):
                for paulis in itertools.product((1, 2, 3), repeat=w):   # X, Z, Y on each qubit
                    total += 1
                    if total > max_candidates:
                        return None
                    v = np.zeros(2 * n, np.uint8)
                    for q, p in zip(supp, paulis):
                        v[q] = p & 1
                        v[n + q] = p >> 1
                    if symplectic_product(self.S, v, n).any():
                        continue
                    if symplectic_product(allL, v, n).any():
                        return w
        return None

    def css_distances(self, max_weight: int = 12, max_candidates: int = 3_000_000):
        """(d_X, d_Z) for a CSS code: lightest X-type logical and lightest Z-type logical."""
        n = self.n
        LX, LZ = self.logicals()
        dX = self._css_min_weight(self.HZ, LZ[:, n:], max_weight, max_candidates)
        dZ = self._css_min_weight(self.HX, LX[:, :n], max_weight, max_candidates)
        return dX, dZ

    @staticmethod
    def _css_min_weight(H_other, L_opp, max_weight, max_candidates):
        """Min weight of v with H_other v = 0 and L_opp v != 0 (vectorised over supports)."""
        n = H_other.shape[1]
        H = H_other.astype(np.int64)
        L = L_opp.astype(np.int64)
        total = 0
        for w in range(1, max_weight + 1):
            combos = itertools.combinations(range(n), w)
            while True:
                chunk = list(itertools.islice(combos, 200_000))
                if not chunk:
                    break
                total += len(chunk)
                if total > max_candidates:
                    return None
                idx = np.array(chunk)
                # column sums mod 2 for each support = H v
                Hv = H[:, idx].sum(axis=2) % 2            # (rows, combos)
                Lv = L[:, idx].sum(axis=2) % 2
                ok = (~Hv.any(axis=0)) & Lv.any(axis=0)
                if ok.any():
                    return w
        return None

    # ---- text summary ------------------------------------------------------
    def params(self) -> str:
        d = self.info.get("d_known")
        return f"[[{self.n},{self.k},{d if d is not None else '?'}]]"

    def summary(self) -> str:
        weights = (self.Sx | self.Sz).sum(1)
        return (f"{self.name}  {self.params()}  family={self.info.get('family', '?')}  "
                f"CSS={self.is_css}  generators={self.m}  weights {weights.min()}-{weights.max()}")
