"""
families.py - every code in the catalogue, written from its GEOMETRY.

What this file is FOR
---------------------
Each function here starts from a picture - qubits on a line, on a grid, on a
triangle, on a pentagon, on a torus - and reads the stabilizers off that
picture.  The check matrices are a *consequence* of the geometry.  That is the
whole point of the project: look at the drawing and recognise the code.

Catalogue (all return a core.StabilizerCode)
--------------------------------------------
    repetition(d)              bit-flip repetition code on a line
    shor()                     [[9,1,3]]   3 x 3 grid, ZZ pairs in rows, X^6 across rows
    steane()                   [[7,1,3]]   triangle with 3 faces (smallest colour code)
    five_qubit()               [[5,1,3]]   pentagon, XZZXI shifted around the ring (not CSS)
    surface_unrotated(d)       [[d^2+(d-1)^2,1,d]]  vertices = X checks, faces = Z checks
    surface_rotated(d)         [[d^2,1,d]]  checkerboard plaquettes, weight-2 boundaries
    toric(L)                   [[2L^2,2,L]] surface code with no boundary (on a torus)
    hypergraph_product(H1,H2)  Tillich-Zemor product of two classical codes (qLDPC family)
    bb(name)                   bivariate-bicycle codes incl. the gross code (Bravyi et al.)

`catalogue()` returns the standard list used by the notebooks, tests and explorer.

Coordinates: (x, y) with y growing downward, 1 unit = one grid step.
"""

from __future__ import annotations

import numpy as np

from .codes import BBCode
from .core import Geometry, StabilizerCode

REF_GOTTESMAN = "Gottesman, PhD thesis (1997), arXiv:quant-ph/9705052"
REF_NC = "Nielsen & Chuang, Quantum Computation and Quantum Information, ch. 10"


def _centroid(points):
    return np.asarray(points, dtype=float).mean(axis=0)


# --------------------------------------------------------------------------
# 1D: repetition code
# --------------------------------------------------------------------------
def repetition(d: int = 5) -> StabilizerCode:
    """
    d qubits on a line; generator i = Z_i Z_{i+1} (sits between the two qubits).
    Protects against bit flips only: X distance d, Z distance 1.
    """
    HZ = np.zeros((d - 1, d), np.uint8)
    for i in range(d - 1):
        HZ[i, [i, i + 1]] = 1
    qpos = np.array([[2.0 * i, 0.0] for i in range(d)])
    cpos = np.array([[2.0 * i + 1, 0.0] for i in range(d - 1)])
    geo = Geometry(qpos, cpos, ["polygon"] * (d - 1))
    return StabilizerCode.css(f"repetition d={d}", np.zeros((0, d), np.uint8), HZ, geo, dict(
        family="repetition", d_known=1, dX=d, dZ=1,
        description="Qubits in a row; each check compares two neighbours (Z_i Z_i+1). "
                    "Catches bit flips (X) only - a single Z flip is a logical error, so as a "
                    "quantum code d = 1, but against X errors the distance is d.",
        references=[REF_NC + " (sec. 10.1)"]))


# --------------------------------------------------------------------------
# small CSS codes
# --------------------------------------------------------------------------
def shor() -> StabilizerCode:
    """
    Shor's [[9,1,3]] code on a 3 x 3 grid: row b = block b.
    Z checks: neighbouring pairs inside a row (6 of them).
    X checks: all six qubits of two neighbouring rows (2 of them).
    """
    q = lambda b, j: 3 * b + j
    HZ, cz = [], []
    for b in range(3):
        for j in range(2):
            row = np.zeros(9, np.uint8); row[[q(b, j), q(b, j + 1)]] = 1
            HZ.append(row); cz.append([2.0 * j + 1, 2.0 * b])
    HX, cx = [], []
    for b in range(2):
        row = np.zeros(9, np.uint8); row[[q(b, j) for j in range(3)] + [q(b + 1, j) for j in range(3)]] = 1
        HX.append(row); cx.append([2.0, 2.0 * b + 1])
    qpos = np.array([[2.0 * j, 2.0 * b] for b in range(3) for j in range(3)])
    geo = Geometry(qpos, np.array(cx + cz), ["polygon"] * 8)
    return StabilizerCode.css("Shor [[9,1,3]]", np.array(HX), np.array(HZ), geo, dict(
        family="small CSS", d_known=3,
        description="Three blocks of three (rows). Inside a row, pairwise Z checks catch a bit flip "
                    "(3-qubit repetition code). Across two rows, a six-qubit X check compares the "
                    "phase of the blocks, catching a phase flip. Degenerate: Z0 and Z1 have the "
                    "same syndrome and either correction works.",
        references=["Shor, Phys. Rev. A 52, R2493 (1995)", REF_NC + " (sec. 10.2)"]))


def steane() -> StabilizerCode:
    """
    Steane's [[7,1,3]] code drawn as the smallest triangular colour code:
    7 qubits on a triangle, 3 four-qubit faces; every face carries BOTH an
    X check and a Z check (so H_X = H_Z = Hamming [7,4,3] parity checks).
    """
    h = np.sqrt(3)
    qpos = np.array([[0, 2 * h], [2, 2 * h], [4, 2 * h],      # bottom: corner, middle, corner
                     [1, h], [2, 4 * h / 3], [3, h],          # left side, centre, right side
                     [2, 0]])                                 # top corner
    faces = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 5, 6]]
    H = np.zeros((3, 7), np.uint8)
    for f, face in enumerate(faces):
        H[f, face] = 1
    cent = [_centroid(qpos[f]) for f in faces]
    cpos = np.array([c + [-0.25, 0] for c in cent] + [c + [0.25, 0] for c in cent])   # X label left, Z right
    geo = Geometry(qpos.astype(float), cpos, ["polygon"] * 6)
    return StabilizerCode.css("Steane [[7,1,3]]", H, H, geo, dict(
        family="small CSS / colour code", d_known=3,
        description="Seven qubits on a triangle with three four-qubit faces. Each face is measured "
                    "twice, once as XXXX and once as ZZZZ. The three faces are the parity checks of "
                    "the classical Hamming [7,4,3] code, so every single error gives a different "
                    "pattern of faces. This is the smallest 2D colour code.",
        references=["Steane, Phys. Rev. Lett. 77, 793 (1996)",
                    "Bombin & Martin-Delgado, Phys. Rev. Lett. 97, 180501 (2006) (colour codes)"]))


def five_qubit() -> StabilizerCode:
    """
    The perfect [[5,1,3]] code: generators XZZXI and its cyclic shifts
    (4 independent ones).  Not CSS - each generator mixes X and Z - so it
    is drawn on a pentagon with the Pauli letter written at each qubit.
    """
    base = "XZZXI"
    n = 5
    Sx, Sz, cpos = [], [], []
    ang = lambda i: -np.pi / 2 + 2 * np.pi * i / n
    qpos = np.array([[2 * np.cos(ang(i)), 2 * np.sin(ang(i))] for i in range(n)]) + 2.5
    for s in range(4):
        word = base[-s:] + base[:-s] if s else base            # cyclic shift by s
        Sx.append([1 if c in "XY" else 0 for c in word])
        Sz.append([1 if c in "ZY" else 0 for c in word])
        supp = [i for i, c in enumerate(word) if c != "I"]
        cpos.append(_centroid(qpos[supp]) * 0.75 + 2.5 * 0.25)   # pulled slightly to the centre
    geo = Geometry(qpos, np.array(cpos), ["polygon"] * 4)
    return StabilizerCode("five-qubit [[5,1,3]]", np.array(Sx), np.array(Sz), geo, dict(
        family="perfect (non-CSS)", d_known=3,
        description="Five qubits on a ring. One stabilizer, XZZXI, is shifted one step around the "
                    "ring four times. It mixes X and Z in the same check, so it is not CSS. It is "
                    "the smallest code that corrects any single-qubit error, and it is 'perfect': "
                    "the 15 single errors plus 'no error' use up all 2^4 = 16 syndromes exactly.",
        references=["Laflamme, Miquel, Paz, Zurek, PRL 77, 198 (1996)",
                    "Bennett, DiVincenzo, Smolin, Wootters, PRA 54, 3824 (1996)", REF_GOTTESMAN]))


# --------------------------------------------------------------------------
# surface codes
# --------------------------------------------------------------------------
def surface_unrotated(d: int = 3) -> StabilizerCode:
    """
    Planar (unrotated) surface code on a (2d-1) x (2d-1) grid:
        data qubits     at (r, c) with r + c even      -> d^2 + (d-1)^2 qubits
        X checks        at (r even, c odd)  ("vertices") touch up/down/left/right
        Z checks        at (r odd,  c even) ("faces")
    Checks on the edge of the grid have only 3 qubits: that is the boundary.
    """
    size = 2 * d - 1
    data = [(r, c) for r in range(size) for c in range(size) if (r + c) % 2 == 0]
    idx = {p: i for i, p in enumerate(data)}
    n = len(data)

    def checks(rpar, cpar):
        rows, pos = [], []
        for r in range(size):
            for c in range(size):
                if r % 2 == rpar and c % 2 == cpar:
                    row = np.zeros(n, np.uint8)
                    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        if (r + dr, c + dc) in idx:
                            row[idx[(r + dr, c + dc)]] = 1
                    rows.append(row); pos.append([float(c), float(r)])
        return np.array(rows), pos

    HX, px = checks(0, 1)
    HZ, pz = checks(1, 0)
    qpos = np.array([[float(c), float(r)] for r, c in data])
    geo = Geometry(qpos, np.array(px + pz), ["polygon"] * (len(px) + len(pz)))
    return StabilizerCode.css(f"surface (unrotated) d={d}", HX, HZ, geo, dict(
        family="surface", d_known=d,
        description="Kitaev's surface code with boundaries. Qubits sit on the edges of a square "
                    "grid. X checks sit on vertices (the star of 4 edges) and Z checks on faces "
                    "(the 4 edges around a square). A logical X is a string of X's crossing the "
                    "grid; any shorter string leaves two lit checks at its ends.",
        references=["Kitaev, Ann. Phys. 303, 2 (2003)", "Bravyi & Kitaev, arXiv:quant-ph/9811052",
                    "Dennis, Kitaev, Landahl, Preskill, J. Math. Phys. 43, 4452 (2002)"]))


def surface_rotated(d: int = 3) -> StabilizerCode:
    """
    Rotated surface code: d x d data qubits at odd grid points (2i+1, 2j+1);
    plaquettes at even points (2i, 2j), coloured like a chessboard (X if i+j even).
    Inside: all 4-qubit plaquettes.  Boundary: only weight-2 plaquettes of the
    right type (X on top/bottom, Z on left/right) - corners have none.
    Same construction as Stim's 'surface_code:rotated_memory_*' (Horsman et al.).
    """
    data = [(2 * i + 1, 2 * j + 1) for i in range(d) for j in range(d)]
    idx = {p: k for k, p in enumerate(data)}
    n = d * d
    HX, HZ, px, pz = [], [], [], []
    for i in range(d + 1):
        for j in range(d + 1):
            kind = "X" if (i + j) % 2 == 0 else "Z"
            interior = 1 <= i <= d - 1 and 1 <= j <= d - 1
            top_bottom = (i in (0, d)) and 1 <= j <= d - 1 and kind == "X"
            left_right = (j in (0, d)) and 1 <= i <= d - 1 and kind == "Z"
            if not (interior or top_bottom or left_right):
                continue
            row = np.zeros(n, np.uint8)
            for di in (-1, 1):
                for dj in (-1, 1):
                    p = (2 * i + di, 2 * j + dj)
                    if p in idx:
                        row[idx[p]] = 1
            (HX if kind == "X" else HZ).append(row)
            (px if kind == "X" else pz).append([2.0 * j, 2.0 * i])
    qpos = np.array([[float(c), float(r)] for r, c in data])
    geo = Geometry(qpos, np.array(px + pz), ["polygon"] * (len(px) + len(pz)))
    return StabilizerCode.css(f"surface (rotated) d={d}", np.array(HX), np.array(HZ), geo, dict(
        family="surface", d_known=d,
        # hook-safe CNOT order (direction signs (dx, dy), y down): found by searching all
        # 24 x 24 orders with circuits.circuit_distance; this is the standard pair
        # "Z-shape for X checks, N-shape for Z checks" (Tomita & Svore 2014, Fig. 2)
        cnot_order={"X": [(-1, -1), (1, -1), (-1, 1), (1, 1)],
                    "Z": [(-1, -1), (-1, 1), (1, -1), (1, 1)]},
        description="The compact surface code used by Google and most experiments. d x d data "
                    "qubits, plaquettes coloured like a chessboard (X or Z). Weight-2 plaquettes on "
                    "the boundary: X ones on two opposite sides, Z ones on the other two. That "
                    "uses about half the qubits of the unrotated code for the same distance.",
        references=["Horsman, Fowler, Devitt, Van Meter, New J. Phys. 14, 123011 (2012)",
                    "Tomita & Svore, Phys. Rev. A 90, 062320 (2014)",
                    "Google Quantum AI, Nature 614, 676 (2023)"]))


def toric(L: int = 4) -> StabilizerCode:
    """Toric code = surface code on a torus (no boundary): BB construction with A = 1+x, B = 1+y."""
    return _from_bb(BBCode(f"toric {L}x{L}", L, L, [(0, 0), (1, 0)], [(0, 0), (0, 1)]), dict(
        family="surface (toric)", d_known=L,
        description="The surface code with opposite edges glued together (a torus), so there is "
                    "no boundary. Every check has 4 qubits. Two logical qubits, from the two "
                    "independent ways a string can wrap around the torus. One X check and one Z "
                    "check are redundant: the product of all X checks is the identity.",
        references=["Kitaev, Ann. Phys. 303, 2 (2003)", "Dennis et al., J. Math. Phys. 43, 4452 (2002)"]),
        shape="polygon")


# --------------------------------------------------------------------------
# qLDPC: hypergraph product and bivariate bicycle
# --------------------------------------------------------------------------
def classical_repetition_H(d: int, closed: bool = False) -> np.ndarray:
    """Parity checks of the classical repetition code (open line or closed ring)."""
    r = d if closed else d - 1
    H = np.zeros((r, d), np.uint8)
    for i in range(r):
        H[i, [i, (i + 1) % d]] = 1
    return H


def hamming_H() -> np.ndarray:
    """[7,4,3] Hamming code: column j = binary of j+1."""
    return np.array([[(j + 1) >> b & 1 for j in range(7)] for b in range(3)], np.uint8)


def hypergraph_product(H1, H2, name: str = "hypergraph product", info: dict | None = None) -> StabilizerCode:
    """
    Tillich-Zemor hypergraph product of two classical codes H1 (r1 x n1), H2 (r2 x n2):

        H_X = [ H1 (x) I_n2  |  I_r1 (x) H2^T ]
        H_Z = [ I_n1 (x) H2  |  H1^T (x) I_r2 ]

    Qubits: n1*n2 "bit x bit" + r1*r2 "check x check".   H_X H_Z^T = 2 H1 (x) H2^T = 0.
    Geometry: bit i of a code sits at coordinate 2i, check j at 2j+1, on both axes.
    The product of two open repetition codes is exactly the unrotated surface code.
    """
    H1, H2 = np.asarray(H1, np.uint8), np.asarray(H2, np.uint8)
    r1, n1 = H1.shape
    r2, n2 = H2.shape
    I = lambda k: np.eye(k, dtype=np.uint8)
    HX = np.hstack([np.kron(H1, I(n2)), np.kron(I(r1), H2.T)]) % 2
    HZ = np.hstack([np.kron(I(n1), H2), np.kron(H1.T, I(r2))]) % 2
    qpos = [[2.0 * b, 2.0 * a] for a in range(n1) for b in range(n2)] + \
           [[2.0 * c + 1, 2.0 * a + 1] for a in range(r1) for c in range(r2)]
    px = [[2.0 * b, 2.0 * a + 1] for a in range(r1) for b in range(n2)]     # X check (c1, v2)
    pz = [[2.0 * c + 1, 2.0 * a] for a in range(n1) for c in range(r2)]     # Z check (v1, c2)
    geo = Geometry(np.array(qpos), np.array(px + pz), ["star"] * (len(px) + len(pz)))
    base = dict(family="hypergraph product (qLDPC)",
                description="Product of two classical codes laid out as a grid: bit x bit qubits and "
                            "check x check qubits. Checks are local along rows and columns only if "
                            "the classical codes are local. Feed in better classical LDPC codes and "
                            "you get more logical qubits at the same distance.",
                references=["Tillich & Zemor, IEEE Trans. Inf. Theory 60, 1193 (2014), arXiv:0903.0566"])
    base.update(info or {})
    return StabilizerCode.css(name, HX, HZ, geo, base)


BB_TABLE = {
    # name: (l, m, A terms, B terms, d)   - Bravyi et al., Nature 627, 778 (2024), Table 3
    "bb72":  (6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 6),
    "bb90":  (15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)], 10),
    "bb108": (9, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 10),
    "gross": (12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)], 12),
    "bb288": (12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)], 18),
}


def bb(name: str = "gross") -> StabilizerCode:
    """A bivariate-bicycle code from Bravyi et al.'s table (bb72, bb90, bb108, gross, bb288)."""
    l, m, A, B, d = BB_TABLE[name]
    label = "gross [[144,12,12]]" if name == "gross" else f"BB [[{2 * l * m},?,{d}]]"
    code = BBCode(label, l, m, A, B)
    from .codes import poly_str
    out = _from_bb(code, dict(
        family="bivariate bicycle (qLDPC)", d_known=d, l=l, m=m,
        A=poly_str(A), B=poly_str(B),
        description=f"Torus of {l} x {m} cells, each holding an L and an R data qubit, an X check "
                    f"and a Z check. A = {poly_str(A)}, B = {poly_str(B)}: each term is one 'arrow' "
                    "from a check to a qubit, the same arrow in every cell. In the toric layout 4 "
                    "arrows reach nearest neighbours and 2 are long-range: the toric code plus "
                    "two long wires per check.",
        references=["Bravyi et al., Nature 627, 778 (2024), arXiv:2308.07915",
                    "Kovalev & Pryadko, PRA 88, 012311 (2013) (generalised bicycle codes)"]),
        shape="star")
    if name != "gross":
        out.name = f"BB [[{out.n},{out.k},{d}]]"
    return out


def _from_bb(code: BBCode, info: dict, shape: str) -> StabilizerCode:
    """Convert a BBCode (codes.py, torus layout) into a StabilizerCode with a wrapping geometry."""
    qrc = code.data_positions().astype(float)                      # (row, col)
    xrc = code.check_positions("X").astype(float)
    zrc = code.check_positions("Z").astype(float)
    rows, cols = code.grid_shape()
    geo = Geometry(qrc[:, ::-1].copy(), np.vstack([xrc, zrc])[:, ::-1].copy(),
                   [shape] * (2 * code.l * code.m), period=(float(cols), float(rows)))
    return StabilizerCode.css(code.name, code.HX, code.HZ, geo, info)


# --------------------------------------------------------------------------
# the standard list
# --------------------------------------------------------------------------
def catalogue() -> list[StabilizerCode]:
    """Codes used across notebooks, tests and the browser explorer (small to large)."""
    hp_surface = hypergraph_product(classical_repetition_H(3), classical_repetition_H(3),
                                    "HGP(rep3, rep3) = surface d=3",
                                    dict(d_known=3, description="The hypergraph product of two "
                                         "3-bit repetition codes. It has exactly the parameters of "
                                         "the unrotated surface code [[13,1,3]]: the surface code "
                                         "IS a hypergraph product."))
    hp_hamming = hypergraph_product(hamming_H(), hamming_H(), "HGP(Hamming, Hamming)",
                                    dict(d_known=3, description="Product of two [7,4,3] Hamming "
                                         "codes: [[58,16,3]]. 16 logical qubits for 58 physical, "
                                         "against 1 logical for 13 in the surface code of the same "
                                         "distance. The price is non-local checks."))
    return [repetition(5), shor(), steane(), five_qubit(),
            surface_unrotated(3), surface_unrotated(5), surface_rotated(3), surface_rotated(5),
            toric(4), hp_surface, hp_hamming,
            bb("bb72"), bb("bb90"), bb("bb108"), bb("gross"), bb("bb288")]
