"""
Tests for qdec/gf2.py and qdec/codes.py.

What these tests are FOR: proving, by computation, that the codes we build
from polynomials are the codes the papers say they are (right n, k, check
weight, commuting checks, valid logical operators), and that our GF(2)
algebra agrees with an independent implementation (ldpc's mod2 module).
"""

import itertools

import numpy as np
import pytest

from qdec import gf2
from qdec.codes import BBCode, bb72, gross_code, poly_matrix, shift_matrices, toric_code

ALL_CODES = [toric_code(3), toric_code(4), toric_code(6), bb72(), gross_code()]


# ---------------------------------------------------------------- GF(2) ----
def test_rank_small_examples():
    assert gf2.rank([[1, 1], [1, 1]]) == 1          # identical rows
    assert gf2.rank([[1, 0], [0, 1]]) == 2
    assert gf2.rank([[1, 1, 0], [0, 1, 1], [1, 0, 1]]) == 2   # row3 = row1 + row2


def test_nullspace_vectors_are_in_kernel_and_complete():
    rng = np.random.default_rng(0)
    for _ in range(20):
        M = rng.integers(0, 2, size=(6, 10))
        N = gf2.nullspace(M)
        assert not np.any((M @ N.T) % 2)                         # really in the kernel
        assert N.shape[0] == M.shape[1] - gf2.rank(M)            # rank-nullity


def test_solve_finds_solution_or_none():
    rng = np.random.default_rng(1)
    for _ in range(20):
        M = rng.integers(0, 2, size=(5, 8))
        x_true = rng.integers(0, 2, size=8)
        b = (M @ x_true) % 2
        x = gf2.solve(M, b)
        assert x is not None and np.array_equal((M @ x) % 2, b)
    # inconsistent system: x0 = 0 and x0 = 1
    assert gf2.solve([[1], [1]], [0, 1]) is None


def test_rank_matches_ldpc():
    mod2 = pytest.importorskip("ldpc.mod2")
    for code in ALL_CODES:
        assert gf2.rank(code.HX) == mod2.rank(code.HX)


# ------------------------------------------------------------- matrices ----
def test_x_and_y_commute_and_are_permutations():
    x, y = shift_matrices(5, 3)
    assert np.array_equal((x @ y) % 2, (y @ x) % 2)
    for M in (x, y):
        assert np.all(M.sum(0) == 1) and np.all(M.sum(1) == 1)


def test_polynomial_matrix_has_one_entry_per_term():
    A = poly_matrix([(3, 0), (0, 1), (0, 2)], 12, 6)
    assert np.all(A.sum(1) == 3)                      # 3 monomials -> weight 3 rows


# ---------------------------------------------------------------- codes ----
@pytest.mark.parametrize("code,n,k,w", [
    (toric_code(3), 18, 2, 4),
    (toric_code(6), 72, 2, 4),
    (bb72(), 72, 12, 6),
    (gross_code(), 144, 12, 12 // 2),          # check weight 6
])
def test_parameters(code, n, k, w):
    assert code.n == n
    assert code.k == k
    assert code.HX.sum(1).max() == w and code.HZ.sum(1).max() == w


@pytest.mark.parametrize("code", ALL_CODES, ids=lambda c: c.name)
def test_checks_commute(code):
    assert code.commutes()


def test_any_polynomials_give_commuting_checks():
    """The construction guarantees H_X H_Z^T = 0 for random A, B too."""
    rng = np.random.default_rng(2)
    for _ in range(10):
        A = [tuple(t) for t in rng.integers(0, 5, size=(3, 2))]
        B = [tuple(t) for t in rng.integers(0, 5, size=(3, 2))]
        assert BBCode("random", 5, 5, A, B).commutes()


@pytest.mark.parametrize("code", ALL_CODES, ids=lambda c: c.name)
def test_logicals_are_valid_and_paired(code):
    LX, LZ = code.logicals()
    assert LX.shape == LZ.shape == (code.k, code.n)
    # X logicals pass every Z check, Z logicals pass every X check
    assert not np.any((code.HZ.astype(int) @ LX.T) % 2)
    assert not np.any((code.HX.astype(int) @ LZ.T) % 2)
    # they are not stabilizers
    for v in LX:
        assert not gf2.in_row_space(code.HX, v)
    for v in LZ:
        assert not gf2.in_row_space(code.HZ, v)
    # canonical pairing: LX_i . LZ_j = delta_ij (mod 2)
    assert np.array_equal((LX.astype(int) @ LZ.T) % 2, np.eye(code.k, dtype=int))


def test_toric_3x3_distance_is_3_by_brute_force():
    """Smallest X logical of the 3x3 toric code has weight 3 (no weight 1 or 2 one exists)."""
    code = toric_code(3)
    found = None
    for w in range(1, 4):
        for support in itertools.combinations(range(code.n), w):
            v = np.zeros(code.n, dtype=np.uint8)
            v[list(support)] = 1
            if not np.any((code.HZ.astype(int) @ v) % 2) and not gf2.in_row_space(code.HX, v):
                found = w
                break
        if found:
            break
    assert found == 3


def test_toric_layout_is_local():
    """In the toric code every check touches only its 4 nearest grid neighbours (on the torus)."""
    code = toric_code(5)
    period = np.array(code.grid_shape())
    dpos = code.data_positions()
    for kind in "XZ":
        cpos = code.check_positions(kind)
        H = code.H(kind)
        for c in range(H.shape[0]):
            for q in np.nonzero(H[c])[0]:
                d = (dpos[q] - cpos[c]) % period
                d = np.minimum(d, period - d)          # shortest way round the torus
                assert d.sum() == 1


def test_positions_are_all_distinct():
    code = gross_code()
    allpos = np.vstack([code.data_positions(), code.check_positions("X"), code.check_positions("Z")])
    assert len({tuple(p) for p in allpos}) == len(allpos) == 4 * code.l * code.m


@pytest.mark.parametrize("code", [bb72(), gross_code()], ids=lambda c: c.name)
def test_bb_toric_layout_has_four_short_edges_per_check(code):
    """Bravyi et al.'s toric layout: every check has 4 nearest-neighbour and 2 long-range edges."""
    from qdec.tanner import torus_displacement
    assert code.layout()["kind"] == "toric"
    period = np.array(code.grid_shape())
    dpos = code.data_positions()
    for kind in "XZ":
        cpos, H = code.check_positions(kind), code.H(kind)
        for c in range(H.shape[0]):
            lens = [np.abs(torus_displacement(cpos[c], dpos[q], period)).sum() for q in np.nonzero(H[c])[0]]
            assert sorted(lens)[:4] == [1, 1, 1, 1] and min(sorted(lens)[4:]) > 1
