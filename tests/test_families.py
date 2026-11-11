"""
Tests for qdec/core.py and qdec/families.py.

What these tests are FOR: every code in the catalogue is checked by
computation, not trusted from the paper - commuting generators, the right
number of qubits n and logical qubits k, logical operators that really are
logical, and (where brute force is affordable) the right distance d.
"""

import numpy as np
import pytest

from qdec.core import StabilizerCode, symplectic_product
from qdec.families import (BB_TABLE, bb, catalogue, classical_repetition_H, five_qubit, hamming_H,
                           hypergraph_product, repetition, shor, steane, surface_rotated,
                           surface_unrotated, toric)

CAT = catalogue()
EXPECTED = {   # name -> (n, k)
    "repetition d=5": (5, 1), "Shor [[9,1,3]]": (9, 1), "Steane [[7,1,3]]": (7, 1),
    "five-qubit [[5,1,3]]": (5, 1), "surface (unrotated) d=3": (13, 1), "surface (unrotated) d=5": (41, 1),
    "surface (rotated) d=3": (9, 1), "surface (rotated) d=5": (25, 1), "toric 4x4": (32, 2),
    "HGP(rep3, rep3) = surface d=3": (13, 1), "HGP(Hamming, Hamming)": (58, 16),
    "BB [[72,12,6]]": (72, 12), "BB [[90,8,10]]": (90, 8), "BB [[108,8,10]]": (108, 8),
    "gross [[144,12,12]]": (144, 12), "BB [[288,12,18]]": (288, 12),
}


@pytest.mark.parametrize("code", CAT, ids=lambda c: c.name)
def test_generators_commute(code):
    assert code.commutes()


@pytest.mark.parametrize("code", CAT, ids=lambda c: c.name)
def test_n_and_k(code):
    assert (code.n, code.k) == EXPECTED[code.name]


@pytest.mark.parametrize("code", CAT, ids=lambda c: c.name)
def test_logicals(code):
    """Logicals commute with all generators, are not stabilizers, and pair up canonically."""
    LX, LZ = code.logicals()
    n, k = code.n, code.k
    assert LX.shape == LZ.shape == (k, 2 * n)
    assert not symplectic_product(code.S, LX, n).any()
    assert not symplectic_product(code.S, LZ, n).any()
    assert np.array_equal(symplectic_product(LX, LZ, n), np.eye(k, dtype=int))
    assert not symplectic_product(LX, LX, n).any() and not symplectic_product(LZ, LZ, n).any()


@pytest.mark.parametrize("code,d", [(shor(), 3), (steane(), 3), (five_qubit(), 3),
                                    (surface_unrotated(3), 3), (surface_unrotated(5), 5),
                                    (surface_rotated(3), 3), (surface_rotated(5), 5),
                                    (toric(4), 4)], ids=lambda x: getattr(x, "name", str(x)))
def test_distance_by_brute_force(code, d):
    assert code.distance(max_weight=d + 1) == d


def test_repetition_protects_against_x_only():
    assert repetition(5).css_distances() == (5, 1)


def test_five_qubit_is_not_css_and_perfect():
    c = five_qubit()
    assert not c.is_css
    # perfect: the 15 single-qubit errors + identity give all 16 syndromes exactly once
    synd = {tuple(np.zeros(4, int))}
    for q in range(5):
        for ex, ez in [(1, 0), (0, 1), (1, 1)]:
            x = np.zeros(5, int); z = np.zeros(5, int); x[q] = ex; z[q] = ez
            synd.add(tuple(c.syndrome(x, z)))
    assert len(synd) == 16


def test_steane_faces_are_hamming_checks():
    c = steane()
    cols = {tuple(col) for col in c.HX.T}
    assert len(cols) == 7 and (0, 0, 0) not in cols          # 7 distinct non-zero columns


def test_hgp_of_open_repetition_is_unrotated_surface():
    hp = hypergraph_product(classical_repetition_H(4), classical_repetition_H(4))
    s = surface_unrotated(4)
    assert (hp.n, hp.k) == (s.n, s.k) == (25, 1)
    assert hp.distance(max_weight=5) == 4


def test_hgp_hamming():
    hp = hypergraph_product(hamming_H(), hamming_H())
    assert (hp.n, hp.k) == (58, 16) and hp.distance(max_weight=4) == 3


@pytest.mark.parametrize("name", list(BB_TABLE))
def test_bb_family_found_logical_not_lighter_than_published_distance(name):
    """Brute force is too slow for BB distances; we at least check our found logicals
    are never lighter than the published d (a lighter one would disprove the table)."""
    c = bb(name)
    LX, LZ = c.logicals()
    n = c.n
    w = min(int((L[:n] | L[n:]).sum()) for L in np.vstack([LX, LZ]))
    assert w >= BB_TABLE[name][4]


def test_syndrome_of_single_x_on_surface_lights_two_or_one_z_checks():
    c = surface_rotated(5)
    for q in range(c.n):
        x = np.zeros(c.n, int); x[q] = 1
        s = c.syndrome(ex=x)
        lit_kinds = {c.kinds[i] for i in np.nonzero(s)[0]}
        assert lit_kinds == {"Z"} and 1 <= s.sum() <= 2


def test_geometry_shapes_match_generators():
    for c in CAT:
        g = c.geometry
        assert g.qubit_pos.shape == (c.n, 2) and g.check_pos.shape == (c.m, 2) and len(g.shape) == c.m
