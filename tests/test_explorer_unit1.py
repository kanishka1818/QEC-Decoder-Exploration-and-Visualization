"""
Tests for explorer unit 1 (8 Nov 2026): Stim-based distances, lighter logical
operators for drawing, logical "strings", and the 'How it's written' forms.

What these tests are FOR: the explorer now shows distances computed by Stim,
draws X-bar / Z-bar as lines and prints each code three ways. Each of those is
checked here against something independent (brute force, the definition of a
logical operator, or the code's own matrices).
"""

import numpy as np
import pytest

from qdec.circuits import code_capacity_circuit, stim_distance
from qdec.core import symplectic_product
from qdec.export import code_to_dict, distance_info, logical_path, reduce_logical, written_forms
from qdec.families import (bb, catalogue, classical_repetition_H, hamming_H, hypergraph_product, repetition,
                           shor, steane, surface_rotated, surface_unrotated, toric)


# ------------------------------------------------------------- distances ----
@pytest.mark.parametrize("d", [3, 5, 9, 15])
def test_stim_distance_surface_exact_at_any_size(d):
    dd, how = stim_distance(surface_rotated(d))
    assert dd == d and how.startswith("exact")


@pytest.mark.parametrize("code", [shor(), steane(), surface_unrotated(3), surface_unrotated(5), toric(4)],
                         ids=lambda c: c.name)
def test_stim_distance_agrees_with_brute_force(code):
    assert stim_distance(code)[0] == code.distance(max_weight=6)


def test_stim_distance_regression_hgp_hamming_is_3_not_1():
    """Stim's shortest_graphlike_error reported d = 1 here (a detector-less decomposition piece);
    stim_distance must use the search for codes whose errors light 3+ checks."""
    d, how = stim_distance(hypergraph_product(hamming_H(), hamming_H()))
    assert d == 3 and how.startswith("upper bound")


def test_stim_distance_bb72_matches_paper():
    d, how = stim_distance(bb("bb72"))
    assert d == 6 and how.startswith("upper bound")


def test_repetition_distance_text():
    info = distance_info(repetition(5))
    assert info["d"] == 1 and "5 against bit flips" in info["how"]


def test_code_capacity_circuit_is_silent_without_errors():
    c = code_capacity_circuit(surface_rotated(5), "Z", p=0.0)
    assert not c.compile_detector_sampler().sample(32, append_observables=True).any()


# --------------------------------------------------------------- logicals ----
@pytest.mark.parametrize("code", [surface_rotated(5), surface_unrotated(5), toric(4), shor(), bb("bb72")],
                         ids=lambda c: c.name)
def test_reduced_logicals_are_still_the_same_logicals(code):
    n = code.n
    LX, LZ = code.logicals()
    S = np.hstack([code.Sx, code.Sz])
    for L in np.vstack([LX, LZ]):
        R = reduce_logical(L, S, n)
        assert int((R[:n] | R[n:]).sum()) <= int((L[:n] | L[n:]).sum())     # never heavier
        assert not symplectic_product(code.S, R, n).any()                  # still passes every check
        # same logical: it anticommutes with exactly the same logicals as before
        allL = np.vstack([LX, LZ])
        assert np.array_equal(symplectic_product(allL, R, n), symplectic_product(allL, L, n))


@pytest.mark.parametrize("code", [surface_rotated(5), surface_unrotated(5), toric(4)], ids=lambda c: c.name)
def test_surface_logicals_draw_as_minimum_weight_strings(code):
    d = code_to_dict(code, dist={"d": None, "how": ""})
    for L in d["logicals"]:
        for t in ("x", "z"):
            assert len(L[t]) == code.info["d_known"]          # reduced all the way down to weight d
            assert L[t + "_path"]["string"]                  # and drawn as a connected line


def test_gross_logical_is_not_a_string():
    p = logical_path(bb("gross"), [3, 9, 12, 15, 21, 24, 39, 45, 48, 51, 57, 60])
    assert p["string"] is False


# ---------------------------------------------------------- written forms ----
@pytest.mark.parametrize("code", catalogue(), ids=lambda c: c.name)
def test_written_forms_exist_and_match_the_matrices(code):
    forms = written_forms(code)
    assert [f["title"] for f in forms][:2] == ["Stabilizer generators", "Check matrices (0/1)"]
    assert len(forms) == 3 and forms[2]["title"] in ("Polynomial form", "Algebraic form")
    # the stabilizer list has one line per generator
    assert len(forms[0]["pre"].splitlines()) == code.m
    # the 0/1 text really is the matrix
    pre = forms[1]["pre"]
    if code.is_css and code.HX.size:
        rows = [ln for ln in pre.splitlines() if ln and set(ln) <= {"0", "1"}]
        assert rows[: code.HX.shape[0]] == ["".join(map(str, r)) for r in code.HX]


def test_bb_polynomial_form_names_the_polynomials():
    pre = written_forms(bb("gross"))[2]["pre"]
    assert "A = x^3 + y + y^2" in pre and "B = y^3 + x + x^2" in pre and "l = 12, m = 6" in pre


def test_five_qubit_polynomial_form():
    from qdec.families import five_qubit
    pre = written_forms(five_qubit())[2]["pre"]
    assert "a(x) = 1 + x³" in pre and "b(x) = x + x²" in pre
