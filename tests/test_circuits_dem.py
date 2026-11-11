"""
Tests for qdec/circuits.py and qdec/dem.py.

What these tests are FOR: the memory circuits we write must be *correct
experiments* (silent when there is no noise, right numbers of detectors and
observables), keep the code's distance when the gate order is hook-safe, and
our DEM reading must give a graph for surface codes and a hypergraph for BB codes.
"""

import numpy as np
import pytest

from qdec.circuits import circuit_distance, circuit_stats, cnot_layers, memory_circuit
from qdec.dem import DecodingProblem, MatchingGraph, fault_degree_histogram, parse_dem
from qdec.families import bb, repetition, shor, steane, surface_rotated, surface_unrotated, toric

CSS = [repetition(5), shor(), steane(), surface_rotated(3), surface_rotated(5), surface_unrotated(3),
       toric(4), bb("bb72"), bb("gross")]


@pytest.mark.parametrize("code", CSS, ids=lambda c: c.name)
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_noiseless_circuit_is_silent(code, basis):
    if code.info["family"] == "repetition" and basis == "X":
        pytest.skip("repetition code has no X checks / X logical readout")
    c = memory_circuit(code, rounds=3, p=0.0, basis=basis)
    s = c.compile_detector_sampler().sample(64, append_observables=True)
    assert not s.any()


@pytest.mark.parametrize("code", CSS, ids=lambda c: c.name)
def test_counts(code):
    r = 3
    c = memory_circuit(code, rounds=r, p=0.0, basis="Z")
    mx, mz = code.HX.shape[0], code.HZ.shape[0]
    assert c.num_observables == code.k
    assert c.num_detectors == mz * r + mx * (r - 1) + mz          # Z every round, X from round 2, final Z
    assert circuit_stats(c)["CX"] == r * int(code.HX.sum() + code.HZ.sum())


@pytest.mark.parametrize("code", [surface_rotated(3), toric(4), bb("gross")], ids=lambda c: c.name)
def test_layers_never_reuse_a_qubit(code):
    for kind in "XZ":
        for layer in cnot_layers(code, kind):
            qs = [q for _, q in layer]
            rs = [r for r, _ in layer]
            assert len(set(qs)) == len(qs) and len(set(rs)) == len(rs)


def test_gross_layers_are_one_per_polynomial_term():
    """Translation invariance: each of the 6 terms of A, B is exactly one full CNOT layer of 72 gates."""
    g = bb("gross")
    for kind in "XZ":
        L = cnot_layers(g, kind)
        assert len(L) == 6 and all(len(layer) == 72 for layer in L)


@pytest.mark.parametrize("d", [3, 5])
def test_hook_safe_order_keeps_distance_plain_order_loses_it(d):
    code = surface_rotated(d)
    for basis in "ZX":
        assert circuit_distance(memory_circuit(code, 3, 0.001, basis)) == d
    assert circuit_distance(memory_circuit(code, 3, 0.001, "X", order="plain")) < d


def test_surface_dem_is_graphlike_after_decomposition():
    c = memory_circuit(surface_rotated(5), 3, 0.001)
    dem = c.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    assert max(fault_degree_histogram(dem, use_parts=True)) <= 2
    assert MatchingGraph(dem).dropped == 0


def test_gross_dem_is_a_hypergraph():
    c = memory_circuit(bb("gross"), 2, 0.001)
    dem = c.detector_error_model(decompose_errors=False, approximate_disjoint_errors=True)
    hist = fault_degree_histogram(dem)
    assert max(hist) >= 3 and sum(v for k, v in hist.items() if k >= 3) > sum(v for k, v in hist.items() if k <= 2)
    with pytest.raises(ValueError):
        c.detector_error_model(decompose_errors=True)


def test_decoding_problem_matches_dem():
    c = memory_circuit(surface_rotated(3), 2, 0.002)
    dem = c.detector_error_model(approximate_disjoint_errors=True)
    pb = DecodingProblem.from_dem(dem)
    faults, nd, no, _ = parse_dem(dem)
    assert pb.H.shape == (nd, len(faults)) and pb.L.shape == (no, len(faults))
    assert np.allclose(pb.priors, [f.p for f in faults])
