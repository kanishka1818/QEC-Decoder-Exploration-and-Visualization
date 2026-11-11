"""
Tests for qdec/decoders.py.

What these tests are FOR: our from-scratch decoders are checked against the
established ones (PyMatching, ldpc) on the same shots, and against exact
answers on problems small enough to enumerate.
"""

import numpy as np
import pytest

from qdec.circuits import memory_circuit
from qdec.decoders import BPOSDDecoder, LookupDecoder, MWPMDecoder, sinter_decoders
from qdec.dem import DecodingProblem, MatchingGraph
from qdec.families import bb, five_qubit, steane, surface_rotated


def _cc(code, p):
    LX, LZ = code.logicals()
    return DecodingProblem.code_capacity(code.HZ, LZ[:, code.n:], p)


def test_lookup_fixes_every_single_error_steane():
    code = steane()
    pb = _cc(code, 0.01)
    dec = LookupDecoder(pb, max_weight=1)
    assert len(dec) == 8                                   # 7 single errors + no error, all distinct
    for q in range(code.n):
        e = np.zeros(code.n, np.uint8); e[q] = 1
        s = (pb.H.astype(int) @ e) % 2
        assert np.array_equal(dec.decode(s), (pb.L.astype(int) @ e) % 2)


@pytest.mark.parametrize("d,p,shots", [(3, 0.005, 1500), (5, 0.005, 400)])
def test_mwpm_agrees_with_pymatching(d, p, shots):
    pymatching = pytest.importorskip("pymatching")
    c = memory_circuit(surface_rotated(d), d, p)
    dem = c.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    ours = MWPMDecoder(MatchingGraph(dem))
    pm = pymatching.Matching.from_detector_error_model(dem)
    det, obs = c.compile_detector_sampler(seed=7).sample(shots, separate_observables=True)
    pred = np.array([ours.decode(s) & 1 for s in det])
    assert np.mean(pred == pm.decode_batch(det)[:, 0]) >= 0.995


def test_mwpm_exact_dp_and_ilp_agree():
    """The bitmask DP and the integer program must find equally good matchings."""
    c = memory_circuit(surface_rotated(5), 5, 0.006)
    dem = c.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    g = MatchingGraph(dem)
    a, b = MWPMDecoder(g, max_exact=30), MWPMDecoder(g, max_exact=0)    # b always uses the ILP
    det = c.compile_detector_sampler(seed=3).sample(150)
    for s in det:
        if 0 < s.sum() <= 14:
            a.decode(s); b.decode(s)
            cost = lambda d: sum(d._dijkstra(u, {v})[v][0] for u, v in d.last_matching)
            assert abs(cost(a) - cost(b)) < 1e-6


def test_bposd_reproduces_syndrome_and_beats_or_matches_ldpc():
    ldpc = pytest.importorskip("ldpc")
    code = bb("gross")
    p = 0.05
    pb = _cc(code, p)
    ours = BPOSDDecoder(pb, osd_order=6)
    lib = ldpc.BpOsdDecoder(code.HZ, error_rate=p, max_iter=50, bp_method="minimum_sum",
                            ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=6)
    rng = np.random.default_rng(11)
    fo = fl = 0
    for _ in range(600):
        e = (rng.random(code.n) < p).astype(np.uint8)
        s = (pb.H.astype(int) @ e) % 2
        corr = ours.decode_error(s)
        assert np.array_equal((pb.H.astype(int) @ corr) % 2, s)       # OSD always explains the syndrome
        true = (pb.L.astype(int) @ e) % 2
        fo += np.any((pb.L.astype(int) @ corr) % 2 != true)
        fl += np.any((pb.L.astype(int) @ lib.decode(s.astype(np.uint8))) % 2 != true)
    assert fo <= fl * 1.3 + 5                                          # within noise of ldpc, or better


def test_lookup_on_non_css_five_qubit_code():
    """Our lookup works on the full symplectic picture too: every single Pauli error is fixed."""
    code = five_qubit()
    n = code.n
    # faults = X, Z, Y on each qubit; detectors = the 4 generators; observables = the 2 logicals
    cols, obs = [], []
    LX, LZ = code.logicals()
    from qdec.core import symplectic_product
    for q in range(n):
        for ex, ez in [(1, 0), (0, 1), (1, 1)]:
            v = np.zeros(2 * n, np.uint8); v[q] = ex; v[n + q] = ez
            cols.append(code.syndrome(v[:n], v[n:]))
            obs.append(symplectic_product(np.vstack([LX, LZ]), v, n)[:, 0])
    pb = DecodingProblem(np.array(cols).T, np.array(obs).T, np.full(len(cols), 0.01), {})
    dec = LookupDecoder(pb, max_weight=1)
    assert len(dec) == 16                                             # perfect code: all 16 syndromes
    for j in range(len(cols)):
        assert np.array_equal(dec.decode(pb.H[:, j]), pb.L[:, j])


def test_sinter_adapters_run():
    sinter = pytest.importorskip("sinter")
    c = memory_circuit(surface_rotated(3), 3, 0.004)
    stats = sinter.collect(num_workers=1, tasks=[sinter.Task(circuit=c, json_metadata={"d": 3})],
                           decoders=["own_mwpm"], custom_decoders=sinter_decoders(), max_shots=300)
    assert stats and stats[0].shots >= 300
