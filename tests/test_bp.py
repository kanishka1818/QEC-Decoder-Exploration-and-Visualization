"""
Tests for qdec/bp.py and qdec/bp_viz.py.

What these tests are FOR: proving our hand-written min-sum BP is a correct
BP (same answers and same beliefs as the independent `ldpc` package), that
its recorded history is self-consistent, that it shows the two textbook
failure behaviours we use in the notebook, and that every picture function runs.
"""

import itertools

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from qdec.bp import MinSumBP, logical_failure
from qdec.bp_viz import bp_frames, qubit_traces, save_bp_gif
from qdec.codes import bb72, gross_code, toric_code


def syndrome_of(H, e):
    return (H.astype(int) @ e) % 2


def test_zero_syndrome_needs_no_iterations():
    code = toric_code(4)
    r = MinSumBP(code.HZ, 0.05).decode(np.zeros(code.HZ.shape[0], dtype=int))
    assert r.converged and r.iterations == 0 and not r.correction.any()


@pytest.mark.parametrize("code", [toric_code(5), bb72(), gross_code()], ids=lambda c: c.name)
def test_every_single_qubit_error_is_corrected(code):
    bp = MinSumBP(code.HZ, 0.03, max_iter=30)
    for q in range(code.n):
        e = np.zeros(code.n, dtype=int); e[q] = 1
        r = bp.decode(syndrome_of(code.HZ, e))
        assert r.converged and np.array_equal(r.correction, e)


def test_history_is_consistent():
    code = gross_code()
    rng = np.random.default_rng(0)
    bp = MinSumBP(code.HZ, 0.04, max_iter=20)
    for _ in range(30):
        e = (rng.random(code.n) < 0.04).astype(int)
        s = syndrome_of(code.HZ, e)
        r = bp.decode(s)
        T = r.iterations + 1
        assert r.llr_history.shape == (T, code.n) and r.guess_history.shape == (T, code.n)
        assert len(r.unsatisfied_history) == T
        assert np.all((r.prob_history >= 0) & (r.prob_history <= 1))
        assert np.array_equal(r.guess_history, (r.llr_history < 0).astype(np.uint8))
        assert r.converged == (r.unsatisfied_history[-1] == 0)
        if r.converged:                                        # converged = syndrome reproduced
            assert np.array_equal(syndrome_of(code.HZ, r.correction), s)


@pytest.mark.parametrize("code,p", [(toric_code(6), 0.05), (bb72(), 0.04), (gross_code(), 0.04)],
                         ids=["toric", "bb72", "gross"])
def test_matches_ldpc_min_sum(code, p):
    """Same corrections, same iteration counts and same final beliefs as ldpc's BpDecoder."""
    ldpc = pytest.importorskip("ldpc")
    lib = ldpc.BpDecoder(code.HZ, error_rate=p, max_iter=30, bp_method="minimum_sum",
                         ms_scaling_factor=0.625, schedule="parallel")
    ours = MinSumBP(code.HZ, p, max_iter=30, alpha=0.625)
    rng = np.random.default_rng(1)
    for _ in range(150):
        e = (rng.random(code.n) < p).astype(np.uint8)
        s = syndrome_of(code.HZ, e)
        r = ours.decode(s)
        d = lib.decode(s.astype(np.uint8))
        assert np.array_equal(r.correction, d)
        assert r.converged == lib.converge
        if s.any():        # for an all-zero syndrome ldpc returns early and leaves old beliefs behind
            assert r.iterations == lib.iter
            assert np.allclose(r.llr_history[-1], np.array(lib.log_prob_ratios), atol=1e-6)


def test_toric_tie_is_never_broken():
    """Diagonal pair of lit checks: two equally short corrections -> BP stays symmetric and never converges."""
    code = toric_code(6)
    e = np.zeros(code.n, dtype=int); e[[15, 56]] = 1
    r = MinSumBP(code.HZ, 0.05, max_iter=40).decode(syndrome_of(code.HZ, e))
    assert not r.converged
    four = r.prob_history[-1][[15, 56, 21, 57]]
    assert np.allclose(four, four[0]) and four[0] < 0.5


def test_all_toric_ties_fail():
    code = toric_code(4)
    H = code.HZ.astype(int)
    seen = {}
    for a, b in itertools.combinations(range(code.n), 2):
        e = np.zeros(code.n, dtype=int); e[[a, b]] = 1
        s = tuple(syndrome_of(H, e))
        if sum(s) == 2:
            seen[s] = seen.get(s, 0) + 1
    tied = [s for s, c in seen.items() if c >= 2]
    bp = MinSumBP(H, 0.05, max_iter=40)
    assert tied and all(not bp.decode(np.array(s)).converged for s in tied)


def test_logical_failure_detects_logicals_not_stabilizers():
    code = toric_code(4)
    LX, LZ = code.logicals()
    zero = np.zeros(code.n, dtype=int)
    stab = code.HX[0].astype(int)                 # an X stabilizer: harmless residual
    assert not logical_failure(stab, zero, LZ)
    assert logical_failure(LX[0].astype(int), zero, LZ)      # an X logical: failure


def test_pictures_run(tmp_path):
    code = toric_code(4)
    e = np.zeros(code.n, dtype=int); e[[0, 5]] = 1
    s = syndrome_of(code.HZ, e)
    r = MinSumBP(code.HZ, 0.05, max_iter=5).decode(s)
    plt.close(bp_frames(code, "Z", r, s, e))
    qubit_traces(r, [0, 5]); plt.close("all")
    path = save_bp_gif(code, "Z", r, s, e, str(tmp_path / "bp.gif"))
    assert (tmp_path / "bp.gif").stat().st_size > 0 and path.endswith(".gif")
