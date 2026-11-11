"""
Tests for qdec/tanner.py.

What these tests are FOR: checking that the graph we draw is exactly the
check matrix (one edge per 1 in H, right node degrees), that the torus
geometry helpers measure distances the short way round, and that the
drawing functions run without error for every code.
"""

import matplotlib

matplotlib.use("Agg")                       # no screen needed

import matplotlib.pyplot as plt
import numpy as np
import pytest

from qdec.codes import bb72, gross_code, toric_code
from qdec.tanner import (central_check, draw_check_support, draw_tanner_torus, edge_lengths,
                         tanner_graph, torus_displacement)

CODES = [toric_code(4), bb72(), gross_code()]


@pytest.mark.parametrize("code", CODES, ids=lambda c: c.name)
def test_graph_matches_check_matrices(code):
    G = tanner_graph(code)
    assert G.number_of_nodes() == code.n + code.HX.shape[0] + code.HZ.shape[0]
    assert G.number_of_edges() == int(code.HX.sum() + code.HZ.sum())
    w = int(code.HX.sum(1).max())
    for c in range(code.HX.shape[0]):
        assert G.degree[("X", c)] == w and G.degree[("Z", c)] == w
    # each data qubit sits in w/2 X checks and w/2 Z checks (one per monomial
    # of A or of B), so its total degree is w
    for q in range(code.n):
        assert G.degree[("q", q)] == w


def test_displacement_takes_short_way_round():
    period = np.array([10, 10])
    assert list(torus_displacement((0, 0), (9, 0), period)) == [-1, 0]   # 1 step back, not 9 forward
    assert list(torus_displacement((0, 0), (3, 4), period)) == [3, 4]


def test_edge_length_counts():
    assert set(edge_lengths(toric_code(5), "X")) == {1}                 # toric: all neighbours
    lens = edge_lengths(gross_code(), "Z")
    assert (lens == 1).sum() == 4 * 72 and (lens > 1).sum() == 2 * 72   # 4 short + 2 long per check


@pytest.mark.parametrize("code", CODES, ids=lambda c: c.name)
def test_drawing_runs(code):
    fig, axs = plt.subplots(1, 2)
    draw_tanner_torus(code, "X", axs[0])
    draw_check_support(code, "Z", central_check(code, "Z"), axs[1])
    plt.close(fig)
