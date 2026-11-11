"""
Tests for qdec/views.py, qdec/export.py and qdec/widgets.py.

What these tests are FOR: every picture function must run on every code
(no code is "special-cased" out of the pictures), the JSON the browser
explorer reads must be complete and consistent with the Python objects, and
the small lookup tables shipped to the browser must give correct answers.
"""

import json

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from qdec.core import symplectic_product
from qdec.decoders import BPOSDDecoder
from qdec.dem import DecodingProblem
from qdec.export import build_atlas_data, build_explorer_html, code_to_dict, default_replays, lookup_table, toric_tie
from qdec.families import bb, catalogue, five_qubit, steane, surface_rotated, toric
from qdec.views import (code_card, draw_check_matrix, draw_lattice, draw_logicals, draw_matching,
                        draw_osd_steps, draw_tanner, draw_tanner_layers)
from qdec.widgets import describe_error, parse_error

CAT = catalogue()


@pytest.mark.parametrize("code", CAT, ids=lambda c: c.name)
def test_every_view_runs_on_every_code(code):
    fig, axs = plt.subplots(1, 5)
    draw_lattice(code, axs[0], highlight=0)
    draw_tanner(code, axs[1])
    draw_tanner_layers(code, axs[2])
    draw_check_matrix(code, axs[3], highlight_row=0)
    draw_logicals(code, axs[4])
    plt.close(fig)


def test_code_card_and_decoding_pictures():
    plt.close(code_card(steane()))
    s = surface_rotated(5)
    e = np.zeros(s.n, np.uint8); e[[6, 12, 13]] = 1
    _, fail = draw_matching(s, e)
    assert fail is False
    e = np.zeros(s.n, np.uint8); e[[0, 5, 10]] = 1
    _, fail = draw_matching(s, e)
    assert fail is True                      # 3 errors down a column beat a d=5 code
    g = bb("gross")
    LX, LZ = g.logicals()
    dec = BPOSDDecoder(DecodingProblem.code_capacity(g.HZ, LZ[:, g.n:], 0.03), osd_order=4)
    e = np.zeros(g.n, np.uint8); e[[18, 79, 94, 117]] = 1
    draw_osd_steps(g, e, dec)
    plt.close("all")


def test_export_is_json_and_consistent():
    data = build_atlas_data(use_stim=False)
    text = json.dumps(data)
    assert len(text) < 5_000_000
    names = {c["name"] for c in data["codes"]}
    for r in data["replays"]:
        assert r["code"] in names                       # every replay points at a code in the list
    for c, obj in zip(data["codes"], CAT):
        assert c["n"] == obj.n and c["k"] == obj.k and len(c["gens"]) == obj.m
        for i, g in enumerate(c["gens"]):
            assert g["support"] == list(map(int, obj.support(i)))


@pytest.mark.parametrize("code", [steane(), five_qubit(), surface_rotated(3)], ids=lambda c: c.name)
def test_lookup_table_corrects_every_single_qubit_error(code):
    table = lookup_table(code, 2)
    LX, LZ = code.logicals()
    L = np.vstack([LX, LZ])
    for q in range(code.n):
        for ex, ez in [(1, 0), (0, 1), (1, 1)]:
            x = np.zeros(code.n, np.uint8); z = np.zeros(code.n, np.uint8); x[q] = ex; z[q] = ez
            corr = table["".join(map(str, code.syndrome(x, z)))]
            for cq, letter in corr.items():
                x[cq] ^= letter in "XY"; z[cq] ^= letter in "ZY"
            assert not code.syndrome(x, z).any()                              # residual is silent
            assert not symplectic_product(L, np.concatenate([x, z]), code.n).any()   # and not logical


def test_replays_tell_the_intended_stories():
    r = {x["title"]: x for x in default_replays()}
    assert r["MWPM fixes a weight-3 chain (surface d=5)"]["logical_error"] is False
    assert r["MWPM loses to 3 errors down a column (surface d=5)"]["logical_error"] is True
    assert r["Three errors that are really one (surface d=5)"]["logical_error"] is False
    assert r["BP stuck on a tie (toric 4x4)"]["converged"] is False
    assert r["BP converges (gross code)"]["converged"] is True
    assert r["BP oscillates (gross code)"]["converged"] is False


def test_toric_tie_has_two_equal_explanations():
    t = toric(4)
    a, b = toric_tie(t)
    e = np.zeros(t.n, int); e[[a, b]] = 1
    s = (t.HZ.astype(int) @ e) % 2
    assert s.sum() == 2


def test_explorer_html_builds(tmp_path):
    out = build_explorer_html("explorer/template.html", str(tmp_path / "x.html"), use_stim=False)
    html = open(out).read()
    assert "/*__ATLAS_DATA__*/" not in html and "Stabilizer Atlas" in html and html.startswith("<!doctype html>")


def test_widget_text_helpers():
    st = steane()
    ex, ez = parse_error("X3 z5 Y2", st.n)
    assert list(np.nonzero(ex)[0]) == [2, 3] and list(np.nonzero(ez)[0]) == [2, 5]
    assert "corrected" in describe_error(st, *parse_error("X3", st.n), lookup_table(st))
    assert "LOGICAL ERROR" in describe_error(st, *parse_error("X0 X1", st.n), lookup_table(st))
