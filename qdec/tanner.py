"""
tanner.py - the Tanner graph of a code, and pictures of it on the torus.

What this file is FOR
---------------------
A decoder never sees "a code"; it sees a graph:
    * one node per data qubit (the unknowns: "was this qubit flipped?"),
    * one node per check (the clues: "did this parity come out odd?"),
    * an edge whenever a check looks at a qubit (a 1 in H).
That is the Tanner graph (Tanner, IEEE Trans. Inf. Theory 27, 533 (1981)).
Belief propagation literally runs on this graph: messages travel along
these edges.  So to "see" decoding, we first have to see this graph.

Because BB codes live on a torus, we draw the graph on the l x m grid of
cells (layout in codes.py).  An edge whose two ends are grid neighbours
is "short"; anything longer is "long-range".  The toric code has only
short edges; the gross code has some long ones - that is exactly what
lets it pack 12 logical qubits into 144 data qubits, and also why it
needs long-range wiring in hardware (Bravyi et al. 2024, Fig. 1).

All drawing is plain matplotlib + networkx (networkx only to hand back a
standard graph object for anyone who wants to analyse it).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .codes import BBCode

# colours used across all pictures (kept in one place so every figure matches)
COLORS = {
    "L": "#4C72B0",       # L data qubits
    "R": "#55A868",       # R data qubits
    "X": "#C44E52",       # X checks
    "Z": "#8172B3",       # Z checks
    "short": "#B0B0B0",   # short edges
    "long": "#DD8452",    # long-range edges
}


# --------------------------------------------------------------------------
# graph object
# --------------------------------------------------------------------------
def tanner_graph(code: BBCode, kinds: str = "XZ") -> nx.Graph:
    """
    Tanner graph as a networkx Graph.

    Node names: ('q', i) for data qubit i, ('X', c) / ('Z', c) for checks.
    Every node carries 'pos' = its (row, col) grid position on the torus.
    """
    G = nx.Graph()
    dpos = code.data_positions()
    for q in range(code.n):
        G.add_node(("q", q), pos=tuple(dpos[q]), block="L" if q < code.n // 2 else "R")
    for kind in kinds:
        cpos = code.check_positions(kind)
        H = code.H(kind)
        for c in range(H.shape[0]):
            G.add_node((kind, c), pos=tuple(cpos[c]))
            for q in np.nonzero(H[c])[0]:
                G.add_edge((kind, c), ("q", int(q)))
    return G


# --------------------------------------------------------------------------
# torus geometry
# --------------------------------------------------------------------------
def torus_period(code: BBCode) -> np.ndarray:
    """Size of the picture grid: the torus wraps after this many (row, col) steps."""
    return np.array(code.grid_shape())


def torus_displacement(p_from, p_to, period) -> np.ndarray:
    """Shortest (row, col) step from p_from to p_to when the grid wraps around."""
    d = (np.asarray(p_to) - np.asarray(p_from)) % period
    return np.where(d > period / 2, d - period, d)      # pick the shorter way round


def edge_lengths(code: BBCode, kind: str) -> np.ndarray:
    """
    Torus distance (|row step| + |col step|) of every check-qubit edge for one
    check type.  1 means "grid neighbour"; anything bigger is long-range.
    """
    period = torus_period(code)
    dpos, cpos, H = code.data_positions(), code.check_positions(kind), code.H(kind)
    out = []
    for c in range(H.shape[0]):
        for q in np.nonzero(H[c])[0]:
            out.append(np.abs(torus_displacement(cpos[c], dpos[q], period)).sum())
    return np.array(out)


# --------------------------------------------------------------------------
# drawing
# --------------------------------------------------------------------------
def draw_nodes(code: BBCode, ax, kinds: str = "XZ", size: float = 60):
    """Scatter all data qubits and checks at their torus grid positions."""
    dpos = code.data_positions()
    half = code.n // 2
    # columns -> x axis, rows -> y axis (row 0 at the top, like a matrix)
    ax.scatter(dpos[:half, 1], dpos[:half, 0], s=size, marker="o", c=COLORS["L"], label="L qubit", zorder=3)
    ax.scatter(dpos[half:, 1], dpos[half:, 0], s=size, marker="o", c=COLORS["R"], label="R qubit", zorder=3)
    for kind in kinds:
        cp = code.check_positions(kind)
        ax.scatter(cp[:, 1], cp[:, 0], s=size, marker="s", c=COLORS[kind], label=f"{kind} check", zorder=3)
    _frame(code, ax)


def _frame(code: BBCode, ax):
    """Common axis styling: torus box, equal aspect, row 0 on top."""
    period = torus_period(code)
    ax.set_xlim(-1, period[1])
    ax.set_ylim(period[0], -1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    # dashed box = one copy of the torus; leaving it on one side = entering on the other
    ax.add_patch(plt.Rectangle((-0.5, -0.5), period[1], period[0], fill=False, ls="--", lw=0.8, ec="0.5"))


def central_check(code: BBCode, kind: str) -> int:
    """Index of the check closest to the middle of the picture (nice for single-check plots)."""
    centre = np.array(code.grid_shape()) / 2
    return int(np.argmin(np.abs(code.check_positions(kind) - centre).sum(1)))


def draw_edges(code: BBCode, ax, kind: str, checks=None, lw: float = 0.8, alpha: float = 0.6,
               mark_ghosts: bool = False):
    """
    Draw the edges of the given checks (default: all checks of this kind).
    Each edge is drawn from the check towards the NEAREST copy of its qubit
    on the torus, so wrap-around edges leave the box instead of crossing it.
    """
    period = torus_period(code)
    dpos, cpos, H = code.data_positions(), code.check_positions(kind), code.H(kind)
    checks = range(H.shape[0]) if checks is None else checks
    for c in checks:
        for q in np.nonzero(H[c])[0]:
            d = torus_displacement(cpos[c], dpos[q], period)
            end = cpos[c] + d
            style = "short" if np.abs(d).sum() == 1 else "long"
            if mark_ghosts and np.any(end != dpos[q]):
                # the shortest route wraps round the torus; for a single-check
                # picture it is clearer to draw a dashed line to the copy of the
                # qubit inside the box (same qubit, the line just "goes round")
                end, ls = dpos[q], "--"
            else:
                ls = "-"
            ax.plot([cpos[c][1], end[1]], [cpos[c][0], end[0]], c=COLORS[style], ls=ls,
                    lw=lw if style == "short" else lw * 1.5, alpha=alpha, zorder=1)


def draw_tanner_torus(code: BBCode, kind: str = "X", ax=None, title: str | None = None):
    """Whole Tanner graph for one check type, drawn on the torus."""
    if ax is None:
        r, c = code.grid_shape()
        _, ax = plt.subplots(figsize=(c * 0.45 + 2, r * 0.45 + 1))
    draw_edges(code, ax, kind, lw=0.6, alpha=0.35)
    draw_nodes(code, ax, kinds=kind, size=25)
    lens = edge_lengths(code, kind)
    ax.set_title(title or f"{code.name}: {kind}-check Tanner graph\n"
                 f"{(lens == 1).sum()} short edges, {(lens > 1).sum()} long-range edges", fontsize=10)
    return ax


def draw_check_support(code: BBCode, kind: str, check: int, ax=None, title: str | None = None):
    """
    Highlight ONE check and the data qubits it looks at.  This is the clearest
    way to see what a polynomial means: each monomial of A or B is one arrow.
    """
    if ax is None:
        r, c = code.grid_shape()
        _, ax = plt.subplots(figsize=(c * 0.45 + 2, r * 0.45 + 1))
    draw_nodes(code, ax, kinds=kind, size=18)
    draw_edges(code, ax, kind, checks=[check], lw=2.0, alpha=0.9, mark_ghosts=True)
    cp = code.check_positions(kind)[check]
    ax.scatter([cp[1]], [cp[0]], s=180, marker="s", facecolors="none", edgecolors="k", lw=2, zorder=4)
    qs = np.nonzero(code.H(kind)[check])[0]
    dp = code.data_positions()[qs]
    ax.scatter(dp[:, 1], dp[:, 0], s=150, facecolors="none", edgecolors="k", lw=2, zorder=4)
    ax.set_title(title or f"{code.name}\n{kind} check #{check} -> qubits {list(map(int, qs))}\n"
                 "grey = neighbour edge, orange = long-range (dashed = wraps round the torus)", fontsize=8)
    return ax


def draw_qubit_set(code: BBCode, qubits, ax=None, title: str = "", color: str = "k"):
    """
    Highlight a set of data qubits (e.g. the support of a logical operator)
    on the torus, with all other qubits drawn faintly.
    """
    if ax is None:
        r, c = code.grid_shape()
        _, ax = plt.subplots(figsize=(c * 0.35 + 1, r * 0.35 + 1))
    dpos = code.data_positions()
    ax.scatter(dpos[:, 1], dpos[:, 0], s=15, c="0.8", zorder=2)
    sel = dpos[np.asarray(qubits, dtype=int)]
    ax.scatter(sel[:, 1], sel[:, 0], s=90, c=color, zorder=3)
    _frame(code, ax)
    ax.set_title(title, fontsize=9)
    return ax
