"""
views.py - the same code seen four ways (matplotlib).

What this file is FOR
---------------------
Recognising a code by looking at it.  For ANY StabilizerCode:

    draw_lattice        the geometric picture: stabilizers as coloured faces
                        (local codes) or as stars with spokes (long-range codes)
    draw_tanner         the Tanner graph drawn on the same geometry
    draw_tanner_layers  the textbook Tanner graph: qubits in one row, checks in another
    draw_check_matrix   the symplectic matrix [S_x | S_z] (red = X, blue = Z, purple = Y)
    draw_logicals       a logical operator drawn on the lattice
    code_card           all of the above in one figure, with the code's description

Colour convention used everywhere (also in the browser explorer):
    X stabilizer = red, Z stabilizer = blue, mixed (non-CSS) = purple.
"""

from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from .core import StabilizerCode

COL = {"X": "#d1495b", "Z": "#2a6fdb", "mixed": "#8a5cc2", "qubit": "#333333",
       "error": "#f2a900", "lit": "#ffd23f", "logX": "#d1495b", "logZ": "#2a6fdb"}


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def _support_points(code: StabilizerCode, i: int) -> np.ndarray:
    """Positions of generator i's qubits, unwrapped around its check position (torus)."""
    g = code.geometry
    c = g.check_pos[i]
    return np.array([g.nearest_image(c, g.qubit_pos[q]) for q in code.support(i)])


def _ordered_polygon(points: np.ndarray, centre: np.ndarray) -> np.ndarray:
    """Sort points by angle around `centre` so they form a simple polygon."""
    ang = np.arctan2(points[:, 1] - centre[1], points[:, 0] - centre[0])
    return points[np.argsort(ang)]


def _face_points(code: StabilizerCode, i: int) -> np.ndarray:
    """
    Corner points of generator i drawn as a face.
    Weight 2: if the check position is off the line between the two qubits
    (rotated-surface boundary), use a triangle through the check position
    (the familiar 'half-moon' boundary plaquette); otherwise a thick line.
    """
    pts = _support_points(code, i)
    c = code.geometry.check_pos[i]
    if len(pts) == 2:
        mid = pts.mean(0)
        if np.linalg.norm(mid - c) > 1e-6:
            return np.array([pts[0], c + (c - mid) * 0.3, pts[1]])
        return pts
    centre = pts.mean(0)
    return _ordered_polygon(pts, centre)


def _setup(ax, code: StabilizerCode, pad: float = 1.0):
    g = code.geometry
    allp = np.vstack([g.qubit_pos, g.check_pos])
    if g.period is not None:
        ax.set_xlim(-pad, g.period[0] + pad - 1)
        ax.set_ylim(g.period[1] + pad - 1, -pad)
        ax.add_patch(plt.Rectangle((-0.5, -0.5), g.period[0], g.period[1], fill=False,
                                   ls="--", lw=0.8, ec="0.6"))
    else:
        lo, hi = allp.min(0) - pad, allp.max(0) + pad
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(hi[1], lo[1])          # y grows downward
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _label(code):
    """'name [[n,k,d]]' without repeating the parameters when the name already has them."""
    return code.name if code.params() in code.name else f"{code.name}  {code.params()}"


def _qubit_size(code):
    return 260 if code.n <= 13 else 110 if code.n <= 60 else 40


# --------------------------------------------------------------------------
# lattice view
# --------------------------------------------------------------------------
def draw_lattice(code: StabilizerCode, ax=None, show: str = "both", highlight=None,
                 error_x=None, error_z=None, labels: bool | None = None, title: str | None = None,
                 alpha: float = 0.35):
    """
    Geometric picture of the code.

    show       : 'both', 'X' or 'Z' - which generators to draw (mixed ones always shown)
    highlight  : generator index (or list) to emphasise, with its Pauli letters written on its qubits
    error_x/z  : 0/1 vectors; the error is marked on the qubits and the lit generators glow yellow
    labels     : write qubit numbers (default: only for small codes)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    g = code.geometry
    kinds = code.kinds
    labels = code.n <= 41 if labels is None else labels
    hl = [] if highlight is None else ([highlight] if np.isscalar(highlight) else list(highlight))
    lit = np.zeros(code.m, bool)
    if error_x is not None or error_z is not None:
        lit = code.syndrome(error_x, error_z).astype(bool)
    # colour-code faces (Steane): an X and a Z generator on the SAME qubits.  Overlapping red
    # and blue would blend to purple and look "mixed", so the Z copy is drawn as a dashed outline.
    x_supports = {tuple(code.support(i)) for i in range(code.m) if kinds[i] == "X"}
    twin = [kinds[i] == "Z" and tuple(code.support(i)) in x_supports and show == "both" for i in range(code.m)]

    for i in range(code.m):
        kind = kinds[i]
        if show != "both" and kind != "mixed" and kind != show:
            continue
        colour = COL[kind]
        emph = i in hl
        a = 0.85 if emph else alpha
        if g.shape[i] == "polygon":
            pts = _face_points(code, i)
            if len(pts) == 2:
                ax.plot(pts[:, 0], pts[:, 1], color=colour, lw=10 if code.n < 30 else 6,
                        alpha=a, solid_capstyle="round", zorder=1)
            elif twin[i] and not emph:
                ax.add_patch(Polygon(pts * 0.9 + pts.mean(0) * 0.1, closed=True, fill=False, ec=colour,
                                     lw=2.2, ls="--", zorder=2))
            else:
                ax.add_patch(Polygon(pts, closed=True, fc=colour, ec=colour if not emph else "k",
                                     lw=2.5 if emph else 0.8, alpha=a, zorder=1))
        else:   # star: spokes from the check node
            c = g.check_pos[i]
            faint = 0.25 if code.n <= 60 else 0.07          # big codes: keep spokes quiet
            for p in _support_points(code, i):
                ax.plot([c[0], p[0]], [c[1], p[1]], color=colour, lw=2.2 if emph else 0.6,
                        alpha=0.9 if emph else faint, zorder=1)
            ax.scatter([c[0]], [c[1]], marker="s", s=60 if emph else 14, color=colour,
                       edgecolors="k" if emph else "none", zorder=3)
        if lit[i]:
            cpos = g.check_pos[i]
            ax.scatter([cpos[0]], [cpos[1]], marker="*", s=260, color=COL["lit"],
                       edgecolors="k", lw=0.8, zorder=6)

    s = _qubit_size(code)
    ax.scatter(g.qubit_pos[:, 0], g.qubit_pos[:, 1], s=s, color="white", edgecolors=COL["qubit"],
               lw=1.2, zorder=4)
    # error markers
    for vec, letter in [(error_x, "X"), (error_z, "Z")]:
        if vec is not None:
            for q in np.nonzero(np.asarray(vec))[0]:
                ax.scatter(*g.qubit_pos[q], s=s, color=COL["error"], edgecolors="k", zorder=5)
    if labels:
        for q, (x, y) in enumerate(g.qubit_pos):
            ax.text(x, y, str(q), ha="center", va="center", fontsize=7 if code.n > 13 else 8, zorder=7)
    for i in hl:
        word = code.generator_string(i)
        for q in code.support(i):
            x, y = g.qubit_pos[q]
            ax.text(x + 0.35, y - 0.35, word[q], fontsize=10, fontweight="bold",
                    color=COL["mixed" if word[q] == "Y" else word[q]], zorder=8)
    _setup(ax, code)
    ax.set_title(title if title is not None else _label(code), fontsize=10)
    return ax


# --------------------------------------------------------------------------
# Tanner graph views
# --------------------------------------------------------------------------
def draw_tanner(code: StabilizerCode, ax=None, show: str = "both", title: str | None = None):
    """Tanner graph on the code's own geometry: qubits (circles) - checks (squares)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    g = code.geometry
    for i, kind in enumerate(code.kinds):
        if show != "both" and kind != "mixed" and kind != show:
            continue
        c = g.check_pos[i]
        for p in _support_points(code, i):
            ax.plot([c[0], p[0]], [c[1], p[1]], color=COL[kind], lw=0.8, alpha=0.6, zorder=1)
        ax.scatter([c[0]], [c[1]], marker="s", s=_qubit_size(code) * 0.5, color=COL[kind], zorder=3)
    ax.scatter(g.qubit_pos[:, 0], g.qubit_pos[:, 1], s=_qubit_size(code) * 0.6, color="white",
               edgecolors=COL["qubit"], zorder=4)
    _setup(ax, code)
    ax.set_title(title if title is not None else f"{code.name}: Tanner graph", fontsize=10)
    return ax


def draw_tanner_layers(code: StabilizerCode, ax=None, title: str | None = None):
    """Textbook layout: checks on top, qubits at the bottom, one edge per non-identity entry."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 3))
    n, m = code.n, code.m
    qx = np.linspace(0, 1, n)
    cx = np.linspace(0, 1, m)
    for i, kind in enumerate(code.kinds):
        for q in code.support(i):
            ax.plot([cx[i], qx[q]], [1, 0], color=COL[kind], lw=0.6, alpha=0.5)
    ax.scatter(cx, np.ones(m), marker="s", c=[COL[k] for k in code.kinds], s=30 if m < 60 else 8, zorder=3)
    ax.scatter(qx, np.zeros(n), c="white", edgecolors="k", s=30 if n < 60 else 8, zorder=3)
    ax.set_ylim(-0.2, 1.2); ax.axis("off")
    ax.text(-0.02, 1, "checks", ha="right", va="center", fontsize=8)
    ax.text(-0.02, 0, "qubits", ha="right", va="center", fontsize=8)
    ax.set_title(title if title is not None else f"{code.name}: Tanner graph (layered)", fontsize=10)
    return ax


# --------------------------------------------------------------------------
# check matrix
# --------------------------------------------------------------------------
def draw_check_matrix(code: StabilizerCode, ax=None, highlight_row=None, title: str | None = None):
    """
    The symplectic matrix: one row per generator, one column per qubit.
    Cell colour: red = X on that qubit, blue = Z, purple = Y, white = identity.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    img = np.ones((code.m, code.n, 3))
    rgb = lambda h: [int(h[k:k + 2], 16) / 255 for k in (1, 3, 5)]
    img[(code.Sx == 1) & (code.Sz == 0)] = rgb(COL["X"])
    img[(code.Sx == 0) & (code.Sz == 1)] = rgb(COL["Z"])
    img[(code.Sx == 1) & (code.Sz == 1)] = rgb(COL["mixed"])
    ax.imshow(img, aspect="auto", interpolation="nearest")
    if highlight_row is not None:
        ax.add_patch(plt.Rectangle((-0.5, highlight_row - 0.5), code.n, 1, fill=False, ec="k", lw=2))
    ax.set_xlabel("qubit"); ax.set_ylabel("generator")
    if code.n <= 25:
        ax.set_xticks(range(code.n))
        ax.tick_params(labelsize=6)
    ax.set_title(title if title is not None else f"{code.name}: stabilizer matrix (red X, blue Z, purple Y)",
                 fontsize=9)
    return ax


# --------------------------------------------------------------------------
# logicals
# --------------------------------------------------------------------------
def draw_logicals(code: StabilizerCode, ax=None, pair: int = 0, title: str | None = None):
    """Draw logical pair `pair`: X-bar qubits in red, Z-bar qubits in blue (shared qubits half/half)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    LX, LZ = code.logicals()
    n = code.n
    draw_lattice(code, ax, alpha=0.12, labels=False, title="")
    g = code.geometry
    s = _qubit_size(code)
    for L, key, off in [(LX[pair], "logX", -0.18), (LZ[pair], "logZ", 0.18)]:
        x, z = L[:n], L[n:]
        word = "".join("IXZY"[int(a) + 2 * int(b)] for a, b in zip(x, z))
        for q in np.nonzero(x | z)[0]:
            px, py = g.qubit_pos[q]
            ax.scatter([px + off], [py], s=s * 0.7, color=COL[key], edgecolors="k", zorder=6)
            if code.n <= 41:
                ax.text(px + off, py, word[q], ha="center", va="center", color="white",
                        fontsize=7, fontweight="bold", zorder=7)
    wx, wz = int((LX[pair][:n] | LX[pair][n:]).sum()), int((LZ[pair][:n] | LZ[pair][n:]).sum())
    ax.set_title(title if title is not None else
                 f"{code.name}: logical pair {pair}  (X-bar weight {wx} red, Z-bar weight {wz} blue)",
                 fontsize=9)
    return ax


# --------------------------------------------------------------------------
# everything on one card
# --------------------------------------------------------------------------
def code_card(code: StabilizerCode, highlight: int | None = None):
    """Lattice + Tanner + matrix + logicals + description, one figure per code."""
    if highlight is None:          # pick a generator near the middle of the picture
        g = code.geometry
        centre = g.qubit_pos.mean(0)
        highlight = int(np.argmin(np.linalg.norm(g.check_pos - centre, axis=1)))
    fig = plt.figure(figsize=(16, 5.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.1, 1])
    draw_lattice(code, fig.add_subplot(gs[0]), highlight=highlight,
                 title=f"lattice: generator {highlight} = {code.generator_string(highlight) if code.n <= 13 else code.kinds[highlight] + '-type, weight ' + str(len(code.support(highlight)))}")
    draw_tanner(code, fig.add_subplot(gs[1]), title="Tanner graph on the lattice")
    draw_check_matrix(code, fig.add_subplot(gs[2]), highlight_row=highlight, title="stabilizer matrix")
    draw_logicals(code, fig.add_subplot(gs[3]), title="logical pair 0 (red X-bar, blue Z-bar)")
    desc = textwrap.fill(code.info.get("description", ""), 230)
    fig.suptitle(f"{_label(code)}   -   {code.info.get('family', '')}\n{desc}",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# decoding pictures on the lattice
# --------------------------------------------------------------------------
def draw_matching(code: StabilizerCode, error_x, ax=None, p: float = 0.05, title: str | None = None):
    """
    Code-capacity MWPM on X errors, drawn on the lattice: lit Z checks (stars),
    matched paths (green dashes through the qubits used), true error (orange),
    correction (green rings).  Returns (ax, logical_error: bool).
    """
    from .decoders import MWPMDecoder
    from .dem import MatchingGraph
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    LX, LZ = code.logicals()
    n = code.n
    g = MatchingGraph.from_check_matrix(code.HZ, LZ[:, n:], p)
    dec = MWPMDecoder(g)
    e = np.asarray(error_x, np.uint8)
    s = (code.HZ.astype(int) @ e) % 2
    dec.decode(s)
    draw_lattice(code, ax, show="Z", error_x=e, labels=code.n <= 41, title="")
    zidx = code.index_of("Z")
    geo = code.geometry
    corr = np.zeros(n, np.uint8)
    for (a, b), path in dec.last_paths.items():
        qs = g.path_faults(path)
        pts = [geo.check_pos[zidx[a]]]
        for q in qs:
            corr[q] ^= 1
            pts.append(geo.nearest_image(pts[-1], geo.qubit_pos[q]))
        if b != g.boundary:
            pts.append(geo.nearest_image(pts[-1], geo.check_pos[zidx[b]]))
        pts = np.array(pts)
        ax.plot(pts[:, 0], pts[:, 1], "--", color="#2f8a4f", lw=3, zorder=6)
    for q in np.nonzero(corr)[0]:
        ax.scatter(*geo.qubit_pos[q], s=_qubit_size(code) * 2.2, facecolors="none", edgecolors="#2f8a4f",
                   lw=2.5, zorder=7)
    fail = bool(((LZ[:, n:].astype(int) @ ((corr + e) % 2)) % 2).any())
    ax.set_title(title if title is not None else
                 f"MWPM on {code.name}: {'LOGICAL ERROR' if fail else 'corrected'}", fontsize=10,
                 color="#c4424f" if fail else "#2f8a4f")
    return ax, fail


def draw_osd_steps(code: StabilizerCode, error_x, decoder, axs=None):
    """
    Three panels for one shot that BP could not finish:
      1. BP's average belief per qubit (what OSD sorts by)
      2. the pivot qubits OSD keeps (most reliable set that can explain any syndrome)
      3. OSD's answer vs the true error
    `decoder` is a decoders.BPOSDDecoder built on this code's H_Z.
    """
    if axs is None:
        _, axs = plt.subplots(1, 3, figsize=(15, 6))
    e = np.asarray(error_x, np.uint8)
    s = (code.HZ.astype(int) @ e) % 2
    corr = decoder.decode_error(s)
    info = decoder.last
    bp = info["bp"]
    avg = bp.llr_history[1:].mean(0) if bp.iterations else bp.llr_history[0]
    prob = 1 / (1 + np.exp(np.clip(avg, -30, 30)))
    geo = code.geometry
    lit = np.zeros(code.m, bool)
    lit[np.array(code.index_of("Z"))[np.nonzero(s)[0]]] = True

    def base(ax, title):
        _setup(ax, code)
        ax.scatter(geo.check_pos[lit, 0], geo.check_pos[lit, 1], marker="*", s=200, color=COL["lit"],
                   edgecolors="k", lw=0.6, zorder=6)
        ax.set_title(title, fontsize=10)

    ax = axs[0]
    sc = ax.scatter(geo.qubit_pos[:, 0], geo.qubit_pos[:, 1], c=prob, cmap="Reds", vmin=0, vmax=1,
                    s=_qubit_size(code), edgecolors="0.5", zorder=4)
    base(ax, f"1. BP's average belief (BP {'converged' if bp.converged else 'did not converge'})")
    plt.colorbar(sc, ax=ax, fraction=0.046, label="P(flipped), averaged over iterations")
    ax = axs[1]
    ax.scatter(geo.qubit_pos[:, 0], geo.qubit_pos[:, 1], s=_qubit_size(code), color="white",
               edgecolors="0.6", zorder=3)
    if info.get("osd"):
        piv = np.asarray(info["order"])[info["pivots"]]
        ax.scatter(geo.qubit_pos[piv, 0], geo.qubit_pos[piv, 1], s=_qubit_size(code), color="#2a62c9",
                   alpha=0.6, zorder=4)
        base(ax, f"2. OSD keeps {len(piv)} most reliable independent qubits (blue)")
    else:
        base(ax, "2. (BP converged - OSD not needed)")
    ax = axs[2]
    ax.scatter(geo.qubit_pos[:, 0], geo.qubit_pos[:, 1], s=_qubit_size(code), color="white",
               edgecolors="0.6", zorder=3)
    eq = np.nonzero(e)[0]; cq = np.nonzero(corr)[0]
    ax.scatter(geo.qubit_pos[eq, 0], geo.qubit_pos[eq, 1], marker="x", s=90, color="k", lw=2.5, zorder=6)
    ax.scatter(geo.qubit_pos[cq, 0], geo.qubit_pos[cq, 1], s=_qubit_size(code) * 2.2, facecolors="none",
               edgecolors="#2f8a4f", lw=2.5, zorder=7)
    LX, LZ = code.logicals()
    fail = bool(((LZ[:, code.n:].astype(int) @ ((corr + e) % 2)) % 2).any())
    base(ax, f"3. answer (green rings) vs true error (x): {'LOGICAL ERROR' if fail else 'corrected'}")
    return axs, fail
