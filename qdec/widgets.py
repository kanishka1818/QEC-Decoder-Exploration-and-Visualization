"""
widgets.py - the same explorer as the browser page, inside Jupyter (ipywidgets).

What this file is FOR
---------------------
Sometimes you want to poke at a code without leaving the notebook.  These
helpers put dropdowns and sliders on top of the matplotlib views:

    code_explorer()   pick a code, a view (lattice / Tanner / matrix / logicals),
                      a stabilizer to highlight, and type an error like "X3 Z7 Y12";
                      the syndrome and (for small codes) the lookup decoder's verdict
                      are printed under the picture.
    bp_player(...)    step through a recorded BP run with a slider.
    mwpm_player(...)  type X-error qubits on a surface/toric code and watch matching.

The browser explorer (explorer/index.html) and these widgets read the same
data, so they always agree.  Widgets only work in a live kernel; in a saved
notebook the last picture stays as a static image.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np

from .core import symplectic_product
from .export import lookup_table
from .families import catalogue
from .views import draw_check_matrix, draw_lattice, draw_logicals, draw_matching, draw_tanner


def parse_error(text: str, n: int):
    """'X3 Z7 Y12' -> (ex, ez) 0/1 vectors.  Ignores anything it cannot read."""
    ex, ez = np.zeros(n, np.uint8), np.zeros(n, np.uint8)
    for letter, q in re.findall(r"([XYZxyz])\s*(\d+)", text or ""):
        q = int(q)
        if q < n:
            L = letter.upper()
            ex[q] ^= L in "XY"
            ez[q] ^= L in "ZY"
    return ex, ez


def describe_error(code, ex, ez, table=None) -> str:
    """Text summary: lit checks, and the lookup decoder's verdict when a table exists."""
    s = code.syndrome(ex, ez)
    lines = [f"{int(s.sum())} of {code.m} checks light up: {list(np.nonzero(s)[0])}"]
    if table is not None:
        corr = table.get("".join(map(str, s)))
        if corr is None:
            lines.append("lookup decoder: syndrome not in the table (too many errors)")
        else:
            cx, cz = np.zeros(code.n, np.uint8), np.zeros(code.n, np.uint8)
            for q, L in corr.items():
                cx[q] = L in "XY"; cz[q] = L in "ZY"
            res = np.concatenate([(ex + cx) % 2, (ez + cz) % 2])
            LX, LZ = code.logicals()
            bad = symplectic_product(np.vstack([LX, LZ]), res, code.n).any()
            lines.append(f"lookup decoder applies {corr or 'nothing'} -> {'LOGICAL ERROR' if bad else 'corrected'}")
    return "\n".join(lines)


def code_explorer(codes=None):
    """Interactive code explorer (needs ipywidgets and a live kernel)."""
    import ipywidgets as W
    from IPython.display import display

    codes = codes or catalogue()
    names = [c.name for c in codes]
    tables = {}
    dd = W.Dropdown(options=names, value=names[6] if len(names) > 6 else names[0], description="code")
    view = W.ToggleButtons(options=["lattice", "Tanner", "matrix", "logicals"], value="lattice")
    gen = W.IntSlider(value=0, min=0, max=codes[0].m - 1, description="stabilizer")
    err = W.Text(value="", placeholder="e.g. X3 Z7 Y12", description="error")
    out = W.Output()

    def code_of(name):
        return codes[names.index(name)]

    def redraw(*_):
        c = code_of(dd.value)
        gen.max = c.m - 1
        ex, ez = parse_error(err.value, c.n)
        if c.n <= 13 and c.name not in tables:
            tables[c.name] = lookup_table(c, 3 if c.n <= 9 else 2)
        with out:
            out.clear_output(wait=True)
            fig, ax = plt.subplots(figsize=(7, 6))
            if view.value == "lattice":
                draw_lattice(c, ax, highlight=gen.value, error_x=ex, error_z=ez)
            elif view.value == "Tanner":
                draw_tanner(c, ax)
            elif view.value == "matrix":
                draw_check_matrix(c, ax, highlight_row=gen.value)
            else:
                draw_logicals(c, ax)
            plt.show()
            print(f"{c.name} {c.params()} - {c.info.get('family', '')}")
            print(f"stabilizer {gen.value}: {c.generator_string(gen.value) if c.n <= 60 else c.kinds[gen.value] + '-type'}")
            if ex.any() or ez.any():
                print(describe_error(c, ex, ez, tables.get(c.name)))

    for w in (dd, view, gen, err):
        w.observe(redraw, names="value")
    display(W.VBox([W.HBox([dd, view]), W.HBox([gen, err]), out]))
    redraw()


def bp_player(replay: dict, codes=None):
    """Slider over the iterations of a recorded BP replay (from export.default_replays())."""
    import ipywidgets as W
    from IPython.display import display

    codes = codes or catalogue()
    code = next(c for c in codes if c.name == replay["code"])
    it = W.IntSlider(value=0, min=0, max=len(replay["probs"]) - 1, description="iteration")
    play = W.Play(value=0, min=0, max=len(replay["probs"]) - 1, interval=400)
    W.jslink((play, "value"), (it, "value"))
    out = W.Output()
    geo = code.geometry

    def redraw(*_):
        with out:
            out.clear_output(wait=True)
            fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.3, 1]})
            draw_lattice(code, ax, show="Z", labels=False, title=f"{replay['title']} - iteration {it.value}")
            ax.scatter(geo.qubit_pos[:, 0], geo.qubit_pos[:, 1], c=replay["probs"][it.value], cmap="Reds",
                       vmin=0, vmax=1, s=120, edgecolors="0.4", zorder=8)
            for q in replay["error"]:
                ax.scatter(*geo.qubit_pos[q], marker="x", s=80, c="k", lw=2, zorder=9)
            ax2.plot(replay["unsat"], "o-", ms=3)
            ax2.axvline(it.value, color="orange")
            ax2.set_xlabel("iteration"); ax2.set_ylabel("unsatisfied checks"); ax2.grid(alpha=0.3)
            plt.show()

    it.observe(redraw, names="value")
    display(W.VBox([W.HBox([play, it]), out]))
    redraw()


def mwpm_player(code=None):
    """Type X-error qubits (e.g. '6 12 13') on a surface/toric code and watch MWPM decode them."""
    import ipywidgets as W
    from IPython.display import display
    from .families import surface_rotated

    code = code or surface_rotated(5)
    txt = W.Text(value="6 12 13", description="X errors on")
    out = W.Output()

    def redraw(*_):
        e = np.zeros(code.n, np.uint8)
        for q in re.findall(r"\d+", txt.value):
            if int(q) < code.n:
                e[int(q)] ^= 1
        with out:
            out.clear_output(wait=True)
            fig, ax = plt.subplots(figsize=(6, 6))
            draw_matching(code, e, ax)
            plt.show()

    txt.observe(redraw, names="value")
    display(W.VBox([txt, out]))
    redraw()
