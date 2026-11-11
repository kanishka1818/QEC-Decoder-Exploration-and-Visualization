"""
bp_viz.py - pictures of belief propagation happening on the torus.

What this file is FOR
---------------------
It turns a BPResult (from bp.py) into pictures you can read:

    * draw_bp_state   one moment of BP: every data qubit coloured by how
                      likely BP thinks it is flipped, lit checks outlined,
                      the true error marked with a cross.
    * bp_frames       a row of those pictures, iteration by iteration, plus
                      the "unsatisfied checks vs iteration" curve.
    * qubit_traces    each chosen qubit's flip probability vs iteration -
                      the clearest way to see a qubit stuck near 50/50.
    * save_bp_gif     the same thing as an animation.

Colour scale: white = "surely fine" (P = 0), dark red = "surely flipped"
(P = 1).  A qubit hovering in the middle is a qubit BP cannot decide on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

from .bp import BPResult
from .codes import BBCode
from .tanner import _frame

CMAP = "Reds"


def draw_bp_state(code: BBCode, check_kind: str, probs, syndrome, true_error=None,
                  guess=None, ax=None, title: str = "", show_colorbar: bool = False):
    """
    One frame.  `check_kind` is the type of check that produced the syndrome
    ('Z' checks see X errors, 'X' checks see Z errors).
    """
    if ax is None:
        r, c = code.grid_shape()
        _, ax = plt.subplots(figsize=(c * 0.35 + 1, r * 0.35 + 1))
    dpos = code.data_positions()
    cpos = code.check_positions(check_kind)
    syndrome = np.asarray(syndrome).astype(bool)
    # checks: small grey squares, lit ones big with a thick border
    ax.scatter(cpos[~syndrome, 1], cpos[~syndrome, 0], s=12, marker="s", c="0.8", zorder=2)
    ax.scatter(cpos[syndrome, 1], cpos[syndrome, 0], s=120, marker="s", facecolors="gold",
               edgecolors="k", lw=1.5, zorder=4, label="lit check")
    # data qubits coloured by P(flipped)
    sc = ax.scatter(dpos[:, 1], dpos[:, 0], s=70, c=probs, cmap=CMAP, vmin=0, vmax=1,
                    edgecolors="0.4", lw=0.5, zorder=3)
    if guess is not None:                       # BP's current hard guess: blue ring
        g = np.asarray(guess).astype(bool)
        ax.scatter(dpos[g, 1], dpos[g, 0], s=170, facecolors="none", edgecolors="#1f77b4",
                   lw=1.8, zorder=5, label="BP guess")
    if true_error is not None:                  # what really happened: black cross
        e = np.asarray(true_error).astype(bool)
        ax.scatter(dpos[e, 1], dpos[e, 0], s=60, marker="x", c="k", lw=2, zorder=6, label="true error")
    _frame(code, ax)
    ax.set_title(title, fontsize=9)
    if show_colorbar:
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02, label="P(flipped)")
    return ax


def bp_frames(code: BBCode, check_kind: str, result: BPResult, syndrome, true_error=None,
              iterations=None, suptitle: str = ""):
    """
    A row of frames at the chosen iterations (default: 0, 1, 2, 3, 5, last)
    plus the unsatisfied-checks curve on the right.
    """
    last = result.iterations
    if iterations is None:
        iterations = sorted({0, 1, 2, 3, 5, last} & set(range(last + 1)))
    r, c = code.grid_shape()
    w = c * 0.28 + 0.6
    fig, axs = plt.subplots(1, len(iterations) + 1,
                            figsize=(w * len(iterations) + 3.5, r * 0.28 + 1.4),
                            gridspec_kw={"width_ratios": [1] * len(iterations) + [1.3]})
    probs = result.prob_history
    for ax, it in zip(axs, iterations):
        draw_bp_state(code, check_kind, probs[it], syndrome, true_error,
                      guess=result.guess_history[it], ax=ax,
                      title=f"iteration {it}\n{result.unsatisfied_history[it]} checks unhappy")
    # unsatisfied-check curve
    ax = axs[-1]
    ax.plot(result.unsatisfied_history, "o-", ms=3)
    ax.set_xlabel("iteration")
    ax.set_ylabel("unsatisfied checks")
    ax.set_title("converged" if result.converged else "did NOT converge", fontsize=9,
                 color="green" if result.converged else "red")
    ax.grid(alpha=0.3)
    ax.set_box_aspect(1)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


def qubit_traces(result: BPResult, qubits, ax=None, labels=None, title: str = ""):
    """P(flipped) of a few qubits against iteration number."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 3.2))
    probs = result.prob_history
    for k, q in enumerate(qubits):
        ax.plot(probs[:, q], "o-", ms=3, label=labels[k] if labels else f"qubit {q}")
    ax.axhline(0.5, ls="--", c="0.5", lw=1)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("iteration")
    ax.set_ylabel("P(flipped)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    return ax


def save_bp_gif(code: BBCode, check_kind: str, result: BPResult, syndrome, true_error,
                path: str, fps: int = 2):
    """Animated GIF: one frame per BP iteration (uses matplotlib's Pillow writer)."""
    r, c = code.grid_shape()
    fig, ax = plt.subplots(figsize=(c * 0.35 + 1.5, r * 0.35 + 1.2))
    probs = result.prob_history

    def frame(it):
        ax.clear()
        draw_bp_state(code, check_kind, probs[it], syndrome, true_error,
                      guess=result.guess_history[it], ax=ax,
                      title=f"{code.name} - iteration {it}: "
                            f"{result.unsatisfied_history[it]} checks unhappy")

    anim = animation.FuncAnimation(fig, frame, frames=result.iterations + 1)
    anim.save(path, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return path
