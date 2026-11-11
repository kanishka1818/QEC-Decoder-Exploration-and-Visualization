"""
export.py - write every code (and some recorded decoding runs) as plain JSON
for the browser explorer and the Jupyter widgets.

What this file is FOR
---------------------
The maths lives in Python; the clicking lives in the browser.  This file is
the bridge: one JSON object per code with positions, generators, logicals and
(for small codes) a full lookup-decoder table, plus "replays" - recorded BP
and MWPM runs that the explorer can step through iteration by iteration.

    code_to_dict(code)        -> dict (JSON-ready)
    lookup_table(code)        -> syndrome -> most likely Pauli correction (small codes)
    bp_replay(...), mwpm_replay(...)  -> recorded decoding runs
    build_atlas_data()        -> {"codes": [...], "replays": [...]} used by explorer/
"""

from __future__ import annotations

import itertools
import json
import os

import numpy as np

from .bp import MinSumBP
from .core import StabilizerCode, symplectic_product
from .decoders import MWPMDecoder
from .dem import MatchingGraph
from .families import bb, catalogue, surface_rotated, toric


def _r(a, nd=3):
    return [round(float(v), nd) for v in a]


def lookup_table(code: StabilizerCode, max_weight: int = 3) -> dict[str, dict]:
    """
    For every reachable syndrome (as a '0101..' string), the lightest Pauli
    error producing it (any of X, Y, Z on up to max_weight qubits) - the
    lookup decoder of Lesson 1, for the code-capacity depolarising picture.
    """
    n = code.n
    table: dict[str, dict] = {"0" * code.m: {}}
    for w in range(1, max_weight + 1):
        for supp in itertools.combinations(range(n), w):
            for letters in itertools.product("XZY", repeat=w):
                ex = np.zeros(n, np.uint8); ez = np.zeros(n, np.uint8)
                for q, L in zip(supp, letters):
                    ex[q] = L in "XY"; ez[q] = L in "ZY"
                key = "".join(map(str, code.syndrome(ex, ez)))
                if key not in table:
                    table[key] = {int(q): L for q, L in zip(supp, letters)}
    return table


# --------------------------------------------------------------------------
# distance, nicer logicals, and "how the code is written"
# --------------------------------------------------------------------------
def distance_info(code: StabilizerCode, use_stim: bool = True) -> dict:
    """
    {'d': int, 'how': str} for display.  Order of preference:
      small codes (n <= 41): exact brute force (and Stim, which agrees);
      CSS codes: Stim (exact for graph-like codes at any size, upper bound for qLDPC);
      otherwise the published value.
    The [[288,12,18]] search takes minutes, so it keeps the paper's value.
    """
    from .circuits import stim_distance
    if code.info.get("family") == "repetition":
        dX, _ = stim_distance_xz(code)
        return {"d": 1, "how": f"1 as a quantum code (no check sees a Z flip); {dX} against bit flips (Stim, exact)"}
    if code.n <= 41 and code.info.get("family") != "repetition":
        d = code.distance(max_weight=8)
        if d is not None:
            return {"d": int(d), "how": "exact (brute force over all small errors)"}
    if use_stim and code.is_css and code.n <= 200:
        d, how = stim_distance(code)
        extra = ""
        if code.info.get("d_known") is not None and how.startswith("upper"):
            extra = f"; matches the published d = {code.info['d_known']}"
        return {"d": int(d), "how": f"Stim, {how}{extra}"}
    return {"d": code.info.get("d_known"), "how": "published value (search too slow here)"}


def stim_distance_xz(code: StabilizerCode):
    """(d against X errors, d against Z errors) for a CSS code, with Stim; None where no checks exist."""
    from .circuits import code_capacity_circuit
    out = []
    for basis, checks in (("Z", code.HZ), ("X", code.HX)):
        if checks.shape[0] == 0:
            out.append(None)
            continue
        c = code_capacity_circuit(code, basis)
        out.append(len(c.shortest_graphlike_error(ignore_ungraphlike_errors=True)))
    return tuple(out)


def reduce_logical(L: np.ndarray, stabs: np.ndarray, n: int) -> np.ndarray:
    """
    Make a logical operator lighter by multiplying in stabilizers while that
    lowers its weight (greedy).  The result is the SAME logical (it differs by
    stabilizers), so pairings and commutation are unchanged - it just draws
    as a clean straight string on surface codes.
    """
    L = L.copy()
    weight = lambda v: int((v[:n] | v[n:]).sum())
    improved = True
    while improved:
        improved = False
        for s in stabs:
            cand = (L ^ s).astype(np.uint8)
            if weight(cand) < weight(L):
                L, improved = cand, True
    return L


def logical_path(code: StabilizerCode, support) -> dict:
    """
    Order a logical's qubits into a chain (nearest neighbour walk) so the
    explorer can draw it as a line.  'string' = True when consecutive qubits
    are close (surface / toric logicals), False for scattered supports (qLDPC),
    which are drawn as dots instead.
    """
    g = code.geometry
    pts = [np.asarray(g.qubit_pos[q], float) for q in support]
    if len(pts) <= 1:
        return {"order": [int(q) for q in support], "string": False}
    # start from the qubit farthest from the centroid, then walk to nearest unvisited
    centre = np.mean(pts, axis=0)
    left = list(range(len(pts)))
    cur = max(left, key=lambda i: np.linalg.norm(pts[i] - centre))
    order, gaps = [cur], []
    left.remove(cur)
    while left:
        dist = [np.linalg.norm(g.nearest_image(pts[cur], pts[j]) - pts[cur]) for j in left]
        j = left[int(np.argmin(dist))]
        gaps.append(min(dist))
        order.append(j)
        left.remove(j)
        cur = j
    step = 2.9 if g.period is None else 2.9
    return {"order": [int(support[i]) for i in order], "string": bool(max(gaps) <= step)}


def written_forms(code: StabilizerCode) -> list[dict]:
    """
    Sections for the explorer's 'How it's written' tab: the same code as
    (1) a list of stabilizers, (2) 0/1 check matrices, (3) its algebraic /
    polynomial recipe.  Each section: {'title', 'text', 'pre' (optional)}.
    """
    n = code.n
    out = []
    lines = []
    for i in range(code.m):
        word = code.generator_string(i)
        prod = " ".join(f"{word[q]}{q}" for q in code.support(i))
        lines.append(f"S{i:<3d} ({code.kinds[i]:5s})  {prod}")
    out.append({"title": "Stabilizer generators",
                "text": "Each line is one check: the Pauli it measures on each qubit (qubit numbers as on the lattice).",
                "pre": "\n".join(lines)})
    as01 = lambda M: "\n".join("".join(str(int(v)) for v in row) for row in M)
    if code.is_css:
        pre = (f"H_X  ({code.HX.shape[0]} x {n})\n" + as01(code.HX) +
               f"\n\nH_Z  ({code.HZ.shape[0]} x {n})\n" + as01(code.HZ)) if code.HX.size else \
              (f"H_Z  ({code.HZ.shape[0]} x {n})\n" + as01(code.HZ))
        text = ("CSS code: X checks and Z checks are separate 0/1 matrices. Row = check, column = qubit, "
                "1 = the check looks at that qubit. They commute exactly when H_X · H_Zᵀ = 0 (mod 2).")
    else:
        pre = f"[ S_x | S_z ]  ({code.m} x {2 * n})\n" + "\n".join(
            "".join(map(str, a)) + " | " + "".join(map(str, b)) for a, b in zip(code.Sx, code.Sz))
        text = ("Not CSS, so each check is a pair of 0/1 rows: its X part and its Z part (both 1 = Y). "
                "Two checks commute when x₁·z₂ + z₁·x₂ = 0 (mod 2).")
    out.append({"title": "Check matrices (0/1)", "text": text, "pre": pre})
    out.append(_algebraic_form(code))
    return out


def _algebraic_form(code: StabilizerCode) -> dict:
    """Family-specific recipe: polynomials for BB/toric/cyclic codes, products for HGP, etc."""
    fam = code.info.get("family", "")
    name = code.name
    if "bivariate" in fam or "toric" in fam:
        l, m = code.info.get("l"), code.info.get("m")
        A, B = code.info.get("A"), code.info.get("B")
        if "toric" in fam:
            L = int(round(np.sqrt(code.n / 2)))
            l = m = L
            A, B = "1 + x", "1 + y"
        terms_note = ("Each term xᵃyᵇ is one arrow: a check touches the qubit a cells away along x "
                      "and b cells away along y (which way round is a labelling choice). The same arrows repeat in every cell, which is why one term is "
                      "one full CNOT layer in the circuit.")
        pre = (f"torus: l = {l}, m = {m}   ({l*m} cells, each with an L qubit, an R qubit, an X check, a Z check)\n"
               f"x = cyclic shift by one cell along the l direction  ({l*m} x {l*m} permutation matrix)\n"
               f"y = cyclic shift by one cell along the m direction\n\n"
               f"A = {A}\nB = {B}\n\n"
               f"H_X = [ A | B ]        (X check at cell c: L qubits c·A, R qubits c·B)\n"
               f"H_Z = [ Bᵀ | Aᵀ ]\n\n"
               f"Why they commute: H_X H_Zᵀ = A B + B A = 0 (mod 2), because x and y commute.")
        return {"title": "Polynomial form", "text": terms_note, "pre": pre}
    if fam == "repetition":
        d = code.n
        return {"title": "Polynomial form",
                "text": "A cyclic-style code: every check is the same pattern (1 1) shifted one qubit along the line.",
                "pre": f"check polynomial h(x) = 1 + x\nrow i of H_Z = xⁱ·(1 + x)   for i = 0 … {d-2}\n"
                       f"(no X checks: a single Z flip goes unseen, so d_Z = 1, while d_X = {d})"}
    if "perfect" in fam:
        return {"title": "Polynomial form",
                "text": "A cyclic code over 5 qubits. Write the X part and Z part of the first check as "
                        "polynomials in x (x = shift by one qubit, x⁵ = 1); the other checks are x, x², x³ times it.",
                "pre": "first check XZZXI:\n  X part a(x) = 1 + x³\n  Z part b(x) = x + x²\n"
                       "check k = ( xᵏ a(x) | xᵏ b(x) )  mod (x⁵ − 1),  k = 0 … 3\n"
                       "(the 5th shift is the product of the other four, so only 4 are independent)"}
    if "colour" in fam:
        return {"title": "Algebraic form",
                "text": "Steane = CSS(Hamming, Hamming): both H_X and H_Z are the [7,4,3] Hamming checks. "
                        "Hamming is a cyclic code with generator g(x) = 1 + x + x³ when the qubits are "
                        "numbered around a cycle (this drawing uses the triangle numbering instead).",
                "pre": "Hamming checks: column j = binary of the face set qubit j belongs to\n"
                       "every non-zero 3-bit pattern appears exactly once -> every single error has its own syndrome"}
    if "small CSS" in fam:
        return {"title": "Algebraic form",
                "text": "Shor's code is a concatenation: an outer 3-qubit phase-flip code whose 'qubits' are "
                        "3-qubit bit-flip blocks. No polynomial; its structure is a product.",
                "pre": "|0_L> = (|000> + |111>)⊗3,  |1_L> = (|000> − |111>)⊗3\n"
                       "inner (per row): Z0Z1, Z1Z2   outer (across rows): X on rows 0+1, X on rows 1+2"}
    if fam == "surface":
        if "unrotated" in name:
            return {"title": "Algebraic form",
                    "text": "The unrotated surface code is the hypergraph product of two repetition codes "
                            "(we check this: same n, k, d).",
                    "pre": "H_rep = repetition checks (d−1) x d\n"
                           "H_X = [ H_rep ⊗ I | I ⊗ H_repᵀ ]\nH_Z = [ I ⊗ H_rep | H_repᵀ ⊗ I ]"}
        return {"title": "Algebraic form",
                "text": "No polynomial: the rotated code is read off a chessboard rule.",
                "pre": "data qubits at (2i+1, 2j+1), plaquettes at (2i, 2j)\n"
                       "plaquette type = X if i + j is even, else Z\n"
                       "keep all inner plaquettes; on the edge keep only X ones (top/bottom) and Z ones (left/right)"}
    if "hypergraph" in fam:
        return {"title": "Algebraic form",
                "text": "Hypergraph product of two classical codes H1 (r1 x n1) and H2 (r2 x n2).",
                "pre": "H_X = [ H1 ⊗ I_n2 | I_r1 ⊗ H2ᵀ ]\nH_Z = [ I_n1 ⊗ H2 | H1ᵀ ⊗ I_r2 ]\n"
                       "n = n1·n2 + r1·r2   (bit×bit qubits, then check×check qubits)"}
    return {"title": "Algebraic form", "text": "No special recipe recorded for this code.", "pre": ""}


def code_to_dict(code: StabilizerCode, with_lookup: bool | None = None, dist: dict | None = None) -> dict:
    g = code.geometry
    LX, LZ = code.logicals()
    n = code.n
    # lighter, nicer-looking representatives for drawing (same logicals up to stabilizers)
    SX = np.hstack([code.Sx, code.Sz])
    LX = np.array([reduce_logical(L, SX, n) for L in LX])
    LZ = np.array([reduce_logical(L, SX, n) for L in LZ])

    def pauli_list(row):
        x, z = row[:n], row[n:]
        return [[int(q), "IXZY"[int(x[q]) + 2 * int(z[q])]] for q in np.nonzero(x | z)[0]]

    gens = []
    for i in range(code.m):
        word = code.generator_string(i)
        gens.append({"kind": code.kinds[i], "support": [int(q) for q in code.support(i)],
                     "letters": "".join(word[q] for q in code.support(i)), "shape": g.shape[i]})
    d = {
        "id": code.name, "name": code.name, "params": code.params(), "n": n, "k": code.k, "m": code.m,
        "d": code.info.get("d_known"), "family": code.info.get("family", ""),
        "description": code.info.get("description", ""), "references": code.info.get("references", []),
        "is_css": code.is_css, "qubits": [_r(p, 3) for p in g.qubit_pos],
        "checks": [_r(p, 3) for p in g.check_pos], "period": list(g.period) if g.period else None,
        "gens": gens,
        "logicals": [{"x": pauli_list(a), "z": pauli_list(b),
                      "x_path": logical_path(code, [int(q) for q in np.nonzero(a[:n] | a[n:])[0]]),
                      "z_path": logical_path(code, [int(q) for q in np.nonzero(b[:n] | b[n:])[0]])}
                     for a, b in zip(LX, LZ)],
        "written": written_forms(code),
    }
    dist = dist or distance_info(code)
    d["d"], d["d_how"] = dist["d"], dist["how"]
    d["params"] = f"[[{n},{code.k},{dist['d']}]]"
    for key in ("A", "B", "l", "m"):
        if key in code.info:
            d["poly_" + key if key in "AB" else key + "_torus"] = code.info[key]
    if with_lookup is None:
        with_lookup = n <= 13
    if with_lookup:
        d["lookup"] = lookup_table(code, max_weight=3 if n <= 9 else 2)
    return d


# --------------------------------------------------------------------------
# recorded decoding runs
# --------------------------------------------------------------------------
def bp_replay(code: StabilizerCode, error_qubits, p: float, title: str, caption: str,
              max_iter: int = 30) -> dict:
    """BP on X errors (seen by Z checks): probabilities, guesses and unsatisfied counts per iteration."""
    e = np.zeros(code.n, np.uint8); e[list(error_qubits)] = 1
    s = (code.HZ.astype(int) @ e) % 2
    res = MinSumBP(code.HZ, p, max_iter=max_iter).decode(s)
    zidx = code.index_of("Z")
    return {"code": code.name, "kind": "bp", "title": title, "caption": caption,
            "error": [int(q) for q in error_qubits],
            "lit": [int(zidx[r]) for r in np.nonzero(s)[0]],
            "probs": [_r(row, 3) for row in res.prob_history],
            "guess": [[int(q) for q in np.nonzero(gu)[0]] for gu in res.guess_history],
            "unsat": [int(u) for u in res.unsatisfied_history], "converged": bool(res.converged)}


def mwpm_replay(code: StabilizerCode, error_qubits, title: str, caption: str, p: float = 0.05) -> dict:
    """Code-capacity MWPM on X errors: lit checks, the chosen pairs and their qubit paths."""
    LX, LZ = code.logicals()
    g = MatchingGraph.from_check_matrix(code.HZ, LZ[:, code.n:], p)
    dec = MWPMDecoder(g)
    e = np.zeros(code.n, np.uint8); e[list(error_qubits)] = 1
    s = (code.HZ.astype(int) @ e) % 2
    pred = dec.decode(s)
    zidx = code.index_of("Z")
    pairs, corr = [], np.zeros(code.n, np.uint8)
    for (a, b), path in dec.last_paths.items():
        qs = g.path_faults(path)
        for q in qs:
            corr[q] ^= 1
        pairs.append({"a": int(zidx[a]), "b": None if b == g.boundary else int(zidx[b]),
                      "qubits": [int(q) for q in qs]})
    residual = (corr + e) % 2
    fail = bool(((LZ[:, code.n:].astype(int) @ residual) % 2).any())
    return {"code": code.name, "kind": "mwpm", "title": title, "caption": caption,
            "error": [int(q) for q in error_qubits], "lit": [int(zidx[r]) for r in np.nonzero(s)[0]],
            "pairs": pairs, "correction": [int(q) for q in np.nonzero(corr)[0]], "logical_error": fail}


def toric_tie(code: StabilizerCode) -> list[int]:
    """
    A weight-2 X error whose two lit Z checks sit on a diagonal, so a second,
    equally short correction exists (picked as close to the picture centre as possible).
    """
    H = code.HZ.astype(int)
    seen: dict[tuple, list] = {}
    for a, b in itertools.combinations(range(code.n), 2):
        e = np.zeros(code.n, int); e[[a, b]] = 1
        s = tuple((H @ e) % 2)
        if sum(s) == 2:
            seen.setdefault(s, []).append((a, b))
    pos = code.geometry.qubit_pos
    centre = pos.mean(0)
    # prefer ties whose two qubits are neighbours in the picture (no wrap-around), easier to read
    ties = [v[0] for v in seen.values() if len(v) >= 2
            and np.linalg.norm(pos[v[0][0]] - pos[v[0][1]]) < 2.5]
    best = min(ties, key=lambda ab: np.linalg.norm(code.geometry.qubit_pos[list(ab)].mean(0) - centre))
    return [int(best[0]), int(best[1])]


def default_replays() -> list[dict]:
    """The decoding stories shown in the explorer (same cases as the notebooks)."""
    out = []
    s5 = surface_rotated(5)
    out.append(mwpm_replay(s5, [6, 12, 13], "MWPM fixes a weight-3 chain (surface d=5)",
                           "Three X errors light Z checks at the two ends of the chain. Matching pairs "
                           "the lit checks along the shortest path, so the correction is the same chain."))
    out.append(mwpm_replay(s5, [0, 5, 10], "MWPM loses to 3 errors down a column (surface d=5)",
                           "Three X errors from the top edge: only one check lights, at the bottom of the "
                           "chain. Going to the far boundary costs 2, back to the near one costs 3, so "
                           "matching picks the far side - error + correction spans the code: logical error. "
                           "A distance-5 code fixes at most 2 errors."))
    out.append(mwpm_replay(s5, [0, 1, 2], "Three errors that are really one (surface d=5)",
                           "X on qubits 0, 1, 2 lights a single check. Qubits 1 and 2 together are a "
                           "weight-2 boundary stabilizer, so this error is equivalent to X on qubit 0 alone. "
                           "Matching corrects qubit 0 and the code is fine: degeneracy you can see."))
    t = toric(4)
    out.append(bp_replay(t, toric_tie(t), 0.05, "BP stuck on a tie (toric 4x4)",
                         "Two lit checks on a diagonal: two equally short corrections. BP gives all "
                         "four qubits the same probability and never chooses (degeneracy)."))
    g = bb("gross")
    rng = np.random.default_rng(2026)
    picked = {"ok": None, "bad": None}
    bp = MinSumBP(g.HZ, 0.03, max_iter=40)
    while picked["ok"] is None or picked["bad"] is None:
        qs = rng.choice(g.n, 4, replace=False)
        e = np.zeros(g.n, np.uint8); e[qs] = 1
        r = bp.decode((g.HZ.astype(int) @ e) % 2)
        key = "ok" if (r.converged and r.iterations >= 3) else "bad" if not r.converged else None
        if key and picked[key] is None:
            picked[key] = sorted(int(q) for q in qs)
    out.append(bp_replay(g, picked["ok"], 0.03, "BP converges (gross code)",
                         "A weight-4 error: the probability gathers on the true error and every check is satisfied.",
                         max_iter=40))
    out.append(bp_replay(g, picked["bad"], 0.03, "BP oscillates (gross code)",
                         "Also weight 4 - well inside the distance 12 - but BP never settles: short loops "
                         "in the Tanner graph make its beliefs flip back and forth. This is why the "
                         "gross code is decoded with BP+OSD or Relay-BP.", max_iter=40))
    return out


def build_atlas_data(use_stim: bool = True) -> dict:
    """use_stim=False skips the slow Stim distance searches (tests use it)."""
    codes = []
    for c in catalogue():
        dist = distance_info(c, use_stim=use_stim) if use_stim else \
            {"d": c.info.get("d_known"), "how": "published / brute-force value"}
        codes.append(code_to_dict(c, dist=dist))
    return {"codes": codes, "replays": default_replays()}


def write_atlas_json(path: str) -> str:
    with open(path, "w") as f:
        json.dump(build_atlas_data(), f, separators=(",", ":"))
    return path


def build_explorer_html(template_path: str, out_path: str, standalone: bool = True,
                        use_stim: bool = True) -> str:
    """
    Inline the atlas data into the explorer template -> one self-contained HTML file.
    standalone=True wraps it in a full HTML document (open it straight from disk);
    False writes the bare fragment used when publishing it as a claude.ai artifact.
    """
    data = json.dumps(build_atlas_data(use_stim), separators=(",", ":"))
    html = open(template_path, encoding="utf-8").read().replace("/*__ATLAS_DATA__*/null", data)
    # the live decoders (explorer/decoders.js, next to the template) go inline too,
    # so the built page stays one self-contained file
    js_path = os.path.join(os.path.dirname(os.path.abspath(template_path)), "decoders.js")
    if "/*__DECODERS_JS__*/" in html and os.path.exists(js_path):
        html = html.replace("/*__DECODERS_JS__*/", open(js_path, encoding="utf-8").read())
    if standalone:
        html = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">'
                '<style>body{margin:0}</style></head><body>\n' + html + "\n</body></html>\n")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    return out_path
