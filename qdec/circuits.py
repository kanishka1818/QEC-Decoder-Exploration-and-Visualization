"""
circuits.py - turn ANY CSS code into a noisy memory experiment (Stim circuit),
written by us from the check matrices.

What this file is FOR
---------------------
Code-capacity pictures assume perfect checks.  Real hardware measures checks
with ancilla qubits and CNOTs, and every gate can fail.  This file writes
that experiment, gate by gate, for any CSS StabilizerCode:

    round r:
      reset ancillas (X-check ancillas in |+>, Z-check ancillas in |0>)
      X checks : CNOT ancilla -> data, one "layer" at a time
      measure X-check ancillas in the X basis
      Z checks : CNOT data -> ancilla, layer by layer
      measure Z-check ancillas
      DETECTORs : compare each check with the same check one round earlier
    finally: measure every data qubit, rebuild the checks of the memory basis
             from the data, and declare the logical observables.

Stim is used only as the *simulator* (it runs the circuit we write).  The
same structure as Kanishka's `gross_from_scratch.py`, generalised from the
gross code's six permutations to any check matrix.

Layers (which CNOTs run at the same time)
-----------------------------------------
Two gates in one layer must not share a qubit.  For each check we sort its
qubits by DIRECTION (angle from the check to the qubit, using the torus
picture if there is one).  On translation-invariant codes - surface, toric,
bivariate bicycle - "the k-th direction" is the same shift for every check,
so it is automatically a clash-free layer: exactly the "one polynomial term =
one CNOT layer" idea.  Elsewhere any clash is split off into an extra layer.

Noise model (uniform circuit noise, strength p)
-----------------------------------------------
    DEPOLARIZE2(p) after every CNOT, DEPOLARIZE1(p) on idle qubits each layer,
    reset and measurement flips with probability p, data start with flips p.
Same family as Stim's `after_clifford_depolarization` / SD6-style models.

Honest limitation: X checks are done before Z checks (always correct, but not
the depth-optimal interleaved order).  The circuit distance can therefore be
below the code distance because of hook errors; `circuit_distance` measures it.

References: Fowler et al., PRA 86, 032324 (2012); Gidney, Quantum 5, 497
(2021) (Stim); Bravyi et al. 2024 (gross code circuits).
"""

from __future__ import annotations

import numpy as np
import stim

from .core import StabilizerCode


# --------------------------------------------------------------------------
# CNOT layers
# --------------------------------------------------------------------------
def cnot_layers(code: StabilizerCode, kind: str, order: str = "auto") -> list[list[tuple[int, int]]]:
    """
    Group the (check, qubit) pairs of one check type into layers with no qubit
    or check used twice.  Returns a list of layers, each a list of
    (row index in H_kind, data qubit).

    order='auto'  : if the code supplies `info['cnot_order'][kind]` - a list of
                    direction signs like (-1,-1) = top-left - use that fixed
                    order (this is how the rotated surface code gets its
                    hook-safe "Z-shape for X checks, N-shape for Z checks").
    order='plain' : ignore it and sort by angle (shows the hook-error problem).
    """
    H = code.HX if kind == "X" else code.HZ
    gen_idx = code.index_of(kind)
    g = code.geometry
    fixed = code.info.get("cnot_order", {}).get(kind) if order == "auto" else None
    if fixed is not None:
        layers = [[] for _ in fixed]
        for r, gi in enumerate(gen_idx):
            c = g.check_pos[gi]
            for q in np.nonzero(H[r])[0]:
                v = tuple(int(s) for s in np.sign(g.nearest_image(c, g.qubit_pos[q]) - c))
                layers[fixed.index(v)].append((r, int(q)))
        return [L for L in layers if L]
    ordered = []                      # per check: qubits sorted by direction
    for r, gi in enumerate(gen_idx):
        c = g.check_pos[gi]
        qs = np.nonzero(H[r])[0]
        vecs = np.array([g.nearest_image(c, g.qubit_pos[q]) - c for q in qs])
        ang = np.round(np.arctan2(vecs[:, 1], vecs[:, 0]), 6)
        dist = np.round(np.linalg.norm(vecs, axis=1), 6)
        order = np.lexsort((dist, ang))          # by angle, then by distance
        ordered.append([(r, int(qs[o])) for o in order])
    layers: list[list[tuple[int, int]]] = []
    depth = max(len(o) for o in ordered) if ordered else 0
    for k in range(depth):                        # k-th direction of every check
        pending = [o[k] for o in ordered if k < len(o)]
        while pending:                            # split clashes into extra layers
            used_q, used_c, layer, rest = set(), set(), [], []
            for r, q in pending:
                if q in used_q or r in used_c:
                    rest.append((r, q))
                else:
                    layer.append((r, q)); used_q.add(q); used_c.add(r)
            layers.append(layer)
            pending = rest
    return layers


# --------------------------------------------------------------------------
# the circuit
# --------------------------------------------------------------------------
def memory_circuit(code: StabilizerCode, rounds: int = 3, p: float = 0.0, basis: str = "Z",
                   order: str = "auto") -> stim.Circuit:
    """
    Memory experiment for a CSS code.

    basis 'Z': data start in |0...0>, Z checks give detectors from round 1,
               the final data readout is in Z, observables = Z logicals.
    basis 'X': the mirror image.
    Qubit numbering: data 0..n-1, then X-check ancillas, then Z-check ancillas.
    QUBIT_COORDS and DETECTOR coordinates (x, y, round) come from the geometry.
    """
    if not code.is_css:
        raise ValueError("memory_circuit needs a CSS code (the 5-qubit code mixes X and Z in one check)")
    n = code.n
    HX, HZ = code.HX, code.HZ
    mx, mz = HX.shape[0], HZ.shape[0]
    XA = list(range(n, n + mx))
    ZA = list(range(n + mx, n + mx + mz))
    DATA = list(range(n))
    ALL = DATA + XA + ZA
    g = code.geometry
    xi, zi = code.index_of("X"), code.index_of("Z")
    lay_x, lay_z = cnot_layers(code, "X", order), cnot_layers(code, "Z", order)

    c = stim.Circuit()
    for q in DATA:
        c.append("QUBIT_COORDS", [q], list(map(float, g.qubit_pos[q])))
    for r, a in enumerate(XA):
        c.append("QUBIT_COORDS", [a], list(map(float, g.check_pos[xi[r]])))
    for r, a in enumerate(ZA):
        c.append("QUBIT_COORDS", [a], list(map(float, g.check_pos[zi[r]])))

    def layer(pairs, x_type):
        flat = []
        for r, q in pairs:
            flat += [XA[r], q] if x_type else [q, ZA[r]]   # X check: ancilla controls
        c.append("CX", flat)
        if p:
            c.append("DEPOLARIZE2", flat, p)
            idle = sorted(set(ALL) - set(flat))
            if idle:
                c.append("DEPOLARIZE1", idle, p)
        c.append("TICK")

    # initial data preparation
    c.append("R" if basis == "Z" else "RX", DATA)
    if p:
        c.append("X_ERROR" if basis == "Z" else "Z_ERROR", DATA, p)
    c.append("TICK")

    first, second = (ZA, XA) if basis == "Z" else (XA, ZA)        # 'first' = deterministic from round 1
    for rnd in range(rounds):
        if mx:
            c.append("RX", XA)
        if mz:
            c.append("R", ZA)
        if p:
            if mx: c.append("Z_ERROR", XA, p)
            if mz: c.append("X_ERROR", ZA, p)
        c.append("TICK")
        for L in lay_x:
            layer(L, True)
        if mx:
            c.append("MX", XA, p)
        for L in lay_z:
            layer(L, False)
        if mz:
            c.append("M", ZA, p)
        # measurement record per round: [X ancillas (mx)] then [Z ancillas (mz)]
        per = mx + mz
        for kind, anc, rows in [("X", XA, xi), ("Z", ZA, zi)]:
            off = 0 if kind == "X" else mx                     # offset inside this round's block
            deterministic_now = (anc is first) or rnd > 0
            if not deterministic_now:
                continue
            for r in range(len(anc)):
                t = [stim.target_rec(-per + off + r)]
                if rnd > 0:
                    t.append(stim.target_rec(-2 * per + off + r))
                x, y = g.check_pos[rows[r]]
                c.append("DETECTOR", t, [float(x), float(y), float(rnd)])

    # final readout of the data in the memory basis
    c.append("M" if basis == "Z" else "MX", DATA, p)
    H_mem, rows_mem, anc_mem = (HZ, zi, ZA) if basis == "Z" else (HX, xi, XA)
    off = mx if basis == "Z" else 0
    per = mx + mz
    for r in range(H_mem.shape[0]):
        t = [stim.target_rec(-n + q) for q in np.nonzero(H_mem[r])[0]]
        t.append(stim.target_rec(-n - per + off + r))           # same check, last round
        x, y = g.check_pos[rows_mem[r]]
        c.append("DETECTOR", t, [float(x), float(y), float(rounds)])
    LX, LZ = code.logicals()
    L = LZ[:, n:] if basis == "Z" else LX[:, :n]
    for k, row in enumerate(L):
        c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-n + q) for q in np.nonzero(row)[0]], k)
    return c


def circuit_stats(circuit: stim.Circuit) -> dict:
    """Counts that matter when comparing schedules."""
    counts = {"CX": 0, "ticks": 0}
    for inst in circuit.flattened():
        if inst.name == "CX":
            counts["CX"] += len(inst.targets_copy()) // 2
        elif inst.name == "TICK":
            counts["ticks"] += 1
    counts["detectors"] = circuit.num_detectors
    counts["observables"] = circuit.num_observables
    counts["qubits"] = circuit.num_qubits
    return counts


def circuit_distance(circuit: stim.Circuit, max_degree: int = 4) -> int:
    """
    Smallest number of circuit faults that flip an observable without lighting
    any detector (Stim's search; graph-like errors use the fast exact method).
    """
    dem = circuit.detector_error_model(decompose_errors=False, approximate_disjoint_errors=True)
    try:
        err = circuit.shortest_graphlike_error(ignore_ungraphlike_errors=False)
        graphlike = len(err)
    except ValueError:
        graphlike = None
    found = circuit.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=max_degree,
        dont_explore_edges_with_degree_above=max_degree,
        dont_explore_edges_increasing_symptom_degree=False,
        canonicalize_circuit_errors=True)
    return min(x for x in [graphlike, len(found)] if x is not None)


# --------------------------------------------------------------------------
# code distance with Stim (no exponential brute force)
# --------------------------------------------------------------------------
def code_capacity_circuit(code: StabilizerCode, basis: str = "Z", p: float = 0.01) -> stim.Circuit:
    """
    One perfect round of checks with independent data errors: X errors (basis 'Z')
    or Z errors (basis 'X').  Its smallest undetectable logical error IS the code's
    X-type (resp. Z-type) distance, so Stim's distance tools apply directly.
    """
    base = memory_circuit(code, rounds=1, p=0.0, basis=basis)
    lines = str(base).splitlines()
    first_tick = lines.index("TICK")               # right after the data reset
    err = "X_ERROR" if basis == "Z" else "Z_ERROR"
    lines.insert(first_tick + 1, f"{err}({p}) " + " ".join(map(str, range(code.n))))
    return stim.Circuit("\n".join(lines))


def stim_distance(code: StabilizerCode, max_degree: int = 7) -> tuple[int, str]:
    """
    Code distance of a CSS code using Stim, returned as (d, how):

      'exact (shortest path)'  every data error lights <= 2 checks (surface, toric,
                               repetition...): the distance is a shortest-path
                               problem on the matching graph - exact and instant
                               at any size (Stim's shortest_graphlike_error).
      'upper bound (search)'   otherwise (qLDPC): Stim searches for the smallest
                               undetectable logical error, ignoring detector sets
                               larger than `max_degree`, so the answer is an upper
                               bound (it matches the published d for BB codes).
    d = min over the X-type and Z-type distances.
    """
    best, how = None, "exact (shortest path)"
    for basis in ("Z", "X"):
        checks = code.HZ if basis == "Z" else code.HX     # the checks that see this error type
        if checks.shape[0] == 0:
            # nothing detects this error type (repetition code has no X checks):
            # a single flip is already a logical error
            best, how = 1, "exact (no checks of this type)"
            continue
        c = code_capacity_circuit(code, basis)
        graphlike = int(checks.sum(axis=0).max()) <= 2      # every single error lights <= 2 checks
        if graphlike:
            d = len(c.shortest_graphlike_error(ignore_ungraphlike_errors=True))
        else:
            # NOT shortest_graphlike_error here: Stim would split 3+-check errors into
            # pieces, and a piece with no detectors at all looks like a weight-1 logical
            # (we hit exactly this on HGP(Hamming): it reported d = 1).
            d = len(c.search_for_undetectable_logical_errors(
                dont_explore_detection_event_sets_with_size_above=max_degree,
                dont_explore_edges_with_degree_above=max_degree,
                dont_explore_edges_increasing_symptom_degree=False))
            if how.startswith("exact"):
                how = "upper bound (search)"
        best = d if best is None else min(best, d)
    return best, how
