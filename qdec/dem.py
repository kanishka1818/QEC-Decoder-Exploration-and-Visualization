"""
dem.py - read Stim's detector error model into OUR OWN graph objects,
and see why the surface code gives a graph but the gross code a hypergraph.

What this file is FOR
---------------------
After simulating a circuit, a decoder never sees qubits again.  It sees:
    * detectors   - "this check changed since last round" (the nodes),
    * fault mechanisms - each possible fault, its probability, which
      detectors it lights and which logical observables it flips.
Stim calls this list the detector error model (DEM).  We parse it ourselves
and build:

    DecodingProblem   H (detectors x faults), L (observables x faults), priors
                      -> the input to BP / OSD / lookup decoders
    MatchingGraph     nodes = detectors (+ one boundary node), edges = faults
                      that light at most 2 detectors -> input to our MWPM
    fault_degree_histogram   how many detectors each fault lights.
                      All <= 2  -> a graph  (matching works)
                      Some >= 3 -> a hypergraph (matching does not apply)

A Y error on a surface-code qubit lights both an X-type and a Z-type
detector pair (4 detectors).  Stim can split such a fault into two graph-like
pieces (`decompose_errors=True`, written with '^'); we keep both views.

Reference: Gidney, "Stim: a fast stabilizer circuit simulator", Quantum 5, 497 (2021);
Stim docs, "The Detector Error Model File Format (.dem)".
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import stim


@dataclass
class Fault:
    p: float
    detectors: tuple[int, ...]
    observables: tuple[int, ...]
    parts: list[tuple[tuple[int, ...], tuple[int, ...]]] = field(default_factory=list)  # '^' pieces


def parse_dem(dem: stim.DetectorErrorModel) -> tuple[list[Fault], int, int, dict]:
    """
    Flatten the DEM and list every fault.  Returns (faults, num_detectors,
    num_observables, detector_coords).  Each fault's `parts` holds the
    graph-like pieces if Stim decomposed it with '^'.
    """
    faults = []
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        parts, cur_d, cur_o = [], [], []
        for t in inst.targets_copy():
            if t.is_separator():
                parts.append((tuple(cur_d), tuple(cur_o))); cur_d, cur_o = [], []
            elif t.is_relative_detector_id():
                cur_d.append(t.val)
            elif t.is_logical_observable_id():
                cur_o.append(t.val)
        parts.append((tuple(cur_d), tuple(cur_o)))
        # whole fault = symmetric difference of the pieces
        dset, oset = Counter(), Counter()
        for dd, oo in parts:
            dset.update(dd); oset.update(oo)
        dets = tuple(sorted(k for k, v in dset.items() if v % 2))
        obs = tuple(sorted(k for k, v in oset.items() if v % 2))
        faults.append(Fault(p, dets, obs, parts if len(parts) > 1 else []))
    coords = dem.get_detector_coordinates()
    return faults, dem.num_detectors, dem.num_observables, coords


# --------------------------------------------------------------------------
# matrix view (for BP / OSD / lookup)
# --------------------------------------------------------------------------
@dataclass
class DecodingProblem:
    H: np.ndarray          # (detectors, faults) 0/1
    L: np.ndarray          # (observables, faults) 0/1
    priors: np.ndarray     # (faults,) probabilities
    coords: dict

    @classmethod
    def from_dem(cls, dem: stim.DetectorErrorModel) -> "DecodingProblem":
        faults, nd, no, coords = parse_dem(dem)
        H = np.zeros((nd, len(faults)), np.uint8)
        L = np.zeros((no, len(faults)), np.uint8)
        for j, f in enumerate(faults):
            H[list(f.detectors), j] = 1
            L[list(f.observables), j] = 1
        return cls(H, L, np.array([f.p for f in faults]), coords)

    @classmethod
    def code_capacity(cls, H, L, p: float) -> "DecodingProblem":
        """Same object for code-capacity noise: faults = data qubits, detectors = checks."""
        H = np.asarray(H, np.uint8)
        return cls(H, np.asarray(L, np.uint8), np.full(H.shape[1], p), {})


def fault_degree_histogram(dem: stim.DetectorErrorModel, use_parts: bool = False) -> Counter:
    """
    How many detectors each fault lights.  use_parts=True counts Stim's
    graph-like pieces instead of whole faults.
    """
    faults, *_ = parse_dem(dem)
    hist = Counter()
    for f in faults:
        if use_parts and f.parts:
            for dd, _ in f.parts:
                hist[len(dd)] += 1
        else:
            hist[len(f.detectors)] += 1
    return hist


# --------------------------------------------------------------------------
# graph view (for our MWPM)
# --------------------------------------------------------------------------
class MatchingGraph:
    """
    Nodes 0..D-1 = detectors, node D = the boundary.
    Edge (u, v) = a fault lighting exactly u and v; a fault lighting only u
    becomes the edge (u, boundary).  Parallel edges are merged (independent
    faults: p = p1(1-p2) + p2(1-p1)), keeping the observable mask of the more
    likely one.  weight = log((1-p)/p): the most likely explanation is the
    lightest set of edges.

    Faults lighting 3+ detectors are split into Stim's '^' pieces if present,
    otherwise counted in `self.dropped` (the decoder cannot use them).
    """

    def __init__(self, dem: stim.DetectorErrorModel | None = None, *, _faults=None):
        if _faults is not None:                      # built by from_check_matrix
            faults, nd, no, coords = _faults
        else:
            faults, nd, no, coords = parse_dem(dem)
        self.num_detectors, self.num_observables = nd, no
        self.boundary = nd
        self.coords = coords
        self.dropped = 0
        edges: dict[tuple[int, int], list] = {}
        self.edge_fault: dict[tuple[int, int], int] = {}   # edge -> index of the fault drawn for it

        def add(dets, obs, p, fault_index):
            if len(dets) == 0:
                return
            u, v = (dets[0], self.boundary) if len(dets) == 1 else (min(dets), max(dets))
            mask = 0
            for o in obs:
                mask ^= 1 << o
            if (u, v) in edges:
                p0, m0 = edges[(u, v)]
                merged = p0 * (1 - p) + p * (1 - p0)
                edges[(u, v)] = [merged, m0 if p0 >= p else mask]
                if p > p0:
                    self.edge_fault[(u, v)] = fault_index
            else:
                edges[(u, v)] = [p, mask]
                self.edge_fault[(u, v)] = fault_index

        for fi, f in enumerate(faults):
            pieces = f.parts if f.parts else [(f.detectors, f.observables)]
            if any(len(dd) > 2 for dd, _ in pieces):
                self.dropped += 1
                continue
            for dd, oo in pieces:
                add(dd, oo, f.p, fi)
        self.edges = {k: (p, m, float(np.log((1 - p) / p))) for k, (p, m) in edges.items()}
        self.adj: list[list[tuple[int, float, int]]] = [[] for _ in range(nd + 1)]
        for (u, v), (p, m, w) in self.edges.items():
            self.adj[u].append((v, w, m))
            self.adj[v].append((u, w, m))

    @classmethod
    def from_check_matrix(cls, H, L, p: float) -> "MatchingGraph":
        """
        Code-capacity matching graph straight from a check matrix: nodes = checks,
        one edge per qubit (column) joining the <= 2 checks it touches; a qubit in
        only one check joins that check to the boundary.  L gives which logical
        observables each qubit flips.  edge_fault[(u, v)] = the qubit index.
        """
        H, L = np.asarray(H, np.uint8), np.asarray(L, np.uint8)
        faults = [Fault(p, tuple(np.nonzero(H[:, j])[0]), tuple(np.nonzero(L[:, j])[0]))
                  for j in range(H.shape[1])]
        return cls(_faults=(faults, H.shape[0], L.shape[0], {}))

    def path_faults(self, path: list[int]) -> list[int]:
        """Fault (or qubit) indices along a node path returned by the MWPM decoder."""
        out = []
        for a, b in zip(path[:-1], path[1:]):
            out.append(self.edge_fault[(min(a, b), max(a, b))])
        return out

    def __repr__(self):
        return (f"MatchingGraph({self.num_detectors} detectors + boundary, {len(self.edges)} edges, "
                f"{self.dropped} faults not graph-like)")


# --------------------------------------------------------------------------
# pictures
# --------------------------------------------------------------------------
def draw_detector_graph_3d(dem: stim.DetectorErrorModel, ax, max_faults: int = 4000,
                           title: str = ""):
    """
    Space-time picture: every detector at (x, y, round).  Faults lighting 2
    detectors are drawn as lines (grey = same round, orange = across rounds);
    faults lighting 3+ detectors are drawn as purple stars from their centre
    (the hyperedges).  Single-detector faults are small ticks to the boundary.
    """
    faults, nd, no, coords = parse_dem(dem)
    P = np.array([coords.get(i, [0, 0, 0])[:3] for i in range(nd)], float)
    for f in faults[:max_faults]:
        d = list(f.detectors)
        if len(d) == 2:
            a, b = P[d[0]], P[d[1]]
            col = "0.55" if a[2] == b[2] else "#e07a1f"
            ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color=col, lw=0.5, alpha=0.5)
        elif len(d) >= 3:
            cen = P[d].mean(0)
            for q in d:
                ax.plot([cen[0], P[q][0]], [cen[1], P[q][1]], [cen[2], P[q][2]],
                        color="#8a5cc2", lw=0.4, alpha=0.25)
    ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=8, c=P[:, 2], cmap="viridis")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("round")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    return ax
