"""
decoders.py - our own decoders, written from scratch, on our own data
structures (dem.DecodingProblem / dem.MatchingGraph).

What this file is FOR
---------------------
Every decoder answers the same question: "given which detectors lit up,
which logical observables were flipped?"  Here are four answers, each one a
step up from the last, each small enough to read in one sitting:

    LookupDecoder     try every fault set of weight 0, 1, 2, ...; store the
                      most likely one per syndrome. Exact for tiny problems,
                      explodes as 2^(number of detectors).            (Lesson 1)
    MWPMDecoder       matching: shortest paths between lit detectors
                      (Dijkstra, ours) + exact minimum-weight perfect matching
                      of the lit detectors (bitmask dynamic programming, ours).
                      Needs a graph: every fault lights <= 2 detectors.
    BPOSDDecoder      our min-sum BP (bp.py) and, if it does not converge,
                      ordered-statistics decoding (OSD-0 / OSD-CS): trust the
                      qubits BP is surest about, solve the rest exactly by
                      GF(2) elimination.  Works on hypergraphs (gross code).
    sinter adapters   wrap the above so `sinter.collect` can run them.

PyMatching and ldpc are NOT imported here - only the tests use them, to check us.

References
----------
* Edmonds, Canad. J. Math. 17, 449 (1965) - matching; Dennis et al.,
  J. Math. Phys. 43, 4452 (2002) - matching decoders for the surface code.
* Dijkstra, Numer. Math. 1, 269 (1959) - shortest paths.
* Held & Karp, J. SIAM 10, 196 (1962) - the subset dynamic-programming idea
  we use for exact matching of few defects.
* Fossorier & Lin, IEEE Trans. Inf. Theory 41, 1379 (1995) - OSD;
  Panteleev & Kalachev, Quantum 5, 585 (2021) - BP+OSD for quantum codes;
  Roffe et al., arXiv:2005.07016 - OSD-CS ("combination sweep").
"""

from __future__ import annotations

import heapq
import itertools
from functools import lru_cache

import numpy as np

from .bp import MinSumBP
from .dem import DecodingProblem, MatchingGraph
from .gf2 import row_reduce


# ==========================================================================
# 1. lookup table
# ==========================================================================
class LookupDecoder:
    """
    Enumerate fault sets in order of weight (number of faults) up to
    `max_weight`; for each syndrome keep the most probable set found.
    Prediction = observables flipped by that set.  Exact maximum-likelihood
    *single-error-set* decoding for small problems.
    """

    def __init__(self, problem: DecodingProblem, max_weight: int = 2):
        H, L, pri = problem.H, problem.L, problem.priors
        self.num_det = H.shape[0]
        llr = np.log((1 - pri) / pri)                       # cost of using each fault
        self.table: dict[bytes, tuple[float, np.ndarray]] = {}
        zero = np.zeros(self.num_det, np.uint8)
        self.table[zero.tobytes()] = (0.0, np.zeros(L.shape[0], np.uint8))
        cols = H.T.astype(np.uint8)
        lcols = L.T.astype(np.uint8)
        for w in range(1, max_weight + 1):
            for combo in itertools.combinations(range(H.shape[1]), w):
                s = np.bitwise_xor.reduce(cols[list(combo)], axis=0)
                key = s.tobytes()
                cost = float(llr[list(combo)].sum())
                if key not in self.table or cost < self.table[key][0]:
                    self.table[key] = (cost, np.bitwise_xor.reduce(lcols[list(combo)], axis=0))

    def decode(self, syndrome) -> np.ndarray:
        """Observable prediction; all-zeros if the syndrome is not in the table."""
        hit = self.table.get(np.asarray(syndrome, np.uint8).tobytes())
        return hit[1] if hit else np.zeros(len(next(iter(self.table.values()))[1]), np.uint8)

    def __len__(self):
        return len(self.table)


# ==========================================================================
# 2. minimum-weight perfect matching (ours)
# ==========================================================================
class MWPMDecoder:
    """
    Matching decoder on a MatchingGraph.

    Step 1 (Dijkstra): from every lit detector, shortest path to every other
            lit detector and to the boundary; remember the observable parity
            along each path (XOR of edge masks).
    Step 2 (exact matching): choose pairs (or "go to boundary") minimising
            total length.  With k lit detectors we solve it exactly by dynamic
            programming over subsets:  best(S) = min over partners of the
            lowest detector in S.  Cost 2^k * k - fine up to ~20 lit detectors.
            Above `max_exact` we solve the same pairing exactly as a small
            integer program with scipy's MILP solver (counted in `self.ilp_calls`);
            only if scipy is missing or fails do we fall back to greedy nearest
            pairing (counted in `self.greedy_calls`: not exact any more).
    Prediction = XOR of the observable parities of the chosen paths.
    """

    def __init__(self, graph: MatchingGraph, max_exact: int = 20):
        self.g = graph
        self.max_exact = max_exact
        self.greedy_calls = 0
        self.last_matching: list[tuple[int, int]] = []        # for pictures
        self.last_paths: dict[tuple[int, int], list[int]] = {}

    def _dijkstra(self, src: int, targets: set[int]):
        """Shortest distance, observable parity and path from src to every target."""
        g = self.g
        dist = {src: 0.0}
        par = {src: 0}
        prev = {src: -1}
        heap = [(0.0, src)]
        remaining = set(targets)
        while heap and remaining:
            d, u = heapq.heappop(heap)
            if d > dist.get(u, np.inf):
                continue
            remaining.discard(u)
            for v, w, m in g.adj[u]:
                nd = d + w
                if nd < dist.get(v, np.inf):
                    dist[v], par[v], prev[v] = nd, par[u] ^ m, u
                    heapq.heappush(heap, (nd, v))
        out = {}
        for t in targets:
            if t in dist:
                path, x = [], t
                while x != -1:
                    path.append(x); x = prev[x]
                out[t] = (dist[t], par[t], path[::-1])
            else:
                out[t] = (np.inf, 0, [])
        return out

    def decode(self, syndrome) -> int:
        """Returns the predicted observable flips as a bit mask (bit k = observable k)."""
        lit = [int(i) for i in np.nonzero(np.asarray(syndrome))[0]]
        self.last_matching, self.last_paths = [], {}
        if not lit:
            return 0
        B = self.g.boundary
        info = {u: self._dijkstra(u, set(lit) | {B}) for u in lit}
        k = len(lit)
        if k > self.max_exact:
            return self._exact_ilp(lit, info)

        @lru_cache(maxsize=None)
        def best(mask: int):
            """(cost, choice list) to match the detectors in `mask`."""
            if mask == 0:
                return 0.0, ()
            i = (mask & -mask).bit_length() - 1                 # lowest remaining detector
            rest = mask & ~(1 << i)
            c0, ch0 = best(rest)                                 # option A: i -> boundary
            cand = (c0 + info[lit[i]][B][0], ((i, -1),) + ch0)
            j_mask = rest
            while j_mask:                                        # option B: i <-> j
                j = (j_mask & -j_mask).bit_length() - 1
                j_mask &= j_mask - 1
                cj, chj = best(rest & ~(1 << j))
                c = cj + info[lit[i]][lit[j]][0]
                if c < cand[0]:
                    cand = (c, ((i, j),) + chj)
            return cand

        _, choice = best((1 << k) - 1)
        return self._collect(lit, info, choice)

    def _collect(self, lit, info, choice):
        B = self.g.boundary
        obs = 0
        for i, j in choice:
            a = lit[i]
            b = B if j == -1 else lit[j]
            _, parity, path = info[a][b]
            obs ^= parity
            self.last_matching.append((a, b))
            self.last_paths[(a, b)] = path
        return obs

    def _exact_ilp(self, lit, info):
        """
        Many lit detectors: the same matching problem written as a small integer
        program (one 0/1 variable per pair and per "go to boundary"; every lit
        detector used exactly once) and solved exactly with scipy's MILP solver
        (HiGHS).  This is the one place we lean on a general-purpose solver; a
        from-scratch blossom algorithm is on the to-do list.  Falls back to
        greedy if scipy is missing.
        """
        try:
            from scipy.optimize import Bounds, LinearConstraint, milp
        except ImportError:
            return self._greedy(lit, info)
        self.ilp_calls = getattr(self, "ilp_calls", 0) + 1
        B = self.g.boundary
        k = len(lit)
        vars_, cost = [], []
        for i in range(k):
            vars_.append((i, -1)); cost.append(info[lit[i]][B][0])
            for j in range(i + 1, k):
                cost.append(info[lit[i]][lit[j]][0]); vars_.append((i, j))
        cost = np.array(cost)
        finite = np.isfinite(cost)
        cost = np.where(finite, cost, 1e6)
        A = np.zeros((k, len(vars_)))
        for v, (i, j) in enumerate(vars_):
            A[i, v] = 1
            if j != -1:
                A[j, v] = 1
        res = milp(cost, constraints=LinearConstraint(A, 1, 1), integrality=np.ones(len(vars_)),
                   bounds=Bounds(0, 1))
        if not res.success:
            return self._greedy(lit, info)
        choice = [vars_[v] for v in np.nonzero(res.x > 0.5)[0]]
        return self._collect(lit, info, choice)

    def _greedy(self, lit, info):
        self.greedy_calls += 1
        B = self.g.boundary
        left = set(range(len(lit)))
        choice = []
        while left:
            i = min(left)
            left.discard(i)
            options = [(info[lit[i]][B][0], -1)] + [(info[lit[i]][lit[j]][0], j) for j in left]
            c, j = min(options)
            if j != -1:
                left.discard(j)
            choice.append((i, j))
        return self._collect(lit, info, choice)


# ==========================================================================
# 3. BP + OSD (ours)
# ==========================================================================
class BPOSDDecoder:
    """
    BP first; if BP's guess does not reproduce the syndrome, OSD:

      1. sort faults from "BP is sure it happened / didn't" ... to least sure
         (by BP's final flip probability, most likely first)
      2. Gaussian-eliminate H with columns in that order; the pivot columns
         are the most reliable set of faults that can explain ANY syndrome
      3. OSD-0: solve H e = s using only pivot columns (others = 0)
      4. OSD-CS (order w): also try flipping each single, and each pair, of
         the first w non-pivot columns; keep the most likely solution
    Works on any H, including the gross code's hypergraph.
    """

    def __init__(self, problem: DecodingProblem, max_iter: int = 50, osd_order: int = 0,
                 alpha: float = 0.625):
        self.pb = problem
        self.H = problem.H.astype(np.uint8)
        self.L = problem.L.astype(np.uint8)
        p = float(np.median(problem.priors))
        self.bp = MinSumBP(self.H, p, max_iter=max_iter, alpha=alpha)
        self.bp_prior = np.log((1 - problem.priors) / problem.priors)   # per-fault priors
        self.order = osd_order
        self.osd_calls = 0
        self.last = {}

    def _bp(self, s):
        """Uniform priors: the plain MinSumBP. Different prior per fault (circuit noise): weighted loop."""
        if np.allclose(self.pb.priors, self.pb.priors[0]):
            return self.bp.decode(s)
        return self._bp_weighted(s)

    def _bp_weighted(self, s):
        """Same min-sum loop as bp.MinSumBP but with a different prior per fault."""
        bp = self.bp
        prior = self.bp_prior
        llr = prior.copy()
        m_vc = np.where(bp.mask, prior[bp.var], 0.0)
        guess = (llr < 0).astype(np.uint8)
        hist = [llr.copy()]
        it, conv = 0, not s.any()
        while not conv and it < bp.max_iter:
            it += 1
            m_cv = bp._check_to_var(m_vc, s)
            llr = prior.copy()
            np.add.at(llr, bp.var[bp.mask], m_cv[bp.mask])
            m_vc = np.where(bp.mask, llr[bp.var] - m_cv, 0.0)
            guess = (llr < 0).astype(np.uint8)
            hist.append(llr.copy())
            conv = not np.any((self.H.astype(np.int64) @ guess) % 2 != s)
        from .bp import BPResult
        return BPResult(guess, conv, it, np.array(hist), (np.array(hist) < 0).astype(np.uint8),
                        np.zeros(len(hist), int))

    def decode_error(self, syndrome) -> np.ndarray:
        """Most likely fault vector found (0/1 per column of H)."""
        s = (np.asarray(syndrome) % 2).astype(np.uint8)
        res = self._bp(s)
        self.last = {"bp": res, "osd": False}
        if res.converged:
            return res.correction
        self.osd_calls += 1
        # rank faults by BP's belief AVERAGED over the iterations: when BP oscillates
        # (exactly the shots that reach OSD) the last iteration is a coin toss,
        # the average is much steadier.  Measured in the notebook: this alone
        # roughly halves OSD's failures on the gross code.
        return self._osd(s, res.llr_history[1:].mean(0) if res.iterations else res.llr_history[0])

    def decode(self, syndrome) -> np.ndarray:
        e = self.decode_error(syndrome)
        return (self.L.astype(np.int64) @ e % 2).astype(np.uint8)

    def _osd(self, s, llr):
        H = self.H
        n = H.shape[1]
        order = np.argsort(llr, kind="stable")     # stable: ties keep qubit order (the explorer JS does the same); most negative LLR = most likely flipped, first
        Hs = H[:, order]
        R, piv = row_reduce(np.hstack([Hs, s[:, None]]))
        if n in piv:                               # syndrome not reachable (should not happen)
            piv = [c for c in piv if c != n]
        piv = [c for c in piv if c < n]
        piv_set = set(piv)
        nonpiv = [c for c in range(n) if c not in piv_set]
        # score of a candidate = -log(probability) up to a constant = sum of the prior
        # weights log((1-p)/p) of the faults it switches on (Hamming weight if all p equal)
        prior_w = self.bp_prior[order]
        cost_of = lambda e: float(np.dot(e.astype(float), prior_w))

        def solve_with(extra):
            """pivot variables from reduced rows, given values of the chosen non-pivot columns."""
            e = np.zeros(n, np.uint8)
            rhs = R[:len(piv), n].copy()
            for c in extra:
                e[c] = 1
                rhs ^= R[:len(piv), c]
            for r, c in enumerate(piv):
                e[c] = rhs[r]
            return e

        best = solve_with([])
        best_cost = cost_of(best)
        if self.order > 0:
            # OSD-CS: every single non-pivot flip, then every pair among the first `order`
            trials = [(c,) for c in nonpiv] + list(itertools.combinations(nonpiv[: self.order], 2))
            for extra in trials:
                e = solve_with(extra)
                c = cost_of(e)
                if c < best_cost:
                    best, best_cost = e, c
        out = np.zeros(n, np.uint8)
        out[order] = best
        self.last.update(osd=True, order=order, pivots=piv, osd_solution=out)
        return out


# ==========================================================================
# 4. sinter adapters (so sinter.collect can sweep OUR decoders)
# ==========================================================================
# sinter runs decoders in worker processes, which must be able to import the
# classes by name - so they live at module level.
try:
    import sinter as _sinter

    class _Compiled(_sinter.CompiledDecoder):
        """Unpacks sinter's bit-packed shots, runs `fn` per shot, packs predictions."""

        def __init__(self, fn, num_det, num_obs):
            self.fn, self.nd, self.no = fn, num_det, num_obs

        def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: np.ndarray) -> np.ndarray:
            dets = np.unpackbits(bit_packed_detection_event_data, axis=1, bitorder="little")[:, :self.nd]
            out = np.zeros((dets.shape[0], self.no), np.uint8)
            for i, s in enumerate(dets):
                out[i] = self.fn(s)
            return np.packbits(out, axis=1, bitorder="little")

    class SinterOwnMWPM(_sinter.Decoder):
        """Our MWPM (needs a graph-like DEM: surface/toric/repetition codes)."""

        def compile_decoder_for_dem(self, *, dem):
            d = MWPMDecoder(MatchingGraph(dem))
            no = dem.num_observables

            def fn(s):
                m = d.decode(s)
                return np.array([(m >> k) & 1 for k in range(no)], np.uint8)
            return _Compiled(fn, dem.num_detectors, no)

    class SinterOwnBPOSD(_sinter.Decoder):
        """Our BP+OSD-CS (any DEM, including hypergraphs)."""

        def __init__(self, osd_order: int = 4):
            self.osd_order = osd_order

        def compile_decoder_for_dem(self, *, dem):
            d = BPOSDDecoder(DecodingProblem.from_dem(dem), osd_order=self.osd_order)
            return _Compiled(d.decode, dem.num_detectors, dem.num_observables)

except ImportError:          # sinter not installed: adapters simply unavailable
    SinterOwnMWPM = SinterOwnBPOSD = None


def sinter_decoders() -> dict:
    """{'own_mwpm': ..., 'own_bposd': ...} for sinter.collect(custom_decoders=...)."""
    if SinterOwnMWPM is None:
        raise ImportError("sinter is not installed")
    return {"own_mwpm": SinterOwnMWPM(), "own_bposd": SinterOwnBPOSD()}
