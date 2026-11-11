/*
 * decoders.js - the explorer's live decoders, in plain JavaScript.
 *
 * What this file is FOR
 * ---------------------
 * The browser explorer decodes the errors YOU click, step by step. These are
 * small, readable re-implementations of the Python decoders in qdec/:
 *
 *   lookup(code, err)            table look-up (codes with <= 13 qubits)
 *   mwpm(code, errSet, kind, p)  minimum-weight matching on the code-capacity graph
 *   bp(code, errSet, kind, ...)  min-sum belief propagation, every iteration recorded
 *   osd0(code, kind, s, llr)     ordered-statistics step when BP does not settle
 *
 * "kind" is the type of check that SEES the error: Z checks see X errors and
 * X checks see Z errors. A Y error is an X part plus a Z part, decoded separately.
 *
 * The tests (tests/test_explorer_decoders.py) run this file in Node and compare
 * it with qdec's Python decoders on random errors, so the explorer cannot drift
 * away from the library.
 *
 * Works in the browser (global QDec) and in Node (module.exports).
 */
(function (root) {
  "use strict";

  // ------------------------------------------------------------------ basics
  /** Indices of the generators of one kind ('X' or 'Z'). */
  function rowsOf(code, kind) {
    const out = [];
    code.gens.forEach((g, i) => { if (g.kind === kind) out.push(i); });
    return out;
  }

  /** Split an error {qubit: 'X'|'Y'|'Z'} into its X part and Z part (sets of qubits). */
  function splitError(err) {
    const x = new Set(), z = new Set();
    for (const q in err) {
      const l = err[q];
      if (l === "X" || l === "Y") x.add(+q);
      if (l === "Z" || l === "Y") z.add(+q);
    }
    return { x, z };
  }

  /** Syndrome bits of the rows of one kind for an error set (of the other Pauli type). */
  function syndromeOf(code, rows, errSet) {
    return rows.map((r) => code.gens[r].support.reduce((s, q) => s ^ (errSet.has(q) ? 1 : 0), 0));
  }

  /** For each qubit, which of these rows contain it (the columns of H). */
  function columns(code, rows) {
    const cols = Array.from({ length: code.n }, () => []);
    rows.forEach((r, k) => code.gens[r].support.forEach((q) => cols[q].push(k)));
    return cols;
  }

  // ------------------------------------------------------------------ lookup
  /** Lookup decoder: the stored lightest correction for this syndrome, or null. */
  function lookup(code, err) {
    if (!code.lookup) return null;
    const syn = code.gens.map((g) => {
      let s = 0;
      g.support.forEach((q, j) => { const e = err[q], l = g.letters[j]; if (e && e !== l) s ^= 1; });
      return s;
    });
    const key = syn.join("");
    const corr = code.lookup[key];
    return { key, rows: Object.keys(code.lookup).length, correction: corr || null };
  }

  // ------------------------------------------------------------------ MWPM
  /** Can matching be used: every qubit sits in at most 2 checks of this kind? */
  function isGraphLike(code, kind) {
    const rows = rowsOf(code, kind);
    return rows.length > 0 && columns(code, rows).every((c) => c.length <= 2);
  }

  /**
   * Minimum-weight perfect matching on the code-capacity graph - the same graph
   * as qdec.dem.MatchingGraph.from_check_matrix and the same search as
   * qdec.decoders.MWPMDecoder (checked by the tests):
   *   nodes: the checks of `kind` plus one boundary node B
   *   edges: one per data qubit, joining the (one or two) checks that contain it.
   *          Two qubits giving the SAME edge (e.g. both only in one boundary check)
   *          are merged: either flipping explains it, p = p1(1-p2) + p2(1-p1),
   *          so that edge is a little cheaper. weight = log((1-p)/p).
   * Shortest paths: Dijkstra. Pairing: exact dynamic programming over subsets
   * of lit checks (up to 18 lit), greedy beyond that (flagged in the result).
   */
  function mwpm(code, errSet, kind, p = 0.05) {
    const rows = rowsOf(code, kind);
    const cols = columns(code, rows);
    if (!rows.length || cols.some((c) => c.length > 2)) return null;
    const B = rows.length;                                   // boundary node id
    const emap = new Map();                                  // "u,v" -> {u, v, p, q (drawn qubit), qubits}
    cols.forEach((c, q) => {
      if (!c.length) return;
      const u = c.length === 1 ? c[0] : Math.min(c[0], c[1]);
      const v = c.length === 1 ? B : Math.max(c[0], c[1]);
      const key = u + "," + v, e = emap.get(key);
      if (e) { e.p = e.p * (1 - p) + p * (1 - e.p); e.qubits.push(q); }   // first qubit stays the drawn one
      else emap.set(key, { u, v, p, q, qubits: [q] });
    });
    const adj = Array.from({ length: B + 1 }, () => []);
    for (const e of emap.values()) {
      e.w = Math.log((1 - e.p) / e.p);
      adj[e.u].push([e.v, e.w, e.q]); adj[e.v].push([e.u, e.w, e.q]);
    }
    const syn = syndromeOf(code, rows, errSet);
    const lit = syn.map((s, k) => (s ? k : -1)).filter((k) => k >= 0);
    // Dijkstra from one lit node: distance to every node + the qubits along the path.
    // Ties are popped smallest node first, like Python's heap of (distance, node).
    const dijkstra = (src) => {
      const dist = new Array(B + 1).fill(Infinity), via = new Array(B + 1).fill(null), done = new Array(B + 1).fill(false);
      dist[src] = 0;
      for (;;) {
        let u = -1;
        for (let x = 0; x <= B; x++) if (!done[x] && dist[x] < Infinity && (u < 0 || dist[x] < dist[u])) u = x;
        if (u < 0) break;
        done[u] = true;
        for (const [v, w, q] of adj[u]) if (dist[u] + w < dist[v]) { dist[v] = dist[u] + w; via[v] = [u, q]; }
      }
      const pathTo = (t) => { const qs = []; let x = t; while (x !== src && via[x]) { qs.push(via[x][1]); x = via[x][0]; } return qs.reverse(); };
      return { dist, pathTo };
    };
    const info = lit.map((u) => dijkstra(u));
    const k = lit.length;
    let choice = [], greedy = false;
    if (k <= 18) {
      const memo = new Map();
      const best = (mask) => {
        if (mask === 0) return [0, []];
        if (memo.has(mask)) return memo.get(mask);
        const i = 31 - Math.clz32(mask & -mask);             // lowest remaining lit node
        const rest = mask & ~(1 << i);
        const [c0, ch0] = best(rest);
        let cand = [c0 + info[i].dist[B], [[i, -1]].concat(ch0)];
        let m = rest;
        while (m) {
          const j = 31 - Math.clz32(m & -m); m &= m - 1;
          const [cj, chj] = best(rest & ~(1 << j));
          const c = cj + info[i].dist[lit[j]];
          if (c < cand[0]) cand = [c, [[i, j]].concat(chj)];
        }
        memo.set(mask, cand);
        return cand;
      };
      choice = best((1 << k) - 1)[1];
    } else {
      greedy = true;
      const left = new Set(lit.map((_, i) => i));
      while (left.size) {
        const i = Math.min(...left); left.delete(i);
        let bj = -1, bc = info[i].dist[B];
        for (const j of left) if (info[i].dist[lit[j]] < bc) { bc = info[i].dist[lit[j]]; bj = j; }
        if (bj >= 0) left.delete(bj);
        choice.push([i, bj]);
      }
    }
    const correction = new Set(), pairs = [];
    let hops = 0, cost = 0;
    choice.forEach(([i, j]) => {
      const target = j === -1 ? B : lit[j];
      const qs = info[i].pathTo(target);
      hops += qs.length; cost += info[i].dist[target];
      qs.forEach((q) => (correction.has(q) ? correction.delete(q) : correction.add(q)));
      pairs.push({ a: rows[lit[i]], b: j === -1 ? null : rows[lit[j]], qubits: qs });
    });
    // the graph, for drawing: one entry per qubit = edge between its checks
    // cost = total log-weight (what was minimised); hops = qubits flipped along the paths
    return { lit: lit.map((x) => rows[x]), pairs, correction, cost, hops, greedy, edges: [...emap.values()].map((e) => ({ ends: [rows[e.u], e.v === B ? null : rows[e.v]], qubits: e.qubits })) };
  }

  // ------------------------------------------------------------------ BP
  /**
   * Normalised min-sum BP, flooding schedule - the same algorithm and the same
   * numbers as qdec.bp.MinSumBP (checked by the tests).
   */
  function bp(code, errSet, kind, p = 0.05, maxIter = 30, alpha = 0.625) {
    const rows = rowsOf(code, kind);
    const supp = rows.map((r) => code.gens[r].support);
    const s = syndromeOf(code, rows, errSet);
    const n = code.n, prior = Math.log((1 - p) / p);
    let mvc = supp.map((qs) => qs.map(() => prior));        // qubit -> check messages
    let llr = new Array(n).fill(prior);
    const guessOf = (l) => l.map((v) => (v < 0 ? 1 : 0));
    const unsatOf = (g) => supp.reduce((acc, qs, c) => acc + (((qs.reduce((t, q) => t ^ g[q], 0)) !== s[c]) ? 1 : 0), 0);
    let guess = guessOf(llr);
    const hist = { llr: [llr.slice()], guess: [guess.slice()], unsat: [unsatOf(guess)] };
    let converged = hist.unsat[0] === 0 && !s.some((v) => v);
    let it = 0;
    while (!converged && it < maxIter) {
      it += 1;
      // check -> qubit: sign = (-1)^s * product of the OTHER signs; size = alpha * min of the OTHER sizes
      const mcv = supp.map((qs, c) => {
        const mags = mvc[c].map(Math.abs);
        let i1 = 0;
        for (let t = 1; t < mags.length; t++) if (mags[t] < mags[i1]) i1 = t;
        let min2 = Infinity;
        for (let t = 0; t < mags.length; t++) if (t !== i1 && mags[t] < min2) min2 = mags[t];
        const signs = mvc[c].map((v) => (v < 0 ? -1 : 1));
        const total = signs.reduce((a, b) => a * b, 1) * (s[c] ? -1 : 1);
        return qs.map((_, t) => alpha * total * signs[t] * (t === i1 ? min2 : mags[i1]));
      });
      llr = new Array(n).fill(prior);
      supp.forEach((qs, c) => qs.forEach((q, t) => { llr[q] += mcv[c][t]; }));
      mvc = supp.map((qs, c) => qs.map((q, t) => llr[q] - mcv[c][t]));
      guess = guessOf(llr);
      const u = unsatOf(guess);
      hist.llr.push(llr.slice()); hist.guess.push(guess.slice()); hist.unsat.push(u);
      converged = u === 0;
    }
    const probs = hist.llr.map((l) => l.map((v) => 1 / (1 + Math.exp(Math.max(-50, Math.min(50, v))))));
    return { rows, syndrome: s, converged, iterations: it, llr: hist.llr, guess: hist.guess, unsat: hist.unsat, probs };
  }

  // ------------------------------------------------------------------ OSD-0
  /**
   * Ordered-statistics decoding, order 0 (as qdec.decoders.BPOSDDecoder with
   * osd_order=0): sort qubits from "BP thinks most likely flipped" to least
   * (by the given LLRs, stable sort), Gaussian-eliminate H in that column order
   * with the syndrome attached, and solve using only the pivot columns.
   */
  function osd0(code, kind, s, llr) {
    const rows = rowsOf(code, kind);
    const n = code.n, m = rows.length;
    const order = [...Array(n).keys()].sort((a, b) => llr[a] - llr[b]);
    const R = rows.map((r, i) => {
      const row = new Array(n + 1).fill(0);
      code.gens[r].support.forEach((q) => { row[order.indexOf(q)] = 1; });
      row[n] = s[i];
      return row;
    });
    const pivots = [];
    let rr = 0;
    for (let c = 0; c < n && rr < m; c++) {
      let pr = -1;
      for (let i = rr; i < m; i++) if (R[i][c]) { pr = i; break; }
      if (pr < 0) continue;
      [R[rr], R[pr]] = [R[pr], R[rr]];
      for (let i = 0; i < m; i++) if (i !== rr && R[i][c]) for (let t = 0; t <= n; t++) R[i][t] ^= R[rr][t];
      pivots.push(c); rr += 1;
    }
    const sorted = new Array(n).fill(0);
    pivots.forEach((c, i) => { sorted[c] = R[i][n]; });
    const correction = new Set();
    sorted.forEach((v, c) => { if (v) correction.add(order[c]); });
    return { correction, pivots: pivots.map((c) => order[c]), order };
  }

  /**
   * BP, then OSD-0 if BP does not settle - exactly qdec's BPOSDDecoder(osd_order=0):
   * OSD ranks qubits by BP's belief AVERAGED over iterations 1..end (steadier
   * than the last iteration when BP oscillates).
   */
  function bposd(code, errSet, kind, p = 0.05, maxIter = 30, alpha = 0.625) {
    const run = bp(code, errSet, kind, p, maxIter, alpha);
    if (run.converged) {
      const correction = new Set();
      run.guess[run.guess.length - 1].forEach((v, q) => { if (v) correction.add(q); });
      return { bp: run, osd: null, avg: null, correction };
    }
    const L = run.llr.slice(1), n = code.n;
    const avg = new Array(n).fill(0).map((_, q) => L.reduce((a, l) => a + l[q], 0) / L.length);
    const osd = osd0(code, kind, run.syndrome, avg);
    return { bp: run, osd, avg, correction: osd.correction };
  }

  // ------------------------------------------------------------------ result
  /**
   * After decoding: residual = error x correction. Returns
   *   'exact'      the decoder found your error exactly
   *   'stabilizer' different from your error, but the difference is harmless (degeneracy)
   *   'logical'    the difference is a logical operator: a logical error
   * plus which logicals (X-bar_k / Z-bar_k) it flips.
   */
  function verdict(code, residual) {
    const anti = (a, b) => a !== "I" && b !== "I" && a !== b;
    const flips = [];
    code.logicals.forEach((L, k) => {
      ["x", "z"].forEach((t) => {
        let s = 0;
        L[t].forEach(([q, l]) => { if (residual[q] && anti(residual[q], l)) s ^= 1; });
        // a residual that anticommutes with Z-bar_k acts as X-bar_k, and vice versa
        if (s) flips.push((t === "z" ? "X̄" : "Z̄") + k);
      });
    });
    const empty = Object.keys(residual).length === 0;
    return { kind: empty ? "exact" : flips.length ? "logical" : "stabilizer", flips };
  }

  const QDec = { rowsOf, splitError, syndromeOf, lookup, isGraphLike, mwpm, bp, osd0, bposd, verdict };
  if (typeof module !== "undefined" && module.exports) module.exports = QDec;
  else root.QDec = QDec;
})(typeof window !== "undefined" ? window : globalThis);
