# Developer reference: `qdec` part 1 (codes, Tanner graphs, BP)

**What this package is for:** everything a decoder needs, built by hand in plain numpy so that every step can be drawn and checked:
codes from polynomials → Tanner graph on the torus → belief propagation with full history. No QEC library is used inside `qdec`.
`ldpc` appears only in the tests and the notebook, as an independent check.

Noise model in part 1: **code capacity**. Data qubits get X errors independently with probability p, and checks are read perfectly.

| File | Purpose |
|---|---|
| `qdec/gf2.py` | Mod-2 linear algebra: row reduction, rank, null space, "is this in the row space", solve. |
| `qdec/codes.py` | `BBCode`: polynomials → H_X, H_Z, n, k, logical operators, torus positions. Presets `toric_code(L)`, `bb72()`, `gross_code()`. |
| `qdec/tanner.py` | Tanner graph object and torus drawings (whole graph, one check's support, a set of qubits). |
| `qdec/bp.py` | `MinSumBP`: normalised min-sum BP that records every iteration; `logical_failure`. |
| `qdec/bp_viz.py` | Pictures of BP: per-iteration heatmaps, unsatisfied-check curve, qubit traces, GIF. |

---

## `gf2.py`

| Function | Parameters | Returns / notes |
|---|---|---|
| `as_gf2(M)` | any array-like | uint8 copy with entries mod 2 |
| `row_reduce(M)` | 0/1 matrix | `(R, pivots)`: reduced row-echelon form and pivot columns |
| `rank(M)` | 0/1 matrix | rank over GF(2) |
| `nullspace(M)` | 0/1 matrix | basis of {v : Mv = 0} as rows; one vector per free column |
| `in_row_space(M, v)` | matrix, vector | `True` if v is a sum of rows of M (checked as "adding v doesn't raise the rank") |
| `solve(M, b)` | matrix, vector | one x with Mx = b (free variables set to 0), or `None`. Will be reused by OSD. |

## `codes.py`

**Construction** (Bravyi et al. 2024): `x = S_l ⊗ I_m`, `y = I_l ⊗ S_m`, with the cyclic shift `S[a, a−1] = 1`.
Then `H_X = [A | B]` and `H_Z = [Bᵀ | Aᵀ]`. The checks commute for any A, B because AB = BA.

| Item | Parameters | What it does |
|---|---|---|
| `cyclic_shift(size)` | int | the shift matrix S |
| `shift_matrices(l, m)` | ints | the commuting matrices x, y |
| `poly_matrix(terms, l, m)` | list of `(a, b)` = xᵃyᵇ | matrix of the polynomial, mod 2 |
| `poly_str(terms)` | list | `"x^3 + y + y^2"` |
| `BBCode(name, l, m, A, B)` | torus size, two term lists | builds `HX`, `HZ` |
| `.n`, `.k` | – | n = 2lm; k = n − rank H_X − rank H_Z |
| `.commutes()` | – | H_X H_Zᵀ = 0 ? |
| `.summary()` | – | one-line description |
| `.logicals()` | – | `(LX, LZ)`, k × n each, paired so LX·LZᵀ = I (mod 2). Cached. Not necessarily minimum weight. |
| `.layout()` | – | drawing layout. `'toric'` when two A terms and two B terms generate the torus (Bravyi's toric layout), otherwise `'plain'` |
| `.grid_shape()` | – | picture grid size (rows, cols) |
| `.data_positions()` | – | (n, 2) grid positions, L qubits then R qubits |
| `.check_positions(kind)` | `'X'` or `'Z'` | (lm, 2) grid positions |
| `.H(kind)` | `'X'` or `'Z'` | `HX` or `HZ` |
| `toric_code(L)` | int | A = 1 + x, B = 1 + y. [[2L², 2, L]] |
| `bb72()` | – | l = m = 6, A = x³ + y + y², B = y³ + x + x². [[72,12,6]] |
| `gross_code()` | – | l = 12, m = 6, same A, B. [[144,12,12]] |

**Toric layout, in words:** the step "A_p − A_q" becomes "one row down" and the step "B_p − B_q" becomes "one column right". If these two steps reach every
cell exactly once, each check gets four neighbour edges, and the remaining terms become long-range edges. For the gross code that's 4 short + 2 long per check.

## `tanner.py`

| Function | Parameters | What it does |
|---|---|---|
| `tanner_graph(code, kinds="XZ")` | code, which checks | networkx graph; nodes `('q', i)`, `('X', c)`, `('Z', c)` with `pos` |
| `torus_period(code)` | – | wrap size of the picture grid |
| `torus_displacement(a, b, period)` | two positions | shortest step from a to b around the torus |
| `edge_lengths(code, kind)` | – | torus distance of every edge (1 = neighbour) |
| `draw_nodes(code, ax, kinds, size)` | – | scatter qubits and checks |
| `draw_edges(code, ax, kind, checks=None, lw, alpha, mark_ghosts=False)` | – | edges. Grey = short, orange = long. With `mark_ghosts=True` a wrapped edge is drawn dashed to the in-box copy of the qubit. |
| `draw_tanner_torus(code, kind, ax, title)` | – | whole graph for one check type |
| `central_check(code, kind)` | – | index of the check nearest the picture centre |
| `draw_check_support(code, kind, check, ax, title)` | – | one check and the qubits it measures |
| `draw_qubit_set(code, qubits, ax, title, color)` | – | highlight a set of qubits (e.g. a logical) |

## `bp.py`

`MinSumBP(H, error_rate, max_iter=50, alpha=0.625)`

| Parameter | Meaning |
|---|---|
| `H` | parity-check matrix. Use `code.HZ` to decode X errors and `code.HX` to decode Z errors. |
| `error_rate` | p. Sets the prior LLR log((1−p)/p). |
| `max_iter` | stop after this many iterations if not converged |
| `alpha` | min-sum scaling factor (ldpc's `ms_scaling_factor`) |

`.decode(syndrome) → BPResult` with `correction`, `converged`, `iterations`, `llr_history` (row 0 = prior), `guess_history`,
`unsatisfied_history`, and the property `prob_history` = 1/(1+e^LLR).

Schedule: flooding/parallel. All check→qubit messages are computed, then all qubit→check messages. The "min over the other qubits" uses the
smallest and second-smallest trick, so one numpy pass updates every edge.

`logical_failure(error, correction, logicals_other)`: the residual error+correction is a logical error if it overlaps an
opposite-type logical an odd number of times. Use `LZ` when decoding X errors.

**Verified:** identical corrections, convergence flags, iteration counts and final LLRs to `ldpc.BpDecoder` (min-sum, parallel, same α)
on toric 6×6, [[72,12,6]] and gross codes.

## `bp_viz.py`

| Function | What it draws |
|---|---|
| `draw_bp_state(code, check_kind, probs, syndrome, true_error, guess, ax, title, show_colorbar)` | One moment: qubits coloured by P(flipped) (white→red), lit checks gold, true error ×, BP guess blue ring. `check_kind` is the type that produced the syndrome (`'Z'` for X errors). |
| `bp_frames(code, check_kind, result, syndrome, true_error, iterations=None, suptitle)` | Row of frames (default iterations 0, 1, 2, 3, 5, last) plus the unsatisfied-checks curve. Returns the figure. |
| `qubit_traces(result, qubits, ax, labels, title)` | P(flipped) against iteration for chosen qubits |
| `save_bp_gif(code, check_kind, result, syndrome, true_error, path, fps=2)` | Animated GIF, one frame per iteration |

## Known limits (part 1)

* Code capacity only. There is no measurement noise and no circuit, so the Stim circuit and detector error model come in a later part.
* Logical operators come from linear algebra, not a minimum-weight search, so they aren't guaranteed to be the shortest.
* Our BP is about 10–30× slower than `ldpc` (Python vs C++). That's fine for pictures; use `ldpc` for large sweeps.
* Only flooding schedule and min-sum. Serial schedule, product-sum and memory terms (Relay-BP) come later.

## References
Bravyi et al., *Nature* 627, 778 (2024), arXiv:2308.07915 · Tanner, IEEE Trans. Inf. Theory 27, 533 (1981) · Gallager, *Low-Density Parity-Check Codes* (1963) ·
Fossorier, Mihaljević, Imai, IEEE Trans. Commun. 47, 673 (1999) · Chen & Fossorier, IEEE Trans. Commun. 50, 406 (2002) · Poulin & Chung, QIC 8, 987 (2008) ·
Roffe et al., arXiv:2005.07016 · Panteleev & Kalachev, *Quantum* 5, 585 (2021).
