# Developer reference: `qdec` part 2 (codes, views, circuits, DEM, decoders, explorer)

Part 1 modules (`gf2`, `codes`, `tanner`, `bp`, `bp_viz`) are documented in `qdec_part1_reference.md`.
Conventions: 0/1 arrays are `uint8`; a Pauli on n qubits is a length-2n row (x | z); geometry uses (x, y) with **y growing downward**.

---

## `core.py`: the common code object

**`Geometry(qubit_pos, check_pos, shape, period=None)`**
| field | meaning |
|---|---|
| `qubit_pos` (n, 2) | where each data qubit is drawn |
| `check_pos` (m, 2) | label spot of each generator (face centre, vertex, ancilla) |
| `shape` list[str] | `'polygon'` (face through its qubits) or `'star'` (node + spokes, for long-range checks) |
| `period` (Px, Py) or None | picture wraps (torus) |
| `.nearest_image(origin, point)` | copy of `point` closest to `origin` on the torus |

**`StabilizerCode(name, Sx, Sz, geometry, info)`**, or **`StabilizerCode.css(name, HX, HZ, geometry, info)`**

| member | returns / meaning |
|---|---|
| `n`, `m`, `k` | qubits, listed generators (may be dependent), logical qubits = n − rank(S) |
| `S` | (m, 2n) `[Sx \| Sz]` |
| `kinds` | `'X'`, `'Z'` or `'mixed'` per generator |
| `is_css`, `HX`, `HZ` | CSS test; rows of the X-type / Z-type generators |
| `index_of(kind)` | generator indices of that kind (HX row r = generator `index_of('X')[r]`) |
| `commutes()` | all generators commute |
| `generator_string(i)`, `support(i)` | `'XZZXI'`-style string; qubits it touches |
| `syndrome(ex, ez)` | which generators the error X^ex Z^ez anticommutes with |
| `logicals()` | `(LX, LZ)`, (k, 2n) each, canonically paired. CSS: X-type / Z-type. Non-CSS: symplectic Gram–Schmidt. Cached. |
| `distance(max_weight=12, max_candidates=3e6)` | brute-force minimum logical weight, or `None` if the search budget runs out. CSS: X and Z parts searched separately and vectorised. |
| `css_distances(...)` | `(d_X, d_Z)` |
| `params()`, `summary()` | `'[[n,k,d]]'` with `d` = `info['d_known']`; one-line description |

`info` keys used across the package: `family`, `description`, `references`, `d_known`, and optionally `dX`, `dZ`, `cnot_order`, `A`, `B`, `l`, `m`.

Helpers: `symplectic_product(a, b, n)` (commutation table, 1 = anticommute) and `pauli_string(x, z)`.

## `families.py`: the catalogue

| function | parameters | builds |
|---|---|---|
| `repetition(d=5)` | length | Z_iZ_{i+1} checks on a line; (d_X, d_Z) = (d, 1) |
| `shor()` | — | [[9,1,3]]: rows = blocks |
| `steane()` | — | [[7,1,3]] triangle, 3 faces × (X and Z) |
| `five_qubit()` | — | [[5,1,3]], XZZXI cyclic, non-CSS |
| `surface_unrotated(d=3)` | distance | [[d²+(d−1)², 1, d]] on a (2d−1)² grid |
| `surface_rotated(d=3)` | distance | [[d², 1, d]], chessboard + weight-2 boundaries; carries the hook-safe `cnot_order` |
| `toric(L=4)` | size | [[2L², 2, L]] via BB with A = 1+x, B = 1+y |
| `classical_repetition_H(d, closed=False)`, `hamming_H()` | — | classical parity-check matrices |
| `hypergraph_product(H1, H2, name, info)` | two classical H | Tillich–Zémor product, grid geometry, `'star'` shapes |
| `bb(name='gross')` | `'bb72'`, `'bb90'`, `'bb108'`, `'gross'`, `'bb288'` | Bravyi et al. Table 3, toric layout |
| `BB_TABLE` | dict | `(l, m, A terms, B terms, d)` per name |
| `catalogue()` | — | the 16 codes used everywhere |

## `views.py`: pictures (matplotlib)

| function | key parameters | draws |
|---|---|---|
| `draw_lattice(code, ax, show, highlight, error_x, error_z, labels, title, alpha)` | `show` = `'both'`/`'X'`/`'Z'`; `highlight` = generator index or list | faces/stars coloured by kind, highlighted generator with Pauli letters, error and lit checks |
| `draw_tanner(code, ax, show, title)` | — | Tanner graph on the geometry |
| `draw_tanner_layers(code, ax, title)` | — | checks in a row on top, qubits below |
| `draw_check_matrix(code, ax, highlight_row, title)` | — | symplectic matrix as colours |
| `draw_logicals(code, ax, pair, title)` | logical pair index | X̄ (red) and Z̄ (blue) supports |
| `code_card(code, highlight)` | — | the four views + description in one figure |
| `draw_matching(code, error_x, ax, p, title)` | X-error vector | code-capacity MWPM on the lattice; returns `(ax, logical_error)` |
| `draw_osd_steps(code, error_x, decoder, axs)` | a `BPOSDDecoder` on the code's H_Z | BP average belief → OSD pivot qubits → answer vs truth; returns `(axs, logical_error)` |

## `circuits.py`: memory experiments (we write, Stim runs)

| function | parameters | notes |
|---|---|---|
| `cnot_layers(code, kind, order='auto')` | `'X'`/`'Z'`; `'auto'` uses `info['cnot_order']` if present, `'plain'` sorts by direction | list of layers of `(row in H_kind, data qubit)`, no qubit or check twice per layer. On BB codes: one layer per polynomial term. |
| `memory_circuit(code, rounds=3, p=0.0, basis='Z', order='auto')` | CSS code only | qubits: data, then X ancillas, then Z ancillas; noise as in the module docstring; detectors carry (x, y, round) |
| `circuit_stats(circuit)` | — | CX count, ticks, detectors, observables, qubits |
| `circuit_distance(circuit, max_degree=4)` | — | fewest faults flipping an observable silently (Stim's searches) |
| `code_capacity_circuit(code, basis='Z', p=0.01)` | — | one perfect round of checks + independent data flips (X for basis Z, Z for basis X); its smallest undetectable logical = the code distance of that type |
| `stim_distance(code, max_degree=7)` | CSS code | `(d, how)`. If every single flip lights ≤ 2 checks: Stim's shortest path, **exact** at any size. Otherwise Stim's pruned search, an **upper bound**. Never uses `shortest_graphlike_error` on non-graph-like codes (it can return a detector-less decomposition piece: HGP(Hamming) gave d = 1). A missing check type gives d = 1 for that error type |

## `dem.py`: detector error model → our objects

| item | parameters | meaning |
|---|---|---|
| `Fault(p, detectors, observables, parts)` | — | one fault; `parts` = Stim's `^` pieces if decomposed |
| `parse_dem(dem)` | `stim.DetectorErrorModel` | `(faults, num_detectors, num_observables, coords)` |
| `DecodingProblem(H, L, priors, coords)` | — | decoder input; `.from_dem(dem)`, `.code_capacity(H, L, p)` |
| `fault_degree_histogram(dem, use_parts=False)` | — | Counter: detectors lit → number of faults |
| `MatchingGraph(dem)` | needs a graph-like (decomposed) DEM | nodes = detectors + boundary; edges merged; `weight = log((1−p)/p)`; `.dropped` counts hyperedges it had to skip |
| `MatchingGraph.from_check_matrix(H, L, p)` | code capacity | one edge per qubit; `.edge_fault[(u, v)]` = qubit |
| `.path_faults(path)` | node path | faults (qubits) along a matched path |
| `draw_detector_graph_3d(dem, ax, max_faults, title)` | 3D axes | detectors at (x, y, round); graph edges and hyperedge stars |

## `decoders.py`: our decoders

| class | constructor | `decode(syndrome)` returns | notes |
|---|---|---|---|
| `LookupDecoder` | `(problem, max_weight=2)` | observable flips (0/1 array) | enumerates fault sets by weight; `len()` = table rows |
| `MWPMDecoder` | `(graph, max_exact=20)` | observable flips as an **int bitmask** | Dijkstra (heapq) + subset DP; above `max_exact` defects, scipy MILP (`ilp_calls` counts them); `last_matching`, `last_paths` for pictures |
| `BPOSDDecoder` | `(problem, max_iter=50, osd_order=0, alpha=0.625)` | observable flips; `decode_error()` gives the fault vector | BP (per-fault priors supported) → OSD-0 / OSD-CS (all single flips + pairs among the first `osd_order` non-pivots); ranks by **average** BP LLR; `last` holds BP result, order and pivots |
| `SinterOwnMWPM`, `SinterOwnBPOSD(osd_order=4)` | — | — | `sinter.Decoder` subclasses; `sinter_decoders()` → `{'own_mwpm', 'own_bposd'}` |

## `export.py`: data for the explorer

| function | does |
|---|---|
| `lookup_table(code, max_weight=3)` | syndrome string → lightest Pauli correction `{qubit: letter}` |
| `distance_info(code, use_stim=True)` | `{'d', 'how'}` for display: brute force (n ≤ 41), else Stim, else the published value |
| `stim_distance_xz(code)` | `(d against X, d against Z)` for graph-like CSS codes (used for the repetition code's text) |
| `reduce_logical(L, stabs, n)` | multiply in stabilizers while the weight drops: same logical, lighter, drawn as a straight string on surface codes |
| `logical_path(code, support)` | nearest-neighbour ordering of a logical's qubits + `string` flag (all steps ≤ 2.9 grid units) |
| `written_forms(code)` | three sections for the explorer: stabilizer list, 0/1 matrices, polynomial/algebraic recipe |
| `code_to_dict(code, with_lookup=None, dist=None)` | JSON-ready dict (lookup included for n ≤ 13); logicals are reduced and carry `x_path`/`z_path`; includes `written`, `d`, `d_how` |
| `bp_replay(code, error_qubits, p, title, caption, max_iter)` | recorded BP run on X errors |
| `mwpm_replay(code, error_qubits, title, caption, p)` | recorded matching (pairs, qubit paths, verdict) |
| `toric_tie(code)` | a weight-2 error with two equally short corrections |
| `default_replays()`, `build_atlas_data(use_stim=True)`, `write_atlas_json(path)` | the explorer's data (`use_stim=False` skips the slow BB searches) |
| `build_explorer_html(template, out, standalone=True, use_stim=True)` | inline data into `explorer/template.html` |

## `widgets.py`: Jupyter

`parse_error(text, n)` reads `'X3 Z7 Y12'` into `(ex, ez)`. `describe_error(code, ex, ez, table)` gives the text verdict.
`code_explorer(codes=None)`, `bp_player(replay, codes=None)` and `mwpm_player(code=None)` build interactive panels (ipywidgets).

## `explorer/`

`template.html` is the page source. `decoders.js` holds the JavaScript decoders. `export.build_explorer_html` pastes the data
(`/*__ATLAS_DATA__*/null`) and `decoders.js` (`/*__DECODERS_JS__*/`) into the template to make `index.html`, one self-contained file.

### `decoders.js` (global `QDec`; `module.exports` in Node)

Functions take the explorer's code dict (`gens[i] = {kind, support, letters}`, `logicals`, `lookup`). `kind` is the check type that
sees the error: $Z$ checks see $X$ errors.

| function | returns | Python counterpart |
|---|---|---|
| `rowsOf(code, kind)` | generator indices of that kind | `StabilizerCode.index_of` |
| `splitError(err)` | `{x: Set, z: Set}`; a Y goes in both | none |
| `syndromeOf(code, rows, set)` | 0/1 per row | `H @ e % 2` |
| `lookup(code, err)` | `{key, rows, correction}`; `null` correction if the syndrome is not in the table | `export.lookup_table` |
| `isGraphLike(code, kind)` | true if every qubit is in at most 2 checks of that kind | the test in `circuits.stim_distance` |
| `mwpm(code, set, kind, p=0.05)` | `{lit, pairs, correction, cost, hops, greedy, edges}`, or `null` if not graph-like | `MatchingGraph.from_check_matrix` + `MWPMDecoder` |
| `bp(code, set, kind, p, maxIter, alpha)` | `{syndrome, converged, iterations, llr[], guess[], unsat[], probs[]}`; index 0 is the prior | `bp.MinSumBP` |
| `osd0(code, kind, s, llr)` | `{correction, pivots, order}` | `BPOSDDecoder._osd`, order 0 |
| `bposd(...)` | `{bp, osd, avg, correction}`; OSD runs only if BP did not converge | `BPOSDDecoder.decode_error` |
| `verdict(code, residual)` | `{kind: 'exact' or 'stabilizer' or 'logical', flips}` | `isLogicalError` |

`mwpm` merges parallel edges, weights them $\log\frac{1-p}{p}$, runs Dijkstra, and pairs exactly for up to 18 lit checks (greedy
beyond that, flagged).

### `template.html`

| piece | what it does |
|---|---|
| `state` | current code, view, errors, selected check, decoder run, theme-independent UI state (`open`: which sections are folded open) |
| `buildList()`, `GROUPS`, `groupOf()` | the code list as four folding groups; `selectCode` opens only the current group |
| `clickQubit(q)` | cycles a qubit's error: I, X, Z, Y, I |
| `render()` | redraws the current view, the legend and the details panel |
| `drawBoard()` | the SVG: checks, error strings, lit checks, ancillas, logicals, qubits, and the current decoding step's overlay |
| `drawLegend()` | the legend on the picture; lists only what is on screen |
| `drawMatrix()`, `drawWritten()` | the Check matrix and How it's written tabs |
| `renderInspector()`, `collapsify(box)` | the details panel; every `div.section` becomes a folding section |
| `buildRun()` | runs the chosen decoder and returns `{dec, steps}`; each step is `{title, text, ov}` |
| `decodePanel()`, `bindDecodePanel()` | the Decode tab's side panel: error, decoder, step, examples |
| `randomError()` | a random error on the current code |
| `decoderOptions()`, `chosenDecoder()` | which decoders can run on this code; "best available" is lookup, then MWPM, then BP + OSD |
| `PRESETS`, `loadPreset(i)` | the eight examples; their results are asserted in the browser test |
| `syncPlayer()`, `drawSpark()` | the BP round slider and the unexplained-checks curve, shown only on a BP step |
| `setTheme(t)` | auto, light or dark, through `data-theme` on `<html>` |

**Step overlays** (`ov`, read by `drawBoard`): `focus` (draw only that check type), `graph` (matching edges), `pairs` and `pairNow`
(matched paths), `bp` (beliefs and guess at `state.iter`), `avgProbs` and `ranked` (OSD), `corr` (correction rings), `residual`
(the last step draws error times correction).
