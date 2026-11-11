# Architecture plan: decoders from scratch

**What this project is for:** understanding quantum codes and decoders through **geometry**, **check matrices** and **pictures**,
by writing everything ourselves. You should be able to look at a drawing and say "that's a rotated surface code" or "that's a gross code",
and watch a decoder think, step by step. It doubles as the code base for a tutorial-review paper on decoders.

Home: `QEC-Ecosystem/from-scratch/decoders/` (personal decoder series, started 30 Oct 2026). Package name: `qdec`.

---

## 1. Principles

1. **Geometry first.** A code is defined by where its qubits and checks sit. H_X and H_Z are *read off* the picture (`families.py`).
2. **One object for every code.** `StabilizerCode` stores generators in symplectic form (S_x | S_z), so CSS and non-CSS codes
   (the five-qubit code) go through the same pictures, syndromes, logicals and decoders.
3. **One object for every decoding problem.** `DecodingProblem` (H, L, priors) is the contract between "the noise" and "the decoder".
   Code-capacity noise and full circuit noise produce the same object.
4. **We write the maths; libraries only simulate or referee.**

   | layer | written by us | external, and why |
   |---|---|---|
   | GF(2) algebra, codes, geometry, logicals, distances | ✔ numpy | — |
   | Tanner graphs, lattice/matrix/logical pictures | ✔ matplotlib / SVG | — |
   | memory circuits (gate by gate) | ✔ | Stim **runs** them (fast stabilizer simulation) |
   | detector error model → graph / hypergraph | ✔ | Stim **produces** the DEM text |
   | decoders: lookup, MWPM, BP, BP+OSD | ✔ | scipy's MILP for MWPM with more than 20 defects only (to be replaced by our own blossom) |
   | sweeps / threshold plots | adapters ✔ | Sinter **schedules** shots |
   | PyMatching, ldpc | — | **referees in tests only**, never imported by `qdec` decoders |
5. **Every claim is checked by computation**: n, k, d by brute force where feasible, decoders against PyMatching/ldpc shot by shot,
   circuits silent without noise, circuit distance measured.
6. **Small logical units**, each with code + tests + developer doc + notebook.

---

## 2. Layers and data flow

```
   families.py ──► StabilizerCode (core.py) ──► views.py / tanner.py ─────► pictures (notebooks, cards)
   (geometry)        Sx, Sz, Geometry, info          │
        │                  │                         └─► export.py ─► atlas JSON ─► explorer/index.html (browser)
        │                  │                                                  └─► widgets.py (Jupyter)
        │                  ├──► circuits.memory_circuit ─► stim.Circuit ─► Stim ─► DEM
        │                  │                                                     │
        │                  │                                    dem.parse_dem ◄──┘
        │                  │                                     ├─► DecodingProblem (H, L, priors)
        │                  │                                     └─► MatchingGraph (graph-like faults)
        │                  └──► DecodingProblem.code_capacity / MatchingGraph.from_check_matrix
        │                                                        │
        │                                       decoders.py ◄────┘  LookupDecoder · MWPMDecoder · BPOSDDecoder (uses bp.MinSumBP)
        │                                                        │
        └── codes.py (BBCode: polynomials, toric layout)          └─► sinter adapters ─► threshold plots
```

## 3. Module map

| module | purpose (plain words) |
|---|---|
| `gf2.py` | Arithmetic where 1 + 1 = 0: rank, null space, solving equations. Used everywhere. |
| `codes.py` | Bivariate-bicycle codes from two polynomials, plus Bravyi's toric layout (part 1). |
| `core.py` | `StabilizerCode` and `Geometry`: the common shape every code takes. |
| `families.py` | The catalogue: repetition, Shor, Steane, five-qubit, surface (unrotated/rotated), toric, hypergraph products, BB family incl. gross. |
| `tanner.py` | Tanner graphs of BB codes on the torus (part 1). |
| `views.py` | Lattice / Tanner / matrix / logicals / matching / OSD pictures for any code. |
| `circuits.py` | Writes the noisy memory experiment for any CSS code; hook-safe gate order; circuit distance. |
| `dem.py` | Reads Stim's detector error model into our `DecodingProblem` and `MatchingGraph`; graph vs hypergraph. |
| `bp.py` | Min-sum belief propagation that records every iteration (part 1). |
| `bp_viz.py` | BP frames, traces and GIFs on the torus (part 1). |
| `decoders.py` | Our lookup, MWPM and BP+OSD decoders, plus Sinter adapters. |
| `export.py` | JSON for the browser explorer: codes, lookup tables, recorded decoding runs. |
| `widgets.py` | ipywidgets versions of the explorer for Jupyter. |
| `explorer/` | `template.html` (source) → `index.html` (self-contained, data inlined). Vanilla JS + SVG, no JS libraries. |

## 4. Key design decisions (and why)

| # | decision | why |
|---|---|---|
| D1 | Symplectic (x \| z) storage for every generator | One commutation formula covers CSS and non-CSS; the five-qubit code needs it. |
| D2 | Geometry stored with the code, `period` for tori | Pictures, circuit coordinates (`QUBIT_COORDS`, detector coords) and the explorer all read the same positions. |
| D3 | Bravyi's **toric layout** for BB codes | Shows "toric code + 2 long edges per check"; also makes each polynomial term one clash-free CNOT layer. |
| D4 | Memory circuit: all X checks, then all Z checks | Always correct for any CSS code (no ordering conflicts between overlapping X and Z checks). Cost: deeper rounds, so a lower threshold (§6). |
| D5 | Fixed hook-safe order for the rotated surface code | Found by searching all 24 × 24 orders: X checks in Z-shape, Z checks in N-shape keep circuit distance = d (3, 5, 7). The plain order drops X-basis distance to about d/2. |
| D6 | `DecodingProblem` as the single decoder input | Code-capacity and circuit-level decoding share code paths; BP/OSD/lookup don't care where H came from. |
| D7 | MWPM = Dijkstra + exact subset DP (≤ 20 defects) + MILP above | Exact, short and readable; agrees with PyMatching on 100% of tested shots. A from-scratch blossom is on the roadmap. |
| D8 | OSD sorts by BP's **average** belief over iterations | Shots that reach OSD are the ones where BP oscillates, and the last iteration is noisy. Measured: this roughly halves OSD failures and matches/beats `ldpc` BP+OSD-CS. |
| D9 | Explorer is plain HTML/SVG/JS with the data inlined | No JS dependencies, works offline, same JSON as the widgets, and can be published as a page. |
| D10 | PyMatching / ldpc only in tests and comparison cells | Keeps "what's happening inside the decoder" fully ours, while proving correctness. |
| D11 | Code distance from Stim (exact shortest path for graph-like codes, pruned search otherwise) | Removes the brute-force size cap; never `shortest_graphlike_error` on non-graph-like codes (decomposition artifact gave d = 1 on HGP(Hamming)). |
| D12 | Logicals reduced by multiplying in stabilizers before drawing | Same logical operator, but the familiar straight string on surface/toric codes. |
| D13 | The explorer decodes live in JavaScript (`explorer/decoders.js`), tested against the Python decoders in Node | A static page cannot call Python, and recorded replays cannot decode *your* errors. A second implementation is a risk, so the tests demand identical answers (pairs, weights, every BP belief, OSD corrections). |
| D14 | One decoding run = a list of steps, each with one sentence and one picture overlay | "Next / Back" is how a person reads an algorithm; every overlay (graph, paths, beliefs, OSD ranking, correction, residual) reuses the lattice drawing. |

## 5. Verified results (Nov 2026)

| check | result |
|---|---|
| n, k for all 16 catalogue codes | match the literature (computed) |
| d by brute force | Shor 3, Steane 3, five-qubit 3, surface unrot. 3/5, rot. 3/5, toric 4, HGP(rep) 3, HGP(Hamming) 3 |
| BB family k | [[72,12]], [[90,8]], [[108,8]], [[144,12]], [[288,12]]; no found logical lighter than the paper's d |
| circuits | silent without noise for every code and both bases |
| circuit distance, rotated surface | hook-safe: 3, 5, 7 · plain order (X basis): 2, 3, 4 |
| DEM structure | surface: all faults split into ≤ 2-detector pieces · gross: faults up to 9 detectors, Stim cannot split them (hypergraph) |
| our MWPM vs PyMatching | 100% identical predictions (d=3: 2000 shots, d=5: 1500, d=7: 1500, circuit noise) |
| our BP vs ldpc BP | identical corrections, iteration counts and LLRs (900 shots) |
| our BP+OSD-CS vs ldpc BP+OSD-CS (gross, code capacity) | equal or better at p = 0.03–0.06 |
| code distance with Stim | rotated surface exact up to d = 25 in seconds; BB [[72]], [[90]], [[108]], gross: 6, 10, 10, 12 (upper bounds, = paper) |
| explorer decoders (JS) vs qdec (Python) | MWPM pairs, BP beliefs, BP+OSD-0 corrections: 1050/1050 random shots identical (notebook 07) |
| tests | 239 passed, 2 skipped |

## 6. Known limits

* Circuits: X-then-Z rounds are twice as deep as an interleaved schedule, so the threshold crossing is ~0.4% vs ~0.5–0.6% for Stim's interleaved circuits (with a lighter noise model).
* Non-CSS codes (five-qubit) have pictures, logicals, lookup decoding and syndromes, but no memory circuit yet.
* MWPM above 20 lit detectors uses scipy's integer-program solver (exact, but not our code). It's slower than PyMatching by orders of magnitude, which is fine for teaching and not for big sweeps.
* BB distances come from Stim's pruned search: an upper bound that matches the paper ([[288,12,18]] keeps the paper's value, its search takes minutes).
* Colour codes beyond Steane, and subsystem codes (Bacon–Shor), are not in the catalogue yet.

## 7. Roadmap (next units)

1. **Own blossom algorithm** (Edmonds) to replace the MILP fallback; time it against PyMatching.
2. **Union-Find decoder** with cluster-growth pictures on the lattice.
3. **Clustering decoders** (BP-LSD, ambiguity clustering) and **Relay-BP**, drawn like OSD.
4. **BP+OSD at circuit level on the gross code** through Sinter (performance work: sparse elimination).
5. **Interleaved X/Z schedules** and **Bravyi's depth-8 gross-code schedule**; compare circuit distance and error rates.
6. More families: 6.6.6 / 4.8.8 colour codes, Bacon–Shor (subsystem), generalised bicycle, lifted products.
7. Explorer: ~~live decoding on your own errors (unit 2)~~ done 8 Nov 2026; Stim circuit panel with Stim's SVG on demand (unit 3), Union-Find (unit 4) and our own blossom (unit 5) with step pictures; space-time (3D) detector view.
8. Paper: one chapter per decoder (theory + our implementation + pictures + comparison).
