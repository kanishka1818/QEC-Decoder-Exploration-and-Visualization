# Architecture

How the project is put together, why, what has been checked, and what is missing. For a file-by-file tour see
[Reading the code](CODE_GUIDE.md).

## Principles

1. **Geometry first.** A code is defined by where its qubits and checks sit. $H_X$ and $H_Z$ are read off the picture.
2. **One object for every code.** `StabilizerCode` stores generators as $(S_x \mid S_z)$, so CSS and non-CSS codes share the same
   pictures, syndromes, logicals and decoders.
3. **One object for every decoding problem.** `DecodingProblem` ($H$, $L$, priors) is the same for code-capacity noise and circuit noise.
4. **We write the maths.** Libraries only simulate or referee (table below).
5. **Every claim is checked by running it.**

| part | written here | external |
|---|---|---|
| mod-2 algebra, codes, geometry, logicals | yes | none |
| pictures | yes | matplotlib, SVG |
| memory circuits | yes | Stim runs them |
| detector error model to graph or hypergraph | yes | Stim produces the model |
| decoders: lookup, MWPM, BP, BP + OSD | yes | scipy's integer-program solver for MWPM above 20 lit detectors |
| threshold sweeps | adapters | Sinter schedules the shots |
| PyMatching, `ldpc` | no | tests only |

## Design decisions

| # | decision | why |
|---|---|---|
| D1 | Store every generator as $(x \mid z)$ | One commutation formula covers CSS and non-CSS codes. |
| D2 | Keep the geometry with the code | Pictures, circuit coordinates and the explorer read the same positions. |
| D3 | Toric layout for bivariate-bicycle codes | Shows "toric code plus 2 long edges per check". Each polynomial term becomes one CNOT layer. |
| D4 | Circuit rounds run all $X$ checks, then all $Z$ checks | Correct for any CSS code. The cost is deeper rounds and a lower threshold. |
| D5 | Fixed hook-safe CNOT order on the rotated surface code | Keeps circuit distance $d$. Found by searching all orders. |
| D6 | `DecodingProblem` is the only decoder input | BP, OSD and lookup do not care where $H$ came from. |
| D7 | MWPM is Dijkstra plus exact pairing | Short and readable. Agrees with PyMatching on every tested shot. |
| D8 | OSD ranks by BP's belief averaged over iterations | The shots that reach OSD are the ones where BP oscillates. Averaging roughly halves OSD failures. |
| D9 | The explorer is plain HTML, SVG and JavaScript with its data inside | No dependencies, works offline, one file. |
| D10 | Reference decoders appear only in tests | What happens inside a decoder stays ours; correctness is still proved. |
| D11 | Distance from Stim | Removes the brute-force size limit. `shortest_graphlike_error` is used only where it is exact. |
| D12 | Logicals are reduced before drawing | Same operator, drawn as the familiar straight string. |
| D13 | The explorer decodes in JavaScript, tested against Python | A static page cannot call Python. The tests demand identical answers. |
| D14 | A decoding run is a list of steps | Each step has one sentence and one overlay on the lattice. |
| D15 | Show little, fold the rest | Code groups, detail sections, legend and display options all fold. |

## Checked results

| check | result |
|---|---|
| $n$, $k$ for all 16 codes | match the literature |
| $d$ by brute force | Shor 3, Steane 3, five-qubit 3, surface 3 and 5, toric 4, hypergraph products 3 |
| bivariate-bicycle $k$ | 12, 8, 8, 12, 12; no logical found lighter than the published $d$ |
| noiseless circuits | no detector fires, for every code and both bases |
| circuit distance, rotated surface | hook-safe order: 3, 5, 7. Plain order, $X$ basis: 2, 3, 4 |
| detector error model | surface: faults split into pieces of at most 2 detectors. Gross: up to 9, cannot be split |
| our MWPM against PyMatching | identical on 5000 shots ($d = 3, 5, 7$, circuit noise) |
| our BP against `ldpc` | identical corrections, iteration counts and beliefs (900 shots) |
| our BP + OSD against `ldpc` | equal or better on the gross code, $p$ from 0.03 to 0.06 |
| distance with Stim | rotated surface exact up to $d = 25$; bivariate bicycle 6, 10, 10, 12, as published |
| explorer decoders against Python | identical on 1050 of 1050 random errors |
| tests | 239 passed, 2 skipped |

## Known limits

- Circuit rounds are twice as deep as an interleaved schedule. The threshold crossing is near 0.4%, against 0.5 to 0.6% for Stim's circuits.
- Non-CSS codes have pictures, syndromes, logicals and lookup decoding, but no memory circuit.
- MWPM above 20 lit detectors uses scipy's solver. It is exact but much slower than PyMatching.
- Bivariate-bicycle distances from Stim are upper bounds. They match the published values.
- No colour codes beyond Steane, and no subsystem codes.

## Roadmap

1. Explorer: the Stim circuit beside the lattice.
2. Union-Find decoder with cluster-growth pictures.
3. Our own blossom algorithm, replacing the scipy fallback.
4. Clustering decoders (BP-LSD, ambiguity clustering) and Relay-BP.
5. BP + OSD at circuit level on the gross code.
6. Interleaved schedules, including the depth-8 gross-code schedule of Bravyi et al.
7. More families: larger colour codes, Bacon-Shor, lifted products.
