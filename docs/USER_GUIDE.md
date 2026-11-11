# User guide: decoders from scratch

**What this folder gives you:** a set of quantum codes you can *look at*, a browser explorer you can *click*, notebooks that tell
the story from codes to decoding, and our own decoders you can run and compare.

## 0. Setup (once)

```bash
cd QEC-Ecosystem/from-scratch/decoders
conda env create -f environment.yml
conda activate qec-decoders
python -m ipykernel install --user --name qec-decoders     # so Jupyter can see the env
python -m pytest -q tests/                                # 239 tests, about 50 s (two need Node.js / Playwright)
```

## 1. The browser explorer (no Python needed)

Open `explorer/index.html` in any browser. It's one file with everything inside.

| you do | you see |
|---|---|
| pick a code on the left | its lattice, n/k/d (and **how d was obtained**: brute force, Stim exact, or Stim upper bound), description, references |
| click a qubit (repeatedly) | its error cycles I → X → Z → Y → I. Your errors appear as chips in the side panel, each with its own × button |
| **Inspect** on (default), then click a face / line / square | that stabilizer highlighted, its Pauli letters on the qubits, its row selected in the matrix. Turn Inspect off to stop face clicks selecting anything |
| **check (ancilla) qubits** | a diamond at every check: the ancilla that measures it, reading 0 or 1 (lit). These readings are the syndrome, the only thing a decoder sees. Our Stim circuits put each ancilla at exactly this position |
| two errors under the same check | an **orange string** joins them *through* that check: the check sees two flips and stays dark, so only the ends of the string light up (amber) |
| no check lights at all | the panel says whether your error is a **stabilizer** (harmless: degeneracy) or a **logical operator** (undetectable failure) |
| small codes (≤ 13 qubits) | the **lookup decoder's** answer for your syndrome, and whether it corrects or causes a logical error. It says so when it picks a *different but equivalent* correction (degeneracy) |
| logical ▾ | X̄ₖ, Z̄ₖ or both. On surface/toric codes they are drawn as **labelled strings** across the patch; on qLDPC codes as rings on scattered qubits |
| **Tanner graph** tab | spokes from every check to its qubits (long-range ones visible on BB codes) |
| **Check matrix** tab | a real table: H_X and H_Z (CSS) or Pauli letters (non-CSS), qubits across, checks down, last column = your syndrome bit, lit rows highlighted; click a row to select a stabilizer, click a column header to cycle that qubit's error |
| **How it's written** tab | the same code three ways: the list of stabilizers, the 0/1 check matrices (copyable), and its **polynomial / algebraic form** (BB and toric: A, B, x, y; five-qubit: a(x), b(x) mod x⁵ − 1; repetition: h(x) = 1 + x; HGP and unrotated surface: the product formula; Steane: Hamming; Shor: concatenation) |
| **Decode** tab | **live decoding of your own errors**, one step at a time (Next / Back / Result, or click a step). Pick a decoder: *lookup table* (≤ 13 qubits), *matching (MWPM)* (codes where every qubit sits in ≤ 2 checks of a type: repetition, surface, toric, HGP(rep)), or *BP + OSD* (any CSS code). Unavailable ones say why. Steps: the syndrome (ancillas reading 1) → the decoder's picture (the matching graph, BP's beliefs with a round slider, or the table row) → its moves (each matched pair, or OSD's ranking) → the correction (green rings with the Pauli applied) → the **result**: what is left after correcting, and whether it is nothing (*exact*), a stabilizer (*degeneracy*, fine) or a logical (*logical error*, with which X̄/Z̄ it flips). X and Z parts of a Y error are decoded separately. Eight **examples** cover success, degeneracy, a too-heavy error, a tie, BP freezing, BP oscillating, and two lookup cases (one non-CSS) |

To rebuild it after changing codes: `python -c "from qdec.export import build_explorer_html as b; b('explorer/template.html','explorer/index.html')"`
(about a minute, mostly the Stim distance searches for the BB codes; pass `use_stim=False` for a fast build with published distances).

## 2. The notebooks (read in order)

| notebook | story |
|---|---|
| `01_codes_tanner_bp` | BB codes from polynomials, Tanner graphs on the torus, our BP watched iteration by iteration |
| `02_code_atlas` | every code in the catalogue, drawn four ways, with n/k/d computed |
| `03_widgets_explorer` | the explorer inside Jupyter (needs a live kernel) |
| `04_circuits_and_detector_graphs` | our memory circuits, hook errors, graph vs hypergraph |
| `05_decoding_pictures` | lookup → MWPM → BP+OSD, each checked against the reference tool |
| `06_distance_logicals_written` | distance with Stim at any size, logicals as strings, the three written forms |
| `07_live_decoding` | the Decode tab's examples in Python, a picture per step; matching's tie and too-heavy failures; BP oscillating and OSD rescuing; the browser decoders checked against Python (1050/1050 identical) |

## 3. Using the package in your own code

```python
from qdec.families import surface_rotated, steane, five_qubit, bb, catalogue
from qdec.views import draw_lattice, code_card

code = surface_rotated(5)
print(code.summary())            # surface (rotated) d=5  [[25,1,5]]  family=surface  CSS=True ...
print(code.n, code.k, code.distance())       # 25 1 5   (brute force)
from qdec.circuits import stim_distance
print(stim_distance(surface_rotated(21)))     # (21, 'exact (shortest path)') in about 2 s
code_card(code)                   # lattice + Tanner + matrix + logicals in one figure
```

**Your own code from check matrices** (any CSS code), with a geometry you choose:

```python
import numpy as np
from qdec.core import StabilizerCode, Geometry
HX = np.array([[1,1,1,1,0,0,0]]); HZ = ...            # your matrices
geo = Geometry(qubit_pos=..., check_pos=..., shape=["polygon"] * (len(HX) + len(HZ)))
mine = StabilizerCode.css("my code", HX, HZ, geo, {"family": "mine"})
```

**Syndromes and logicals**

```python
s = code.syndrome(ex=x_error_vector, ez=z_error_vector)   # 1 = generator anticommutes
LX, LZ = code.logicals()                                  # (k, 2n) symplectic rows, paired
```

**A noisy circuit and its detector graph**

```python
from qdec.circuits import memory_circuit, circuit_distance
from qdec.dem import DecodingProblem, MatchingGraph, fault_degree_histogram
c = memory_circuit(surface_rotated(5), rounds=5, p=0.003, basis="Z")
print(circuit_distance(c))                                  # 5
dem = c.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
graph = MatchingGraph(dem)                                  # for MWPM
problem = DecodingProblem.from_dem(dem)                     # for BP / OSD / lookup
```

**Decoding with our decoders**

```python
from qdec.decoders import MWPMDecoder, BPOSDDecoder, LookupDecoder
det, obs = c.compile_detector_sampler().sample(1000, separate_observables=True)
mwpm = MWPMDecoder(graph)
pred = [mwpm.decode(s) & 1 for s in det]                    # bit k = observable k
```

**Sweeps with Sinter, running our decoders**

```python
import sinter
from qdec.decoders import sinter_decoders
stats = sinter.collect(num_workers=4, tasks=[sinter.Task(circuit=c, json_metadata={"d": 5})],
                       decoders=["own_mwpm"], custom_decoders=sinter_decoders(), max_shots=2000)
```

## 4. Colour key (pictures, explorer, widgets)

red = X stabilizer · blue = Z stabilizer · purple = mixed (non-CSS) · dashed blue outline = Z check on the same face as an X check (colour code) ·
amber = lit check (syndrome) · × = true error · green = matching path / correction (letter = Pauli applied) · thin grey lines = matching graph · white→red = BP's flip probability · dashed blue ring = BP's current guess · numbered squares = OSD's top-8 ranking.
