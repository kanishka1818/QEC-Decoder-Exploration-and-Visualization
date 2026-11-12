# QEC Decoder Exploration and Visualization

This project is my exploration of quantum error-correcting codes and their decoders. It started from one thought: any code can be
written down from its check matrix (for many families, just two polynomials), and that matrix is also the best starting point for
visualizing it. From there I build the codes, their circuits and my own decoders, and draw every step.

![The explorer decoding an error on the rotated surface code](docs/figures/49_explorer_decode_new_layout.png)

**[Open the explorer](https://kanishka1818.github.io/QEC-Decoder-Exploration-and-Visualization/explorer/index.html)** ·
**[Read the docs](https://kanishka1818.github.io/QEC-Decoder-Exploration-and-Visualization/)**

## What is here

- **16 codes** built from their check matrices: repetition, Shor, Steane, five-qubit, surface (rotated and unrotated), toric,
  hypergraph products, and the bivariate-bicycle family up to the gross code $[[144,12,12]]$.
- **Four views of each code:** lattice, Tanner graph, check matrix, and its written form (stabilizers, matrices, polynomials).
- **Decoders written from scratch:** lookup table, matching (MWPM), belief propagation, and BP + OSD.
- **An explorer in the browser:** click qubits to add errors and watch a decoder correct them step by step.
- **Seven notebooks** that go from codes to circuits to decoding, with a figure at each step.

Stim runs the circuits. PyMatching and `ldpc` are used only in the tests, to check my decoders against them.

## Run it

```bash
git clone https://github.com/kanishka1818/QEC-Decoder-Exploration-and-Visualization.git
cd QEC-Decoder-Exploration-and-Visualization
conda env create -f environment.yml
conda activate qec-decoders
python -m pytest -q tests/        # about a minute
jupyter lab notebooks/            # start with 01
```

The explorer needs no install: open `explorer/index.html` in a browser.

```python
from qdec.families import surface_rotated
from qdec.views import code_card

code_card(surface_rotated(5))     # lattice, Tanner graph, check matrix and logicals in one figure
```

## Checked results

| claim | how it was checked |
|---|---|
| My MWPM gives the same answers as PyMatching | every shot, rotated surface code, $d = 3, 5, 7$, circuit noise |
| My BP matches `ldpc`; my BP + OSD is at least as good | gross code, code-capacity noise |
| The Z/N CNOT order avoids hook errors | circuit distance stays $d = 3, 5, 7$; a plain order drops it to 2, 3, 4 |
| The explorer's decoders match the Python ones | 1050 of 1050 random errors |

## Where to read more

| | |
|---|---|
| [Getting started](docs/GETTING_STARTED.md) | install, tests, notebooks, building the docs site |
| [User guide](docs/USER_GUIDE.md) | the explorer and the package |
| [Reading the code](docs/CODE_GUIDE.md) | what each file is for and the order to read them |
| [Code catalogue](docs/CODE_CATALOGUE.md) | the theory of each code family |
| [Architecture](docs/ARCHITECTURE.md) | design decisions, limits, roadmap |

## Next

The Stim circuit beside the lattice, Union-Find decoding, and my own blossom algorithm.

## Licence

MIT, see [LICENSE](LICENSE). References are listed in the [code catalogue](docs/CODE_CATALOGUE.md).
