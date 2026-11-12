# Getting started

Everything you need to install, run and publish this project, in order. The commands are the same on Windows (Anaconda Prompt),
macOS and Linux unless a step says otherwise.

## 0. What you need

| tool | why | required? |
|---|---|---|
| [conda](https://docs.conda.io/) (Anaconda or Miniforge) | builds the Python environment from `environment.yml` in one command | yes |
| a browser | the explorer (`explorer/index.html`) | yes |
| [Node.js](https://nodejs.org/) | the 19 tests in `test_explorer_decoders.py` run the explorer's JavaScript decoders and compare them with the Python ones | optional (those tests skip without it) |
| Playwright + Chromium | 1 test clicks through the explorer in a headless browser | optional (skips without it) |
| `requirements-docs.txt` | building the documentation website | only for the website |

## 1. Get the code and build the environment

```bash
git clone https://github.com/kanishka1818/QEC-Decoder-Exploration-and-Visualization.git
cd QEC-Decoder-Exploration-and-Visualization
conda env create -f environment.yml      # creates the env "qec-decoders" (Python 3.12), takes a few minutes
conda activate qec-decoders
python -m ipykernel install --user --name qec-decoders   # so Jupyter lists this env as a kernel
```

> **Windows:** run these in *Anaconda Prompt*. If your path has spaces, quote it: `cd "C:\Users\me\My Projects\qdec"`.

Updating an existing environment after `environment.yml` changes: `conda env update -f environment.yml --prune`.

## 2. Check that everything works

```bash
python -m pytest -q tests/
```

Expected: **239 passed, 2 skipped** in about a minute. The two skips are by design. Without Node.js or Playwright you get more skips, not failures.

Optional extras for the full suite:

```bash
conda install -c conda-forge nodejs      # JavaScript-vs-Python decoder tests
pip install playwright
playwright install chromium              # the headless-browser test of the Decode tab
```

Useful variations:

```bash
python -m pytest -q tests/test_decoders.py          # one file
python -m pytest -q -k "mwpm"                       # tests whose name contains "mwpm"
python -m pytest -q -x --lf                         # stop at the first failure; rerun only last failures
```

## 3. Open the explorer

Double-click `explorer/index.html` (or drag it into a browser). It is one self-contained file and needs no server and no internet.

The full tour is in the [user guide](USER_GUIDE.md).

## 4. Run the notebooks

```bash
jupyter lab notebooks/
```

Read them in order; each one tells one part of the story and draws a figure at every step.

| notebook | about | run time |
|---|---|---|
| `01_codes_tanner_bp` | BB codes from polynomials, Tanner graph on the torus, BP watched iteration by iteration | ~15 s |
| `02_code_atlas` | all 16 codes drawn four ways, n/k/d computed | ~15 s |
| `03_widgets_explorer` | the explorer inside Jupyter (needs a live kernel to be interactive) | ~5 s |
| `04_circuits_and_detector_graphs` | our memory circuits, hook errors, graph vs hypergraph | ~50 s |
| `05_decoding_pictures` | lookup → MWPM → BP + OSD, each checked against the reference tool | ~45 s |
| `06_distance_logicals_written` | distance with Stim at any size, logicals as strings, the three written forms | ~80 s |
| `07_live_decoding` | the explorer's Decode tab in Python; JavaScript vs Python agreement (uses Node.js if present) | ~10 s |

Running all of them without opening Jupyter:

```bash
cd notebooks
jupyter nbconvert --to notebook --execute --inplace 0*.ipynb
```

## 5. Use the package in your own code

```python
from qdec.families import catalogue, surface_rotated
from qdec.views import code_card
from qdec.circuits import memory_circuit, stim_distance

for code in catalogue():
    print(code.summary())

code = surface_rotated(5)
code_card(code)                                   # one figure: lattice, Tanner, matrix, logicals
print(stim_distance(surface_rotated(21)))         # (21, 'exact (shortest path)')
circuit = memory_circuit(code, rounds=5, p=0.003, basis="Z")   # a stim.Circuit
```

More examples: [user guide](USER_GUIDE.md#the-package).

## 6. Rebuild the explorer after changing code

The explorer is generated: `explorer/template.html` (page) + `explorer/decoders.js` (JavaScript decoders) + data computed by `qdec.export`
→ `explorer/index.html`. After changing a code, the catalogue, the template or `decoders.js`:

```bash
python -c "from qdec.export import build_explorer_html as b; b('explorer/template.html', 'explorer/index.html')"
```

This takes about a minute (Stim's distance searches for the BB codes). Add `use_stim=False` for a 5-second build that uses the published
distances instead. Then re-run `python -m pytest -q tests/test_explorer_decoders.py tests/test_explorer_unit1.py`.

## 7. The documentation website

The `docs/` folder plus `mkdocs.yml` form a website built with MkDocs Material. On top of the Markdown pages it adds:

* the seven **notebooks**, rendered with their figures (not re-executed when the site is built);
* an **API reference** generated from the docstrings in `qdec/`;
* the **explorer** itself, so anyone can click it online.

The notebooks and the explorer live outside `docs/`. A small build hook (`docs_hooks.py`) copies them into `docs/` when the site is built;
those copies are git-ignored.

```bash
pip install -r requirements-docs.txt      # mkdocs-material, mkdocs-jupyter, mkdocstrings-python
mkdocs serve                              # live preview at http://127.0.0.1:8000, reloads on save
mkdocs build --strict                     # the static site in site/; --strict turns warnings into errors
```

### Publishing on GitHub Pages

The workflow `.github/workflows/docs.yml` builds the site on every push to `main` and publishes it. One-time setup:

1. Push the repository to GitHub (public, or a paid plan for private Pages).
2. On GitHub: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
3. Push. The **Actions** tab shows the run; the site appears at `https://kanishka1818.github.io/QEC-Decoder-Exploration-and-Visualization/`.

If your default branch is not `main`, change `branches: [main]` at the top of the workflow.

## 8. Troubleshooting

| symptom | fix |
|---|---|
| `ModuleNotFoundError: qdec` in a notebook | start Jupyter from the repository root or `notebooks/`, and pick the `qec-decoders` kernel. The notebooks add the repository root to `sys.path` themselves. |
| many tests skipped | Node.js or Playwright missing: see step 2 (optional). |
| a notebook "hangs" in `nbconvert` | almost always the Jupyter kernel failing to start; rerun it. All seven run in under 2 minutes each. |
| explorer shows old data | rebuild it (step 6) and hard-refresh the browser (Ctrl+F5). |
| `mkdocs build --strict` fails on a link | a link in a `.md` file points at a file that does not exist; the message names the page. |
