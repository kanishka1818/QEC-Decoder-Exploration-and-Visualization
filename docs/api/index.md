# API reference

Generated from the docstrings in `qdec/` (with [mkdocstrings](https://mkdocstrings.github.io/)). For the *why* behind each module, and the
order to read them in, see [Reading the code](../CODE_GUIDE.md); for compact tables, the [developer reference](../DEVELOPER_REFERENCE.md).

| module | what it is for |
|---|---|
| [`qdec.gf2`](gf2.md) | Mod-2 linear algebra: row reduction, rank, null space, solving H·e = s. |
| [`qdec.core`](core.md) | `StabilizerCode` and `Geometry`: the one object every code becomes. |
| [`qdec.families`](families.md) | The 16-code catalogue and the constructors behind it. |
| [`qdec.codes`](codes.md) | Bivariate-bicycle codes from polynomials, and Bravyi's toric layout. |
| [`qdec.views`](views.md) | Matplotlib pictures: lattice, Tanner graph, check matrix, logicals, matching, OSD. |
| [`qdec.circuits`](circuits.md) | Noisy memory circuits for any CSS code; distances with Stim. |
| [`qdec.dem`](dem.md) | Detector error model → `DecodingProblem` (any decoder) and `MatchingGraph` (matching). |
| [`qdec.bp`](bp.md) | Min-sum belief propagation that records every iteration. |
| [`qdec.decoders`](decoders.md) | Our decoders: lookup, MWPM, BP + OSD, and Sinter adapters. |
| [`qdec.export`](export.md) | JSON for the browser explorer, lookup tables, reduced logicals, written forms, the explorer build. |
| [`qdec.tanner`](tanner.md) | Tanner graphs of BB codes drawn on the torus. |
| [`qdec.bp_viz`](bp_viz.md) | Pictures and GIFs of BP running on BB codes. |
| [`qdec.widgets`](widgets.md) | The explorer inside Jupyter (ipywidgets). |
