# User guide

## The explorer

Open `explorer/index.html` in a browser, or the [live copy](explorer/index.html) on this site.

**Layout.** Codes on the left, in four folding groups. The picture in the middle. Details on the right, in folding sections; the
**Details** button hides the whole panel. **Theme** switches between your system's setting, light and dark.

| to do this | do this |
|---|---|
| add an error | click a qubit: X, then Z, then Y, then none |
| read a check | keep **Inspect checks** on and click a face, line or square |
| see the ancilla readings | **Show**, then *Ancilla readings (0 / 1)* |
| see a logical operator | pick one under **Logical** |
| hide a check type, error strings or qubit numbers | **Show** |
| move and zoom | drag, and scroll or pinch |

**Tabs**

- **Lattice:** the code drawn from its geometry. Lit checks turn amber.
- **Tanner graph:** every check joined to its qubits.
- **Check matrix:** $H_X$ and $H_Z$ as tables. The last column is your syndrome. Click a column header to add an error.
- **How it's written:** the stabilizer list, the 0/1 matrices, and the polynomial or algebraic form.
- **Decode:** a decoder corrects your error one step at a time.

**Decode tab**

1. Click qubits, or press **Random error on this code**.
2. Pick a decoder, or leave it on *Best available*. Options that cannot run on this code say why.
3. Press **Next**. The steps go from the syndrome, through the decoder's working, to the correction and the result.

The result is one of three:

- **Exact:** the decoder found your error.
- **Degeneracy:** it applied a different correction, but the difference is a stabilizer. Still corrected.
- **Logical error:** what is left is a logical operator. The decoder failed.

The **Examples** list has ready-made cases. Examples on other codes say so, and open that code.

**Colours.** Red is $X$, blue is $Z$, violet is $Y$ or a mixed check, amber is a lit check, green is the correction.

## The notebooks

Read them in order.

| notebook | about |
|---|---|
| `01_codes_tanner_bp` | bivariate-bicycle codes from polynomials; BP watched round by round |
| `02_code_atlas` | all 16 codes, drawn four ways |
| `03_widgets_explorer` | the explorer inside Jupyter |
| `04_circuits_and_detector_graphs` | noisy circuits, hook errors, graph or hypergraph |
| `05_decoding_pictures` | lookup, matching, BP + OSD, each checked against a reference |
| `06_distance_logicals_written` | distance with Stim, logicals as strings, written forms |
| `07_live_decoding` | the Decode tab's examples in Python |

## The package

```python
from qdec.families import surface_rotated, catalogue
from qdec.views import code_card
from qdec.circuits import stim_distance

code = surface_rotated(5)
print(code.summary())
code_card(code)                               # four views in one figure
print(stim_distance(surface_rotated(21)))     # (21, 'exact (shortest path)')
```

**Your own code from check matrices**

```python
import numpy as np
from qdec.core import StabilizerCode, Geometry

HX = np.array([[1, 1, 1, 1, 0, 0, 0]])        # your matrices
HZ = ...
geo = Geometry(qubit_pos=..., check_pos=..., shape=["polygon"] * (len(HX) + len(HZ)))
mine = StabilizerCode.css("my code", HX, HZ, geo, {"family": "mine"})
```

**A noisy circuit, decoded by our matching**

```python
from qdec.circuits import memory_circuit
from qdec.dem import MatchingGraph
from qdec.decoders import MWPMDecoder

circuit = memory_circuit(surface_rotated(5), rounds=5, p=0.003, basis="Z")
dem = circuit.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
decoder = MWPMDecoder(MatchingGraph(dem))
det, obs = circuit.compile_detector_sampler().sample(1000, separate_observables=True)
predictions = [decoder.decode(s) & 1 for s in det]
```

For codes where matching does not apply, use `DecodingProblem.from_dem(dem)` with `BPOSDDecoder`.

**A Sinter sweep with our decoders**

```python
import sinter
from qdec.decoders import sinter_decoders

stats = sinter.collect(num_workers=4, tasks=[sinter.Task(circuit=circuit, json_metadata={"d": 5})],
                       decoders=["own_mwpm"], custom_decoders=sinter_decoders(), max_shots=2000)
```

Every function is listed in the [API reference](api/index.md).
