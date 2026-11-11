# Code catalogue: theory notes for every code in `qdec`

**What this doc is for:** the reading companion to the atlas (notebook 02, explorer). For each code: the geometry it is built from, how the
stabilizers come out of that geometry, its parameters (computed by `qdec`), what to notice in the pictures, and where it comes from.
Written to become the "codes" chapter of the decoder paper/tutorial.

Notation: [[n, k, d]] = physical qubits, logical qubits, distance. A **CSS** code has generators that are all-X or all-Z, so H_X and H_Z can be listed
separately. Two Paulis commute iff x₁·z₂ + z₁·x₂ = 0 (mod 2), the **symplectic product** every check below uses.

---

## Repetition code [[d, 1, 1]] (d_X = d, d_Z = 1)
* **Geometry:** d qubits on a line; a check between each neighbouring pair.
* **Stabilizers:** Z_i Z_{i+1}. H_Z is the (d−1) × d "two ones per row, shifted" matrix.
* **Notice:** a single X flip lights the 1 or 2 checks next to it, while a Z flip lights nothing. It's a classical code inside a quantum one.
  Every error is an *edge* between two checks, or between a check and the end of the line (the boundary). This is the first matching graph.
* **Ref:** Nielsen & Chuang §10.1.

## Shor code [[9, 1, 3]]
* **Geometry:** 3 × 3 grid, row b = block b.
* **Stabilizers:** six Z-pairs inside the rows (bit-flip repetition per block) and two X⁶ checks spanning two neighbouring rows (phase comparison between blocks).
* **Notice:** **degenerate**. Z₀ and Z₁ give the same syndrome and differ by the stabilizer Z₀Z₁, so either correction works.
  Logicals: the X-type one is X on a whole row, the Z-type one is Z on one qubit of each row. Which of them is called "logical X" depends on how |0_L⟩ is defined: with |0_L⟩ = (|000⟩+|111⟩)^⊗3 the X-row flips |+_L⟩ ↔ |−_L⟩ (Lesson 1 follow-up).
* **Ref:** Shor, PRA 52, R2493 (1995).

## Steane code [[7, 1, 3]], the smallest colour code
* **Geometry:** a triangle with 7 qubits (corners, edge midpoints, centre) and three 4-qubit faces.
* **Stabilizers:** each face twice, as XXXX and as ZZZZ. The faces are the parity checks of the Hamming [7,4,3] code: qubit j sits in a distinct non-empty set of faces.
* **Notice:** the same face carries both colours (drawn as a filled red face with a dashed blue outline). Every single X (or Z) error has its own syndrome.
  Transversal Clifford gates follow from the face structure.
* **Ref:** Steane, PRL 77, 793 (1996); Bombin & Martin-Delgado, PRL 97, 180501 (2006).

## Five-qubit code [[5, 1, 3]], perfect and not CSS
* **Geometry:** a ring (pentagon).
* **Stabilizers:** XZZXI and its cyclic shifts (four are independent). Each mixes X and Z, drawn purple.
* **Notice:** the 15 single-qubit Paulis plus identity give exactly 16 = 2⁴ syndromes, the quantum Hamming bound with equality ("perfect").
  It is the smallest code correcting any single-qubit error. Our symplectic core treats it exactly like the CSS codes.
* **Ref:** Laflamme et al., PRL 77, 198 (1996); Bennett et al., PRA 54, 3824 (1996).

## Surface code, unrotated [[d² + (d−1)², 1, d]]
* **Geometry:** (2d−1) × (2d−1) grid; data on sites with r + c even, X checks on (even, odd) "vertices", Z checks on (odd, even) "faces".
* **Stabilizers:** each check = its 4 grid neighbours (3 on the boundary).
* **Notice:** a logical operator is a string of X (or Z) crossing the patch from one boundary to the opposite one; the X-type and Z-type strings run in perpendicular directions. Any shorter string leaves two lit checks at its ends.
  It equals the hypergraph product of two open repetition codes (checked).
* **Ref:** Kitaev, Ann. Phys. 303, 2 (2003); Bravyi & Kitaev, quant-ph/9811052; Dennis et al., J. Math. Phys. 43, 4452 (2002).

## Surface code, rotated [[d², 1, d]]
* **Geometry:** d × d data qubits; plaquettes on a chessboard (X if i + j even), weight-2 "half-moon" plaquettes on the boundary: X type on top and bottom, Z type on left and right.
* **Notice:** about half the qubits of the unrotated code at the same d. For circuits, the CNOT order matters. **Z-shape for X checks and N-shape for Z checks**
  keeps circuit distance = d, while other orders lose up to half of it (hook errors; measured in notebook 04).
* **Ref:** Horsman et al., NJP 14, 123011 (2012); Tomita & Svore, PRA 90, 062320 (2014); Google Quantum AI, Nature 614, 676 (2023).

## Toric code [[2L², 2, L]]
* **Geometry:** the surface code on a torus, with opposite edges glued. In `qdec` it is the BB code with A = 1 + x, B = 1 + y.
* **Notice:** two logical qubits (two independent ways to wrap around). The product of all X checks is the identity, so one check is redundant (m = 2L² listed, rank 2L² − 2).
  BP gets stuck on **diagonal ties** here: two equally short corrections (notebook 01, explorer replay).
* **Ref:** Kitaev 2003; Dennis et al. 2002.

## Hypergraph-product codes (qLDPC family 1)
* **Construction:** classical codes H₁ (r₁ × n₁), H₂ (r₂ × n₂):
  H_X = [H₁ ⊗ I | I ⊗ H₂ᵀ], H_Z = [I ⊗ H₂ | H₁ᵀ ⊗ I]; n = n₁n₂ + r₁r₂.
* **Geometry:** a grid with bits at even and checks at odd coordinates on both axes. Checks reach along their row and column, which is local only if the classical codes are.
* **In the catalogue:** HGP(rep₃, rep₃) = [[13,1,3]] (the surface code) and HGP(Hamming, Hamming) = [[58,16,3]] (16 logical qubits for 58 physical, vs 1 for 13).
* **Ref:** Tillich & Zémor, IEEE TIT 60, 1193 (2014), arXiv:0903.0566.

## Bivariate-bicycle codes, including the gross code (qLDPC family 2)
* **Construction:** on an l × m torus, x and y are the two cyclic shifts; A, B are three-term polynomials; H_X = [A | B], H_Z = [Bᵀ | Aᵀ].
  Because AB = BA, the checks commute for *any* A, B.
* **Geometry:** each cell holds L and R data qubits, an X check and a Z check. In Bravyi's **toric layout** 4 of each check's 6 arrows are nearest-neighbour and 2 are long-range:
  "the toric code plus two long wires per check". Each polynomial term is one full layer of CNOTs.
* **The family (Bravyi et al. Table 3; k computed, d quoted):**

| code | l × m | A | B | n/k |
|---|---|---|---|---|
| [[72,12,6]] | 6 × 6 | x³+y+y² | y³+x+x² | 6 |
| [[90,8,10]] | 15 × 3 | x⁹+y+y² | 1+x²+x⁷ | 11 |
| [[108,8,10]] | 9 × 6 | x³+y+y² | y³+x+x² | 14 |
| **[[144,12,12]] gross** | 12 × 6 | x³+y+y² | y³+x+x² | 12 |
| [[288,12,18]] | 12 × 12 | x³+y²+y⁷ | y³+x+x² | 24 |

* **Notice for decoding:** a single circuit fault lights up to 9 detectors, so the decoding problem is a **hypergraph** and matching doesn't apply.
  Plain BP oscillates on some weight-4 errors (short loops + degeneracy). BP+OSD fixes them (notebooks 01, 04, 05).
* **Ref:** Bravyi et al., Nature 627, 778 (2024), arXiv:2308.07915; Kovalev & Pryadko, PRA 88, 012311 (2013).

---

### How to tell them apart at a glance
| picture | code |
|---|---|
| line of qubits | repetition |
| separate red and blue faces | CSS (surface, toric, Shor, HGP, BB) |
| one face, both colours | colour code (Steane) |
| purple checks | non-CSS (five-qubit) |
| chessboard + half-moons | rotated surface |
| everything wraps, no edges | toric / BB |
| stars with long spokes | qLDPC (HGP, BB, gross) |
