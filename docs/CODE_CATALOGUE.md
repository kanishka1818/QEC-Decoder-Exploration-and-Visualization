# Code catalogue

The theory behind each code in `qdec`: the geometry, the stabilizers, what to look for in the picture, and the source.

**Notation.** $[[n, k, d]]$ means $n$ physical qubits, $k$ logical qubits, distance $d$. A Pauli on $n$ qubits is a pair of bit
vectors $(x \mid z)$. Two Paulis commute exactly when their symplectic product vanishes:

$$
x_1 \cdot z_2 + z_1 \cdot x_2 = 0 \pmod 2 .
$$

A **CSS** code has checks that are all-$X$ or all-$Z$, so it has two matrices $H_X$ and $H_Z$ with $H_X H_Z^{\mathsf T} = 0$.

## Repetition code

$[[d, 1, 1]]$, with distance $d$ against bit flips only.

- **Geometry:** $d$ qubits on a line, one check between each neighbouring pair.
- **Stabilizers:** $Z_i Z_{i+1}$.
- **Look for:** an $X$ flip lights the one or two checks next to it. A $Z$ flip lights nothing. Each error is an edge between two
  checks, or between a check and the end of the line. This is the first matching graph.
- **Source:** Nielsen & Chuang, §10.1.

## Shor code

$[[9, 1, 3]]$

- **Geometry:** a $3 \times 3$ grid; each row is one block.
- **Stabilizers:** six $ZZ$ pairs inside the rows, and two weight-6 $X$ checks, each spanning two neighbouring rows.
- **Look for:** degeneracy. $Z_0$ and $Z_1$ give the same syndrome and differ by the stabilizer $Z_0 Z_1$, so either correction works.
- **Source:** Shor, PRA 52, R2493 (1995).

## Steane code

$[[7, 1, 3]]$, the smallest colour code.

- **Geometry:** a triangle with 7 qubits and three 4-qubit faces.
- **Stabilizers:** each face twice, once as $XXXX$ and once as $ZZZZ$. The faces are the checks of the Hamming $[7,4,3]$ code.
- **Look for:** one face carries both colours. Every single-qubit error has its own syndrome.
- **Source:** Steane, PRL 77, 793 (1996); Bombin & Martin-Delgado, PRL 97, 180501 (2006).

## Five-qubit code

$[[5, 1, 3]]$, not CSS.

- **Geometry:** a ring of five qubits.
- **Stabilizers:** $XZZXI$ and its cyclic shifts; four of them are independent.
- **Look for:** the 15 single-qubit errors plus "no error" give exactly $16 = 2^4$ syndromes. No syndrome is wasted, which is why the
  code is called perfect.
- **Source:** Laflamme et al., PRL 77, 198 (1996); Bennett et al., PRA 54, 3824 (1996).

## Surface code, unrotated

$[[d^2 + (d-1)^2,\ 1,\ d]]$

- **Geometry:** a $(2d-1) \times (2d-1)$ grid. Data qubits where row + column is even; $X$ checks and $Z$ checks on the other sites.
- **Stabilizers:** each check acts on its 4 grid neighbours (3 at the boundary).
- **Look for:** a logical operator is a string from one boundary to the opposite one. A shorter string lights a check at each end.
- **Source:** Kitaev, Ann. Phys. 303, 2 (2003); Dennis et al., J. Math. Phys. 43, 4452 (2002).

## Surface code, rotated

$[[d^2, 1, d]]$

- **Geometry:** $d \times d$ data qubits. Square checks in a chessboard pattern, weight-2 checks on the boundary.
- **Look for:** about half the qubits of the unrotated code at the same $d$.
- **Circuits:** the CNOT order matters. A Z-shaped order for $X$ checks and an N-shaped order for $Z$ checks keeps the circuit
  distance at $d$. Other orders lose up to half of it (hook errors, measured in notebook 04).
- **Source:** Horsman et al., NJP 14, 123011 (2012); Tomita & Svore, PRA 90, 062320 (2014).

## Toric code

$[[2L^2, 2, L]]$

- **Geometry:** the surface code on a torus. In `qdec` it is the bivariate-bicycle code with $A = 1 + x$ and $B = 1 + y$.
- **Look for:** two logical qubits, one for each way of wrapping around. BP gets stuck here when two corrections are equally short.
- **Source:** Kitaev (2003); Dennis et al. (2002).

## Hypergraph-product codes

Built from two classical codes $H_1$ ($r_1 \times n_1$) and $H_2$ ($r_2 \times n_2$):

$$
H_X = \big[\, H_1 \otimes I \ \big|\ I \otimes H_2^{\mathsf T} \,\big], \qquad
H_Z = \big[\, I \otimes H_2 \ \big|\ H_1^{\mathsf T} \otimes I \,\big], \qquad
n = n_1 n_2 + r_1 r_2 .
$$

- **In the catalogue:** two repetition codes give the surface code $[[13,1,3]]$. Two Hamming codes give $[[58,16,3]]$: 16 logical
  qubits in 58.
- **Source:** Tillich & Zémor, IEEE TIT 60, 1193 (2014), arXiv:0903.0566.

## Bivariate-bicycle codes and the gross code

On an $l \times m$ torus, let $x$ and $y$ be the two cyclic shifts. Pick two polynomials $A$ and $B$ in $x$ and $y$:

$$
H_X = [\, A \mid B \,], \qquad H_Z = [\, B^{\mathsf T} \mid A^{\mathsf T} \,] .
$$

The checks commute for any choice, because $AB = BA$.

- **Geometry:** each cell holds two data qubits, one $X$ check and one $Z$ check. In the toric layout, 4 of each check's 6 edges go
  to nearest neighbours and 2 are long-range.

| code | $l \times m$ | $A$ | $B$ |
|---|---|---|---|
| $[[72,12,6]]$ | $6 \times 6$ | $x^3 + y + y^2$ | $y^3 + x + x^2$ |
| $[[90,8,10]]$ | $15 \times 3$ | $x^9 + y + y^2$ | $1 + x^2 + x^7$ |
| $[[108,8,10]]$ | $9 \times 6$ | $x^3 + y + y^2$ | $y^3 + x + x^2$ |
| $[[144,12,12]]$ (gross) | $12 \times 6$ | $x^3 + y + y^2$ | $y^3 + x + x^2$ |
| $[[288,12,18]]$ | $12 \times 12$ | $x^3 + y^2 + y^7$ | $y^3 + x + x^2$ |

- **Decoding:** one circuit fault can light up to 9 detectors, so the problem is a hypergraph and matching does not apply. BP alone
  oscillates on some weight-4 errors. BP + OSD fixes them.
- **Source:** Bravyi et al., Nature 627, 778 (2024), arXiv:2308.07915; Kovalev & Pryadko, PRA 88, 012311 (2013).

## Telling them apart by the picture

| picture | code |
|---|---|
| a line of qubits | repetition |
| separate red and blue faces | CSS: surface, toric, Shor |
| one face with both colours | colour code (Steane) |
| purple checks | non-CSS (five-qubit) |
| chessboard with small boundary checks | rotated surface |
| everything wraps around | toric, bivariate bicycle |
| checks with long spokes | qLDPC: hypergraph product, bivariate bicycle |
