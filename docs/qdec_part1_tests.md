# Test-suite reference: `qdec` part 1

Run with `python -m pytest -q tests/` from the repo root. **46 tests**, about 2 seconds. Tests that compare against `ldpc` are skipped if it isn't installed.

## `tests/test_codes.py` (26 tests): are the codes the ones the papers describe?

| Test | What it proves |
|---|---|
| `test_rank_small_examples` | rank on hand-checkable matrices (dependent rows counted once) |
| `test_nullspace_vectors_are_in_kernel_and_complete` | every null-space vector satisfies Mv = 0 and the count obeys rank–nullity (20 random matrices) |
| `test_solve_finds_solution_or_none` | `solve` reproduces b, and returns `None` for x₀ = 0 and x₀ = 1 |
| `test_rank_matches_ldpc` | our GF(2) rank equals `ldpc.mod2.rank` on every code |
| `test_x_and_y_commute_and_are_permutations` | the shift matrices commute and are permutations |
| `test_polynomial_matrix_has_one_entry_per_term` | 3 monomials give rows of weight 3 |
| `test_parameters` (×4) | n, k and check weight: toric 3×3 [[18,2]], toric 6×6 [[72,2]], [[72,12]], gross [[144,12]] |
| `test_checks_commute` (×5) | H_X H_Zᵀ = 0 for every code |
| `test_any_polynomials_give_commuting_checks` | the same holds for 10 random A, B (the construction guarantees it) |
| `test_logicals_are_valid_and_paired` (×5) | logicals pass the opposite checks, aren't stabilizers, and are paired as LX·LZᵀ = I |
| `test_toric_3x3_distance_is_3_by_brute_force` | no X logical of weight 1 or 2 exists, and one of weight 3 does |
| `test_toric_layout_is_local` | every toric-code edge goes to a grid neighbour |
| `test_positions_are_all_distinct` | no two objects share a grid point |
| `test_bb_toric_layout_has_four_short_edges_per_check` (×2) | Bravyi's layout: exactly 4 neighbour edges + 2 long edges per check |

## `tests/test_tanner.py` (8 tests): is the drawn graph the matrix?

| Test | What it proves |
|---|---|
| `test_graph_matches_check_matrices` (×3) | node count, one edge per 1 in H, check degree w, qubit degree w |
| `test_displacement_takes_short_way_round` | torus distances go the short way |
| `test_edge_length_counts` | toric: all length 1; gross: 288 short and 144 long Z edges |
| `test_drawing_runs` (×3) | the drawing functions run for every code |

## `tests/test_bp.py` (12 tests): is our BP a correct BP?

| Test | What it proves |
|---|---|
| `test_zero_syndrome_needs_no_iterations` | an empty syndrome gives an empty correction and 0 iterations |
| `test_every_single_qubit_error_is_corrected` (×3) | all n single-qubit errors are fixed exactly, on all three codes |
| `test_history_is_consistent` | history shapes, probabilities in [0,1], guesses = sign of LLR, converged ⇔ 0 unsatisfied ⇔ syndrome reproduced |
| `test_matches_ldpc_min_sum` (×3) | 150 random shots per code: same correction, convergence, iteration count and final LLRs (to 1e-6) as `ldpc.BpDecoder` |
| `test_toric_tie_is_never_broken` | the notebook's diagonal example: never converges, and the four tied qubits keep identical probabilities below ½ |
| `test_all_toric_ties_fail` | on the 4×4 toric code, every syndrome with two weight-2 explanations defeats BP |
| `test_logical_failure_detects_logicals_not_stabilizers` | residual = stabilizer → not a failure; residual = logical → failure |
| `test_pictures_run` | frames, traces and GIF are produced |
