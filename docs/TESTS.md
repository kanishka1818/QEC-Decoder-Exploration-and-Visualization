# Tests

Run `python -m pytest -q tests/` from the folder root: **239 passed, 2 skipped**, in about a minute. The two skips are the repetition
code's missing $X$ checks. The first three files (`test_codes.py`, `test_tanner.py`, `test_bp.py`) are described in
[Part 1 tests](qdec_part1_tests.md).

## `test_families.py` (68): is every code what the papers say?
| test | proves |
|---|---|
| `test_generators_commute` ×16 | all generators commute, for every catalogue code |
| `test_n_and_k` ×16 | n and k match the literature (e.g. HGP(Hamming) [[58,16]], BB [[90,8]]) |
| `test_logicals` ×16 | logicals commute with all generators, pair canonically, commute within X and within Z sets |
| `test_distance_by_brute_force` ×8 | d = 3 (Shor, Steane, 5-qubit, surfaces), 5 (surfaces d=5), 4 (toric 4×4) |
| `test_repetition_protects_against_x_only` | (d_X, d_Z) = (5, 1) |
| `test_five_qubit_is_not_css_and_perfect` | 15 single errors + identity → 16 distinct syndromes |
| `test_steane_faces_are_hamming_checks` | 7 distinct non-zero columns |
| `test_hgp_of_open_repetition_is_unrotated_surface` | HGP(rep4, rep4) = [[25,1,4]] = unrotated surface d=4 |
| `test_hgp_hamming` | [[58,16,3]] |
| `test_bb_family_found_logical_not_lighter_than_published_distance` ×5 | no found logical is lighter than the published d |
| `test_syndrome_of_single_x_on_surface_lights_two_or_one_z_checks` | every single X lights 1–2 Z checks only |
| `test_geometry_shapes_match_generators` | geometry arrays have the right sizes |

## `test_circuits_dem.py` (36): are the circuits correct experiments?
| test | proves |
|---|---|
| `test_noiseless_circuit_is_silent` ×18 | no detector or observable fires without noise (9 codes × 2 bases) |
| `test_counts` ×9 | detectors = mz·r + mx·(r−1) + mz; observables = k; CX count = r × (edges of H_X + H_Z) |
| `test_layers_never_reuse_a_qubit` ×3 | each CNOT layer is a matching |
| `test_gross_layers_are_one_per_polynomial_term` | 6 layers of 72 CNOTs per check type |
| `test_hook_safe_order_keeps_distance_plain_order_loses_it` ×2 | circuit distance = d with the hook-safe order; less with the plain order |
| `test_surface_dem_is_graphlike_after_decomposition` | all pieces have ≤ 2 detectors; nothing dropped |
| `test_gross_dem_is_a_hypergraph` | most faults light ≥ 3 detectors; Stim refuses to decompose |
| `test_decoding_problem_matches_dem` | H/L shapes and priors match the parsed DEM |

## `test_decoders.py` (7): do our decoders give the right answers?
| test | proves |
|---|---|
| `test_lookup_fixes_every_single_error_steane` | 8-row table, every single X error fixed |
| `test_mwpm_agrees_with_pymatching` ×2 | ≥ 99.5% identical predictions at d = 3 and 5 (circuit noise) |
| `test_mwpm_exact_dp_and_ilp_agree` | subset DP and the MILP find matchings of equal weight |
| `test_bposd_reproduces_syndrome_and_beats_or_matches_ldpc` | OSD output always reproduces the syndrome; failures within noise of `ldpc` or fewer |
| `test_lookup_on_non_css_five_qubit_code` | full-Pauli lookup on the five-qubit code: 16 syndromes, all single errors fixed |
| `test_sinter_adapters_run` | `sinter.collect` runs our MWPM |

## `test_views_export.py` (25): do the pictures and the explorer data hold up?
| test | proves |
|---|---|
| `test_every_view_runs_on_every_code` ×16 | lattice, Tanner (both), matrix and logicals draw for every code |
| `test_code_card_and_decoding_pictures` | the matching picture's verdicts (corrected / logical error) and OSD panels |
| `test_export_is_json_and_consistent` | JSON size, every replay's code exists, supports match the Python objects |
| `test_lookup_table_corrects_every_single_qubit_error` ×3 | shipped tables fix every single X/Y/Z (Steane, 5-qubit, rotated d=3) |
| `test_replays_tell_the_intended_stories` | each recorded run has the outcome its caption claims |
| `test_toric_tie_has_two_equal_explanations` | the tie example lights exactly 2 checks |
| `test_explorer_html_builds` | data inlined, full HTML document |
| `test_widget_text_helpers` | error parsing and verdict text |

## `test_explorer_unit1.py` (40): explorer update 1
| test | proves |
|---|---|
| `test_stim_distance_surface_exact_at_any_size` ×4 | Stim gives d = 3, 5, 9, 15 exactly for rotated surface codes |
| `test_stim_distance_agrees_with_brute_force` ×5 | Shor, Steane, unrotated surface d=3/5, toric: Stim = brute force |
| `test_stim_distance_regression_hgp_hamming_is_3_not_1` | guards the Stim decomposition trap (detector-less piece reported as d = 1) |
| `test_stim_distance_bb72_matches_paper` | upper-bound search gives 6 |
| `test_repetition_distance_text` | d = 1 as a quantum code, 5 against bit flips |
| `test_code_capacity_circuit_is_silent_without_errors` | the distance circuit is a valid experiment |
| `test_reduced_logicals_are_still_the_same_logicals` ×5 | weight never grows, still passes every check, same commutation with all logicals |
| `test_surface_logicals_draw_as_minimum_weight_strings` ×3 | rotated/unrotated surface and toric: reduced logicals have weight d and draw as strings |
| `test_gross_logical_is_not_a_string` | scattered supports are drawn as rings, not lines |
| `test_written_forms_exist_and_match_the_matrices` ×16 | three sections for every code; stabilizer list has m lines; the 0/1 text is exactly H_X |
| `test_bb_polynomial_form_names_the_polynomials` | gross: A, B, l, m printed |
| `test_five_qubit_polynomial_form` | a(x) = 1 + x³, b(x) = x + x² |

## `test_explorer_decoders.py` (19): explorer update 2, live decoding
The browser decoders (`explorer/decoders.js`) are run in **Node.js** on random errors and compared with the Python decoders. Skipped without Node.

| test | proves |
|---|---|
| `test_mwpm_matches_python` ×8 (4 codes × X/Z) | rotated d=5, unrotated d=5, toric 4×4, repetition: same matched pairs and same total log-weight (1e-9) as `MWPMDecoder` on `MatchingGraph.from_check_matrix`, 40 random errors each; the JS correction always clears the syndrome. (Repetition × X skips: no X checks.) |
| `test_mwpm_refuses_hypergraph_codes` | Steane: `isGraphLike` is false and `mwpm` returns null instead of a wrong answer |
| `test_bp_beliefs_match_python` ×4 | rotated d=5, toric, BB72, gross: same LLR on every qubit at every iteration as `MinSumBP` (1e-9), same convergence and iteration count, 15 errors each |
| `test_bposd_matches_python` ×4 | toric, BB72, gross, HGP(Hamming): same correction as `BPOSDDecoder(osd_order=0)`, 30 errors each; the gross run includes the explorer's oscillation example so the OSD branch is exercised |
| `test_template_inlines_decoders` | the built page carries decoders.js inline (one self-contained file) |
| `test_decode_tab_in_browser` | headless Chromium: every example in the Decode tab gives its stated decoder and verdict (exact / degeneracy / logical); every one of the 16 codes steps through all steps (BP slider moved) without a JavaScript error; **Random error on this code** gives an error and steps on all 16 codes without changing the code; the theme button cycles light, dark, auto; a folded details section stays folded after a redraw |

To make OSD deterministic in both languages, `BPOSDDecoder._osd` now sorts with `np.argsort(..., kind="stable")` (ties keep qubit order, as JavaScript's sort does).
