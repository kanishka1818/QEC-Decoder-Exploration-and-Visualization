"""
qdec - quantum codes and decoders written from scratch, to SEE what happens.

    gf2        mod-2 linear algebra
    core       StabilizerCode + Geometry (one shape for every code, CSS or not)
    families   repetition, Shor, Steane, five-qubit, surface (unrotated/rotated), toric,
               hypergraph products, bivariate-bicycle family incl. the gross code
    codes      BB codes from polynomials + Bravyi's toric layout
    views      lattice / Tanner / matrix / logical / matching / OSD pictures
    tanner     Tanner graphs of BB codes on the torus
    circuits   noisy memory circuits for any CSS code (Stim runs them)
    dem        detector error model -> DecodingProblem / MatchingGraph (graph vs hypergraph)
    bp, bp_viz min-sum belief propagation that records every iteration, and its pictures
    decoders   lookup, MWPM, BP+OSD (ours) + Sinter adapters
    export     JSON for the browser explorer (explorer/index.html)
    widgets    the explorer inside Jupyter (ipywidgets)
"""
