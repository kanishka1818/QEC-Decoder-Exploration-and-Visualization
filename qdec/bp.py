"""
bp.py - belief propagation (normalised min-sum), written from scratch,
that REMEMBERS every iteration so we can draw it.

What this file is FOR
---------------------
Library decoders (ldpc, Relay-BP) hand back only the final answer.  To
*see* why BP works on some syndromes and gets stuck on others, we need the
whole history: each qubit's belief after every round of messages.  This
file computes exactly the same thing a library min-sum BP does, but keeps
that history.

The problem being solved (code-capacity noise)
----------------------------------------------
Each data qubit flipped independently with probability p.  We see the
syndrome s = H e (mod 2) and want a guess e_hat with H e_hat = s that is
likely.  BP works with log-likelihood ratios (LLRs):

    LLR(qubit) = log( P(not flipped) / P(flipped) )
    LLR > 0  -> "probably fine",   LLR < 0 -> "probably flipped"
    P(flipped) = 1 / (1 + exp(LLR))

The algorithm (one "flooding" iteration = all checks, then all qubits)
---------------------------------------------------------------------
    start:         every qubit tells every check it touches: prior = log((1-p)/p)
    check -> qubit: "given what the OTHER qubits told me, and my syndrome
                     bit, here is what I think about you":
                     sign = (-1)^s_c * product of the others' signs
                     size = alpha * (smallest |message| among the others)
    qubit belief:   LLR_v = prior + sum of all incoming check messages
    qubit -> check: LLR_v minus what that same check just said
                    (never echo a check's own message back to it)
    guess:          e_hat_v = 1 if LLR_v < 0
    stop when:      H e_hat = s  (all checks satisfied), or max_iter reached

`alpha` (< 1) shrinks the min-sum messages to correct for min-sum
over-confidence ("normalised min-sum", Chen & Fossorier 2002).  The ldpc
package calls this `ms_scaling_factor`.

References
----------
* Gallager, "Low-Density Parity-Check Codes", MIT Press (1963) - BP on LDPC.
* Fossorier, Mihaljevic, Imai, IEEE Trans. Commun. 47, 673 (1999) - min-sum.
* Chen & Fossorier, IEEE Trans. Commun. 50, 406 (2002) - normalised min-sum.
* Poulin & Chung, Quantum Inf. Comput. 8, 987 (2008) - why BP struggles
  on quantum codes (degeneracy, short loops).
* Roffe et al., Phys. Rev. Research 2, 043423 (2020), arXiv:2005.07016 -
  BP+OSD and the `ldpc` package we compare against.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class BPResult:
    """Everything BP did, iteration by iteration."""

    correction: np.ndarray                 # final guess e_hat (0/1 per qubit)
    converged: bool                        # did e_hat satisfy every check?
    iterations: int                        # how many iterations were run
    llr_history: np.ndarray                # (iterations+1, n) beliefs; row 0 = prior
    guess_history: np.ndarray              # (iterations+1, n) hard decisions
    unsatisfied_history: np.ndarray        # (iterations+1,) number of checks still violated
    prior_llr: float = field(default=0.0)

    @property
    def prob_history(self) -> np.ndarray:
        """P(flipped) for every qubit after every iteration (from the LLRs)."""
        return 1.0 / (1.0 + np.exp(np.clip(self.llr_history, -50, 50)))


class MinSumBP:
    """
    Normalised min-sum BP on a parity-check matrix H (rows = checks, cols = qubits).

    Messages are stored on the edges of the Tanner graph.  We keep them in a
    padded (n_checks, max_check_degree) table so one numpy line updates a
    whole layer of messages at once.
    """

    def __init__(self, H, error_rate: float, max_iter: int = 50, alpha: float = 0.625):
        self.H = (np.asarray(H) % 2).astype(np.uint8)
        self.n_checks, self.n = self.H.shape
        self.p = error_rate
        self.max_iter = max_iter
        self.alpha = alpha
        # --- edge tables -------------------------------------------------
        deg = self.H.sum(1)
        self.dmax = int(deg.max())
        self.var = np.zeros((self.n_checks, self.dmax), dtype=np.int64)   # qubit on each edge slot
        self.mask = np.zeros((self.n_checks, self.dmax), dtype=bool)      # is this slot a real edge?
        for c in range(self.n_checks):
            qs = np.nonzero(self.H[c])[0]
            self.var[c, :len(qs)] = qs
            self.mask[c, :len(qs)] = True

    # ------------------------------------------------------------------
    def _check_to_var(self, m_vc: np.ndarray, syndrome: np.ndarray) -> np.ndarray:
        """
        Min-sum check update for every edge at once.
        For edge (c, v): sign = (-1)^s_c * prod of other signs,
                         size = alpha * min of other |messages|.
        Trick: "min over the others" is the overall min, unless v IS the
        minimum, in which case it is the second-smallest.
        """
        big = 1e9
        mags = np.where(self.mask, np.abs(m_vc), big)          # padding never the min
        signs = np.where(self.mask & (m_vc < 0), -1.0, 1.0)    # padding counts as +
        order = np.argsort(mags, axis=1)
        min1 = mags[np.arange(self.n_checks), order[:, 0]]
        min2 = mags[np.arange(self.n_checks), order[:, 1]]
        is_min = np.zeros_like(self.mask)
        is_min[np.arange(self.n_checks), order[:, 0]] = True
        other_min = np.where(is_min, min2[:, None], min1[:, None])
        total_sign = np.prod(signs, axis=1) * np.where(syndrome == 1, -1.0, 1.0)
        other_sign = total_sign[:, None] * signs               # dividing by own sign = multiplying (+-1)
        return np.where(self.mask, self.alpha * other_sign * other_min, 0.0)

    def decode(self, syndrome) -> BPResult:
        """Run BP on one syndrome and return the full history."""
        s = (np.asarray(syndrome) % 2).astype(np.uint8)
        prior = np.log((1 - self.p) / self.p)
        llr = np.full(self.n, prior)
        m_vc = np.where(self.mask, prior, 0.0)                 # start: qubits send the prior

        guess = (llr < 0).astype(np.uint8)
        llrs, guesses = [llr.copy()], [guess.copy()]
        unsat = [int(np.sum((self.H.astype(np.int64) @ guess) % 2 != s))]
        converged = unsat[0] == 0 and not s.any()

        it = 0
        while not converged and it < self.max_iter:
            it += 1
            m_cv = self._check_to_var(m_vc, s)                 # checks talk
            llr = np.full(self.n, prior)
            np.add.at(llr, self.var[self.mask], m_cv[self.mask])   # qubits add up what they heard
            m_vc = np.where(self.mask, llr[self.var] - m_cv, 0.0)  # reply, minus the echo
            guess = (llr < 0).astype(np.uint8)
            n_bad = int(np.sum((self.H.astype(np.int64) @ guess) % 2 != s))
            llrs.append(llr.copy())
            guesses.append(guess.copy())
            unsat.append(n_bad)
            converged = n_bad == 0

        return BPResult(correction=guess, converged=converged, iterations=it,
                        llr_history=np.array(llrs), guess_history=np.array(guesses),
                        unsatisfied_history=np.array(unsat), prior_llr=prior)


def logical_failure(error, correction, logicals_other) -> bool:
    """
    After correcting, did a LOGICAL error happen?
    residual = error + correction (mod 2).  It passes all checks if BP
    converged; it is harmless only if it is a stabilizer, i.e. if it
    overlaps every opposite-type logical an even number of times.
    `logicals_other` = Z logicals when decoding X errors (and vice versa).
    """
    residual = (np.asarray(error) + np.asarray(correction)) % 2
    return bool(np.any((np.asarray(logicals_other, dtype=np.int64) @ residual) % 2))
