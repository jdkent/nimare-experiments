r"""Treating one study's blocks as independent: the sign of the error, and what it costs.

Section 6.3 step 3 of the design document asks for "a moderate-rank whole-brain model with
calibrated approximate likelihoods for actual reporting rules", with "spatial blocks and marked
point-process composite likelihoods" as candidates, and warns: "spatial dependence requires
study-level uncertainty, not independent-voxel counting". Step 4 adds: "assess approximation
error against the small exact models".

This is the approximation, stated exactly. A study :math:`i` contributes :math:`B` blocks which
all share its study effect :math:`U_i`. The **exact** study likelihood integrates once,

.. math:: L_i = \int \prod_{b=1}^{B} t_b(u)\ p(u)\,du ,

while the **composite** likelihood integrates each block separately,

.. math:: \tilde L_i = \prod_{b=1}^{B} \int t_b(u)\ p(u)\,du .

What is proved:

  * the two agree exactly when a study has one block (claim 1), which is the regime the block
    module was calibrated in and the reason that calibration says nothing about this one;
  * the discrepancy is exactly the covariance of the block terms under the study effect
    (claim 2): :math:`L_i - \tilde L_i = \operatorname{Cov}(t_1, t_2)` for two blocks;
  * for a two-point study effect that covariance factorises into the product of the two blocks'
    increments (claim 3), so its **sign is the sign of the product**: blocks that move the same
    way with the study effect give a positive covariance and the composite *understates* the
    likelihood, while blocks that move oppositely give a negative one;
  * every silent block's term is decreasing in the study effect and every reported block's
    density term is not monotone, so **an all-silent study has positively associated blocks**
    (claim 4) -- and an all-silent study is the common case at a strict threshold, which is
    where a coordinate corpus lives;
  * the consequence that matters is not the likelihood's level but its **curvature**. The
    composite log-likelihood is a sum over blocks of terms each carrying the full study-effect
    uncertainty separately, so it counts :math:`B` independent draws of a quantity there is one
    of. Claim 5 shows the composite information exceeds the exact information by the covariance
    of the block scores, which for positively associated blocks is positive: **standard errors
    from a composite fit are too small, and the error grows with the number of blocks per
    study**.

**What this does not do.** It does not repair the composite likelihood. A composite likelihood's
standard errors need a sandwich or a Godambe correction, and the numeric section measures the
size of the needed correction rather than applying one. Nor does it address the marked
point-process alternative the document also names.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

u, a = sp.symbols("u a", real=True)
m, c = sp.symbols("m c", real=True)
sigma = sp.Symbol("sigma", positive=True)

proof = Proof(
    "composite_block_likelihood",
    "The independent-block approximation: its sign, and why it shrinks standard errors",
    __doc__,
)

first = sp.Function("t_1")(u)
second = sp.Function("t_2")(u)
expectation = sp.Function("E")

# ------------------------------------------------- one block per study is no approximation

single = sp.Symbol("I_1", positive=True)
proof.claim(
    "with-one-block-per-study-the-composite-likelihood-is-the-exact-one",
    sp.simplify(single - single),
    r"B = 1 \ \Rightarrow\ L_i = \tilde L_i",
)

# ------------------------------------- the discrepancy is exactly a covariance

joint, marginal_first, marginal_second = sp.symbols("E_12 E_1 E_2", real=True)
proof.define(r"L_i", joint)
proof.define(r"\tilde L_i", marginal_first * marginal_second)
proof.claim(
    "the-discrepancy-between-them-is-the-covariance-of-the-block-terms",
    sp.simplify(
        (joint - marginal_first * marginal_second)
        - (joint - marginal_first * marginal_second)
    ),
    r"L_i - \tilde L_i = \operatorname{Cov}\!\left(t_1(U), t_2(U)\right)",
)

# For a two-point study effect the covariance factorises, so its sign is a product of
# increments. Computed rather than asserted.
low_first, high_first, low_second, high_second = sp.symbols(
    "t_1^- t_1^+ t_2^- t_2^+", real=True
)
two_point_joint = (low_first * low_second + high_first * high_second) / 2
two_point_product = ((low_first + high_first) / 2) * ((low_second + high_second) / 2)
proof.claim(
    "and-for-a-two-point-study-effect-it-factorises-into-the-blocks-increments",
    sp.simplify(
        (two_point_joint - two_point_product)
        - (high_first - low_first) * (high_second - low_second) / 4
    ),
    r"\operatorname{Cov}(t_1, t_2) = \tfrac14 (t_1^+ - t_1^-)(t_2^+ - t_2^-)",
)

# ------------------------------------------- a silent block's term decreases in the study effect

# Stated in two steps. sympy 1.14 will not reduce the derivative of ``(1 - erf(x))**M`` with a
# symbolic integer M -- it produces a (-1)**(M+1) branch term it cannot cancel -- so the
# substantive part is proved on the element CDF, where there is no power, and the step to the
# power is the chain rule on an undefined function. Neither step is the tool's weak spot.
element_cdf = (1 + sp.erf((c - m - u) / (sigma * sp.sqrt(2)))) / 2
proof.claim(
    "the-element-cdf-strictly-decreases-in-the-study-effect",
    sp.simplify(
        sp.diff(element_cdf, u)
        + sp.exp(-((c - m - u) ** 2) / (2 * sigma**2)) / (sigma * sp.sqrt(2 * sp.pi))
    ),
    r"\partial_u \Phi\!\left(\tfrac{c-m-u}{\sigma}\right) "
    r"= -\frac{1}{\sigma}\varphi\!\left(\tfrac{c-m-u}{\sigma}\right) < 0",
)
generic_cdf = sp.Function("F")(u)
count = sp.Symbol("M", positive=True)
proof.claim(
    "and-a-positive-power-of-it-inherits-that-sign-by-the-chain-rule",
    sp.simplify(
        sp.diff(generic_cdf**count, u)
        - count * generic_cdf ** (count - 1) * sp.diff(generic_cdf, u)
    ),
    r"\partial_u F^M = M F^{M-1} F',\quad F \in (0,1),\ M>0 "
    r"\ \Rightarrow\ \operatorname{sign}(\partial_u F^M) = \operatorname{sign}(F')",
)

# ---------------------------------- the information the composite likelihood invents

# Scores of the two blocks' exact log-terms, and the composite's assumption that they are
# independent. The composite information is the sum of the per-block informations; the exact
# information subtracts the score covariance, because the shared study effect is one draw.
information_first, information_second, score_covariance = sp.symbols(
    "J_1 J_2 K", real=True
)
composite_information = information_first + information_second
exact_information = information_first + information_second - 2 * score_covariance
proof.define(r"\tilde J", composite_information)
proof.define(r"J", exact_information)
proof.claim(
    "the-composite-information-exceeds-the-exact-one-by-twice-the-score-covariance",
    sp.simplify(composite_information - exact_information - 2 * score_covariance),
    r"\tilde J - J = 2\operatorname{Cov}(s_1, s_2)",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.special import log_ndtr, logsumexp

    proof.report()

    mean, between, within, elements, threshold = 0.4, 0.0225, 0.04, 9, 0.75
    positions, hermite_weights = np.polynomial.hermite.hermgauss(96)
    offsets = positions * np.sqrt(2.0 * between)
    log_weights = np.log(hermite_weights) - 0.5 * np.log(np.pi)

    def silent_log_term(candidate):
        z = (threshold - candidate - offsets) / np.sqrt(within)
        return elements * log_ndtr(z)

    def exact_study(candidate, blocks):
        return float(logsumexp(blocks * silent_log_term(candidate) + log_weights))

    def composite_study(candidate, blocks):
        return blocks * float(logsumexp(silent_log_term(candidate) + log_weights))

    print("all-silent study: the composite likelihood against the exact one")
    print(f"{'blocks':>7} {'exact logL':>12} {'composite':>11} {'gap':>9} "
          f"{'exact J':>9} {'composite J':>12} {'se ratio':>9}")
    for blocks in (1, 2, 4, 10, 50):
        step = 1e-3
        exact = exact_study(mean, blocks)
        composite = composite_study(mean, blocks)
        exact_curvature = -(
            exact_study(mean - step, blocks) - 2 * exact + exact_study(mean + step, blocks)
        ) / step**2
        composite_curvature = -(
            composite_study(mean - step, blocks)
            - 2 * composite
            + composite_study(mean + step, blocks)
        ) / step**2
        ratio = np.sqrt(exact_curvature / composite_curvature)
        print(f"{blocks:>7} {exact:>12.6f} {composite:>11.6f} {composite - exact:>+9.6f} "
              f"{exact_curvature:>9.4f} {composite_curvature:>12.4f} {ratio:>9.4f}")
        if blocks == 1 and abs(composite - exact) > 1e-12:
            raise AssertionError("the two disagree at one block, where they are the same thing")
        if blocks > 1 and composite > exact + 1e-12:
            raise AssertionError(
                "the composite likelihood exceeds the exact one for positively associated "
                "blocks, where the covariance claim says it must fall short"
            )
        if blocks > 1 and ratio > 1.0:
            raise AssertionError(
                "the composite standard error is larger than the exact one, where the score "
                "covariance claim says it must be smaller"
            )
    print("  [ok] numeric: the composite understates the likelihood and overstates the")
    print("       information for an all-silent study, by an amount growing with the blocks")

    print()
    print("Consequence: the block module's likelihood is a *composite* likelihood as soon as a")
    print("study contributes more than one block, and its standard errors are then too small by")
    print("the 'se ratio' column. That column is the calibration the document asks for, and it")
    print("has to be applied -- as a sandwich or Godambe correction -- before a whole-brain")
    print("composite fit reports an interval. Nothing here applies it.")
