r"""One study does not become many, however many voxels it has. A bound, not a caution.

Section 6.1 of the design document warns: "One study with thousands of voxels can help identify
a few spatial coefficients under such assumptions. It does not become thousands of independent
studies or reveal between-laboratory variation." That is the right instinct and it is stated as
a caution. It is a theorem, with a number attached, and the number decides whether the
moderate-rank spatial model the same section proposes can be fitted at all.

The structure that makes it true. Write a study's own level as

.. math:: \eta_i = m + U_i, \qquad U_i \sim N(0, \tau^2),

so :math:`\eta_i \sim N(m, \tau^2)`. Every observation the study produces -- every block's
maximum, every silence, at every voxel -- depends on the parameters *only* through
:math:`\eta_i`, because the within-study noise and the reporting rule do not involve
:math:`(m, \tau^2)`. So :math:`\eta_i` is **sufficient**, and the study's likelihood is

.. math:: L_i(m,\tau^2) = \int \underbrace{\prod_b t_b(\eta)}_{\text{parameter-free}}
    \ N(\eta; m, \tau^2)\,d\eta .

What follows:

  * the block terms carry no dependence on the parameters once written in :math:`\eta`
    (claim 1), which is the factorisation above;
  * the study's data are therefore a **noisy measurement of one draw** from
    :math:`N(m,\tau^2)`, and by the data-processing inequality for Fisher information no amount
    of within-study data can carry more information about :math:`(m,\tau^2)` than the draw
    itself. That ceiling is the information of a single normal observation, which is
    :math:`\operatorname{diag}(1/\tau^2,\ 1/(2\tau^4))` (claims 2-4);
  * so the per-study information about the heterogeneity is at most :math:`1/(2\tau^4)`
    **whatever the voxel count**, and :math:`K` studies give at most :math:`K/(2\tau^4)`.
    Inverting (claim 5), a relative standard error :math:`r` on :math:`\tau^2` needs

    .. math:: K \ \ge\ \frac{2}{r^2}

    studies **with images** -- 200 for 10%, 800 for 5%. No spatial model, no basis, no rank
    choice and no quantity of coordinate tables moves that number, because it is a statement
    about how many draws from the between-study distribution exist;
  * and the composite likelihood violates the bound outright (claim 6): treating a study's
    :math:`B` blocks as independent studies makes its information :math:`B` times a single
    block's, which grows without limit. That is a second and independent proof that the
    independent-block approximation is wrong, reached from sufficiency rather than from the
    covariance of the block terms.

**What this does not say.** It bounds information about :math:`(m, \tau^2)` -- the between-study
parameters. Within-study spatial structure is not bounded this way: a single study's thousands of
voxels do inform the *shape* of its own map, which is exactly what section 6.1's basis
coefficients are for. The bound says that shape cannot be converted into knowledge of
between-laboratory variation, which is the conversion the caution warns against.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

eta, m = sp.symbols("eta m", real=True)
tau2 = sp.Symbol("tau^2", positive=True)
threshold, element_sd = sp.symbols("c sigma", positive=True)
count = sp.Symbol("M", positive=True, integer=True)

proof = Proof(
    "information_per_study_is_bounded",
    "Per-study information about the between-study parameters has a ceiling the voxel count cannot lift",
    __doc__,
)

# --------------------------------- the block terms are parameter-free in the study's own level

# A block's silence probability written as a function of the study's level eta rather than of
# its departure U: the parameters have gone.
block_term = ((1 + sp.erf((threshold - eta) / (element_sd * sp.sqrt(2)))) / 2) ** count
proof.define(r"t_b(\eta)", block_term)
proof.claim(
    "a-blocks-term-does-not-depend-on-the-mean-once-written-in-the-studys-own-level",
    sp.simplify(sp.diff(block_term, m)),
    r"\partial_m t_b(\eta) = 0",
)
proof.claim(
    "nor-on-the-heterogeneity",
    sp.simplify(sp.diff(block_term, tau2)),
    r"\partial_{\tau^2} t_b(\eta) = 0",
)

# ------------------------------------- the information in one draw, which is the ceiling

log_density = sp.log(
    sp.exp(-((eta - m) ** 2) / (2 * tau2)) / sp.sqrt(2 * sp.pi * tau2)
)
proof.define(r"\log N(\eta; m, \tau^2)", log_density)

# Expected squared score in the mean. E[(eta - m)^2] = tau^2, so the information is 1/tau^2.
score_mean = sp.diff(log_density, m)
proof.claim(
    "the-score-in-the-mean-is-the-residual-over-the-heterogeneity",
    sp.simplify(score_mean - (eta - m) / tau2),
    r"\partial_m \log N = \frac{\eta - m}{\tau^2}",
)
proof.claim(
    "so-one-draw-carries-exactly-the-reciprocal-heterogeneity-about-the-mean",
    sp.simplify(tau2 / tau2**2 - 1 / tau2),
    r"I_{mm} = \frac{\mathbb{E}(\eta-m)^2}{\tau^4} = \frac{1}{\tau^2}",
)

# And in the heterogeneity: the second derivative's expectation. E[(eta-m)^2] = tau^2 and
# E[(eta-m)^4] = 3 tau^4, so the information is 1/(2 tau^4).
score_variance = sp.diff(log_density, tau2)
proof.claim(
    "the-score-in-the-heterogeneity-is-the-standardised-residual-squared-less-one",
    sp.simplify(score_variance - ((eta - m) ** 2 / tau2 - 1) / (2 * tau2)),
    r"\partial_{\tau^2}\log N = \frac{1}{2\tau^2}\left[\frac{(\eta-m)^2}{\tau^2} - 1\right]",
)
# E[((eta-m)^2/tau^2 - 1)^2] = (3 tau^4 - 2 tau^4 + tau^4)/tau^4 = 2.
proof.claim(
    "whose-expected-square-is-two-so-the-information-is-one-over-twice-the-squared-heterogeneity",
    sp.simplify((3 * tau2**2 - 2 * tau2**2 + tau2**2) / tau2**2 / (4 * tau2**2) - 1 / (2 * tau2**2)),
    r"I_{\tau^2\tau^2} = \frac{2}{4\tau^4} = \frac{1}{2\tau^4}",
)

# ----------------------------------------- the study count the ceiling implies

studies, relative = sp.symbols("K r", positive=True)
standard_error = sp.sqrt(2 * tau2**2 / studies)
proof.define(r"\operatorname{se}(\hat\tau^2)", standard_error)
proof.claim(
    "so-a-relative-precision-on-the-heterogeneity-needs-two-over-its-square-studies",
    sp.simplify(
        sp.solve(sp.Eq(standard_error / tau2, relative), studies)[0] - 2 / relative**2
    ),
    r"\frac{\operatorname{se}(\hat\tau^2)}{\tau^2} = r \iff K = \frac{2}{r^2}",
)

# ------------------------------- the composite likelihood breaks the ceiling

blocks = sp.Symbol("B", positive=True)
single_block_information = sp.Symbol("J_1", positive=True)
# Stated as the reciprocal vanishing. "Equals infinity" cannot be written as an expression
# that reduces to zero -- ``limit(B*J_1, B, oo) - oo`` is ``nan``, not zero, and a claim that
# evaluates to nan is not a claim.
proof.claim(
    "treating-a-studys-blocks-as-studies-makes-its-information-grow-without-limit",
    sp.simplify(sp.limit(1 / (blocks * single_block_information), blocks, sp.oo)),
    r"\lim_{B\to\infty} \frac{1}{B\,J_1} = 0,\ \text{so}\ B J_1 \to \infty "
    r"> \frac{1}{2\tau^4}",
)

if __name__ == "__main__":
    import numpy as np

    proof.report()

    print("studies with images needed for a stated relative precision on the heterogeneity")
    print(f"{'relative se':>12} {'studies':>9}")
    for relative_value in (0.50, 0.25, 0.10, 0.05):
        needed = 2.0 / relative_value**2
        print(f"{relative_value:>12.2f} {needed:>9.0f}")

    # The ceiling, checked numerically. Fisher information is an *expectation* over outcomes,
    # not the curvature at one of them: a first attempt took the observed curvature of an
    # all-silent study's log-likelihood and got -110, a negative "information", because that
    # outcome's likelihood is monotone rather than peaked in the heterogeneity. The sum below is
    # the real thing -- every possible report count, weighted by its probability.
    print()
    print("per-study Fisher information about tau^2, as a study contributes more blocks")
    from scipy.special import log_ndtr, logsumexp
    from scipy.stats import binom

    mean_value, between, within, elements, cut = 0.4, 0.0225, 0.04, 9, 0.75
    ceiling = 1.0 / (2.0 * between**2)

    def outcome_probabilities(candidate_between, block_count, node_count):
        """P(r blocks reported) for r = 0..B, integrating over the study's own level."""
        nodes, quad_weights = np.polynomial.hermite.hermgauss(node_count)
        # numpy's own hermgauss underflows above roughly 400 nodes: at 600 it returns 162 zero
        # weights and warns internally, and the resulting rule produces nan. A convergence check
        # run against such a rule condemns every row. Refuse the rule rather than use it.
        if np.any(quad_weights <= 0):
            raise ValueError(
                f"hermgauss({node_count}) returned {int((quad_weights <= 0).sum())} zero "
                "weights; the rule has underflowed and cannot be used"
            )
        log_quad = np.log(quad_weights) - 0.5 * np.log(np.pi)
        offsets = nodes * np.sqrt(2.0 * candidate_between)
        silent = np.exp(elements * log_ndtr((cut - mean_value - offsets) / np.sqrt(within)))
        report = np.clip(1.0 - silent, 1e-300, 1.0 - 1e-16)
        counts = np.arange(block_count + 1)
        log_binomial = binom.logpmf(counts[None, :], block_count, report[:, None])
        return np.exp(logsumexp(log_binomial + log_quad[:, None], axis=0))

    def information_at(block_count, node_count, step_scale=1e-3):
        step = between * step_scale
        high = outcome_probabilities(between + step, block_count, node_count)
        low = outcome_probabilities(between - step, block_count, node_count)
        centre = outcome_probabilities(between, block_count, node_count)
        slope = (high - low) / (2 * step)
        usable = centre > 1e-10 * centre.max()
        return float(np.sum(slope[usable] ** 2 / centre[usable])), float(centre[usable].sum())

    print(f"{'blocks':>7} {'information':>12} {'ceiling':>9} {'share':>7} "
          f"{'mass kept':>10} {'node check':>11}")
    previous = -1.0
    for block_count in (1, 2, 5, 20, 100, 1000):
        information, retained = information_at(block_count, 200)
        # Two independent refinements. A row is usable only if it survives both: more quadrature
        # nodes, and a smaller differencing step. The second is what catches a large block count,
        # where the outcome distribution has a thousand cells and the tail ones contribute
        # slope-squared-over-probability terms that are pure noise.
        finer, _ = information_at(block_count, 300)
        shorter, _ = information_at(block_count, 200, step_scale=5e-4)
        node_gap = abs(finer - information) / max(abs(information), 1.0)
        step_gap = abs(shorter - information) / max(abs(information), 1.0)
        agreement = max(node_gap, step_gap)
        trustworthy = retained > 1 - 1e-9 and agreement < 0.01
        verdict = f"{agreement:.1e}" if trustworthy else "UNUSABLE"
        print(f"{block_count:>7} {information:>12.2f} {ceiling:>9.2f} "
              f"{information / ceiling:>7.3f} {retained:>10.7f} {verdict:>11}")
        if not trustworthy:
            print(f"          the discrete outcome sum is not reliable at {block_count} blocks "
                  f"(nodes differ by {node_gap:.0%}, step by {step_gap:.0%}); "
                  "not evidence either way")
            continue
        if information > ceiling * (1 + 1e-3):
            raise AssertionError(
                f"a study with {block_count} blocks carries {information:.2f} about the "
                f"heterogeneity, above the ceiling {ceiling:.2f} that sufficiency forbids"
            )
        if information < previous - 1e-6:
            raise AssertionError("information fell as blocks were added")
        previous = information
    print("  [ok] numeric: over every block count the computation can resolve, the information")
    print("       rises and stays under the ceiling one exactly known draw would give")

    print()
    print("Consequence for the moderate-rank spatial model of section 6.1: its spatial")
    print("coefficients can be informed by one study's voxels, and its between-study covariance")
    print("cannot. Fitting the second needs studies, about 200 with images for a tenth relative")
    print("precision on tau^2, and no basis, rank or table count substitutes. Any collection")
    print("with a handful of images must borrow that covariance from outside and propagate the")
    print("borrowing's uncertainty, exactly as the document says -- but the size of what has to")
    print("be borrowed is now a number rather than a caveat.")
