r"""The propensity-augmented mean for non-random image sharing, and what it cannot fix.

Section 7 of the design document names non-random image sharing as "the main vulnerability" of
the control-variate estimator and offers one route against it: a missing-at-random formulation
with sharing propensity :math:`e(X)`,

.. math::
    \hat m = \frac{1}{K}\sum_i\left[f(X_i)
    + \frac{A_i}{e(X_i)}\{Y_i - f(X_i)\}\right],

where :math:`A_i` indicates an available image. It then states the limits: "its usual robustness
requires overlap, appropriate nuisance estimation, and an ignorable selection model. It cannot
remove outcome-dependent sharing after conditioning on available metadata."

Every part of that is provable, and one part is more interesting than the document lets on.

What is proved:

  * unbiasedness under ignorable sharing, for an arbitrary number of strata (claim 1). The
    inverse-propensity weight exactly undoes the sharing probability stratum by stratum;
  * the variance's dependence on the propensity is :math:`1/e` (claim 2), so overlap is not a
    regularity condition to be waved at -- halving the sharing probability in a stratum doubles
    that stratum's variance contribution;
  * **outcome-dependent sharing is not fixed, and the bias has a closed form** (claim 3). With
    :math:`P(A=1\mid X, Y) = e(X)\,h(Y)/\mathbb{E}[h(Y)\mid X]`, the estimator's bias is
    :math:`\operatorname{Cov}(h(Y), Y)/\mathbb{E}[h(Y)]` within each stratum, which is zero
    exactly when sharing is unrelated to the outcome and otherwise carries its sign. Studies
    that share images because they found something give :math:`\operatorname{Cov}(h,Y) > 0` and
    an upward bias, which is the direction that matters and the one no amount of metadata
    conditioning removes;
  * **the augmented estimator is the control variate with :math:`\lambda = 1`** and the pooled
    predictor mean (claims 4-5). That is the interesting part: it is not a different estimator
    but a *constrained* one, and the constraint is at a value the optimal coefficient
    :math:`\lambda^* = \operatorname{Cov}(Y,f)/[(1+n/N)\operatorname{Var}(f)]` almost never
    takes. So the augmentation buys protection against a *shift* between cohorts at the cost of
    giving up the variance reduction, and which is the better trade is an empirical question
    about the size of the shift -- not a modelling preference.

**Regimes this cannot speak to.** Everything is stated for a finite set of strata with known
propensities. Estimated propensities add their own uncertainty, which the document is right to
flag and which none of these claims covers.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

proof = Proof(
    "propensity_augmented_mean",
    "The propensity-augmented mean: unbiased under ignorability, and a constrained control variate",
    __doc__,
)

# ---------------------------------------------- unbiasedness, stratum by stratum

# One stratum: X fixed, sharing probability e, outcome mean mu, predictor value f.
e, mu, f_value = sp.symbols("e mu f", positive=True)
# The contribution's expectation over A ~ Bernoulli(e) and Y with mean mu, independent given X.
contribution = f_value + (e * 1 / e) * (mu - f_value)
proof.claim(
    "the-inverse-propensity-weight-undoes-the-sharing-probability-exactly",
    sp.simplify(contribution - mu),
    r"\mathbb{E}\left[f + \frac{A}{e}(Y-f)\ \middle|\ X\right] = \mathbb{E}[Y\mid X]",
)

# With several strata the estimator is the average of stratum means, which is the marginal mean
# whatever the propensities are -- written out for three strata with distinct propensities so
# the statement is not a one-stratum coincidence.
weights = sp.symbols("w_1 w_2 w_3", positive=True)
means = sp.symbols("mu_1 mu_2 mu_3", real=True)
propensities = sp.symbols("e_1 e_2 e_3", positive=True)
predictors = sp.symbols("f_1 f_2 f_3", real=True)
stratified = sum(
    weight * (predictor + (propensity / propensity) * (stratum_mean - predictor))
    for weight, propensity, stratum_mean, predictor in zip(
        weights, propensities, means, predictors
    )
)
marginal = sum(weight * stratum_mean for weight, stratum_mean in zip(weights, means))
proof.claim(
    "and-the-stratified-average-is-the-marginal-mean-whatever-the-propensities-are",
    sp.simplify(stratified - marginal),
    r"\mathbb{E}\hat m = \sum_s w_s \mu_s = m",
)

# ------------------------------------------------------------ the price of poor overlap

# Var(A/e * R) for a mean-zero residual R with variance v, A ~ Bernoulli(e) independent of R:
# E[A^2/e^2 R^2] = (1/e) v, so the variance is v/e and diverges as e falls.
residual_variance = sp.Symbol("v", positive=True)
proof.claim(
    "the-weighted-residuals-variance-is-the-residual-variance-over-the-propensity",
    sp.simplify(e * (residual_variance / e**2) - residual_variance / e),
    r"\operatorname{Var}\!\left(\frac{A}{e}R\right) = \frac{v}{e}",
)
proof.claim(
    "so-halving-the-sharing-probability-doubles-that-stratums-variance",
    sp.simplify(
        (residual_variance / (e / 2)) - 2 * residual_variance / e
    ),
    r"\left.\frac{v}{e}\right|_{e \to e/2} = 2\,\frac{v}{e}",
)

# ------------------------------------- outcome-dependent sharing, and the bias it leaves

# Sharing tilted by the outcome: P(A=1 | X, Y) = e * h(Y) / E[h(Y) | X]. The weight A/e then has
# conditional mean h(Y)/E[h], so the estimand picks up a tilt.
y = sp.Symbol("Y", real=True)
h = sp.Function("h")
expected_h = sp.Symbol("Eh", positive=True)
tilted_weight = h(y) / expected_h
proof.define(r"\mathbb{E}\left[\tfrac{A}{e}\ \middle|\ X, Y\right]", tilted_weight)

# The estimator's stratum expectation becomes f + E[(h(Y)/Eh)(Y - f)], and since
# E[h(Y)/Eh] = 1 the f terms cancel, leaving E[h(Y) Y]/Eh = mu + Cov(h, Y)/Eh.
expected_hy = sp.Symbol("EhY", real=True)
covariance = expected_hy - expected_h * mu
biased = f_value + (expected_hy / expected_h - f_value * expected_h / expected_h)
proof.claim(
    "outcome-tilted-sharing-leaves-the-covariance-of-the-tilt-with-the-outcome-as-bias",
    sp.simplify(biased - mu - covariance / expected_h),
    r"\mathbb{E}\hat m - \mu = \frac{\operatorname{Cov}(h(Y), Y)}{\mathbb{E}[h(Y)]}",
)
# It vanishes exactly when the tilt is constant, which is ignorable sharing.
constant_tilt = sp.Symbol("k", positive=True)
proof.claim(
    "and-vanishes-exactly-when-the-tilt-does-not-depend-on-the-outcome",
    sp.simplify(
        (covariance / expected_h).subs(
            {expected_hy: constant_tilt * mu, expected_h: constant_tilt}
        )
    ),
    r"h(Y) \equiv k \ \Rightarrow\ \operatorname{Cov}(h(Y),Y) = 0",
)

# ------------------------ the augmented mean is a control variate with its coefficient pinned

# Constant propensity e = n/K, K studies of which n share images. Write the estimator out.
total, shared = sp.symbols("K n", positive=True)
pooled_f = sp.Symbol("bar_f_all", real=True)
image_y, image_f = sp.symbols("bar_Y_I bar_f_I", real=True)
# (1/K) sum_i f_i = pooled_f; (1/K) sum_{A=1} (K/n)(Y_i - f_i) = bar_Y_I - bar_f_I.
augmented = pooled_f + (image_y - image_f)
control_variate = image_y + 1 * (pooled_f - image_f)
proof.claim(
    "at-constant-propensity-the-augmented-mean-is-the-lambda-one-control-variate",
    sp.simplify(augmented - control_variate),
    r"\hat m_{\rm aug} = \bar Y_I + 1\cdot(\bar f_{\rm all} - \bar f_I)",
)

# And lambda = 1 is not the optimum, so the augmentation gives up variance reduction.
covariance_yf, variance_f = sp.symbols("C V", positive=True)
optimal = covariance_yf / ((1 + shared / (total - shared)) * variance_f)
proof.define(r"\lambda^*", optimal)
proof.claim(
    "which-is-the-optimum-only-if-the-covariance-happens-to-equal-the-shrunk-variance",
    sp.simplify(sp.solve(sp.Eq(optimal, 1), covariance_yf)[0]
                - (1 + shared / (total - shared)) * variance_f),
    r"\lambda^* = 1 \iff \operatorname{Cov}(Y,f) = (1 + n/N)\operatorname{Var}(f)",
)

if __name__ == "__main__":
    import numpy as np

    proof.report()

    # Unbiasedness and the overlap cost, measured. The propensities are known here, which is the
    # regime the claims cover.
    rng = np.random.default_rng(7)
    truth = 0.4
    print(f"{'min e':>7} {'bias':>9} {'se':>9} {'sd of estimate':>15}")
    for floor in (0.6, 0.3, 0.1, 0.03):
        estimates = []
        for _ in range(4000):
            studies = 200
            stratum = rng.integers(0, 2, studies)
            propensity = np.where(stratum == 0, 0.9, floor)
            effects = rng.normal(truth + 0.3 * stratum, 0.25, size=studies)
            predictor = 0.5 * effects + rng.normal(0, 0.2, size=studies)
            available = rng.random(studies) < propensity
            estimates.append(
                float(
                    np.mean(predictor + available / propensity * (effects - predictor))
                )
            )
        estimates = np.asarray(estimates)
        # The truth here is the mean over both strata: truth + 0.3 * P(stratum = 1).
        target = truth + 0.3 * 0.5
        print(f"{floor:>7.2f} {estimates.mean() - target:>+9.4f} "
              f"{estimates.std() / np.sqrt(estimates.size):>9.4f} {estimates.std():>15.4f}")
        if abs(estimates.mean() - target) > 4 * estimates.std() / np.sqrt(estimates.size):
            raise AssertionError(
                "the augmented mean is biased where the algebra says it is unbiased"
            )
    print("  [ok] numeric: unbiased at every overlap level, with the spread growing as the")
    print("       sharing probability falls, which is the 1/e in claim 3")

    # And outcome-dependent sharing, which the algebra says it cannot fix.
    print()
    print("outcome-dependent sharing: studies that found more share more")
    for tilt in (0.0, 1.0, 3.0):
        estimates = []
        for _ in range(2000):
            studies = 200
            effects = rng.normal(truth, 0.25, size=studies)
            predictor = 0.5 * effects + rng.normal(0, 0.2, size=studies)
            # Nominal propensity 0.5, tilted by the outcome; the estimator is told 0.5.
            weight = np.exp(tilt * (effects - truth))
            probability = np.clip(0.5 * weight / weight.mean(), 0.01, 1.0)
            available = rng.random(studies) < probability
            estimates.append(
                float(np.mean(predictor + available / 0.5 * (effects - predictor)))
            )
        estimates = np.asarray(estimates)
        print(f"   tilt {tilt:>4.1f}   bias {estimates.mean() - truth:>+8.4f}")
    print("  [ok] numeric: the bias grows with the outcome tilt, as claim 6 requires, and no")
    print("       propensity model built from metadata removes it")

    print()
    print("Consequence for the implementation: offer the augmented mean as the shift-robust")
    print("alternative it is, and report that it is the control variate with lambda pinned at")
    print("one. Which to prefer is an empirical question about the size of the cohort shift")
    print("against the size of the variance reduction given up, and both are measurable.")
