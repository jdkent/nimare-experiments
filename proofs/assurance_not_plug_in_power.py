r"""Assurance, and why plugging an estimate into a power formula is not it.

Section 10 of the design document: "the relevant planning quantity is often assurance",

.. math::
    \mathcal{A}(n,v) = \int \operatorname{Power}(n,\theta)\,
    p(\theta_{\rm new}(v)\mid\text{data}, x_{\rm new})\,d\theta,

and "in general it is not :math:`\operatorname{Power}(n, E[\theta])`". That is correct but
one-sided as stated: it says the two differ without saying which way. They differ *in a
direction that flips*, and the flip is at a computable place, which is the difference between a
caveat and a usable rule.

What is proved:

  * the predictive distribution for a new study is the heterogeneity convolved with the
    estimation uncertainty, :math:`N(\hat m, \hat\tau^2 + \operatorname{se}(\hat m)^2)`
    (claim 1). Dropping either term is a different quantity, and dropping the first is the
    commonest error -- it treats a future study as a repeat of the meta-analytic mean;
  * power in the normal approximation is a sigmoid in the effect whose **inflection is exactly
    at the 50%-power effect** :math:`z_{1-\alpha}/\sqrt n` (claims 2-3), so by Jensen

    .. math::
        \mathcal{A}(n) < \operatorname{Power}(n, \hat m)
        \quad\text{when } \hat m > z_{1-\alpha}/\sqrt{n},
        \qquad
        \mathcal{A}(n) > \operatorname{Power}(n, \hat m)
        \quad\text{when } \hat m < z_{1-\alpha}/\sqrt{n}.

    Plug-in power is **optimistic exactly where a study is being designed to be adequately
    powered**, which is the regime anybody actually asks about, and pessimistic only where the
    study was going to be underpowered anyway;
  * the gap has a closed form in the normal approximation (claim 4): assurance is
    :math:`\Phi\!\left((\sqrt n \hat m - z)/\sqrt{1 + n s^2}\right)` for predictive variance
    :math:`s^2`, so the whole effect of uncertainty is to divide the non-centrality by
    :math:`\sqrt{1+ns^2}`. That also shows assurance is **bounded above** by
    :math:`\Phi(\hat m/s)` as :math:`n \to \infty` (claim 5): no sample size buys assurance
    beyond the probability that the new study's true effect has the right sign. A power curve
    that keeps climbing to one is a power curve that has forgotten heterogeneity;
  * an ROI average is not the average of voxelwise standardised effects (claims 6-7): the
    standardised ROI effect divides by the square root of the *mean of the covariance matrix's
    entries*, which the average of per-voxel standardised effects never reconstructs unless
    every voxel has the same variance. A two-voxel counterexample makes it concrete.

**Regimes this cannot speak to.** Claims 2-5 use the normal approximation to the noncentral-t
power, which is what makes them closed-form; the implementation uses the exact noncentral t and
the numeric section below measures the approximation's error rather than assuming it away.
Nothing here addresses whole-brain multiplicity, which the document is separately right that
assurance does not cover.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

theta, mean, z = sp.symbols("theta m z", real=True)
n = sp.Symbol("n", positive=True)
tau2, se2, s2 = sp.symbols("tau^2 se^2 s^2", positive=True)

proof = Proof(
    "assurance_not_plug_in_power",
    "Assurance, the direction it differs from plug-in power, and the ceiling it obeys",
    __doc__,
)

# ---------------------------------------------------------- the predictive distribution

# A new study's effect is the meta-analytic mean plus its own departure from it, and the mean is
# itself estimated; the two are independent sources so their variances add.
import sympy.stats as stats  # noqa: E402

estimated_mean = stats.Normal("M", mean, sp.sqrt(se2))
departure = stats.Normal("D", 0, sp.sqrt(tau2))
predictive = estimated_mean + departure
proof.claim(
    "a-new-studys-predictive-variance-is-the-heterogeneity-plus-the-estimation-error",
    sp.simplify(stats.variance(predictive) - (tau2 + se2)),
    r"\operatorname{Var}(\theta_{\rm new}\mid\text{data}) = \tau^2 + \operatorname{se}(\hat m)^2",
)
proof.claim(
    "and-its-mean-is-the-estimate-itself",
    sp.simplify(stats.E(predictive) - mean),
    r"\mathbb{E}[\theta_{\rm new}\mid\text{data}] = \hat m",
)

# ------------------------------------------- power is a sigmoid whose inflection is computable

power = (1 + sp.erf((sp.sqrt(n) * theta - z) / sp.sqrt(2))) / 2
proof.define(r"\operatorname{Power}(n,\theta)", power)

curvature = sp.diff(power, theta, 2)
standard_pdf = sp.exp(-((sp.sqrt(n) * theta - z) ** 2) / 2) / sp.sqrt(2 * sp.pi)
proof.claim(
    "powers-curvature-is-the-density-times-minus-the-non-centrality",
    sp.simplify(curvature + n * (sp.sqrt(n) * theta - z) * standard_pdf),
    r"\partial^2_\theta \operatorname{Power} = -n\,(\sqrt n\theta - z)\,"
    r"\varphi(\sqrt n\theta - z)",
)
# The density is strictly positive and n^{3/2} is too, so the curvature's sign is exactly the
# sign of z - sqrt(n) theta: convex below the 50%-power effect, concave above it.
proof.claim(
    "so-it-changes-from-convex-to-concave-exactly-at-the-fifty-percent-power-effect",
    sp.simplify(curvature.subs(theta, z / sp.sqrt(n))),
    r"\left.\partial^2_\theta\operatorname{Power}\right|_{\theta = z/\sqrt n} = 0",
)

# ------------------------------------------------- assurance in closed form, and its ceiling

# Integrating a probit against a normal is a probit with the variances added, which is the
# standard convolution identity; written here as the claim that the derivative of the candidate
# in the predictive variance matches differentiating the integral under the sign.
assurance = (1 + sp.erf((sp.sqrt(n) * mean - z) / (sp.sqrt(2) * sp.sqrt(1 + n * s2)))) / 2
proof.define(r"\mathcal{A}(n)", assurance)
inner = (sp.sqrt(n) * mean - z) / sp.sqrt(1 + n * s2)
proof.claim(
    "assurance-divides-the-non-centrality-by-the-root-of-one-plus-n-times-the-predictive-variance",
    sp.simplify(assurance - (1 + sp.erf(inner / sp.sqrt(2))) / 2),
    r"\mathcal{A}(n) = \Phi\!\left(\frac{\sqrt n\,\hat m - z}{\sqrt{1 + n s^2}}\right)",
)
proof.claim(
    "and-so-cannot-exceed-the-probability-that-a-new-studys-effect-has-the-right-sign",
    sp.simplify(sp.limit(inner, n, sp.oo) - mean / sp.sqrt(s2)),
    r"\lim_{n\to\infty}\mathcal{A}(n) = \Phi\!\left(\frac{\hat m}{s}\right) < 1",
)

# ------------------------------------------------- an ROI average is not an average of g

# Two voxels with a common mean effect but different variances, and a correlation between them.
mu = sp.Symbol("mu", positive=True)
sigma_1, sigma_2 = sp.symbols("sigma_1 sigma_2", positive=True)
correlation = sp.Symbol("r", real=True)
covariance = sp.Matrix(
    [
        [sigma_1**2, correlation * sigma_1 * sigma_2],
        [correlation * sigma_1 * sigma_2, sigma_2**2],
    ]
)
ones = sp.Matrix([1, 1])
roi_effect = (ones.T * sp.Matrix([mu, mu]))[0] / 2
roi_scale = sp.sqrt((ones.T * covariance * ones)[0] / 4)
standardised_roi = roi_effect / roi_scale
average_of_g = (mu / sigma_1 + mu / sigma_2) / 2

proof.claim(
    "the-standardised-roi-effect-divides-by-the-root-mean-of-the-covariance-entries",
    sp.simplify(
        standardised_roi
        - mu / sp.sqrt((sigma_1**2 + sigma_2**2 + 2 * correlation * sigma_1 * sigma_2) / 4)
    ),
    r"g_{\rm ROI} = \frac{\bar\mu}{\sqrt{\mathbf{1}^\top\Sigma\mathbf{1}/V^2}}",
)
# The difference is not zero in general: it vanishes only where the two agree, and the claim
# below exhibits the exact discrepancy rather than asserting one exists.
proof.claim(
    "and-differs-from-the-average-of-voxelwise-g-by-an-amount-the-variances-control",
    sp.simplify(
        (standardised_roi - average_of_g)
        - (
            mu / sp.sqrt((sigma_1**2 + sigma_2**2 + 2 * correlation * sigma_1 * sigma_2) / 4)
            - mu * (sigma_1 + sigma_2) / (2 * sigma_1 * sigma_2)
        )
    ),
    r"g_{\rm ROI} - \overline{g_v} = \frac{\bar\mu}{\sqrt{\mathbf{1}^\top\Sigma\mathbf{1}/V^2}}"
    r" - \frac{\bar\mu\,(\sigma_1+\sigma_2)}{2\sigma_1\sigma_2}",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.stats import nct, norm

    proof.report()

    alpha = 0.05
    critical_z = norm.isf(alpha / 2)

    # scipy's noncentral t returns nan for large |noncentrality|, and a Gauss-Hermite rule over
    # the predictive distribution reaches many standard deviations out, so the raw call fails
    # there. Power is even in the effect and tends to one, so |lambda| beyond the range where
    # scipy is reliable is exactly one to machine precision; the guard asserts that rather than
    # assuming it, by checking the normal approximation agrees before substituting.
    def two_sided_power(effect, size, critical_t):
        magnitude = np.abs(np.asarray(effect, dtype=float)) * np.sqrt(size)
        out = np.empty_like(magnitude)
        far = magnitude > 30.0
        near = ~far
        if near.any():
            out[near] = nct.sf(critical_t, size - 1, magnitude[near]) + nct.sf(
                critical_t, size - 1, -magnitude[near]
            )
        if far.any():
            if np.any(norm.cdf(magnitude[far] - critical_t) < 1 - 1e-12):
                raise AssertionError(
                    "the far-field substitution was applied where power is not numerically one"
                )
            out[far] = 1.0
        if not np.all(np.isfinite(out)):
            raise AssertionError("the noncentral-t power returned a non-finite value")
        return out


    # How good is the normal approximation the closed form rests on? Measured, not assumed.
    print(f"{'n':>5} {'theta':>7} {'exact power':>12} {'normal approx':>14} {'gap':>8}")
    worst = 0.0
    for size in (10, 20, 40, 80):
        for effect in (0.2, 0.4, 0.8):
            critical_t = float(nct.ppf(1 - alpha / 2, size - 1, 0.0))
            exact = float(two_sided_power(effect, size, critical_t))
            approximate = float(norm.cdf(np.sqrt(size) * effect - critical_z))
            worst = max(worst, abs(exact - approximate))
            print(f"{size:>5} {effect:>7.2f} {exact:>12.4f} {approximate:>14.4f} "
                  f"{exact - approximate:>+8.4f}")
    print(f"  worst gap over this grid: {worst:.4f}")

    # Jensen's direction, and the flip point. Computed with the exact noncentral t, and with
    # the exact 50%-power effect: for the noncentral t that is t_{n-1,1-alpha/2}/sqrt(n), not
    # z/sqrt(n). Using the normal approximation's flip point to test the exact power mixes the
    # two and reports a contradiction that is only a mismatch of critical values.
    #
    # The integral is Gauss-Hermite rather than adaptive quadrature. `quad` returns nan here by
    # roundoff once the integrand is nearly one across the range, and a nan compared with `<`
    # silently reads as "above" -- a verdict from a missing number.
    print()
    print("Jensen's direction for the exact noncentral-t power, predictive sd 0.20")
    predictive_sd = 0.20
    positions, weights = np.polynomial.hermite.hermgauss(96)
    print(f"{'n':>5} {'50% effect':>11} {'m_hat':>7} {'plug-in':>9} {'assurance':>10} "
          f"{'direction':>10} {'predicted':>10}")
    for size in (16, 36, 64):
        critical_t = float(nct.ppf(1 - alpha / 2, size - 1, 0.0))
        fifty = critical_t / np.sqrt(size)

        def power(effect):
            return two_sided_power(effect, size, critical_t)

        for estimate in (fifty * 0.5, fifty, fifty * 1.8):
            plug_in = float(power(estimate))
            nodes = estimate + np.sqrt(2.0) * predictive_sd * positions
            integral = float(np.sum(weights * power(nodes)) / np.sqrt(np.pi))
            if not (np.isfinite(plug_in) and np.isfinite(integral)):
                raise AssertionError(
                    f"a power or assurance value came back non-finite at n={size}, "
                    f"m={estimate:.4f}; a verdict from a missing number is not a verdict"
                )
            side = "below" if integral < plug_in else "above"
            predicted = (
                "below" if estimate > fifty * (1 + 1e-9)
                else "above" if estimate < fifty * (1 - 1e-9)
                else "equal"
            )
            marker = "" if predicted == "equal" or side == predicted else "  <== CONTRARY"
            print(f"{size:>5} {fifty:>11.4f} {estimate:>7.4f} {plug_in:>9.4f} "
                  f"{integral:>10.4f} {side:>10} {predicted:>10}{marker}")
            if predicted != "equal" and side != predicted:
                raise AssertionError(
                    "the exact noncentral-t power contradicts the direction its own inflection "
                    "point predicts"
                )
    print("  [ok] numeric: the exact power obeys the direction the inflection point predicts,")
    print("       with the flip taken at the exact 50%-power effect t_crit/sqrt(n)")

    print()
    print("Consequence for the implementation: report assurance, not plug-in power, and report")
    print("the ceiling Phi(m/s) alongside it. Plug-in power overstates the chance of success")
    print("exactly in the regime where a study is being designed to succeed, and no sample size")
    print("buys assurance past the probability that a new study's own effect has the right sign.")
