r"""Retention: what a reporting probability does to the censored likelihood, and to ignoring it.

The design document's section 8.1 separates two models of the same data. A coordinate study
whose estimate clears a directional threshold :math:`c` is *retained* -- actually printed in the
table -- only with probability :math:`\rho`. The document's "correct retention model" has
:math:`\rho` in the likelihood; its "ignoring retention" model does not, and its intervals
collapse to 0% coverage at 500 coordinate studies while its RMSE stays middling. That is the
whole point of the experiment: a wrong reporting model becomes more confidently wrong as tables
accumulate.

This file derives the retention likelihood before it is implemented, and derives what ignoring
it costs, so that neither is fitted to a simulation.

With :math:`S(m,\tau^2) = 1 - \Phi\!\left(\frac{c-m}{\sigma}\right)` and
:math:`\sigma^2 = s^2 + \tau^2`:

.. math::
    P(\text{reported}) = \rho S, \qquad P(\text{not reported}) = 1 - \rho S.

Claims 1-4 give the exact scores in :math:`(m, \tau^2)` for both outcomes. Claim 5 identifies
the misspecified model exactly: ignoring retention **is** the :math:`\rho = 1` boundary case, so
it is not a different family but a fixed point on the edge of this one.

Claims 6-7 are the useful part. The misspecified estimating equation has a population root that
satisfies

.. math:: S(\hat m) = \rho\,S(m_0),

so a model with no retention in it sets the *exceedance* probability equal to the *reported*
rate, and therefore reads a retention shortfall as a smaller effect. Solving,
:math:`\hat m = c - \sigma\,\Phi^{-1}(1 - \rho S_0)`, which is strictly below :math:`m_0` for
every :math:`\rho < 1`. The bias is downward, and it is computable.

**Calibration, and it is a sharp one.** At the document's own settings -- :math:`m_0 = .4`,
:math:`\sigma = .25`, :math:`c = .6`, :math:`\rho = .5` -- claim 7 predicts an asymptotic bias of
:math:`-.1122`. The document measures an ignoring-retention RMSE of **.114** at 1 image + 500
coordinate studies, where the variance is negligible and RMSE is essentially |bias|. The
prediction is analytic and was not fitted to that number. Agreement to .002 confirms both the
algebra and my reading of their setup.

Claim 8 is the limit on all of this: from report rates alone, :math:`m` and :math:`\rho` are
**not separately identified**. The Bernoulli information in :math:`(m, \rho)` is an outer
product, so studies sharing a threshold and a precision repeatedly estimate one number. The
cross-determinant is proportional to :math:`\rho(\varphi_1 S_2/\sigma_1 - \varphi_2 S_1/\sigma_2)`,
which vanishes exactly when the two studies' hazards agree. This is why a freely estimated
per-study detection curve can absorb the effect itself, and why retention has to come from
images, from threshold and precision spread, or from external calibration -- not from the tables
it is applied to.

**Regimes this cannot speak to.** A single scalar retention probability, independent across
studies, applied to a one-sided threshold. Real retention depends on cluster extent, on
proximity to other peaks, and on table length, and none of that is here.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

m, c = sp.symbols("m c", real=True)
s2, tau2 = sp.symbols("s^2 tau^2", positive=True)
rho = sp.Symbol("rho", positive=True)

total_sd = sp.sqrt(s2 + tau2)
z = (c - m) / total_sd
standard_pdf = sp.exp(-z**2 / 2) / sp.sqrt(2 * sp.pi)
survival = (1 - sp.erf(z / sp.sqrt(2))) / 2  # 1 - Phi(z)

proof = Proof(
    "retention_in_the_reporting_model",
    "A retention probability in the reporting model, and the cost of leaving it out",
    __doc__,
)

reported = rho * survival
unreported = 1 - rho * survival
proof.define(r"P(\text{reported}) = \rho S", reported)
proof.define(r"P(\text{not reported}) = 1 - \rho S", unreported)

# ------------------------------------------------------------------- the four exact scores

proof.claim(
    "a-reported-study-scores-the-hazard-of-its-own-threshold",
    sp.simplify(sp.diff(sp.log(reported), m) - standard_pdf / (total_sd * survival)),
    r"\partial_m \log(\rho S) = \frac{\varphi(z)}{\sigma S}",
)
proof.claim(
    "an-unreported-one-scores-the-same-hazard-negated-and-weighted-by-retention",
    sp.simplify(
        sp.diff(sp.log(unreported), m) + rho * standard_pdf / (total_sd * unreported)
    ),
    r"\partial_m \log(1-\rho S) = -\frac{\rho\,\varphi(z)}{\sigma\,(1-\rho S)}",
)
proof.claim(
    "and-in-the-heterogeneity-each-carries-its-standardised-threshold-as-a-factor",
    sp.simplify(
        sp.diff(sp.log(reported), tau2) - z * standard_pdf / (2 * (s2 + tau2) * survival)
    ),
    r"\partial_{\tau^2} \log(\rho S) = \frac{z\,\varphi(z)}{2\sigma^2 S}",
)
proof.claim(
    "with-the-unreported-term-again-the-negated-retention-weighted-version",
    sp.simplify(
        sp.diff(sp.log(unreported), tau2)
        + rho * z * standard_pdf / (2 * (s2 + tau2) * unreported)
    ),
    r"\partial_{\tau^2} \log(1-\rho S) = -\frac{\rho\,z\,\varphi(z)}{2\sigma^2(1-\rho S)}",
)

# ------------------------------- ignoring retention is the rho = 1 edge of the same family

proof.claim(
    "ignoring-retention-is-exactly-the-model-that-assumes-every-exceedance-is-printed",
    sp.simplify(unreported.subs(rho, 1) - (1 + sp.erf(z / sp.sqrt(2))) / 2),
    r"\left.1-\rho S\right|_{\rho=1} = \Phi(z)",
)

# -------------------------------------- where the misspecified estimating equation lands

# In the large-coordinate limit a fraction rho*S_0 of studies report. The rho = 1 score, summed
# in expectation at a candidate mean whose survival is S_hat, is proportional to
#   rho*S_0 / S_hat - (1 - rho*S_0) / (1 - S_hat),
# and the root of that is what the misspecified fit converges to.
survival_true, survival_hat = sp.symbols("S_0 S", positive=True)
expected_score = survival_true * rho / survival_hat - (1 - rho * survival_true) / (
    1 - survival_hat
)
proof.define(r"\mathbb{E}\,\partial_m\ell_{\rho=1}", expected_score)
proof.claim(
    "the-misspecified-fit-sets-the-exceedance-probability-equal-to-the-reported-rate",
    sp.simplify(sp.solve(sp.Eq(expected_score, 0), survival_hat)[0] - rho * survival_true),
    r"S(\hat m) = \rho\,S(m_0)",
)

# Inverting the survival function gives the bias in closed form, and it is downward because the
# survival function decreases in its argument while rho < 1 shrinks the right-hand side.
probit = sp.Function("Phi_inv")
bias_solution = c - total_sd * probit(1 - rho * survival_true)
proof.define(r"\hat m", bias_solution)
proof.claim(
    "so-the-shortfall-is-read-as-a-smaller-effect-and-the-bias-grows-as-retention-falls",
    sp.simplify(
        sp.diff(rho * survival_true, rho) - survival_true
    ),
    r"\partial_\rho\left[\rho S_0\right] = S_0 > 0 "
    r"\ \Rightarrow\ S(\hat m)\ \text{increases with}\ \rho,\ \hat m\ \text{with it}",
)

# ------------------------------------------- what report rates alone cannot separate

# A Bernoulli report indicator with success probability q = rho*S has information
# grad(q) grad(q)^T / (q(1-q)), an outer product. Two studies' gradients in (m, rho) are
# (rho*phi/sigma, S); they are parallel exactly when the hazards agree.
sigma_1, sigma_2 = sp.symbols("sigma_1 sigma_2", positive=True)
pdf_1, pdf_2 = sp.symbols("varphi_1 varphi_2", positive=True)
S_1, S_2 = sp.symbols("S_1 S_2", positive=True)
gradient_1 = sp.Matrix([rho * pdf_1 / sigma_1, S_1])
gradient_2 = sp.Matrix([rho * pdf_2 / sigma_2, S_2])
cross = gradient_1[0] * gradient_2[1] - gradient_1[1] * gradient_2[0]
proof.claim(
    "magnitude-and-retention-separate-only-where-two-studies-hazards-differ",
    sp.simplify(cross - rho * (pdf_1 * S_2 / sigma_1 - pdf_2 * S_1 / sigma_2)),
    r"a_{11}a_{22}-a_{12}a_{21} = \rho\left(\frac{\varphi_1 S_2}{\sigma_1}"
    r"-\frac{\varphi_2 S_1}{\sigma_2}\right)",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.stats import norm

    proof.report()

    # Calibration against the document's own measured number, computed analytically from claim 7.
    true_mean, sd, threshold, retention = 0.4, 0.25, 0.6, 0.5
    survival_0 = norm.sf((threshold - true_mean) / sd)
    predicted_mean = threshold - sd * norm.isf(retention * survival_0)
    predicted_bias = predicted_mean - true_mean
    document_rmse = 0.114
    print(f"  calibration: S_0 = {survival_0:.4f}, rho*S_0 = {retention * survival_0:.4f}")
    print(f"               claim 7 predicts m_hat = {predicted_mean:.4f}, "
          f"bias {predicted_bias:.4f}")
    print(f"               the document measures RMSE {document_rmse:.3f} at 1+500, where the")
    print(f"               variance is negligible, so |bias| should equal it")
    if abs(abs(predicted_bias) - document_rmse) > 0.005:
        raise AssertionError(
            f"the analytic bias {predicted_bias:.4f} disagrees with the document's measured "
            f"{document_rmse:.3f}; one of the two is wrong and guessing which is not allowed"
        )
    print(f"  [ok] calibration: analytic bias matches the independently measured value to "
          f"{abs(abs(predicted_bias) - document_rmse):.4f}")

    # The bias is monotone in retention and vanishes only at rho = 1.
    grid = np.linspace(0.05, 1.0, 96)
    biases = threshold - sd * norm.isf(grid * survival_0) - true_mean
    if not (np.all(np.diff(biases) > 0) and abs(biases[-1]) < 1e-9 and np.all(biases[:-1] < 0)):
        raise AssertionError("the bias is not monotone in retention, or does not vanish at rho=1")
    print("  [ok] numeric: the bias is negative for every rho < 1, increasing in rho, and")
    print("       exactly zero at rho = 1")

    print()
    print("Consequence for the implementation: retention enters as one extra parameter with the")
    print("scores above, and leaving it out is not a neutral simplification but a fixed choice")
    print("of rho = 1 whose bias is computable in advance. It cannot be estimated from report")
    print("rates alone (claim 8), so it has to come from images, from spread in thresholds or")
    print("precisions, or from external calibration.")
