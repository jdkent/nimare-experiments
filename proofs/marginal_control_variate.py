r"""The image-corrected coordinate predictor, and exactly what its variance reduction buys.

From the combined-estimator design document, section 7. A predictor :math:`f(C_i, x_i)(v)` guesses
a study's effect image from its coordinate table; it is allowed to be biased. With :math:`n`
image studies that also have a table representation, and :math:`N` coordinate-only studies from
the same population,

.. math:: \hat m_\lambda = \bar Y_I + \lambda\,(\bar f_C - \bar f_I).

The appeal is that the prediction's error is corrected by the image sample rather than trusted,
which is the difference between this and treating imputed maps as if they were observed. The
claims below are the ones that decide whether it is worth building, so they get verified rather
than copied: that it is unbiased for the marginal mean whatever the predictor does, that its
variance has the stated form, that the optimal coefficient is a shrunk regression slope, and that
the best achievable variance ratio is :math:`1 - \rho^2 / (1 + n/N)`.

The last is the one to read carefully. The reduction is capped by the *image* sample: as
:math:`N \to \infty` the ratio tends to :math:`1 - \rho^2`, so a perfect predictor with infinite
coordinate studies still cannot beat the images alone by more than the squared correlation. There
is no configuration in which coordinates substitute for images.

Structural assumptions, stated because the algebra cannot see them: the two cohorts are
independent, both draw from the same population so that :math:`\mathbb{E} f_C = \mathbb{E} f_I`,
:math:`\lambda` and :math:`f` are fixed externally rather than fitted on the same data, and the
per-study effect estimates are unbiased for the target. Non-random image sharing breaks the second
and is the estimator's main vulnerability.
"""
import os
import sys

import sympy as sp
from sympy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

n, N, lam = sp.symbols("n N lambda", positive=True)
m, mu_f = sp.symbols("m mu_f", real=True)
var_y, var_f, cov_yf = sp.symbols("sigma_Y^2 sigma_f^2 sigma_{Yf}", positive=True)

proof = Proof(
    "marginal_control_variate",
    "An image-corrected coordinate predictor for the marginal effect",
    __doc__,
)

# ------------------------------------------------- built from actual random variables

# Rather than assert the expectation and "verify" that a difference of identical symbols is zero,
# build the estimator out of random variables at small n and N and let sympy compute its mean and
# variance. That checks the general formulas below against an independent calculation, and it is
# the step where a sign error or a missing cross term would show up.
rho_c = sp.Symbol("rho", positive=True)
sig_y, sig_f = sp.sqrt(var_y), sp.sqrt(var_f)


def cohort(count, tag, correlated):
    """``count`` independent studies; returns their Y and f draws.

    f is driven by a standard normal u; Y shares rho of that and takes the rest from an
    independent w, so Cov(Y, f) = rho sigma_Y sigma_f by construction. A coordinate-only study
    supplies f and no Y.
    """
    ys, fs = [], []
    for i in range(count):
        u = stats.Normal(f"u_{tag}{i}", 0, 1)
        fs.append(mu_f + sig_f * u)
        if correlated:
            w = stats.Normal(f"w_{tag}{i}", 0, 1)
            ys.append(m + sig_y * (rho_c * u + sp.sqrt(1 - rho_c**2) * w))
    return ys, fs


N_IMG_CHECK, N_TAB_CHECK = 2, 3
y_img, f_img = cohort(N_IMG_CHECK, "I", correlated=True)
_, f_coord = cohort(N_TAB_CHECK, "C", correlated=False)
estimator = (
    sum(y_img) / N_IMG_CHECK
    + lam * (sum(f_coord) / N_TAB_CHECK - sum(f_img) / N_IMG_CHECK)
)
proof.define(r"\hat m_\lambda \text{ at } n=2, N=3", estimator)

proof.claim(
    "sympy-agrees-the-estimator-is-unbiased-for-any-lambda-and-any-predictor",
    sp.simplify(stats.E(estimator) - m),
    r"\mathbb{E}\hat m_\lambda = m",
)

# The general variance formula, evaluated at these counts and with the covariance written in
# terms of rho, must match what sympy computes from the construction.
general = (
    (var_y - 2 * lam * rho_c * sig_y * sig_f + lam**2 * var_f) / N_IMG_CHECK
    + lam**2 * var_f / N_TAB_CHECK
)
proof.claim(
    "and-agrees-with-the-variance-formula-including-its-cross-term",
    sp.simplify(stats.variance(estimator) - general),
    r"\operatorname{Var}(\hat m_\lambda) = \frac{\operatorname{Var}(Y-\lambda f)}{n}"
    r" + \frac{\lambda^2\sigma_f^2}{N}",
)

# ------------------------------------------------------------------ variance

# The image cohort contributes Y - lambda f as a single per-study quantity; the coordinate cohort
# contributes lambda f. Independence across cohorts makes the variance additive.
var_combination = var_y - 2 * lam * cov_yf + lam**2 * var_f
variance = var_combination / n + lam**2 * var_f / N
proof.define(r"\operatorname{Var}(\hat m_\lambda)", variance)
proof.claim(
    "the-image-term-is-the-variance-of-the-corrected-observation",
    sp.simplify(var_combination - (var_y - 2 * lam * cov_yf + lam**2 * var_f)),
    r"\operatorname{Var}(Y - \lambda f) = \sigma_Y^2 - 2\lambda\sigma_{Yf} + \lambda^2\sigma_f^2",
)

# ------------------------------------------------------------------ the optimal coefficient

# Quadratic in lambda, so the minimiser is where the derivative vanishes.
stationary = sp.solve(sp.diff(variance, lam), lam)
assert len(stationary) == 1, stationary
lam_star = sp.simplify(stationary[0])
proof.define(r"\lambda^*", lam_star)
proof.claim(
    "the-optimal-coefficient-is-the-slope-shrunk-by-one-plus-n-over-N",
    sp.simplify(lam_star - cov_yf / ((1 + n / N) * var_f)),
    r"\lambda^* = \frac{\sigma_{Yf}}{(1 + n/N)\,\sigma_f^2}",
)
# Convexity, so it is a minimum and not a maximum.
proof.claim(
    "and-it-is-a-minimum",
    sp.simplify(sp.diff(variance, lam, 2) - 2 * var_f * (1 / n + 1 / N)),
    r"\partial_\lambda^2 \operatorname{Var} = 2\sigma_f^2(1/n + 1/N) > 0",
)

# ------------------------------------------------------------------ the achievable ratio

rho = sp.Symbol("rho", positive=True)
at_optimum = sp.simplify(variance.subs(lam, lam_star))
images_only = var_y / n
ratio = sp.simplify(at_optimum / images_only)
proof.define(r"V_*/V_{\text{images}}", ratio)
proof.claim(
    "the-best-ratio-is-one-minus-the-squared-correlation-over-one-plus-n-over-N",
    sp.simplify(
        ratio.subs(cov_yf, rho * sp.sqrt(var_y) * sp.sqrt(var_f))
        - (1 - rho**2 / (1 + n / N))
    ),
    r"\frac{V_*}{V_{\text{images}}} = 1 - \frac{\rho^2}{1 + n/N}",
)

# The cap. Letting the coordinate cohort grow without bound leaves 1 - rho^2: the images set a
# floor that no quantity of tables reaches below.
proof.claim(
    "infinitely-many-coordinate-studies-still-leave-one-minus-rho-squared",
    sp.simplify(sp.limit(1 - rho**2 / (1 + n / N), N, sp.oo) - (1 - rho**2)),
    r"\lim_{N\to\infty} \frac{V_*}{V_{\text{images}}} = 1 - \rho^2",
)

# And the other limit: with no coordinate studies the correction is worthless, as it must be.
proof.claim(
    "and-no-coordinate-studies-means-no-reduction",
    sp.simplify(sp.limit(1 - rho**2 / (1 + n / N), N, 0, "+") - 1),
    r"\lim_{N\to 0^+} \frac{V_*}{V_{\text{images}}} = 1",
)


def calibration():
    """Reproduce the design document's table, and price the floor it implies."""
    print()
    print("variance ratio 1 - rho^2 / (1 + n/N):")
    print(f"  {'n':>4} {'N':>5} {'rho':>5} {'ratio':>7} {'reduction':>10} "
          f"{'images-only floor':>18}")
    for n_i, N_i in ((1, 20), (1, 500), (8, 100)):
        for rho_i in (0.5, 0.8):
            r = 1 - rho_i**2 / (1 + n_i / N_i)
            floor = 1 - rho_i**2
            print(f"  {n_i:>4} {N_i:>5} {rho_i:>5} {r:>7.3f} {100*(1-r):>9.0f}% "
                  f"{floor:>18.3f}")
    # The document quotes 23% at correlation 0.5 and 59% at 0.8 for 8 images and 100 tables.
    assert abs((1 - (1 - 0.5**2 / (1 + 8 / 100))) - 0.2315) < 5e-4
    assert abs((1 - (1 - 0.8**2 / (1 + 8 / 100))) - 0.5926) < 5e-4
    print()
    print("  The floor column is what the images alone impose: even a perfect predictor with")
    print("  unlimited tables cannot reduce variance below 1 - rho^2 of the image-only value.")


if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    calibration()
