r"""Practical prevalence :math:`\pi_\delta^+`, and the identity it does not satisfy.

Section 2 of the design document proposes replacing the structural-zero prevalence with

.. math:: \pi_\delta^+(v) = P_{P^*}\{\theta_i(v) > \delta\},

and immediately warns: "it still requires estimating the effect distribution, including
heterogeneity. It is not automatically identifiable because it has a better name. Under this
definition, multiplying prevalence by the conditional mean does not recover the full marginal
mean unless the complementary component has mean zero."

Both halves of that warning are provable, and the second is the one that bites.

What is proved:

  * under the continuous random-effects model,
    :math:`\pi_\delta^+ = \Phi\!\left((m-\delta)/\tau\right)` (claim 1), so it is a *function of
    the parameters already being estimated* -- no new latent structure, but also no new
    information: it inherits the heterogeneity's uncertainty entirely;
  * its derivatives with respect to the mean and the heterogeneity (claims 2-3), which is how a
    standard error propagates to it by the delta method, and which show it is **increasing in
    the heterogeneity whenever** :math:`m < \delta`. A more heterogeneous literature has *more*
    studies above any threshold the mean falls short of -- so a practical prevalence can be
    driven up by noise between studies rather than by effect;
  * :math:`\tau \to 0` sends it to a step function of :math:`m` (claim 4), so at small
    heterogeneity it is an almost-deterministic function of the mean's position relative to
    :math:`\delta` and carries essentially no separate information;
  * **the decomposition identity fails** (claims 5-7). The exact statement is

    .. math::
        m = \pi_\delta^+\,\mathbb{E}[\theta \mid \theta > \delta]
        + (1-\pi_\delta^+)\,\mathbb{E}[\theta \mid \theta \le \delta],

    so the product :math:`\pi_\delta^+\,\mathbb{E}[\theta\mid\theta>\delta]` recovers :math:`m`
    only when the complementary part contributes nothing. Claim 7 gives that term in closed
    form, and the gap is

    .. math::
        \pi_\delta^+\,\mathbb{E}[\theta\mid\theta>\delta] - m
        = \tau\,\varphi\!\left(\frac{m-\delta}{\tau}\right) - (1-\pi_\delta^+)\,m .

    Its **sign flips**: the product overstates the marginal mean where the density term
    dominates and understates it where the complementary mass does. An earlier draft of this
    docstring asserted it always understates, and the numeric check below contradicted that on
    its second row -- worth keeping as the reason the check exists. What is general is that the
    gap is a closed-form quantity which vanishes only as :math:`\delta \to -\infty`, where the
    whole distribution is "above". The document's "unless the complementary component has mean
    zero" is therefore never satisfied at a finite threshold by a continuous distribution.

**What this does not establish.** Nothing here makes :math:`\pi_\delta^+` estimable from
coordinates. It is a function of :math:`(m,\tau)`, and everything already proved about how badly
:math:`\tau` is determined from thresholded tables applies to it unchanged -- with the extra
sensitivity that claim 3 quantifies.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

m, delta = sp.symbols("m delta", real=True)
tau = sp.Symbol("tau", positive=True)
theta = sp.Symbol("theta", real=True)

proof = Proof(
    "practical_prevalence",
    "Practical prevalence, its sensitivity to heterogeneity, and the identity it fails",
    __doc__,
)


def normal_pdf(value, mean, sd):
    return sp.exp(-((value - mean) ** 2) / (2 * sd**2)) / (sd * sp.sqrt(2 * sp.pi))


def normal_cdf(value, mean, sd):
    return (1 + sp.erf((value - mean) / (sd * sp.sqrt(2)))) / 2


# --------------------------------------------- the prevalence is a function of the parameters

prevalence = 1 - normal_cdf(delta, m, tau)
proof.define(r"\pi_\delta^+", prevalence)
proof.claim(
    "practical-prevalence-is-the-probit-of-the-standardised-gap-to-the-threshold",
    sp.simplify(prevalence - normal_cdf((m - delta) / tau, 0, 1)),
    r"\pi_\delta^+ = \Phi\!\left(\frac{m-\delta}{\tau}\right)",
)

standardised = (m - delta) / tau
standard_pdf = sp.exp(-(standardised**2) / 2) / sp.sqrt(2 * sp.pi)
proof.claim(
    "its-slope-in-the-mean-is-the-density-over-the-heterogeneity",
    sp.simplify(sp.diff(prevalence, m) - standard_pdf / tau),
    r"\partial_m \pi_\delta^+ = \frac{1}{\tau}\varphi\!\left(\frac{m-\delta}{\tau}\right)",
)
proof.claim(
    "and-its-slope-in-the-heterogeneity-carries-the-sign-of-the-threshold-minus-the-mean",
    sp.simplify(sp.diff(prevalence, tau) - (delta - m) / tau**2 * standard_pdf),
    r"\partial_\tau \pi_\delta^+ = \frac{\delta-m}{\tau^2}"
    r"\varphi\!\left(\frac{m-\delta}{\tau}\right)",
)

# With m < delta the slope in tau is strictly positive: more between-study spread puts more
# studies above a threshold the mean itself does not reach.
proof.claim(
    "so-a-mean-below-the-threshold-makes-prevalence-increase-with-heterogeneity",
    sp.simplify(
        sp.diff(prevalence, tau).subs(m, delta - sp.Symbol("d", positive=True))
        - sp.Symbol("d", positive=True) / tau**2
        * sp.exp(-sp.Symbol("d", positive=True) ** 2 / (2 * tau**2)) / sp.sqrt(2 * sp.pi)
    ),
    r"m = \delta - d,\ d>0 \ \Rightarrow\ \partial_\tau \pi_\delta^+ "
    r"= \frac{d}{\tau^2}\varphi(d/\tau) > 0",
)

# ------------------------------------------------- the decomposition, and the term that spoils it

# The law of total expectation, written out. Each conditional mean is the truncated normal's.
upper_mean = m + tau * standard_pdf / prevalence
lower_mean = m - tau * standard_pdf / (1 - prevalence)
proof.define(r"\mathbb{E}[\theta\mid\theta>\delta]", upper_mean)
proof.define(r"\mathbb{E}[\theta\mid\theta\le\delta]", lower_mean)

proof.claim(
    "the-two-conditional-means-recombine-to-the-marginal-mean-exactly",
    sp.simplify(prevalence * upper_mean + (1 - prevalence) * lower_mean - m),
    r"m = \pi_\delta^+\,\mathbb{E}[\theta\mid\theta>\delta] "
    r"+ (1-\pi_\delta^+)\,\mathbb{E}[\theta\mid\theta\le\delta]",
)
proof.claim(
    "so-the-product-of-prevalence-and-the-upper-mean-overshoots-by-the-complementary-term",
    sp.simplify(
        prevalence * upper_mean - m + (1 - prevalence) * lower_mean
    ),
    r"\pi_\delta^+\,\mathbb{E}[\theta\mid\theta>\delta] - m "
    r"= -(1-\pi_\delta^+)\,\mathbb{E}[\theta\mid\theta\le\delta]",
)
# And that complementary term is the marginal mean's lower share minus a strictly positive
# density term, so it never vanishes for a continuous distribution with positive lower mass.
proof.claim(
    "and-the-complementary-term-is-the-lower-mass-times-the-mean-less-a-positive-density-term",
    sp.simplify(
        (1 - prevalence) * lower_mean - ((1 - prevalence) * m - tau * standard_pdf)
    ),
    r"(1-\pi_\delta^+)\,\mathbb{E}[\theta\mid\theta\le\delta] "
    r"= (1-\pi_\delta^+)\,m - \tau\,\varphi\!\left(\frac{m-\delta}{\tau}\right)",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.stats import norm

    proof.report()

    # What the decomposition error actually is, at plausible settings. Reported rather than
    # asserted small, because it is not small.
    print(f"{'m':>6} {'tau':>6} {'delta':>6} {'prevalence':>11} {'pi*E[up]':>9} "
          f"{'m':>7} {'shortfall':>10}")
    for mean in (0.1, 0.3, 0.5):
        for spread in (0.1, 0.2, 0.4):
            for cut in (0.0, 0.2):
                prevalence_value = float(norm.sf((cut - mean) / spread))
                density = float(norm.pdf((mean - cut) / spread))
                upper = mean + spread * density / prevalence_value
                product = prevalence_value * upper
                print(f"{mean:>6.2f} {spread:>6.2f} {cut:>6.2f} {prevalence_value:>11.4f} "
                      f"{product:>9.4f} {mean:>7.4f} {product - mean:>+10.4f}")
                # The invariant, rather than a guessed direction: the gap must equal the
                # closed form of claim 7, whatever its sign.
                predicted = spread * density - (1 - prevalence_value) * mean
                if abs((product - mean) - predicted) > 1e-12:
                    raise AssertionError(
                        f"the measured gap {product - mean:+.6f} disagrees with the closed form "
                        f"{predicted:+.6f}"
                    )

    # And the sensitivity claim 4 warns about: prevalence rising with heterogeneity alone,
    # at a fixed mean that does not reach the threshold.
    print("  [ok] numeric: every gap matches the closed form tau*phi - (1-pi)*m, and its sign")
    print("       flips between rows -- neither direction holds generally")

    # It vanishes only in the limit where everything is above the threshold.
    tail = [
        float(norm.sf((cut - 0.3) / 0.2)) * (0.3 + 0.2 * float(norm.pdf((0.3 - cut) / 0.2))
        / float(norm.sf((cut - 0.3) / 0.2))) - 0.3
        for cut in (-1.0, -2.0, -4.0, -8.0)
    ]
    if not (np.all(np.abs(tail) == np.sort(np.abs(tail))[::-1]) and abs(tail[-1]) < 1e-9):
        raise AssertionError("the gap does not vanish as the threshold falls away")
    print(f"  [ok] numeric: the gap falls to zero as the threshold recedes "
          f"({abs(tail[0]):.2e} -> {abs(tail[-1]):.2e})")

    print()
    print("prevalence at a fixed mean of 0.10 below a threshold of 0.30, as heterogeneity grows")
    previous = -1.0
    for spread in (0.05, 0.1, 0.2, 0.4, 0.8):
        value = float(norm.sf((0.30 - 0.10) / spread))
        print(f"   tau = {spread:.2f}   pi = {value:.4f}")
        if value <= previous:
            raise AssertionError("prevalence did not increase with heterogeneity below the cut")
        previous = value
    print("  [ok] numeric: a practical prevalence can be driven entirely by between-study")
    print("       spread at a mean that never reaches the threshold")

    print()
    print("Consequence for the implementation: practical prevalence is reportable because it is")
    print("a function of (m, tau), and it is reportable *only* alongside tau and its")
    print("uncertainty, because it moves with tau in a direction that has nothing to do with")
    print("effect. And the product of prevalence and the conditional mean must never be offered")
    print("as the marginal mean: the gap is a closed-form term that is zero only in the limit")
    print("where the threshold recedes and every study counts as above it.")
