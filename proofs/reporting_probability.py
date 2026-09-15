"""Why the smooth-field maxima density has the sign backwards at a signal peak.

The censoring term uses `P(|g| >= c)` as the probability a study reports a voxel. A paper reports
a voxel only if it cleared `c` *and* the value there was a local maximum, so the model's
probability is too large. Measured, the over-statement runs 1.00, 1.78, 1.08 and 1.13 at true `g`
of 0.2, 0.4, 0.6 and 0.8 -- non-monotone, so no constant reweighting absorbs it.

The obvious correction is the random-field expected-maxima density, and it was tried. It made
every focus worse, and no constant `q` helped either. The numerics said no without saying why.
This says why, and the reason is a term the standard density does not contain.

Write the observed field as `Z(x) = m(x) + e(x)` with `e` a smooth stationary Gaussian field of
mean zero and `m` the signal. A local maximum at the origin requires `Z'(0) = 0` and
`Z''(0) < 0`. At a *signal peak* the mean is stationary, `m'(0) = 0`, so the first condition is a
statement about the noise alone -- but the second is not: `Z''(0) = m''(0) + e''(0)`, and at a
peak `m''(0) < 0`. The signal's own curvature therefore makes a local maximum MORE likely
exactly where the effect is, and the expected-maxima density of a zero-mean field is the special
case `m'' = 0`. Applying it at a signal peak understates the probability of a maximum, so the
correction it supplies is too large, in the direction the measurements found.
"""
import sympy as sp

from proofs.latex import Proof

PREAMBLE = r"""Let $Z=m+e$ on $\mathbb{R}$, with $e$ a smooth stationary Gaussian field of mean
zero, $\operatorname{sd}(e''(0))=\sigma_2$, and $m$ the signal. At a signal peak $m'(0)=0$ and
$\kappa:=-m''(0)>0$. A local maximum at $0$ needs $Z'(0)=0$ and $Z''(0)<0$."""


def normal_cdf(x):
    return (1 + sp.erf(x / sp.sqrt(2))) / 2


def run():
    proof = Proof("reporting_probability", "The maxima density at a signal peak",
                  preamble=PREAMBLE)
    kappa, sigma2 = sp.symbols("kappa sigma_2", positive=True)

    # --- 1. the curvature condition, conditional on the signal's own curvature -------------
    # Z''(0) = m''(0) + e''(0) = -kappa + e''(0), so
    #     P(Z''(0) < 0) = P(e''(0) < kappa) = Phi(kappa / sigma_2).
    curvature_probability = normal_cdf(kappa / sigma2)
    proof.define(r"P\big(Z''(0)<0\big)=\Phi(\kappa/\sigma_2)", curvature_probability)

    # The zero-mean field is the case kappa = 0, where the curvature condition is a coin flip.
    proof.claim(
        "a zero-mean field puts the curvature condition at exactly one half",
        sp.simplify(curvature_probability.subs(kappa, 0) - sp.Rational(1, 2)),
        r"\Phi(0)-\tfrac12",
    )

    # --- 2. and it is strictly increasing in the peak's sharpness --------------------------
    # So the standard density is not merely different at a signal peak, it is a LOWER BOUND.
    derivative = sp.diff(curvature_probability, kappa)
    proof.claim(
        "the probability increases in kappa: the derivative is a positive density over sigma_2",
        sp.simplify(derivative - sp.exp(-kappa**2 / (2 * sigma2**2))
                    / (sigma2 * sp.sqrt(2 * sp.pi))),
        r"\partial_\kappa \Phi(\kappa/\sigma_2)-\phi(\kappa/\sigma_2)/\sigma_2",
    )
    assert sp.simplify(derivative.subs({kappa: 1, sigma2: 1})) > 0, (
        "the derivative must be strictly positive, or claim 2 is vacuous")

    # --- 3. how wrong the zero-mean density is, as a factor --------------------------------
    # The ratio of the true curvature probability to the one a zero-mean density supplies is
    # 2 Phi(kappa / sigma_2), which is 1 at kappa = 0 and tends to 2 as the peak sharpens. So a
    # null-field density can understate the chance of a local maximum by up to a factor of two,
    # and it understates it most where the signal is sharpest. That gives the *direction* of the
    # error and a bound on its size; it does NOT by itself produce the non-monotone pattern in
    # the measurements (1.00, 1.78, 1.08, 1.13), which would need the interaction with the
    # height threshold as well, and is not derived here.
    ratio = curvature_probability / sp.Rational(1, 2)
    proof.claim(
        "the understatement factor is 2 Phi(kappa/sigma_2), which is 1 at kappa = 0",
        sp.simplify(ratio.subs(kappa, 0) - 1),
        r"2\Phi(\kappa/\sigma_2)\big|_{\kappa=0}-1",
    )
    proof.claim(
        "and tends to 2 as the peak sharpens",
        sp.limit(ratio, kappa, sp.oo) - 2,
        r"\lim_{\kappa\to\infty} 2\Phi(\kappa/\sigma_2)-2",
    )

    # The gradient condition is left in the preamble rather than made a claim. At a signal peak
    # Z'(0) = m'(0) + e'(0) = e'(0), so the first-order condition is the same statement about
    # the noise as it is under the null -- but that is true *by the definition* of a peak, and
    # asserting `0 == 0` in sympy would look like a proof while being one. It matters for the
    # diagnosis all the same: exactly one of the two conditions defining a local maximum picks
    # up the signal, and it is the one the zero-mean density pins at one half.
    return proof


if __name__ == "__main__":
    run().report()
