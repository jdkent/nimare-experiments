"""CBES's implemented derivatives: are the closed forms in the code the right ones?

`censoring.py` proves things about the *model*. This proves things about the *implementation*.
Three hand-derived closed forms in `nimare/meta/cbma/effectsize.py` have never been checked
against symbolic differentiation, and one wrong limit has already shipped from this exact
neighbourhood:

  1. `_censoring_terms` returns `prob`, `score` and `d2_over_prob` for a reporting indicator,
     with a `sign` telling the silent limb from the reported one;
  2. `_mu_derivatives` assembles the score and curvature of the weighted log-likelihood in `mu`,
     using `d2_over_prob - score**2` for the censored part;
  3. `_observed_information` builds the three blocks of the 2x2 observed information from a
     responsibility `r`, a per-observation score and a per-observation Hessian.

Each is stated here as the implementation writes it and checked against sympy's derivative of
the thing it claims to differentiate. A failure is a bug in the estimator, not in the proof.
"""
import sympy as sp

from proofs.latex import Proof

PREAMBLE = r"""Let $S(\mu)=\Phi(u)-\Phi(l)$ with $u=(c-\mu)/\sigma$ and $l=(-c-\mu)/\sigma$ be
the probability a study is silent, and let $\varepsilon=\pm1$ be the indicator's sign: $+1$ for
a silent pair, whose probability is $S$, and $-1$ for a pair that reported, whose probability is
$1-S$. Write $p_\varepsilon$ for the probability of the event actually observed."""


def normal_pdf(x):
    return sp.exp(-x**2 / 2) / sp.sqrt(2 * sp.pi)


def normal_cdf(x):
    return (1 + sp.erf(x / sp.sqrt(2))) / 2


def run():
    proof = Proof("censored_information", "CBES's implemented derivatives", preamble=PREAMBLE)
    mu, sigma, c = sp.symbols("mu sigma c", positive=True)
    pi, g, prec, w = sp.symbols("pi g tau w", positive=True)

    upper = (c - mu) / sigma
    lower = (-c - mu) / sigma
    silent = normal_cdf(upper) - normal_cdf(lower)
    proof.define(r"S(\mu)=\Phi\big(\tfrac{c-\mu}{\sigma}\big)"
                 r"-\Phi\big(\tfrac{-c-\mu}{\sigma}\big)", silent)

    # --- 1. _censoring_terms' `first` and `second` -----------------------------------------
    # The code computes, before any normalisation:
    #     first  = -(phi(upper) - phi(lower)) / sigma
    #     second = -(upper phi(upper) - lower phi(lower)) / sigma^2
    # and claims these are dS/dmu and d2S/dmu2. Note the second is *not* the derivative of the
    # first by the chain rule on its face -- the sign works out only because d(phi)/dx = -x phi,
    # which is exactly the kind of step worth checking rather than trusting.
    first = -(normal_pdf(upper) - normal_pdf(lower)) / sigma
    second = -(upper * normal_pdf(upper) - lower * normal_pdf(lower)) / sigma**2
    proof.claim("`first` is dS/dmu", sp.simplify(sp.diff(silent, mu) - first),
                r"\partial_\mu S - \texttt{first}")
    proof.claim("`second` is d2S/dmu2", sp.simplify(sp.diff(silent, mu, 2) - second),
                r"\partial^2_\mu S - \texttt{second}")

    # --- 2. the sign generalisation is the complement's derivative, not a fudge -------------
    # The code writes, for either limb,
    #     prob = S if sign > 0 else 1 - S
    #     score = first * sign / prob
    #     d2_over_prob = second * sign / prob
    # For the reported limb this must be the derivative of log(1 - S), and the single `sign`
    # factor has to serve both derivatives. That is true because 1 - S differentiates to -S' and
    # -S'' alike, so one sign flip covers both -- but it is asserted in the code and proved here.
    for label, sign, prob in [("silent limb", 1, silent), ("reported limb", -1, 1 - silent)]:
        proof.claim(
            f"{label}: score is d/dmu log p",
            sp.simplify(sp.diff(sp.log(prob), mu) - first * sign / prob),
            r"\partial_\mu \log p_\varepsilon - \varepsilon\,\texttt{first}/p_\varepsilon",
        )
        # `d2_over_prob - score**2` is what _mu_derivatives adds to the curvature, and it must
        # be the second derivative of the log. This is the identity
        #     d2/dmu2 log p = p''/p - (p'/p)^2
        # with p'' carrying the same single sign factor as p'.
        proof.claim(
            f"{label}: `d2_over_prob - score^2` is d2/dmu2 log p",
            sp.simplify(sp.diff(sp.log(prob), mu, 2)
                        - (second * sign / prob - (first * sign / prob) ** 2)),
            r"\partial^2_\mu \log p_\varepsilon"
            r" - \big(\varepsilon\,\texttt{second}/p_\varepsilon - \texttt{score}^2\big)",
        )

    # --- 3. the image limb's score and curvature -------------------------------------------
    # _mu_derivatives uses (g - mu) * precision and -precision for a value-bearing pair, which
    # are the derivatives of a normal log-density with precision held fixed.
    log_density = -prec * (g - mu) ** 2 / 2
    proof.claim("an image pair's score is (g - mu) * precision",
                sp.simplify(sp.diff(log_density, mu) - (g - mu) * prec),
                r"\partial_\mu \log\phi - (g-\mu)\tau")
    proof.claim("an image pair's curvature is minus the precision",
                sp.simplify(sp.diff(log_density, mu, 2) + prec),
                r"\partial^2_\mu \log\phi + \tau")

    # --- 4. the three observed-information blocks ------------------------------------------
    # _observed_information computes, per observation of mixture density f = pi a + (1-pi) b,
    # with r = pi a / f the responsibility, `score` = d/dmu log a and `hessian` = d2/dmu2 log a:
    #     i_mu   = -w (r * hessian + r(1-r) * score^2)
    #     cross  = -w * score * r(1-r) / (pi (1-pi))
    #     i_pi   =  w (r/pi - (1-r)/(1-pi))^2
    # Each must equal minus the corresponding second derivative of w log f. The `a` limb is
    # given an explicit mu-dependence so the derivatives are not vacuous; `b` has none, which
    # is the structural fact the whole estimator turns on.
    a = sp.Function("a")(mu)
    b = sp.Symbol("b", positive=True)
    f = pi * a + (1 - pi) * b
    ell = w * sp.log(f)
    r_resp = pi * a / f
    a_score = sp.diff(sp.log(a), mu)
    a_hess = sp.diff(sp.log(a), mu, 2)

    i_mu = -w * (r_resp * a_hess + r_resp * (1 - r_resp) * a_score**2)
    proof.claim(
        "i_mu equals minus d2/dmu2 of the weighted log-mixture",
        sp.simplify(-sp.diff(ell, mu, 2) - i_mu),
        r"-\partial^2_\mu \ell - i_{\mu\mu}",
    )
    cross = -w * a_score * r_resp * (1 - r_resp) / (pi * (1 - pi))
    proof.claim(
        "the cross block equals minus the mixed partial",
        sp.simplify(-sp.diff(ell, mu, pi) - cross),
        r"-\partial^2_{\mu\pi}\ell - i_{\mu\pi}",
    )
    i_pi = w * (r_resp / pi - (1 - r_resp) / (1 - pi)) ** 2
    proof.claim(
        "i_pi equals minus d2/dpi2, exactly, because f is affine in pi",
        sp.simplify(-sp.diff(ell, pi, 2) - i_pi),
        r"-\partial^2_\pi \ell - i_{\pi\pi}",
    )
    # Anti-vacuity: for a generic mu-dependent `a` the blocks must not collapse to zero, or
    # the three claims above hold only because both sides are empty.
    concrete = sp.exp(-mu)
    assert sp.simplify(i_mu.subs(a, concrete).doit()) != 0, "i_mu must not vanish"
    assert sp.simplify(cross.subs(a, concrete).doit()) != 0, "the cross block must not vanish"
    assert sp.simplify(i_pi.subs(a, concrete).doit()) != 0, "i_pi must not vanish"
    return proof


if __name__ == "__main__":
    run().report()
