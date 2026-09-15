"""The zero-inflated censored likelihood: what a silence can and cannot identify.

Five claims, each of which a simulation in this program either reached the long way or got
wrong. Together they say why CBES's `se` is not a bound on `mu`, which kind of study
heterogeneity makes a magnitude recoverable, and where one shipped bug came from.

The model. Study `i` either has the effect (probability `pi`) or has none. Its contribution to
the likelihood at one voxel is a two-component mixture

    f_i = pi * a_i + (1 - pi) * b_i,

with `a_i` the density or probability of what was observed under an effect of `mu` and `b_i` the
same under an effect of exactly zero. For an image `a_i` is a normal density at `mu`; for a
*silence* it is `S(mu) = P(|g| < c | mu)`; for a *report* it is `1 - S(mu)`. Only `a_i` carries
`mu`, and that single fact drives everything below.
"""
import sympy as sp

from proofs.latex import Proof

PREAMBLE = r"""Let $f=\pi a(\mu)+(1-\pi)b$ be one study's mixture contribution at a voxel,
with $b$ the density or probability under an effect of exactly zero. A silence has
$a(\mu)=S(\mu)=\Phi\!\big(\tfrac{c-\mu}{\sigma}\big)-\Phi\!\big(\tfrac{-c-\mu}{\sigma}\big)$
and a report has $a(\mu)=1-S(\mu)$. Write $\ell=\log f$."""


def normal_cdf(x):
    """Phi, written through erf so sympy can differentiate and take limits of it."""
    return (1 + sp.erf(x / sp.sqrt(2))) / 2


def silent_probability(mu, sigma, c):
    """S(mu) = P(|g| < c) for g ~ N(mu, sigma^2)."""
    return normal_cdf((c - mu) / sigma) - normal_cdf((-c - mu) / sigma)


def run():
    proof = Proof("censoring", "The censored mixture's identifiability", preamble=PREAMBLE)
    pi, mu, sigma, c, eps = sp.symbols("pi mu sigma c epsilon", positive=True)
    a, b = sp.symbols("a b", positive=True)
    n, z = sp.symbols("n z", positive=True)

    # --- 1. the pi block is exact, not an approximation ------------------------------------
    # f is affine in pi, so the second derivative of the log is exactly minus the square of the
    # first. This is why an outer product of scores can sit in the same information matrix as a
    # Hessian-based mu block without mixing two different estimators of information.
    f = pi * a + (1 - pi) * b
    ell = sp.log(f)
    proof.claim(
        "d2/dpi2 log f = -(d/dpi log f)^2",
        sp.simplify(sp.diff(ell, pi, 2) + sp.diff(ell, pi) ** 2),
        r"\partial^2_\pi \ell + (\partial_\pi \ell)^2",
    )
    # The same expression as a negative square is concavity in pi at fixed mu, which licenses
    # the fixed-point inner solve in the profile interval as a global maximiser.
    proof.claim(
        "log f is concave in pi: the curvature is minus a square",
        sp.simplify(sp.diff(ell, pi, 2) + ((a - b) / f) ** 2),
        r"\partial^2_\pi \ell + \big((a-b)/f\big)^2",
    )

    # --- 2. one study configuration cannot separate pi from mu -----------------------------
    # Every pair sharing (sigma, c) has the SAME score direction, so any number of them gives a
    # rank-one information matrix. That singularity is the ridge, and it is exact rather than
    # numerical.
    S_mu = silent_probability(mu, sigma, c)
    S_null = silent_probability(0, sigma, c)
    P = pi * S_mu + (1 - pi) * S_null
    score = sp.Matrix([sp.diff(sp.log(P), mu), sp.diff(sp.log(P), pi)])
    # The content is not that one score direction spans a line -- that is trivial -- but that a
    # SILENT and a REPORTED pair at the same (sigma, c) give *collinear* scores. Both depend on
    # (pi, mu) only through P, so their gradients are both multiples of grad P and no mixture of
    # observed outcomes at one configuration escapes the line. This is the ridge, and it is the
    # invariant measured constant to 0.4% along the profile while pi*mu swung by hundreds.
    reported = sp.Matrix([sp.diff(sp.log(1 - P), mu), sp.diff(sp.log(1 - P), pi)])
    proof.claim(
        "a silence and a report at one configuration have collinear scores",
        sp.simplify(score[0] * reported[1] - score[1] * reported[0]),
        r"s_\mu r_\pi - s_\pi r_\mu\ \textrm{for}\ s=\nabla\log P,\ r=\nabla\log(1-P)",
    )
    k = sp.Symbol("k", positive=True, integer=True)   # any number of identical studies
    information = k * score * score.T + reported * reported.T
    proof.claim(
        "so any number of them, silent or reporting, is still singular",
        sp.simplify(information.det()),
        r"\det\big(k\,s s^\top + r r^\top\big)",
    )
    # Anti-vacuity: two *distinct* configurations must not be singular, or the claim says
    # nothing about heterogeneity helping.
    sigma2, c2 = sp.symbols("sigma_2 c_2", positive=True)
    P2 = pi * silent_probability(mu, sigma2, c2) + (1 - pi) * silent_probability(0, sigma2, c2)
    score2 = sp.Matrix([sp.diff(sp.log(P2), mu), sp.diff(sp.log(P2), pi)])
    two = score * score.T + score2 * score2.T
    numeric = two.det().subs(
        {mu: sp.Rational(1, 2), pi: sp.Rational(3, 5), sigma: sp.Rational(1, 5),
         c: sp.Rational(2, 5), sigma2: sp.Rational(1, 10), c2: sp.Rational(9, 10)})
    assert abs(float(numeric.evalf())) > 1e-12, (
        "two distinct configurations must give a nonsingular information matrix, "
        "or claim 2 is vacuous")

    # --- 3. the profile likelihood has a horizontal asymptote in mu ------------------------
    # A silence's and an image's mu-dependence both vanish as |mu| grows; a report's does not,
    # tending to 1. So the limit is free of mu but still a function of pi -- which is exactly
    # the boundedness criterion the estimator now documents.
    proof.claim(
        "a silence loses all information about mu in the limit",
        sp.limit(S_mu, mu, sp.oo),
        r"\lim_{\mu\to\infty} S(\mu)",
    )
    proof.claim(
        "a report's probability tends to one, so pi survives the limit",
        sp.limit(1 - S_mu, mu, sp.oo) - 1,
        r"\lim_{\mu\to\infty}\big(1-S(\mu)\big) - 1",
    )
    # Hence the silent limb's contribution tends to log((1 - pi) S(0)), with no mu in it.
    silent_term = sp.log(pi * S_mu + (1 - pi) * S_null)
    proof.claim(
        "the silent limb's limit is free of mu",
        sp.limit(silent_term, mu, sp.oo) - sp.log((1 - pi) * S_null),
        r"\lim_{\mu\to\infty}\log\big(\pi S(\mu)+(1-\pi)S(0)\big)-\log\big((1-\pi)S(0)\big)",
    )

    # --- 4. the responsibility limit that a shipped patch got wrong ------------------------
    # Holding pi at 1 was claimed to make the cross block vanish because r -> 1 and r(1-r) -> 0.
    # With pi clamped at 1 - eps the limit is a RATIO OF DENSITIES, order one, not zero: the
    # Schur complement would have gone on subtracting a term no parameter earned.
    r = (1 - eps) * a / ((1 - eps) * a + eps * b)
    proof.claim(
        "(1 - r)/(1 - pi) tends to b/a, not to zero",
        sp.simplify(sp.limit((1 - r) / eps, eps, 0) - b / a),
        r"\lim_{\epsilon\to0}\frac{1-r}{\epsilon}-\frac{b}{a}",
    )
    assert sp.simplify(b / a) != 0, "b/a is not zero, which is the whole point of claim 4"

    # --- 5. why sample size identifies and the threshold does not --------------------------
    # The cutoff measured in sampling standard deviations is the reported statistic again. With
    # the null sampling sd sigma = 1/sqrt(n) and a cutoff of z/sqrt(n) in effect-size units,
    # c/sigma = z exactly -- free of n. So the NULL component's silence probability depends on
    # the threshold alone, while the active one keeps mu*sqrt(n).
    sigma_n = 1 / sp.sqrt(n)
    c_n = z / sp.sqrt(n)
    proof.claim(
        "the cutoff in sampling-sd units is the reported statistic, independent of n",
        sp.simplify(c_n / sigma_n - z),
        r"c/\sigma - z\ \textrm{at}\ \sigma=n^{-1/2},\ c=zn^{-1/2}",
    )
    S_null_n = silent_probability(0, sigma_n, c_n)
    proof.claim(
        "the null component's silence probability does not depend on n",
        sp.simplify(sp.diff(S_null_n, n)),
        r"\partial_n S(0;\,n)",
    )
    S_mu_n = silent_probability(mu, sigma_n, c_n)
    proof.claim(
        "the active component's does, and only through mu*sqrt(n)",
        sp.simplify(S_mu_n - (normal_cdf(z - mu * sp.sqrt(n))
                              - normal_cdf(-z - mu * sp.sqrt(n)))),
        r"S(\mu;n)-\big[\Phi(z-\mu\sqrt n)-\Phi(-z-\mu\sqrt n)\big]",
    )
    # Anti-vacuity: that derivative must actually be nonzero, or "n identifies" is empty.
    d_active = sp.diff(S_mu_n, n).subs({mu: sp.Rational(1, 2), z: 3, n: 30})
    assert abs(float(d_active.evalf())) > 1e-12, (
        "the active component must vary with n, or claim 5 is vacuous")

    # --- 6. Hedges' variance weights against the effect it measures ------------------------
    # v(g) = 1/n + g^2/(2n) is increasing in |g|, so an inverse-variance weight is smaller for
    # the studies reporting larger effects and the pooled mean is pulled toward zero. This is
    # the -0.02 to -0.09 focus bias present in every arm table, textbook arm included.
    g = sp.Symbol("g", positive=True)
    v = 1 / n + g**2 / (2 * n)
    proof.claim(
        "the inverse-variance weight decreases in the observed magnitude",
        sp.simplify(sp.diff(1 / v, g) + (g / n) / v**2),
        r"\partial_g v^{-1} + (g/n)\,v^{-2}",
    )
    assert sp.simplify(sp.diff(1 / v, g).subs({g: 1, n: 30})) < 0, (
        "the weight must strictly decrease, or claim 6 is vacuous")
    return proof


if __name__ == "__main__":
    run().report()
