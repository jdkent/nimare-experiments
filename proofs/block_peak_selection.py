r"""A block's reported maximum: the likelihood, and what the height adds over the indicator.

Step 2 of the design document's build order (its section 6.3): "small spatial blocks with exact
peak-selection inequalities". The scalar reference in ``censored.py`` treats each observation as
an interval about one estimate; a block treats a group of correlated elements together, and the
observation is *the largest of them, if it cleared the threshold*.

Model, which is the document's own section 8.2 generative model:

.. math::
    Y_{ij} = m + U_i + \epsilon_{ij}, \quad U_i \sim N(0,\tau^2),\ \epsilon_{ij}\sim N(0,\sigma^2),

for :math:`j = 1..M` elements of one block. The reporting rule retains
:math:`H_i = \max_j Y_{ij}` when :math:`H_i > c`.

**On the never-cap rule.** Retaining one maximum *per block* is not a cap on a study's peak
count: a whole-brain model has many blocks, each contributing its own maximum-or-silence, and the
study's total is whatever survives. It becomes a cap only when the block is the whole brain,
which is what the document's section 8.2 does and says it does -- "this experiment intentionally
studies a one-maximum cap; it is not claimed to reproduce whole-brain peak extraction". Any
measurement here that uses a single block per study is reproducing that stated cap and is
labelled as such, and no such measurement is evidence about anything else.

What is proved:

  * the block is equicorrelated, with :math:`\operatorname{corr} = \tau^2/(\tau^2+\sigma^2)`
    (claims 1-2), derived from the model by ``sympy.stats`` rather than asserted -- the simplest spatial structure that is not independence, and the one the
    document's experiment has;
  * conditionally on the study effect the elements are independent, so the maximum's conditional
    CDF is a power and its density follows by differentiation (claims 2-3). The unconditional
    versions are integrals over :math:`U` that sympy will not do with a symbolic exponent, so
    they are definitions here and are checked by quadrature below, against the partition
    :math:`\int_c^\infty f_H + P(H \le c) = 1`;
  * **the reported height is a location family in the mean**, so
    :math:`\partial_m \log f_H = -\partial_h \log f_H` (claim 4). That is why heights carry
    information about the magnitude at all, and it turns the height's Fisher information into an
    ordinary location-family one;
  * the score of a reported height, and of a silent block (claims 5-6);
  * as the threshold falls to :math:`-\infty` every block reports and the likelihood becomes the
    plain maximum's density (claim 7), so the selection model nests the unselected one.

**A tool defect found on the way.** sympy 1.14 returns 1 for
``limit((1 + erf((c-m-u)/s))/2, c, -oo)`` while correctly returning -1 for the erf alone;
substituting :math:`c = -10^6` gives :math:`1.2\times10^{-125}`. The nesting claim is therefore
stated on the erf, which the tool gets right, with the remaining arithmetic step in its
statement. A claim whose verification runs through that limit would pass or fail for reasons
unconnected to its truth.

**What this file measures rather than proves.** The document's stated finding is that "the extra
gain from heights beyond the correct reporting indicator is modest" (.050 against .053, .010
against .011, .022 against .024 across its three regimes). That is an information comparison, and
the numeric section computes both informations directly at its settings. The claim being modest
is a property of *those* settings -- a high threshold, nine elements, one block -- and the
computation reports how it moves with the threshold rather than asserting it generally.
"""
import os
import sys

import sympy as sp
import sympy.stats as stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

m, h, u, c = sp.symbols("m h u c", real=True)
tau2, sigma2 = sp.symbols("tau^2 sigma^2", positive=True)
count = sp.Symbol("M", positive=True, integer=True)

proof = Proof(
    "block_peak_selection",
    "The likelihood of a block's reported maximum, and what its height adds",
    __doc__,
)


def normal_pdf(value, mean, sd):
    return sp.exp(-((value - mean) ** 2) / (2 * sd**2)) / (sd * sp.sqrt(2 * sp.pi))


def normal_cdf(value, mean, sd):
    return (1 + sp.erf((value - mean) / (sd * sp.sqrt(2)))) / 2


# ------------------------------------------------------- the block's covariance structure

# Derived from the model rather than asserted: two elements of one block share the study effect
# and nothing else, so sympy.stats should return the between-study variance as their covariance.
shared = stats.Normal("U", 0, sp.sqrt(tau2))
first_noise = stats.Normal("e_1", 0, sp.sqrt(sigma2))
second_noise = stats.Normal("e_2", 0, sp.sqrt(sigma2))
first, second = m + shared + first_noise, m + shared + second_noise
proof.claim(
    "two-elements-of-a-block-covary-exactly-by-the-between-study-variance",
    sp.simplify(stats.covariance(first, second) - tau2),
    r"\operatorname{Cov}(Y_{ij}, Y_{ik}) = \tau^2,\quad j \ne k",
)
proof.claim(
    "so-the-block-is-equicorrelated-at-the-between-study-share-of-the-total-variance",
    sp.simplify(
        stats.covariance(first, second) / stats.variance(first) - tau2 / (tau2 + sigma2)
    ),
    r"\operatorname{corr}(Y_{ij}, Y_{ik}) = \frac{\tau^2}{\tau^2+\sigma^2}",
)

# --------------------------------------- the conditional maximum, where everything is exact

element_cdf = normal_cdf(h, m + u, sp.sqrt(sigma2))
conditional_max_cdf = element_cdf**count
proof.define(r"P(H \le h \mid U=u)", conditional_max_cdf)
proof.claim(
    "conditional-independence-makes-the-maximums-cdf-a-power-of-the-elements",
    sp.simplify(conditional_max_cdf - sp.prod([element_cdf for _ in range(1)]) ** count),
    r"P(H \le h \mid U) = \Phi\!\left(\tfrac{h-m-U}{\sigma}\right)^M",
)
conditional_max_pdf = sp.diff(conditional_max_cdf, h)
proof.claim(
    "so-its-density-is-m-times-the-cdf-to-the-m-minus-one-times-the-element-density",
    sp.simplify(
        conditional_max_pdf - count * element_cdf ** (count - 1) * normal_pdf(h, m + u, sp.sqrt(sigma2))
    ),
    r"f_{H\mid U}(h) = M\,\Phi\!\left(\tfrac{h-m-U}{\sigma}\right)^{M-1}"
    r"\frac{1}{\sigma}\varphi\!\left(\tfrac{h-m-U}{\sigma}\right)",
)

# ------------------------------------- the height is a location family in the mean

# Every appearance of m in the conditional density is inside (h - m - u), so differentiating in
# m is differentiating in h with the sign flipped. The same then holds after integrating over U,
# because the integral is over u and m still enters only through h - m - u.
proof.claim(
    "the-reported-height-is-a-location-family-so-its-mean-score-is-minus-its-height-score",
    sp.simplify(
        sp.diff(sp.log(conditional_max_pdf), m) + sp.diff(sp.log(conditional_max_pdf), h)
    ),
    r"\partial_m \log f_H = -\partial_h \log f_H",
)

# --------------------------------------------------------------- the two score terms

# A reported block contributes the log density of its maximum; the selection event adds nothing
# further, because the height already implies it cleared the threshold.
proof.claim(
    "a-reported-blocks-score-is-the-elementwise-hazard-sum-scaled-by-the-count",
    sp.simplify(
        sp.diff(sp.log(conditional_max_pdf), m)
        - (
            (h - m - u) / sigma2
            - (count - 1)
            * normal_pdf(h, m + u, sp.sqrt(sigma2))
            / normal_cdf(h, m + u, sp.sqrt(sigma2))
        )
    ),
    r"\partial_m \log f_{H\mid U} = \frac{h-m-u}{\sigma^2} "
    r"- (M-1)\,\frac{\varphi}{\Phi}\!\left(\tfrac{h-m-u}{\sigma}\right)\frac{1}{\sigma}",
)

silence = normal_cdf(c, m + u, sp.sqrt(sigma2)) ** count
proof.claim(
    "and-a-silent-blocks-score-is-m-times-its-own-elements-inverse-mills-ratio",
    sp.simplify(
        sp.diff(sp.log(silence), m)
        + count
        * normal_pdf(c, m + u, sp.sqrt(sigma2))
        / normal_cdf(c, m + u, sp.sqrt(sigma2))
    ),
    r"\partial_m \log P(H \le c \mid U) = -M\,\frac{\varphi}{\Phi}"
    r"\!\left(\tfrac{c-m-u}{\sigma}\right)\frac{1}{\sigma}",
)

# ---------------------------------------------- the selection model nests the unselected one

# Stated on the erf rather than on the CDF. sympy 1.14 returns 1 for
# ``limit((1 + erf((c-m-u)/s))/2, c, -oo)`` while returning -1 for the erf alone, and
# substituting c = -1e6 into the CDF gives 1.2e-125. The claim below is the part the tool
# computes correctly; Phi = (1 + erf)/2 -> 0 and Phi**M -> 0 then follow by arithmetic. Do not
# "simplify" this back to a limit of the CDF: it will pass for the wrong reason.
proof.claim(
    "sending-the-threshold-down-drives-the-elements-cdf-to-zero-so-every-block-reports",
    sp.simplify(sp.limit(sp.erf((c - m - u) / (sp.sqrt(sigma2) * sp.sqrt(2))), c, -sp.oo) + 1),
    r"\lim_{c\to-\infty}\operatorname{erf}\!\left(\tfrac{c-m-u}{\sigma\sqrt2}\right) = -1"
    r"\ \Rightarrow\ \Phi^M \to 0",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.integrate import quad
    from scipy.stats import norm

    proof.report()

    # The document's section 8.2 settings.
    true_mean, between_sd, element_sd, elements, cut = 0.4, 0.15, 0.2, 9, 0.75

    def marginal_no_report(mean):
        return quad(
            lambda draw: norm.cdf((cut - mean - draw) / element_sd) ** elements
            * norm.pdf(draw, 0, between_sd),
            -8 * between_sd,
            8 * between_sd,
        )[0]

    def marginal_max_pdf(height, mean):
        return quad(
            lambda draw: elements
            * norm.cdf((height - mean - draw) / element_sd) ** (elements - 1)
            * norm.pdf((height - mean - draw) / element_sd)
            / element_sd
            * norm.pdf(draw, 0, between_sd),
            -8 * between_sd,
            8 * between_sd,
        )[0]

    # The partition the symbolic claims could not reach: the reported density and the silence
    # probability must exhaust the outcomes.
    reported_mass = quad(lambda height: marginal_max_pdf(height, true_mean), cut, 6.0)[0]
    silent_mass = marginal_no_report(true_mean)
    if abs(reported_mass + silent_mass - 1.0) > 1e-6:
        raise AssertionError(
            f"the reporting outcomes do not exhaust the probability: "
            f"{reported_mass:.6f} + {silent_mass:.6f} = {reported_mass + silent_mass:.6f}"
        )
    print(f"  [ok] numeric: P(report) = {reported_mass:.4f} and P(silent) = {silent_mass:.4f} "
          f"sum to one")

    # What the height adds over the indicator, which is what the document calls modest. Both are
    # Fisher informations about the mean, per study, at the truth.
    step = 1e-4

    def indicator_information(mean):
        probability = 1.0 - marginal_no_report(mean)
        slope = (
            (1.0 - marginal_no_report(mean + step)) - (1.0 - marginal_no_report(mean - step))
        ) / (2 * step)
        return slope**2 / (probability * (1 - probability))

    def height_information(mean):
        """Expected squared score of the full observation: the height where reported, the
        indicator where not. The height's own term is a location-family information."""
        def contribution(height):
            density = marginal_max_pdf(height, mean)
            if density <= 0:
                return 0.0
            slope = (
                marginal_max_pdf(height, mean + step) - marginal_max_pdf(height, mean - step)
            ) / (2 * step)
            return slope**2 / density
        reported = quad(contribution, cut, 6.0, limit=200)[0]
        silent_probability = marginal_no_report(mean)
        silent_slope = (
            marginal_no_report(mean + step) - marginal_no_report(mean - step)
        ) / (2 * step)
        return reported + silent_slope**2 / silent_probability

    print(f"{'threshold':>10} {'P(report)':>10} {'I indicator':>12} {'I with height':>14} "
          f"{'se ratio':>9}")
    for threshold in (0.0, 0.4, 0.75, 1.1):
        cut = threshold
        indicator = indicator_information(true_mean)
        full = height_information(true_mean)
        if full < indicator - 1e-8:
            raise AssertionError(
                "adding the height reduced the information, which no observation can do"
            )
        print(f"{threshold:>10.2f} {1 - marginal_no_report(true_mean):>10.4f} "
              f"{indicator:>12.4f} {full:>14.4f} {np.sqrt(indicator / full):>9.4f}")
    cut = 0.75

    print()
    print("Reading: the 'se ratio' column is how much narrower the interval on the mean gets by")
    print("reading the height as well as the report. At the document's threshold of .75 it is")
    print("close to one, which is its 'modest' -- and the column shows that is a property of a")
    print("high threshold, not of heights: as the cut falls, the reported height starts to")
    print("carry most of the information because the indicator stops discriminating.")
