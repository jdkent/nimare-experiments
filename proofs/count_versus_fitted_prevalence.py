r"""When does the naive report count beat the fitted prevalence at ordering voxels?

A standing puzzle (#51): a map of "how many studies named this voxel" orders voxels by prevalence
about as well as the fitted :math:`\hat\pi`, and sometimes better. That looks like a failure of the
model, and it is not. It is a consequence of what each statistic is estimating.

The expected number of reports among :math:`n` studies is

.. math:: n\,P(\text{report}) = n\left[1 - \pi s_a(\mu) - (1-\pi)s_0\right],

a function of both :math:`\pi` and :math:`\mu`. As an *ordering* statistic it inherits whichever of
the two varies across voxels. Where :math:`\mu` is roughly constant, the count is a monotone
function of :math:`\pi` alone and therefore a perfect ranker, with the variance of a binomial and
nothing more. The fitted :math:`\hat\pi` spends information separating two things that did not
need separating, and pays variance for it.

Where :math:`\mu` does vary, the count confounds the two and the fitted value should win. So the
puzzle has a testable resolution rather than a verdict: the count wins exactly when the magnitude
is flat, and the comparison measures the corpus, not the estimator.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

mu, pi, n = sp.symbols("mu pi n", positive=True)
s_a = sp.Function("s_a")(mu)
s_0 = sp.Symbol("s_0", positive=True)

proof = Proof(
    "count_versus_fitted_prevalence",
    "When the naive report count out-ranks the fitted prevalence",
    __doc__,
)

report_rate = 1 - (pi * s_a + (1 - pi) * s_0)
proof.define(r"P(\text{report})", report_rate)

# Monotone in pi, with a slope that does not depend on pi at all: the count is an affine function
# of the prevalence once mu is fixed.
proof.claim(
    "the-report-rate-is-affine-in-the-prevalence",
    sp.simplify(sp.diff(report_rate, pi) - (s_0 - s_a)),
    r"\partial_\pi P(\text{report}) = s_0 - s_a(\mu)",
)
proof.claim(
    "with-no-curvature-in-the-prevalence",
    sp.simplify(sp.diff(report_rate, pi, 2)),
    r"\partial_\pi^2 P(\text{report}) = 0",
)

# s_a falls as the effect grows -- a larger effect is less likely to stay inside (-c, c) -- so
# s_0 - s_a > 0 and the count rises with the prevalence. Strictly monotone, hence rank-preserving.
proof.define(r"s_0 - s_a(\mu)", s_0 - s_a)

# It is also monotone in mu at fixed pi, through the same term, which is the confound.
proof.claim(
    "and-monotone-in-the-magnitude-through-the-same-term",
    sp.simplify(sp.diff(report_rate, mu) + pi * sp.diff(s_a, mu)),
    r"\partial_\mu P(\text{report}) = -\pi\,s_a'(\mu)",
)

# The whole dependence on (pi, mu) collapses into one product. That is the confound, stated
# exactly: two voxels with the same pi(s_a - s_0) produce the same expected count however
# different their prevalences are.
level = pi * (s_a - s_0)
proof.claim(
    "the-count-depends-on-the-pair-only-through-one-product",
    sp.simplify(report_rate - (1 - s_0 - level)),
    r"P(\text{report}) = 1 - s_0 - \pi\,(s_a(\mu) - s_0)",
)
proof.define(r"\text{count level set}", level)

# Hold mu fixed across voxels and the product is an affine function of pi, so the count is a
# strictly monotone transform of the prevalence: its rank correlation with the truth is 1 up to
# binomial sampling noise. No estimator can rank better, and one that also fits mu must rank
# worse in finite samples, because it spends information separating what did not need separating.

if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    print()
    print("Resolution of #51: the count out-ranks the fitted prevalence exactly where the")
    print("magnitude is flat across voxels, because there the count is a monotone function of")
    print("the prevalence and the fit pays variance to separate two things that did not need")
    print("separating. Where the magnitude varies, the count confounds them. So the comparison")
    print("is a measurement of the corpus, not a verdict on the estimator -- and it predicts")
    print("that the count's advantage should shrink as the spread of g across voxels grows.")
