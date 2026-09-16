"""When does a biased channel stop paying for itself?

Two channels estimate the same `mu`. The images are unbiased with variance `v_i`; the coordinates
carry a fixed bias `b` and variance `v_c`. Inverse-variance weighting puts weight `w` on the
coordinates. The question the pain image-count sweep raised is when the combination beats the
images alone -- and the answer has to be free of `w`, because `w` is chosen by the same
information that sets the variances, so a condition still containing it would say nothing.

Verified here rather than asserted, because the cancellation is the whole content of the result:
every term in the information ratio drops out and what survives is a statement about the variance
of the *difference* between the channels, which is what makes the criterion a Hausman statistic
and computable without knowing the truth.

The last two claims are the ones that matter in practice: that `Q < 1` is exactly the condition,
and that `Q` can be rebuilt from what the estimator reports -- the full estimate, the images-only
estimate, the reported standard error and the information share -- with no access to either
channel in isolation.
"""
import sympy as sp

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof

mu, b, v_i, v_c, w, T = sp.symbols("mu b v_i v_c w T", positive=True)
g_i, g_c = sp.symbols("g_i g_c", real=True)

proof = Proof(
    "channel_bias_crossing",
    "When a biased channel stops paying for itself",
    r"An unbiased channel of variance $v_i$ against a biased one of variance $v_c$ and bias $b$.",
)

# The inverse-variance weight on the coordinate channel.
w_star = proof.define(r"w^\*", v_i / (v_i + v_c))

# Mean squared error of the combination `(1 - w) g_i + w g_c`, where `g_i` is unbiased with
# variance `v_i` and `g_c` has mean `mu + b` and variance `v_c`.
mse_combined = proof.define(
    r"\operatorname{MSE}(w)", (w * b) ** 2 + (1 - w) ** 2 * v_i + w**2 * v_c
)
mse_images = proof.define(r"\operatorname{MSE}_{\text{img}}", v_i)

proof.claim(
    "weight-is-the-minimiser-when-unbiased",
    sp.simplify(sp.solve(sp.diff(mse_combined.subs(b, 0), w), w)[0] - w_star),
    r"At $b = 0$ the MSE-minimising weight is $w^\* = v_i/(v_i + v_c)$.",
)

# The excess MSE at the inverse-variance weight, which is what the estimator actually uses.
excess = sp.simplify(mse_combined.subs(w, w_star) - mse_images)

# Every term in the information ratio cancels: the excess is negative exactly when
# `b**2 < v_i + v_c`, so the criterion does not involve the weight at all.
proof.claim(
    "excess-factorises",
    sp.simplify(sp.factor(excess) - v_i**2 * (b**2 - v_i - v_c) / (v_i + v_c) ** 2),
    r"$\operatorname{MSE}(w^\*) - v_i = v_i^2\,(b^2 - v_i - v_c)/(v_i + v_c)^2$.",
)

# The prefactor is positive for all positive variances, so the sign of the excess is the sign of
# `b**2 - v_i - v_c`. Checked as a property of the symbols rather than at a chosen point: a claim
# that substitutes numbers and then verifies a sign proves nothing about the general case.
prefactor = v_i**2 / (v_i + v_c) ** 2
assert prefactor.is_positive is True, "the prefactor's sign is not settled by positivity alone"
proof.claim(
    "criterion-is-weight-free",
    sp.simplify(sp.factor(excess) - prefactor * (b**2 - v_i - v_c)),
    r"$\operatorname{MSE}(w^\*) - v_i$ is a positive multiple of $b^2 - v_i - v_c$, so the "
    r"combination helps iff $b^2 < v_i + v_c$ -- the variance of the channels' difference, "
    r"with no $w$ in it.",
)

# Hausman form: `b` is the mean of `g_c - g_i`, whose variance is `v_i + v_c`.
q = proof.define(r"Q", (g_c - g_i) ** 2 / (v_i + v_c))
proof.claim(
    "q-is-the-criterion",
    sp.simplify(q.subs(g_c - g_i, b) - b**2 / (v_i + v_c)),
    r"With $g_c - g_i$ at its mean $b$, $Q = b^2/(v_i+v_c)$, so $Q < 1$ is the criterion.",
)

# What the estimator reports: the combined estimate, the images-only estimate, the total
# information `T = 1/se**2`, and the information share `w`. Note `w` is the share of the
# *information*, which for inverse-variance weighting is the same number as the weight.
delta = proof.define(r"\Delta", w * (g_c - g_i))
proof.claim(
    "q-from-reported-maps",
    sp.simplify(
        (delta**2 * T * (1 - w) / w)
        - q.subs({v_i: 1 / ((1 - w) * T), v_c: 1 / (w * T)})
    ),
    r"$\Delta^2\,T\,(1-w)/w = Q$ when $v_i = 1/((1-w)T)$ and $v_c = 1/(wT)$, so $Q$ follows "
    r"from $g$, $g$ with selection off, $se$ and the information share alone.",
)

# The crossing image count, from `v_i = 1/(k n)`: solving `b**2 = v_i + v_c` for `k`.
k, n = sp.symbols("k n", positive=True)
proof.claim(
    "crossing-image-count",
    sp.simplify(sp.solve(sp.Eq(b**2, 1 / (k * n) + v_c), k)[0] - 1 / (n * (b**2 - v_c))),
    r"$k^\* = 1/\!\left(n\,(b^2 - v_c)\right)$, which is infinite -- the coordinates never "
    r"hurt -- whenever $v_c \ge b^2$.",
)

# Anti-vacuity: the excess must actually change sign, or every claim above is about nothing.
assert sp.simplify(excess.subs({b: 2, v_i: 1, v_c: 1})) > 0, "excess never turns positive"
assert sp.simplify(excess.subs({b: 0, v_i: 1, v_c: 1})) < 0, "excess never turns negative"
# And the reported-map form must not be trivially zero.
assert sp.simplify(delta.subs({w: sp.Rational(1, 2), g_c: 1, g_i: 0})) != 0, "delta is vacuous"

if __name__ == "__main__":
    proof.report()
