r"""Can retention be estimated at all, and from what kind of spread?

``censored.py`` takes retention as a supplied parameter, and
``retention_in_the_reporting_model.py`` shows why: from report rates at a *common* threshold and
precision, magnitude and retention have rank-one information and cannot separate. The design
document says the same thing from the other direction -- "a free per-study detection curve is
too flexible to learn from one table and can absorb the effect itself".

But a real corpus does not have one precision. This file asks what spread is needed to identify
the full triple :math:`(m, \tau^2, \rho)` from report indicators, and the answer has a clean
geometric form.

Write :math:`\sigma_i^2 = s_i^2 + \tau^2`, :math:`z_i = (c_i-m)/\sigma_i`,
:math:`S_i = 1-\Phi(z_i)`, and :math:`q_i = \rho S_i`. A Bernoulli report indicator contributes
:math:`\nabla q_i \nabla q_i^\top / [q_i(1-q_i)]`, so everything turns on the *direction* of

.. math:: \nabla q_i = \left(\frac{\rho\varphi_i}{\sigma_i},\
    \frac{\rho z_i \varphi_i}{2\sigma_i^2},\ S_i\right).

Claim 1 is the structural fact that makes the rest work: the first two components are
proportional, with factor :math:`z_i/(2\sigma_i)`. So the gradient is
:math:`A_i\,(1,\ u_i,\ v_i)` with :math:`A_i = \rho\varphi_i/\sigma_i`,
:math:`u_i = (c_i-m)/(2\sigma_i^2)` and :math:`v_i = \sigma_i/(\rho\,h_i)` for the hazard
:math:`h_i = \varphi_i/S_i`. A gradient's *direction* therefore depends on the study only through
its threshold and its precision -- not on anything else about it.

Claims 2-3 then give the identification condition exactly. The determinant of three studies'
gradients factors as :math:`A_1A_2A_3` times the determinant of

.. math:: \begin{pmatrix} 1 & 1 & 1 \\ u_1 & u_2 & u_3 \\ v_1 & v_2 & v_3\end{pmatrix},

which vanishes precisely when the three points :math:`(u_i, v_i)` are **collinear in the plane**.
So:

  * studies sharing a threshold *and* a precision give identical directions -- one point, rank
    one, no separation however many of them there are;
  * two distinct precisions give two points, hence rank at most two: enough for two parameters,
    never for three;
  * three parameters need three non-collinear points, so at least three distinct
    (threshold, precision) pairs, and collinearity is a real degeneracy rather than a
    measure-zero curiosity -- it is a single scalar equation on the design.

**What this deliberately does not claim.** Non-collinearity is local identification and nothing
more. The document is explicit that this "is not a guarantee of adequate finite-sample
precision", and the numeric section below measures the conditioning at realistic settings rather
than declaring victory: at a plausible spread of sample sizes the triple is formally identified
and the smallest eigenvalue is still four orders of magnitude below the largest, which is a
statement about how many studies it would take, not about whether it is possible.

**Not a restatement of the rank-one theorem.** ``mixture_identification.py`` proves the
two-dimensional determinant identity. This file does not reuse it: it computes a particular
three-dimensional determinant and factors it, which is a different statement in a different
dimension about a different parameter triple.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

m, c = sp.symbols("m c", real=True)
tau2 = sp.Symbol("tau^2", nonnegative=True)
s2 = sp.Symbol("s^2", positive=True)
rho = sp.Symbol("rho", positive=True)

sigma = sp.sqrt(s2 + tau2)
z = (c - m) / sigma
pdf = sp.exp(-z**2 / 2) / sp.sqrt(2 * sp.pi)
survival = (1 - sp.erf(z / sp.sqrt(2))) / 2
report_probability = rho * survival

proof = Proof(
    "what_identifies_retention",
    "What spread of thresholds and precisions identifies magnitude, heterogeneity and retention",
    __doc__,
)

gradient = sp.Matrix(
    [
        sp.diff(report_probability, m),
        sp.diff(report_probability, tau2),
        sp.diff(report_probability, rho),
    ]
)
proof.define(r"\nabla q", gradient)

# --------------------------------------------------- the gradient's three components

proof.claim(
    "the-magnitude-component-is-the-retention-weighted-density-over-the-scale",
    sp.simplify(gradient[0] - rho * pdf / sigma),
    r"\partial_m q = \frac{\rho\,\varphi(z)}{\sigma}",
)
proof.claim(
    "the-retention-component-is-the-survival-itself",
    sp.simplify(gradient[2] - survival),
    r"\partial_\rho q = S",
)
proof.claim(
    "and-the-heterogeneity-component-is-the-magnitude-one-scaled-by-z-over-two-sigma",
    sp.simplify(gradient[1] - gradient[0] * z / (2 * sigma)),
    r"\partial_{\tau^2} q = \frac{z}{2\sigma}\,\partial_m q",
)

# ------------------------------------- so the direction depends only on threshold and precision

amplitude = rho * pdf / sigma
hazard = pdf / survival
direction = sp.Matrix([1, (c - m) / (2 * (s2 + tau2)), sigma / (rho * hazard)])
proof.claim(
    "so-a-studys-gradient-is-one-amplitude-times-a-direction-fixed-by-its-threshold-and-precision",
    sp.simplify(gradient - amplitude * direction),
    r"\nabla q = A\,(1,\ u,\ v)^\top,\quad A = \frac{\rho\varphi}{\sigma},\ "
    r"u = \frac{c-m}{2\sigma^2},\ v = \frac{\sigma}{\rho h}",
)

# ------------------------------------------------- the three-study determinant, factored

A_1, A_2, A_3 = sp.symbols("A_1 A_2 A_3", positive=True)
u_1, u_2, u_3 = sp.symbols("u_1 u_2 u_3", real=True)
v_1, v_2, v_3 = sp.symbols("v_1 v_2 v_3", real=True)
stacked = sp.Matrix(
    [[A_1, A_2, A_3], [A_1 * u_1, A_2 * u_2, A_3 * u_3], [A_1 * v_1, A_2 * v_2, A_3 * v_3]]
)
shape = sp.Matrix([[1, 1, 1], [u_1, u_2, u_3], [v_1, v_2, v_3]])
proof.claim(
    "the-amplitudes-factor-straight-out-of-the-three-study-determinant",
    sp.simplify(stacked.det() - A_1 * A_2 * A_3 * shape.det()),
    r"\det[\nabla q_1\ \nabla q_2\ \nabla q_3] = A_1A_2A_3\,"
    r"\det\begin{pmatrix}1&1&1\\u_1&u_2&u_3\\v_1&v_2&v_3\end{pmatrix}",
)

# That shape determinant is twice the signed area of the triangle on the three (u, v) points, so
# it vanishes exactly on collinearity. Stating it as the cross product of two edges makes the
# geometry explicit rather than asserted.
edge_one = sp.Matrix([u_2 - u_1, v_2 - v_1])
edge_two = sp.Matrix([u_3 - u_1, v_3 - v_1])
proof.claim(
    "and-what-is-left-is-twice-the-triangle-area-so-it-vanishes-exactly-on-collinearity",
    sp.simplify(shape.det() - (edge_one[0] * edge_two[1] - edge_one[1] * edge_two[0])),
    r"\det(\cdot) = (u_2-u_1)(v_3-v_1) - (v_2-v_1)(u_3-u_1)",
)

# Two studies can never do it: a 3x3 matrix with a repeated column is singular, and two distinct
# directions span a plane. Written as the rank statement rather than a determinant.
repeated = sp.Matrix(
    [[A_1, A_2, A_2], [A_1 * u_1, A_2 * u_2, A_2 * u_2], [A_1 * v_1, A_2 * v_2, A_2 * v_2]]
)
proof.claim(
    "two-distinct-precisions-can-carry-at-most-two-of-the-three-parameters",
    sp.simplify(repeated.det()),
    r"\det[\nabla q_1\ \nabla q_2\ \nabla q_2] = 0",
)

if __name__ == "__main__":
    import numpy as np
    from scipy.stats import norm

    proof.report()

    # Conditioning at realistic settings: the document's own scalar parameters, with sample
    # sizes spread over a range a real corpus plausibly shows, and one common threshold.
    true_mean, between, retention, threshold = 0.4, 0.0225, 0.5, 0.6
    sample_sizes = np.array([12, 16, 20, 28, 40, 60, 90, 140])
    within = 1.0 / sample_sizes  # sampling variance of a one-sample standardised mean, roughly
    scales = np.sqrt(within + between)
    standardised = (threshold - true_mean) / scales
    densities = norm.pdf(standardised)
    survivals = norm.sf(standardised)
    probabilities = retention * survivals

    gradients = np.stack(
        [
            retention * densities / scales,
            retention * standardised * densities / (2 * scales**2),
            survivals,
        ],
        axis=1,
    )
    weights = 1.0 / (probabilities * (1 - probabilities))
    information = (gradients * weights[:, None]).T @ gradients
    eigenvalues = np.linalg.eigvalsh(information)

    print(f"  report probabilities run {probabilities.min():.3f} to {probabilities.max():.3f} "
          f"across n = {sample_sizes.min()}-{sample_sizes.max()}")
    print(f"  information eigenvalues per study: "
          + ", ".join(f"{value:.3g}" for value in eigenvalues))
    print(f"  condition number {eigenvalues[-1] / eigenvalues[0]:.3g}")
    if eigenvalues[0] <= 0:
        raise AssertionError(
            "the triple is not locally identified even with eight distinct precisions, which "
            "contradicts the collinearity claim and means one of the two is wrong"
        )
    print("  [ok] numeric: eight distinct precisions at a common threshold do identify the")
    print("       triple, so the (u, v) points are not collinear there")

    # And the degenerate case, as a check that the diagnostic can detect a failure at all.
    common = np.full(8, scales[3])
    standardised_common = (threshold - true_mean) / common
    flat = np.stack(
        [
            retention * norm.pdf(standardised_common) / common,
            retention * standardised_common * norm.pdf(standardised_common) / (2 * common**2),
            norm.sf(standardised_common),
        ],
        axis=1,
    )
    flat_information = flat.T @ flat
    flat_eigenvalues = np.linalg.eigvalsh(flat_information)
    if flat_eigenvalues[0] > 1e-10 * flat_eigenvalues[-1]:
        raise AssertionError("a common precision should give rank one and does not")
    print(f"  [ok] numeric: a common precision collapses to rank one "
          f"(smallest eigenvalue {flat_eigenvalues[0]:.1e} against "
          f"largest {flat_eigenvalues[-1]:.3g})")

    print()
    print("Consequence: retention is estimable in principle from precision spread alone, needing")
    print("at least three distinct (threshold, precision) pairs whose (u, v) points are not")
    print("collinear -- but the conditioning above says how expensive that is, and it is the")
    print(f"condition number {eigenvalues[-1] / eigenvalues[0]:.3g}, not the rank, that decides")
    print("whether a real corpus can do it. That is a simulation question, and it is the next")
    print("step rather than a conclusion drawn here.")
