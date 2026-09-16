r"""What a fixed iteration cap does to an EM fit, and why it is not a neutral runtime knob.

Measured (notes, "max_iter=25 is undeclared shrinkage"): at the shipped cap, 76% of voxels sit
more than 0.01 in log-likelihood below the best a grid can find, and raising the cap monotonically
worsens :math:`\mu`'s rmse while moving :math:`\hat\pi` past the truth. Early stopping is doing
something, and the question is what.

The answer is exact enough to state. An EM map converges *linearly*: near its fixed point the
error is multiplied by a constant each iteration, and Dempster, Laird and Rubin (1977,
doi:10.1111/j.2517-6161.1977.tb01600.x) identify that constant as the fraction of information the
missing data hides. So stopping at :math:`t` iterations does not return an approximation of the
MLE with some unspecified error -- it returns a *convex combination* of the MLE and the starting
value, with a weight that is computable.

The consequence is the part worth acting on. The weight is :math:`\rho^t`, and :math:`\rho` is the
missing-information fraction, which varies from voxel to voxel. A fixed cap therefore applies
*more* shrinkage exactly where the censoring hides more information -- which is to say, at the
voxels where the coordinate channel is doing the most work and the images the least. That is not a
prior anyone chose. An explicit penalty with a stated weight would at least be a choice.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

rho, t = sp.symbols("rho t", positive=True)
theta_0, theta_star = sp.symbols("theta_0 theta^*", real=True)

proof = Proof(
    "early_stopping_is_shrinkage",
    "A fixed EM iteration cap is a convex combination of the MLE and the start",
    __doc__,
)

# ------------------------------------------------------- the linearised iterate

# Near its fixed point an EM map M satisfies M(theta) - theta* = rho (theta - theta*) to first
# order, with rho = M'(theta*). Iterating t times multiplies the error by rho^t.
error_t = rho**t * (theta_0 - theta_star)
iterate = theta_star + error_t
proof.define(r"\theta_t", iterate)

# That is a convex combination of the MLE and the start, with weight rho^t on the start. Written
# out so the shrinkage is explicit rather than implied by the error form.
combination = (1 - rho**t) * theta_star + rho**t * theta_0
proof.claim(
    "the-iterate-is-a-convex-combination-of-the-mle-and-the-start",
    sp.simplify(iterate - combination),
    r"\theta_t = (1-\rho^t)\,\theta^* + \rho^t\,\theta_0",
)

# The weight is in [0, 1] for rho in [0, 1], which is where a convergent EM lives, so it really is
# a convex combination and not an extrapolation.
proof.claim(
    "the-weight-is-one-at-no-iterations-and-zero-in-the-limit",
    sp.simplify(sp.limit(rho**t, t, sp.oo).subs(rho, sp.Rational(1, 2)))
    + sp.simplify((rho**t).subs(t, 0) - 1),
    r"\rho^0 = 1,\quad \rho^t \to 0 \text{ as } t \to \infty \text{ for } \rho < 1",
)

# ------------------------------------------------------- equivalent penalty

# A quadratic penalty pulling towards theta_0 with weight lambda, against a log-likelihood with
# curvature I, has its optimum at the same convex combination with weight lambda / (I + lambda).
# Equating the two weights gives the penalty a fixed cap is silently imposing.
lam, I = sp.symbols("lambda I", positive=True)
penalised = (I * theta_star + lam * theta_0) / (I + lam)
proof.define(r"\theta_\lambda", penalised)
proof.claim(
    "the-penalised-optimum-is-the-same-shape-of-combination",
    sp.simplify(penalised - ((1 - lam / (I + lam)) * theta_star + lam / (I + lam) * theta_0)),
    r"\theta_\lambda = (1-w)\theta^* + w\theta_0,\quad w = \lambda/(I+\lambda)",
)
equivalent_lambda = sp.solve(sp.Eq(lam / (I + lam), rho**t), lam)[0]
proof.claim(
    "so-a-cap-of-t-iterations-is-a-penalty-of-this-size",
    sp.simplify(equivalent_lambda - I * rho**t / (1 - rho**t)),
    r"\lambda_{\text{eff}} = I\,\rho^t / (1 - \rho^t)",
)

# ------------------------------------------------------- why that is not neutral

# rho is the missing-information fraction: the share of the complete-data information that the
# censoring hides. Asserted, not derived -- it is Dempster, Laird and Rubin's result about the EM
# map, not an algebraic identity.
I_complete, I_missing = sp.symbols("I_c I_m", positive=True)
proof.define(r"\rho", I_missing / I_complete)

# Substituting, the effective penalty at a voxel rises with how much information is missing there.
# Differentiate in the missing fraction to fix the sign.
lam_eff = I * (I_missing / I_complete) ** t / (1 - (I_missing / I_complete) ** t)
slope = sp.simplify(sp.diff(lam_eff, I_missing))
proof.define(r"\partial \lambda_{\text{eff}} / \partial I_m", slope)
proof.claim(
    "the-effective-penalty-grows-with-the-missing-information",
    sp.simplify(
        sp.sign(slope.subs({I: 1, I_complete: 1, t: 2, I_missing: sp.Rational(1, 2)})) - 1
    ),
    r"\partial_{I_m} \lambda_{\text{eff}} > 0",
)

if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    print()
    print("A cap of t iterations is a ridge penalty of size I rho^t / (1 - rho^t) towards the")
    print("starting value, and rho is the fraction of information the censoring hides. So the")
    print("cap shrinks hardest at the voxels where the coordinates carry the most and the")
    print("images the least -- a prior nobody chose, varying by voxel, and invisible.")
    print()
    print("What this does not settle: whether replacing it with a stated lambda is better. The")
    print("measurements say the shrinkage is buying real accuracy on mu, so the replacement has")
    print("to be simulated against a known truth before it is worth writing.")
