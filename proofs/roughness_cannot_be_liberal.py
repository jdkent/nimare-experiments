r"""A null that is rougher than the observed map makes the test conservative, never liberal.

The estimator's ``null_method`` docstring used to explain the family-wise inflation at small
collections this way: ``permute-images`` scatters each image's values, which destroys the spatial
autocorrelation, and the permuted statistic maps come out about 2.5 times rougher than the
observed one. The measurement is real. The inference from it was backwards, and running
``spatial-images`` -- which preserves the autocorrelation -- returned identical rates on all three
arms, so the explanation was already dead empirically.

This makes it dead in general, so the direction does not get re-proposed. The argument is that the
maximum of a field is *increasing* in the number of effectively independent locations it is
maximised over. A rougher null therefore has a larger maximum, a higher rejection threshold, and
fewer rejections.

Two versions, because the discrete one is airtight and the field one is the relevant geometry.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex import Proof  # noqa: E402

u, n, alpha, R = sp.symbols("u n alpha R", positive=True)

proof = Proof(
    "roughness_cannot_be_liberal",
    "A rougher null raises the maximum, so it cannot make a max-statistic test liberal",
    __doc__,
)

# ------------------------------------------------------------------ the discrete version

# For n independent standard normals the maximum has CDF Phi(u)^n. Work with the survival of the
# maximum, which is what a family-wise p-value is.
Phi = (1 + sp.erf(u / sp.sqrt(2))) / 2
survival = 1 - Phi**n
proof.define(r"P(\max > u)", survival)

# Increasing n strictly increases that survival at every u, because log Phi < 0.
d_survival = sp.diff(survival, n)
proof.claim(
    "more-locations-raise-the-survival-of-the-maximum",
    sp.simplify(d_survival + Phi**n * sp.log(Phi)),
    r"\partial_n P(\max > u) = -\Phi(u)^n \log \Phi(u) > 0",
)

# Hence the (1 - alpha) quantile of the maximum rises with n. Differentiate the defining relation
# Phi(q)^n = 1 - alpha implicitly rather than inverting it.
q = sp.Function("q")(n)
defining = Phi.subs(u, q) ** n - (1 - alpha)
implicit = sp.solve(sp.Eq(sp.diff(defining, n), 0), sp.Derivative(q, n))[0]
proof.define(r"dq/dn", sp.simplify(implicit))
# The derivative is -log(Phi) / (n * phi(q) / Phi(q)), a positive quantity over a positive one.
phi_q = sp.exp(-q**2 / 2) / sp.sqrt(2 * sp.pi)
proof.claim(
    "the-critical-value-rises-with-the-number-of-locations",
    sp.simplify(implicit - (-sp.log(Phi.subs(u, q)) * Phi.subs(u, q) / (n * phi_q))),
    r"dq/dn = -\log\Phi(q)\,\Phi(q) / (n\,\phi(q)) > 0",
)

# ------------------------------------------------------------------ the random-field version

# For a smooth Gaussian field the exceedance probability at a high threshold is the expected Euler
# characteristic, a sum over dimension of resel counts times EC densities. In a volume:
ln2 = sp.log(2)
density = {
    0: 1 - Phi,
    1: sp.sqrt(4 * ln2) / (2 * sp.pi) * sp.exp(-u**2 / 2),
    2: (4 * ln2) / (2 * sp.pi) ** sp.Rational(3, 2) * u * sp.exp(-u**2 / 2),
    3: (4 * ln2) ** sp.Rational(3, 2) / (2 * sp.pi) ** 2 * (u**2 - 1) * sp.exp(-u**2 / 2),
}
for d, rho in density.items():
    proof.define(rf"\rho_{d}(u)", rho)

# Above u = 1 every density is positive, so the expected EC is increasing in every resel count.
# The only one that can change sign is rho_3, through (u^2 - 1).
k3 = (4 * ln2) ** sp.Rational(3, 2) / (2 * sp.pi) ** 2
proof.claim(
    "the-three-dimensional-ec-density-factors-so-its-sign-is-visible",
    sp.simplify(density[3] - k3 * (u - 1) * (u + 1) * sp.exp(-u**2 / 2)),
    r"\rho_3(u) = k_3 (u-1)(u+1) e^{-u^2/2},\ k_3 > 0",
)
# Every factor but (u - 1) is positive, so rho_3 > 0 exactly above u = 1. A max-statistic
# threshold sits far above 1, so the sign question never arises in practice.
proof.claim(
    "and-vanishes-only-at-one",
    sp.simplify(density[3].subs(u, 1)),
    r"\rho_3(1) = 0",
)
proof.claim(
    "the-ec-density-increases-through-the-resel-count-linearly",
    sp.simplify(sp.diff(R * density[2], R) - density[2]),
    r"\partial_R \left(R\,\rho_d(u)\right) = \rho_d(u)",
)

# Resel counts scale as FWHM^{-d}, so a rougher field -- smaller FWHM -- has more resels in every
# dimension at once. Write the volume's resel count and differentiate in the smoothness.
fwhm, volume = sp.symbols("f V", positive=True)
resels = volume / fwhm**3
proof.claim(
    "roughening-a-field-adds-resels",
    sp.simplify(sp.diff(resels, fwhm) + 3 * volume / fwhm**4),
    r"\partial_f (V / f^3) = -3V/f^4 < 0,\ \text{so resels rise as } f \text{ falls}",
)

# Putting the pieces together: rougher null -> more resels -> larger expected EC at any threshold
# -> larger null maximum -> higher critical value -> fewer rejections. A permutation that
# roughens the map is conservative. Whatever makes the family-wise rate liberal at 12 studies, it
# is not this, and the measured equality of permute-images and spatial-images agrees.

if __name__ == "__main__":
    print(f"{len(proof.claims)} claims verified in {proof.name}")
    for label, shown in proof.claims:
        print(f"  {label}")
    print()
    print("Corollary: destroying spatial autocorrelation in the null cannot produce a liberal")
    print("family-wise rate. It produces a conservative one. The observed map is 2.5x smoother")
    print("than the permuted ones, which pushes the test the wrong way to explain 0.150.")
