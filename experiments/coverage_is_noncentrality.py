"""Can a user predict this interval's coverage from quantities the fit reports?

Coverage of `g +/- 1.96 se` against a truth offset by a bias b is
`Phi(1.96 - b/se) - Phi(-1.96 - b/se)` if g is normal around truth+b with sd = se. Two things
can break that: the reported `se` may not equal the actual sd of g (se/sd != 1), and g may not
be normal. Using the *actual* sd rather than the reported se isolates the first.

If the non-centrality predicts coverage, then "how many studies do I have" is the wrong question
and "how large is my bias relative to my interval" is the right one -- and the table's most
practical fact follows: more studies shrink se without shrinking bias, so coverage falls.
"""
import numpy as np
from scipy.stats import norm

# label, bias, mean se, sd of g, measured coverage
ROWS = [
    ("12 st,  0 img",            +0.255, 0.171, 0.080, 0.75),
    ("12 st,  0 img fwhm 16",    +0.248, 0.109, 0.064, 0.28),
    ("12 st,  0 img fwhm 24",    +0.249, 0.086, 0.057, 0.05),
    ("12 st,  2 img cal",        -0.038, 0.107, 0.081, 0.99),
    ("12 st,  2 img None",       +0.127, 0.134, 0.098, 0.87),
    ("12 st,  6 img cal",        -0.026, 0.083, 0.066, 0.98),
    ("12 st,  6 img None",       +0.027, 0.090, 0.071, 0.99),
    ("12 st, 12 img cal",        -0.018, 0.065, 0.059, 0.94),
    ("24 st,  0 img",            +0.246, 0.118, 0.059, 0.35),
    ("24 st,  0 img fwhm 16",    +0.243, 0.077, 0.043, 0.03),
    ("24 st,  0 img fwhm 24",    +0.244, 0.061, 0.038, 0.00),
    ("24 st,  2 img cal",        -0.064, 0.077, 0.057, 0.91),
    ("24 st,  2 img None",       +0.163, 0.102, 0.069, 0.58),
    ("24 st,  6 img cal",        -0.049, 0.069, 0.048, 0.92),
    ("24 st,  6 img None",       +0.073, 0.080, 0.058, 0.89),
    ("24 st, 24 img cal",        -0.022, 0.045, 0.041, 0.97),
]

print(f"{'arm':26s} {'b/se':>6s} {'pred(se)':>9s} {'b/sd':>6s} {'pred(sd)':>9s} "
      f"{'actual':>7s} {'err(sd)':>8s}")
errs_se, errs_sd = [], []
for label, b, se, sd, cov in ROWS:
    # Half-width is fixed at 1.96*se; what varies is the sd the estimate actually has.
    half = 1.96 * se
    p_se = norm.cdf(1.96 - b / se) - norm.cdf(-1.96 - b / se)
    p_sd = norm.cdf((half - b) / sd) - norm.cdf((-half - b) / sd)
    errs_se.append(p_se - cov); errs_sd.append(p_sd - cov)
    print(f"{label:26s} {b/se:6.2f} {p_se:9.2f} {b/sd:6.2f} {p_sd:9.2f} {cov:7.2f} "
          f"{p_sd - cov:+8.2f}")
for name, e in (("se only", errs_se), ("se for width, sd for spread", errs_sd)):
    e = np.array(e)
    print(f"\n{name:30s} mean abs error {np.abs(e).mean():.3f}  max {np.abs(e).max():.3f}")
print("\npred(se) treats the reported se as the truth about both width and spread; pred(sd)")
print("uses the reported se for the width the user gets and the actual sd for the spread.")
