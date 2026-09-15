"""An out-of-sample test of the prediction that a larger collection covers worse.

Coverage across all sixteen arms of the calibrated table is predicted, to a mean absolute error
of 0.034 and a maximum of 0.124, by nothing more than a shifted normal:

    coverage = Phi((1.96 se - b) / sd) - Phi((-1.96 se - b) / sd)

with `b` the bias and `sd` the estimator's actual spread (`experiments/coverage_is_noncentrality`).
There is no residual pathology to explain -- the interval fails exactly when and as much as the
bias-to-width ratio says it should.

That has an uncomfortable consequence, because `se` falls roughly as 1/sqrt(studies) while the
bias does not fall at all. Coverage therefore *decreases* as a collection grows, in every
configuration measured:

    coordinates only       0.75 at 12 studies  ->  0.35 at 24
    2 images, calibrated   0.99 at 12 studies  ->  0.91 at 24
    6 images, calibrated   0.98 at 12 studies  ->  0.92 at 24

So the usual reassurance -- more studies, better inference -- is false for this interval, and the
calibrated configuration's good coverage at twelve studies may be a small-collection property
rather than a validated operating point. The extrapolation says 48 studies at the same image
fraction as the 2-of-24 arm should land near 0.85, and 96 near 0.75.

This arm tests that directly. The image fraction is held at 1/12 so the weight share does not
move, and only the collection grows: 12 studies with 1 image, 24 with 2, 48 with 4. If coverage
tracks the prediction, the honest statement about the interval is conditional on collection size,
and the PR should say so. If coverage holds up instead, the extrapolation is wrong and the
normal-shift model breaks somewhere between 24 and 48 studies -- also worth knowing, because it
is the model everything else here is being read through.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from joblib import Parallel, delayed
from scipy.stats import norm
import interval_coverage as bed

N = int(os.environ.get("NSIMS", 40))
CAL = dict(peak_bias="per-study", peak_bias_scale="images")

if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; image fraction held at 1/12; {N} replications")
    print(f"{'studies':>8s} {'images':>7s} {'mean g':>8s} {'bias':>8s} {'mean se':>8s} "
          f"{'sd':>7s} {'se/sd':>6s} {'cover':>6s} {'predicted':>10s}")
    for n_studies in (12, 24, 48):
        n_image = n_studies // 12
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(bed.one)(s, n_studies, n_image, 0.0, CAL["peak_bias"], 10.0,
                             CAL["peak_bias_scale"]) for s in range(N)) if r is not None]
        if not rows:
            print(f"{n_studies:8d} {n_image:7d}   no usable fits"); continue
        g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
        ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
        b = g[ok].mean() - bed.TRUE_G
        sd = g[ok].std(ddof=1)
        half = 1.96 * se[ok].mean()
        cover = np.mean((g[ok] - 1.96*se[ok] <= bed.TRUE_G)
                        & (bed.TRUE_G <= g[ok] + 1.96*se[ok]))
        pred = norm.cdf((half - b) / sd) - norm.cdf((-half - b) / sd)
        print(f"{n_studies:8d} {n_image:7d} {g[ok].mean():8.3f} {b:+8.3f} "
              f"{se[ok].mean():8.3f} {sd:7.3f} {se[ok].mean()/sd:6.2f} {cover:6.2f} "
              f"{pred:10.2f}  (n={ok.sum()})", flush=True)
    print("\nThe prediction column is the shifted normal evaluated on this arm's own bias and")
    print("spread, so it checks the model's shape rather than its extrapolation. The claim being")
    print("tested is the trend across rows: coverage should fall as the collection grows.")
