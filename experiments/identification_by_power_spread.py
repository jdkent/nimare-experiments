"""Does separating prevalence from magnitude depend on the spread of study power? (E3)

Section 18 of the first-principles notes recasts the coordinate observation as an occupancy
record with imperfect detection: a study reports near a voxel with probability
`pi * D(mu, n, u)`. In that model occupancy and detection are separately identified only through
repeat visits or through covariates that move detection without moving occupancy. Studies are the
repeat visits; sample size and reporting threshold are the detection covariates.

So the prediction is specific. With a homogeneous roster -- every study the same size at the same
threshold -- only the product `pi * mu` is identified and the split into factors is set by the
likelihood's shape rather than by the data. As the spread of `n` and `u` grows, the factors should
become estimable and the fitted prevalence should start tracking the truth.

Three rosters over the same truth, with true prevalence swept so that "tracks the truth" can be
measured as a slope rather than asserted from one point:

  fixed      n = 30 for every study, u = 3.2905 for every study
  n varies   n ~ U(15, 120), u = 3.2905
  both vary  n ~ U(15, 120), u drawn from the thresholds papers actually use

Scored by the slope of fitted prevalence on true prevalence (1.0 would be unbiased tracking, 0.0
no information), and by the same slope for the product. If the theory holds the prevalence slope
rises across the three rosters while the product's slope stays roughly constant -- the product
being identified all along.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

SITE = (0.0, 0.0, 0.0)
TRUE_MU = 0.6
TRUE_PI_SWEEP = (0.25, 0.50, 0.75, 1.00)
N_SIMS = int(os.environ.get("NSIMS", 24))
N_STUDIES = 24
LOCALISATION_SD = 4.0
#: Thresholds in real use: FDR-ish, p<0.005, p<0.001, and a stricter corrected height.
THRESHOLDS = (2.5758, 3.0902, 3.2905, 3.7190)
ROSTERS = ("fixed", "n varies", "both vary")


def build_mask():
    shape, step = (21, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def draw_roster(rng, roster):
    if roster == "fixed":
        return [(30, 3.2905)] * N_STUDIES
    if roster == "n varies":
        return [(int(rng.integers(15, 121)), 3.2905) for _ in range(N_STUDIES)]
    return [(int(rng.integers(15, 121)), float(rng.choice(THRESHOLDS)))
            for _ in range(N_STUDIES)]


def one(seed, roster, true_pi):
    rng = np.random.default_rng(seed)
    mask = build_mask()
    studies = []
    for k, (n, u) in enumerate(draw_roster(rng, roster)):
        meta = {"sample_sizes": [n]}
        points = []
        if rng.random() < true_pi:
            var = 1.0 / n + TRUE_MU**2 / (2.0 * n)
            z = rng.normal(TRUE_MU, np.sqrt(var)) * np.sqrt(n)
            if abs(z) >= u:
                loc = np.asarray(SITE) + rng.normal(0, LOCALISATION_SD, 3)
                points.append((loc, z))
        # One noise focus, just clearing this study's own cut with the overshoot of a
        # smooth field, so the roster's threshold spread shows up in the table as it would.
        loc = rng.uniform(-36, 36, 3)
        points.append((loc, (u + rng.exponential(1.0 / u)) * rng.choice([-1.0, 1.0])))
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(c) for c in loc],
                 "values": [{"kind": "Z", "value": float(z)}]} for loc, z in points]}]})
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none")
    res = est.fit(Studyset({"id": "i", "name": "i", "studies": studies},
                           target=None, mask=mask))
    at = mm2vox(np.asarray([SITE]), mask.affine)[0]
    pos = int(np.ravel_multi_index(tuple(at), mask.shape))
    g = float(res.get_map("g", return_type="array").ravel()[pos])
    pi = float(res.get_map("prevalence", return_type="array").ravel()[pos])
    return g, pi, g * pi


if __name__ == "__main__":
    print(f"true mu {TRUE_MU}, {N_STUDIES} studies, {N_SIMS} replications per cell")
    print(f"true marginal at each prevalence: "
          + " ".join(f"{TRUE_MU*p:.3f}" for p in TRUE_PI_SWEEP) + "\n")
    for roster in ROSTERS:
        pis, gs, margs = [], [], []
        for true_pi in TRUE_PI_SWEEP:
            rows = [r for r in Parallel(n_jobs=8)(
                delayed(one)(s, roster, true_pi) for s in range(N_SIMS)) if r is not None]
            a = np.array(rows, dtype=float)
            ok = np.isfinite(a).all(axis=1)
            gs.append(a[ok, 0].mean()); pis.append(a[ok, 1].mean()); margs.append(a[ok, 2].mean())
        truth = np.array(TRUE_PI_SWEEP)
        pi_slope = stats.linregress(truth, pis).slope
        marg_slope = stats.linregress(truth * TRUE_MU, margs).slope
        print(f"--- roster: {roster} ---")
        print("   true prevalence  " + " ".join(f"{v:7.2f}" for v in truth))
        print("  fitted prevalence " + " ".join(f"{v:7.3f}" for v in pis)
              + f"   slope on truth {pi_slope:+.3f}")
        print("  fitted g          " + " ".join(f"{v:7.3f}" for v in gs))
        print("  fitted g_marginal " + " ".join(f"{v:7.3f}" for v in margs)
              + f"   slope on truth {marg_slope:+.3f}", flush=True)
        print()
    print("Prediction: the prevalence slope rises across the three rosters as the spread of")
    print("power grows, while the marginal slope stays roughly constant -- the product being")
    print("identified all along and the factors only becoming so with detection covariates.")
