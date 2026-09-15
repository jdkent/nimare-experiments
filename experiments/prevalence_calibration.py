"""What does `prevalence` converge to, and is its compression invertible?

The map is documented as compressed toward the middle of its range -- a true 0.25 reads 0.49 to
0.60 -- and read ordinally for that reason. But "compressed" is not a characterisation. If the
relationship between the truth and the estimate is affine with stable coefficients, it can be
inverted and the map becomes usable as a fraction, which matters because the decision-relevant
quantity for planning a study is prevalence on a real scale: power under zero inflation is capped
at `pi + (1 - pi) * alpha`, so below about 0.8 no sample size reaches 80%.

Three things are separated here that the earlier measurement did not.

**Which truth.** A study can have an effect and still fail to report it, so the fraction of
studies *with* an effect at a voxel and the fraction that *report* near it are different
quantities, the second strictly smaller. Both are recorded, so which one the estimate tracks is
readable rather than assumed.

**The functional form.** A fine grid of true prevalence rather than two points, so an affine fit
is identifiable and its residuals visible.

**Stability.** The coefficients are estimated at two sample-size ranges and two effect sizes. A
correction is only worth having if it does not depend on things the analyst cannot know.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
N_SIMS = 12
GRID = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


def build_mask():
    shape, step = (21, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, prevalence, sample_size, effect):
    mask = build_mask()
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=effect, n_studies=24, sample_size=sample_size,
        tau=0.0, prevalence=prevalence, seed=seed, n_noise_foci=3, noise_extent=40.0)
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none")
    result = est.fit(studyset)
    at = mm2vox(np.asarray([TRUTH]), mask.affine)[0]
    flat = int(np.ravel_multi_index(tuple(at), mask.shape))
    pi_hat = float(result.get_map("prevalence", return_type="array").ravel()[flat])
    # The fraction of analyses that actually placed a focus within 10 mm of the truth: what a
    # convergence measure sees, and strictly at most the fraction that have an effect.
    coords = studyset.coordinates
    xyz = coords[["x", "y", "z"]].values.astype(float)
    near = np.linalg.norm(xyz - np.asarray(TRUTH), axis=1) <= 10.0
    n_analyses = len(set(coords["id"]))
    reporting = len(set(coords.loc[near, "id"])) / max(n_analyses, 1)
    return pi_hat, reporting


print(f"{N_SIMS} simulations per cell, 24 studies, prevalence read at the truth voxel\n")
print(f"  {'N range':>10} {'effect':>7} {'true pi':>8} {'reporting':>10} {'pi_hat':>8} "
      f"{'sd':>6}")
rows = {}
for sample_size, effect in (((20, 40), 0.8), ((20, 40), 0.4), ((10, 200), 0.8)):
    for prevalence in GRID:
        out = Parallel(n_jobs=4)(
            delayed(one)(seed, prevalence, sample_size, effect) for seed in range(N_SIMS))
        pis = np.array([o[0] for o in out], dtype=float)
        reps = np.array([o[1] for o in out], dtype=float)
        rows.setdefault((sample_size, effect), []).append(
            (prevalence, np.nanmean(reps), np.nanmean(pis)))
        print(f"  {str(sample_size):>10} {effect:7.1f} {prevalence:8.2f} "
              f"{np.nanmean(reps):10.3f} {np.nanmean(pis):8.3f} {np.nanstd(pis):6.3f}",
              flush=True)
    print()

print(f"  {'N range':>10} {'effect':>7} {'vs true pi':>28} {'vs reporting':>28}")
print(f"  {'':>10} {'':>7} {'slope':>9} {'intercept':>10} {'r':>7} "
      f"{'slope':>9} {'intercept':>10} {'r':>7}")
for key, vals in rows.items():
    truth_pi = np.array([v[0] for v in vals])
    reporting = np.array([v[1] for v in vals])
    pi_hat = np.array([v[2] for v in vals])
    a = stats.linregress(truth_pi, pi_hat)
    b = stats.linregress(reporting, pi_hat)
    print(f"  {str(key[0]):>10} {key[1]:7.1f} {a.slope:9.3f} {a.intercept:10.3f} "
          f"{a.rvalue:7.3f} {b.slope:9.3f} {b.intercept:10.3f} {b.rvalue:7.3f}")
print("\n  A stable slope and intercept against one of the two truths would make the"
      "\n  compression invertible, and prevalence usable as a fraction rather than ordinally.")
