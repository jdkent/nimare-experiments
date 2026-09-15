"""Why is the familywise rate 0.18 when studies report two foci each?

The corrected error-rate run found the within-analysis null anticonservative in exactly the
regime cluster-extent reporting produces: 20 studies, 2 foci apiece, voxel-level familywise rate
0.18 against a nominal 0.05, with the degenerate-collection guard silent because 2^20
arrangements clears its 10^4 threshold.

Two mechanisms could do that and they call for different fixes.

**A narrow null.** Swapping two values inside a study barely perturbs the map, so the permutation
distribution of the maximum is concentrated, its 95th percentile sits close to the typical
observed maximum, and too many maps clear it. The fix is a guard on the spread of the attained
maxima rather than on the arrangement count.

**The tail fit.** Corrected p-values are read off a generalized Pareto fitted to the upper tail of
the maximum-statistic null. Fitted to a null that is nearly discrete -- few distinct attained
maxima -- the extrapolation can run well below the empirical tail and return p-values that are
too small. The fix is to refuse the fit, not the collection.

They are separated by rerunning the same simulations with ``tail_approximation`` off, which falls
back to the empirical tail. If the rate drops to nominal the tail fit did it; if it stays near
0.18 the null itself is too narrow; if it lands between, both contribute.

Reports the spread of the null maxima alongside, since that is the quantity a better guard would
have to threshold.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

N_SIMS = 60
N_ITERS = 200
ALPHA = 0.05
TRUTH = (0.0, 0.0, 0.0)


def build_mask():
    shape, step = (25, 25, 25), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, n_noise_foci, tail):
    mask = build_mask()
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=0.0, n_studies=20, sample_size=(10, 1000),
        tau=0.0, prevalence=0.0, seed=seed, n_noise_foci=n_noise_foci, noise_extent=40.0)
    est = CBES(fwhm=10.0, mask=mask, peak_bias="per-study", n_iters=N_ITERS,
               seed=seed, cluster_threshold=None)
    result = est.fit(studyset)
    p = result.get_map("p", return_type="array").ravel()
    if np.all(p == 1.0):
        return np.nan, np.nan, np.nan
    maps, _, _ = est.correct_fwe_montecarlo(result, vfwe_only=True,
                                            tail_approximation=tail)
    logp = maps["logp_level-voxel"].ravel()
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    rejected = float((logp[covered] >= -np.log10(ALPHA)).any())
    # The null maxima themselves: how much does shuffling actually move the maximum?
    null = est.null_distributions_.get("values_level-voxel")
    if null is None:
        for key, value in est.null_distributions_.items():
            if "voxel" in key and np.ndim(value) == 1 and np.size(value) > 10:
                null = value
                break
    if null is None:
        return rejected, np.nan, np.nan
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    spread = float(np.std(null) / max(np.mean(null), 1e-9))
    distinct = float(np.unique(np.round(null, 4)).size)
    return rejected, spread, distinct


print(f"{N_SIMS} simulations, 20 studies, N 10-1000, global null, nominal {ALPHA}")
print(f"binomial standard error is about {np.sqrt(ALPHA*(1-ALPHA)/N_SIMS):.3f}\n")
print(f"  {'foci/study':>11} {'tail fit':>9} {'voxel FWE':>10} "
      f"{'null max CV':>12} {'distinct maxima':>16}")
for n_foci in (2, 6):
    for tail in (True, False):
        rows = Parallel(n_jobs=4)(
            delayed(one)(seed, n_foci, tail) for seed in range(N_SIMS))
        rate = np.nanmean([r[0] for r in rows])
        cv = np.nanmean([r[1] for r in rows])
        distinct = np.nanmean([r[2] for r in rows])
        print(f"  {n_foci:>11} {str(tail):>9} {rate:10.3f} {cv:12.3f} {distinct:16.0f}",
              flush=True)
print("\nIf turning the tail fit off restores nominal, the extrapolation is the fault and the"
      "\nremedy is to refuse the fit. If it does not, the null itself is too narrow and the"
      "\nguard has to threshold the spread of the attained maxima.")
