"""Does the max-statistic guard actually fire in the cell it was built for?

`why_fwe_fails.py` established the diagnosis: twenty studies reporting two foci each give a
voxel-level familywise rate of 0.150 to 0.180 against a nominal 0.05, the generalized Pareto tail
is not implicated (identical rates with it off), and the cause is that across 200 permutations the
maximum statistic attains **six distinct values** with a coefficient of variation of 0.032.
Swapping two values inside a study does not move the maximum.

The guard was then rewritten to threshold what actually matters -- the number of distinct attained
maxima and their spread -- rather than the arrangement count, which waves 2^20 through. The unit
test pins the behaviour on a constructed null. What has *not* been measured is whether the guard
fires on the simulated cell that exposed the defect, which is the only evidence that matters for
the claim in the PR description.

Three numbers per cell, and they answer different questions:

  guard fired    the fraction of fits where voxel-level FWE was withheld. This is the fix
                 working as designed.
  voxel FWE      the familywise rejection rate over all fits. Withholding sets logp to 0, so a
                 withheld fit cannot reject and this should fall to near zero in the bad cell.
  rejected | kept   the rate among fits the guard let through. This is the number that decides
                 whether the guard is *sufficient*: if the fits it passes still over-reject, the
                 threshold is in the wrong place even though the mechanism is right.

The six-foci cell is the control. The guard must not fire there, and the rate must stay near
nominal -- a guard that refuses everything is not a fix.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

N_SIMS = int(os.environ.get("NSIMS", 60))
N_ITERS = int(os.environ.get("NITERS", 200))
ALPHA = 0.05
TRUTH = (0.0, 0.0, 0.0)


def build_mask():
    shape, step = (25, 25, 25), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, n_noise_foci):
    mask = build_mask()
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=0.0, n_studies=20, sample_size=(10, 1000),
        tau=0.0, prevalence=0.0, seed=seed, n_noise_foci=n_noise_foci, noise_extent=40.0)
    est = CBES(fwhm=10.0, mask=mask, peak_bias="per-study", n_iters=N_ITERS,
               seed=seed, cluster_threshold=None)
    result = est.fit(studyset)
    p = result.get_map("p", return_type="array").ravel()
    if np.all(p == 1.0):
        return None
    maps, _, _ = est.correct_fwe_montecarlo(result, vfwe_only=True)
    logp = maps["logp_level-voxel"].ravel()
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    # Withholding zeroes the map, so an all-zero logp over covered voxels is the guard's
    # signature. Checked rather than read off a flag so the measurement does not depend on an
    # attribute name staying put.
    withheld = bool(np.all(logp[covered] == 0.0))
    rejected = float((logp[covered] >= -np.log10(ALPHA)).any())
    return rejected, withheld


if __name__ == "__main__":
    se = np.sqrt(ALPHA * (1 - ALPHA) / N_SIMS)
    print(f"{N_SIMS} simulations, 20 studies, N 10-1000, global null, nominal {ALPHA}")
    print(f"binomial standard error about {se:.3f}; {N_ITERS} permutations\n")
    print(f"  {'foci/study':>11} {'guard fired':>12} {'voxel FWE':>10} "
          f"{'rejected|kept':>14} {'n kept':>7}")
    for n_foci in (2, 6):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(seed, n_foci) for seed in range(N_SIMS)) if r is not None]
        rej = np.array([r[0] for r in rows])
        wh = np.array([r[1] for r in rows], dtype=bool)
        kept = ~wh
        among = rej[kept].mean() if kept.any() else np.nan
        print(f"  {n_foci:>11} {wh.mean():12.3f} {rej.mean():10.3f} {among:14.3f} "
              f"{int(kept.sum()):7d}", flush=True)
    print("\nBefore the guard the two-foci cell rejected at 0.150 to 0.180. The fix works if that")
    print("cell's 'voxel FWE' is now near zero, and is sufficient only if 'rejected|kept' is")
    print("also near nominal. The six-foci control must show the guard not firing.")
