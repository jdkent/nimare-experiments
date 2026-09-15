"""Error rates under the within-analysis null, at the condition that broke the old one.

The across-table shuffle rejected at 96.7% against a nominal 5% once sample sizes ranged
widely: a value carries its study's precision, so moving it to another study's voxel built a
null whose spread came from the roster rather than from the observed map. The within-analysis
shuffle keeps each value inside its own study. Two things to establish:

  * the false-positive rate is back to nominal under heterogeneous N, not merely better;
  * what the restriction costs in power, since a study reporting one focus now contributes no
    randomness and even a 6-focus study contributes only 720 arrangements.

A global null (effect 0, prevalence 0) makes every rejection a false one, and the familywise
rate is the share of *maps* carrying one. Power is read at the
truth voxel of a focal effect. Run on a 25-voxel box at 4mm so 100 simulations per cell is
affordable; the null is per-voxel, so a smaller brain changes the multiplicity but not the
per-voxel rate.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

N_SIMS = 100
N_ITERS = 200
ALPHA = 0.05
TRUTH = (0.0, 0.0, 0.0)


def build_mask():
    shape, step = (25, 25, 25), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, sample_size, n_studies, effect_size, n_noise_foci):
    mask = build_mask()
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=effect_size, n_studies=n_studies, sample_size=sample_size,
        tau=0.0, prevalence=0.0 if effect_size == 0.0 else 1.0, seed=seed,
        n_noise_foci=n_noise_foci, noise_extent=40.0)
    estimator = CBES(fwhm=10.0, mask=mask, peak_bias="per-study", n_iters=N_ITERS,
                     seed=seed, cluster_threshold=None)
    result = estimator.fit(studyset)
    p = result.get_map("p", return_type="array").ravel()
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    if not covered.any() or np.all(p == 1.0):
        # p == 1 everywhere means the collection admitted too few arrangements to test; that
        # is a refusal, not a result, and has to be counted separately rather than as 0% FPR.
        return np.nan, np.nan, np.nan, np.nan, bool(np.all(p == 1.0))
    maps, _, _ = estimator.correct_fwe_montecarlo(result, vfwe_only=True)
    logp = maps["logp_level-voxel"].ravel()
    reject_fwe = logp >= -np.log10(ALPHA)
    if effect_size == 0.0:
        # The familywise rate is the probability of *any* false rejection, not the average
        # share of voxels rejected. An earlier version returned the share, which sits near zero
        # even in a map that does reject somewhere, so it could not show whether the test held
        # its level -- it read 0.0000 where the familywise rate was about 0.10.
        return (float(np.mean(p[covered] <= ALPHA)), float(reject_fwe[covered].any()),
                np.nan, np.nan, False)
    # Power is read at the truth, which is the voxel a reader cares about.
    at = np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), mask.affine)[0], mask.shape)
    flat = np.zeros(mask.shape, dtype=bool)
    flat.ravel()[at] = True
    idx = int(np.flatnonzero(flat.ravel()[np.ravel(np.ones(mask.shape, bool))])[0])
    return (np.nan, np.nan, float(p[idx] <= ALPHA), float(reject_fwe[idx]), False)


print(f"{N_SIMS} simulations per cell, {N_ITERS} permutations, nominal alpha {ALPHA}")
print(f"binomial standard error on a rate across {N_SIMS} draws is about "
      f"{np.sqrt(ALPHA * (1 - ALPHA) / N_SIMS):.3f}\n")

CELLS = [
    ("global null, N 20-40",      (20, 40),   20, 0.0, 6),
    ("global null, N 10-1000",    (10, 1000), 20, 0.0, 6),
    ("global null, N 10-1000",    (10, 1000), 12, 0.0, 6),
    ("global null, 2 foci/study", (10, 1000), 20, 0.0, 2),
    ("power g=0.8, N 20-40",      (20, 40),   30, 0.8, 6),
    ("power g=0.8, 2 foci/study", (20, 40),   30, 0.8, 2),
]
print(f"{'cell':>28s} {'studies':>7s} {'uncorrected':>14s} {'voxel FWE':>12s} {'refused':>8s}")
for label, sample_size, n_studies, effect, n_noise in CELLS:
    rows = Parallel(n_jobs=4)(
        delayed(one)(seed, sample_size, n_studies, effect, n_noise) for seed in range(N_SIMS))
    refused = float(np.mean([r[4] for r in rows]))
    if effect == 0.0:
        unc = np.array([r[0] for r in rows], dtype=float)
        fwe = np.array([r[1] for r in rows], dtype=float)
        unc_s = f"{np.nanmean(unc):.4f}" if np.isfinite(unc).any() else "n/a"
        fwe_s = f"{np.nanmean(fwe):.4f}" if np.isfinite(fwe).any() else "n/a"
    else:
        unc = np.array([r[2] for r in rows], dtype=float)
        fwe = np.array([r[3] for r in rows], dtype=float)
        unc_s = f"{np.nanmean(unc):.3f} pw" if np.isfinite(unc).any() else "n/a"
        fwe_s = f"{np.nanmean(fwe):.3f} pw" if np.isfinite(fwe).any() else "n/a"
    print(f"{label:>28s} {n_studies:7d} {unc_s:>14s} {fwe_s:>12s} {refused:>8.2f}", flush=True)
