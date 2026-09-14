"""What a mis-inferred reporting threshold actually costs the fit.

cluster_extent_reporting.py established the error: inverting the order statistic recovers the
height threshold to +0.06 z when peaks really were selected by height, but overshoots by +0.41
at a 10-voxel extent threshold and +1.10 at 50 voxels, because extent thresholding keeps broad
clusters (whose peaks are higher) and discards isolated low ones. What that error is worth is a
separate question, and the one that matters: the threshold feeds both rho_k, which divides the
reported effect sizes, and the censoring cutoff, which decides how surprising silence is.

Measured by handing the estimator a threshold that is deliberately wrong by a known amount, on
simulated data with a known effect size and prevalence, so both outputs can be scored.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
TRUE_U = 3.2905
TRUE_EFFECT = 0.8
TRUE_PREVALENCE = 0.6
N_REPLICATES = 10


def focus_index(masker):
    ijk = mm2vox(np.array([TRUTH]), masker.mask_img.affine)[0]
    mask = np.asarray(masker.mask_img.dataobj).astype(bool)
    flat = np.full(mask.shape, -1, dtype=np.int64)
    flat[mask] = np.arange(mask.sum())
    return int(flat[tuple(ijk)])


print(f"true threshold z = {TRUE_U}, effect {TRUE_EFFECT}, prevalence {TRUE_PREVALENCE}, "
      f"{N_REPLICATES} replicates\n")
print(f"{'threshold given':>17s} {'error':>7s} {'g at focus':>11s} {'g/true':>8s} "
      f"{'pi at focus':>12s}")
for offset in (0.0, 0.41, 0.80, 1.10, -0.40):
    g_vals, pi_vals = [], []
    for seed in range(N_REPLICATES):
        studyset = create_effect_size_coordinate_studyset(
            [TRUTH],
            effect_sizes=TRUE_EFFECT,
            n_studies=30,
            sample_size=(20, 40),
            tau=0.1,
            prevalence=TRUE_PREVALENCE,
            threshold_z=TRUE_U,
            seed=seed,
            n_noise_foci=2,
            noise_extent=40.0,
        )
        estimator = CBES(
            fwhm=10.0,
            null_method="none",
            peak_bias="per-study",
            threshold=float(TRUE_U + offset),
        )
        result = estimator.fit(studyset)
        index = focus_index(studyset.masker)
        if index < 0:
            continue
        g_vals.append(float(result.get_map("g", return_type="array").ravel()[index]))
        pi_vals.append(float(result.get_map("prevalence", return_type="array").ravel()[index]))
    g_mean = float(np.nanmean(g_vals))
    print(f"{TRUE_U + offset:17.2f} {offset:+7.2f} {g_mean:11.3f} "
          f"{g_mean / TRUE_EFFECT:8.2f} {np.nanmean(pi_vals):12.3f}", flush=True)
