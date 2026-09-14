"""Finding 7, measured on a simulator that actually has peak-height inflation.

An earlier version of this used the default point-based simulator, where a study's reported
value is drawn at a location rather than being the local maximum of a smooth field. That
simulator has no peak-height inflation at all, so it cannot test a peak-height correction --
it can only show that correcting for an absent bias hurts, which is close to tautological.
``simulate_field=True`` reports genuine local maxima of smooth fields, and its reported heights
run about +0.2 z above the null peak expectation, so there is something for rho to recover.

Three questions:

  1. With peak_bias=None, is `g` biased *up* once the reported values are real maxima? That is
     what decides whether any peak-height correction is warranted, over and above the censoring
     the likelihood already does.
  2. Does peak_bias="per-study" reduce that bias, or overshoot it?
  3. Does it help or hurt MSE, separately on collections where rho varies little and a lot?
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
TRUE_G = 0.8
N_SIMS = 100

shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
MASK = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
AT = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))


def one(seed, sample_size, peak_bias):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=TRUE_G, n_studies=24, sample_size=sample_size, tau=0.1,
        prevalence=1.0, seed=seed, simulate_field=True, smoothness_fwhm=10.0,
        field_zooms=4.0)
    estimator = CBES(fwhm=12.0, mask=MASK, peak_bias=peak_bias, null_method="none")
    result = estimator.fit(studyset)
    g = result.get_map("g", return_type="array").ravel()[AT]
    rho = getattr(estimator, "_peak_bias_", None)
    spread = float(rho.max() / rho.min()) if rho is not None else 1.0
    excess = estimator.peak_information_["excess_z"]
    return float(g), spread, float(excess)


def summarise(label, sample_size, peak_bias):
    rows = [one(seed, sample_size, peak_bias) for seed in range(N_SIMS)]
    g = np.array([r[0] for r in rows])
    keep = np.isfinite(g)
    g = g[keep]
    bias = g.mean() - TRUE_G
    mse = np.mean((g - TRUE_G) ** 2)
    # Monte Carlo error on each, so a difference can be read against its own precision.
    se_bias = g.std(ddof=1) / np.sqrt(g.size)
    se_mse = np.std((g - TRUE_G) ** 2, ddof=1) / np.sqrt(g.size)
    print(f"{label:>18} {str(peak_bias):>10} {np.mean([r[1] for r in rows]):6.2f} "
          f"{np.mean([r[2] for r in rows]):+7.3f} {g.mean():7.3f} "
          f"{bias:+7.3f}+-{se_bias:.3f} {mse:7.4f}+-{se_mse:.4f}", flush=True)


print(f"true g = {TRUE_G}, {N_SIMS} simulations of 24 studies, smooth fields, local maxima\n")
print(f"{'collection':>18} {'peak_bias':>10} {'rho sp':>6} {'excess':>7} {'mean g':>7} "
      f"{'bias':>14} {'MSE':>15}")
for label, sample_size in (("N 20-40", (20, 40)), ("N 15-400", (15, 400)), ("N=30 fixed", (30, 30))):
    for peak_bias in (None, "per-study"):
        summarise(label, sample_size, peak_bias)
