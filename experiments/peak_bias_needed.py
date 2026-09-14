"""Does the censored likelihood already correct the winner's curse, leaving peak_bias idle?

The censored likelihood observes both outcomes -- a reported value, or a silence -- and those
two probabilities sum to one over the sample space, so threshold selection is already fully
modelled. If that works, `g` should be near the truth across sample sizes with peak_bias=None,
and any rho on top of it is a second charge for the same selection.

Two questions, measured against simulator truth:

  1. Homogeneous N, peak_bias=None: is `g` unbiased as N varies? This is the censored
     likelihood on its own. Note that per-study rho normalises to the median study, so on a
     homogeneous collection it is exactly 1.0 and there is nothing to compare against.
  2. Heterogeneous N: bias and MSE of `g` with peak_bias None against "per-study", which is
     where rho actually varies between studies and where the review found the damage.
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
N_SIMS = 30

shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
MASK = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
AT = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))


def one(seed, sample_size, peak_bias):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=TRUE_G, n_studies=24, sample_size=sample_size, tau=0.1,
        prevalence=1.0, seed=seed, n_noise_foci=3, noise_extent=30.0, spatial_sd=5.0)
    estimator = CBES(fwhm=12.0, mask=MASK, peak_bias=peak_bias, null_method="none")
    result = estimator.fit(studyset)
    g = result.get_map("g", return_type="array").ravel()[AT]
    rho = getattr(estimator, "_peak_bias_", None)
    spread = (rho.max() / rho.min()) if rho is not None else 1.0
    return float(g), float(spread)


def summarise(label, sample_size, peak_bias):
    rows = [one(seed, sample_size, peak_bias) for seed in range(N_SIMS)]
    g = np.array([r[0] for r in rows])
    spread = np.mean([r[1] for r in rows])
    bias = np.mean(g) - TRUE_G
    mse = np.mean((g - TRUE_G) ** 2)
    print(f"{label:>26} {str(peak_bias):>10} {spread:7.2f} {np.mean(g):8.3f} "
          f"{bias:+8.3f} {mse:8.4f}", flush=True)


print(f"true g = {TRUE_G} at the focus, {N_SIMS} simulations of 24 studies each\n")
print(f"{'collection':>26} {'peak_bias':>10} {'rho sp':>7} {'mean g':>8} {'bias':>8} {'MSE':>8}")
print("\n-- 1. the censored likelihood alone, as N varies --")
for n in (20, 40, 80, 200, 400):
    summarise(f"N = {n} (homogeneous)", (n, n), None)

print("\n-- 2. heterogeneous N, where rho actually varies --")
for label, sample_size in (("N 15-400", (15, 400)), ("N 20-40", (20, 40))):
    for peak_bias in (None, "per-study"):
        summarise(label, sample_size, peak_bias)
