"""Is the common peak-height scale a constant, or does it depend on the effect it is scaling?

Reporting density turned out to move the required common scale only a little: 0.796 at 5 peaks
per study, 0.818 at 15, 0.837 at 28, all at a true g of 0.8. If that were the only axis, a
default near 0.81 would be defensible.

The axis that decides it is the true effect size. A peak-height correction derived at mu = 0
assumes the reported height is noise-dominated. A large true effect clears the threshold
easily and is barely selected, so it needs almost no correction; a marginal one is heavily
selected and needs a lot. If the required scale swings with the truth, then no fixed default
is right and the constant genuinely cannot be set without external information -- which is
what `peak_bias_scale` already says, and would mean the honest answer is to leave it alone.

Fixed at 5 peaks per study, the realistic end, and reported as the scale that would put the
estimate exactly on the truth.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox
from peak_bias_realistic_peaks import cap_peaks

TRUTH = (0.0, 0.0, 0.0)
N_SIMS = 100
MAX_PEAKS = 5

shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
MASK = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
AT = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))


def one(seed, true_g):
    """One simulation, or NaN where the effect is too small to be reported at all.

    At a small true effect nothing clears the threshold inside the mask and the fit raises.
    That is a fact about coordinate meta-analysis rather than an error, so it is counted.
    """
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=true_g, n_studies=24, sample_size=(20, 40), tau=0.1,
        prevalence=1.0, seed=seed, simulate_field=True, smoothness_fwhm=10.0, field_zooms=4.0)
    studyset = cap_peaks(studyset, MAX_PEAKS)
    estimator = CBES(fwhm=12.0, mask=MASK, peak_bias=None, null_method="none")
    try:
        result = estimator.fit(studyset)
    except ValueError:
        return np.nan, np.nan
    return (float(result.get_map("g", return_type="array").ravel()[AT]),
            float(estimator.peak_information_["excess_z"]))


print(f"{N_SIMS} simulations, 24 studies, N 20-40, {MAX_PEAKS} peaks/study, peak_bias=None\n")
print(f"{'true g':>7} {'usable':>7} {'excess z':>9} {'mean g':>8} {'bias':>8} {'ratio':>7} "
      f"{'scale needed':>13}")
for true_g in (0.2, 0.4, 0.6, 0.8, 1.2, 1.6):
    rows = [one(seed, true_g) for seed in range(N_SIMS)]
    g = np.array([r[0] for r in rows])
    excess = np.array([r[1] for r in rows])
    g, excess = g[np.isfinite(g)], excess[np.isfinite(excess)]
    if g.size < 10:
        print(f"{true_g:7.1f} {g.size:7d}   too few studies report anything to estimate")
        continue
    se = g.std(ddof=1) / np.sqrt(g.size)
    print(f"{true_g:7.1f} {g.size:7d} {excess.mean():+9.3f} {g.mean():8.3f} "
          f"{g.mean() - true_g:+8.3f} {g.mean() / true_g:7.3f} "
          f"{true_g / g.mean():8.3f}+-{true_g * se / g.mean() ** 2:.3f}", flush=True)
