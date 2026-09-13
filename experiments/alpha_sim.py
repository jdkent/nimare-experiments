"""Does the flatter kernel discount survive on simulated data with a known answer?

The pain sweep favours w**0.25 to w**0.5 over the current w, but that is 21 studies of one
contrast scored against a criterion (the images) that is itself smooth -- exactly the kind of
thing a blurrier estimate flatters. Simulation has a known truth and a known null, so it can
say whether the gain is real and whether it costs false positive control.

Two questions, both of which have to come out right before the default moves:
  recovery  -- g at the true focus, and the prevalence there, with peaks selected from a field
  null      -- the uncorrected rejection rate when nothing is true anywhere
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

ALPHAS = (1.0, 0.5, 0.25)
N_REC, N_NULL = 6, 8
TRUE_G = 0.6

affine = np.array([[4., 0, 0, -40.], [0, 4., 0, -40.], [0, 0, 4., -40.], [0, 0, 0, 1.]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)
centre = tuple(mm2vox(np.array([[0., 0., 0.]]), affine)[0])

original = CBES._kernel_support


def with_alpha(alpha):
    def patched(self, sample_size=None, _a=alpha):
        offsets, values = original(self, sample_size=sample_size)
        return offsets, values**_a
    return patched


print(f"recovery: known g = {TRUE_G}, field-simulated peaks, {N_REC} sims")
print(f"{'alpha':>7s} {'g at truth':>11s} {'g_marginal':>11s} {'prevalence':>11s}", flush=True)
for alpha in ALPHAS:
    CBES._kernel_support = with_alpha(alpha)
    try:
        gs, gms, pis = [], [], []
        for seed in range(N_REC):
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=TRUE_G, n_studies=30, sample_size=(20, 40),
                n_noise_foci=2, noise_extent=36.0, seed=seed, simulate_field=True,
            )
            res = CBES(fwhm=12.0, mask=mask, null_method="none",
                       peak_bias="per-study").fit(ss)
            gs.append(res.get_map("g").get_fdata()[centre])
            gms.append(res.get_map("g_marginal").get_fdata()[centre])
            pis.append(res.get_map("prevalence").get_fdata()[centre])
        print(f"{alpha:7.2f} {np.mean(gs):11.3f} {np.mean(gms):11.3f} {np.mean(pis):11.3f}",
              flush=True)
    finally:
        CBES._kernel_support = original

print(f"\nglobal null, {N_NULL} sims (nominal 0.050 / 0.0010)")
print(f"{'alpha':>7s} {'unc<.05':>9s} {'unc<.001':>10s}", flush=True)
for alpha in ALPHAS:
    CBES._kernel_support = with_alpha(alpha)
    try:
        a, b = [], []
        for seed in range(N_NULL):
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
                prevalence=0.0, n_noise_foci=8, noise_extent=36.0, seed=seed,
                threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
            )
            res = CBES(fwhm=12.0, mask=mask, n_iters=200, seed=seed,
                       threshold="study-min", peak_bias="per-study").fit(ss)
            p = res.get_map("p", return_type="array")
            a.append(np.mean(p < 0.05))
            b.append(np.mean(p < 0.001))
        print(f"{alpha:7.2f} {np.mean(a):9.3f} {np.mean(b):10.4f}", flush=True)
    finally:
        CBES._kernel_support = original
