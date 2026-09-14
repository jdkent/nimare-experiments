"""Does peak_bias leave z, p and prevalence alone, as the parameter docstring claims?

The claim matters for how much the peak_bias validity problem costs. If the correction only
moves the magnitude map, a reader who follows the advice to read `g_relative` is untouched by
it. If it moves `z`, it moves the inference.

A scalar rho and a per-study rho are different cases. A scalar multiplies every study's g by
the same factor and its variance by the square, so every inverse-variance weight is scaled
identically and the pooled z is unchanged. A per-study rho scales each study's variance by its
own factor, which reweights the studies -- so it can move the pooled estimate, its precision,
and the prevalence/magnitude split. Measured on a collection whose sample sizes span the range
that makes rho vary at all.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
mask = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)

for label, sample_size in (("N 15-400 (wide)", (15, 400)), ("N 20-40 (narrow)", (20, 40))):
    studyset = create_effect_size_coordinate_studyset(
        [(0.0, 0.0, 0.0)], effect_sizes=0.8, n_studies=24, sample_size=sample_size,
        tau=0.1, seed=11, n_noise_foci=3, noise_extent=30.0, spatial_sd=5.0)
    out = {}
    for name, peak_bias in (("none", None), ("per-study", "per-study"), ("scalar0.5", 0.5)):
        estimator = CBES(fwhm=12.0, mask=mask, peak_bias=peak_bias, null_method="none")
        result = estimator.fit(studyset)
        out[name] = {k: result.get_map(k, return_type="array").ravel()
                     for k in ("g", "z", "prevalence", "n_studies")}
        if peak_bias == "per-study":
            rho = estimator._peak_bias_
            spread = rho.max() / rho.min()
            print(f"\n{label}: rho_k min {rho.min():.3f} median {np.median(rho):.3f} "
                  f"max {rho.max():.3f}, spread {spread:.1f}x")
    covered = out["none"]["n_studies"] > 0
    print(f"{covered.sum()} covered voxels")
    print(f"{'variant':>10} {'map':>11} {'median ratio':>13} {'max rel diff':>13}")
    for name in ("per-study", "scalar0.5"):
        for key in ("g", "z", "prevalence"):
            a, b = out["none"][key][covered], out[name][key][covered]
            keep = np.isfinite(a) & np.isfinite(b) & (np.abs(a) > 1e-6)
            print(f"{name:>10} {key:>11} {np.median(b[keep] / a[keep]):13.4f} "
                  f"{np.max(np.abs(b[keep] - a[keep]) / np.abs(a[keep])):13.4f}", flush=True)
