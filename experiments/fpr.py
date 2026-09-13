"""False positive rate of CBES under a global null: foci are noise, everywhere.

Every study reports peaks, but none of them marks a real effect. A valid estimator should
flag ~5% of voxels at an uncorrected p < .05, and should produce *any* surviving voxel in
<=5% of whole simulations once corrected.
"""
import sys, warnings
import numpy as np
import nibabel as nib
from nimare.correct import FDRCorrector, FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

warnings.simplefilter("ignore")

N_SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 50
MODELS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["none", "tobit", "zero-inflated"]

affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)

def null_studyset(seed):
    """30 studies, 8 noise foci each, no true effect anywhere."""
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
        prevalence=0.0, n_noise_foci=8, noise_extent=36.0, seed=seed,
    )

print(f"global null, {N_SIMS} simulations, 30 studies x 8 noise foci\n")
print(f"{'model':16s} {'unc p<.05':>10s} {'unc p<.001':>11s} {'FWE-bonf':>9s} "
      f"{'FDR q=.05':>10s} {'FWE-mc':>8s}")

for model in MODELS:
    frac_05, frac_001 = [], []
    any_bonf, any_fdr, any_mc = [], [], []
    for seed in range(N_SIMS):
        ss = null_studyset(seed)
        est = CBES(fwhm=12.0, mask=mask, selection_model=model)
        res = est.fit(ss)
        p = res.get_map("p", return_type="array")
        frac_05.append(np.mean(p < 0.05))
        frac_001.append(np.mean(p < 0.001))

        bonf = FWECorrector(method="bonferroni").transform(res)
        any_bonf.append(
            np.any(bonf.get_map("p_corr-FWE_method-bonferroni", return_type="array") < 0.05))
        fdr = FDRCorrector(method="indep", alpha=0.05).transform(res)
        any_fdr.append(
            np.any(fdr.get_map("p_corr-FDR_method-indep", return_type="array") < 0.05))

        if seed < min(N_SIMS, 20):  # montecarlo is the expensive one
            mc = FWECorrector(method="montecarlo", n_iters=100).transform(res)
            any_mc.append(
                np.any(mc.get_map("p_corr-FWE_method-montecarlo", return_type="array") < 0.05))

    print(f"{model:16s} {np.mean(frac_05):10.3f} {np.mean(frac_001):11.4f} "
          f"{np.mean(any_bonf):9.2f} {np.mean(any_fdr):10.2f} {np.mean(any_mc):8.2f}")
