"""False positive rate of CBES under a global null, on the code paths added since it was
last calibrated: per-study rho, the new threshold rules, and a heterogeneous literature.

Nothing here has a real effect anywhere. A valid estimator should flag ~5% of voxels at an
uncorrected p < .05, and should produce *any* surviving voxel in <=5% of whole simulations
once corrected. The threshold varies from study to study, as it would in a literature search,
which is precisely the regime the new options were built for -- and a regime that could break
calibration, because rho rescales each study by a different factor.
"""
import sys, warnings
import numpy as np
import nibabel as nib
from nimare.correct import FDRCorrector, FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

warnings.simplefilter("ignore")
N_SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 20

affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)

def null_studyset(seed):
    """30 studies, 8 noise foci each, thresholds all over the place, no effect anywhere."""
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
        prevalence=0.0, n_noise_foci=8, noise_extent=36.0, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
    )

CONFIGS = [
    ("pooled-min, no rho",        dict(threshold="pooled-min",          peak_bias=None)),
    ("study-min",       dict(threshold="study-min", peak_bias=None)),
    ("+ per-study rho",           dict(threshold="study-min", peak_bias="per-study")),
    ("oracle threshold + rho",    dict(threshold="reporting_threshold", peak_bias="per-study")),
]

print(f"global null, {N_SIMS} simulations, 30 studies x 8 noise foci, mixed thresholds")
print("(nominal: 0.050 / 0.0010 / 0.05 / 0.05 / 0.05)\n")
print(f"{'configuration':>24s} {'unc p<.05':>10s} {'unc p<.001':>11s} {'FWE-bonf':>9s} "
      f"{'FDR q=.05':>10s} {'FWE-mc':>8s}")

for label, kwargs in CONFIGS:
    frac_05, frac_001, any_bonf, any_fdr, any_mc = [], [], [], [], []
    for seed in range(N_SIMS):
        ss = null_studyset(seed)
        res = CBES(fwhm=12.0, mask=mask, n_iters=200, seed=seed, **kwargs).fit(ss)
        p = res.get_map("p", return_type="array")
        frac_05.append(np.mean(p < 0.05))
        frac_001.append(np.mean(p < 0.001))

        bonf = FWECorrector(method="bonferroni").transform(res)
        any_bonf.append(
            np.any(bonf.get_map("p_corr-FWE_method-bonferroni", return_type="array") < 0.05))
        fdr = FDRCorrector(method="indep", alpha=0.05).transform(res)
        any_fdr.append(
            np.any(fdr.get_map("p_corr-FDR_method-indep", return_type="array") < 0.05))
        mc = FWECorrector(method="montecarlo", n_iters=200).transform(res)
        # The montecarlo corrector emits logp at voxel level, not a plain p map.
        logp = mc.get_map("logp_level-voxel_corr-FWE_method-montecarlo", return_type="array")
        any_mc.append(bool(np.any(logp > -np.log10(0.05))))
        # Emit a running total every simulation, so a timeout cannot lose the whole run.
        print(f"{label:>24s} {np.mean(frac_05):10.3f} {np.mean(frac_001):11.4f} "
              f"{np.mean(any_bonf):9.2f} {np.mean(any_fdr):10.2f} {np.mean(any_mc):8.2f}"
              f"   [{seed + 1}/{N_SIMS}]", flush=True)
