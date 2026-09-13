"""Which of the choices actually earn their place: validity, error rate, cost.

Two regimes. Under a global null nothing is real, so a valid estimator flags ~5% of voxels at
p < .05 and produces any surviving voxel in <=5% of whole simulations once corrected. With a
known effect present, what matters is recovering it. Runtime is measured in both.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.correct import FDRCorrector, FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 25
N_SIGNAL = int(sys.argv[2]) if len(sys.argv) > 2 else 8

affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)
centre = tuple(int(c) for c in mm2vox(np.array([[0.,0.,0.]]), affine)[0])
THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]

CONFIGS = [
    ("pointwise + montecarlo",  dict(censoring="pointwise", null_method="montecarlo")),
    ("rft + montecarlo",        dict(censoring="rft",       null_method="montecarlo")),
    ("pointwise + approximate", dict(censoring="pointwise", null_method="approximate")),
    ("rft + approximate",       dict(censoring="rft",       null_method="approximate")),
]
BASE = dict(fwhm=10.0, mask=mask, threshold="study-min-corrected", peak_bias="per-study",
            smoothness_fwhm=12.0, n_iters=200)

def null_set(seed):
    return create_effect_size_coordinate_studyset(
        [(0,0,0)], effect_sizes=0.0, n_studies=30, sample_size=(15,45), prevalence=0.0,
        n_noise_foci=8, noise_extent=36.0, seed=seed, threshold_z=THRESHOLDS)

def signal_set(seed):
    return create_effect_size_coordinate_studyset(
        [(0,0,0)], effect_sizes=0.6, n_studies=30, sample_size=(15,45), tau=0.1,
        n_noise_foci=2, noise_extent=30.0, seed=1000+seed, threshold_z=THRESHOLDS)

print(f"global null: {N_NULL} sims, 30 studies x 8 noise foci, mixed thresholds")
print("nominal:  unc .050 / .0010 ;  corrected .05\n")
print(f"{'configuration':>24s} {'unc<.05':>8s} {'unc<.001':>9s} {'FDR':>5s} {'FWEmc':>6s} "
      f"{'s/fit':>7s}")
results = {}
for label, kwargs in CONFIGS:
    f05, f001, fdr_hit, mc_hit, secs = [], [], [], [], []
    for seed in range(N_NULL):
        ss = null_set(seed)
        start = time.perf_counter()
        res = CBES(seed=seed, **BASE, **kwargs).fit(ss)
        secs.append(time.perf_counter() - start)
        p = res.get_map("p", return_type="array")
        f05.append(np.mean(p < 0.05)); f001.append(np.mean(p < 0.001))
        fdr = FDRCorrector(method="indep", alpha=0.05).transform(res)
        fdr_hit.append(bool(np.any(
            fdr.get_map("p_corr-FDR_method-indep", return_type="array") < 0.05)))
        mc = FWECorrector(method="montecarlo", n_iters=200).transform(res)
        mc_hit.append(bool(np.any(
            mc.get_map("logp_level-voxel_corr-FWE_method-montecarlo",
                       return_type="array") > -np.log10(0.05))))
    results[label] = dict(f05=np.mean(f05), secs=np.mean(secs))
    print(f"{label:>24s} {np.mean(f05):8.3f} {np.mean(f001):9.4f} {np.mean(fdr_hit):5.2f} "
          f"{np.mean(mc_hit):6.2f} {np.mean(secs):7.1f}", flush=True)

print(f"\nknown effect g = 0.6 at the origin, {N_SIGNAL} sims")
print(f"{'configuration':>24s} {'g at truth':>11s} {'prevalence':>11s} {'s/fit':>7s}")
for label, kwargs in CONFIGS[:2]:
    got, prev, secs = [], [], []
    for seed in range(N_SIGNAL):
        ss = signal_set(seed)
        start = time.perf_counter()
        res = CBES(seed=seed, **{**BASE, "null_method": "parametric"},
                   **{k: v for k, v in kwargs.items() if k != "null_method"}).fit(ss)
        secs.append(time.perf_counter() - start)
        got.append(float(abs(res.get_map("g", return_type="image").get_fdata()[centre])))
        prev.append(float(res.get_map("prevalence", return_type="image").get_fdata()[centre]))
    print(f"{label:>24s} {np.mean(got):11.3f} {np.mean(prev):11.3f} {np.mean(secs):7.1f}",
          flush=True)
