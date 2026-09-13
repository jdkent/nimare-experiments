"""Which censoring and null choices earn their place: validity and cost.

Scoped to the measurement that is well powered here. The uncorrected false positive rate is a
mean over ~9000 voxels per simulation, so 15 simulations pin it tightly; the corrected rate is
one Bernoulli per simulation and needs ~100, which is a separate exercise. Prints per
simulation so a timeout cannot lose the run.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 15
N_SIG = int(sys.argv[2]) if len(sys.argv) > 2 else 6
affine = np.array([[4., 0, 0, -40.], [0, 4., 0, -40.], [0, 0, 4., -40.], [0, 0, 0, 1.]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)
centre = tuple(int(c) for c in mm2vox(np.array([[0., 0., 0.]]), affine)[0])
THR = [2.3263, 3.0902, 3.2905, 4.2649]
BASE = dict(fwhm=10.0, mask=mask, threshold="study-min", peak_bias="per-study",
            smoothness_fwhm=12.0, n_iters=100)
CONFIGS = [
    ("pointwise + montecarlo", dict(censoring="pointwise", null_method="montecarlo")),
    ("rft + montecarlo", dict(censoring="rft", null_method="montecarlo")),
    ("pointwise + approximate", dict(censoring="pointwise", null_method="approximate")),
    ("rft + approximate", dict(censoring="rft", null_method="approximate")),
]

print(f"global null, {N_NULL} sims; nominal .050 / .0010\n", flush=True)
print(f"{'configuration':>24s} {'unc<.05':>8s} {'unc<.001':>9s} {'s/fit':>7s}", flush=True)
for label, kw in CONFIGS:
    f05, f001, secs = [], [], []
    for seed in range(N_NULL):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(15, 45), prevalence=0.0,
            n_noise_foci=8, noise_extent=36.0, seed=seed, threshold_z=THR)
        t0 = time.perf_counter()
        res = CBES(seed=seed, **BASE, **kw).fit(ss)
        secs.append(time.perf_counter() - t0)
        p = res.get_map("p", return_type="array")
        f05.append(np.mean(p < 0.05))
        f001.append(np.mean(p < 0.001))
        print(f"{label:>24s} {np.mean(f05):8.3f} {np.mean(f001):9.4f} {np.mean(secs):7.1f}"
              f"   [{seed + 1}/{N_NULL}]", flush=True)

print(f"\nknown g = 0.6, field-simulated so peaks are selected, {N_SIG} sims", flush=True)
print(f"{'configuration':>24s} {'g at truth':>11s} {'prevalence':>11s} {'s/fit':>7s}", flush=True)
for label, kw in CONFIGS[:2]:
    got, prev, secs = [], [], []
    for seed in range(N_SIG):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.6, n_studies=25, sample_size=25, seed=100 + seed,
            simulate_field=True, noise_extent=40.0, threshold_z=3.2905, smoothness_fwhm=12.0)
        t0 = time.perf_counter()
        res = CBES(seed=seed, **{**BASE, "null_method": "none"},
                   **{k: v for k, v in kw.items() if k != "null_method"}).fit(ss)
        secs.append(time.perf_counter() - t0)
        got.append(float(abs(res.get_map("g", return_type="image").get_fdata()[centre])))
        prev.append(float(res.get_map("prevalence", return_type="image").get_fdata()[centre]))
        print(f"{label:>24s} {np.mean(got):11.3f} {np.mean(prev):11.3f} {np.mean(secs):7.1f}"
              f"   [{seed + 1}/{N_SIG}]", flush=True)
