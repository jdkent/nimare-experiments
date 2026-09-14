"""E. Does the uniform relocation null survive foci that clump for non-effect reasons?

The single most dangerous assumption CBES inherits. correct_fwe_montecarlo relocates every
focus to a uniformly drawn in-mask voxel. Real foci are not uniform: they concentrate in gray
matter and in a handful of much-studied regions, for reasons that have nothing to do with the
effect under study. If the data clump and the null does not, the null is too diffuse, the
observed clumping looks extraordinary, and the false positive rate goes up.

Here there is no true effect anywhere at any setting. The only thing that changes between rows
is how tightly the noise foci are scattered. A method whose error rate depends on that is
reporting the clumping, not the effect.
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), AFFINE)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12

print(f"no true effect anywhere; {N} simulations per row\n", flush=True)
print(f"{'null foci spread':>18s} {'unc p<.05':>10s} {'unc p<.001':>11s} {'FWE .05':>9s} "
      f"{'verdict':>9s}", flush=True)
rows = []
for extent, label in ((36.0, "diffuse 36mm"), (18.0, "clumped 18mm"), (10.0, "tight 10mm")):
    u5, u1, fwe = [], [], []
    for seed in range(N):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
            prevalence=0.0, n_noise_foci=8, noise_extent=extent, seed=seed,
            threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])
        est = CBES(fwhm=12.0, mask=MASK, null_method="montecarlo", n_iters=150,
                   threshold="study-min", peak_bias="per-study", seed=seed)
        res = est.fit(ss)
        p = res.get_map("p", return_type="array")
        u5.append(float(np.mean(p < 0.05)))
        u1.append(float(np.mean(p < 0.001)))
        maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
        fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
    rate = float(np.mean(fwe))
    ok = rate <= 0.05 + 3 * np.sqrt(0.05 * 0.95 / N)
    rows.append((label, ok, rate))
    print(f"{label:>18s} {np.mean(u5):10.4f} {np.mean(u1):11.4f} {rate:9.2f} "
          f"{'PASS' if ok else 'FAIL':>9s}", flush=True)

print("\n(nominal 0.050 / 0.0010 / 0.05)")
worst = max(r for _, _, r in rows)
best = min(r for _, _, r in rows)
print(f"FWE rate varies {best:.2f} to {worst:.2f} across clumping conditions")
print("E VERDICT:", "PASS" if all(ok for _, ok, _ in rows) else
      "FAIL -- the error rate depends on how the null foci clump")
