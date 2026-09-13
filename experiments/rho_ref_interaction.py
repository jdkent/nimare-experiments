"""Do per-study rho and the reference calibration cancel, compound, or act independently?

They use the same sample sizes on entirely different grounds. rho is mechanical: a peak at a
given threshold implies a larger effect when N is small, so small-N studies are discounted
harder. The reference is sociological: researchers who ran few subjects were studying larger
effects, so small-N collections are expected to show more. They push opposite ways.

If they cancel, applying both leaves the answer where neither did, and the pair is doing less
than it appears. Measured by varying N at a fixed truth, which separates them: the truth does
not move, so anything that tracks N is the machinery rather than the effect.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset as gen
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import null_peak_mean_g, reference_magnitude
from nimare.utils import mm2vox

affine = np.array([[4., 0, 0, -40.], [0, 4., 0, -40.], [0, 0, 4., -40.], [0, 0, 0, 1.]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)
centre = tuple(int(c) for c in mm2vox(np.array([[0., 0., 0.]]), affine)[0])
TRUE_G, U = 0.8, 3.2905

print("Truth is fixed at g = 0.8 throughout; only the sample size varies.\n")
print(f"{'N':>5s} {'rho(u,N)':>9s} {'reference':>10s} {'raw g':>7s} {'rho only':>9s} "
      f"{'rho+ref':>8s} {'scale':>7s}")
rows = []
for n in (12, 16, 25, 40, 60, 100):
    raw, rho_only, both, scales = [], [], [], []
    for seed in range(3):
        ss = gen([(0, 0, 0)], effect_sizes=TRUE_G, n_studies=25, sample_size=n, seed=400 + seed,
                 simulate_field=True, noise_extent=40.0, threshold_z=U, smoothness_fwhm=10.0)
        common = dict(fwhm=8.0, mask=mask, null_method="none", threshold="reporting_threshold")
        def value(est):
            return float(abs(est.fit(ss).get_map("g", return_type="image").get_fdata()[centre]))
        raw.append(value(CBES(peak_bias=None, **common)))
        rho_only.append(value(CBES(peak_bias="per-study", peak_bias_scale=1.0, **common)))
        e = CBES(peak_bias="per-study", peak_bias_scale="reference", **common)
        both.append(value(e))
        scales.append(e._peak_bias_scale_)
    # rho for a single-N collection is 1.0 by construction (it is normalised to the median
    # study), so report the quantity it is built from instead.
    rows.append((n, null_peak_mean_g(U, n), reference_magnitude([n] * 25),
                 np.mean(raw), np.mean(rho_only), np.mean(both), np.mean(scales)))
    print(f"{n:5d} {rows[-1][1]:9.3f} {rows[-1][2]:10.3f} {np.mean(raw):7.3f} "
          f"{np.mean(rho_only):9.3f} {np.mean(both):8.3f} {np.mean(scales):7.3f}", flush=True)

arr = np.array([(r[0], r[3], r[4], r[5]) for r in rows], dtype=float)
print("\nHow much does each estimate drift with N, when the truth does not move?")
print("(slope of log estimate on log N; 0 means the machinery has removed the N dependence)")
for i, name in ((1, "raw (no correction)"), (2, "rho only"), (3, "rho + reference")):
    slope = np.polyfit(np.log(arr[:, 0]), np.log(np.clip(arr[:, i], 1e-9, None)), 1)[0]
    spread = arr[:, i].max() / arr[:, i].min()
    print(f"  {name:>22s}  slope {slope:+.3f}   max/min across N {spread:.2f}")
