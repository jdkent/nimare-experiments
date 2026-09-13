"""Does rate-based calibration recover the true effect size, with no images?

The simulator generates a known g, thresholds each study at its own level, and reports only
coordinates. A calibrated estimator should return that g; an uncalibrated one returns it times
an arbitrary constant. The comparison is against the truth, not against another estimator.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

affine = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)
centre = tuple(int(c) for c in mm2vox(np.array([[0.0, 0.0, 0.0]]), affine)[0])

print(f"{'true g':>7s} {'studies':>8s} {'uncalibrated':>13s} {'rate-calibrated':>16s} "
      f"{'scale':>7s} {'error':>8s}")
for true_g in (0.3, 0.5, 0.8):
    for n_studies in (30, 60):
        raw, calibrated, scales = [], [], []
        for seed in range(5):
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=true_g, n_studies=n_studies,
                sample_size=(12, 60), tau=0.1, seed=seed, spatial_sd=4.0,
                n_noise_foci=2, noise_extent=30.0,
                threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
            )
            common = dict(fwhm=8.0, mask=mask, null_method="parametric",
                          threshold="reporting_threshold", peak_bias="per-study")
            a = CBES(peak_bias_scale=1.0, **common).fit(ss)
            b_est = CBES(peak_bias_scale="rates", **common)
            b = b_est.fit(ss)
            def grab(result):
                volume = result.get_map("g", return_type="image").get_fdata()
                return float(np.abs(volume[centre]))
            raw.append(grab(a))
            calibrated.append(grab(b))
            scales.append(b_est._peak_bias_scale_)
        raw, calibrated = np.mean(raw), np.mean(calibrated)
        print(f"{true_g:7.2f} {n_studies:8d} {raw:13.3f} {calibrated:16.3f} "
              f"{np.mean(scales):7.3f} {calibrated - true_g:+8.3f}", flush=True)
