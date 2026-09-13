"""How many images does it take to calibrate the scale, against a known truth?

Earlier figures for this came from a fixture that read a value at the ground-truth location, so
it had no peak-height inflation and nothing to calibrate. The field simulator selects local
maxima, so the inflation is real and the question is answerable.

Also fixes the rho test: rho is normalised to the median study, so it does nothing when every
study shares a sample size. Heterogeneous N within a collection is what it exists for.
"""
import warnings; warnings.simplefilter("ignore")
import os, tempfile
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset as gen
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

affine = np.array([[4., 0, 0, -40.], [0, 4., 0, -40.], [0, 0, 4., -40.], [0, 0, 0, 1.]])
shape = (21, 21, 21)
mask_img = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
centre = tuple(int(c) for c in mm2vox(np.array([[0., 0., 0.]]), affine)[0])
TRUE_G, U = 0.8, 3.2905
grid = np.stack(np.indices(shape), -1) * 4.0 - 40.0
blob_sigma = 10.0 / (2 * np.sqrt(2 * np.log(2)))
truth_field = TRUE_G * np.exp(-((grid ** 2).sum(-1)) / (2 * blob_sigma ** 2))

print("Part 1: rho with heterogeneous sample sizes, which is what it corrects.\n")
print(f"{'sample sizes':>22s} {'raw g':>7s} {'rho only':>9s} {'truth':>7s}")
for label, spec in (("uniform N = 25", 25), ("mixed N = 12 to 60", (12, 60)),
                    ("mixed N = 10 to 100", (10, 100))):
    raw, rho = [], []
    for seed in range(3):
        ss = gen([(0, 0, 0)], effect_sizes=TRUE_G, n_studies=25, sample_size=spec,
                 seed=500 + seed, simulate_field=True, noise_extent=40.0,
                 threshold_z=[2.3263, 3.0902, 3.2905, 4.2649], smoothness_fwhm=10.0)
        common = dict(fwhm=8.0, mask=mask_img, null_method="none",
                      threshold="reporting_threshold")
        f = lambda e: float(abs(e.fit(ss).get_map("g", return_type="image").get_fdata()[centre]))
        raw.append(f(CBES(peak_bias=None, **common)))
        rho.append(f(CBES(peak_bias="per-study", peak_bias_scale=1.0, **common)))
    print(f"{label:>22s} {np.mean(raw):7.3f} {np.mean(rho):9.3f} {TRUE_G:7.3f}", flush=True)

print("\nPart 2: how many images does peak_bias_scale='auto' need?\n")
directory = tempfile.mkdtemp()
print(f"{'images':>7s} {'coord studies':>14s} {'calibrated g':>13s} {'scale':>7s} {'err x':>7s}")
for n_images in (0, 2, 5, 10, 15):
    got, scales = [], []
    for seed in range(3):
        rng = np.random.default_rng(600 + seed)
        ss = gen([(0, 0, 0)], effect_sizes=TRUE_G, n_studies=25, sample_size=(15, 45),
                 seed=600 + seed, simulate_field=True, noise_extent=40.0,
                 threshold_z=U, smoothness_fwhm=10.0)
        # Add image studies drawn from the same truth, with honest sampling noise.
        studies = list(ss.to_dict()["studies"])
        for k in range(n_images):
            n = int(rng.integers(15, 45))
            obs = truth_field + rng.normal(0, 1 / np.sqrt(n), shape)
            gp = f"{directory}/s{seed}_{k}_g.nii.gz"
            vp = f"{directory}/s{seed}_{k}_v.nii.gz"
            nib.save(nib.Nifti1Image(obs.astype(np.float32), affine), gp)
            nib.save(nib.Nifti1Image(np.full(shape, 1.0 / n, np.float32), affine), vp)
            meta = {"sample_sizes": [n], "reporting_threshold": U}
            studies.append({"id": f"img{k}", "name": f"img{k}", "metadata": meta, "analyses": [{
                "id": f"img{k}-1", "name": "1", "metadata": meta,
                "points": [{"space": "MNI", "coordinates": [0., 0., 0.],
                            "values": [{"kind": "Z", "value": float(U + 0.5)}]}],
                "images": [{"url": gp, "filename": "g", "space": "MNI", "value_type": "g"},
                           {"url": vp, "filename": "v", "space": "MNI", "value_type": "g_var"}]}]})
        mixed = Studyset({"id": "m", "name": "m", "studies": studies}, target=None,
                         mask=mask_img)
        est = CBES(fwhm=8.0, mask=mask_img, null_method="none", peak_bias="per-study",
                   threshold="reporting_threshold",
                   peak_bias_scale=("auto" if n_images else 1.0), use_images=bool(n_images))
        got.append(float(abs(est.fit(mixed).get_map("g", return_type="image").get_fdata()[centre])))
        scales.append(est._peak_bias_scale_)
    g = np.mean(got)
    print(f"{n_images:7d} {25:14d} {g:13.3f} {np.mean(scales):7.3f} "
          f"{max(g, TRUE_G) / min(g, TRUE_G):7.2f}", flush=True)
