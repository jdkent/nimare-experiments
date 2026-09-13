"""Does image-calibrated peak_bias recover a KNOWN effect field from coordinates alone?"""
import json, os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/cmp")
import numpy as np, nibabel as nib
from scipy import stats
from field_sim import SHAPE, true_field, study_image, report_peaks
from nimare.studyset import Studyset
from nimare.meta.cbma import CBES
from nimare.transforms import d_to_g, t_to_z

OUT = "/tmp/claude-0/cmp/simimg"
os.makedirs(OUT, exist_ok=True)
AFFINE = np.diag([3.0, 3.0, 3.0, 1.0]); AFFINE[:3, 3] = -43.5
THRESH_Z = 3.2905

def build(n_studies=20, peak_g=0.6, smooth=2.0, seed=0):
    rng = np.random.default_rng(seed)
    truth = true_field(peak_g=peak_g)
    nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE), f"{OUT}/mask.nii.gz")
    studies = []
    for k in range(n_studies):
        n = int(rng.integers(15, 40))
        img = study_image(truth, n, smooth, rng)
        nib.save(nib.Nifti1Image(img.astype(np.float32), AFFINE), f"{OUT}/s{k}_g.nii.gz")
        nib.save(nib.Nifti1Image(np.full(SHAPE, 1.0 / n, np.float32), AFFINE),
                 f"{OUT}/s{k}_gvar.nii.gz")
        hits, values = report_peaks(img, n, THRESH_Z)
        points = []
        for ijk, g in zip(hits, values):
            xyz = nib.affines.apply_affine(AFFINE, ijk)
            bias = 1.0 - 3.0 / (4.0 * (n - 1) - 1)
            z = float(t_to_z(np.array([g / bias * np.sqrt(n)]), n - 1)[0])
            points.append({"space": "MNI", "coordinates": [float(c) for c in xyz],
                           "values": [{"kind": "Z", "value": z}]})
        studies.append({
            "id": f"s{k}", "name": f"s{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"s{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                          "points": points,
                          "images": [
                              {"url": f"{OUT}/s{k}_g.nii.gz", "filename": f"s{k}_g.nii.gz",
                               "space": "MNI", "value_type": "g"},
                              {"url": f"{OUT}/s{k}_gvar.nii.gz", "filename": f"s{k}_gvar.nii.gz",
                               "space": "MNI", "value_type": "g_var"}]}]})
    return truth, Studyset({"id": "sim", "name": "sim", "studies": studies},
                           target=None, mask=f"{OUT}/mask.nii.gz")

truth, ss = build()
masker = ss.masker
truth_vec = masker.transform(nib.Nifti1Image(truth.astype(np.float32), AFFINE)).ravel()
ids = np.array(list(ss.ids))
print(f"{len(ids)} studies, {len(ss.coordinates)} reported peaks, "
      f"true peak g = {truth_vec.max():.3f}, true mean g = {truth_vec.mean():.3f}")

def fit(studyset, use_images, peak_bias=None):
    est = CBES(fwhm=9.0, null_method="parametric", use_images=use_images, peak_bias=peak_bias)
    return est.fit(studyset).get_map("g", return_type="array").ravel()

def restrict(s, keep): return s.slice(analyses=[a for a in s.ids if a in set(keep)])

# Calibrate rho on half the studies, from their images vs their coordinates.
rng = np.random.default_rng(3)
perm = rng.permutation(len(ids)); calib, test = ids[perm[:10]], ids[perm[10:]]
cal = restrict(ss, calib)
im, co = fit(cal, True), fit(cal, False)
ok = (im != 0) & (co != 0)
rho = float(im[ok].mean() / co[ok].mean())
print(f"calibrated rho (from 10 calibration studies' images) = {rho:.3f}\n")

hot = truth_vec > 0.1
print(f"{'estimate':34s} {'mean g (true blob)':>19s} {'r with truth':>13s}")
print(f"{'TRUE FIELD':34s} {truth_vec[hot].mean():19.3f} {'':>13s}")
test_set = restrict(ss, test)
for label, kw in (("held-out coords, uncorrected", dict(use_images=False)),
                  ("held-out coords, peak_bias=rho", dict(use_images=False, peak_bias=rho)),
                  ("held-out images (upper bound)", dict(use_images=True))):
    g = fit(test_set, **kw)
    cov = g != 0
    print(f"{label:34s} {g[hot & cov].mean():19.3f} "
          f"{stats.pearsonr(g[cov], truth_vec[cov])[0]:13.3f}")
