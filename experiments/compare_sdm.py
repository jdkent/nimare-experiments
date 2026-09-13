"""Does SDM-PSI recover the images better than CBES does, from the same coordinates?

Both methods are given exactly the same peaks, derived from the same 21 NIDM pain images at the
same threshold, so nothing turns on curation. The images are the truth: a fixed-effect mean of
the per-study Hedges' g maps is what a meta-analysis of the full images would report, and it is
what each coordinate-only method is trying to recover without them.

Fairness notes, both of which cost CBES an unearned advantage:

* ``use_images=False``. The studyset still carries the images the coordinates came from, and
  CBES would happily use them -- comparing an image meta-analysis against a coordinate one.
* CBES is handed the same reporting threshold SDM gets in its ``t_thr`` column rather than
  inferring it. ``"study-min"`` is also run, as the realistic case where nobody tells you.

Pattern and magnitude are scored separately because they are separate claims: a method can put
the effect in the right place and be wrong about its size by a factor of two, which is exactly
the failure mode this estimator has been fighting. Magnitude uses a paired median ratio over the
reference's top quartile, not a ratio of means -- the reference's whole-brain mean is near zero,
and dividing by it manufactures ratios that describe the denominator (notes section 17).
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
SDM_DIR = "/tmp/claude-0/sdm_input"
SDM_COEFF = f"{SDM_DIR}/analysis_MyMean/mi/coeff.nii.gz"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

# Fixed-effect mean of the per-study g maps, in masker space.
total = weight = None
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    g = np.where(ok, g, 0.0)
    total = w * g if total is None else total + w * g
    weight = w if weight is None else weight + w
truth_vec = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
truth_img = masker.inverse_transform(truth_vec)

sdm_img = nib.load(SDM_COEFF)
sdm = sdm_img.get_fdata()
truth = resample_to_img(truth_img, sdm_img, interpolation="continuous").get_fdata()

sdm_mask = nib.load(f"{SDM_DIR}/sdm_mask.nii.gz")
in_sdm = resample_to_img(sdm_mask, sdm_img, interpolation="nearest").get_fdata() > 0
valid = in_sdm & np.isfinite(truth) & (truth != 0) & np.isfinite(sdm)
print(f"scoring over {valid.sum()} voxels of the SDM gray-matter mask\n", flush=True)


def score(name, estimate):
    est, ref = estimate[valid], truth[valid]
    good = np.isfinite(est) & np.isfinite(ref)
    est, ref = est[good], ref[good]
    r = float(np.corrcoef(est, ref)[0, 1])
    top = np.abs(ref) >= np.percentile(np.abs(ref), 75)
    ratio = float(np.median(np.abs(est[top]) / np.maximum(np.abs(ref[top]), 1e-6)))
    # Correlation restricted to where the reference says there is signal, so a good score
    # cannot come from agreeing about empty brain.
    r_top = float(np.corrcoef(est[top], ref[top])[0, 1])
    print(f"{name:>24s} {r:8.3f} {r_top:8.3f} {ratio:9.2f} {np.abs(est[top]).mean():9.3f} "
          f"{np.abs(ref[top]).mean():9.3f}", flush=True)


print(f"{'method':>24s} {'r all':>8s} {'r top':>8s} {'mag':>9s} {'|est|':>9s} {'|ref|':>9s}",
      flush=True)
score("SDM-PSI (MLE coeff)", sdm)

coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss)

for label, kwargs in [
    ("CBES, known u", dict(threshold=U, peak_bias=None)),
    ("CBES, known u + rho", dict(threshold=U, peak_bias="per-study")),
    ("CBES, study-min + rho", dict(threshold="study-min", peak_bias="per-study")),
]:
    res = CBES(fwhm=None, mask=masker.mask_img, use_images=False, null_method="none",
               **kwargs).fit(coords)
    arr = resample_to_img(res.get_map("g"), sdm_img, interpolation="continuous").get_fdata()
    score(label, arr)
