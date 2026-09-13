"""Does flattening the kernel's distance profile recover the images better?

weighting4 showed the defect is the 1/w(d) variance inflation: a coordinate's actual error is
nearly flat across the kernel's support (0.063 -> 0.082 over 13 mm) while the model's claimed
variance grows 34x. A scalar multiplier cannot fix a shape, so the question is what shape is
right -- and the images are a criterion that can answer it.

``w -> w**alpha`` spans the hypotheses with one parameter, leaving the *support* untouched so
coverage, the censoring term and n_studies are all unchanged:

    alpha = 1    what the estimator does now: precision w/var
    alpha = 0.5  square-root discounting
    alpha = 0    the kernel decides only *which* voxels a study speaks to, not how precisely

Scored against the same image-based truth the SDM comparison uses, on pattern and magnitude
separately, coordinate-only (so the images are purely the criterion and never an input).
"""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as fx
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
CACHE = "/tmp/claude-0/cmp/pain_coords_ss.pkl"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

total = weight = None
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)

if os.path.exists(CACHE):
    coords = pickle.load(open(CACHE, "rb"))
else:
    coords = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
    ).transform(ss)
    pickle.dump(coords, open(CACHE, "wb"))

signal = np.abs(truth) >= np.percentile(np.abs(truth[truth != 0]), 75)
print(f"truth: {int((truth != 0).sum())} voxels, {int(signal.sum())} in the top quartile\n",
      flush=True)
print(f"{'alpha':>7s} {'r all':>8s} {'r top':>8s} {'mag':>8s} {'|est| top':>10s}", flush=True)

original = CBES._kernel_support


def score(label, est):
    good = np.isfinite(est) & (truth != 0)
    r = float(np.corrcoef(est[good], truth[good])[0, 1])
    s = good & signal
    r_top = float(np.corrcoef(est[s], truth[s])[0, 1])
    mag = float(np.median(np.abs(est[s]) / np.maximum(np.abs(truth[s]), 1e-6)))
    print(f"{label:>7s} {r:8.3f} {r_top:8.3f} {mag:8.2f} {np.abs(est[s]).mean():10.3f}",
          flush=True)


for alpha in (1.0, 0.75, 0.5, 0.25, 0.0):
    def patched(self, sample_size=None, _a=alpha):
        offsets, values = original(self, sample_size=sample_size)
        return offsets, values**_a
    CBES._kernel_support = patched
    try:
        res = CBES(fwhm=None, mask=masker.mask_img, use_images=False, null_method="none",
                   threshold=U, peak_bias="per-study").fit(coords)
        score(f"{alpha:.2f}", res.get_map("g", return_type="array"))
    finally:
        CBES._kernel_support = original
