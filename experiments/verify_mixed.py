"""Confirm the mixed test really withheld images: count them, and check the all-images ceiling."""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
total = weight = None
per = {}
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    per[sid] = 1
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
signal = (truth != 0) & (np.abs(truth) >= np.percentile(np.abs(truth[truth != 0]), 75))

coords = pickle.load(open("/tmp/claude-0/cmp/pain_coords_ss.pkl", "rb"))
image_ids = sorted(per)[:5]
original_load = CBES._load_image_studies

def only_selected(self, dataset):
    loaded = original_load(self, dataset)
    kept = {k: v for k, v in (loaded or {}).items() if str(k) in image_ids}
    print(f"  [filter] loaded {len(loaded or {})} -> kept {len(kept)}", flush=True)
    return kept

def report(name, est):
    good = np.isfinite(est) & (truth != 0)
    s = good & signal
    print(f"{name:>22s}  r_all {np.corrcoef(est[good], truth[good])[0,1]:.3f}  "
          f"r_top {np.corrcoef(est[s], truth[s])[0,1]:.3f}  "
          f"mag {np.median(np.abs(est[s])/np.maximum(np.abs(truth[s]),1e-6)):.2f}", flush=True)

print("A) all 21 images -- the ceiling, should be ~1.0 if the harness is sound")
res = CBES(fwhm=None, mask=masker.mask_img, use_images=True, null_method="none",
           threshold=3.2905, peak_bias="per-study").fit(coords)
report("all images", res.get_map("g_marginal", return_type="array"))
print(f"   _image_studies_ = {len(res.estimator._image_studies_)}")

print("\nB) only 5 images, rest coordinates")
CBES._load_image_studies = only_selected
try:
    res = CBES(fwhm=None, mask=masker.mask_img, use_images=True, null_method="none",
               threshold=3.2905, peak_bias="per-study").fit(coords)
    report("5 images + coords", res.get_map("g_marginal", return_type="array"))
    print(f"   _image_studies_ = {len(res.estimator._image_studies_)}")
finally:
    CBES._load_image_studies = original_load
