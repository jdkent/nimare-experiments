"""The mixed case: do under-weighted coordinates actually get swamped by images?

weighting4 measured coordinates at 0.035x the weight their errors justify. The alpha sweep then
showed that reshaping the distance profile changes neither recovery nor the null, because the
badly-weighted contributions are the far ones that carry little weight anyway. So the
mis-weighting can only matter where coordinates compete with images -- which is the open issue.

The test gives CBES the real images for a handful of studies and coordinates only for the rest,
and scores against the truth built from *all* 21 images. If coordinates are swamped, the fit sits
at the few-image baseline however many coordinate studies are added, and up-weighting them should
move it toward the full-image truth. If it does not move, they are not being swamped and the
measured mis-weighting has no practical consequence.

``scale`` multiplies every coordinate study's kernel weight; images always enter at 1. Because
tau^2 and mu are invariant to a uniform rescaling of all weights, this changes nothing except the
balance between the two kinds of study.
"""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
N_IMAGE = 5
CACHE = "/tmp/claude-0/cmp/pain_coords_ss.pkl"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

per_study = {}
total = weight = None
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    per_study[sid] = (np.where(ok, g, 0.0), w)
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
signal = (truth != 0) & (np.abs(truth) >= np.percentile(np.abs(truth[truth != 0]), 75))

image_ids = sorted(per_study)[:N_IMAGE]
print(f"{len(per_study)} studies; images given for {N_IMAGE}, coordinates for the rest")
print(f"truth from all 21 images; {int(signal.sum())} voxels in the top quartile\n", flush=True)


def score(name, est):
    good = np.isfinite(est) & (truth != 0)
    r = float(np.corrcoef(est[good], truth[good])[0, 1])
    s = good & signal
    r_top = float(np.corrcoef(est[s], truth[s])[0, 1])
    mag = float(np.median(np.abs(est[s]) / np.maximum(np.abs(truth[s]), 1e-6)))
    print(f"{name:>28s} {r:8.3f} {r_top:8.3f} {mag:8.2f}", flush=True)


print(f"{'configuration':>28s} {'r all':>8s} {'r top':>8s} {'mag':>8s}", flush=True)

# Baseline: the few images on their own, with no coordinates at all.
tot = wgt = None
for sid in image_ids:
    g, w = per_study[sid]
    tot = w * g if tot is None else tot + w * g
    wgt = w if wgt is None else wgt + w
score(f"{N_IMAGE} images alone", np.divide(tot, wgt, out=np.zeros_like(tot), where=wgt > 0))

if os.path.exists(CACHE):
    coords = pickle.load(open(CACHE, "rb"))
else:
    coords = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
    ).transform(ss)
    pickle.dump(coords, open(CACHE, "wb"))

original_load = CBES._load_image_studies
original_kernel = CBES._kernel_support


def only_selected(self, dataset):
    loaded = original_load(self, dataset)
    return {k: v for k, v in (loaded or {}).items() if str(k) in image_ids}


for scale in (1.0, 5.0, 29.0):
    def scaled(self, sample_size=None, _s=scale):
        offsets, values = original_kernel(self, sample_size=sample_size)
        return offsets, values * _s
    CBES._load_image_studies = only_selected
    CBES._kernel_support = scaled
    try:
        res = CBES(fwhm=None, mask=masker.mask_img, use_images=True, null_method="none",
                   threshold=U, peak_bias="per-study", peak_bias_scale=1.0).fit(coords)
        score(f"+ coords, weight x{scale:g}", res.get_map("g_marginal", return_type="array"))
    finally:
        CBES._load_image_studies = original_load
        CBES._kernel_support = original_kernel
