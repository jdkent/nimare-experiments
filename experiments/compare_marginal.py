"""Score the *marginal* CBES map against the images, not the conditional one.

The zero-inflated model separates prevalence pi from magnitude mu: ``g`` is the effect
conditional on a study having a non-zero effect at that voxel, and ``g_marginal = g * pi`` is
the population average. The image truth is an inverse-variance mean over *all* studies -- a
marginal quantity, which includes the studies that have no effect there.

Comparing ``g`` against it therefore overstates the estimate by 1/pi by construction, and every
"CBES runs about 2x high" statement made from the ``g`` map needs checking against this before
it can be believed.
"""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
SDM_COEFF = "/tmp/claude-0/sdm_input/analysis_MyMean/mi/coeff.nii.gz"
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

signal = (truth != 0) & (np.abs(truth) >= np.percentile(np.abs(truth[truth != 0]), 75))
print(f"truth: {int((truth != 0).sum())} voxels, {int(signal.sum())} in the top quartile\n",
      flush=True)


def score(name, est):
    good = np.isfinite(est) & (truth != 0)
    r = float(np.corrcoef(est[good], truth[good])[0, 1])
    s = good & signal
    r_top = float(np.corrcoef(est[s], truth[s])[0, 1])
    mag = float(np.median(np.abs(est[s]) / np.maximum(np.abs(truth[s]), 1e-6)))
    print(f"{name:>26s} {r:8.3f} {r_top:8.3f} {mag:8.2f} {np.abs(est[s]).mean():10.3f}",
          flush=True)


print(f"{'map':>26s} {'r all':>8s} {'r top':>8s} {'mag':>8s} {'|est| top':>10s}", flush=True)

# SDM on the same mask, resampled into masker space so every row is scored identically.
sdm_vec = masker.transform(
    resample_to_img(nib.load(SDM_COEFF), masker.mask_img, interpolation="continuous")
).ravel()
score("SDM-PSI (MLE coeff)", sdm_vec)

if os.path.exists(CACHE):
    coords = pickle.load(open(CACHE, "rb"))
else:
    coords = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
    ).transform(ss)
    pickle.dump(coords, open(CACHE, "wb"))

res = CBES(fwhm=None, mask=masker.mask_img, use_images=False, null_method="none",
           threshold=U, peak_bias="per-study").fit(coords)
g = res.get_map("g", return_type="array")
gm = res.get_map("g_marginal", return_type="array")
pi = res.get_map("prevalence", return_type="array")
score("CBES g (conditional)", g)
score("CBES g_marginal", gm)
print(f"\nprevalence over the top quartile: median {np.median(pi[signal]):.3f}, "
      f"mean {pi[signal].mean():.3f}")
print(f"conditional/marginal magnitude ratio implied by pi: {1/np.median(pi[signal]):.2f}")
