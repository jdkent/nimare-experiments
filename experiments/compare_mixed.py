"""Head to head on the same mixture: 5 images + 16 coordinate studies, scored against all 21.

Every row is scored identically -- masker space, the same inverse-variance truth built from all
21 images, the same top-quartile signal mask -- so the only thing that differs between rows is
the procedure. The coordinate-only rows are carried over for reference, which is what makes the
mixture's contribution readable: how much does each method gain from five real images?
"""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer

N_IMAGE = 5
SDM_COORD = "/tmp/claude-0/sdm_input/analysis_MyMean/mi/coeff.nii.gz"
SDM_MIXED = "/tmp/claude-0/sdm_mixed/analysis_MyMean/mi/coeff.nii.gz"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

per = {}
total = weight = None
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    per[sid] = (np.where(ok, g, 0.0), w)
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
signal = (truth != 0) & (np.abs(truth) >= np.percentile(np.abs(truth[truth != 0]), 75))
image_ids = sorted(per)[:N_IMAGE]

print(f"truth from all 21 images; {int(signal.sum())} voxels in the top quartile")
print(f"mixture: images for {', '.join(i.split('.')[0] for i in image_ids)}\n", flush=True)
print(f"{'procedure':>30s} {'r all':>8s} {'r top':>8s} {'mag':>8s}", flush=True)


def score(name, est):
    good = np.isfinite(est) & (truth != 0)
    s = good & signal
    r = float(np.corrcoef(est[good], truth[good])[0, 1])
    r_top = float(np.corrcoef(est[s], truth[s])[0, 1])
    mag = float(np.median(np.abs(est[s]) / np.maximum(np.abs(truth[s]), 1e-6)))
    print(f"{name:>30s} {r:8.3f} {r_top:8.3f} {mag:8.2f}", flush=True)


def from_nifti(path):
    return masker.transform(
        resample_to_img(nib.load(path), masker.mask_img, interpolation="continuous")
    ).ravel()


tot = wgt = None
for sid in image_ids:
    g, w = per[sid]
    tot = w * g if tot is None else tot + w * g
    wgt = w if wgt is None else wgt + w
score(f"{N_IMAGE} images alone", np.divide(tot, wgt, out=np.zeros_like(tot), where=wgt > 0))

if os.path.exists(SDM_COORD):
    score("SDM, coordinates only", from_nifti(SDM_COORD))
if os.path.exists(SDM_MIXED):
    score("SDM, 5 images + coords", from_nifti(SDM_MIXED))

coords = pickle.load(open("/tmp/claude-0/cmp/pain_coords_ss.pkl", "rb"))
original_load = CBES._load_image_studies


def only_five(self, dataset):
    loaded = original_load(self, dataset)
    return {k: v for k, v in (loaded or {}).items() if str(k) in image_ids}


res = CBES(fwhm=None, mask=masker.mask_img, use_images=False, null_method="none",
           threshold=3.2905, peak_bias="per-study").fit(coords)
score("CBES, coordinates only", res.get_map("g_marginal", return_type="array"))

CBES._load_image_studies = only_five
try:
    res = CBES(fwhm=None, mask=masker.mask_img, use_images=True, null_method="none",
               threshold=3.2905, peak_bias="per-study").fit(coords)
    assert len(res.estimator._image_studies_) == N_IMAGE, "image filter did not apply"
    score("CBES, 5 images + coords", res.get_map("g_marginal", return_type="array"))
finally:
    CBES._load_image_studies = original_load
