"""Does ``peak_bias='per-study'`` flatten the threshold-driven level, judged on held-out data?

The estimator documents ``peak_bias`` as the remedy for exactly the failure the cross-dataset
run measured: a level set by ``(threshold, sample size)`` rather than by the effect. It divides
each study's values by the mean height its own threshold and sample size would produce under
pure noise, which is the right shape of correction if that is what is wrong. It has never been
scored against a reference the peaks did not produce.

Same split-half protocol as ``cross_dataset_floor``: half the studies become coordinates, the
other half is pooled by inverse variance for the truth. Three settings are compared on the same
splits -- no correction, the per-study correction, and the marginal ``pi * mu`` -- so the
question "does dividing out the noise height fix the flatness" gets a direct answer.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.transforms import ImageTransformer
from load_pain import load_pain
from cross_dataset_floor import (MAX_PEAKS, MIN_HALF, mask_img, masker, affine,
                                 peaks_of, pooled_truth)
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

CUT = 3.2905
N_SPLITS = 6
rng = np.random.default_rng(0)


def build(maps, sizes):
    studies = []
    for k, (z, n) in enumerate(zip(maps, sizes)):
        found = peaks_of(z, CUT, MAX_PEAKS)
        if len(found) < 2:
            continue
        meta = {"sample_sizes": [int(n)]}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": v}]} for ijk, v in found]}]})
    return studies


ss = ImageTransformer(target="z").transform(load_pain())
maps, sizes = [], []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    path = getattr(row, "z", None)
    if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
        continue
    img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                          force_resample=True, copy_header=True)
    maps.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
    sizes.append(float(n))
sizes = np.asarray(sizes, dtype=float)
n = len(maps)
print(f"NIDM pain: {n} studies with usable z maps, mean N {sizes.mean():.0f}\n", flush=True)

SETTINGS = (("peak_bias=None", None), ("peak_bias='per-study'", "per-study"))
out = {name: {"ratio": [], "r": [], "strata": []} for name, _ in SETTINGS}
out["pi*mu (peak_bias=None)"] = {"ratio": [], "r": [], "strata": []}

for _ in range(N_SPLITS):
    order = rng.permutation(n)
    lo, hi = order[: n // 2], order[n // 2:]
    truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
    studies = build([maps[i] for i in lo], sizes[lo])
    if len(studies) < MIN_HALF:
        continue
    bands = [(truth >= np.percentile(truth, a)) &
             (truth < np.percentile(truth, b) if b < 100 else np.ones_like(truth, bool))
             for a, b in ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100))]
    for name, setting in SETTINGS:
        est = CBES(fwhm=10.0, mask=masker, peak_bias=setting, null_method="none",
                   threshold="study-min")
        res = est.fit(Studyset({"id": "p", "name": "p", "studies": studies},
                               target=None, mask=mask_img))
        g = np.abs(res.get_map("g", return_type="array").ravel())
        use = (res.get_map("n_studies", return_type="array").ravel() > 0) & np.isfinite(g) & (g > 0)
        if use.sum() < 100:
            continue
        for label, value in ((name, g),) + ((("pi*mu (peak_bias=None)",
             g * res.get_map("prevalence", return_type="array").ravel()),) if setting is None else ()):
            out[label]["ratio"].append(value[use].mean() / max(truth[use].mean(), 1e-9))
            out[label]["r"].append(stats.pearsonr(value[use], truth[use])[0])
            out[label]["strata"].append([
                value[use & b].mean() / max(truth[use & b].mean(), 1e-9)
                if (use & b).sum() >= 30 else np.nan for b in bands])

print(f"{'setting':>24} {'overall ratio':>14} {'r':>7}   "
      + " ".join(f"{s:>8}" for s in ("0-50%", "50-75%", "75-90%", "90-99%", "99-100%")))
for label, rec in out.items():
    if not rec["ratio"]:
        print(f"{label:>24}   no usable split")
        continue
    strata = np.nanmean(np.array(rec["strata"], dtype=float), axis=0)
    print(f"{label:>24} {np.mean(rec['ratio']):14.2f} {np.mean(rec['r']):+7.3f}   "
          + " ".join(f"{s:8.2f}" for s in strata))
print("\nA ratio flat across the strata means the map is reading the threshold; one falling"
      "\ntoward 1.00 at the right means it is reading the effect.")
