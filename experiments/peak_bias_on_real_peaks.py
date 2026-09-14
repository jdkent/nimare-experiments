"""Finding 7 on real peaks, with an image-based truth, independent of any simulator.

Every finding-7 number so far comes from the field simulator. The pain collection is the
independent check: its peaks are genuine local maxima of genuine statistic images, extracted
at a known threshold, and the images themselves give a truth that owes nothing to the
simulator's assumptions about smoothness, blob shape or effect distribution.

Truth is the inverse-variance pooled g across all images. Coordinates are extracted from the
same images and capped at a realistic number of peaks per study, strongest-first. Scored where
the truth is largest, since that is where a magnitude claim is made, and on the voxels CBES
covers.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates

# Image truth: inverse-variance pooled g over every study that supplies both maps.
g_maps, w_maps = [], []
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g = masker.transform(str(row.g)).ravel()
    v = masker.transform(str(row.g_var)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_maps.append(np.where(ok, g, 0.0))
    w_maps.append(np.where(ok, 1.0 / np.maximum(v, 1e-12), 0.0))
G, W = np.array(g_maps), np.array(w_maps)
truth = np.sum(G * W, axis=0) / np.maximum(np.sum(W, axis=0), 1e-12)
print(f"image truth from {len(g_maps)} studies; |truth| median "
      f"{np.median(np.abs(truth)):.3f}, top-quartile mean "
      f"{np.abs(truth)[np.abs(truth) >= np.percentile(np.abs(truth), 75)].mean():.3f}\n")

ids = sorted(set(str(i) for i in coords["id"]))


def build(max_peaks):
    studies = []
    for sid in ids:
        sub = coords[coords["id"].astype(str) == sid].copy()
        sub["_abs"] = np.abs(sub["z_stat"].astype(float))
        sub = sub.sort_values("_abs", ascending=False).head(max_peaks)
        meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
        studies.append({
            "id": sid, "name": sid, "metadata": meta,
            "analyses": [{"id": sid, "name": "1", "metadata": meta,
                          "points": [{"space": "MNI",
                                      "coordinates": [float(r.x), float(r.y), float(r.z)],
                                      "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                                     for r in sub.itertuples()]}],
        })
    return Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img)


hot = np.abs(truth) >= np.percentile(np.abs(truth), 75)
print(f"{'peaks/study':>12} {'peak_bias':>10} {'mean |g| hot':>13} {'truth hot':>10} "
      f"{'ratio':>7} {'r (all cov)':>12}")
for max_peaks in (5, 15, 40):
    for peak_bias in (None, "per-study"):
        estimator = CBES(fwhm=10.0, mask=masker, peak_bias=peak_bias, null_method="none",
                         threshold="reporting_threshold")
        result = estimator.fit(build(max_peaks))
        g = result.get_map("g", return_type="array").ravel()
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        both = covered & hot & np.isfinite(g)
        r = stats.pearsonr(np.abs(g[covered]), np.abs(truth[covered]))[0] if covered.sum() > 10 else np.nan
        print(f"{max_peaks:12d} {str(peak_bias):>10} {np.abs(g[both]).mean():13.3f} "
              f"{np.abs(truth[both]).mean():10.3f} "
              f"{np.abs(g[both]).mean() / np.abs(truth[both]).mean():7.3f} {r:12.3f}",
              flush=True)
