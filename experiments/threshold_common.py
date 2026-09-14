"""Does "study-min" cost anything when every study really did share a threshold?

radius_threshold.py showed "study-min" beating "pooled-min" decisively when studies threshold
differently, which is the case it exists for. Before changing a default, the other case has to
be checked: with one common threshold, "pooled-min" is very nearly correct by construction,
while "study-min" infers each study's cut separately and so pays estimation noise for nothing.

Same scoring as radius_threshold.py. The "true" row hands over the threshold that actually
generated the peaks, which is the ceiling here rather than an approximation.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.studyset import Studyset
from nimare.transforms import ImagesToCoordinates, ImageTransformer

U = 3.2905
ss_full = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss_full.masker

total = weight = None
for gp, vp in zip(ss_full.images["g"], ss_full.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))

sizes = {r.study_id: float(np.mean(r.sample_sizes))
         for r in ss_full.metadata.itertuples() if r.sample_sizes}
es.CBES._load_image_studies = lambda self, dataset: {}

table = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                            remove_subpeaks=True).transform(ss_full).coordinates


def build(cap):
    studies = []
    for sid, rows in table.groupby("study_id"):
        rows = rows.assign(_abs=rows["z_stat"].astype(float).abs()).nlargest(cap, "_abs")
        if not len(rows):
            continue
        n = sizes.get(sid, 20.0)
        studies.append({"id": sid, "name": sid, "metadata": {"sample_sizes": [n]},
                        "analyses": [{"id": f"{sid}-1", "name": "1",
                                      "metadata": {"sample_sizes": [n]},
                                      "points": [
                                          {"space": "MNI",
                                           "coordinates": [float(r.x), float(r.y), float(r.z)],
                                           "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                                          for r in rows.itertuples()],
                                      "images": []}]})
    return Studyset({"id": "p", "name": "p", "studies": studies})


def score(studyset, threshold):
    result = CBES(fwhm=10.0, mask=masker.mask_img, null_method="none",
                  peak_bias="per-study", threshold=threshold).fit(studyset)
    g = result.get_map("g", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    prev = result.get_map("prevalence", return_type="array").ravel()
    covered = k > 0
    both = covered & signal & np.isfinite(g)
    ratios = np.abs(g[both]) / np.maximum(np.abs(truth[both]), 1e-9)
    return (covered.mean(),
            stats.spearmanr(np.abs(g[both]), np.abs(truth[both])).statistic,
            float(np.median(ratios)), float(np.median(prev[covered])))


print(f"one common threshold at z = {U}\n")
print(f"{'peaks':>6s} {'rule':>12s} {'covered':>8s} {'r_top':>7s} {'mag':>7s} {'prev':>7s}")
for cap in (10, 20, 10_000):
    studyset = build(cap)
    label = "all" if cap > 1000 else str(cap)
    for name, rule in (("pooled-min", "pooled-min"), ("study-min", "study-min"), ("true", U)):
        cov, r, mag, prev = score(studyset, rule)
        print(f"{label:>6s} {name:>12s} {100 * cov:7.1f}% {r:7.3f} {mag:7.2f} {prev:7.3f}",
              flush=True)
