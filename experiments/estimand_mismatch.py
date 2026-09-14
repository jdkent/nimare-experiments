"""Is the pain 'overestimate' an estimand mismatch rather than a bias?

CBES's ``g`` is mu(v): the effect among the studies that have an effect at v. The
inverse-variance pooled image map averages every study, including those with nothing there, so
it estimates the marginal pi(v) * mu(v). Comparing one against the other should show a ratio
near 1/pi, not near 1.

Tested by taking CBES's own prevalence map and asking whether ``g * prevalence`` -- the
marginal CBES implies -- lines up with the pooled image truth, where plain ``g`` did not.
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
coords = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                             remove_subpeaks=True).transform(ss).coordinates

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
pooled = np.sum(G * W, axis=0) / np.maximum(np.sum(W, axis=0), 1e-12)

ids = sorted(set(str(i) for i in coords["id"]))
for max_peaks in (10, 40):
    studies = []
    for sid in ids:
        sub = coords[coords["id"].astype(str) == sid].copy()
        sub["_a"] = np.abs(sub["z_stat"].astype(float))
        sub = sub.sort_values("_a", ascending=False).head(max_peaks)
        meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [
            {"id": sid, "name": "1", "metadata": meta,
             "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                         "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                        for r in sub.itertuples()]}]})
    studyset = Studyset({"id": "pain", "name": "pain", "studies": studies},
                        target=None, mask=masker.mask_img)
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="reporting_threshold")
    r = est.fit(studyset)
    g = r.get_map("g", return_type="array").ravel()
    pi = r.get_map("prevalence", return_type="array").ravel()
    cov = r.get_map("n_studies", return_type="array").ravel() > 0
    hot = np.abs(pooled) >= np.percentile(np.abs(pooled), 75)
    both = cov & hot & np.isfinite(g)
    marginal = g * pi
    print(f"\n--- {max_peaks} peaks/study, {cov.sum()} covered voxels ---")
    print(f"  pooled image truth (hot):      {np.abs(pooled[both]).mean():.3f}")
    print(f"  CBES g (hot):                  {np.abs(g[both]).mean():.3f}  "
          f"ratio {np.abs(g[both]).mean()/np.abs(pooled[both]).mean():.2f}")
    print(f"  CBES prevalence (hot):         {pi[both].mean():.3f}  "
          f"implied 1/pi {1.0/max(pi[both].mean(), 1e-6):.2f}")
    print(f"  CBES g * prevalence (hot):     {np.abs(marginal[both]).mean():.3f}  "
          f"ratio {np.abs(marginal[both]).mean()/np.abs(pooled[both]).mean():.2f}")
    print(f"  r(|g|, |truth|) over covered:            "
          f"{stats.pearsonr(np.abs(g[cov]), np.abs(pooled[cov]))[0]:.3f}")
    print(f"  r(|g*prevalence|, |truth|) over covered: "
          f"{stats.pearsonr(np.abs(marginal[cov]), np.abs(pooled[cov]))[0]:.3f}")
