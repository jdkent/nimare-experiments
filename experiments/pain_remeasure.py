"""Re-measure the pain inflation with a metric that cannot divide by near-zero.

The figure quoted all session -- 2.65x -- is the ratio of mean |g| over every covered voxel to
the image reference over the same voxels. On simulated data where the truth is known that same
metric returns 80 to 182x, because most covered voxels hold no effect and the denominator
collapses. It measures how much true-zero is in the covered set, not how inflated the estimate
is.

Scored instead where the reference says there is something to estimate.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))

# Image reference: inverse-variance pooled g over all 21 studies.
gs, ws = [], []
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    gs.append(np.where(ok, g, 0.0)); ws.append(np.where(ok, 1.0 / v, 0.0))
gs, ws = np.array(gs), np.array(ws)
reference = (gs * ws).sum(0) / np.clip(ws.sum(0), 1e-12, None)

coords = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                             remove_subpeaks=True).transform(ss).coordinates
studies = []
for sid, sub in coords.groupby("id"):
    sid = str(sid)
    meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
    studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [{
        "id": sid, "name": "1", "metadata": meta,
        "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                    "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                   for r in sub.itertuples()]}]})
cs = Studyset({"id": "pain", "name": "pain", "studies": studies}, target=None,
              mask=masker.mask_img)

fitted = CBES(fwhm=10.0, null_method="none", peak_bias="per-study",
              threshold="reporting_threshold", mask=masker).fit(cs)
g = fitted.get_map("g", return_type="array").ravel()
covered = g != 0
ref = np.abs(reference)

print(f"covered voxels: {covered.mean():.2f} of the brain")
print(f"{'scored over':>40s} {'voxels':>8s} {'estimate':>9s} {'reference':>10s} {'ratio':>7s}")
rows = [("every covered voxel (the quoted metric)", covered),
        ("reference above its median", covered & (ref >= np.median(ref[covered]))),
        ("reference in its top quartile", covered & (ref >= np.percentile(ref[covered], 75))),
        ("reference in its top decile", covered & (ref >= np.percentile(ref[covered], 90)))]
for label, sel in rows:
    if sel.sum() < 50:
        continue
    est, tru = np.abs(g[sel]).mean(), ref[sel].mean()
    print(f"{label:>40s} {sel.sum():8d} {est:9.3f} {tru:10.3f} {est/tru:7.2f}")

# A ratio of means is fragile; a paired median ratio is not.
strong = covered & (ref >= np.percentile(ref[covered], 75))
print(f"\nmedian of the per-voxel ratio, top quartile: "
      f"{np.median(np.abs(g[strong]) / np.clip(ref[strong], 1e-9, None)):.2f}")
