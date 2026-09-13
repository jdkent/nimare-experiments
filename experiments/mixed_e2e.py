"""End-to-end: does the threshold rule matter, on a collection that thresholded unevenly?

Each of the 21 NIDM pain studies gets its own randomly assigned threshold, so the collection
looks like a literature search rather than a consortium. A studyset is rebuilt from the peaks
each study would then have reported, carrying its true threshold as metadata. The images are
held out as the reference: the inverse-variance pooled g map over all 21, which is what the
coordinate fit is trying to recover.

  r      -- correlation with the image reference, over covered voxels (is the shape right?)
  ratio  -- mean |g| coordinates / mean |g| images (is the scale right? 1.0 is correct)
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import pandas as pd
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
rng = np.random.default_rng(0)
assigned = {str(sid): float(rng.choice(THRESHOLDS)) for sid in ss.ids}
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))

# --- image reference: inverse-variance pooled g over all 21 studies -------------------
g_stack, w_stack = [], []
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_stack.append(np.where(ok, g, 0.0))
    w_stack.append(np.where(ok, 1.0 / v, 0.0))
g_stack, w_stack = np.array(g_stack), np.array(w_stack)
reference = (g_stack * w_stack).sum(0) / np.clip(w_stack.sum(0), 1e-12, None)

# --- mixed-threshold coordinates -----------------------------------------------------
frames = []
for u in THRESHOLDS:
    df = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=u, two_sided=True, remove_subpeaks=True
    ).transform(ss).coordinates
    frames.append(df[[assigned.get(str(s)) == u for s in df["id"]]])
coords = pd.concat(frames, ignore_index=True)

studies = []
for sid, sub in coords.groupby("id"):
    sid = str(sid)
    meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": assigned[sid]}
    studies.append({
        "id": sid, "name": sid, "metadata": meta,
        "analyses": [{
            "id": f"{sid}-1", "name": "1", "metadata": meta,
            "points": [
                {"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                 "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                for r in sub.itertuples()
            ],
        }],
    })
mixed = Studyset({"id": "mixed", "name": "mixed", "studies": studies})
print(f"mixed collection: {len(studies)} studies, {len(coords)} foci, "
      f"{coords.groupby('id').size().median():.0f} foci/study (median); "
      f"thresholds {sorted(set(assigned.values()))}")

def evaluate(threshold, peak_bias):
    est = CBES(fwhm=10.0, null_method="parametric", use_images=False,
               threshold=threshold, peak_bias=peak_bias, mask=masker)
    fitted = est.fit(mixed).get_map("g", return_type="array").ravel()
    ok = (fitted != 0) & (reference != 0) & np.isfinite(fitted) & np.isfinite(reference)
    cutoffs = est._cutoffs_z_
    shared = [i for i in cutoffs.index if str(i) in assigned]
    u_err = np.mean([abs(cutoffs[i] - assigned[str(i)]) for i in shared])
    return (stats.pearsonr(fitted[ok], reference[ok])[0],
            np.abs(fitted[ok]).mean() / np.abs(reference[ok]).mean(), u_err)

print(f"\n{'threshold rule':>22s} {'peak_bias':>12s} {'r vs images':>12s} "
      f"{'scale ratio':>12s} {'mean |u_hat-u|':>15s}")
for threshold in ("pooled-min", "study-min", "study-min-corrected", "reporting_threshold"):
    for peak_bias in (None, "per-study"):
        try:
            r, ratio, u_err = evaluate(threshold, peak_bias)
            print(f"{threshold:>22s} {str(peak_bias):>12s} {r:12.3f} {ratio:12.2f} "
                  f"{u_err:15.3f}")
        except Exception as exc:  # noqa: BLE001
            print(f"{threshold:>22s} {str(peak_bias):>12s}  failed: {type(exc).__name__}: {exc}")
