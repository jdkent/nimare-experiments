"""How does accuracy depend on peaks per study? The realistic range is 3-20.

Every magnitude figure so far used the pain collection's full peak tables -- about 130 peaks
per study, which no real paper reports. A typical coordinate table has 3 to 20. That matters
for more than sample size: threshold="study-min" infers each study's reporting threshold from
its smallest reported peak, so with three peaks the order-statistic correction is working from
almost nothing, and the per-study correction rho_k is a function of that threshold.

Peaks are kept strongest-first, so a low cap is what a paper reporting only its top peaks
gives. Scored as before: inverse-variance truth from all 21 images, magnitude on the truth's
top quartile, restricted to the voxels CBES covers.
"""
import json, os, sys, warnings; warnings.simplefilter("ignore")
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
total = weight = None
for gp, vp in zip(ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel(); v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    gg = np.where(ok, g, 0.0)
    total = w * gg if total is None else total + w * gg
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))

full = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                           remove_subpeaks=True).transform(ss)
coords = full.coordinates.copy()
coords["_abs"] = coords["z_stat"].abs()
per_study = coords.groupby("id").size()
print(f"full tables: {len(coords)} peaks, {coords['id'].nunique()} studies, "
      f"median {int(per_study.median())} per study\n")

base = json.loads(json.dumps(full.to_dict()))
print(f"{'peaks/study':>12s} {'total':>6s} {'covered':>8s} {'rho':>7s} {'mag ratio':>10s}")
for cap in (3, 5, 10, 20, 10**6):
    keep = (coords.sort_values("_abs", ascending=False)
            .groupby("id", sort=False).head(cap).index)
    trimmed = coords.loc[sorted(keep)]
    by_analysis = {k: v for k, v in trimmed.groupby("id")}
    obj = json.loads(json.dumps(base))
    for study in obj["studies"]:
        for analysis in study["analyses"]:
            rows = by_analysis.get(analysis["id"], by_analysis.get(study["id"]))
            analysis["points"] = [] if rows is None else [
                {"space": "MNI",
                 "coordinates": [float(r["x"]), float(r["y"]), float(r["z"])],
                 "values": [{"kind": "Z", "value": float(r["z_stat"])}]}
                for _, r in rows.iterrows()]
            analysis["images"] = []
    est = CBES(fwhm=10.0, null_method="none", threshold="study-min", peak_bias="per-study",
               use_images=False).fit(Studyset(obj))
    g_est = est.get_map("g", return_type="array").ravel()
    n = est.get_map("n_studies", return_type="array").ravel()
    scored = signal & (n > 0)
    if scored.sum() < 50:
        print(f"{str(cap):>12s} {len(trimmed):6d} {(n>0).mean():8.1%} {'-':>7s} {'-':>10s}")
        continue
    rho = stats.spearmanr(np.abs(g_est[scored]), np.abs(truth[scored])).statistic
    mag = float(np.median(np.abs(g_est[scored]) / np.maximum(np.abs(truth[scored]), 1e-6)))
    print(f"{('all' if cap > 10**5 else str(cap)):>12s} {len(trimmed):6d} "
          f"{(n>0).mean():8.1%} {rho:7.3f} {mag:10.3f}")
