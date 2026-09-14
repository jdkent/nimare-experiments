"""With only one image, does it matter *which* one?

One image fixes the scale well -- 2.04x down to 1.42x on the pain collection -- because the
scale is a single scalar and one map supplies thousands of constraints on it. But it is one
study's scalar. If that study's effect is atypical, every coordinate study is rescaled to its
idiosyncrasy, and the caller has no way to know: they contributed the image they happened to
have.

Each of the 21 pain studies takes its turn as the sole image, with the other 20 contributing
coordinates only. The spread across those 21 fits is the risk a one-image caller is running.
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
base = json.loads(json.dumps(full.to_dict()))
ids = [st["id"] for st in base["studies"]]
print(f"{len(ids)} studies; each takes its turn as the sole image supplier\n")
print(f"{'image study':>16s} {'r_top':>7s} {'mag ratio':>10s}")

ratios = []
for keep_id in ids:
    obj = json.loads(json.dumps(base))
    for study in obj["studies"]:
        for analysis in study["analyses"]:
            if study["id"] != keep_id:
                analysis["images"] = []
    try:
        est = CBES(fwhm=10.0, null_method="none", threshold="study-min",
                   peak_bias="per-study", peak_bias_scale="auto", use_images=True,
                   ).fit(Studyset(obj))
    except Exception as exc:
        print(f"{str(keep_id)[:16]:>16s} {'-':>7s} {type(exc).__name__:>10s}")
        continue
    g_est = est.get_map("g", return_type="array").ravel()
    n = est.get_map("n_studies", return_type="array").ravel()
    scored = signal & (n > 0)
    rho = stats.spearmanr(np.abs(g_est[scored]), np.abs(truth[scored])).statistic
    mag = float(np.median(np.abs(g_est[scored]) / np.maximum(np.abs(truth[scored]), 1e-6)))
    ratios.append(mag)
    print(f"{str(keep_id)[:16]:>16s} {rho:7.3f} {mag:10.3f}")

if ratios:
    a = np.array(ratios)
    print(f"\nmag ratio across choice of image: median {np.median(a):.3f}, "
          f"range {a.min():.3f}-{a.max():.3f}, spread factor {a.max()/a.min():.2f}x")
