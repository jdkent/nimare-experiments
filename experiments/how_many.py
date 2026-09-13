"""How many studies with images are needed to calibrate peak_bias?

Because the fitted map scales exactly linearly in rho, the relative sampling error of
rho IS the relative error of the reported effect sizes. So the question reduces to: how
precisely is rho estimated from k calibration studies?
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ids = np.array(list(cs.ids))

def restrict(s, keep): return s.slice(analyses=[a for a in s.ids if a in set(keep)])
def fit(s, use_images):
    return CBES(fwhm=10.0, null_method="parametric",
                use_images=use_images).fit(s).get_map("g", return_type="array").ravel()

def rho_from(subset):
    sub = restrict(cs, subset)
    im, co = fit(sub, True), fit(sub, False)
    ok = (im != 0) & (co != 0)
    return float(im[ok].mean() / co[ok].mean()) if ok.sum() > 100 else np.nan

KS = [2, 3, 4, 5, 8, 12, 16]
REPS = 12
rng = np.random.default_rng(0)
results = {k: [] for k in KS}
for rep in range(REPS):
    pool = rng.permutation(ids)
    for k in KS:
        r = rho_from(pool[:k])
        if np.isfinite(r):
            results[k].append(r)
    print(f"  rep {rep} done", flush=True)

print(f"\n{'images used':>11s} {'mean rho':>9s} {'sd':>7s} {'rel. sd':>8s} "
      f"{'implied error in g':>19s}")
for k in KS:
    v = np.array(results[k])
    if len(v) < 3:
        continue
    print(f"{k:11d} {v.mean():9.3f} {v.std(ddof=1):7.3f} {v.std(ddof=1)/v.mean()*100:7.1f}% "
          f"{'+-' + format(v.std(ddof=1)/v.mean()*100, '.0f') + '%':>19s}")
