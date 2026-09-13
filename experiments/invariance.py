"""Invariance test: swapping a study from coordinates to images should not move the answer.

An estimator whose result depends on how many studies happen to have shared images is
incoherent, whatever its correlation with the truth. This sweeps k = images used and asks
whether the pooled effect is flat in k.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ref = np.load("/tmp/claude-0/cmp/reference_g.npy")
ids = list(cs.ids)
rng = np.random.default_rng(0)
order = list(rng.permutation(len(ids)))

def fit(k, peak_bias):
    keep = {ids[i] for i in order[:k]}
    est = CBES(fwhm=10.0, null_method="parametric", use_images=k > 0, peak_bias=peak_bias)
    original = est._load_image_studies
    est._load_image_studies = lambda ds: {a: b for a, b in original(ds).items() if a in keep}
    g = est.fit(cs).get_map("g", return_type="array").ravel()
    cov = g != 0
    sel = cov & (ref > 0.2)
    return g[sel].mean(), stats.pearsonr(g[cov], ref[cov])[0]

for peak_bias in (None, 0.479):
    label = "uncorrected" if peak_bias is None else f"peak_bias={peak_bias}"
    print(f"\n{label}   (reference mean where ref>0.2 = {ref[ref>0.2].mean():.3f})")
    print(f"{'images used':>12s} {'mean g':>8s} {'r':>7s}")
    means = []
    for k in (0, 2, 5, 10, 15, 21):
        m, r = fit(k, peak_bias)
        means.append(m)
        print(f"{k:12d} {m:8.3f} {r:7.3f}", flush=True)
    print(f"  spread across k: {max(means)-min(means):.3f}")
