"""Does shrinking the prevalence help on a real collection, or only in simulation?

The use case: 21 NIDM pain studies whose full images are available, so the answer is known.
CBES is run on peaks thresholded out of those same images -- the standing validation setup --
and scored against the inverse-variance map from all 21, on the truth's top quartile. This is
the identical scoring used for the SDM comparison, so the numbers are comparable to it.

The documented failure being tested is the magnitude: CBES reaches rho = 0.84 against these
images but overestimates the effect about twofold. If shrinking the prevalence does not move
that, and does not improve the ranking, then it is not useful here and should be discarded.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
total = weight = None
for gp, vp in zip(ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    gg = np.where(ok, g, 0.0)
    total = w * gg if total is None else total + w * gg
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))
print(f"truth from 21 images; {int(signal.sum())} voxels in the top quartile")

peaks = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                            remove_subpeaks=True).transform(ss)
print(f"\n{'prior':>7s} {'rho(top)':>9s} {'mag ratio g':>12s} {'mag ratio marg':>15s} "
      f"{'median pi':>10s}")
print(f"{'(want)':>7s} {'high':>9s} {'1.000':>12s} {'1.000':>15s} {'':>10s}")
for prior in (0.0, 0.1, 0.25, 0.5, 1.0):
    est = CBES(fwhm=10.0, null_method="none", threshold="study-min", peak_bias="per-study",
               use_images=False, prevalence_prior=prior, max_iter=200).fit(peaks)
    g = est.get_map("g", return_type="array").ravel()
    m = est.get_map("g_marginal", return_type="array").ravel()
    pi = est.get_map("prevalence", return_type="array").ravel()
    rho = stats.spearmanr(np.abs(g[signal]), np.abs(truth[signal])).statistic
    mag = float(np.median(np.abs(g[signal]) / np.maximum(np.abs(truth[signal]), 1e-6)))
    magm = float(np.median(np.abs(m[signal]) / np.maximum(np.abs(truth[signal]), 1e-6)))
    print(f"{prior:7.2f} {rho:9.3f} {mag:12.3f} {magm:15.3f} "
          f"{np.median(pi[pi > 0]):10.3f}")
