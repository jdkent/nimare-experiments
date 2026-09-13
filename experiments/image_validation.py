"""Can CBES recover, from thresholded peaks alone, the effect-size map you get from the images?

Reference: a voxelwise random-effects pooling of per-study Hedges' g, computed from the full
t images. Test: threshold those same images, keep only the peak coordinates and their z values,
and run CBES on that. Anything CBES recovers is recovered from ~1% of the information.
"""
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import nimare
from nimare.meta.cbma import ALE, CBES, MKDADensity
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.tests.utils import get_test_data_path
from nimare.transforms import ImagesToCoordinates, ImageTransformer

warnings.simplefilter("ignore")

dset = nimare.dataset.Dataset(os.path.join(get_test_data_path(), "nidm_pain_dset.json"))
dset.update_path("/root/.nimare/nidm_21pain")
ss = nimare.studyset.normalize_collection(dset)
ss = ImageTransformer(target="z").transform(ss)

masker = ss.masker
n_studies = len(ss.ids)
sample_sizes = np.asarray(ss.sample_sizes(), dtype=float)
print(f"{n_studies} studies, sample sizes {np.nanmin(sample_sizes):.0f}-{np.nanmax(sample_sizes):.0f}")

# ---------------------------------------------------------------- reference from full images
images = ss.images
t_col = [c for c in images.columns if c == "t"][0]
g_stack, var_stack = [], []
for i, path in enumerate(images[t_col].values):
    if path is None or not os.path.isfile(str(path)):
        continue
    t_data = masker.transform(str(path)).ravel()
    g, var_g = peak_stat_to_hedges_g(t_data, np.full(t_data.shape, sample_sizes[i]),
                                     stat_type="t", design="one-sample")
    g_stack.append(g)
    var_stack.append(var_g)

g_stack = np.vstack(g_stack)
var_stack = np.vstack(var_stack)
print(f"reference built from {g_stack.shape[0]} t images over {g_stack.shape[1]} voxels")

# Voxelwise DerSimonian-Laird pooling of the per-study g images.
w = 1.0 / var_stack
fixed = (w * g_stack).sum(0) / w.sum(0)
q = (w * (g_stack - fixed) ** 2).sum(0)
c = w.sum(0) - (w**2).sum(0) / w.sum(0)
tau2 = np.clip((q - (g_stack.shape[0] - 1)) / np.where(c > 0, c, np.nan), 0, None)
tau2 = np.nan_to_num(tau2)
w_re = 1.0 / (var_stack + tau2)
reference_g = (w_re * g_stack).sum(0) / w_re.sum(0)
reference_z = reference_g * np.sqrt(w_re.sum(0))
print(f"reference g: mean {reference_g.mean():+.3f}, max {reference_g.max():.3f}")

# ------------------------------------------------------- coordinates from thresholded images
coord_ss = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=3.2905, two_sided=True, remove_subpeaks=True
).transform(ss)
coords = coord_ss.coordinates
print(f"\nthresholded peaks: {len(coords)} foci from {coords['id'].nunique()} studies "
      f"({len(coords) / max(g_stack.shape[1], 1) * 100:.2f}% of voxels retained)")

results = {}
support = None
for name, est in [
    ("CBES none", CBES(fwhm=10.0, selection_model="none", null_method="parametric")),
    ("CBES tobit", CBES(fwhm=10.0, selection_model="tobit", null_method="parametric")),
    ("CBES zero-inflated", CBES(fwhm=10.0, selection_model="zero-inflated", null_method="parametric")),
    ("ALE", ALE()),
    ("MKDADensity", MKDADensity()),
]:
    res = est.fit(coord_ss)
    if name.startswith("CBES"):
        results[name] = res.get_map("g", return_type="array").ravel()
        if "g_marginal" in res.maps:
            results[name + " (marginal)"] = res.get_map("g_marginal", return_type="array").ravel()
        if support is None:
            support = res.get_map("n_studies", return_type="array").ravel()
    else:
        results[name] = res.get_map("z", return_type="array").ravel()

well_supported = support >= 3
print(f"voxels with >=3 contributing studies: {well_supported.sum()} of {support.size}")

def slope(est_vals, ref):
    """Regression of reference on estimate: 1.0 means calibrated, <1 means inflated."""
    ok = np.isfinite(est_vals) & np.isfinite(ref)
    if ok.sum() < 10 or np.allclose(est_vals[ok], 0):
        return np.nan
    return np.polyfit(est_vals[ok], ref[ok], 1)[0]

print(f"\n{'estimator':30s} {'r':>7s} {'rho':>7s} {'r|sup':>7s} {'rho|sup':>8s} "
      f"{'slope':>7s} {'mean|ref>.2':>12s}")
big = reference_g > 0.2
for name, arr in results.items():
    r = stats.pearsonr(arr, reference_g)[0]
    rho = stats.spearmanr(arr, reference_g)[0]
    r_s = stats.pearsonr(arr[well_supported], reference_g[well_supported])[0]
    rho_s = stats.spearmanr(arr[well_supported], reference_g[well_supported])[0]
    sl = slope(arr[well_supported], reference_g[well_supported])
    sel = big & well_supported
    mean_here = arr[sel].mean() if sel.sum() else np.nan
    print(f"{name:30s} {r:7.3f} {rho:7.3f} {r_s:7.3f} {rho_s:8.3f} {sl:7.3f} {mean_here:12.3f}")
print(f"{'(reference)':30s} {'':7s} {'':7s} {'':7s} {'':8s} {'':7s} "
      f"{reference_g[big & well_supported].mean():12.3f}")
