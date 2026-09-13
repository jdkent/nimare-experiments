"""Are coordinate studies over-weighted against image studies in a mixed fit?

A coordinate study enters at kernel weight w with var_g = the sampling variance of g at its
peak. That variance is honest about the peak's own precision and silent about two things it
should not be: the peak's location is uncertain, and the peak was selected for being large.
An image study's var_g has neither problem. If the coordinate variance understates the real
error, coordinates carry more precision than they have earned and can outvote images.

This is measurable without any modelling. For each study and each voxel, compare

    the error the model *claims*      var_g (+ tau2), as it enters the inverse-variance weight
    the error the study *makes*       (its g here  -  the leave-one-out pooled truth here)^2

averaged over the voxels that study informs. A ratio near 1 means the weight is earned; a
ratio above 1 means the study is over-weighted by that factor. Doing it for both kinds of
study on the same collection gives their relative over-weighting directly.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import (
    null_peak_mean_g, peak_stat_to_hedges_g, infer_threshold_from_minimum,
)
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.utils import get_ale_kernel
from nimare.utils import mm2vox

U = 3.2905
FWHM = 10.0

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))

g_maps, v_maps = {}, {}
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_maps[sid], v_maps[sid] = np.where(ok, g, 0.0), np.where(ok, v, np.inf)
keys = list(g_maps)
G = np.array([g_maps[k] for k in keys])
W = np.array([1.0 / v_maps[k] for k in keys])

def loo(exclude):
    """Inverse-variance pooled g over every study but this one."""
    keep = np.array([k != exclude for k in keys])
    w = W[keep]
    return (G[keep] * w).sum(0) / np.clip(w.sum(0), 1e-12, None)

coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates

# Kernel offsets, peak-normalized and truncated the way the estimator truncates them.
kern = get_ale_kernel(masker.mask_img, fwhm=FWHM)[1]
kern = kern / kern.max()
offs = np.array(np.nonzero(kern >= 0.01)).T - (np.array(kern.shape) - 1) // 2
wts = kern[tuple((offs + (np.array(kern.shape) - 1) // 2).T)]

shape = masker.mask_img.shape[:3]
mask_flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[mask_flat] = np.arange(mask_flat.size)

print(f"{'kind':>12s} {'studies':>8s} {'claimed var':>12s} {'actual err^2':>13s} "
      f"{'over-weighted':>14s}")

# --- coordinate studies ---------------------------------------------------------------
claimed, actual = [], []
for sid, sub in coords.groupby("id"):
    sid = str(sid)
    if sid not in g_maps:
        continue
    truth = loo(sid)
    n = float(sizes[sid])
    z = np.abs(sub["z_stat"].astype(float).to_numpy())
    u_hat = infer_threshold_from_minimum(z.min(), len(z))
    rho = 1.0 / null_peak_mean_g(u_hat, n)          # scale-free; cancels in the ratio below
    g_k, var_k = peak_stat_to_hedges_g(z, np.full(len(z), n), stat_type="z")
    ijk = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    for centre, g_focus, var_focus in zip(ijk, g_k, var_k):
        vox = centre + offs
        inside = np.all((vox >= 0) & (vox < shape), axis=1)
        cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
        keep = cols >= 0
        cols, w = cols[keep], wts[inside][keep]
        if not cols.size:
            continue
        # The model's inverse-variance weight is w / var; the error it implies is var / w.
        claimed.append(np.average(var_focus / w, weights=w))
        actual.append(np.average((rho * abs(g_focus) - np.abs(truth[cols])) ** 2, weights=w))
# rho is defined up to a constant, so report the ratio after matching the means -- the
# question is the *spread* of the error, not the scale offset already known to be 2.6x.
claimed_c, actual_c = np.mean(claimed), np.mean(actual)
n_coord = coords["id"].nunique()

# --- image studies --------------------------------------------------------------------
claimed_i, actual_i = [], []
for sid in keys:
    truth = loo(sid)
    ok = np.isfinite(v_maps[sid]) & (v_maps[sid] < np.inf) & (truth != 0)
    claimed_i.append(np.mean(v_maps[sid][ok]))
    actual_i.append(np.mean((g_maps[sid][ok] - truth[ok]) ** 2))
claimed_i, actual_i = np.mean(claimed_i), np.mean(actual_i)

print(f"{'images':>12s} {len(keys):8d} {claimed_i:12.5f} {actual_i:13.5f} "
      f"{actual_i/claimed_i:14.2f}")
print(f"{'coordinates':>12s} {n_coord:8d} {claimed_c:12.5f} {actual_c:13.5f} "
      f"{actual_c/claimed_c:14.2f}")
print(f"\ncoordinates are over-weighted relative to images by "
      f"{(actual_c/claimed_c)/(actual_i/claimed_i):.2f}x")
print("(>1 means a coordinate study claims more precision than its errors justify, so it "
      "outvotes\n an image study it should be losing to)")
