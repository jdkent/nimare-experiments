"""Are coordinate studies over- or under-weighted against image studies in a mixed fit?

v2. The first attempt was confounded two ways and its number meant nothing:

  * the coordinate "actual error" used rho = 1/null_peak_mean_g, which is only defined up to
    the very constant the coordinates cannot identify, so it was dominated by a known scale
    offset rather than by the spread being measured;
  * "claimed variance" omitted tau^2, which matters because real between-study heterogeneity
    is part of the error a study makes, and it is the same for both kinds of study.

Both are fixed here. rho is calibrated by matching means against the truth before any error is
taken, so only the residual spread is measured; tau^2 is estimated from the images and added to
both sides. The comparison is then like for like:

    claimed   (var_k + tau^2) / w        the effective variance behind the model's weight
    actual    (g_k - truth_LOO)^2        the error the study actually makes there

ratio > 1 means the study claims more precision than its errors justify.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.meta.utils import get_ale_kernel
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import mm2vox

U, FWHM = 3.2905, 10.0
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
V = np.array([v_maps[k] for k in keys])
W = 1.0 / V

def loo(exclude):
    keep = np.array([k != exclude for k in keys])
    w = W[keep]
    return (G[keep] * w).sum(0) / np.clip(w.sum(0), 1e-12, None)

truths = {sid: loo(sid) for sid in keys}
pooled = (G * W).sum(0) / np.clip(W.sum(0), 1e-12, None)

# tau^2: DerSimonian-Laird at each voxel, from the images. Part of the error either kind of
# study makes, and identical for both, so it belongs on both sides of the comparison.
q = (W * (G - pooled) ** 2).sum(0)
denom = W.sum(0) - (W**2).sum(0) / np.clip(W.sum(0), 1e-12, None)
tau2 = np.clip((q - (len(keys) - 1)) / np.clip(denom, 1e-12, None), 0.0, None)
print(f"tau^2 from the images: median {np.median(tau2[tau2 > 0]):.4f}")

coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates

kern = get_ale_kernel(masker.mask_img, fwhm=FWHM)[1]
kern = kern / kern.max()
half = (np.array(kern.shape) - 1) // 2
sel = kern >= 0.01
offs = np.array(np.nonzero(sel)).T - half
wts = kern[sel]

shape = masker.mask_img.shape[:3]
mask_flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[mask_flat] = np.arange(mask_flat.size)

# Pass 1: collect every (focus, voxel) contribution, uncalibrated.
rows = []
for sid, sub in coords.groupby("id"):
    sid = str(sid)
    if sid not in g_maps:
        continue
    truth = truths[sid]
    n = float(sizes[sid])
    z = np.abs(sub["z_stat"].astype(float).to_numpy())
    g_k, var_k = peak_stat_to_hedges_g(z, np.full(len(z), n), stat_type="z")
    ijk = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    for centre, g_focus, var_focus in zip(ijk, np.abs(g_k), var_k):
        vox = centre + offs
        inside = np.all((vox >= 0) & (vox < shape), axis=1)
        cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
        keep = cols >= 0
        cols, w = cols[keep], wts[inside][keep]
        if cols.size:
            rows.append((g_focus, var_focus, w, np.abs(truth[cols]), tau2[cols]))

# Calibrate rho by matching means, so what follows measures spread and not the known offset.
num = sum(np.sum(w * t) for _, _, w, t, _ in rows)
den = sum(g * np.sum(w) for g, _, w, _, _ in rows)
rho = float(num / den)
print(f"calibrated rho (means matched): {rho:.3f}")

claimed_c = sum(np.sum(w * ((var + t2) / w)) for _, var, w, _, t2 in rows)
actual_c = sum(np.sum(w * (rho * g - t) ** 2) for g, _, w, t, t2 in rows)
weight_c = sum(np.sum(w) for _, _, w, _, _ in rows)

claimed_i, actual_i, weight_i = 0.0, 0.0, 0.0
for position, sid in enumerate(keys):
    truth = truths[sid]
    ok = np.isfinite(V[position]) & (V[position] < np.inf) & (truth != 0)
    claimed_i += np.sum(V[position][ok] + tau2[ok])
    actual_i += np.sum((G[position][ok] - truth[ok]) ** 2)
    weight_i += ok.sum()

ratio_i = (actual_i / weight_i) / (claimed_i / weight_i)
ratio_c = (actual_c / weight_c) / (claimed_c / weight_c)
print(f"\n{'kind':>12s} {'claimed var':>12s} {'actual err^2':>13s} {'over-weighted':>14s}")
print(f"{'images':>12s} {claimed_i/weight_i:12.5f} {actual_i/weight_i:13.5f} {ratio_i:14.2f}")
print(f"{'coordinates':>12s} {claimed_c/weight_c:12.5f} {actual_c/weight_c:13.5f} {ratio_c:14.2f}")
print(f"\ncoordinates over-weighted relative to images: {ratio_c/ratio_i:.2f}x")
