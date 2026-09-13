"""Are coordinate studies correctly weighted against image studies? v3.

Two corrections to v2, both of which mattered:

* **rho^2 on the variance.** ``_apply_peak_bias`` rescales ``g`` by rho *and* ``var_g`` by
  rho^2. v2 shrank the value and not the variance, so it credited coordinates with ~20x more
  claimed variance than the estimator actually uses, and the d = 0 ratio of 0.227 was an
  artefact of that alone.
* **the image side must not be scored on voxels selected by its own peak.** v2 evaluated study
  k's image error at voxels where study k reported, where G_k is large by construction -- the
  image's own winner's curse. It made images look 4.5x miscalibrated at d = 0. Here the image
  baseline at those voxels comes from the *other* studies, which were not selected on.

What is left is a like-for-like comparison at the same voxels:

    claimed_c(d) = (rho^2 var_k + tau^2) / w(d)     what the model says a coordinate knows
    actual_c(d)  = (rho g_peak - truth_k(v))^2      what it actually gets wrong
    claimed_i    = mean_{j != k} (var_j(v) + tau^2)
    actual_i     = mean_{j != k} (G_j(v) - truth_j(v))^2

ratio > 1 means the study claims more precision than its errors justify. The mis-weighting is
ratio_c / ratio_i; below 1 means coordinates are under-weighted against images.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import pandas as pd
from load_pain import load_pain
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.meta.utils import get_ale_kernel
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import mm2vox

U, FWHM, MIN_W = 3.2905, 10.0, 0.01
CACHE = "/tmp/claude-0/cmp/pain_coords_cache.csv"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = {str(i): float(n) for i, n in zip(ss.ids, ss.sample_sizes())}

g_maps, v_maps = {}, {}
for sid, gp, vp in zip(ss.images["id"].astype(str), ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_maps[sid], v_maps[sid] = np.where(ok, g, 0.0), np.where(ok, v, np.inf)

keys = list(g_maps)
index = {k: i for i, k in enumerate(keys)}
G = np.array([g_maps[k] for k in keys])
V = np.array([v_maps[k] for k in keys])
OK = np.isfinite(V) & (V < np.inf)
W = np.where(OK, 1.0 / np.where(OK, V, 1.0), 0.0)
pooled = (G * W).sum(0) / np.clip(W.sum(0), 1e-12, None)

q = (W * (G - pooled) ** 2).sum(0)
den = W.sum(0) - (W**2).sum(0) / np.clip(W.sum(0), 1e-12, None)
tau2 = np.clip((q - (len(keys) - 1)) / np.clip(den, 1e-12, None), 0.0, None)

# Leave-one-out truth for every study, as a matrix.
TRUTH = np.empty_like(G)
for k, i in index.items():
    keep = np.arange(len(keys)) != i
    w = W[keep]
    TRUTH[i] = (G[keep] * w).sum(0) / np.clip(w.sum(0), 1e-12, None)

# Per-study image error and claimed variance, zeroed where the study has no usable value there.
ERR_I = np.where(OK, (G - TRUTH) ** 2, 0.0)
CLAIM_I = np.where(OK, np.where(OK, V, 0.0) + tau2[None, :], 0.0)
COUNT_I = OK.astype(float)

if os.path.exists(CACHE):
    coords = pd.read_csv(CACHE)
else:
    coords = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
    ).transform(ss).coordinates
    coords.to_csv(CACHE, index=False)
coords["id"] = coords["id"].astype(str)

kern = get_ale_kernel(masker.mask_img, fwhm=FWHM)[1]
kern = kern / kern.max()
half = (np.array(kern.shape) - 1) // 2
sel = kern >= MIN_W
offs = np.array(np.nonzero(sel)).T - half
wts = kern[sel]
vox_mm = np.abs(np.diag(masker.mask_img.affine)[:3])
dist = np.sqrt(((offs * vox_mm) ** 2).sum(1))

EDGES = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 26.0])
which = np.clip(np.digitize(dist, EDGES) - 1, 0, len(EDGES) - 2)
N_BINS = len(EDGES) - 1

shape = masker.mask_img.shape[:3]
mask_flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[mask_flat] = np.arange(mask_flat.size)

foci, peak_est, peak_truth = [], [], []
for sid, sub in coords.groupby("id"):
    if sid not in g_maps:
        continue
    n = sizes[sid]
    z = sub["z_stat"].astype(float).to_numpy()
    g_k, var_k = peak_stat_to_hedges_g(np.abs(z), np.full(len(z), n), stat_type="z")
    g_k = np.sign(z) * g_k
    ijk = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    for centre, gf, vf in zip(ijk, g_k, var_k):
        col = lookup[np.ravel_multi_index(np.clip(centre, 0, np.array(shape) - 1), shape)]
        if col < 0:
            continue
        foci.append((sid, centre, gf, vf))
        peak_est.append(abs(gf))
        peak_truth.append(abs(TRUTH[index[sid]][col]))

rho = float(np.mean(peak_truth) / np.mean(peak_est))
print(f"{len(foci)} foci; peak inflation rho = {rho:.3f} "
      f"(reported mean |g| {np.mean(peak_est):.3f} vs truth at the peak "
      f"{np.mean(peak_truth):.3f})\n", flush=True)

acc = {k: np.zeros(N_BINS) for k in ("n", "w", "claim_c", "act_c", "claim_i", "act_i")}
for sid, centre, gf, vf in foci:
    i = index[sid]
    t = TRUTH[i]
    vox = centre + offs
    inside = np.all((vox >= 0) & (vox < shape), axis=1)
    cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
    keep = cols >= 0
    cols, w, b = cols[keep], wts[inside][keep], which[inside][keep]
    if not cols.size:
        continue

    # Image baseline from the other studies only: never selected on its own peak.
    cnt = COUNT_I[:, cols].sum(0) - COUNT_I[i, cols]
    good = cnt > 0
    if not good.any():
        continue
    cols, w, b, cnt = cols[good], w[good], b[good], cnt[good]
    err_i = (ERR_I[:, cols].sum(0) - ERR_I[i, cols]) / cnt
    cl_i = (CLAIM_I[:, cols].sum(0) - CLAIM_I[i, cols]) / cnt

    np.add.at(acc["n"], b, 1.0)
    np.add.at(acc["w"], b, w)
    np.add.at(acc["claim_c"], b, (rho**2 * vf + tau2[cols]) / w)
    np.add.at(acc["act_c"], b, (rho * gf - t[cols]) ** 2)
    np.add.at(acc["claim_i"], b, cl_i)
    np.add.at(acc["act_i"], b, err_i)

n = np.maximum(acc["n"], 1.0)
print(f"{'d (mm)':>10s} {'pairs':>9s} {'w(d)':>7s} {'claim_c':>9s} {'act_c':>8s} {'c ratio':>8s} "
      f"{'claim_i':>9s} {'act_i':>8s} {'i ratio':>8s} {'mis-wt':>8s}")
for b in range(N_BINS):
    if acc["n"][b] == 0:
        continue
    cc, ac = acc["claim_c"][b] / n[b], acc["act_c"][b] / n[b]
    ci, ai = acc["claim_i"][b] / n[b], acc["act_i"][b] / n[b]
    rc, ri = ac / cc, ai / ci
    print(f"{EDGES[b]:4.0f}-{EDGES[b+1]:<5.0f} {acc['n'][b]:9.0f} {acc['w'][b]/n[b]:7.3f} "
          f"{cc:9.4f} {ac:8.4f} {rc:8.3f} {ci:9.4f} {ai:8.4f} {ri:8.3f} {rc/ri:8.2f}")

tc = acc["act_c"].sum() / acc["claim_c"].sum()
ti = acc["act_i"].sum() / acc["claim_i"].sum()
print(f"\noverall: coordinates {tc:.3f}, images {ti:.3f}, "
      f"coordinates weighted {tc/ti:.3f}x what they deserve")
print("(<1 = under-weighted: they claim more variance than their errors justify)")
