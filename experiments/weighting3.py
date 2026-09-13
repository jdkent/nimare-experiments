"""What is a coordinate study's actual error as a function of distance from its reported peak?

The estimator gives an image study weight 1 at every voxel and a coordinate study the kernel
weight w(d) at distance d, with precision w/var in both cases. A coordinate contribution
therefore has effective variance var/w(d) against an image's var. Whether coordinates are
correctly weighted against images reduces to one question: is 1/w(d) the right variance
inflation?

That is measurable. Both kinds of study are scored on *the same voxels* -- every voxel within
kernel support of a reported peak -- against the same leave-one-out truth. The earlier version
of this diagnostic scored images over the whole brain and coordinates only near peaks, so its
43x was partly a comparison of two different voxel populations.

    claimed_c(d) = (var_k + tau^2) / w(d)      what the model says a coordinate knows at d
    actual_c(d)  = (rho g_peak - truth(v))^2   what it actually gets wrong there
    claimed_i    = var_k(v) + tau^2            the same for the image, at the same voxels
    actual_i     = (G_k(v) - truth(v))^2

If 1/w(d) is right, actual_c/claimed_c is flat in d and equal to actual_i/claimed_i. The shape
of actual_c(d) is the deliverable: it *is* the correct variance model, whatever the kernel says.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.meta.utils import get_ale_kernel
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import mm2vox

U, FWHM, MIN_W = 3.2905, 10.0, 0.01
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
G = np.array([g_maps[k] for k in keys])
V = np.array([v_maps[k] for k in keys])
W = 1.0 / V
pooled = (G * W).sum(0) / np.clip(W.sum(0), 1e-12, None)

# tau^2 by DerSimonian-Laird at each voxel, from the images. The same for both kinds of study.
q = (W * (G - pooled) ** 2).sum(0)
den = W.sum(0) - (W**2).sum(0) / np.clip(W.sum(0), 1e-12, None)
tau2 = np.clip((q - (len(keys) - 1)) / np.clip(den, 1e-12, None), 0.0, None)

index = {k: i for i, k in enumerate(keys)}
def loo(sid):
    keep = np.arange(len(keys)) != index[sid]
    w = W[keep]
    return (G[keep] * w).sum(0) / np.clip(w.sum(0), 1e-12, None)
truths = {sid: loo(sid) for sid in keys}

coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates
coords["id"] = coords["id"].astype(str)

# Kernel support, with each offset's distance in mm, so contributions can be binned by distance.
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

# Pass 1: rho, from the peak voxels alone -- how inflated a reported peak is, with no distance
# effect mixed in.
peak_est, peak_truth = [], []
foci = []
for sid, sub in coords.groupby("id"):
    if sid not in g_maps:
        continue
    n = sizes[sid]
    z = sub["z_stat"].astype(float).to_numpy()
    g_k, var_k = peak_stat_to_hedges_g(np.abs(z), np.full(len(z), n), stat_type="z")
    g_k = np.sign(z) * g_k
    ijk = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    truth = truths[sid]
    for centre, g_focus, var_focus in zip(ijk, g_k, var_k):
        flat = np.ravel_multi_index(np.clip(centre, 0, np.array(shape) - 1), shape)
        col = lookup[flat]
        if col < 0:
            continue
        foci.append((sid, centre, g_focus, var_focus))
        peak_est.append(abs(g_focus))
        peak_truth.append(abs(truth[col]))

rho = float(np.mean(peak_truth) / np.mean(peak_est))
print(f"{len(foci)} foci; peak inflation rho = {rho:.3f} "
      f"(mean |g| reported {np.mean(peak_est):.3f} vs mean |truth| at the peak "
      f"{np.mean(peak_truth):.3f})\n", flush=True)

# Pass 2: accumulate claimed and actual, by distance bin, for both kinds on the same voxels.
acc = {name: np.zeros(N_BINS) for name in
       ("n", "w", "claim_c", "act_c", "claim_i", "act_i")}
for sid, centre, g_focus, var_focus in foci:
    truth = truths[sid]
    gk, vk = G[index[sid]], V[index[sid]]
    vox = centre + offs
    inside = np.all((vox >= 0) & (vox < shape), axis=1)
    cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
    keep = cols >= 0
    cols = cols[keep]
    if not cols.size:
        continue
    w = wts[inside][keep]
    b = which[inside][keep]
    finite = np.isfinite(vk[cols]) & (vk[cols] < np.inf)
    cols, w, b = cols[finite], w[finite], b[finite]
    if not cols.size:
        continue

    t, t2 = truth[cols], tau2[cols]
    np.add.at(acc["n"], b, 1.0)
    np.add.at(acc["w"], b, w)
    np.add.at(acc["claim_c"], b, (var_focus + t2) / w)
    np.add.at(acc["act_c"], b, (rho * g_focus - t) ** 2)
    np.add.at(acc["claim_i"], b, vk[cols] + t2)
    np.add.at(acc["act_i"], b, (gk[cols] - t) ** 2)

n = np.maximum(acc["n"], 1.0)
print(f"{'d (mm)':>10s} {'pairs':>9s} {'w(d)':>7s} {'claim_c':>9s} {'act_c':>8s} "
      f"{'c ratio':>8s} {'claim_i':>9s} {'act_i':>8s} {'i ratio':>8s} {'mis-wt':>8s}")
for b in range(N_BINS):
    if acc["n"][b] == 0:
        continue
    cc, ac = acc["claim_c"][b] / n[b], acc["act_c"][b] / n[b]
    ci, ai = acc["claim_i"][b] / n[b], acc["act_i"][b] / n[b]
    rc, ri = ac / cc, ai / ci
    print(f"{EDGES[b]:4.0f}-{EDGES[b+1]:<5.0f} {acc['n'][b]:9.0f} {acc['w'][b]/n[b]:7.3f} "
          f"{cc:9.4f} {ac:8.4f} {rc:8.3f} {ci:9.4f} {ai:8.4f} {ri:8.3f} {rc/ri:8.2f}")

tot_c = acc["act_c"].sum() / acc["claim_c"].sum()
tot_i = acc["act_i"].sum() / acc["claim_i"].sum()
print(f"\noverall: coordinates {tot_c:.3f}, images {tot_i:.3f}, "
      f"coordinates over-weighted by {tot_c/tot_i:.2f}x")
print("(<1 means coordinates claim more variance than their errors justify: under-weighted)")
