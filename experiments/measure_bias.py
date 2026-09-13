"""What does a reported peak actually tell us about the true effect nearby?

For each reported peak of study k, compare the value study k reports against a
leave-one-out pooled estimate of the true effect at that location and at increasing
distances from it. Leave-one-out matters: pooling that included study k would inherit
the very selection we are trying to measure.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy.spatial import cKDTree
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import _mask_img_to_bool, vox2mm

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
mask_img = masker.mask_img
world = vox2mm(np.vstack(np.where(_mask_img_to_bool(mask_img))).T, mask_img.affine)
tree = cKDTree(world)
ids = list(ss.ids)
sizes = np.asarray(ss.sample_sizes(), dtype=float)

G, V = [], []
for gp, vp in zip(ss.images["g"], ss.images["g_var"]):
    G.append(masker.transform(str(gp)).ravel())
    V.append(masker.transform(str(vp)).ravel())
G, V = np.vstack(G), np.vstack(V)
W = 1.0 / np.clip(V, 1e-12, None)

cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
coords = cs.coordinates

BINS = [0, 2, 4, 6, 8, 10, 12, 16, 20]
own_at_peak, loo = {b: [] for b in BINS}, {b: [] for b in BINS}
peak_n = []
for k, sid in enumerate(ids):
    sub = coords[coords["id"] == sid]
    if not len(sub):
        continue
    others = [j for j in range(len(ids)) if j != k]
    wl = W[others]
    pooled = (wl * G[others]).sum(0) / wl.sum(0)      # leave-one-out pooled effect
    for point in sub[["x", "y", "z"]].to_numpy(float):
        peak_idx = tree.query(point)[1]
        own = abs(G[k][peak_idx])
        peak_n.append(sizes[k])
        for b in BINS:
            shell = tree.query_ball_point(point, r=b + 1.0)
            if b > 0:
                inner = set(tree.query_ball_point(point, r=b - 1.0))
                shell = [i for i in shell if i not in inner]
            if len(shell) < 3:
                continue
            own_at_peak[b].append(own)
            loo[b].append(np.abs(pooled[shell]).mean())

print(f"{len(peak_n)} reported peaks, mean N = {np.mean(peak_n):.1f}\n")
print(f"{'distance (mm)':>14s} {'reported |g|':>13s} {'true |g| (LOO)':>15s} {'ratio':>7s} {'n':>6s}")
for b in BINS:
    if len(loo[b]) < 20:
        continue
    a, t = np.array(own_at_peak[b]), np.array(loo[b])
    print(f"{b:14d} {a.mean():13.3f} {t.mean():15.3f} {t.mean()/a.mean():7.3f} {len(t):6d}")
