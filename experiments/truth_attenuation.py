"""Does the image-pooled truth understate the per-study effect, and by how much?

CBES scored against an image-pooled truth on the pain collection came out 4x high, while the
corrected simulator says it runs slightly *low*. One candidate reconciliation does not involve
the estimator at all: the reference. Averaging unthresholded maps across studies whose
activations sit at slightly different coordinates blurs a focal effect and shrinks its peak,
so the pooled map understates what any single study measured at its own peak.

Measured without CBES. For each study, take the magnitude of its own g at its own strongest
locations, and compare with the pooled truth at those same voxels. A pooled reference that is
a faithful estimate of the common effect would sit near the average of the per-study values;
one that is attenuated by misalignment will sit well below it.

The dispersion of the studies' peak locations is reported alongside, since that is what drives
the attenuation and is what would let the effect be predicted on another collection.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from load_pain import load_pain
from nimare.transforms import ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

g_maps, w_maps = [], []
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g = masker.transform(str(row.g)).ravel()
    v = masker.transform(str(row.g_var)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g_maps.append(np.where(ok, g, 0.0))
    w_maps.append(np.where(ok, 1.0 / np.maximum(v, 1e-12), 0.0))
G, W = np.array(g_maps), np.array(w_maps)
pooled = np.sum(G * W, axis=0) / np.maximum(np.sum(W, axis=0), 1e-12)
print(f"{len(G)} studies\n")

# Each study's own strongest voxels, and what the pooled map says there.
print(f"{'top %':>6} {'per-study |g| at its own top':>29} {'pooled |g| there':>17} {'ratio':>7}")
for pct in (99.9, 99.5, 99.0, 95.0):
    own, there = [], []
    for g in G:
        cut = np.percentile(np.abs(g), pct)
        hot = np.abs(g) >= cut
        own.append(np.abs(g[hot]).mean())
        there.append(np.abs(pooled[hot]).mean())
    own, there = np.mean(own), np.mean(there)
    print(f"{pct:6.1f} {own:29.3f} {there:17.3f} {own / there:7.2f}")

# How far apart are the studies' strongest locations? That is what drives the attenuation.
coords = []
for g in G:
    coords.append(int(np.argmax(np.abs(g))))
xyz = masker.inverse_transform(np.eye(1, G.shape[1], 0)).affine
import nibabel as nib
mask_idx = np.argwhere(np.asarray(masker.mask_img_.get_fdata()) > 0)
points = np.array([nib.affines.apply_affine(masker.mask_img_.affine, mask_idx[c].astype(float))
                   for c in coords])
centroid = points.mean(axis=0)
spread = np.sqrt(((points - centroid) ** 2).sum(axis=1))
print(f"\nglobal-maximum locations: median distance from their centroid "
      f"{np.median(spread):.1f} mm, range {spread.min():.1f}-{spread.max():.1f} mm")
print(f"pooled map's own max |g|: {np.abs(pooled).max():.3f}; "
      f"mean of per-study max |g|: {np.mean([np.abs(g).max() for g in G]):.3f}; "
      f"ratio {np.mean([np.abs(g).max() for g in G]) / np.abs(pooled).max():.2f}")
