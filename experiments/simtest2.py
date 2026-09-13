"""Does spatial similarity predict effect magnitude, across 258 real reference maps?"""
import json, warnings; warnings.simplefilter("ignore")
import numpy as np
from scipy import stats

X = np.load("/tmp/claude-0/corpus/g_vectors.npy").astype(np.float64)
meta = json.load(open("/tmp/claude-0/corpus/g_meta.json"))
n_maps = len(X)
print(f"{n_maps} maps, {X.shape[1]} voxels\n")

# Magnitude: mean |g| in each map's own top decile. A whole-brain mean would mostly measure how
# much of the brain a contrast lights up, not how strong the effect is where it exists.
absX = np.abs(X)
cut = np.percentile(absX, 90, axis=1, keepdims=True)
magnitude = np.array([a[a >= c].mean() for a, c in zip(absX, cut)])
print(f"magnitude (mean |g| in top decile): median {np.median(magnitude):.3f}, "
      f"5-95 pct {np.percentile(magnitude,5):.3f}-{np.percentile(magnitude,95):.3f}")

# Similarity: spatial correlation between maps.
centred = X - X.mean(axis=1, keepdims=True)
norm = np.linalg.norm(centred, axis=1, keepdims=True)
norm[norm == 0] = np.inf
unit = centred / norm
S = unit @ unit.T
iu = np.triu_indices(n_maps, k=1)
sim = S[iu]
log_mag = np.log(np.clip(magnitude, 1e-6, None))
gap = np.abs(log_mag[iu[0]] - log_mag[iu[1]])       # magnitude dissimilarity, log scale
print(f"pairwise similarity: median {np.median(sim):+.3f}, "
      f"5-95 pct {np.percentile(sim,5):+.3f}-{np.percentile(sim,95):+.3f}")

r, p = stats.pearsonr(sim, gap)
print(f"\ncorr(spatial similarity, |log magnitude difference|) = {r:+.4f}  (p={p:.2g})")
print(f"  variance explained: {r**2:.2%}")
print("  negative r is the hypothesis: more similar maps -> smaller magnitude gap")

# The decisive comparison: does similarity weighting beat the corpus mean, held out?
print(f"\n{'predictor':>34s} {'median |log error|':>19s} {'vs corpus mean':>15s}")
base_err, knn_err, kern_err = [], [], []
for i in range(n_maps):
    others = np.delete(np.arange(n_maps), i)
    base_err.append(abs(log_mag[i] - np.median(log_mag[others])))
    order = others[np.argsort(-S[i, others])]
    knn_err.append(abs(log_mag[i] - np.median(log_mag[order[:15]])))
    w = np.clip(S[i, others], 0, None) ** 4          # kernel weighting on similarity
    kern_err.append(abs(log_mag[i] - (np.sum(w * log_mag[others]) / max(w.sum(), 1e-12))))
for name, err in (("corpus median (no similarity)", base_err),
                  ("15 most similar maps", knn_err),
                  ("similarity-kernel weighted", kern_err)):
    med = np.median(err)
    rel = med / np.median(base_err)
    print(f"{name:>34s} {med:19.4f} {rel:15.2f}")
print("\nratio < 1 means similarity helps; ~1 means it collapses to the corpus average.")

# Is any relationship spatial, or just 'same paradigm'?
para = [m.get("paradigm") for m in meta]
same = np.array([para[a] is not None and para[a] == para[b] for a, b in zip(*iu)])
if same.sum() > 20:
    print(f"\nsame cognitive paradigm ({same.sum()} pairs): median |log gap| "
          f"{np.median(gap[same]):.3f} vs {np.median(gap[~same]):.3f} for different")
