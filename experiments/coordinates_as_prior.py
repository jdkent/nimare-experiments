"""Use coordinates for *where*, images for *how much*, and see whether that beats images alone.

Pooling coordinates and images at the same voxel is measurably worse than the images by
themselves, because a reported peak value is a biased observation of the quantity the image
measures without bias, and a shared bias does not average away. But the coordinate *locations*
are informative: a count of studies reporting nearby correlates +0.54 to +0.60 with held-out
truth under realistic reporting, better than the magnitude map manages.

That suggests a division of labour rather than a pooling. The image estimate carries the scale
and is unbiased but noisy when only a few studies shared maps. The coordinate density says where
the literature has repeatedly found something. Shrinking the image estimate toward zero where
nothing was ever reported, and leaving it alone where many studies were, uses each for what it is
good at and never averages a biased value into an unbiased one.

Four estimators, all scored against studies neither arm saw:

  * **images only** -- inverse-variance pooling of the `k` image studies;
  * **shrunk by coordinate density** -- the same, multiplied by `c / (c + half)` where `c` counts
    studies reporting within 10 mm, a soft gate with one parameter;
  * **hard coordinate gate** -- the same, zeroed where no study reported within 10 mm;
  * **ORACLE gate** -- zeroed where the held-out truth is below its median, which no real
    procedure can do. It bounds what any gating of the image estimate could achieve, so a
    shrinkage that lands far short of it has room, and one that matches it is done.

If the shrunk row beats images alone, coordinates have a job in a collection that has images.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
RADIUS_MM = 10.0
HALF = 2.0            # studies at which the soft gate reaches one half
N_SPLITS = 10
N_IMAGES = (1, 2, 3, 5)

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape = mask_img.shape
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)

_r = int(np.ceil(RADIUS_MM / zooms.min()))
_grid = np.stack(np.meshgrid(*[np.arange(-_r, _r + 1)] * 3, indexing="ij"), -1)
BALL = _grid[np.linalg.norm(_grid * zooms, axis=-1) <= RADIUS_MM]


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g = to_g(maps[i], sizes[i])
        var = 1.0 / sizes[i] + g**2 / (2.0 * sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(var, 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def study_count(members, maps, sizes):
    """How many of these studies reported a focus within the radius of each voxel."""
    total = np.zeros(shape, dtype=np.int32)
    for i in members:
        foci, _ = report_peaks(maps[i], mask_bool, shape, zooms, scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        seen = np.zeros(shape, dtype=bool)
        for ijk, _ in foci:
            pts = np.asarray(ijk) + BALL
            ok = np.all((pts >= 0) & (pts < np.array(shape)), axis=1)
            pts = pts[ok]
            seen[(pts[:, 0], pts[:, 1], pts[:, 2])] = True
        total += seen
    return total[mask_bool].astype(float)


def score(est, truth, use):
    e, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(e)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
           if n_pos and n_neg else np.nan)
    return (stats.pearsonr(e, t)[0], float(stats.spearmanr(e, t)[0]), float(auc),
            float(np.sqrt(np.mean((e - t) ** 2))))


ss = ImageTransformer(target="z").transform(load_pain())
maps, sizes = [], []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    path = getattr(row, "z", None)
    if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
        continue
    img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                          force_resample=True, copy_header=True)
    maps.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
    sizes.append(float(n))
sizes = np.asarray(sizes, dtype=float)
n = len(maps)
print(f"NIDM pain: {n} studies, {SCHEME}/{FOCUS} coordinates, {RADIUS_MM:.0f} mm radius, "
      f"{N_SPLITS} splits\n", flush=True)

acc = {}
rng = np.random.default_rng(0)
for _ in range(N_SPLITS):
    order = rng.permutation(n)
    work, hold = order[: n // 2], order[n // 2:]
    truth = np.abs(pooled(hold, maps, sizes))
    for k in N_IMAGES:
        if k >= len(work):
            continue
        images, tables = list(work[:k]), list(work[k:])
        base = np.abs(pooled(images, maps, sizes))
        count = study_count(tables, maps, sizes)
        use = np.isfinite(base) & np.isfinite(truth)
        if use.sum() < 100:
            continue
        variants = {
            "images only": base,
            "shrunk by coordinate density": base * (count / (count + HALF)),
            "hard coordinate gate": base * (count > 0),
            "ORACLE gate on the truth": base * (truth >= np.median(truth[use])),
        }
        for name, est in variants.items():
            acc.setdefault((k, name), []).append(score(est, truth, use))

print(f"  {'images':>7} {'estimator':>30} {'r':>7} {'rank r':>8} {'AUC':>7} {'rmse':>7}")
for k in N_IMAGES:
    for name in ("images only", "shrunk by coordinate density", "hard coordinate gate",
                 "ORACLE gate on the truth"):
        rows = acc.get((k, name))
        if not rows:
            continue
        r, rho, auc, rmse = np.nanmean(np.array(rows), axis=0)
        print(f"  {k:>7} {name:>30} {r:+7.3f} {rho:+8.3f} {auc:7.3f} {rmse:7.3f}")
    print()
print("The oracle row is not attainable; it bounds what any gate on the image estimate could do.")
