"""Given a few images, do the coordinate-only studies add anything or cost something?

The magnitude from coordinates is compressed and its level depends on the reporting convention.
Images are unbiased and pin the scale. A real collection has a handful of images and many more
studies that only ever published a table, so the design question is not which is better but
whether the tables are worth including at all.

Three estimates are scored against the same held-out truth:

  * **images only** -- inverse-variance pooling of the ``k`` image studies, which is what a user
    would do instead of this estimator;
  * **images plus coordinates** -- CBES with those ``k`` as images and the remaining studies of
    the same half entering as thresholded coordinate tables;
  * **coordinates only** -- CBES on the whole half as tables, with no images at all.

If the middle beats the first, the coordinate pathway earns its place and the estimator should
keep it. If it does not, the compressed tables are diluting good data and the estimator should
say so, or refuse them.

The mixed fit is run under two settings, because the first version of this script got it wrong.
With ``peak_bias=None`` the coordinates enter on their inflated scale -- a mean g near 2.0 on
this collection -- while the images enter on the true one near 0.5, and nothing reconciles them.
The estimator ships ``peak_bias="per-study"`` with ``peak_bias_scale="images"`` for exactly this
case, and its own docstring warns that a mismatched constant makes the two kinds of study
disagree about the same voxel. Judging mixing without it tests a configuration the code tells
you not to use, and would explain a mix scoring worse than either of its parts.

Truth is the inverse-variance pooling of the *other* half, so no study contributes to both
sides. Coordinates come from :mod:`reporting`: multiplicity-corrected, whole surviving clusters,
one focus each, nothing capped.
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
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_IMAGES = (1, 2, 3, 5)
N_SPLITS = 8

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def g_and_var(z, n):
    g = to_g(z, n)
    return g, 1.0 / n + g**2 / (2.0 * n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g, var = g_and_var(maps[i], sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(var, 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def write_image(values, path):
    volume = np.zeros(shape, dtype=np.float32)
    volume[mask_bool] = values
    nib.save(nib.Nifti1Image(volume, affine), path)
    return path


def cbes_fit(image_members, coord_members, maps, sizes, workdir, calibrate=False):
    studies = []
    for i in image_members:
        g, var = g_and_var(maps[i], sizes[i])
        meta = {"sample_sizes": [int(sizes[i])]}
        studies.append({"id": f"i{i}", "name": f"i{i}", "metadata": meta, "analyses": [
            {"id": f"i{i}", "name": "1", "metadata": meta, "points": [], "images": [
                {"url": write_image(g, f"{workdir}/g_{i}.nii.gz"),
                 "filename": f"g_{i}.nii.gz", "value_type": "g", "space": "MNI"},
                {"url": write_image(var, f"{workdir}/gvar_{i}.nii.gz"),
                 "filename": f"gvar_{i}.nii.gz", "value_type": "g_var", "space": "MNI"}]}]})
    for i in coord_members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                    scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(sizes[i])], "reporting_threshold": float(height)}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": v}]} for ijk, v in foci]}]})
    if len(studies) < 2:
        return None
    calibrated = calibrate and image_members and coord_members
    est = CBES(
        fwhm=10.0,
        mask=masker,
        peak_bias="per-study" if calibrated else None,
        peak_bias_scale="images" if calibrated else 1.0,
        null_method="none",
        threshold="reporting_threshold" if coord_members else "study-min",
    )
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    g = res.get_map("g", return_type="array").ravel()
    covered = res.get_map("n_studies", return_type="array").ravel() > 0
    return np.abs(g), covered


def score(estimate, use, truth):
    """Magnitude accuracy and localisation, which are different questions.

    The magnitude columns ask whether the number is right. The localisation columns ask only
    whether the map ranks voxels correctly and finds the strong ones -- which is what a reader
    uses a meta-analytic map for, and which a biased but monotone estimate can do well. Scoring
    only the first is what made the earlier version of this experiment answer half the question.
    """
    est, tru = estimate[use], truth[use]
    top = tru >= np.percentile(tru, 90)
    # AUC for "is this voxel in the truth's top decile", from the rank-sum identity.
    order = stats.rankdata(est)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
           if n_pos and n_neg else np.nan)
    return (stats.pearsonr(est, tru)[0],
            est.mean() / max(tru.mean(), 1e-9),
            float(np.sqrt(np.mean((est - tru) ** 2))),
            float(stats.spearmanr(est, tru)[0]),
            float(auc))


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
workdir = "/tmp/claude-0/coordhelp"
os.makedirs(workdir, exist_ok=True)
print(f"NIDM pain: {n} studies, {SCHEME}/{FOCUS} reporting, "
      f"{n // 2} per half over {N_SPLITS} splits\n", flush=True)

rows = {}
rng = np.random.default_rng(0)
for _ in range(N_SPLITS):
    order = rng.permutation(n)
    work, hold = order[: n // 2], order[n // 2:]
    truth = np.abs(pooled(hold, maps, sizes))
    for k in N_IMAGES:
        if k >= len(work):
            continue
        images, tables = list(work[:k]), list(work[k:])
        got_both = cbes_fit(images, tables, maps, sizes, workdir)
        got_cal = cbes_fit(images, tables, maps, sizes, workdir, calibrate=True)
        got_coord = cbes_fit([], list(work), maps, sizes, workdir)
        if got_both is None or got_coord is None or got_cal is None:
            continue
        # Images-only reference, scored on the voxels the mixed fit covers so the three
        # estimates are judged on the same set.
        only = np.abs(pooled(images, maps, sizes))
        use = (got_both[1] & got_coord[1] & got_cal[1]
               & np.isfinite(got_both[0]) & np.isfinite(got_cal[0]) & np.isfinite(only))
        if use.sum() < 100:
            continue
        rows.setdefault(k, {"images only": [], "mixed, uncalibrated": [],
                            "mixed, scaled to images": [], "coordinates only": []})
        rows[k]["images only"].append(score(only, use, truth))
        rows[k]["mixed, uncalibrated"].append(score(got_both[0], use, truth))
        rows[k]["mixed, scaled to images"].append(score(got_cal[0], use, truth))
        rows[k]["coordinates only"].append(score(got_coord[0], use, truth))

print(f"  {'images':>7} {'estimate':>24} {'r':>7} {'ratio':>7} {'rmse':>7} "
      f"{'rank r':>7} {'AUC top':>8}")
for k in N_IMAGES:
    if k not in rows:
        continue
    for name, values in rows[k].items():
        r, ratio, rmse, rho, auc = np.mean(np.array(values), axis=0)
        print(f"  {k:>7} {name:>24} {r:+7.3f} {ratio:7.2f} {rmse:7.3f} "
              f"{rho:+7.3f} {auc:8.3f}")
    print()
print("Magnitude accuracy is r, ratio and rmse; localisation is the rank correlation and the"
      "\narea under the curve for recovering the truth's top decile. Coordinates can lose the"
      "\nfirst and win the second, and only the first has been measured until now.")
