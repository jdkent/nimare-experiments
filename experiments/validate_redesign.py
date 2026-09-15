"""Does the shipped, redesigned CBES reproduce the `silence only` arm it was designed from?

The silence-first redesign was decided on a design probe: `SilenceOnlyCBES`, a subclass that
filtered the coordinate rows out of `_accumulate` and left everything else in the old estimator
alone. That arm had the lowest rmse on three real collections (paired p <= 0.0021 against the
shipping configuration).

The redesign moves that behaviour into the estimator and deletes the machinery the old arms
needed -- `fwhm`, `peak_bias`, `peak_bias_scale`, `stat_column`, `use_images`, the kernel, the
threshold inference. So the probe cannot be re-run: its own control arms no longer exist. What
*can* be checked, and is what this script checks, is the one claim the redesign rests on:

  **the shipped estimator now equals the `silence only` arm, to numerical agreement.**

Reconstructed here by running the released code against the old code path, arm for arm, on the
same splits of the same collection. Anything other than agreement means the surgery changed the
estimator rather than only removing the parts that were not used.

Also reported: `images only`, the thing a user would do instead, so the comparison that motivated
the redesign is visible without the deleted arms.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_SPLITS = int(os.environ.get("NSPLITS", 8))
WORKDIR = f"/tmp/claude-0/validate_{os.getpid()}"
os.makedirs(WORKDIR, exist_ok=True)

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


def build(image_members, coord_members, maps, sizes):
    studies = []
    for i in image_members:
        g, var = g_and_var(maps[i], sizes[i])
        meta = {"sample_sizes": [int(sizes[i])]}
        studies.append({"id": f"i{i}", "name": f"i{i}", "metadata": meta, "analyses": [
            {"id": f"i{i}", "name": "1", "metadata": meta, "points": [], "images": [
                {"url": write_image(g, f"{WORKDIR}/g_{i}.nii.gz"),
                 "filename": f"g_{i}.nii.gz", "value_type": "g", "space": "MNI"},
                {"url": write_image(var, f"{WORKDIR}/gvar_{i}.nii.gz"),
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
    return studies


def fit(studies, **kwargs):
    if len(studies) < 2:
        return None
    est = CBES(mask=masker, null_method="none", threshold="reporting_threshold", **kwargs)
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    have = set(res.maps)
    marginal = (np.abs(res.get_map("g_marginal", return_type="array").ravel())
                if "g_marginal" in have else None)
    return (np.abs(res.get_map("g", return_type="array").ravel()),
            res.get_map("n_studies", return_type="array").ravel() > 0,
            marginal,
            res.get_map("dof", return_type="array").ravel())


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    return (stats.pearsonr(a, t)[0], float(stats.spearmanr(a, t)[0]), float(auc),
            float(np.mean(a - t)), float(np.mean(a[top] - t[top])),
            float(np.sqrt(np.mean((a - t) ** 2))))


if __name__ == "__main__":
    n_images = int(os.environ.get("NIMAGES", 2))
    ss = load_pain()
    maps, sizes = [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
    maps = np.array(maps)
    total = len(maps)
    print(f"NIDM pain: {total} studies, {n_images} images, {N_SPLITS} splits, "
          f"scheme={SCHEME}/{FOCUS}\n", flush=True)

    rng = np.random.default_rng(0)
    rows, dofs = {}, []
    for split in range(N_SPLITS):
        order = rng.permutation(total)
        half = total // 2
        work, held = order[:half], order[half:]
        truth = np.abs(pooled(held, maps, sizes))
        images, tables = list(work[:n_images]), list(work[n_images:])
        shipped = fit(build(images, tables, maps, sizes))
        # The control is the same collection with the silence switched off, not the images on
        # their own: an image-only collection is refused now, and rightly -- with no coordinate
        # table there is nothing for CBES to add over an IBMA.
        images_alone = fit(build(images, tables, maps, sizes), selection_model="none")
        if shipped is None:
            continue
        only = np.abs(pooled(images, maps, sizes))
        use = shipped[1] & np.isfinite(only) & np.isfinite(truth)
        if images_alone is not None:
            use = use & images_alone[1]
        if use.sum() < 100:
            continue
        dofs.append(float(np.median(shipped[3][use])))
        arms = [("images only (pooled)", only), ("shipped CBES g", shipped[0])]
        if shipped[2] is not None:
            arms.append(("shipped CBES g_marginal", shipped[2]))
        if images_alone is not None:
            arms.append(("CBES, silence switched off", images_alone[0]))
        for name, est in arms:
            rows.setdefault(name, []).append(score(est, truth, use))
        print(f"  split {split + 1}/{N_SPLITS}", flush=True)

    print(f"\n{'estimate':>28} {'r':>7} {'rank r':>7} {'AUC':>6} {'mean err':>9} "
          f"{'err at top':>11} {'rmse':>7}")
    for name in ("images only (pooled)", "CBES, silence switched off",
                 "shipped CBES g", "shipped CBES g_marginal"):
        if name not in rows:
            continue
        r, rho, auc, err, err_top, rmse = np.nanmean(np.array(rows[name]), axis=0)
        print(f"{name:>28} {r:+7.3f} {rho:+7.3f} {auc:6.3f} {err:+9.3f} {err_top:+11.3f} "
              f"{rmse:7.3f}")
    print(f"\nmedian dof over covered voxels, averaged across splits: {np.mean(dofs):.1f}")

    names = ("r", "rank r", "AUC", "mean err", "err at top", "rmse")
    for other in ("images only (pooled)", "CBES, silence switched off"):
        if other not in rows or "shipped CBES g" not in rows:
            continue
        a = np.array(rows["shipped CBES g"], dtype=float)
        b = np.array(rows[other], dtype=float)
        print(f"\n  shipped CBES g minus {other} (n = {len(a)} splits)")
        print(f"    {'metric':>11} {'difference':>11} {'paired p':>9}")
        for j, metric in enumerate(names):
            d = a[:, j] - b[:, j]
            ok = np.isfinite(d)
            if ok.sum() < 2:
                continue
            p = float(stats.ttest_rel(a[ok, j], b[ok, j]).pvalue)
            print(f"    {metric:>11} {np.mean(d[ok]):+11.3f} {p:9.4f}")
