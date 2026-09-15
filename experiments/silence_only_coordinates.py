"""Should a coordinate table contribute only its silence, and not its peak heights?

The stratified real-data result says the coordinate corpus helps a two-image collection through
its **silence** and hurts through its **values**. The bias reduction is 54% at voxels no focus
reaches and 6% where two studies report, shrinking monotonically as more report -- the opposite of
what the magnitude channel would predict. A study that reported nothing near a voxel still enters
the censoring term, and that silence is evidence the effect is small.

So the design question: keep the coordinate tables in the censoring term, drop their magnitudes
from the pooled mean entirely, and see whether that captures the benefit without the cost.

**This is not the flattening idea that was measured and withdrawn.** Flattening replaced every
peak height with a constant in a *coordinates-only* fit, where the magnitudes are the only source
of localisation -- so it took `r` from 0.230 to 0.04 and the top-decile AUC from 0.653 to chance.
Here the images supply the pattern and the coordinates are asked for nothing but presence and
absence. Whether that distinction rescues the idea is exactly what is unmeasured.

The variant is implemented as an experiment-local subclass rather than a change to the estimator,
because it is a design probe and should not be shipped on the strength of a hypothesis. It works
because the two channels are already separate inside the fit: `_accumulate` builds the pooled mean
from the focus table, while `_apply_selection_model` takes the roster from `sample_sizes` and the
silence geometry from the same table. Filtering the table inside `_accumulate` alone therefore
removes the coordinate magnitudes and leaves every study's silence intact.

Four arms, scored against the same held-out half:

  images only          inverse-variance pooling of the two images, what a user would do instead
  mixed                CBES as it ships, documented configuration
  silence only         the variant: coordinate tables in the censoring term, magnitudes dropped
  coordinates only     CBES with no images, for the floor

Reported with signed error beside absolute, because the absolute one can be gamed by shrinkage,
and with pattern metrics beside both, because those are what the flattening retraction turned on.
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
N_IMAGES = int(os.environ.get("NIMAGES", 2))
N_SPLITS = int(os.environ.get("NSPLITS", 8))
WORKDIR = f"/tmp/claude-0/silence_{os.getpid()}"
os.makedirs(WORKDIR, exist_ok=True)

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


class SilenceOnlyCBES(CBES):
    """CBES that takes presence and absence from coordinate tables but not their magnitudes.

    Only ``_accumulate`` is touched, and only to drop the coordinate rows from the pooled mean.
    ``_apply_selection_model`` is called by ``_statistic`` with the unfiltered table, so every
    study's silence geometry and the full roster survive -- which is the whole point.
    """

    def _accumulate(self, table, image_studies=None):
        # Only filter when there are images to carry the mean. `_resolve_peak_bias_scale` calls
        # `_statistic` with `image_studies=None` to fit the coordinates alone, and filtering
        # there leaves nothing at all -- which is how the first version of this subclass raised
        # "No study contributed any in-mask voxels" from inside the calibration sequence.
        #
        # That failure is informative rather than incidental: with the coordinate magnitudes
        # discarded there is nothing for a peak-height scale to correct, so this arm is run with
        # `peak_bias=None` and never reaches the calibration path at all. The guard stays so the
        # subclass cannot silently misbehave if that changes.
        keep = set(image_studies or ())
        if keep and len(table):
            table = table[table["id"].isin(keep)]
        return super()._accumulate(table, image_studies)


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


def fit(cls, studies, calibrate):
    if len(studies) < 2:
        return None
    est = cls(fwhm=10.0, mask=masker, null_method="none",
              peak_bias="per-study" if calibrate else None,
              peak_bias_scale="images" if calibrate else 1.0,
              threshold="reporting_threshold")
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    return (np.abs(res.get_map("g", return_type="array").ravel()),
            res.get_map("n_studies", return_type="array").ravel() > 0)


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    strong = t >= np.percentile(t, 90)
    return (stats.pearsonr(a, t)[0], float(stats.spearmanr(a, t)[0]), float(auc),
            float(np.mean(a - t)), float(np.mean(a[strong] - t[strong])),
            float(np.sqrt(np.mean((a - t) ** 2))))


if __name__ == "__main__":
    ss = load_pain()
    maps, sizes = [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
    maps = np.array(maps)
    total = len(maps)
    print(f"NIDM pain: {total} studies, {N_IMAGES} images held fixed, {N_SPLITS} splits\n")

    rng = np.random.default_rng(0)
    rows = {}
    for split in range(N_SPLITS):
        order = rng.permutation(total)
        half = total // 2
        work, held = order[:half], order[half:]
        truth = np.abs(pooled(held, maps, sizes))
        images, tables = list(work[:N_IMAGES]), list(work[N_IMAGES:])
        studies = build(images, tables, maps, sizes)
        mixed = fit(CBES, studies, True)
        # peak_bias off: there are no coordinate magnitudes left for a scale to act on.
        silence = fit(SilenceOnlyCBES, studies, False)
        coords = fit(CBES, build([], list(work), maps, sizes), False)
        if mixed is None or silence is None or coords is None:
            continue
        only = np.abs(pooled(images, maps, sizes))
        use = mixed[1] & silence[1] & coords[1] & np.isfinite(only) & np.isfinite(truth)
        if use.sum() < 100:
            continue
        for name, est in (("images only", only), ("mixed (ships)", mixed[0]),
                          ("silence only", silence[0]), ("coordinates only", coords[0])):
            rows.setdefault(name, []).append(score(est, truth, use))
        print(f"  split {split + 1}/{N_SPLITS}", flush=True)

    print(f"\n{'estimate':>18} {'r':>7} {'rank r':>7} {'AUC':>6} {'mean err':>9} "
          f"{'err at top':>11} {'rmse':>7}")
    for name in ("images only", "mixed (ships)", "silence only", "coordinates only"):
        if name not in rows:
            continue
        r, rho, auc, err, err_top, rmse = np.nanmean(np.array(rows[name]), axis=0)
        print(f"{name:>18} {r:+7.3f} {rho:+7.3f} {auc:6.3f} {err:+9.3f} {err_top:+11.3f} "
              f"{rmse:7.3f}")
    # Paired tests across splits, because "wins on every metric" read off eight-split means is
    # eyeballing. The splits are paired -- every arm sees the same studies in the same split --
    # so the paired t on the per-split differences is the right test, and it is the difference
    # from the shipping configuration that decides whether this is worth proposing.
    names = ("r", "rank r", "AUC", "mean err", "err at top", "rmse")
    print("\nPaired across splits, silence-only minus each comparator "
          f"(n = {len(rows['silence only'])} splits):")
    for other in ("mixed (ships)", "images only"):
        a = np.array(rows["silence only"], dtype=float)
        b = np.array(rows[other], dtype=float)
        print(f"\n  vs {other}")
        print(f"    {'metric':>11} {'difference':>11} {'paired p':>9}")
        for j, metric in enumerate(names):
            d = a[:, j] - b[:, j]
            ok = np.isfinite(d)
            if ok.sum() < 2:
                continue
            p = float(stats.ttest_rel(a[ok, j], b[ok, j]).pvalue)
            # For the error columns a difference toward zero is the improvement, so say which.
            better = ""
            if metric in ("mean err", "err at top"):
                better = " (closer to 0)" if abs(np.mean(a[ok, j])) < abs(np.mean(b[ok, j])) else ""
            elif metric == "rmse":
                better = " (lower)" if np.mean(d) < 0 else ""
            else:
                better = " (higher)" if np.mean(d) > 0 else ""
            print(f"    {metric:>11} {np.mean(d):+11.3f} {p:9.4f}{better}")
    print("\nThe variant earns its place only if it holds pattern near images-only AND beats it")
    print("on the error at the top decile. Pattern collapsing toward the coordinates-only floor")
    print("would mean the silence geometry distorts the map, which the flattening retraction")
    print("would have predicted and this design was built to escape.")
