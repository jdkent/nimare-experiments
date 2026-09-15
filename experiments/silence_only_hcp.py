"""Does silence-only replicate on HCP, where the truth comes from held-out *subjects*?

The silence-only result -- keep coordinate tables in the censoring term, drop their magnitudes
from the pooled mean -- was measured on one collection, the 21 NIDM pain studies, with the truth
from a held-out half of the same studies. Two things make that a weak base. Pain is a strong,
spatially consistent, high-prevalence effect, the friendliest case there is. And a held-out half
of the same literature shares its conventions, its populations and its smoothness.

HCP collection 4337 has per-subject maps, which fixes the second problem properly: synthetic
studies are cut from one set of subjects and the truth is measured on a **disjoint set of
subjects**, so the selection that produces the peaks is statistically independent of the quantity
compared against. And the construct is different -- a motor contrast rather than pain.

Design, following `hcp_heldout_validation.py`:

  * subjects split into a used part and a held part, disjoint;
  * the used part cut into ``n_studies`` synthetic studies of ``n_per_study`` subjects;
  * each study's own one-sample t mapped to z on ``n - 1`` df -- the reporting model CBES
    inverts -- then `reporting.report_peaks` for the coordinates it would tabulate, corrected
    for multiplicity, whole surviving clusters, one focus each, **nothing capped**;
  * each study's Hedges' g and its variance for the image arm, so an image study and a
    coordinate study describe the same underlying data on one convention;
  * the truth is Hedges' g over the held-out subjects.

The four arms match the pain run: images only (inverse-variance pooling of the ``k`` images),
mixed as CBES ships, silence only, and coordinates only for the floor.

``ratio`` is reported alongside the error columns because the open question this bears on is
whether ``g`` can be read as an *absolute* magnitude with a single image. A ratio near 1 is what
that claim needs; the error columns cannot answer it because they are not scale-free.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z
from reporting import report_peaks

CONTRAST = os.environ.get("CONTRAST", "MOTOR_LH")
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"
N_IMAGES = int(os.environ.get("NIMAGES", 1))
N_PER_STUDY = int(os.environ.get("NPER", 20))
N_STUDIES = int(os.environ.get("NSTUDIES", 12))
N_SPLITS = int(os.environ.get("NSPLITS", 8))
SCHEME, FOCUS = "cluster", "max"
WORKDIR = f"/tmp/claude-0/soh_{os.getpid()}"
os.makedirs(WORKDIR, exist_ok=True)

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
ZOOMS = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


class SilenceOnlyCBES(CBES):
    """Presence and absence from the coordinate tables, but not their magnitudes.

    Only ``_accumulate`` is touched. ``_apply_selection_model`` is called by ``_statistic`` with
    the unfiltered table, so every study's silence geometry and the full roster survive. The
    guard on an empty ``image_studies`` is there because ``_resolve_peak_bias_scale`` fits the
    coordinates alone, which this arm never needs -- with the magnitudes gone there is no scale.
    """

    def _accumulate(self, table, image_studies=None):
        keep = set(image_studies or ())
        if keep and len(table):
            table = table[table["id"].isin(keep)]
        return super()._accumulate(table, image_studies)


def hedges(mean, sd, n):
    correction = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
    return correction * mean / np.maximum(sd, 1e-9)


def write_image(values, path):
    volume = np.zeros(shape, dtype=np.float32)
    volume[mask_bool] = values
    nib.save(nib.Nifti1Image(volume, affine), path)
    return path


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    return (stats.pearsonr(a, t)[0], float(stats.spearmanr(a, t)[0]), float(auc),
            float(a.mean() / max(t.mean(), 1e-9)),
            float(np.mean(a - t)), float(np.mean(a[top] - t[top])),
            float(np.sqrt(np.mean((a - t) ** 2))))


if __name__ == "__main__":
    S = np.load(CACHE)
    print(f"{CONTRAST}: {S.shape[0]} subjects, {S.shape[1]} voxels at 4 mm")
    print(f"{N_STUDIES} synthetic studies of {N_PER_STUDY} subjects, {N_IMAGES} as images, "
          f"{N_SPLITS} splits; truth from the held-out subjects\n")
    rng = np.random.default_rng(0)
    rows, foci_seen = {}, []
    for split in range(N_SPLITS):
        need = N_PER_STUDY * N_STUDIES
        order = rng.permutation(S.shape[0])
        used, held = order[:need], order[need:]
        truth = np.abs(hedges(S[held].mean(axis=0), S[held].std(axis=0, ddof=1), len(held)))

        studies, image_ids, n_foci = [], [], 0
        for k in range(N_STUDIES):
            block = S[used[k * N_PER_STUDY:(k + 1) * N_PER_STUDY]]
            mean, sd = block.mean(axis=0), block.std(axis=0, ddof=1)
            t = mean / np.maximum(sd / np.sqrt(N_PER_STUDY), 1e-9)
            z = np.nan_to_num(t_to_z(t, dof=N_PER_STUDY - 1))
            g = hedges(mean, sd, N_PER_STUDY)
            var = 1.0 / N_PER_STUDY + g**2 / (2.0 * N_PER_STUDY)
            meta = {"sample_sizes": [N_PER_STUDY]}
            if k < N_IMAGES:
                image_ids.append(f"s{k}")
                studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
                    {"id": f"s{k}", "name": "1", "metadata": meta, "points": [], "images": [
                        {"url": write_image(g, f"{WORKDIR}/g_{split}_{k}.nii.gz"),
                         "filename": "g.nii.gz", "value_type": "g", "space": "MNI"},
                        {"url": write_image(var, f"{WORKDIR}/v_{split}_{k}.nii.gz"),
                         "filename": "v.nii.gz", "value_type": "g_var", "space": "MNI"}]}]})
                continue
            found, height = report_peaks(z, mask_bool, shape, ZOOMS,
                                         scheme=SCHEME, focus=FOCUS)
            if not found:
                continue
            n_foci += len(found)
            meta = dict(meta, reporting_threshold=float(height))
            studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
                {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                    {"space": "MNI",
                     "coordinates": [float(c) for c in nib.affines.apply_affine(
                         affine, np.asarray(ijk, dtype=float))],
                     "values": [{"kind": "Z", "value": v}]} for ijk, v in found]}]})
        coord_count = len(studies) - len(image_ids)
        foci_seen.append((coord_count, n_foci))
        if coord_count < 2:
            print(f"  split {split + 1}: only {coord_count} coordinate studies reported, skipped")
            continue

        def run(cls, calibrate):
            est = cls(fwhm=10.0, mask=masker, null_method="none",
                      peak_bias="per-study" if calibrate else None,
                      peak_bias_scale="images" if calibrate else 1.0,
                      threshold="reporting_threshold")
            res = est.fit(Studyset({"id": "h", "name": "h", "studies": studies},
                                   target=None, mask=mask_img))
            return (np.abs(res.get_map("g", return_type="array").ravel()),
                    res.get_map("n_studies", return_type="array").ravel() > 0)

        mixed = run(CBES, True)
        silence = run(SilenceOnlyCBES, False)
        # Images-only reference: inverse-variance pooling of the image studies, which is what a
        # user would do instead of this estimator.
        gs, ws = [], []
        for k in range(N_IMAGES):
            block = S[used[k * N_PER_STUDY:(k + 1) * N_PER_STUDY]]
            g = hedges(block.mean(axis=0), block.std(axis=0, ddof=1), N_PER_STUDY)
            var = 1.0 / N_PER_STUDY + g**2 / (2.0 * N_PER_STUDY)
            gs.append(g)
            ws.append(1.0 / np.maximum(var, 1e-9))
        gs, ws = np.array(gs), np.array(ws)
        only = np.abs((gs * ws).sum(axis=0) / np.maximum(ws.sum(axis=0), 1e-12))
        use = mixed[1] & silence[1] & np.isfinite(only) & np.isfinite(truth)
        if use.sum() < 100:
            continue
        for name, est in (("images only", only), ("mixed (ships)", mixed[0]),
                          ("silence only", silence[0])):
            rows.setdefault(name, []).append(score(est, truth, use))
        print(f"  split {split + 1}/{N_SPLITS}: {coord_count} coordinate studies, "
              f"{n_foci / max(coord_count, 1):.1f} foci each", flush=True)

    print(f"\n{'estimate':>16} {'r':>7} {'rank r':>7} {'AUC':>6} {'ratio':>6} "
          f"{'mean err':>9} {'err at top':>11} {'rmse':>7}")
    for name in ("images only", "mixed (ships)", "silence only"):
        if name not in rows:
            continue
        r, rho, auc, ratio, err, err_top, rmse = np.nanmean(np.array(rows[name]), axis=0)
        print(f"{name:>16} {r:+7.3f} {rho:+7.3f} {auc:6.3f} {ratio:6.2f} {err:+9.3f} "
              f"{err_top:+11.3f} {rmse:7.3f}")
    names = ("r", "rank r", "AUC", "ratio", "mean err", "err at top", "rmse")
    if "silence only" in rows:
        print(f"\nPaired across splits (n = {len(rows['silence only'])}):")
        for other in ("mixed (ships)", "images only"):
            if other not in rows:
                continue
            a, b = np.array(rows["silence only"]), np.array(rows[other])
            print(f"\n  silence-only minus {other}")
            for j, metric in enumerate(names):
                d = a[:, j] - b[:, j]
                ok = np.isfinite(d)
                if ok.sum() < 2:
                    continue
                p = float(stats.ttest_rel(a[ok, j], b[ok, j]).pvalue)
                print(f"    {metric:>11} {np.mean(d):+9.3f}  p {p:.4f}")
    print("\nThe question this was built for: does the pain conclusion hold where the truth comes")
    print("from different subjects and the construct is motor rather than pain. And whether the")
    print("ratio column is near enough to 1 for a single image to support an absolute magnitude.")
