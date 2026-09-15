"""What does `coverage_radius` buy, and what does it cost?

The convergence comparison turned on a number: CBES covers 9.2% of the brain and 34.1% of the
truth's top decile, so on its own support `g_marginal` localises better than any convergence
statistic (+0.111 AUC) and over the whole brain it loses, because two-thirds of the signal lies
where CBES returns exactly zero.

**The lever is `fwhm`, not `coverage_radius`** -- which a first version of this script got wrong.
A voxel gets an estimate when some study's focus reaches it through the *pooling kernel*, whose
support `fwhm` sets; `coverage_radius` decides something else entirely, namely how far from a
focus a study is taken to have been *silent* rather than uninformative, which changes the
estimates without changing where they exist. Sweeping `coverage_radius` from 8 mm to 45 mm leaves
the covered share of the brain at 0.089 to three decimals at every radius, which is how the
error showed itself. So both are swept here, for their own questions.

Two scores, because the trade-off has two sides:

  on the covered voxels   does the estimate localise well where it is offered
  on the whole mask       does the map localise well over the brain a reader looks at

A radius that is too small wins the first and loses the second; one too large should do the
reverse, and somewhere in between is a defensible default. Reported alongside the share of the
brain covered and the share of the truth's top decile captured, so the choice can be read rather
than inferred.
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
from nimare.meta.cbma import CBES, MKDADensity
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_SPLITS = int(os.environ.get("NSPLITS", 8))
#: Kernel widths. This is what sets how much of the brain gets an estimate at all.
FWHMS = (6.0, 10.0, 16.0, 24.0)
#: Silence radii, at the default kernel. This changes the estimates, not their extent.
RADII = (8.0, 20.0, 30.0)                  # 20 mm is the default, being 2 * fwhm

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled(members, maps, sizes):
    stack = np.array([to_g(maps[i], sizes[i]) for i in members])
    w = np.array([1.0 / np.maximum(1.0 / sizes[i] + to_g(maps[i], sizes[i]) ** 2
                                   / (2 * sizes[i]), 1e-9) for i in members])
    return np.sum(stack * w, axis=0) / np.maximum(w.sum(axis=0), 1e-12)


def build(members, maps, sizes):
    studies = []
    for i in members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms, SCHEME, FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(sizes[i])], "reporting_threshold": float(height)}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(c) for c in
                 nib.affines.apply_affine(affine, np.asarray(ijk, float))],
                 "values": [{"kind": "Z", "value": float(v)}]} for ijk, v in foci]}]})
    return (Studyset({"id": "r", "name": "r", "studies": studies}, target=None, mask=mask_img)
            if len(studies) >= 2 else None)


def auc(est, truth, use):
    a, t = est[use], truth[use]
    if not np.isfinite(a).all() or np.allclose(a, a[0]):
        return np.nan
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    return (float((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
            if n_pos and n_neg else np.nan)


if __name__ == "__main__":
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
    sizes = np.asarray(sizes, float)
    total = len(maps)
    print(f"NIDM pain: {total} studies, {SCHEME}/{FOCUS} reporting, {total//2} per half, "
          f"{N_SPLITS} splits")
    print("The default coverage_radius is 2 * fwhm = 20 mm.\n")

    by_fwhm = {f: [] for f in FWHMS}
    by_radius = {r: [] for r in RADII}
    mkda = []
    rng = np.random.default_rng(0)
    for _ in range(N_SPLITS):
        order = rng.permutation(total)
        work, hold = order[: total // 2], order[total // 2:]
        truth = np.abs(pooled(hold, maps, sizes))
        studyset = build(list(work), maps, sizes)
        if studyset is None:
            continue
        whole = np.isfinite(truth)
        top = truth >= np.percentile(truth[whole], 90)
        out = MKDADensity(null_method="approximate", mask=masker).fit(studyset)
        mkda.append((auc(out.get_map("stat", return_type="array").ravel(), truth, whole),))

        def score(est):
            res = est.fit(studyset)
            g = np.abs(res.get_map("g", return_type="array").ravel())
            pi = (res.get_map("prevalence", return_type="array").ravel()
                  if "prevalence" in res.maps else np.ones_like(g))
            marginal = g * pi
            covered = res.get_map("n_studies", return_type="array").ravel() > 0
            return (float(covered.mean()), float(covered[top].mean()),
                    auc(marginal, truth, covered & whole), auc(marginal, truth, whole))

        for width in FWHMS:
            by_fwhm[width].append(score(CBES(
                fwhm=width, mask=masker, peak_bias=None, null_method="none",
                threshold="reporting_threshold")))
        for radius in RADII:
            by_radius[radius].append(score(CBES(
                fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
                threshold="reporting_threshold", coverage_radius=radius)))

    header = (f"{'setting':>14s} {'brain covered':>14s} {'top decile covered':>19s} "
              f"{'AUC on covered':>15s} {'AUC whole mask':>15s}")

    print("--- kernel width, which sets how much of the brain gets an estimate ---")
    print(header)
    for width in FWHMS:
        a = np.array(by_fwhm[width], dtype=float)
        if not a.size:
            continue
        tag = "  <- default" if width == 10.0 else ""
        print(f"{f'fwhm {width:.0f} mm':>14s} {np.nanmean(a[:,0]):14.3f} "
              f"{np.nanmean(a[:,1]):19.3f} {np.nanmean(a[:,2]):15.3f} "
              f"{np.nanmean(a[:,3]):15.3f}{tag}")
    m = np.array(mkda, dtype=float)
    print(f"{'MKDA':>14s} {'1.000':>14s} {'1.000':>19s} {'--':>15s} "
          f"{np.nanmean(m[:,0]):15.3f}")

    print("\n--- silence radius at fwhm 10, which changes the estimates not their extent ---")
    print(header)
    for radius in RADII:
        a = np.array(by_radius[radius], dtype=float)
        if not a.size:
            continue
        tag = "  <- default" if radius == 20.0 else ""
        print(f"{f'radius {radius:.0f} mm':>14s} {np.nanmean(a[:,0]):14.3f} "
              f"{np.nanmean(a[:,1]):19.3f} {np.nanmean(a[:,2]):15.3f} "
              f"{np.nanmean(a[:,3]):15.3f}{tag}")

    print("\nThe whole-mask column is the one the convergence comparison turned on: CBES lost")
    print("there because it returns zero where no kernel reaches. If a wider kernel closes that")
    print("gap without costing much on the covered voxels, the default is simply too narrow.")
