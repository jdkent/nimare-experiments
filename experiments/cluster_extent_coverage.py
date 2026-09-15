"""Is the compression caused by reading a study's own significant territory as silence?

``silence_misread`` established the fact: under cluster-extent reporting, 46% of the voxels a
study found significant sit outside the coverage radius of the one focus its table gives for
that cluster, and the estimator enters every one of them as the study having been *silent*
there -- positive evidence against an effect. The fraction tracks cluster size, so it is largest
where the true effect is largest, which is the shape a mechanism for the compression would have.

Fact is not cause. This tests it by removing the error and nothing else. Each study is given an
analysis mask -- the estimator's existing device for "this study never looked here" -- covering
everything except its own significant-but-uncovered territory. Silence there stops being read as
evidence; every value, weight, threshold and focus stays exactly as it was.

If the compression is caused by the misreading, the stratum ratios should fall toward one as the
truth rises. If it barely moves, the flatness is something else and the coverage model is not
worth changing.

The mask here is built from each study's full image, which a real meta-analysis does not have.
That is the point: it measures the size of the prize before anything is built to claim part of
it from a reported cluster extent.
"""
import logging, os, shutil, sys, tempfile, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks, cluster_extent_threshold, CLUSTER_FORMING_Z

SCHEME, FOCUS = "cluster", "max"
COVERAGE_RADIUS_MM = 20.0
N_SPLITS = 6
MIN_HALF = 5

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
grid = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing="ij"), axis=-1)
workdir = tempfile.mkdtemp(prefix="extent-coverage-")


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled_truth(maps, sizes):
    sizes = np.asarray(sizes, dtype=float)[:, None]
    stack = np.array([to_g(z, float(n)) for z, n in zip(maps, sizes.ravel())])
    weights = 1.0 / np.maximum(1.0 / sizes + stack**2 / (2.0 * sizes), 1e-9)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def significant_set(z, height):
    volume = np.zeros(shape); volume[mask_bool] = z
    supra = (np.abs(volume) >= height) & mask_bool
    labels, n = ndimage.label(supra)
    if not n:
        return supra
    sizes = np.bincount(labels.ravel())
    critical = cluster_extent_threshold(volume, mask_bool, zooms, CLUSTER_FORMING_Z)
    keep = np.zeros(sizes.size, bool)
    keep[sizes >= critical] = True
    keep[0] = False
    return supra & keep[labels]


def build(maps, sizes, with_mask):
    """The studyset, optionally telling each study not to claim silence in its own clusters."""
    studies = []
    for k, (z, n) in enumerate(zip(maps, sizes)):
        foci, height = report_peaks(z, mask_bool, shape, zooms, scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(n)], "reporting_threshold": float(height)}
        analysis = {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
            {"space": "MNI",
             "coordinates": [float(c) for c in nib.affines.apply_affine(
                 affine, np.asarray(ijk, dtype=float))],
             "values": [{"kind": "Z", "value": v}]} for ijk, v in foci]}
        if with_mask:
            covered = np.zeros(shape, dtype=bool)
            for ijk, _ in foci:
                covered |= np.linalg.norm(
                    (grid - np.asarray(ijk)) * zooms, axis=-1) <= COVERAGE_RADIUS_MM
            examined = mask_bool & ~(significant_set(z, height) & ~covered)
            path = os.path.join(workdir, f"examined_{k}.nii.gz")
            nib.save(nib.Nifti1Image(examined.astype(np.int16), affine), path)
            analysis["images"] = [{"url": path, "filename": os.path.basename(path),
                                   "value_type": "analysis_mask", "space": "MNI"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    return studies


def run(maps, sizes, with_mask):
    studies = build(maps, sizes, with_mask)
    if len(studies) < MIN_HALF:
        return None
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="reporting_threshold",
               **({"analysis_mask": "analysis_mask"} if with_mask else {}))
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    g = np.abs(res.get_map("g", return_type="array").ravel())
    pi = (res.get_map("prevalence", return_type="array").ravel()
          if "prevalence" in res.maps else np.ones_like(g))
    covered = res.get_map("n_studies", return_type="array").ravel() > 0
    return g, pi, covered


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
print(f"NIDM pain: {n} studies, scheme={SCHEME}, focus={FOCUS}\n", flush=True)

BANDS = ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100))
acc = {False: [], True: []}
corr = {False: [], True: []}
rng = np.random.default_rng(0)
for _ in range(N_SPLITS):
    order = rng.permutation(n)
    lo, hi = order[: n // 2], order[n // 2:]
    truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
    for with_mask in (False, True):
        got = run([maps[i] for i in lo], sizes[lo], with_mask)
        if got is None:
            continue
        g, pi, cov = got
        use = cov & np.isfinite(g) & (g > 0)
        if use.sum() < 100:
            continue
        corr[with_mask].append(stats.pearsonr(g[use], truth[use])[0])
        cells = []
        for a, b in BANDS:
            band = (truth >= np.percentile(truth, a)) & (
                truth < np.percentile(truth, b) if b < 100 else np.ones_like(truth, bool))
            pick = use & band
            ok = pick.sum() >= 30
            cells.append((truth[pick].mean() if ok else np.nan,
                          g[pick].mean() if ok else np.nan))
        acc[with_mask].append(cells)

print(f"  {'truth stratum':>16} {'truth g':>9} "
      f"{'g as now':>9} {'ratio':>7}   {'g, silence fixed':>17} {'ratio':>7}")
blocks = {k: np.array(v, dtype=float) for k, v in acc.items() if v}
for j, name in enumerate(("0-50%", "50-75%", "75-90%", "90-99%", "99-100%")):
    t = np.nanmean(blocks[False][:, j, 0])
    a = np.nanmean(blocks[False][:, j, 1])
    b = np.nanmean(blocks[True][:, j, 1]) if True in blocks else np.nan
    print(f"  {name:>16} {t:9.3f} {a:9.3f} {a / max(t, 1e-9):7.2f}   "
          f"{b:17.3f} {b / max(t, 1e-9):7.2f}")
print(f"\n  correlation with the truth: as now {np.mean(corr[False]):+.3f}, "
      f"silence fixed {np.mean(corr[True]):+.3f}")
print("\nRatios falling toward 1.00 as the truth rises means the misreading was the cause.")
shutil.rmtree(workdir, ignore_errors=True)
