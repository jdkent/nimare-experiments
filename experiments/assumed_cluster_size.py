"""Does assuming a cluster size around each focus reduce the compression?

Cancelling the misread silence bought about 1%, because a voxel deep inside a cluster has no
focus within kernel reach either and is decided by other studies. The stronger version of the
same idea is to let the assumed cluster carry the study's reported *value*, not merely cancel
its silence -- so the focus speaks for the whole cluster rather than for a 13 mm neighbourhood
of its peak.

The estimator already has the knob: ``fwhm`` sets how far a reported value is spread, and the
coverage radius follows it at twice its size. Sweeping it is exactly "assume every cluster is
this big", applied to the value and the silence together, and it needs no reported extent --
which is the point, since papers do not give one reliably.

The statistical cost is real and worth naming. A 25 mm kernel asserts the effect is roughly
constant across a region that large, which for a 7000-voxel cluster is defensible and for the
median 131-voxel one is not. If the compression does not move across the sweep, the cost need
not be argued about.

Scored the same way as everything else: pain studies split in half, coordinates from one half,
inverse-variance pooled truth from the other, cluster-extent reporting with a family-wise
extent test.
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
FWHMS = (10.0, 15.0, 20.0, 25.0)
N_SPLITS = 6
MIN_HALF = 5
BANDS = ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100))

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled_truth(maps, sizes):
    sizes = np.asarray(sizes, dtype=float)[:, None]
    stack = np.array([to_g(z, float(n)) for z, n in zip(maps, sizes.ravel())])
    weights = 1.0 / np.maximum(1.0 / sizes + stack**2 / (2.0 * sizes), 1e-9)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def build(maps, sizes):
    studies = []
    for k, (z, n) in enumerate(zip(maps, sizes)):
        foci, height = report_peaks(z, mask_bool, shape, zooms, scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(n)], "reporting_threshold": float(height)}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": v}]} for ijk, v in foci]}]})
    return studies


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
print(f"NIDM pain: {n} studies, {SCHEME} reporting, focus at the {FOCUS}\n", flush=True)

acc = {f: {"strata": [], "r": []} for f in FWHMS}
rng = np.random.default_rng(0)
for _ in range(N_SPLITS):
    order = rng.permutation(n)
    lo, hi = order[: n // 2], order[n // 2:]
    truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
    studies = build([maps[i] for i in lo], sizes[lo])
    if len(studies) < MIN_HALF:
        continue
    bands = [(truth >= np.percentile(truth, a)) &
             (truth < np.percentile(truth, b) if b < 100 else np.ones_like(truth, bool))
             for a, b in BANDS]
    for fwhm in FWHMS:
        est = CBES(fwhm=fwhm, mask=masker, peak_bias=None, null_method="none",
                   threshold="reporting_threshold")
        res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                               target=None, mask=mask_img))
        g = np.abs(res.get_map("g", return_type="array").ravel())
        use = (res.get_map("n_studies", return_type="array").ravel() > 0) \
            & np.isfinite(g) & (g > 0)
        if use.sum() < 100:
            continue
        acc[fwhm]["r"].append(stats.pearsonr(g[use], truth[use])[0])
        acc[fwhm]["strata"].append([
            (truth[use & b].mean(), g[use & b].mean()) if (use & b).sum() >= 30
            else (np.nan, np.nan) for b in bands])

print(f"  {'fwhm':>6} {'coverage':>9} {'r':>7}   "
      + " ".join(f"{s:>9}" for s in ("0-50%", "50-75%", "75-90%", "90-99%", "99-100%")))
print(f"  {'':>6} {'':>9} {'':>7}   " + " ".join(f"{'ratio':>9}" for _ in BANDS))
for fwhm in FWHMS:
    rec = acc[fwhm]
    if not rec["strata"]:
        print(f"  {fwhm:6.0f}   no usable split")
        continue
    block = np.array(rec["strata"], dtype=float)
    ratios = [np.nanmean(block[:, j, 1]) / max(np.nanmean(block[:, j, 0]), 1e-9)
              for j in range(len(BANDS))]
    print(f"  {fwhm:6.0f} {2 * fwhm:8.0f}mm {np.mean(rec['r']):+7.3f}   "
          + " ".join(f"{r:9.2f}" for r in ratios))
truths = np.nanmean(np.array(acc[FWHMS[0]]["strata"], dtype=float)[:, :, 0], axis=0)
print(f"\n  held-out truth in each stratum: "
      + " ".join(f"{t:.3f}" for t in truths)
      + f"  ({truths[-1] / max(truths[0], 1e-9):.0f}-fold range)")
print("\nRatios converging toward 1.00 across the strata would mean the assumption helps.")
