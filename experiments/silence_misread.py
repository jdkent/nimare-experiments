"""How much of a study's significant territory does the model read as silence?

With one focus per surviving cluster -- which is what a table gives -- a 40 mm cluster is
represented by a single coordinate. The estimator decides a study was *silent* at a voxel by
asking whether any of its foci lie within ``coverage_radius``, so everything in that cluster
beyond the radius is entered as evidence that the study had no effect there, when in fact the
study found a significant one.

That error is not uniform. A bigger true effect makes a bigger cluster, and a bigger cluster has
proportionally more territory outside the radius of its one focus. So the misreading grows with
the effect, and it pushes the estimate down hardest exactly where the effect is largest -- which
is a mechanism for the compression measured against every held-out reference so far.

This quantifies it without touching the estimator: take each study's real thresholded map, take
the coordinates a table would carry, and compare the significant set against the covered set.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.transforms import ImageTransformer
from load_pain import load_pain
from reporting import report_peaks, CLUSTER_FORMING_Z, cluster_extent_threshold

COVERAGE_RADIUS_MM = 20.0     # the estimator's default: 2 x fwhm
mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape = mask_img.shape
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)

grid = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing="ij"), axis=-1)


def covered_mask(foci):
    """Voxels within the coverage radius of any focus, as the estimator computes silence."""
    out = np.zeros(shape, dtype=bool)
    for ijk, _ in foci:
        out |= np.linalg.norm((grid - np.asarray(ijk)) * zooms, axis=-1) <= COVERAGE_RADIUS_MM
    return out


ss = ImageTransformer(target="z").transform(load_pain())
rows = []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    path = getattr(row, "z", None)
    if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
        continue
    img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                          force_resample=True, copy_header=True)
    rows.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
print(f"NIDM pain: {len(rows)} studies\n")

for scheme in ("fdr", "fwe", "cluster"):
    print(f"--- {scheme} ---")
    print(f"  {'study':>6} {'foci':>5} {'significant vox':>16} {'covered':>8} "
          f"{'read as silent':>15} {'largest cluster vox':>20}")
    totals = []
    for k, z in enumerate(rows):
        foci, height = report_peaks(z, mask_bool, shape, zooms, scheme=scheme, focus="max")
        if not foci:
            continue
        volume = np.zeros(shape); volume[mask_bool] = z
        supra = (np.abs(volume) >= height) & mask_bool
        if scheme == "cluster":
            labels, n_labels = ndimage.label(supra)
            sizes = np.bincount(labels.ravel())
            k_crit = cluster_extent_threshold(volume, mask_bool, zooms, CLUSTER_FORMING_Z)
            keep = np.zeros(sizes.size, bool)
            keep[sizes >= k_crit] = True
            keep[0] = False
            supra &= keep[labels]
        if not supra.any():
            continue
        cov = covered_mask(foci)
        missed = int((supra & ~cov).sum())
        biggest = int(np.bincount(ndimage.label(supra)[0].ravel())[1:].max())
        totals.append(missed / supra.sum())
        if k < 8:
            print(f"  {k:>6} {len(foci):>5} {int(supra.sum()):>16} "
                  f"{int((supra & cov).sum()):>8} {missed / supra.sum():>14.1%} {biggest:>20}")
    if totals:
        print(f"  mean across {len(totals)} studies: "
              f"{np.mean(totals):.1%} of significant territory entered as silence\n", flush=True)
