"""The matched comparison again, stratified by how much signal a voxel actually carries.

The reference-bias calibration showed the conditional reference is inflated about sevenfold at
voxels where the true effect is near zero: conditioning a study on clearing its own threshold
at a null voxel selects a pure noise excursion, and the mean |g| among those is large. Most
voxels in a brain are like that, so an unrestricted comparison is dominated by them -- and
CBES, fed peaks that are themselves selected, may be estimating the same selected noise. Two
inflated numbers agreeing says nothing.

So the comparison is stratified. The pooled image map over 21 studies is a far better indicator
of where signal is than any single study, and it is used only to *choose* strata, never as the
quantity compared against; the comparison within each stratum is still CBES's g against the
matched conditional reference. A second stratification by how many studies are active at the
voxel gives a consensus criterion that does not depend on the pooled map at all.

The threshold is varied rather than the peak count capped: capping fixes the count and lets the
cut float, and that cut is then each study's Nth peak, which tracks its signal strength.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from scipy.ndimage import maximum_filter
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer

THRESHOLDS = (3.2905, 4.0)
MIN_VOXELS = 100

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
mask_bool = np.asarray(masker.mask_img.get_fdata() > 0)
shape = mask_bool.shape
affine = masker.mask_img.affine

g_rows, z_rows, w_rows, ids = [], [], [], []
for row in ss.images.itertuples():
    if row.g is None or row.g_var is None:
        continue
    g = masker.transform(str(row.g)).ravel()
    v = masker.transform(str(row.g_var)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    g = np.where(ok, g, 0.0)
    g_rows.append(g)
    z_rows.append(np.where(ok, g / np.sqrt(np.where(ok, v, np.inf)), 0.0))
    w_rows.append(np.where(ok, 1.0 / np.maximum(v, 1e-12), 0.0))
    ids.append(str(row.id))
G, Z, W = np.array(g_rows), np.array(z_rows), np.array(w_rows)
pooled = np.sum(G * W, axis=0) / np.maximum(W.sum(axis=0), 1e-12)
print(f"{len(G)} studies, {G.shape[1]} voxels; pooled |g| median "
      f"{np.median(np.abs(pooled)):.3f}\n")


def peaks_above(z_masked, threshold):
    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z_masked
    magnitude = np.abs(volume)
    is_peak = (magnitude == maximum_filter(magnitude, size=3)) & (magnitude >= threshold) \
        & mask_bool
    idx = np.argwhere(is_peak)
    if not len(idx):
        return []
    values = volume[tuple(idx.T)]
    order = np.argsort(-np.abs(values))
    return [(idx[i], float(values[i])) for i in order]


for threshold in THRESHOLDS:
    studies, kept = [], []
    for pos, sid in enumerate(ids):
        found = peaks_above(Z[pos], threshold)
        if len(found) < 2:
            continue
        kept.append(pos)
        meta = {"sample_sizes": [int(sizes.get(sid, 30))]}
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [
            {"id": sid, "name": "1", "metadata": meta,
             "points": [{"space": "MNI",
                         "coordinates": [float(c) for c in nib.affines.apply_affine(
                             affine, np.asarray(ijk, dtype=float))],
                         "values": [{"kind": "Z", "value": value}]} for ijk, value in found]}]})
    if len(studies) < 5:
        print(f"threshold {threshold}: only {len(studies)} studies report, skipped")
        continue
    rows = np.array(kept)
    active = np.abs(Z[rows]) >= threshold
    n_active = active.sum(axis=0)
    with np.errstate(invalid="ignore"):
        mu = np.where(n_active > 0,
                      (np.abs(G[rows]) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="study-min")
    result = est.fit(Studyset({"id": "p", "name": "p", "studies": studies}, target=None,
                              mask=masker.mask_img))
    g = np.abs(result.get_map("g", return_type="array").ravel())
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    inferred = float(np.median(np.abs(est._cutoffs_z_.values)))
    base = covered & (n_active >= 2) & np.isfinite(mu) & (g > 0)
    print(f"--- threshold {threshold:.2f} (inferred {inferred:.2f}), {len(studies)} studies, "
          f"{int(base.sum())} comparable voxels ---")

    absolute = np.abs(pooled)
    edges = [0, 50, 75, 90, 99, 100]
    print(f"  {'pooled-truth stratum':>22} {'vox':>7} {'|pooled|':>9} {'mu_hat':>8} "
          f"{'CBES g':>8} {'ratio':>7} {'r':>7}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        band = (absolute >= np.percentile(absolute, lo)) & (absolute < np.percentile(absolute, hi)) \
            if hi < 100 else (absolute >= np.percentile(absolute, lo))
        use = base & band
        if use.sum() < MIN_VOXELS:
            print(f"  {f'{lo}-{hi}%':>22} {int(use.sum()):>7}   too few")
            continue
        r = stats.pearsonr(g[use], mu[use])[0]
        print(f"  {f'{lo}-{hi}%':>22} {int(use.sum()):>7} {absolute[use].mean():9.3f} "
              f"{mu[use].mean():8.3f} {g[use].mean():8.3f} "
              f"{g[use].mean()/mu[use].mean():7.3f} {r:7.3f}")

    print(f"  {'consensus stratum':>22} {'vox':>7} {'|pooled|':>9} {'mu_hat':>8} "
          f"{'CBES g':>8} {'ratio':>7} {'r':>7}")
    for k in (2, 3, 5, 8):
        use = base & (n_active >= k)
        if use.sum() < MIN_VOXELS:
            print(f"  {f'>= {k} studies active':>22} {int(use.sum()):>7}   too few")
            continue
        r = stats.pearsonr(g[use], mu[use])[0]
        print(f"  {f'>= {k} studies active':>22} {int(use.sum()):>7} {absolute[use].mean():9.3f} "
              f"{mu[use].mean():8.3f} {g[use].mean():8.3f} "
              f"{g[use].mean()/mu[use].mean():7.3f} {r:7.3f}", flush=True)
    print(flush=True)
