"""Turn a group z map into the coordinate table a paper would actually print.

Everything measured so far extracted peaks by fixing an uncorrected height cut and keeping the
strongest ten. That is not what anyone reports, and it gets the two things that matter wrong:
the cut is not uncorrected, and the count is not chosen by the author. Capping at ten fixes the
count and lets the effective cut float to whatever the tenth peak happened to be, which makes
the cut a function of the signal -- the confound that invalidated the earlier threshold fit.

What a paper does instead:

  * **Correct for multiplicity.** Either a voxelwise FDR at q = 0.05, or a cluster-forming cut
    at p < 0.001 uncorrected kept only where the cluster is large enough to survive, or a
    voxelwise family-wise cut. Which one is used varies; all three are represented here.
  * **Keep whole clusters.** Peaks come only from clusters that survived, so a study with no
    surviving cluster reports nothing at all -- a case the capped extraction could never
    produce, and one the selection model's silence term exists to handle.
  * **Separate the peaks.** SPM and FSL tabulate local maxima at least 8 mm apart, so adjacent
    voxels are not separate rows in a table. No limit is placed on how many come back. Their
    conventional "three per cluster" is a cap, and a cap re-introduces the confound this module
    exists to remove: dropping the fourth peak of a big cluster raises the smallest reported
    value, so the effective cut moves with the signal again.

The height threshold that comes back is the real one for that map, so it can be handed to the
estimator as metadata instead of being inferred from the smallest reported value.
"""
import numpy as np
from scipy import ndimage, stats

MIN_SEPARATION_MM = 8.0


def fdr_height(z, q=0.05):
    """Two-sided Benjamini-Hochberg height threshold on |z|, or inf if nothing survives."""
    p = 2.0 * stats.norm.sf(np.abs(z[np.isfinite(z)]))
    if not p.size:
        return np.inf
    p = np.sort(p)
    passed = np.flatnonzero(p <= q * np.arange(1, p.size + 1) / p.size)
    if not passed.size:
        return np.inf
    return float(stats.norm.isf(p[passed[-1]] / 2.0))


def bonferroni_height(z, alpha=0.05):
    n = int(np.isfinite(z).sum())
    return float(stats.norm.isf(alpha / (2.0 * max(n, 1))))


def report_peaks(z_masked, mask_bool, shape, zooms, scheme="fdr",
                 cluster_forming=3.0902, min_cluster_voxels=10, q=0.05):
    """Return ``(peaks, height)``: the table a paper would print, and the cut behind it.

    ``peaks`` is a list of ``(ijk, signed z)``, every local maximum that survived, however many
    that is. ``height`` is the height threshold actually applied, which is what the estimator
    should be told rather than left to infer.

    ``scheme`` is one of ``'fdr'`` (voxelwise Benjamini-Hochberg), ``'fwe'`` (voxelwise
    Bonferroni) or ``'cluster'`` (p < 0.001 forming cut, clusters of at least
    ``min_cluster_voxels``). An empty list is a real outcome, not a failure: it is a study whose
    map had nothing to report.
    """
    if scheme == "fdr":
        height = fdr_height(z_masked, q)
    elif scheme == "fwe":
        height = bonferroni_height(z_masked, q)
    elif scheme == "cluster":
        height = float(cluster_forming)
    else:
        raise ValueError(f"unknown scheme {scheme!r}")
    if not np.isfinite(height):
        return [], height

    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z_masked
    magnitude = np.abs(volume)
    supra = (magnitude >= height) & mask_bool
    if not supra.any():
        return [], height

    labels, n_labels = ndimage.label(supra)
    if scheme == "cluster":
        sizes = np.bincount(labels.ravel())
        keep = np.zeros(sizes.size, dtype=bool)
        keep[min_cluster_voxels <= sizes] = True
        keep[0] = False
        supra &= keep[labels]
        if not supra.any():
            return [], height
        labels, n_labels = ndimage.label(supra)

    # Local maxima of |z| inside the surviving clusters, strongest first.
    local = (magnitude == ndimage.maximum_filter(magnitude, size=3)) & supra
    idx = np.argwhere(local)
    if not len(idx):
        return [], height
    strength = magnitude[tuple(idx.T)]
    order = np.argsort(-strength)

    zooms = np.asarray(zooms, dtype=float)
    kept = []
    for i in order:
        ijk = idx[i]
        if kept and np.any(np.linalg.norm(
                (np.array([k for k, _ in kept]) - ijk) * zooms, axis=1) < MIN_SEPARATION_MM):
            continue
        kept.append((ijk, float(volume[tuple(ijk)])))
    return kept, height
