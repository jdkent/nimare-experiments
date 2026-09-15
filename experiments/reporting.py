"""Turn a group z map into the coordinate table a paper would actually print.

Three things decide what ends up in a table, and getting any of them wrong makes the extraction
a different experiment from the one intended.

**The threshold is corrected.** Nobody tabulates an uncorrected height. Either a voxelwise
Benjamini-Hochberg FDR, a voxelwise family-wise cut, or -- most commonly -- a cluster-forming
cut at p < 0.001 with only those clusters large enough to survive a family-wise cluster-extent
test. The last is computed here from random field theory rather than a made-up minimum size, so
the surviving set is the one SPM or FSL would have kept.

**Only surviving clusters are reported.** A study whose map has nothing left reports nothing at
all and drops out of the meta-analysis. That is a real outcome, and it is the one the selection
model's silence term exists to handle.

**One focus per cluster.** A table gives a cluster's location as a single row: either the voxel
holding the maximum statistic, or the cluster's centre of mass, with the statistic there. Both
are in use, they differ substantially -- the centre of mass sits well below the peak -- and both
are available here through ``focus``.

Nothing caps the number of rows. A cap fixes the count and lets the effective threshold float to
whatever the last kept peak happened to be, which makes the cut a function of the signal and
contaminates everything downstream of it. The count is whatever survives the correction.
"""
import numpy as np
from scipy import ndimage, stats
from scipy.special import gamma as gamma_function

#: Cluster-forming height for the cluster-extent scheme, p < 0.001 one-sided.
CLUSTER_FORMING_Z = 3.0902

#: Family-wise error rate the cluster-extent and voxelwise family-wise schemes control.
ALPHA = 0.05


def fdr_height(z, q=ALPHA):
    """Two-sided Benjamini-Hochberg height threshold on ``|z|``, or inf if nothing survives."""
    p = 2.0 * stats.norm.sf(np.abs(z[np.isfinite(z)]))
    if not p.size:
        return np.inf
    p = np.sort(p)
    passed = np.flatnonzero(p <= q * np.arange(1, p.size + 1) / p.size)
    if not passed.size:
        return np.inf
    return float(stats.norm.isf(p[passed[-1]] / 2.0))


def bonferroni_height(z, alpha=ALPHA):
    """Two-sided voxelwise family-wise height threshold."""
    return float(stats.norm.isf(alpha / (2.0 * max(int(np.isfinite(z).sum()), 1))))


def smoothness_fwhm(volume, mask_bool, zooms):
    """Estimated FWHM in mm per axis, from the variance of the map's own gradient.

    The standard estimator: for a Gaussian field, ``Var(dz/dx) / Var(z) = 4 ln 2 / FWHM^2`` in
    voxel units. Computed on the map rather than on residuals, which are not available when all
    that survives of a study is its group statistic image, so this is an approximation -- but
    the alternative is to invent a minimum cluster size, which is worse.
    """
    z = np.where(mask_bool, volume, np.nan)
    variance = np.nanvar(z)
    fwhm = []
    for axis, spacing in enumerate(zooms):
        difference = np.diff(z, axis=axis)
        lam = np.nanvar(difference) / (2.0 * max(variance, 1e-12))
        lam = min(max(lam, 1e-6), 0.5)
        fwhm.append(float(spacing * np.sqrt(4.0 * np.log(2.0) / (-2.0 * np.log(1.0 - lam)))
                          if lam < 1.0 else spacing))
    return np.array(fwhm)


def cluster_extent_threshold(volume, mask_bool, zooms, height=CLUSTER_FORMING_Z, alpha=ALPHA):
    """Smallest cluster, in voxels, that survives a family-wise extent test at ``alpha``.

    Random field theory, as SPM computes it. With ``R`` resels in the search volume and a
    Gaussian field thresholded at ``u``, the expected number of clusters is ``R`` times the
    three-dimensional Euler characteristic density, the expected number of supra-threshold
    voxels is the search volume times the tail probability, and a cluster's size is
    approximately ``exp(-beta k^(2/3))`` distributed. Inverting the family-wise probability at
    ``alpha`` gives the critical size.
    """
    voxel_volume = float(np.prod(zooms))
    n_voxels = int(mask_bool.sum())
    if not n_voxels:
        return np.inf
    fwhm = smoothness_fwhm(volume, mask_bool, zooms)
    resels = n_voxels * voxel_volume / float(np.prod(fwhm))

    u = float(height)
    # Euler characteristic density for a 3D Gaussian field, two-sided via the factor of two.
    ec_density = (4.0 * np.log(2.0)) ** 1.5 / (2.0 * np.pi) ** 2 * (u * u - 1.0) \
        * np.exp(-u * u / 2.0)
    expected_clusters = 2.0 * resels * ec_density
    if expected_clusters <= 0:
        return np.inf
    expected_supra = n_voxels * 2.0 * stats.norm.sf(u)
    if expected_supra <= 0:
        return np.inf
    expected_size = expected_supra / expected_clusters

    beta = (gamma_function(1.0 + 3.0 / 2.0) / max(expected_size, 1e-9)) ** (2.0 / 3.0)
    target = -np.log1p(-alpha) / expected_clusters
    if not (0.0 < target < 1.0):
        return np.inf
    return float((-np.log(target) / beta) ** 1.5)


def report_peaks(z_masked, mask_bool, shape, zooms, scheme="cluster", focus="max"):
    """Return ``(foci, height)`` -- the rows a paper would print, and the cut behind them.

    ``foci`` is a list of ``(ijk, signed z)``, one entry per surviving cluster. ``focus`` picks
    which voxel represents the cluster: ``'max'`` for the strongest statistic in it, ``'com'``
    for its centre of mass. ``height`` is the height threshold actually applied, which the
    estimator should be told rather than left to infer.

    ``scheme`` is ``'fdr'`` (voxelwise Benjamini-Hochberg), ``'fwe'`` (voxelwise Bonferroni) or
    ``'cluster'`` (p < 0.001 forming cut, family-wise extent test). An empty list means the
    study had nothing to report, which is a real outcome rather than a failure.
    """
    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z_masked
    magnitude = np.abs(volume)

    if scheme == "fdr":
        height, min_size = fdr_height(z_masked), 1
    elif scheme == "fwe":
        height, min_size = bonferroni_height(z_masked), 1
    elif scheme == "cluster":
        height = CLUSTER_FORMING_Z
        min_size = cluster_extent_threshold(volume, mask_bool, zooms, height)
    else:
        raise ValueError(f"unknown scheme {scheme!r}")
    if not np.isfinite(height) or not np.isfinite(min_size):
        return [], height

    supra = (magnitude >= height) & mask_bool
    if not supra.any():
        return [], height
    labels, n_labels = ndimage.label(supra)
    if not n_labels:
        return [], height

    sizes = np.bincount(labels.ravel())
    survivors = [c for c in range(1, n_labels + 1) if sizes[c] >= min_size]
    if not survivors:
        return [], height

    foci = []
    for cluster in survivors:
        member = labels == cluster
        if focus == "max":
            ijk = np.unravel_index(np.argmax(np.where(member, magnitude, -np.inf)), shape)
            ijk = np.array(ijk, dtype=np.int64)
        elif focus == "com":
            ijk = np.rint(ndimage.center_of_mass(member)).astype(np.int64)
            ijk = np.clip(ijk, 0, np.array(shape) - 1)
            if not member[tuple(ijk)]:
                # A horseshoe-shaped cluster's centre of mass can fall outside it; a table
                # would give a voxel in the cluster, so take the member nearest to it.
                inside = np.argwhere(member)
                ijk = inside[np.argmin(np.linalg.norm((inside - ijk) * zooms, axis=1))]
        else:
            raise ValueError(f"unknown focus {focus!r}")
        foci.append((np.asarray(ijk, dtype=np.int64), float(volume[tuple(ijk)])))
    return foci, height


def study_t_field(truth_g, n_subjects, smooth_vox, rng, shape=None):
    """One study's observed t map: a genuine t statistic, not a known-variance z.

    Every bed in this repository used to build a statistic as ``(truth + noise / sqrt(n)) *
    sqrt(n)``, which is a *normal* statistic with known variance -- exactly ``d * sqrt(n)``. The
    estimator, correctly for real data, reads a reported z as a p-value-preserving image of a t
    on ``n - 1`` degrees of freedom and maps it back before converting to an effect size. In the
    far tail, where every reported peak lives, that map is strongly expansive: at ``z = 5.33``
    with ``n = 30`` it turns ``d = 0.98`` into ``g = 1.25``. Feeding a known-variance z into it
    inflates the recovered effect size by about a third, and that inflation was mistaken for a
    property of the estimator across a whole session's measurements.

    So the statistic is built the way a one-sample group analysis produces one: an effect plus
    smooth noise, over an independently drawn residual standard deviation on ``n - 1`` degrees
    of freedom.

    The variance field is made smooth by a probability integral transform of a smooth Gaussian
    field, *not* by smoothing a chi-square field. Smoothing chi-square draws averages
    independent variates, so the denominator becomes nearly constant and the t collapses back
    into a z -- which is the very confusion this function exists to avoid, and which a first
    version of it walked straight into. The transform keeps the marginal exactly
    chi-square on ``n - 1`` degrees of freedom while giving it spatial structure.
    """
    shape = tuple(shape) if shape is not None else np.asarray(truth_g).shape
    df = n_subjects - 1

    numerator = ndimage.gaussian_filter(rng.standard_normal(shape), smooth_vox)
    numerator *= 1.0 / (numerator.std() + 1e-12)
    numerator = np.asarray(truth_g) * np.sqrt(n_subjects) + numerator

    latent = ndimage.gaussian_filter(rng.standard_normal(shape), smooth_vox)
    latent *= 1.0 / (latent.std() + 1e-12)
    chi = stats.chi2.ppf(np.clip(stats.norm.cdf(latent), 1e-9, 1 - 1e-9), df)
    return numerator / np.sqrt(chi / df)


def assert_statistic_convention(stat, n_subjects, kind, threshold=3.0, tolerance=0.35):
    """Fail loudly if a bed's statistic is not on the scale its ``kind`` claims.

    The check every bed should run once, on a null field, before reading any estimate out of it.

    It compares *tail* rates rather than variances. A t on 29 degrees of freedom has variance
    1.074 against a standard normal's 1.000 -- a 7% difference that sampling noise hides, and a
    first version of this check duly passed a field whose t had collapsed into a z. The tails
    are not close: ``P(|T_29| > 3)`` is 0.0055 against ``P(|Z| > 3)`` of 0.0027, a factor of two,
    and the gap widens further out. Since reported peaks live in exactly that tail, it is also
    the discrepancy that matters.
    """
    values = np.abs(np.asarray(stat).ravel())
    observed = float(np.mean(values > threshold))
    df = n_subjects - 1
    if kind.lower().startswith("t"):
        expected = float(2.0 * stats.t.sf(threshold, df))
    else:
        expected = float(2.0 * stats.norm.sf(threshold))
    # Binomial noise on the observed rate, so a small field is not failed for being small.
    noise = 3.0 * np.sqrt(max(expected, 1e-9) * (1 - expected) / max(values.size, 1))
    if abs(observed - expected) > tolerance * expected + noise:
        raise AssertionError(
            f"Statistic declared as {kind!r} exceeds {threshold} at rate {observed:.5f} against "
            f"{expected:.5f} expected on {df} degrees of freedom. A bed that reports a "
            f"known-variance z as a t -- or a t as a z -- shifts the recovered effect size by "
            f"about a third, because the estimator maps a reported z back to a t on n - 1 "
            f"degrees of freedom before converting it."
        )
