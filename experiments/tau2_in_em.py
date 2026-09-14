"""Item 13: does re-estimating tau2 about the fitted mu change the answer?

tau2 is estimated once by kernel-weighted DerSimonian-Laird on the *naive* weighted mean, then
held fixed while the EM fits mu and pi. The naive mean is winner's-curse inflated and knows
nothing about censoring, so the residuals heterogeneity is measured from are not the residuals
about the value finally reported. Whether that matters has never been checked.

The comparison has to change one thing only. An earlier version of this script recomputed tau2
with an ad-hoc moment estimator, which collapsed it to zero nearly everywhere and so measured
"DL versus no heterogeneity at all" rather than "DL centred here versus there" -- the drift it
showed was the wrong drift. Here the *same* DerSimonian-Laird expression is used, with Q taken
about a supplied centre instead of about the a-weighted mean, and it is asserted to reproduce
the estimator in the package exactly when that centre is the a-weighted mean.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.meta.cbma.effectsize import _local_dersimonian_laird
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)


def moment_sums(contributions, n_voxels, tau2_for_weights=None):
    """The DL moment sums, per voxel, optionally with tau2 already in the weights."""
    out = {k: np.zeros(n_voxels) for k in ("w", "a", "a2", "ag", "ag2", "w2_over_s2", "n")}
    for _, cols, weights, g, var_g in contributions:
        s2 = var_g if tau2_for_weights is None else var_g + tau2_for_weights[cols]
        a = weights / s2
        out["w"] += np.bincount(cols, weights=weights, minlength=n_voxels)
        out["a"] += np.bincount(cols, weights=a, minlength=n_voxels)
        out["a2"] += np.bincount(cols, weights=a**2, minlength=n_voxels)
        out["ag"] += np.bincount(cols, weights=a * g, minlength=n_voxels)
        out["ag2"] += np.bincount(cols, weights=a * g * g, minlength=n_voxels)
        out["w2_over_s2"] += np.bincount(cols, weights=weights**2 / s2, minlength=n_voxels)
        out["n"] += np.bincount(cols, minlength=n_voxels)
    return out


def dl_about(sums, centre=None):
    """DerSimonian-Laird with Q taken about ``centre``, or the a-weighted mean when None."""
    sum_a = sums["a"]
    tau2 = np.zeros_like(sum_a)
    usable = (sums["n"] >= 2) & (sum_a > 0)
    if not np.any(usable):
        return tau2
    a_u = sum_a[usable]
    if centre is None:
        q = sums["ag2"][usable] - sums["ag"][usable] ** 2 / a_u
    else:
        c = centre[usable]
        q = sums["ag2"][usable] - 2.0 * c * sums["ag"][usable] + c**2 * a_u
    expected = sums["w"][usable] - sums["w2_over_s2"][usable] / a_u
    scale = a_u - sums["a2"][usable] / a_u
    out = np.zeros_like(a_u)
    positive = scale > 0
    out[positive] = (q[positive] - expected[positive]) / scale[positive]
    tau2[usable] = np.maximum(out, 0.0)
    return tau2


def focus_index(masker):
    ijk = mm2vox(np.array([TRUTH]), masker.mask_img.affine)[0]
    mask = np.asarray(masker.mask_img.dataobj).astype(bool)
    lookup = np.full(mask.shape, -1, dtype=np.int64)
    lookup[mask] = np.arange(mask.sum())
    return int(flat := lookup[tuple(ijk)]) if lookup[tuple(ijk)] >= 0 else -1


print("Same DL expression, Q re-centred on the fitted mu, three extra rounds.\n")
for true_tau in (0.0, 0.15, 0.35):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH],
        effect_sizes=0.8,
        n_studies=30,
        sample_size=(20, 40),
        tau=true_tau,
        prevalence=1.0,
        seed=1,
        n_noise_foci=2,
        noise_extent=40.0,
    )
    index = focus_index(studyset.masker)

    rounds = []
    centre = {"g": None}
    original_pool = es.CBES._pool

    def pool_recentred(self, table, image_studies=None):
        fit = original_pool(self, table, image_studies)
        if centre["g"] is None:
            return fit  # round 0 is the estimator as shipped
        contributions, n_voxels = fit["contributions"], fit["n_voxels"]
        sums = moment_sums(contributions, n_voxels)
        # Oracle: with the a-weighted mean as the centre this must be the shipped estimator.
        naive = np.divide(sums["ag"], sums["a"], out=np.zeros(n_voxels), where=sums["a"] > 0)
        assert np.allclose(
            dl_about(sums, None),
            _local_dersimonian_laird(
                sums["w"], sums["a"], sums["a2"], sums["ag"], sums["ag2"],
                sums["w2_over_s2"], sums["n"],
            ),
        ), "re-implementation does not reproduce the packaged estimator"
        assert np.allclose(dl_about(sums, naive), dl_about(sums, None), atol=1e-10)

        tau2 = dl_about(sums, centre["g"])
        fit["tau2"] = tau2
        numer = np.zeros(n_voxels)
        denom = np.zeros(n_voxels)
        var_numer = np.zeros(n_voxels)
        for _, cols, weights, g, var_g in contributions:
            total_var = var_g + tau2[cols]
            w = weights / total_var
            numer += np.bincount(cols, weights=w * g, minlength=n_voxels)
            denom += np.bincount(cols, weights=w, minlength=n_voxels)
            var_numer += np.bincount(cols, weights=weights**2 / total_var, minlength=n_voxels)
        covered = denom > 0
        fit["g"] = np.zeros(n_voxels)
        fit["g"][covered] = numer[covered] / denom[covered]
        fit["se"] = np.full(n_voxels, np.inf)
        fit["se"][covered] = np.sqrt(var_numer[covered]) / denom[covered]
        fit["covered"] = covered
        return fit

    es.CBES._pool = pool_recentred
    try:
        for r in range(4):
            result = CBES(fwhm=10.0, null_method="none", peak_bias="per-study").fit(studyset)
            g = result.get_map("g", return_type="array").ravel()
            tau2 = result.get_map("tau2", return_type="array").ravel()
            rounds.append((g, tau2))
            centre["g"] = g
    finally:
        es.CBES._pool = original_pool

    covered = rounds[0][0] != 0
    print(f"true tau = {true_tau:.2f} (tau^2 = {true_tau**2:.4f}); truth at focus g = 0.8")
    for r, (g, tau2) in enumerate(rounds):
        shift = "-" if r == 0 else f"{np.max(np.abs(g[covered] - rounds[r-1][0][covered])):.4f}"
        print(f"  round {r}: g at focus {g[index]:6.3f}   tau2 at focus {tau2[index]:8.5f}   "
              f"tau2 median {np.median(tau2[covered]):8.5f}   max |dg| {shift}", flush=True)
    print()
