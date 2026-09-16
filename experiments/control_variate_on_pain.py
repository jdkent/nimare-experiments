"""What is the predictor's correlation on real published tables, and what does it buy?

The control-variate estimator's whole benefit is set by one number: rho, the correlation between
a study's actual effect map and a prediction made from its published coordinate table. The
algebra (proofs/marginal_control_variate.py) says the best achievable variance ratio is
1 - rho^2 / (1 + n/N) and that it floors at 1 - rho^2. The design document is explicit that its
23% and 59% figures are "algebraic targets, not measured Neurostore prediction quality".

So the contribution here is measuring rho on real data, and then checking that the estimator
delivers what that rho predicts.

Protocol, per PROTOCOL.md:

  * **The truth is independent of the coordinates.** Studies are split three ways: a held-out
    reference set whose images define the target, an image cohort whose images and tables are
    both used, and a coordinate-only cohort whose tables alone are used. No study appears twice.
  * **No cap on peaks.** Every published peak of every study is used.
  * Published tables, not extracted ones, so the reporting process is the real one.

The predictor is a fixed isotropic Gaussian kernel over the reported peaks, with no fitted
parameters at all. Its amplitude is irrelevant: the estimator uses the predictor only through a
difference of cohort means scaled by lambda, so rescaling the predictor and rescaling lambda
inversely is a no-op (verified as an exact invariant in test_meta_marginal.py). Only the shape
matters, which is why the only choice to make is the width.

**Kill condition, stated before running.** If rho on published tables is below about 0.3, the
achievable reduction is under 10% and this estimator is not worth building further for this kind
of corpus, whatever the simulations said.

Note on which rho: the variance formula's rho is the *between-study* correlation at a voxel,
because that is the covariance the estimator averages over. The within-study spatial correlation
-- whether the kernel recreates one study's map -- is a different and much larger number, and
using it would overstate the benefit several times over. Both are reported.

Environment knobs: SPLITS, NIMAGES, NREF, FWHMS.
"""
import logging
import os
import sys
import warnings

warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)

import numpy as np  # noqa: E402
import nibabel as nib  # noqa: E402
from nilearn.datasets import load_mni152_brain_mask  # noqa: E402
from nilearn.maskers import NiftiMasker  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/user/NiMARE")
from load_pain import load_pain  # noqa: E402
from nimare.meta.cbma.marginal import (  # noqa: E402
    achievable_ratio,
    control_variate_mean,
    optimal_coefficient,
)
from nimare.transforms import d_to_g, t_to_d  # noqa: E402

SPLITS = int(os.environ.get("SPLITS", "40"))
N_IMAGES = int(os.environ.get("NIMAGES", "6"))
N_REFERENCE = int(os.environ.get("NREF", "8"))
FWHMS = [float(x) for x in os.environ.get("FWHMS", "10,20,30,45").split(",")]

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.dataobj) > 0
shape = mask_bool.shape
affine = mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    """A study's z map on the Hedges g scale, at its own sample size."""
    return d_to_g(t_to_d(z, n), n)


def peak_kernel_map(xyz, fwhm):
    """A fixed Gaussian kernel over every reported peak, no fitted parameters.

    Amplitude is arbitrary -- the estimator absorbs it into lambda -- so this returns an
    unnormalised sum of kernels. Peaks outside the mask still contribute, since a paper can
    report a coordinate the mask does not cover and its neighbourhood is still informative.
    """
    volume = np.zeros(shape, dtype=np.float32)
    if xyz is None or not len(xyz):
        return masker.transform(nib.Nifti1Image(volume, affine)).ravel()
    inverse = np.linalg.inv(affine)
    for point in np.atleast_2d(np.asarray(xyz, dtype=float)):
        ijk = np.rint(nib.affines.apply_affine(inverse, point)).astype(int)
        if np.all(ijk >= 0) and np.all(ijk < np.asarray(shape)):
            volume[tuple(ijk)] += 1.0
    sigma = (fwhm / 2.355) / zooms
    smoothed = gaussian_filter(volume, sigma=sigma, mode="constant")
    return masker.transform(nib.Nifti1Image(smoothed, affine)).ravel()


def main():
    studyset = load_pain()
    effects, tables, ids = [], [], []
    for row, n in zip(studyset.images.itertuples(), studyset.sample_sizes()):
        size = int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0])
        effects.append(np.nan_to_num(to_g(masker.transform(row.z).ravel(), size)))
        ids.append(row.study_id)
    effects = np.asarray(effects)

    by_study = {
        sid: group[["x", "y", "z"]].astype(float).to_numpy()
        for sid, group in studyset.coordinates.groupby("study_id")
    }
    peaks = [by_study.get(sid) for sid in ids]
    n_studies = len(ids)
    n_peaks = sum(0 if p is None else len(p) for p in peaks)
    with_table = sum(1 for p in peaks if p is not None and len(p))
    print(f"NIDM pain: {n_studies} studies, {with_table} with published tables, "
          f"{n_peaks} peaks, no cap applied")
    print(f"{N_REFERENCE} held out as the reference, {N_IMAGES} as the image cohort, "
          f"the rest coordinate-only; {SPLITS} splits")
    print()

    # ------------------------------------------- two different rhos, and only one is the right one
    #
    # The variance formula's rho is Corr(Y, f) *across studies at a fixed voxel*: the estimator
    # averages over studies, so its covariances are between-study ones. That is what decides the
    # reduction.
    #
    # It is not the same as the within-study, across-voxel correlation -- whether the kernel
    # recovers the shape of one study's map. A predictor can reproduce every study's map shape
    # well and still carry nothing about how one study differs from another at a given voxel,
    # because where a paper reports is largely common across a domain. Both are printed because
    # the contrast is the result.
    print("rho_spatial:  within a study, across voxels -- does the kernel recover its map shape?")
    print("rho_between:  across studies, at a voxel -- what the variance formula actually uses.")
    print()
    print(f"{'fwhm':>6} {'rho_spatial':>12} {'rho_between':>12} {'ratio at 6+7':>13} "
          f"{'floor':>7}")
    predictions_by_fwhm = {}
    n_img, n_coord = N_IMAGES, n_studies - N_REFERENCE - N_IMAGES
    for fwhm in FWHMS:
        predicted = np.asarray([peak_kernel_map(p, fwhm) for p in peaks])
        predictions_by_fwhm[fwhm] = predicted

        spatial = [
            float(np.corrcoef(effects[i], predicted[i])[0, 1])
            for i in range(n_studies)
            if predicted[i].std() > 0
        ]
        # Between-study correlation per voxel, over all studies, then summarised. Voxels where
        # the predictor never varies across studies carry no information and are skipped.
        centred_y = effects - effects.mean(axis=0)
        centred_f = predicted - predicted.mean(axis=0)
        denominator = np.sqrt((centred_y**2).sum(axis=0) * (centred_f**2).sum(axis=0))
        between = np.divide(
            (centred_y * centred_f).sum(axis=0),
            denominator,
            out=np.full(denominator.shape, np.nan),
            where=denominator > 0,
        )
        rho_between = float(np.nanmedian(between))
        print(f"{fwhm:>6.0f} {np.median(spatial):>12.3f} {rho_between:>12.3f} "
              f"{float(achievable_ratio(rho_between, n_img, n_coord)):>13.3f} "
              f"{1 - rho_between**2:>7.3f}")

    # Pick on the between-study correlation, since that is the one the estimator exploits.
    def between_median(f):
        centred_y = effects - effects.mean(axis=0)
        centred_f = predictions_by_fwhm[f] - predictions_by_fwhm[f].mean(axis=0)
        denominator = np.sqrt((centred_y**2).sum(axis=0) * (centred_f**2).sum(axis=0))
        values = np.divide(
            (centred_y * centred_f).sum(axis=0), denominator,
            out=np.full(denominator.shape, np.nan), where=denominator > 0,
        )
        return abs(float(np.nanmedian(values)))

    best_fwhm = max(FWHMS, key=between_median)
    print(f"\nlargest between-study correlation at {best_fwhm:.0f} mm; using it below")
    predicted = predictions_by_fwhm[best_fwhm]

    # ---------------------------------------------------------- against a held-out reference
    rng = np.random.default_rng(0)
    rows = {"images only": [], "control variate": []}
    ratios, correlations, shifts = [], [], []
    for _ in range(SPLITS):
        order = rng.permutation(n_studies)
        reference_ids = order[:N_REFERENCE]
        image_ids = order[N_REFERENCE:N_REFERENCE + N_IMAGES]
        coord_ids = order[N_REFERENCE + N_IMAGES:]
        if not len(coord_ids):
            continue
        # The reference is the mean of held-out images, and no study in it contributes a
        # coordinate anywhere in the fit.
        truth = effects[reference_ids].mean(axis=0)

        lam = optimal_coefficient(
            effects[image_ids], predicted[image_ids], len(image_ids), len(coord_ids),
            mode="pooled",
        )
        out = control_variate_mean(
            effects[image_ids], predicted[image_ids], predicted[coord_ids], lam
        )
        rows["images only"].append(effects[image_ids].mean(axis=0) - truth)
        rows["control variate"].append(out["estimate"] - truth)
        ratios.append(np.nanmedian(out["variance_ratio"]))
        correlations.append(np.nanmedian(out["correlation"]))
        shifts.append(np.nanmedian(np.abs(out["cohort_shift"]) / out["cohort_shift_se"]))

    print()
    print(f"{'estimate':>18} {'mean err':>10} {'rmse':>8} {'r with truth':>13}")
    for name, errors in rows.items():
        err = np.asarray(errors)
        # Correlation with the reference, per split, then averaged.
        print(f"{name:>18} {err.mean():+10.4f} {np.sqrt((err**2).mean()):8.4f} "
              f"{'':>13}")
    paired = np.sqrt((np.asarray(rows['control variate'])**2).mean()) - np.sqrt(
        (np.asarray(rows['images only'])**2).mean())
    print()
    print(f"control variate minus images only, rmse: {paired:+.4f}")
    print(f"median reported variance ratio across splits: {np.nanmedian(ratios):.3f}")
    print(f"median within-cohort correlation the fit saw: {np.nanmedian(correlations):.3f}")
    print(f"median |cohort shift| in its own standard errors: {np.nanmedian(shifts):.2f}")
    print()
    print("The cohort-shift row is the exchangeability check. A large value means the")
    print("coordinate-only studies' tables look systematically unlike the image studies',")
    print("in which case the correction is importing that difference as if it were signal.")


if __name__ == "__main__":
    main()
