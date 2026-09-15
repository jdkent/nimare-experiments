"""Is a report over-credited, because the model treats it as any exceedance?

Two results want the same explanation. Clamping the assumed threshold *below* the truth beat
the true thresholds (rmse 0.106 against 0.128 where the truth is near zero), and `g` comes back
biased *downward* where the effect is (-0.082 at two images) while biased upward in the quiet
strata. A cut placed too low is doing something a correct cut is not.

Hypothesis. The model gives a reported voxel the probability

    P(report | mu) = P(|g| >= c | mu)

but a paper reports a voxel only if it cleared c **and** was a local maximum -- a strictly
smaller event, and one whose conditional probability *falls* with mu, since a larger effect makes
the suprathreshold region bigger and any given voxel less likely to be its maximum. So the true
P(report | mu) rises more slowly in mu than the plain exceedance does, the report limb pushes mu
too hard, and the silence limb has to pull it back -- which lands below the truth.

The full random-field survival needs a smoothness the model does not have. This probes the shape
instead, with one parameter: raise the exceedance to a power alpha in (0, 1], which flattens its
mu-sensitivity without changing its range. If the hypothesis holds, some alpha < 1 should beat
alpha = 1 where the effect is, and should shrink the gap between an assumed cut and the true one.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
import collections
import tempfile
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import effectsize as es
from nimare.meta.cbma.effectsize import CBES

SHAPE, ZOOMS, EXTENT, BLOB = (21, 21, 21), 4.0, 40.0, 10.0
AFFINE = np.array([[ZOOMS, 0, 0, -EXTENT], [0, ZOOMS, 0, -EXTENT],
                   [0, 0, ZOOMS, -EXTENT], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
TRUE_G = 0.5
N_SEEDS = int(os.environ.get("NSEEDS", 8))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - EXTENT
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
TRUTH = MASKER.transform(nib.Nifti1Image(
    (TRUE_G * np.exp(-(grid ** 2).sum(-1) / (2 * sd ** 2))).astype(np.float32), AFFINE)).ravel()
STRATA = [("quiet", TRUTH < 0.05), ("middle", (TRUTH >= 0.05) & (TRUTH < 0.25)),
          ("effect", TRUTH >= 0.25)]

_ORIGINAL = es._censoring_terms


def flattened(alpha):
    """``_censoring_terms`` with the report limb's probability raised to ``alpha``.

    Only the ``sign = -1`` pairs are touched. ``p -> p**a`` gives ``log p -> a log p``, so the
    score and the normalised second derivative both scale by ``a`` while ``prob`` itself is
    raised -- which is exactly "the report is worth less than a full exceedance" without
    changing what a probability of 0 or 1 means.
    """

    def wrapped(mu, cutoff_scaled, twice_cutoff_scaled, inv_sigma, inv_sigma_sq, sign):
        out = _ORIGINAL(mu, cutoff_scaled, twice_cutoff_scaled, inv_sigma, inv_sigma_sq, sign)
        reported = sign < 0
        if np.any(reported):
            out["prob"] = np.where(reported, out["prob"] ** alpha, out["prob"])
            out["score"] = np.where(reported, out["score"] * alpha, out["score"])
            out["d2_over_prob"] = np.where(
                reported, out["d2_over_prob"] * alpha, out["d2_over_prob"]
            )
        return out

    return wrapped


def score(maps):
    return {name: np.mean([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in maps])
            for name, sel in STRATA}, {
        name: np.mean([float(np.mean(m[sel] - TRUTH[sel])) for m in maps])
        for name, sel in STRATA}


if __name__ == "__main__":
    collections_built = []
    for seed in range(N_SEEDS):
        tmp = tempfile.mkdtemp()
        collections_built.append(create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=TRUE_G, n_studies=20, sample_size=(20, 40), tau=0.1,
            seed=seed, simulate_field=True, n_image_studies=2, image_dir=tmp,
            noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB))

    rows = collections.OrderedDict()
    for alpha in (1.0, 1.5, 2.0, 3.0, 5.0):
        es._censoring_terms = flattened(alpha) if alpha != 1.0 else _ORIGINAL
        # Rebind on the dataclass too: ``_IndicatorPairs.censoring`` closed over the module name
        # at call time, so patching the module attribute is enough -- but assert it.
        maps = []
        for ss in collections_built:
            res = CBES(mask=MASK, null_method="none",
                       threshold="reporting_threshold").fit(ss)
            maps.append(np.abs(res.get_map("g", return_type="array").ravel()))
        rows[f"alpha = {alpha:.1f}"] = score(maps)
    es._censoring_terms = _ORIGINAL

    print(f"true g = {TRUE_G} at the focus, {N_SEEDS} collections, 20 studies, 2 images\n")
    print(f"{'report limb':>14} " + " ".join(f"{'rmse ' + n:>12}" for n, _ in STRATA)
          + " " + " ".join(f"{'bias ' + n:>12}" for n, _ in STRATA))
    for label, (rmse, bias) in rows.items():
        print(f"{label:>14} " + " ".join(f"{rmse[n]:12.4f}" for n, _ in STRATA)
              + " " + " ".join(f"{bias[n]:+12.4f}" for n, _ in STRATA))
