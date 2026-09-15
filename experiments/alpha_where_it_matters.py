"""Re-run the report-limb probe in the regime the first one could not reach.

The first alpha sweep scored a single-focus bed in three strata, and the middle stratum came back
**exactly invariant** to alpha -- 0.0737 at every value from 0.1 to 5.0. That is not a null
result, it is a dead test: in that bed the middle stratum is the blob skirt, where almost nothing
is reported, so there are no `sign = -1` pairs for alpha to act on. The probe never exercised the
limb it was probing.

Where it would matter is the *window of detectability*: truth near the cut, so the chance of
reporting responds strongly to mu. The four-foci bed spans it -- at cutoffs near 0.6 g the
reporting probability runs about 1.5%, 14%, 50% and 86% at true g of 0.2, 0.4, 0.6 and 0.8 -- and
the measured per-focus error is non-monotone in exactly the suspicious way:

    truth      0.2     0.4     0.6     0.8
    g err%    -6.0   -13.8    +0.8    -0.3

The worst error sits where the likelihood is most sensitive to the reporting model, which is what
a mis-specified reporting model would look like. The known mis-specification: a reported peak is
`|g| >= c` **and** a local maximum, a strictly smaller event, so P(report | mu) is overstated.

So: sweep alpha again, scoring each focus separately. If the report limb is the cause, some
alpha < 1 should lift the 0.4 focus specifically without wrecking 0.6 and 0.8.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
import tempfile
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import effectsize as es
from nimare.meta.cbma.effectsize import CBES

ZOOMS, EXTENT, BLOB = 4.0, 48.0, 10.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
FOCI = [(-28, -28, 0), (28, -28, 0), (-28, 28, 0), (28, 28, 0)]
EFFECTS = [0.2, 0.4, 0.6, 0.8]
N_SEEDS = int(os.environ.get("NSEEDS", 6))
ALPHAS = [float(a) for a in os.environ.get("ALPHAS", "1.0,0.7,0.5,0.3").split(",")]

LOOKUP = MASKER.transform(nib.Nifti1Image(
    np.arange(np.prod(SHAPE)).reshape(SHAPE).astype(np.float32), AFFINE)).ravel()
POSITIONS = []
for focus in FOCI:
    ijk = np.rint(np.linalg.solve(AFFINE, np.array([*focus, 1.0]))[:3]).astype(int)
    POSITIONS.append(int(np.flatnonzero(LOOKUP == np.ravel_multi_index(tuple(ijk), SHAPE))[0]))

_ORIGINAL = es._censoring_terms


def flattened(alpha):
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


if __name__ == "__main__":
    built = [create_effect_size_coordinate_studyset(
        FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1, seed=seed,
        simulate_field=True, n_image_studies=2, image_dir=tempfile.mkdtemp(),
        noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB) for seed in range(N_SEEDS)]

    # How much of the report limb is even present at each focus, so a null result can be told
    # apart from a dead test this time.
    est = CBES(mask=MASK, null_method="none", threshold="reporting_threshold")
    est.fit(built[0])
    table, roster = est._focus_table_, list(est._sample_sizes_.index)
    active = np.arange(int(np.asarray(MASK.dataobj).astype(bool).sum()))
    col, pos, sign = est._indicator_entries(
        table, roster, active, active.size, image_ids=tuple(est._image_studies_))
    ind = np.zeros((len(roster), active.size))
    ind[pos, col] = sign
    print(f"four foci at {EFFECTS}, {N_SEEDS} collections, 20 studies, 2 images")
    print("indicator pairs at each focus (seed 0):")
    for effect, p in zip(EFFECTS, POSITIONS):
        print(f"   truth {effect}: {int((ind[:, p] == -1).sum())} reported, "
              f"{int((ind[:, p] == 1).sum())} silent, {int((ind[:, p] == 0).sum())} mute")

    print(f"\n{'alpha':>7} " + " ".join(f"{'g@' + str(e):>14}" for e in EFFECTS))
    for alpha in ALPHAS:
        es._censoring_terms = _ORIGINAL if alpha == 1.0 else flattened(alpha)
        peaks = []
        for ss in built:
            res = CBES(mask=MASK, null_method="none", threshold="reporting_threshold").fit(ss)
            g = np.abs(res.get_map("g", return_type="array").ravel())
            peaks.append([g[p] for p in POSITIONS])
        mean = np.array(peaks).mean(axis=0)
        cells = [f"{v:.3f} ({100 * (v - e) / e:+.0f}%)" for v, e in zip(mean, EFFECTS)]
        print(f"{alpha:7.1f} " + " ".join(f"{c:>14}" for c in cells))
    es._censoring_terms = _ORIGINAL
