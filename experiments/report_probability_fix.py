"""Does correcting `P(report)` fix the magnitude? Test the idea before engineering it.

The diagnosis: the likelihood gives a reported voxel `P(|g| >= c | mu)`, but a paper reports a
voxel only if it cleared `c` **and** was a local maximum -- a strictly smaller event. Measured,
the model's reporting rate exceeds the observed one by 1.00, 1.78, 1.08 and 1.13 at true g of
0.2, 0.4, 0.6 and 0.8. This is implicated in four separate deficits: the 18% shortfall at the 0.4
focus, `se/sd` of 1.55 to 3.74, the HCP magnitude at 0.63, and the prevalence absorbing
heterogeneity (0.924 to 0.617 as tau goes 0 to 0.6 at a true prevalence of 1).

The fix is not a reweighting -- that was tried and made everything worse. It is that **both limbs
must be complementary probabilities of the same event**. Writing `E = P(|g| >= c | mu)` for the
exceedance and `q` for the chance that a study which exceeded at this voxel actually *names* it:

    reported here      q E            instead of   E
    silent here        1 - q E        instead of   1 - E

Both still sum to 1, so the model stays a proper likelihood, but it now describes reporting
rather than exceeding. And the arithmetic says exactly where the benefit comes from:

    report limb   score = (qE)'/(qE) = E'/E          -- q cancels, unchanged
    silent limb   score = -qE'/(1 - qE)              -- scaled by roughly q

So `q` leaves the report limb alone and **weakens the silence pull by about q**, which is
precisely the over-shrinkage that was diagnosed. Nothing here needs a field smoothness.

This fixes `q` at a plausible constant to see whether the deficit closes at all. If it does, `q`
is worth estimating properly (a global scalar profiled out of the likelihood, or a per-study
moment condition on the focus count); if it does not, the diagnosis is wrong.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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
N_SEEDS = int(os.environ.get("NSEEDS", 8))
QS = [float(x) for x in os.environ.get("QS", "1.0,0.75,0.56,0.4,0.25").split(",")]

LOOKUP = MASKER.transform(nib.Nifti1Image(
    np.arange(np.prod(SHAPE)).reshape(SHAPE).astype(np.float32), AFFINE)).ravel()
POSITIONS = [
    int(np.flatnonzero(LOOKUP == np.ravel_multi_index(
        tuple(np.rint(np.linalg.solve(AFFINE, np.array([*f, 1.0]))[:3]).astype(int)), SHAPE))[0])
    for f in FOCI
]
_ORIGINAL = es._censoring_terms


def corrected(q):
    """`_censoring_terms` where both limbs describe *reporting* rather than exceeding."""

    def wrapped(mu, cutoff_scaled, twice_cutoff_scaled, inv_sigma, inv_sigma_sq, sign):
        raw = _ORIGINAL(mu, cutoff_scaled, twice_cutoff_scaled, inv_sigma, inv_sigma_sq,
                        np.ones_like(sign))
        # `raw` is the silence limb: prob = P(|g| < c), score = P'/P, d2 = P''/P.
        silent = raw["prob"]
        first = raw["score"] * silent          # d/dmu P(|g| < c)
        second = raw["d2_over_prob"] * silent  # d2/dmu2 P(|g| < c)
        exceed = 1.0 - silent
        e1, e2 = -first, -second               # derivatives of the exceedance

        reported = sign < 0
        prob = np.where(reported, q * exceed, 1.0 - q * exceed)
        prob = np.clip(prob, es._PROBABILITY_FLOOR, None)
        score = np.where(reported, q * e1, -q * e1) / prob
        d2_over_prob = np.where(reported, q * e2, -q * e2) / prob
        return {"prob": prob, "score": score, "d2_over_prob": d2_over_prob}

    return wrapped


if __name__ == "__main__":
    built = [create_effect_size_coordinate_studyset(
        FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1, seed=seed,
        simulate_field=True, n_image_studies=2, image_dir=tempfile.mkdtemp(),
        noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB) for seed in range(N_SEEDS)]

    print(f"four foci at {EFFECTS}, {N_SEEDS} collections, 20 studies, 2 images")
    print("q = 1.0 is the shipped model; the measured over-statement was about 1.78 at the "
          "0.4 focus, so q ~ 0.56 there.\n")
    print(f"{'q':>6} " + " ".join(f"{'g@' + str(e):>14}" for e in EFFECTS) +
          f" {'mean |err|':>11}")
    for q in QS:
        es._censoring_terms = _ORIGINAL if q == 1.0 else corrected(q)
        peaks = []
        for ss in built:
            res = CBES(mask=MASK, null_method="none", threshold="reporting_threshold").fit(ss)
            g = np.abs(res.get_map("g", return_type="array").ravel())
            peaks.append([g[p] for p in POSITIONS])
        mean = np.array(peaks).mean(axis=0)
        err = np.mean(np.abs(100 * (mean - np.array(EFFECTS)) / np.array(EFFECTS)))
        cells = [f"{v:.3f} ({100 * (v - e) / e:+.0f}%)" for v, e in zip(mean, EFFECTS)]
        print(f"{q:6.2f} " + " ".join(f"{c:>14}" for c in cells) + f" {err:10.1f}%")
    es._censoring_terms = _ORIGINAL
