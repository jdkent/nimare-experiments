"""Why is `g` 18% low at the 0.4 focus, and only there?

Arithmetic first, because it points somewhere specific. At the 0.4 focus, cutoffs near 0.6 g and
sigma = sqrt(1/30) = 0.183, the model's reporting probability is

    P(report | mu = 0.4) = 1 - Phi((0.6 - 0.4) / 0.183) = 0.138

so 17 coordinate studies should produce about 2.3 reports. The representative collection showed
**1**. Inverting the model on that count gives

    1 - Phi((0.6 - mu) / 0.183) = 1/17   =>   mu = 0.315

against the 0.329 measured at 24 seeds. So the estimate is not mis-fitting: it is faithfully
reporting a count that came in below what the truth predicts, and the indicator channel -- 17
observations against 2 image values -- dominates.

Two candidate reasons the count runs low *systematically* rather than just noisily:

  displacement  a study whose peak lands one voxel from the focus is not a report *at* the focus,
                so the observed count understates the true detection rate. At 0.4 the count is
                small, so each lost report costs a lot of mu. Measured on a single-focus bed,
                widening the report radius was catastrophic -- but that bed had no voxel in the
                window, so it could not see this. Re-test on the four-foci bed, per focus.

  Jensen        mu-hat solves observed count = expected count, and that map is strongly
                non-linear near the cut, so E[f^-1(count)] != f^-1(E[count]) even with an
                unbiased count. Diagnosed by comparing the fitted mu against the mu the
                *observed* count implies, focus by focus: if they agree, the inversion is
                faithful and the count is what is biased.

Reports both, per focus, so the two can be told apart rather than argued about.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
import tempfile
import numpy as np
import nibabel as nib
from scipy.stats import norm
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES, reporting_cutoff_to_g

ZOOMS, EXTENT, BLOB = 4.0, 48.0, 10.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
FOCI = [(-28, -28, 0), (28, -28, 0), (-28, 28, 0), (28, 28, 0)]
EFFECTS = [0.2, 0.4, 0.6, 0.8]
N_SEEDS = int(os.environ.get("NSEEDS", 16))
RADII = [float(r) for r in os.environ.get("RADII", "0,4,8").split(",")]

LOOKUP = MASKER.transform(nib.Nifti1Image(
    np.arange(np.prod(SHAPE)).reshape(SHAPE).astype(np.float32), AFFINE)).ravel()
POSITIONS = [
    int(np.flatnonzero(LOOKUP == np.ravel_multi_index(
        tuple(np.rint(np.linalg.solve(AFFINE, np.array([*f, 1.0]))[:3]).astype(int)), SHAPE))[0])
    for f in FOCI
]


def mu_from_count(reported, total, cutoff_g, sigma):
    """Invert the model's reporting probability on an observed count, for the diagnosis."""
    rate = np.clip(reported / max(total, 1), 1e-6, 1 - 1e-6)
    # P(|g| >= c | mu) = 1 - (Phi((c-mu)/s) - Phi((-c-mu)/s)); the lower tail is negligible here.
    return float(cutoff_g - sigma * norm.isf(rate))


if __name__ == "__main__":
    print(f"four foci at {EFFECTS}, {N_SEEDS} collections, 20 studies, 2 images\n")
    built = [create_effect_size_coordinate_studyset(
        FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1, seed=seed,
        simulate_field=True, n_image_studies=2, image_dir=tempfile.mkdtemp(),
        noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB) for seed in range(N_SEEDS)]

    # Diagnosis: is the fit faithful to the count it was given?
    print("Is the fit faithful to the observed count? (radius 0, the shipped behaviour)")
    print(f"{'truth':>6} {'reported':>9} {'silent':>7} {'fitted mu':>10} "
          f"{'mu from count':>14} {'true rate':>10}")
    rows = {j: {"rep": [], "sil": [], "mu": [], "inv": [], "true_rate": []}
            for j in range(len(FOCI))}
    for ss in built:
        est = CBES(mask=MASK, null_method="none", threshold="reporting_threshold")
        res = est.fit(ss)
        g = np.abs(res.get_map("g", return_type="array").ravel())
        roster = list(est._sample_sizes_.index)
        active = np.arange(int(np.asarray(MASK.dataobj).astype(bool).sum()))
        col, pos, sign = est._indicator_entries(
            est._focus_table_, roster, active, active.size,
            image_ids=tuple(est._image_studies_))
        ind = np.zeros((len(roster), active.size))
        ind[pos, col] = sign
        cutoff_g = float(np.median(est._thresholds_.values))
        sigma = float(np.median(np.sqrt(1.0 / est._sample_sizes_.values)))
        for j, (effect, p) in enumerate(zip(EFFECTS, POSITIONS)):
            reported = int((ind[:, p] == -1).sum())
            silent = int((ind[:, p] == 1).sum())
            rows[j]["rep"].append(reported)
            rows[j]["sil"].append(silent)
            rows[j]["mu"].append(float(g[p]))
            rows[j]["inv"].append(mu_from_count(reported, reported + silent, cutoff_g, sigma))
            rows[j]["true_rate"].append(
                1.0 - (norm.cdf((cutoff_g - effect) / sigma) - norm.cdf((-cutoff_g - effect) / sigma))
            )
    print(f"\n{'truth':>6} {'reports':>9} {'total':>7} {'obs rate':>10} {'model rate':>11} "
          f"{'q = obs/model':>15}")
    for j, effect in enumerate(EFFECTS):
        r = rows[j]
        # The ratio of two rates built from a handful of reports. Its uncertainty is dominated
        # by the Poisson count in the numerator, so report that rather than a bare point
        # estimate -- a "1.78 over-statement" resting on 14 reports is not a measurement.
        total_reports = float(np.sum(r["rep"]))
        total_pairs = float(np.sum(r["rep"]) + np.sum(r["sil"]))
        obs = total_reports / max(total_pairs, 1.0)
        model = float(np.mean(r["true_rate"]))
        rel = 1.0 / np.sqrt(max(total_reports, 1e-9)) if total_reports else np.inf
        q = obs / model if model > 0 else np.nan
        lo, hi = q * (1 - 1.96 * rel), q * (1 + 1.96 * rel)
        print(f"{effect:6.1f} {total_reports:9.0f} {total_pairs:7.0f} {obs:10.4f} "
              f"{model:11.4f} {q:8.2f} [{max(lo, 0):.2f}, {hi:.2f}]")

    print(f"\nDoes widening the report radius help *at the window*? (per focus)")
    print(f"{'radius':>8} " + " ".join(f"{'g@' + str(e):>14}" for e in EFFECTS))
    for radius in RADII:
        peaks = []
        for ss in built:
            est = CBES(mask=MASK, null_method="none", threshold="reporting_threshold")
            est._report_radius = radius
            res = est.fit(ss)
            g = np.abs(res.get_map("g", return_type="array").ravel())
            peaks.append([g[p] for p in POSITIONS])
        mean = np.array(peaks).mean(axis=0)
        cells = [f"{v:.3f} ({100 * (v - e) / e:+.0f}%)" for v, e in zip(mean, EFFECTS)]
        label = "named" if radius == 0 else f"{radius:g} mm"
        print(f"{label:>8} " + " ".join(f"{c:>14}" for c in cells))
