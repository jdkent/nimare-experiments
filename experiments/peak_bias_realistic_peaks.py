"""Finding 7 again, with a realistic number of reported peaks per study.

The field simulator reports every local maximum that clears the threshold, which comes to
28-29 per study. A real coordinate table has 3 to 20, and the ones it keeps are the *strongest*
-- a paper reports its top peaks. That matters twice over for a peak-height correction:

  * top-k selection is more selective than "everything above u", so capping should *increase*
    the winner's curse rather than leave it alone;
  * threshold="study-min" infers each study's threshold from its smallest reported peak, and
    with 28 peaks that minimum sits almost at the true cut while with 4 the order-statistic
    correction is doing real work.

So the uncapped measurement is the easy end of the range, not a neutral one. Capped
strongest-first here, at 5 and 15 peaks, against the uncapped 28 for reference.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
TRUE_G = 0.8
N_SIMS = 100

shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
MASK = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
AT = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))


def cap_peaks(studyset, max_peaks):
    """Rebuild the studyset keeping each analysis's strongest ``max_peaks`` foci.

    Strongest-first, because a paper reporting only some of its peaks reports its top ones.
    Keyed on ``(study_id, contrast_id)``: the coordinate table's ``id`` is those two
    concatenated, which does not match ``analysis.id``.
    """
    if max_peaks is None:
        return studyset
    coords = studyset.coordinates.copy()
    coords["_abs"] = np.abs(coords["z_stat"].astype(float))
    keep = (
        coords.sort_values("_abs", ascending=False)
        .groupby(["study_id", "contrast_id"], sort=False)
        .head(max_peaks)
    )
    meta = {
        (str(study.id), str(analysis.id)): dict(analysis.metadata)
        for study in studyset.studies
        for analysis in study.analyses
    }
    studies = []
    for (study_id, contrast_id), sub in keep.groupby(["study_id", "contrast_id"], sort=False):
        info = meta[(str(study_id), str(contrast_id))]
        studies.append({
            "id": str(study_id), "name": str(study_id), "metadata": info,
            "analyses": [{
                "id": str(contrast_id), "name": "1", "metadata": info,
                "points": [{"space": "MNI",
                            "coordinates": [float(r.x), float(r.y), float(r.z)],
                            "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                           for r in sub.itertuples()],
            }],
        })
    return Studyset({"id": "capped", "name": "capped", "studies": studies})


def one(seed, max_peaks, peak_bias):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=TRUE_G, n_studies=24, sample_size=(20, 40), tau=0.1,
        prevalence=1.0, seed=seed, simulate_field=True, smoothness_fwhm=10.0, field_zooms=4.0)
    studyset = cap_peaks(studyset, max_peaks)
    n_per = studyset.coordinates.groupby("id").size().median()
    estimator = CBES(fwhm=12.0, mask=MASK, peak_bias=peak_bias, null_method="none")
    result = estimator.fit(studyset)
    g = result.get_map("g", return_type="array").ravel()[AT]
    return float(g), float(n_per), float(estimator.peak_information_["excess_z"])


def summarise(max_peaks, peak_bias):
    rows = [one(seed, max_peaks, peak_bias) for seed in range(N_SIMS)]
    g = np.array([r[0] for r in rows])
    g = g[np.isfinite(g)]
    mse = np.mean((g - TRUE_G) ** 2)
    print(f"{str(max_peaks):>10} {str(peak_bias):>10} {np.mean([r[1] for r in rows]):6.1f} "
          f"{np.mean([r[2] for r in rows]):+7.3f} {g.mean():7.3f} "
          f"{g.mean() - TRUE_G:+7.3f}+-{g.std(ddof=1)/np.sqrt(g.size):.3f} "
          f"{mse:7.4f}+-{np.std((g - TRUE_G) ** 2, ddof=1)/np.sqrt(g.size):.4f}", flush=True)


if __name__ == "__main__":
    print(f"true g = {TRUE_G}, {N_SIMS} simulations of 24 studies, N 20-40, smooth fields\n")
    print(f"{'cap':>10} {'peak_bias':>10} {'foci/st':>6} {'excess':>7} {'mean g':>7} "
          f"{'bias':>14} {'MSE':>15}")
    for max_peaks in (5, 15, None):
        for peak_bias in (None, "per-study"):
            summarise(max_peaks, peak_bias)
