"""Is the large-effect inflation caused by selection, or by the pooling that follows it?

The conversion is exact (g -> z -> g round-trips to 1.0000 up to g = 2.4), so the inflation
at large true effects is either in what gets reported or in what the fit does with it. Two
discriminators, each removing one suspect:

  * ``selection_model="none"`` drops the censored mixture entirely. If the inflation survives,
    the zero-inflated fit is not causing it.
  * **Truth-substituted values.** Each reported focus keeps its location but its statistic is
    replaced by the one the *true* g at that location would have produced. Selection still
    decided which locations appear, but the values carry no winner's curse at all. If the
    estimate still runs high, the inflation is in the pooling, not in the reported magnitudes.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z
from nimare.utils import mm2vox
from peak_bias_realistic_peaks import cap_peaks

TRUTH = (0.0, 0.0, 0.0)
N_SIMS = 40
shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2
MASK = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
AT = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))


def substitute_truth(studyset):
    """Rebuild with each focus's statistic set to what its own true g would produce."""
    coords = studyset.coordinates
    sizes = {}
    for study in studyset.studies:
        for analysis in study.analyses:
            sizes[(str(study.id), str(analysis.id))] = dict(analysis.metadata)
    studies = []
    for (study_id, contrast_id), sub in coords.groupby(["study_id", "contrast_id"], sort=False):
        info = sizes[(str(study_id), str(contrast_id))]
        n = float(info["sample_sizes"][0])
        j = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
        points = []
        for r in sub.itertuples():
            g_here = float(r.value_trueg)
            z = float(t_to_z(np.array([g_here / j * np.sqrt(n)]), dof=n - 1)[0])
            points.append({"space": "MNI",
                           "coordinates": [float(r.x), float(r.y), float(r.z)],
                           "values": [{"kind": "Z", "value": z}]})
        studies.append({"id": str(study_id), "name": str(study_id), "metadata": info,
                        "analyses": [{"id": str(contrast_id), "name": "1",
                                      "metadata": info, "points": points}]})
    return Studyset({"id": "truthsub", "name": "truthsub", "studies": studies})


def run(true_g, variant):
    values = []
    for seed in range(N_SIMS):
        studyset = create_effect_size_coordinate_studyset(
            [TRUTH], effect_sizes=true_g, n_studies=24, sample_size=(20, 40), tau=0.1,
            prevalence=1.0, seed=seed, simulate_field=True, smoothness_fwhm=10.0,
            field_zooms=4.0)
        studyset = cap_peaks(studyset, 5)
        selection = "zero-inflated"
        if variant == "truth-substituted":
            studyset = substitute_truth(studyset)
        elif variant == "selection=none":
            selection = "none"
        estimator = CBES(fwhm=12.0, mask=MASK, peak_bias=None, null_method="none",
                         selection_model=selection)
        try:
            result = estimator.fit(studyset)
        except ValueError:
            continue
        values.append(float(result.get_map("g", return_type="array").ravel()[AT]))
    g = np.array([v for v in values if np.isfinite(v)])
    if g.size < 5:
        print(f"{true_g:7.1f} {variant:>20} too few usable", flush=True)
        return
    print(f"{true_g:7.1f} {variant:>20} {g.size:6d} {g.mean():8.3f} "
          f"{g.mean() / true_g:8.3f}+-{g.std(ddof=1)/np.sqrt(g.size)/true_g:.3f}", flush=True)


print(f"{N_SIMS} simulations, 24 studies, N 20-40, 5 peaks/study, peak_bias=None\n")
print(f"{'true g':>7} {'variant':>20} {'n':>6} {'mean g':>8} {'ratio to truth':>15}")
for true_g in (0.8, 1.2, 1.6):
    for variant in ("as reported", "selection=none", "truth-substituted"):
        run(true_g, variant)
