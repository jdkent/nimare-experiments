"""Why does a large true effect come back inflated, when selection predicts the opposite?

Scoring CBES against simulator truth at 5 peaks per study gave a ratio of estimate to truth
that *grows* with the true effect: 1.26 at g = 0.8, 1.50 at 1.2, 2.23 at 1.6. The winner's
curse predicts the reverse -- a large effect clears the reporting threshold easily, so it is
barely selected and needs little correction.

So something other than selection is acting at the high end. Three candidates, separated here
before any of it is allowed to drive a recommendation:

  1. the simulator: are the reported peaks themselves inflated relative to the true g at their
     own location? The coordinate table carries ``value_trueg``, so this is a direct read.
  2. the conversion: does ``peak_stat_to_hedges_g`` recover the true g from the reported z at
     the sample size the study had, or does the z->t->d->g chain drift at large z?
  3. the estimator: given reported values whose conversion is sound, does the fit inflate?

Reported per true effect size, with no CBES fit involved in 1 and 2.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from peak_bias_realistic_peaks import cap_peaks

print("24 studies, N 20-40, smooth fields, 5 peaks/study, 20 seeds\n")
print(f"{'true g':>7} {'reported z':>11} {'true g at':>10} {'converted':>10} "
      f"{'conv/true':>10} {'peak/true':>10}")
print(f"{'':>7} {'(mean)':>11} {'the peak':>10} {'g (mean)':>10} {'at peak':>10} {'at peak':>10}")
for true_g in (0.4, 0.8, 1.2, 1.6):
    z_all, trueg_all, conv_all = [], [], []
    for seed in range(20):
        studyset = create_effect_size_coordinate_studyset(
            [(0.0, 0.0, 0.0)], effect_sizes=true_g, n_studies=24, sample_size=(20, 40),
            tau=0.1, prevalence=1.0, seed=seed, simulate_field=True, smoothness_fwhm=10.0,
            field_zooms=4.0)
        studyset = cap_peaks(studyset, 5)
        coords = studyset.coordinates
        sizes = {}
        for study in studyset.studies:
            for analysis in study.analyses:
                sizes[(str(study.id), str(analysis.id))] = float(
                    dict(analysis.metadata)["sample_sizes"][0])
        n = np.array([sizes[(str(r.study_id), str(r.contrast_id))] for r in coords.itertuples()])
        z = np.abs(coords["z_stat"].astype(float).values)
        truth_here = np.abs(coords["value_trueg"].astype(float).values)
        converted, _ = peak_stat_to_hedges_g(z, n, stat_type="z", design="one-sample")
        # Only peaks that actually sit on signal: a pure-noise peak has a true g of zero and
        # a ratio against it is meaningless.
        on_signal = truth_here > 0.05 * true_g
        z_all.append(z[on_signal])
        trueg_all.append(truth_here[on_signal])
        conv_all.append(np.abs(converted)[on_signal])
    z = np.concatenate(z_all)
    truth_here = np.concatenate(trueg_all)
    converted = np.concatenate(conv_all)
    print(f"{true_g:7.1f} {z.mean():11.3f} {truth_here.mean():10.3f} {converted.mean():10.3f} "
          f"{np.mean(converted / truth_here):10.3f} {truth_here.mean() / true_g:10.3f}",
          flush=True)
