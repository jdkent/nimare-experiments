"""Items 4 and 12: is `prevalence` recovering anything, and what radius recovers it best?

`prevalence` has been the least validated output in the estimator. Its absolute level swings
0.165 to 0.438 across coverage_radius and 0.095 to 0.635 across the threshold rule, and nothing
has ever checked it against a known value -- yet it is scale-free and so the best-identified
interesting quantity the model produces, which is exactly why it matters that it be right.

The simulator takes `prevalence` directly: below 1.0 the studies are a genuine mixture, some
with an effect of the stated size and the rest with none. So a truth exists. Recovered
prevalence is read at the ground-truth focus, where the mixture actually lives; away from it the
true prevalence is zero and the estimate is unconstrained.

Swept over the true prevalence, coverage_radius and the threshold rule, so that the two knobs
that move prevalence most can be judged against the value they are supposed to recover.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

TRUTH = (0.0, 0.0, 0.0)
N_REPLICATES = 8


def recovered(studyset, radius, rule):
    estimator = CBES(
        fwhm=10.0,
        null_method="none",
        peak_bias="per-study",
        coverage_radius=radius,
        threshold=rule,
    )
    result = estimator.fit(studyset)
    prevalence = result.get_map("prevalence", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    masker = studyset.masker
    focus = masker.transform(
        masker.inverse_transform(np.ones((1, prevalence.size), dtype=np.float32))
    )
    # Index of the ground-truth focus, via the mask's own affine.
    from nimare.utils import mm2vox

    ijk = mm2vox(np.array([TRUTH]), masker.mask_img.affine)[0]
    mask = np.asarray(masker.mask_img.dataobj).astype(bool)
    flat = np.zeros(mask.shape, dtype=np.int64) - 1
    flat[mask] = np.arange(mask.sum())
    index = int(flat[tuple(ijk)])
    if index < 0:
        return np.nan, np.nan, np.nan
    return float(prevalence[index]), float(k[index]), float(np.median(prevalence[k > 0]))


print(f"{N_REPLICATES} replicates per cell, 30 studies each, effect 0.8, tau 0.1\n")
print(f"{'true pi':>8s} {'radius':>7s} {'rule':>11s} {'pi at focus':>12s} {'k':>5s} "
      f"{'pi median':>10s}")
for true_pi in (0.25, 0.5, 0.75, 1.0):
    for radius in (8.0, 14.0, 20.0, 28.0):
        for rule in ("study-min", "pooled-min"):
            at_focus, counts, median = [], [], []
            for seed in range(N_REPLICATES):
                studyset = create_effect_size_coordinate_studyset(
                    [TRUTH],
                    effect_sizes=0.8,
                    n_studies=30,
                    sample_size=(20, 40),
                    tau=0.1,
                    prevalence=true_pi,
                    seed=seed,
                    n_noise_foci=2,
                    noise_extent=40.0,
                )
                a, kk, m = recovered(studyset, radius, rule)
                at_focus.append(a)
                counts.append(kk)
                median.append(m)
            print(f"{true_pi:8.2f} {radius:7.0f} {rule:>11s} "
                  f"{np.nanmean(at_focus):12.3f} {np.nanmean(counts):5.1f} "
                  f"{np.nanmean(median):10.3f}", flush=True)
