"""Is selection_model="none" actually calibrated, on the simulator that was validated?

The configuration sweep put this cell at 0.078-0.081 against a 0.059 coordinate-only baseline,
which looked like a 1.4x inflation. But that sweep used a hand-rolled generator whose own
baseline sits at 0.059 where the validated simulator gives 0.041, so it inflates everything and
cannot convict anything.

Here the same question is asked with create_effect_size_coordinate_studyset -- the simulator the
false positive work was done on -- with zero-inflated as the side-by-side control. Both models
see identical data, so any difference between the rows is the selection model and nothing else.

The option's premise is that you accept a known bias in exchange for not reading silence as
evidence. That trade is only worth making if its *false positive rate* is still honest.
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.correct import FDRCorrector, FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), AFFINE)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 25


def null_set(seed):
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
        prevalence=0.0, n_noise_foci=8, noise_extent=36.0, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])


print(f"{N} global-null simulations on the validated simulator.")
print("nominal: unc 0.050 / 0.0010, corrected 0.05\n", flush=True)
print(f"{'selection model':>16s} {'null':>12s} {'unc p<.05':>10s} {'unc p<.001':>11s} "
      f"{'FWE':>6s} {'FDR':>6s}", flush=True)

for model in ("zero-inflated", "none"):
    for null in ("approximate", "montecarlo"):
        u5, u1, fwe, fdr = [], [], [], []
        for seed in range(N):
            kw = dict(n_iters=150) if null == "montecarlo" else {}
            est = CBES(fwhm=12.0, mask=MASK, null_method=null, selection_model=model,
                       threshold="study-min", peak_bias="per-study", seed=seed, **kw)
            res = est.fit(null_set(seed))
            p = res.get_map("p", return_type="array")
            u5.append(float(np.mean(p < 0.05)))
            u1.append(float(np.mean(p < 0.001)))
            fdr.append(bool(np.any(FDRCorrector(method="indep", alpha=0.05).transform(res)
                                   .get_map("p_corr-FDR_method-indep", return_type="array")
                                   < 0.05)))
            if null == "montecarlo":
                maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
                fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
        f = f"{np.mean(fwe):6.2f}" if fwe else f"{'-':>6s}"
        print(f"{model:>16s} {null:>12s} {np.mean(u5):10.4f} {np.mean(u1):11.5f} {f} "
              f"{np.mean(fdr):6.2f}", flush=True)
