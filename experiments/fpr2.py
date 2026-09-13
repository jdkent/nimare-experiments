"""FPR under a global null, with the Monte Carlo null powering the uncorrected p-values."""
import sys, warnings
import numpy as np, nibabel as nib
from nimare.correct import FDRCorrector, FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
warnings.simplefilter("ignore")

N_SIMS = int(sys.argv[1]); N_ITERS = int(sys.argv[2])
affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)

print(f"global null: 30 studies x 8 noise foci, {N_SIMS} sims, {N_ITERS} null iters", flush=True)
print(f"{'model':16s} {'null':12s} {'p<.05':>7s} {'p<.01':>7s} {'FWE-bonf':>9s} "
      f"{'FDR':>6s} {'FWE-mc':>7s}", flush=True)

for model in ("none", "zero-inflated"):
    for null_method in ("parametric", "montecarlo"):
        f05, f01, bonf, fdr, mc = [], [], [], [], []
        for seed in range(N_SIMS):
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
                prevalence=0.0, n_noise_foci=8, noise_extent=36.0, seed=seed)
            est = CBES(fwhm=12.0, mask=mask, selection_model=model,
                       null_method=null_method, n_iters=N_ITERS, seed=1000 * seed)
            res = est.fit(ss)
            p = res.get_map("p", return_type="array")
            f05.append(np.mean(p < 0.05)); f01.append(np.mean(p < 0.01))
            b = FWECorrector(method="bonferroni").transform(res)
            bonf.append(np.any(b.maps["p_corr-FWE_method-bonferroni"] < 0.05))
            f = FDRCorrector(method="indep", alpha=0.05).transform(res)
            fdr.append(np.any(f.maps["p_corr-FDR_method-indep"] < 0.05))
            if null_method == "montecarlo":
                m = FWECorrector(method="montecarlo", n_iters=N_ITERS).transform(res)
                mc.append(np.any(m.maps["p_corr-FWE_method-montecarlo"] < 0.05))
        print(f"{model:16s} {null_method:12s} {np.mean(f05):7.3f} {np.mean(f01):7.3f} "
              f"{np.mean(bonf):9.2f} {np.mean(fdr):6.2f} "
              f"{(np.mean(mc) if mc else float('nan')):7.2f}", flush=True)
