"""Item 11: the false-positive rate, resolved properly, and under heterogeneous sample sizes.

The rate on record came from 20 simulations at 200 permutations, which cannot distinguish 0.03
from 0.05 -- the binomial standard error on 20 draws is about 0.05 itself. And every simulation
drew N from (20, 40), so the permutation null's exchangeability was never tested against the
range a real collection spans, where g_k carries very different variances between studies.

Two things measured per configuration: the uncorrected rate at nominal .05 over voxels the fit
covers, and the voxel-level FWE rate, which is the one a reader acts on. Under a global null
every rejection is a false one, so the rate is read directly.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

N_SIMS = 150
N_ITERS = 200
ALPHA = 0.05


def one(seed, sample_size, n_studies):
    """A global null: no true effect anywhere, so every rejection is a false positive."""
    studyset = create_effect_size_coordinate_studyset(
        [(0.0, 0.0, 0.0)],
        effect_sizes=0.0,
        n_studies=n_studies,
        sample_size=sample_size,
        tau=0.0,
        prevalence=0.0,
        seed=seed,
        n_noise_foci=6,
        noise_extent=60.0,
    )
    estimator = CBES(
        fwhm=10.0, peak_bias="per-study", n_iters=N_ITERS, seed=seed, cluster_threshold=None
    )
    result = estimator.fit(studyset)
    p = result.get_map("p", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    covered = k > 0
    if not covered.any():
        return np.nan, np.nan
    uncorrected = float(np.mean(p[covered] <= ALPHA))

    maps, _, _ = estimator.correct_fwe_montecarlo(result, vfwe_only=True)
    logp = maps["logp_level-voxel"].ravel()
    fwe = float(np.mean(logp[covered] >= -np.log10(ALPHA)))
    return uncorrected, fwe


print(f"{N_SIMS} simulations per configuration, {N_ITERS} permutations, nominal alpha {ALPHA}")
print("binomial standard error on the FWE rate across simulations is about "
      f"{np.sqrt(ALPHA * (1 - ALPHA) / N_SIMS):.3f}\n")
print(f"{'sample size':>16s} {'studies':>8s} {'uncorrected':>22s} {'voxel FWE':>22s}")
for label, sample_size in (("20-40 (as before)", (20, 40)),
                           ("10-500 (wide)", (10, 500)),
                           ("10-1000 (extreme)", (10, 1000))):
    for n_studies in (12, 30):
        rows = [one(seed, sample_size, n_studies) for seed in range(N_SIMS)]
        unc = np.array([r[0] for r in rows], dtype=float)
        fwe = np.array([r[1] for r in rows], dtype=float)
        # Mean rejection fraction per map, and the share of maps with any FWE rejection.
        print(f"{label:>16s} {n_studies:8d} "
              f"{f'{np.nanmean(unc):.4f} +- {np.nanstd(unc) / np.sqrt(N_SIMS):.4f}':>22s} "
              f"{f'{np.nanmean(fwe):.4f} (any: {np.nanmean(fwe > 0):.3f})':>22s}", flush=True)
