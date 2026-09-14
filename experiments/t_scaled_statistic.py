"""Would a t-scaled statistic buy anything in the permutation null?

`se` is miscalibrated against a normal reference (73-94% coverage of nominal 95%, fixed to
91-97% by a t on n_eff - 1), and |z| = |g/se| is also the statistic the null ranks. So the
question is whether ranking by a dof-aware statistic detects effects that |z| misses.

One prediction constrains where it could possibly help. The uncorrected p-value is a *per-voxel*
permutation rank: positions never move, so a voxel is reached by the same studies in every
iteration and its dof is constant across the null. Any transform that depends only on |z| and
dof is therefore monotone *within* a voxel, and cannot change that voxel's rank -- so the
uncorrected map must come out bit-identical. It is asserted below rather than assumed.

That leaves familywise correction, where the maximum is taken *across* voxels with different
dof. There a low-dof voxel with a large |z| can dominate the max-statistic null and raise the
bar for everywhere else; putting voxels on a common scale should stop that. Measured as FWE
power at a known focal effect, and FWE error under a global null.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
N_SIMS = 60
N_ITERS = 200
ALPHA = 0.05

plain_statistic = es.CBES._statistic


def t_scaled(self, table, sample_sizes, thresholds, image_studies=None):
    """|g/se| referred to a t on n_eff - 1, re-expressed as an equivalent normal deviate.

    The re-expression matters only in that it keeps the statistic on a familiar scale; what
    changes is that a given |g/se| counts for less where fewer studies contribute.
    """
    fit, z_values = plain_statistic(self, table, sample_sizes, thresholds, image_studies)
    dof = np.clip(fit["n_eff"] - 1.0, 0.0, None)
    usable = (dof > 0) & np.isfinite(z_values) & (z_values != 0)
    scaled = np.zeros_like(z_values)
    if np.any(usable):
        tail = stats.t.sf(np.abs(z_values[usable]), dof[usable])
        tail = np.clip(tail, 1e-300, 0.5)
        scaled[usable] = np.sign(z_values[usable]) * stats.norm.isf(tail)
    return fit, scaled


def focus_index(masker):
    ijk = mm2vox(np.array([TRUTH]), masker.mask_img.affine)[0]
    mask = np.asarray(masker.mask_img.dataobj).astype(bool)
    lookup = np.full(mask.shape, -1, dtype=np.int64)
    lookup[mask] = np.arange(mask.sum())
    return int(lookup[tuple(ijk)])


def run(seed, effect, statistic):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH],
        effect_sizes=effect,
        n_studies=20,
        sample_size=(20, 40),
        tau=0.1,
        prevalence=1.0 if effect else 0.0,
        seed=seed,
        n_noise_foci=4,
        noise_extent=50.0,
    )
    es.CBES._statistic = statistic
    try:
        estimator = CBES(
            fwhm=10.0, peak_bias="per-study", n_iters=N_ITERS, seed=seed, cluster_threshold=None
        )
        result = estimator.fit(studyset)
        uncorrected = result.get_map("p", return_type="array").ravel().copy()
        maps, _, _ = estimator.correct_fwe_montecarlo(result, vfwe_only=True)
        logp = maps["logp_level-voxel"].ravel()
    finally:
        es.CBES._statistic = plain_statistic
    k = result.get_map("n_studies", return_type="array").ravel()
    covered = k > 0
    index = focus_index(studyset.masker)
    detected = bool(index >= 0 and logp[index] >= -np.log10(ALPHA))
    return uncorrected, covered, float(np.mean(logp[covered] >= -np.log10(ALPHA))), detected


print(f"{N_SIMS} simulations, 20 studies, {N_ITERS} permutations, alpha {ALPHA}\n")

# The prediction: the uncorrected map is a within-voxel rank, so it cannot move.
same = 0
for seed in range(6):
    p_plain, covered, _, _ = run(seed, 0.8, plain_statistic)
    p_scaled, _, _, _ = run(seed, 0.8, t_scaled)
    same += int(np.array_equal(p_plain[covered], p_scaled[covered]))
print(f"uncorrected p identical under the transform: {same} of 6 seeds")
print("(predicted: all of them -- a within-voxel monotone transform preserves ranks)\n")

print(f"{'statistic':>12s} {'condition':>13s} {'FWE rejection rate':>20s} {'focus detected':>16s}")
for label, statistic in (("|z|", plain_statistic), ("t-scaled", t_scaled)):
    for condition, effect in (("global null", 0.0), ("effect 0.8", 0.8)):
        rows = [run(seed, effect, statistic) for seed in range(N_SIMS)]
        rate = float(np.mean([r[2] for r in rows]))
        power = float(np.mean([r[3] for r in rows]))
        print(f"{label:>12s} {condition:>13s} {rate:20.4f} {power:16.3f}", flush=True)
