"""Does cluster-level FWE keep the power that the voxel maximum loses under the new null?

Voxel-level FWE power at a focal g = 0.8 across 30 studies came out at 0.030, against 0.710
uncorrected. The reason is visible in the null maxima: the observed max |z| sits only about
19% above the null median, because a within-analysis shuffle keeps every large value in the
map and merely scatters it, and a maximum over thousands of voxels picks up whichever scattered
voxel drew well. The maximum discards the thing the alternative actually produces, which is
*concentration* -- many studies large at the same place.

Cluster size and cluster mass are statistics of concentration, so they should survive the same
shuffle much better. Measured here against the same global null (for validity) and the same
focal effect (for power), with the cluster-forming threshold at .01 -- .001 is unreachable at
200 permutations, where the uncorrected p floors at 1/201.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

N_SIMS = 60
N_ITERS = 200
ALPHA = 0.05
CLUSTER_P = 0.01
TRUTH = (0.0, 0.0, 0.0)

shape, step = (25, 25, 25), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2


def one(seed, effect_size):
    mask = nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)
    at = int(np.ravel_multi_index(mm2vox(np.asarray([TRUTH]), affine)[0], shape))
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=effect_size, n_studies=30, sample_size=(20, 40), tau=0.0,
        prevalence=0.0 if effect_size == 0.0 else 1.0, seed=seed,
        n_noise_foci=6, noise_extent=40.0)
    estimator = CBES(fwhm=10.0, mask=mask, peak_bias="per-study", n_iters=N_ITERS,
                     seed=seed, cluster_threshold=CLUSTER_P)
    result = estimator.fit(studyset)
    try:
        maps, _, _ = estimator.correct_fwe_montecarlo(result, voxel_thresh=CLUSTER_P)
    except Exception:
        return None
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    out = {}
    for label, key in (("voxel", "logp_level-voxel"),
                       ("size", "logp_desc-size_level-cluster"),
                       ("mass", "logp_desc-mass_level-cluster")):
        if key not in maps:
            out[label] = np.nan
            continue
        logp = np.asarray(maps[key]).ravel()
        reject = logp >= -np.log10(ALPHA)
        out[label] = float(reject[at]) if effect_size else float(np.mean(reject[covered]) > 0)
    return out


for label, effect in (("global null (any rejection)", 0.0), ("power at g=0.8 (at the truth)", 0.8)):
    rows = [r for r in Parallel(n_jobs=4)(
        delayed(one)(seed, effect) for seed in range(N_SIMS)) if r]
    print(f"\n{label}: {len(rows)}/{N_SIMS} usable, cluster-forming p={CLUSTER_P}")
    for key in ("voxel", "size", "mass"):
        vals = np.array([r[key] for r in rows], dtype=float)
        vals = vals[np.isfinite(vals)]
        if not vals.size:
            print(f"  {key:>5}: n/a")
            continue
        print(f"  {key:>5}: {vals.mean():.3f} +- {vals.std(ddof=1)/np.sqrt(vals.size):.3f}",
              flush=True)
