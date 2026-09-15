"""Does prevalence order *voxels within one map*, which is what the docstring claims?

`prevalence_calibration` varied prevalence across collections at a single voxel and found a floor
near 0.2 at a true zero, unstable affine coefficients, and non-monotonicity at weak effects. None
of that bears directly on the claim the estimator actually makes, which is that "ordering
survives, so comparing voxels within one map is sound" -- an ordering across *voxels* inside one
fit, not across collections.

So this builds a single collection containing several sites whose prevalences differ by
construction, and asks whether the fitted map ranks those sites correctly. Two effect sizes,
because the across-collection test found the weak case close to noise, and the rank correlation is
taken over sites rather than voxels so that spatial smoothness cannot manufacture it.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

# Four sites, well separated, with prevalences the fit has to put in order.
SITES = [(-32.0, 0.0, 0.0), (-10.0, 0.0, 0.0), (12.0, 0.0, 0.0), (34.0, 0.0, 0.0)]
SITE_PREVALENCE = [0.25, 0.50, 0.75, 1.00]
N_SIMS = 12
N_STUDIES = 24


def build_mask():
    shape, step = (25, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, effect):
    """One collection whose studies each carry a random subset of the sites."""
    local = np.random.default_rng(seed)
    mask = build_mask()
    # Build the collection by generating a single-site studyset per site and merging the
    # coordinates, so each site's prevalence is controlled exactly.
    frames, sizes = [], {}
    for site, prevalence in zip(SITES, SITE_PREVALENCE):
        ss = create_effect_size_coordinate_studyset(
            [site], effect_sizes=effect, n_studies=N_STUDIES,
            sample_size=(20, 40), tau=0.0, prevalence=prevalence,
            seed=int(local.integers(1 << 30)), n_noise_foci=1, noise_extent=40.0)
        frames.append(ss.coordinates)
        for study_id, n in zip(ss.ids, ss.sample_sizes()):
            sizes[str(study_id)] = n
    import pandas as pd
    merged = pd.concat(frames, ignore_index=True)
    # Reuse one studyset's skeleton and replace its coordinates with the merged table.
    ss = create_effect_size_coordinate_studyset(
        [SITES[0]], effect_sizes=effect, n_studies=N_STUDIES, sample_size=(20, 40),
        tau=0.0, prevalence=1.0, seed=seed, n_noise_foci=1, noise_extent=40.0)
    keep = merged[merged["id"].astype(str).isin(set(ss.coordinates["id"].astype(str)))]
    if keep.empty:
        return None
    ss.coordinates = keep.reset_index(drop=True)
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none")
    result = est.fit(ss)
    pi_map = result.get_map("prevalence", return_type="array").ravel()
    out = []
    for site in SITES:
        at = mm2vox(np.asarray([site]), mask.affine)[0]
        out.append(float(pi_map[int(np.ravel_multi_index(tuple(at), mask.shape))]))
    return out


for effect in (0.8, 0.4):
    rows = Parallel(n_jobs=4)(delayed(one)(seed, effect) for seed in range(N_SIMS))
    rows = [r for r in rows if r is not None]
    if not rows:
        print(f"effect {effect}: no usable simulation\n")
        continue
    block = np.array(rows, dtype=float)
    print(f"effect {effect}, {len(rows)} simulations, {N_STUDIES} studies")
    print(f"  {'site prevalence':>16} " + " ".join(f"{p:8.2f}" for p in SITE_PREVALENCE))
    print(f"  {'fitted mean':>16} " + " ".join(f"{v:8.3f}" for v in block.mean(axis=0)))
    print(f"  {'fitted sd':>16} " + " ".join(f"{v:8.3f}" for v in block.std(axis=0)))
    # Within-map ordering: per simulation, does the fit rank the four sites correctly?
    rhos = [stats.spearmanr(row, SITE_PREVALENCE)[0] for row in block]
    exact = np.mean([list(np.argsort(row)) == list(np.argsort(SITE_PREVALENCE)) for row in block])
    print(f"  rank correlation within a map: {np.mean(rhos):+.3f} "
          f"(sd {np.std(rhos):.3f}); exact ordering in {exact:.0%} of maps\n", flush=True)
print("The docstring claims ordering within one map survives. A rank correlation near 1 supports")
print("it; one near zero would mean the ordinal reading is unsupported too.")
