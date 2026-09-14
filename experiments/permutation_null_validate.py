"""Does the magnitude-permutation null hold up where the relocation null does not?

Three claims to check, the first two of which are asserted in the estimator's own docstrings:

  1. calibration    on a global null it should return the nominal rate, on the validated
                    simulator.

  2. mask-fill      the relocation null that used to be here was anticonservative when the foci
     robustness     occupied much less of the mask than the mask itself: 0.1337 at foci spanning
                    +/-12 mm of a +/-22 mm mask. Permutation never moves a focus, so the mask
                    fill cannot reach it -- that is the claim, and it is the one worth breaking.

  3. power          a null that is merely harder to beat is no use. Measured at a true focal
                    effect where the magnitude at the site really is elevated, which is the
                    configuration permutation is entitled to detect.

  4. few images     one image study enters at weight 1 at every voxel while the coordinate
                    studies reach only their own peaks, so it can dominate the map. Under the
                    approximate null that was removed, one image read 0.0702 against a 0.056
                    baseline while two read 0.0546.

The p floor is 1 / (1 + n_iters), so n_iters is kept high enough that 0.05 is well resolved.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/permval"); OUT.mkdir(parents=True, exist_ok=True)
SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 15
ITERS = int(sys.argv[2]) if len(sys.argv) > 2 else 200

BIG = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
BIG_MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), BIG)
SMALL = np.diag([4.0, 4.0, 4.0, 1.0]); SMALL[:3, 3] = -22.0
SMALL_MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones((12, 12, 12), np.int32), SMALL), SMALL_MASK)


def validated_null(seed):
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40), prevalence=0.0,
        n_noise_foci=8, noise_extent=36.0, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])


def validated_effect(seed):
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.8, n_studies=30, sample_size=(20, 40), prevalence=1.0,
        n_noise_foci=6, noise_extent=36.0, spatial_sd=5.0, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])


def confined(seed, span=12.0, n_coord=15):
    """The configuration that breaks relocation: foci in the middle of a larger mask."""
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        pts = [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-span, span, 3)],
                "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
               for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


print(f"{SIMS} simulations per cell, {ITERS} iterations each (p floor {1/(1+ITERS):.4f}).\n",
      flush=True)
print(f"{'condition':>26s} {'null':>20s} {'unc p<.05':>10s} {'FWE .05':>8s}", flush=True)

CASES = (
    ("global null, validated sim", validated_null, BIG_MASK, 12.0, "none", False),
    ("confined foci, wide mask", confined, SMALL_MASK, 10.0, "none", False),
    ("true effect, validated sim", validated_effect, BIG_MASK, 12.0, "zero-inflated", True),
)
for label, build, mask, fwhm, selection, is_power in CASES:
    for null in ("permute-magnitudes",):
        unc, fwe = [], []
        for seed in range(SIMS):
            est = CBES(fwhm=fwhm, mask=mask, null_method=null, selection_model=selection,
                       threshold="study-min", peak_bias="per-study", n_iters=ITERS, seed=seed)
            res = est.fit(build(seed))
            p = res.get_map("p", return_type="array")
            unc.append(float(np.mean(p < 0.05)))
            maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
            fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
        print(f"{label:>26s} {null:>20s} {np.mean(unc):10.4f} {np.mean(fwe):8.2f}", flush=True)
    print(flush=True)
