"""Item 8: does treating a contrast as a study break anything?

The estimator groups on `table["id"]`, which NiMARE builds as study_id + contrast_id. A paper
contributing two contrasts therefore enters as two observations that the model treats as
independent, although they share subjects. Every test collection to date has exactly one
contrast per study, so this has never been exercised.

Two predictions worth separating:

  * `tau2` should be *under*-estimated, because two contrasts of one study agree more than two
    genuine studies would, and `se` should shrink as though there were more information than
    there is. Both inflate `z`.
  * The false-positive rate may survive anyway. The permutation null holds positions and study
    membership fixed, so a voxel reached by two contrasts of one study is reached by two in
    every permutation too -- the null inherits the same inflation and may absorb it.

Measured by splitting each simulated study into two contrasts that share its study-level effect
and differ only in sampling noise, against the same collection left whole, at matched totals.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

TRUTH = (0.0, 0.0, 0.0)
N_SIMS = 40
N_ITERS = 200
ALPHA = 0.05


def rebuild(studyset, split):
    """Re-emit the studyset, optionally as two contrasts per study sharing its foci."""
    coords = studyset.coordinates
    sizes = {
        r.study_id: float(np.mean(r.sample_sizes))
        for r in studyset.metadata.itertuples()
        if r.sample_sizes
    }
    rng = np.random.default_rng(abs(hash(tuple(coords["z_stat"].astype(float).values[:4]))) % 2**32)
    studies = []
    for sid, rows in coords.groupby("study_id"):
        n = sizes.get(sid, 30.0)
        n_contrasts = 2 if split else 1
        analyses = []
        for c in range(n_contrasts):
            points = []
            for r in rows.itertuples():
                z = float(r.z_stat)
                if c:
                    # Same study-level effect, fresh sampling noise: correlated, not identical.
                    z = z + float(rng.normal(0.0, 0.35))
                points.append(
                    {
                        "space": "MNI",
                        "coordinates": [float(r.x), float(r.y), float(r.z)],
                        "values": [{"kind": "Z", "value": z}],
                    }
                )
            analyses.append(
                {
                    "id": f"{sid}-{c + 1}",
                    "name": str(c + 1),
                    "metadata": {"sample_sizes": [n]},
                    "points": points,
                    "images": [],
                }
            )
        studies.append(
            {"id": sid, "name": sid, "metadata": {"sample_sizes": [n]}, "analyses": analyses}
        )
    return Studyset({"id": "c", "name": "c", "studies": studies})


def focus_index(masker):
    ijk = mm2vox(np.array([TRUTH]), masker.mask_img.affine)[0]
    mask = np.asarray(masker.mask_img.dataobj).astype(bool)
    flat = np.full(mask.shape, -1, dtype=np.int64)
    flat[mask] = np.arange(mask.sum())
    return int(flat[tuple(ijk)])


def measure(split, effect, seed):
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH],
        effect_sizes=effect,
        n_studies=20,
        sample_size=(20, 40),
        tau=0.15,
        prevalence=1.0 if effect else 0.0,
        seed=seed,
        n_noise_foci=4,
        noise_extent=50.0,
    )
    rebuilt = rebuild(studyset, split)
    estimator = CBES(
        fwhm=10.0, peak_bias="per-study", n_iters=N_ITERS, seed=seed, cluster_threshold=None
    )
    result = estimator.fit(rebuilt)
    p = result.get_map("p", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    tau2 = result.get_map("tau2", return_type="array").ravel()
    se = result.get_map("se", return_type="array").ravel()
    covered = k > 0
    index = focus_index(rebuilt.masker)
    return dict(
        rate=float(np.mean(p[covered] <= ALPHA)),
        k=float(np.mean(k[covered])),
        tau2=float(np.mean(tau2[covered])),
        se=float(np.mean(se[covered & np.isfinite(se)])),
        p_focus=float(p[index]) if index >= 0 else np.nan,
    )


print(f"{N_SIMS} simulations, 20 studies, {N_ITERS} permutations, nominal alpha {ALPHA}\n")
for label, effect in (("global null", 0.0), ("effect 0.8", 0.8)):
    print(f"--- {label} ---")
    print(f"{'contrasts/study':>16s} {'mean k':>8s} {'mean tau2':>10s} {'mean se':>9s} "
          f"{'p<=.05 rate':>12s} {'p at focus':>11s}")
    for split in (False, True):
        rows = [measure(split, effect, seed) for seed in range(N_SIMS)]
        print(f"{2 if split else 1:16d} "
              f"{np.mean([r['k'] for r in rows]):8.2f} "
              f"{np.mean([r['tau2'] for r in rows]):10.4f} "
              f"{np.mean([r['se'] for r in rows]):9.4f} "
              f"{np.mean([r['rate'] for r in rows]):12.4f} "
              f"{np.nanmean([r['p_focus'] for r in rows]):11.4f}", flush=True)
    print()
