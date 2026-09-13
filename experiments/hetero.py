"""What do study-specific reporting thresholds do to a coordinate effect-size map?

Each paper in the literature picks its own threshold. Here each of the 21 pain studies gets
one drawn from a realistic set, coordinates are generated at that study's own threshold, and
we ask (a) how much the reported effect size depends on the threshold rather than the brain,
and (b) which `threshold=` setting copes.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import null_peak_overshoot

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ids = list(ss.ids)
sizes = np.asarray(ss.sample_sizes(), float)
ref = np.load("/tmp/claude-0/cmp/reference_g.npy")

# Realistic spread of published thresholds.
LEVELS = {"p<.01 unc": 2.3263, "p<.001 unc": 3.0902, "p<.001 (2t)": 3.2905, "FWE-ish": 4.2649}
rng = np.random.default_rng(0)
assigned = {sid: list(LEVELS.values())[rng.integers(0, len(LEVELS))] for sid in ids}

# Threshold the whole collection once per level, then pick each study from the level it used.
# (Slicing before ImagesToCoordinates trips a metadata-sizing bug in Studyset.with_metadata.)
LEVEL_SETS = {
    z: ImagesToCoordinates(merge_strategy="demolish", z_threshold=z, two_sided=True,
                           remove_subpeaks=True).transform(ss).to_dict()
    for z in sorted(set(LEVELS.values()))
}

def coords_at(threshold_map):
    """One studyset in which each study is thresholded at its own level."""
    picked = []
    for sid in ids:
        source = LEVEL_SETS[threshold_map[sid]]
        for study in source["studies"]:
            if any(a["id"] == sid for a in study["analyses"]):
                picked.append(study)
                break
    base = dict(LEVEL_SETS[sorted(set(LEVELS.values()))[0]])
    return Studyset({**base, "studies": picked})

# (a) Does the reported effect size track the threshold rather than the brain?
print("reported |g| by the threshold the study happened to use:")
print(f"{'threshold':>12s} {'studies':>8s} {'mean reported |g|':>18s} {'null-peak prediction':>21s}")
homog = {name: Studyset(LEVEL_SETS[z]) for name, z in LEVELS.items()}
for name, z in LEVELS.items():
    df = homog[name].coordinates
    g = []
    for sid in ids:
        sub = df[df["id"] == sid]
        if not len(sub): continue
        n = sizes[ids.index(sid)]
        from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
        gg, _ = peak_stat_to_hedges_g(sub["z_stat"].astype(float).to_numpy(),
                                      np.full(len(sub), n), stat_type="z")
        g.extend(np.abs(gg))
    # What a pure-noise peak at this threshold would imply, at the median N.
    n_med = np.median(sizes)
    pred, _ = peak_stat_to_hedges_g(np.array([null_peak_overshoot(z)]),
                                    np.array([n_med]), stat_type="z")
    print(f"{name:>12s} {df['id'].nunique():8d} {np.mean(g):18.3f} {abs(pred[0]):21.3f}")

# (b) Which threshold= setting copes with a mixed literature?
mixed = coords_at(assigned)
print(f"\nmixed-threshold studyset: {len(mixed.coordinates)} foci")
print(f"{'threshold setting':>20s} {'mean g':>8s} {'r with image truth':>19s}")
for setting in ("pooled-min", "study-min"):
    g = CBES(fwhm=10.0, null_method="parametric", use_images=False,
             threshold=setting).fit(mixed).get_map("g", return_type="array").ravel()
    cov = g != 0
    print(f"{setting:>20s} {g[cov & (ref > 0.2)].mean():8.3f} "
          f"{stats.pearsonr(g[cov], ref[cov])[0]:19.3f}")
print(f"{'(image reference)':>20s} {ref[ref > 0.2].mean():8.3f}")
