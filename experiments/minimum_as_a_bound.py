"""Can the smallest reported statistic *bound* the reporting threshold, rather than estimate it?

jdkent's suggestion: where a coordinate table carries statistics, use them to inform the
threshold. Estimating the threshold from the minimum was measured and rejected -- it undoes an
order statistic it cannot identify, returned z = 4.0 against a true cluster-forming cut of 3.1,
and produced 0.201 of prevalence error against 0.008 for a plausible constant.

But there is a strictly weaker use that has not been tested. Anything a study reported *cleared*
its cut, so

    true cutoff <= min |reported statistic|

is a hard inequality, not an inference. It cannot overshoot the way the old rule did. So it can
**clamp** an assumed constant rather than replace it: use min(assumed, study minimum). That only
ever moves the cutoff down, and only for a study whose table shows the assumption was too high.

When does it bite? When the assumed constant is above what a study actually applied -- a paper
reporting at p < .01 in a collection assumed to be at p < .001, or a collection where the
assumption is simply wrong. So the bed has to include studies whose thresholds disagree with the
assumption; if every study matches, the clamp is a no-op and says nothing.

Three arms against a known truth:

  assumed         one constant for every study, what ships by default
  clamped         min(constant, that study's smallest reported |g|)
  oracle          each study's real threshold, the ceiling nothing can beat
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collections
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES, reporting_cutoff_to_g
from nimare.transforms import d_to_g, t_to_d, z_to_t

SHAPE, ZOOMS, EXTENT, BLOB = (21, 21, 21), 4.0, 40.0, 10.0
AFFINE = np.array([[ZOOMS, 0, 0, -EXTENT], [0, ZOOMS, 0, -EXTENT],
                   [0, 0, ZOOMS, -EXTENT], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
TRUE_G = 0.5
# The assumption CBES makes by default, and thresholds that straddle it so the clamp can bite.
ASSUMED_Z = 3.2905267314919255
SPREAD_Z = (2.4, 2.8, 3.29, 3.8)

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - EXTENT
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
TRUTH = MASKER.transform(nib.Nifti1Image(
    (TRUE_G * np.exp(-(grid ** 2).sum(-1) / (2 * sd ** 2))).astype(np.float32), AFFINE)).ravel()
STRATA = [("quiet", TRUTH < 0.05), ("middle", (TRUTH >= 0.05) & (TRUTH < 0.25)),
          ("effect", TRUTH >= 0.25)]


class ClampedCBES(CBES):
    """CBES whose assumed cutoff is clamped by each study's own smallest reported value.

    Only ``_study_cutoffs_z`` is touched. The clamp is a hard inequality -- anything reported
    cleared the cut -- so unlike the retired inference it can only move a cutoff downward, and
    only for a study whose table contradicts the assumption.
    """

    def _study_cutoffs_z(self, table, sample_sizes):
        cutoffs = super()._study_cutoffs_z(table, sample_sizes)
        coords = self.inputs_["coordinates"]
        column = "z_stat" if "z_stat" in coords.columns else None
        if column is None:
            return cutoffs
        smallest = coords.assign(_abs=np.abs(coords[column].astype(float))).groupby("id")[
            "_abs"
        ].min()
        bound = smallest.reindex(cutoffs.index).astype(float)
        return cutoffs.where(~np.isfinite(bound), np.minimum(cutoffs, bound))


def build(seed, tmp):
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=TRUE_G, n_studies=20, sample_size=(20, 40), tau=0.1,
        seed=seed, simulate_field=True, n_image_studies=2, image_dir=tmp,
        noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB, threshold_z=list(SPREAD_Z))


def run(cls, collection, **kwargs):
    est = cls(mask=MASK, null_method="none", **kwargs)
    res = est.fit(collection)
    return (np.abs(res.get_map("g", return_type="array").ravel()), est)


def report(title, spread, assumed):
    import tempfile
    global SPREAD_Z
    SPREAD_Z = spread
    rows = collections.defaultdict(list)
    biting = []
    for seed in range(8):
        tmp = tempfile.mkdtemp()
        ss = build(seed, tmp)
        assumed_map, est_a = run(CBES, ss, threshold=assumed)
        clamped, est_c = run(ClampedCBES, ss, threshold=assumed)
        oracle, _ = run(CBES, ss, threshold="reporting_threshold")
        moved = int(np.sum(est_c._cutoffs_z_.values < est_a._cutoffs_z_.values - 1e-9))
        biting.append((moved, len(est_a._cutoffs_z_)))
        for label, values in (("assumed", assumed_map), ("clamped", clamped), ("oracle", oracle)):
            rows[label].append(values)

    print(f"\n=== {title}: thresholds {spread} z, assumption {assumed:.2f} z")
    print(f"studies whose cutoff the clamp moved: "
          f"{np.mean([m for m, _ in biting]):.1f} of {biting[0][1]} per collection\n")
    print(f"{'threshold rule':>16} " + " ".join(f"{'rmse ' + n:>13}" for n, _ in STRATA)
          + " " + " ".join(f"{'bias ' + n:>13}" for n, _ in STRATA))
    for label in ("assumed", "clamped", "oracle"):
        maps = rows[label]
        cells = [np.mean([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in maps])
                 for _, sel in STRATA]
        bias = [np.mean([float(np.mean(m[sel] - TRUTH[sel])) for m in maps])
                for _, sel in STRATA]
        print(f"{label:>16} " + " ".join(f"{c:13.4f}" for c in cells)
              + " " + " ".join(f"{b:+13.4f}" for b in bias))

    from scipy import stats
    print("\nPaired across seeds, clamped minus assumed, rmse per stratum:")
    for name, sel in STRATA:
        a = np.array([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in rows["clamped"]])
        b = np.array([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in rows["assumed"]])
        pv = stats.ttest_rel(a, b).pvalue if np.any(a != b) else 1.0
        print(f"  {name:>8} {np.mean(a - b):+8.4f}  paired p {pv:.4f}")


if __name__ == "__main__":
    # Where the clamp should bite: the assumption sits above what half the studies applied.
    report("assumption too high for some", (2.4, 2.8, 3.29, 3.8), ASSUMED_Z)
    # Where it must be a no-op: every study really did apply the assumed cut.
    report("assumption correct for all", (ASSUMED_Z,), ASSUMED_Z)
    # Where the tables are thin: papers reporting only strong peaks, so the minimum is far
    # above the true cut and the bound is vacuous.
    report("assumption too low for all", (4.5,), ASSUMED_Z)
