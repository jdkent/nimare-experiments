"""Is the residual excess width caused by kernel weights that depend on the data?

`why_se_exceeds_sd.py` located most of the over-wide interval in the censoring term: switching it
off takes `se/sd` from 2.28 to 1.35 coordinates-only, from 1.67 to 0.93 under heterogeneity and
from 1.33 to 1.12 with six donors. The 0.93 and the 1.12 are about where a calibrated *conditional*
`se` belongs. The 1.35 is not, and it is specific to the coordinates-only homogeneous case.

One candidate was cheap to rule out by reading the code: `_apply_peak_bias` scales `g` by `rho_k`
and `var_g` by `rho_k**2`, consistently, so the correction is not putting the estimate and its
variance on different scales.

What is left is that **the kernel weights are random and correlated with the values.** A
coordinate study contributes at the kernel weight of whichever of its foci is nearest the voxel,
and both which focus that is and how far away it lands are functions of that study's noise. The
pooled variance treats those weights as known constants. A study that drew a high peak close to
the voxel contributes a large value *at a large weight*; the formula does not know the weight was
chosen partly by the same noise that set the value.

The dilution pattern fits. With donors, weight-1 image contributions dominate the weighted sum and
the effect is diluted: 1.35 -> 1.12. With real heterogeneity the between-study variance dominates
the weighting noise: 1.35 -> 0.93.

**This is a mechanism probe, not a reporting scenario, and it deliberately breaks realism in one
specific way.** The foci are placed at *fixed* positions -- the five site centres -- with only the
magnitudes drawn, so the kernel weights are identical in every replication and the only randomness
left is in the values. Nothing is capped: every study reports at every site. If `se/sd` falls to 1
when the weights stop varying, the variance formula's treatment of them is the residual. If it
stays near 1.35, the weights are exonerated and the residual is elsewhere.

The comparison arm keeps the same fixed positions but draws the *number* of reporting sites per
study from the same reporting process, so position-randomness is removed while
which-studies-report randomness is kept. That separates "the weights vary" from "the roster
varies", which the pooled variance also treats as fixed.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.studyset import Studyset
import interval_coverage as bed

N = int(os.environ.get("NSIMS", 120))


def one(seed, mode):
    """``mode`` is 'fixed' (every study reports every site) or 'roster' (sites drop out)."""
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(12):
        n = int(rng.integers(20, 41))
        meta = {"sample_sizes": [n]}
        points = []
        for (i, j, l), weight in zip(bed.SITE_IJK, bed.SITE_WEIGHT):
            if mode == "roster" and rng.random() < 0.4:
                continue
            truth = bed.PEAK_G * weight
            # One draw per site, on the t scale the estimator expects, with no spatial field --
            # so the position is fixed by construction and only the magnitude is random.
            t_value = truth * np.sqrt(n) + rng.standard_normal()
            if t_value <= 0:
                continue
            points.append({"space": "MNI",
                           "coordinates": [float(v) for v in nib.affines.apply_affine(
                               bed.AFF, (i, j, l))],
                           "values": [{"kind": "T", "value": float(t_value)}]})
        if not points:
            continue
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}]})
    if len(studies) < 2:
        return None
    ss = Studyset({"id": "w", "name": "w", "studies": studies}, target=None, mask=bed.MASK)
    # The censoring term is off throughout: it is already known to carry most of the excess, and
    # leaving it on would hide whatever the weights do.
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=False,
               peak_bias=None, selection_model="none")
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]),
            float(res.get_map("n_eff", return_type="array").ravel()[pos]))


if __name__ == "__main__":
    print("selection_model='none' throughout, peak_bias=None, foci at fixed site centres")
    print(f"{N} replications; the reference is se/sd, which should be at or just below 1\n")
    print(f"{'mode':>8s} {'mean g':>8s} {'mean se':>8s} {'sd':>7s} {'se/sd':>6s} "
          f"{'n_eff':>6s} {'sd(n_eff)':>10s}")
    for mode in ("fixed", "roster"):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, mode) for s in range(N)) if r is not None]
        if not rows:
            print(f"{mode:>8s}   no usable fits"); continue
        g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
        ne = np.array([r[2] for r in rows])
        ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
        sd = g[ok].std(ddof=1)
        print(f"{mode:>8s} {g[ok].mean():8.3f} {se[ok].mean():8.3f} {sd:7.3f} "
              f"{se[ok].mean()/max(sd,1e-9):6.2f} {ne[ok].mean():6.2f} "
              f"{ne[ok].std(ddof=1):10.3f}  (n={ok.sum()})", flush=True)
    print("\n'fixed' has constant weights and a constant roster: se/sd near 1 there and near 1.35")
    print("in the spatial bed would put the residual in the data-dependent weights. 'roster' adds")
    print("back only which-studies-report randomness, so the gap between the two rows is the part")
    print("the variance formula loses by treating the roster as fixed.")
