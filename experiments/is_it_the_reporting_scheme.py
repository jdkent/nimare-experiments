"""Is the excess interval width misspecification of the reporting process?

`se/sd` runs about 2.4 in the coordinates-only arm: the reported error is more than twice the
estimator's own spread across replications. The `tau2` explanation is falsified (fitted `tau2` is
exactly 0.0000 and turning DerSimonian-Laird off changes nothing) and the censoring term is what
carries it -- `selection_model="none"` takes the `se` from 0.189 to 0.100.

That leaves two readings, and one piece of evidence already separates them. The oracle bed -- a
brute-force MLE against the same censored mixture with *known* variances -- covers 94.5% to 98.4%,
close to nominal. So the observed information is about right **where the model is true**. Which
makes the excess in this bed more likely to be misspecification of the *reporting process* than of
the likelihood.

There is an obvious candidate. The censoring term asks for the probability that a study would be
silent near a voxel given a **height** threshold. The bed extracts foci with a cluster-forming cut
and a minimum extent, so silence is a *cluster-level* event: a study can have suprathreshold
voxels and report nothing because the blob was too small. Under a pure voxelwise height threshold
-- FDR or Bonferroni -- silence really is the voxelwise event the model assumes.

So: hold everything else fixed and change only how foci are extracted. If `se/sd` falls toward 1
under a height scheme, the excess is the reporting process and the likelihood is fine. If it stays
near 2.4 under every scheme, the reporting process is exonerated and the information itself is
overstated.

The schemes are not equally strict, so the number of foci and the fraction of studies reporting
will differ between arms -- both are printed, because an arm where two studies report is not
comparable to one where twelve do, and reading `se/sd` across arms with different effective
collection sizes is the mistake this file exists to avoid making silently.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from scipy.stats import t as student_t
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import interval_coverage as bed
import reporting

N = int(os.environ.get("NSIMS", 80))


def one(seed, scheme, focus, selection_model="zero-inflated"):
    rng = np.random.default_rng(seed)
    studies, n_foci, reported = [], 0, 0
    for k in range(12):
        n = int(rng.integers(20, 41))
        _, _, t_map = bed.study_fields(rng, n, 0.0)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme=scheme, focus=focus)
        n_foci += len(foci)
        reported += bool(foci)
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}]})
    if reported < 2:
        return None
    ss = Studyset({"id": "rs", "name": "rs", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=False,
               peak_bias="per-study", selection_model=selection_model)
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    dof = (float(res.get_map("dof", return_type="array").ravel()[pos])
           if "dof" in res.maps else np.nan)
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]),
            n_foci / 12.0, reported / 12.0, dof)


if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; 12 studies, coordinates only; {N} replications\n")
    print(f"{'scheme':>18s} {'foci/study':>10s} {'report':>7s} {'mean g':>7s} {'bias':>7s} "
          f"{'mean se':>8s} {'sd':>7s} {'se/sd':>6s} {'dof':>5s} {'cov(t)':>6s}")
    for scheme, focus, selection in (
        ("cluster", "max", "zero-inflated"),
        ("fdr", "max", "zero-inflated"),
        ("fwe", "max", "zero-inflated"),
        # The control: the same cluster extraction with the censoring term switched off, so the
        # scheme rows can be read against the size of the effect that term has at all.
        ("cluster", "max", "none"),
    ):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, scheme, focus, selection) for s in range(N)) if r is not None]
        if not rows:
            print(f"{scheme + '/' + selection[:4]:>18s}   no usable fits"); continue
        g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
        dof = np.array([r[4] for r in rows], dtype=float)
        ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
        sd = g[ok].std(ddof=1)
        crit = np.where(np.isfinite(dof[ok]) & (dof[ok] > 0),
                        student_t.ppf(0.975, np.maximum(dof[ok], 1.0)), 1.96)
        cov_t = np.mean((g[ok] - crit*se[ok] <= bed.TRUE_G)
                        & (bed.TRUE_G <= g[ok] + crit*se[ok]))
        label = scheme if selection != "none" else f"{scheme} (no censor)"
        print(f"{label:>18s} {np.mean([r[2] for r in rows]):10.2f} "
              f"{np.mean([r[3] for r in rows]):7.2f} {g[ok].mean():7.3f} "
              f"{g[ok].mean()-bed.TRUE_G:+7.3f} {se[ok].mean():8.3f} {sd:7.3f} "
              f"{se[ok].mean()/max(sd,1e-9):6.2f} {np.nanmedian(dof):5.1f} {cov_t:6.2f}"
              f"  (n={ok.sum()})", flush=True)
    print("\nse/sd falling toward 1 under a height scheme puts the excess in the reporting")
    print("process. Staying near 2.4 under all three exonerates it and indicts the information.")
    print("Compare foci/study and report across rows before believing any of it.")
