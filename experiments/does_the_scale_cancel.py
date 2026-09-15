"""Does the unidentified scale cancel in a ratio? Asked without a normalizer this time.

`relative_interval_coverage.py` tried to test this by comparing `g_relative` against
`mu / P95(|mu|)` and got coverage 0.00 in every arm. That was my test design, not the estimator:
the fitted P95 of |g| came out 1.212 while the truth's P95 of |mu| is 0.175, because the truth
field is mostly zero over a 30^3 volume while the fitted map is kernel-smoothed, censored and
inflated over its covered support. Two quantiles of differently shaped distributions do not
cancel a common factor. That is precisely the "it depends entirely on the voxel set" mistake
recorded as task #38, made in a script whose own comment warns against it.

The question does not need a normalizer. If the peak-height inflation is one common
multiplicative factor `c`, then g = c*mu everywhere and the *ratio between two voxels* is free of
it exactly:

    g(v1) / g(v2)  ->  mu(v1) / mu(v2)

That is the claim "coordinates identify the pattern but not the scale", stated so it can be
falsified. Two sites of the bed differ only in magnitude -- 0.800 and 0.600 -- so the truth is
0.800/0.600 = 1.333 and nothing about the voxel set enters.

The interval comes from the delta method on a log ratio of two independent estimates,

    se(log(g1/g2))^2 = (se1/g1)^2 + (se2/g2)^2

which treats the two voxels as independent. They are not -- both are 10 mm-kernel fits and the
sites are far apart but share studies -- so if that assumption is what fails, `se/sd` comes out
below 1 and this test says so rather than hiding it.

Three arms: coordinates only at twelve and twenty-four studies, where the interval on `g` covers
0.75 and 0.35, and two donors at twelve, where it covers 0.99. If the ratio covers coordinates-only
then the scale really is the whole problem and a ratio-valued output is worth having. If it does
not, the inflation is not a common factor and "coordinates identify the pattern" is too generous.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import interval_coverage as bed
import reporting

N = int(os.environ.get("NSIMS", 100))
#: Two sites of the bed that differ only in magnitude: SITE_WEIGHT 1.0 and 0.75 on PEAK_G.
SITE_A, SITE_B = bed.SITE_IJK[0], bed.SITE_IJK[2]
TRUE_A = bed.PEAK_G * bed.SITE_WEIGHT[0]
TRUE_B = bed.PEAK_G * bed.SITE_WEIGHT[2]
TRUE_RATIO = TRUE_A / TRUE_B


def one(seed, n_studies, n_image):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = bed.study_fields(rng, n, 0.0)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            tag = f"rt{seed}_{n_studies}_{n_image}_{k}"
            gp, vp = bed.OUT / f"{tag}_g.nii.gz", bed.OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), bed.AFF), gp)
            nib.save(nib.Nifti1Image(var_map.astype(np.float32), bed.AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    ss = Studyset({"id": "rt", "name": "rt", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=n_image > 0,
               peak_bias="per-study", peak_bias_scale="images" if n_image else 1.0)
    res = est.fit(ss)
    g = res.get_map("g", return_type="array").ravel()
    se = res.get_map("se", return_type="array").ravel()
    pa = int(np.ravel_multi_index(SITE_A, bed.SHAPE))
    pb = int(np.ravel_multi_index(SITE_B, bed.SHAPE))
    ga, gb, sa, sb = float(g[pa]), float(g[pb]), float(se[pa]), float(se[pb])
    if not (np.isfinite(ga) and np.isfinite(gb) and gb > 0 and ga > 0):
        return None
    log_ratio = np.log(ga / gb)
    se_log = np.sqrt((sa / ga) ** 2 + (sb / gb) ** 2)
    return ga, gb, log_ratio, se_log


if __name__ == "__main__":
    print(f"site A truth {TRUE_A:.3f}, site B truth {TRUE_B:.3f}, ratio {TRUE_RATIO:.3f}; "
          f"{N} replications\n")
    print(f"{'arm':22s} {'mean gA':>8s} {'mean gB':>8s} {'ratio':>7s} {'bias':>7s} "
          f"{'se/sd':>6s} {'cover':>6s}")
    target = np.log(TRUE_RATIO)
    for n_studies, n_image in ((12, 0), (24, 0), (12, 2)):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, n_studies, n_image) for s in range(N)) if r is not None]
        if not rows:
            print(f"{n_studies} studies, {n_image} img: no usable fits"); continue
        lr = np.array([r[2] for r in rows]); sl = np.array([r[3] for r in rows])
        ok = np.isfinite(lr) & np.isfinite(sl) & (sl > 0)
        sd = lr[ok].std(ddof=1)
        cover = np.mean((lr[ok] - 1.96*sl[ok] <= target) & (target <= lr[ok] + 1.96*sl[ok]))
        label = f"{n_studies} studies, {n_image} img"
        print(f"{label:22s} {np.mean([r[0] for r in rows]):8.3f} "
              f"{np.mean([r[1] for r in rows]):8.3f} {np.exp(lr[ok].mean()):7.3f} "
              f"{np.exp(lr[ok].mean())-TRUE_RATIO:+7.3f} {sl[ok].mean()/sd:6.2f} {cover:6.2f}"
              f"  (n={ok.sum()})", flush=True)
    print("\nThe ratio needs no normalizer and no voxel set, so a failure here is the estimator")
    print("rather than the comparison. se/sd below 1 would indict the independence assumption")
    print("in the delta method instead.")
