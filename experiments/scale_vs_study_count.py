"""Why does the image-calibrated scale overshoot further when there are more coordinate studies?

The weight-share model, fitted on the `peak_bias=None` rows at twelve studies, predicts the
twenty-four-study row at a new image fraction almost exactly:

    b(f) = b0 (1-f) / ((1-f) + r f),  b0 = 0.273 (excess over the all-image floor), r = 4.7
    f = 0.083 (2 of 24)  ->  predicted +0.174, measured +0.163

The same form fitted to the *calibrated* rows does not transfer. It says 2 of 24 should overshoot
by -0.027 and the measurement is -0.046, nearly twice that. Dilution is therefore not the whole
story for the calibrated configuration: something about the scale estimate itself moves with the
number of coordinate studies, even though the number of donors it is read from is fixed at two.

The candidate mechanism is the voxel set. The scale is read off by comparing, at shared voxels,
what a coordinate-only fit says against what each donor's image says. More coordinate studies
means more foci, means the coordinate fit is non-zero over more voxels, means the comparison
runs over a larger and differently-composed set -- and voxels far from any peak have a different
ratio than voxels at a peak. Task #38 already found that a conclusion about this estimator
"depends entirely on the voxel set"; this would be the same thing biting the calibration.

So: hold the donors at two, vary the coordinate studies, and read `scale_` directly. If the scale
rises with study count, the mechanism is identified and the fix is to pin the comparison to a
voxel set that does not grow. If it is flat, the overshoot is coming from somewhere else and this
note should say so rather than leave the tidy explanation standing.
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

N = int(os.environ.get("NSIMS", 24))


def one(seed, n_studies, n_image=2):
    """A fit that reports the calibration internals rather than only the estimate."""
    rng = np.random.default_rng(seed)
    studies, n_foci = [], 0
    for k in range(n_studies):
        n = int(rng.integers(20, 41))
        gmap, var_map, t_map = bed.study_fields(rng, n, 0.0)
        foci, _ = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                         bed.ZOOMS, scheme="cluster", focus="max")
        n_foci += len(foci)
        meta = {"sample_sizes": [n]}
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}
        if k < n_image:
            tag = f"sc{seed}_{n_studies}_{k}"
            gp = bed.OUT / f"{tag}_g.nii.gz"
            vp = bed.OUT / f"{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), bed.AFF), gp)
            nib.save(nib.Nifti1Image(var_map.astype(np.float32), bed.AFF), vp)
            analysis["images"] = [
                {"url": str(gp), "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": str(vp), "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    ss = Studyset({"id": "sc", "name": "sc", "studies": studies}, target=None, mask=bed.MASK)
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=True,
               peak_bias="per-study", peak_bias_scale="images")
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    g = float(res.get_map("g", return_type="array").ravel()[pos])
    scale = float(getattr(est, "_peak_bias_scale_", np.nan))
    return g, scale, n_foci, getattr(est, "scale_source_", "?")


if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; two image donors throughout; {N} replications")
    print(f"{'studies':>8s} {'mean g':>8s} {'bias':>8s} {'mean scale':>11s} {'sd scale':>9s} "
          f"{'foci/study':>11s} {'source':>10s}")
    for ns in (6, 12, 24, 36):
        rows = [r for r in Parallel(n_jobs=4)(delayed(one)(s, ns) for s in range(N))
                if r is not None]
        g = np.array([r[0] for r in rows], dtype=float)
        sc = np.array([r[1] for r in rows], dtype=float)
        foci = np.array([r[2] for r in rows], dtype=float) / ns
        okg = np.isfinite(g)
        oks = np.isfinite(sc)
        print(f"{ns:8d} {g[okg].mean():8.3f} {g[okg].mean()-bed.TRUE_G:+8.3f} "
              f"{sc[oks].mean():11.3f} {sc[oks].std(ddof=1) if oks.sum() > 1 else 0:9.3f} "
              f"{foci.mean():11.2f} {rows[0][3]:>10s}", flush=True)
    print("\nA scale that rises with the coordinate-study count identifies the voxel set as the")
    print("mechanism. A flat scale says the overshoot is elsewhere and this hypothesis is wrong.")
