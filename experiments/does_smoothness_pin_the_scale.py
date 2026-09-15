"""Could a smoothness metadata field pin the scale that coordinates cannot identify?

Established, and it narrows the problem sharply. The coordinates-only bias is spatial maximum
selection: a reported focus is the maximum of a smooth field over a blob, which exceeds the value
at the blob's centre. It is **not** about how selectively peaks are reported -- going from 3.58 to
47.92 foci per study, at a much lower threshold, moves the bias from +0.249 to +0.240. And it is
not the threshold: handing the estimator each study's true cut is worth about 0.044.

`peak_bias="per-study"` is supposed to correct exactly this, from the random-field peak-height
distribution at each study's own threshold and sample size. It moves the bias by 0.004. The
reason it cannot work is that the peak-height distribution depends on the field's **smoothness**
-- how many independent resels the blob spans decides how far the maximum over it exceeds the
centre -- and the estimator is never told the smoothness. So `rho_k` is wrong by a factor common
to every study, which is precisely the shape of "an unidentified common scale".

That is a hopeful diagnosis rather than a hopeless one, because **papers report smoothness.** An
estimated FWHM is standard in an fMRI methods section, and where it is absent the applied
smoothing kernel usually is. If the residual really is a smoothness-dependent factor, a smoothness
metadata field could pin the scale with no image donor at all -- which is what four measured
remedies failed to do.

The prerequisite is that the bias actually tracks smoothness, predictably, in the direction and by
roughly the amount random-field theory says. Sweeping the bed's own smoothness tests that:

  * if the bias rises with smoothness (a smoother field means the blob spans fewer independent
    resels, so its maximum is a maximum over fewer effective samples and should be *less*
    inflated -- so the prediction is that bias *falls* as smoothness rises), the factor is
    identifiable in principle and worth modelling;
  * if the bias is flat in smoothness, the residual is something else and a smoothness field
    would buy nothing.

The direction is stated in advance, and so is the assumption it rests on, because without that
assumption it is not clean. **Holding the region fixed**, a smoother field spans fewer independent
resels, so its maximum is a maximum over fewer effective competitors and the winner's curse is
smaller -- the bias should fall as smoothness rises. But the region is not fixed here: at a fixed
cluster-forming cut a smoother field yields *larger* clusters, so the resel count inside a cluster
need not fall at all, and the two effects oppose each other. A third effect pushes the same way as
neither: with a very smooth field the observed maximum can sit further from the true peak, so the
reported focus describes a location where the truth is lower.

So a flat column would be genuinely ambiguous -- it could mean smoothness does not matter, or that
two effects of similar size cancel -- and only a monotone column in either direction is
interpretable. That is a weakness of this design and it is better to know it before reading the
output than to invent an explanation for whichever shape appears. If the column is flat the next
step is to hold the cluster region fixed by construction rather than by threshold.

`fwhm_mm` is reported alongside, from the bed's own estimate, so the sweep can be read against
the smoothness a paper would print rather than against a simulator parameter.
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

N = int(os.environ.get("NSIMS", 60))


def one(seed, smooth_vox):
    """A coordinates-only fit whose studies all share one smoothness."""
    rng = np.random.default_rng(seed)
    studies, n_foci, reported = [], 0, 0
    for k in range(12):
        n = int(rng.integers(20, 41))
        t_map = reporting.study_t_field(bed.TRUTH, n, smooth_vox, rng, shape=bed.SHAPE)
        foci, height = reporting.report_peaks(t_map[bed.MASK_BOOL], bed.MASK_BOOL, bed.SHAPE,
                                              bed.ZOOMS, scheme="cluster", focus="max")
        n_foci += len(foci)
        reported += bool(foci)
        meta = {"sample_sizes": [n]}
        if np.isfinite(height):
            meta["reporting_threshold"] = float(height)
        points = [{"space": "MNI",
                   "coordinates": [float(v) for v in nib.affines.apply_affine(bed.AFF, ijk)],
                   "values": [{"kind": "T", "value": float(zv)}]} for ijk, zv in foci]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}]})
    if reported < 2:
        return None
    ss = Studyset({"id": "sm", "name": "sm", "studies": studies}, target=None, mask=bed.MASK)
    # The threshold is handed over so the sweep isolates smoothness: `study-min` mis-infers a
    # cluster-forming cut by about 1 z, and by an amount that would itself move with smoothness.
    est = CBES(fwhm=10.0, mask=bed.MASK, null_method="none", use_images=False,
               peak_bias="per-study", threshold="reporting_threshold")
    res = est.fit(ss)
    pos = int(np.ravel_multi_index(bed.READ_AT, bed.SHAPE))
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]),
            n_foci / 12.0, reported / 12.0)


if __name__ == "__main__":
    print(f"truth {bed.TRUE_G:.3f}; 12 studies, coordinates only, true threshold supplied; "
          f"{N} replications")
    print(f"the coverage table's regime is smooth_vox {bed.SMOOTH_VOX} "
          f"(~{3.3 * bed.SMOOTH_VOX * 4.0:.1f} mm FWHM)\n")
    print(f"{'smooth_vox':>10s} {'FWHM mm':>8s} {'foci/study':>10s} {'report':>7s} "
          f"{'mean g':>7s} {'bias':>7s} {'mean se':>8s}")
    for smooth_vox in (0.5, 0.8, 1.2, 1.8):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, smooth_vox) for s in range(N)) if r is not None]
        if not rows:
            print(f"{smooth_vox:10.1f}   no usable fits"); continue
        g = np.array([r[0] for r in rows]); se = np.array([r[1] for r in rows])
        ok = np.isfinite(g) & np.isfinite(se) & (se > 0)
        # Smoothing white noise with sd sigma gives FWHM ~ 3.3 sigma, not 2.355 sigma -- the
        # error that made the first version of the coverage bed a 27 mm field at 4 mm voxels.
        print(f"{smooth_vox:10.1f} {3.3 * smooth_vox * 4.0:8.1f} "
              f"{np.mean([r[2] for r in rows]):10.2f} {np.mean([r[3] for r in rows]):7.2f} "
              f"{g[ok].mean():7.3f} {g[ok].mean()-bed.TRUE_G:+7.3f} {se[ok].mean():8.3f}"
              f"  (n={ok.sum()})", flush=True)
    print("\nA bias that falls as smoothness rises is the random-field prediction and makes the")
    print("factor identifiable from a number papers already print. A flat column means the")
    print("residual is not smoothness and a metadata field would buy nothing.")
