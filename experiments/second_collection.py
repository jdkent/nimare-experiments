"""Does the magnitude overestimate replicate on a second, independent collection?

Every magnitude figure quoted for CBES comes from the 21 NIDM pain studies. One collection
cannot tell a property of the method from a property of that literature -- pain is a strong,
spatially consistent effect with high prevalence, which is the friendliest case there is.

Here the same procedure runs on the NeuroVault "animal" set: 11 studies, movie-watching,
as-Animal contrast, unthresholded effect and variance maps. Scoring is identical to the pain
run and to the SDM comparison before it -- inverse-variance truth from all the images, CBES on
peaks thresholded out of those same images, magnitude on the truth's top quartile.

What this does *not* fix is the circularity: the coordinates are still manufactured by
thresholding the images with the reporting model CBES assumes. It tests whether the bias is
general or collection-specific, not whether the reporting model is right.
"""
import glob, os, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
from scipy import stats
from nilearn.image import resample_to_img
from nimare.meta.cbma import CBES
from nimare.utils import get_masker, get_template
from nimare.transforms import ImagesToCoordinates
from nimare.studyset import Studyset

U = 3.2905
HOME = os.path.expanduser("~/.nimare")
studies = sorted({os.path.basename(p).split("-")[1] for p in
                  glob.glob(f"{HOME}/study-*-animal_2.0x2.0x2.0_g.nii.gz")})
print(f"{len(studies)} studies: {' '.join(studies)}")

masker = get_masker(get_template(space="mni152_2mm", mask="brain"))
records, total, weight = [], None, None
for sid in studies:
    gp = f"{HOME}/study-{sid}-animal_2.0x2.0x2.0_g.nii.gz"
    vp = f"{HOME}/study-{sid}-animal_2.0x2.0x2.0_g_var.nii.gz"
    zp = f"{HOME}/study-{sid}-animal_2.0x2.0x2.0_z.nii.gz"
    if not all(os.path.exists(q) for q in (gp, vp, zp)):
        continue
    g = masker.transform(gp).ravel(); v = masker.transform(vp).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    gg = np.where(ok, g, 0.0)
    total = w * gg if total is None else total + w * gg
    weight = w if weight is None else weight + w
    records.append((sid, gp, vp, zp))
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))
print(f"truth from {len(records)} images; {int(signal.sum())} voxels in the top quartile")

# a studyset carrying the images, so ImagesToCoordinates can make the peaks
entries = []
for sid, gp, vp, zp in records:
    entries.append({"id": sid, "name": sid, "metadata": {"sample_sizes": [30]},
        "analyses": [{"id": f"{sid}-1", "name": "1", "metadata": {"sample_sizes": [30]},
            "points": [], "images": [
                {"url": gp, "filename": os.path.basename(gp), "space": "MNI",
                 "value_type": "g"},
                {"url": vp, "filename": os.path.basename(vp), "space": "MNI",
                 "value_type": "g_var"},
                {"url": zp, "filename": os.path.basename(zp), "space": "MNI",
                 "value_type": "z"}]}]})
ss = Studyset({"id": "animal", "name": "animal", "studies": entries})
peaks = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                            remove_subpeaks=True).transform(ss)
n_foci = 0 if peaks.coordinates is None else len(peaks.coordinates)
print(f"{n_foci} peaks extracted at z > {U}\n")
if not n_foci:
    sys.exit("no peaks: nothing to compare")

print("scored only where CBES makes an estimate: 35 peaks over a 2 mm whole brain reach")
print("very little of the truth's top quartile, and g is exactly 0 elsewhere.\n")
print(f"{'configuration':>26s} {'covered':>8s} {'rho':>7s} {'mag ratio':>10s} "
      f"{'n vox scored':>13s}")
for label, kw in (("default", {}),
                  ("peak_bias=per-study", dict(peak_bias="per-study")),
                  ("scale=reference", dict(peak_bias="per-study",
                                           peak_bias_scale="reference")),
                  ("scale=auto", dict(peak_bias="per-study", peak_bias_scale="auto"))):
    est = CBES(fwhm=10.0, null_method="none", threshold="study-min", use_images=False,
               **kw).fit(peaks)
    g = est.get_map("g", return_type="array").ravel()
    n_studies = est.get_map("n_studies", return_type="array").ravel()
    covered = n_studies > 0
    scored = signal & covered
    if scored.sum() < 20:
        print(f"{label:>26s} {covered.mean():8.3%} {'-':>7s} {'-':>10s} {scored.sum():13d}")
        continue
    rho = stats.spearmanr(np.abs(g[scored]), np.abs(truth[scored])).statistic
    mag = float(np.median(np.abs(g[scored]) / np.maximum(np.abs(truth[scored]), 1e-6)))
    print(f"{label:>26s} {covered.mean():8.3%} {rho:7.3f} {mag:10.3f} {scored.sum():13d}")
