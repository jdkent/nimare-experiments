"""Do a couple of images still fix the scale when only 3-20 coordinates are reported?

The reassuring numbers -- 2.04x down to 1.42x with one image, 1.19x with two -- came from the
pain collection's full peak tables, a median of 137 peaks per study. No paper reports that. In
the realistic 3-20 range the coordinate side reaches much less of the brain, and image
calibration works by taking the ratio of the image fit to the coordinate fit *over voxels both
cover*. Fewer shared voxels means a noisier ratio, and below _MIN_CALIBRATION_VOXELS the
selection widens to every shared voxel instead. So the two factors have to be varied together.

The ``images = 0`` column is the coordinates-only sweep, which is why both questions are asked
by one script: those rows say how much of the map survives at 3 peaks per study at all, and
the rest say whether an image rescues it.

Built by constructing the studyset explicitly from the peak table rather than by round-tripping
to_dict() and mutating it, which silently produced studysets with no coordinates at all.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker

# --- inverse-variance truth from the studies that never supply an image.
# Pooling all 21 and then handing some of them back as inputs measures how much one image
# resembles the average of twenty-one, not whether an image fixes the coordinate scale. The
# first N_DONOR studies are the only ones eligible to donate an image, and the truth is built
# from the rest, so no image is ever inside its own target.
N_DONOR = 5
donor_ids = {str(i) for i in list(ss.ids)[:N_DONOR]}
total = weight = None
for sid, gp, vp in zip(ss.images["id"], ss.images["g"], ss.images["g_var"]):
    if gp is None or vp is None or str(sid) in donor_ids:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    gg = np.where(ok, g, 0.0)
    total = w * gg if total is None else total + w * gg
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))

# --- peaks, sample sizes and image paths, all keyed the same way
full = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                           remove_subpeaks=True).transform(ss)
coords = full.coordinates.copy()
coords["id"] = coords["id"].astype(str)
coords["_abs"] = coords["z_stat"].abs()
sizes = {str(i): int(n) for i, n in zip(ss.ids, ss.sample_sizes())}
images = {str(r.id): (str(r.g), str(r.g_var)) for r in ss.images.itertuples()
          if r.g is not None and r.g_var is not None}

study_ids = sorted(coords["id"].unique())
matched_sizes = sum(1 for s in study_ids if s in sizes)
matched_images = sum(1 for s in study_ids if s in images)
print(f"{len(study_ids)} studies with peaks; {matched_sizes} matched to a sample size, "
      f"{matched_images} to images")
if not matched_sizes or not matched_images:
    sys.exit(f"id mismatch -- peaks {study_ids[:2]}, sizes {list(sizes)[:2]}, "
             f"images {list(images)[:2]}")
print(f"{len(coords)} peaks, median {int(coords.groupby('id').size().median())}/study")
print(f"truth from the {21 - len(donor_ids)} non-donor studies; images drawn only from the "
      f"{len(donor_ids)} donors, so no image is inside the truth\n")


def build(cap, n_images):
    """A studyset with at most ``cap`` peaks per study and images for the first few."""
    keep = (coords.sort_values("_abs", ascending=False)
            .groupby("id", sort=False).head(cap).index)
    trimmed = coords.loc[sorted(keep)]
    with_images = set(study_ids[:n_images]) & donor_ids
    studies = []
    for sid, rows in trimmed.groupby("id"):
        n = sizes[sid]
        entry = {"id": sid, "name": sid, "metadata": {"sample_sizes": [n]},
                 "analyses": [{"id": f"{sid}-1", "name": "1",
                               "metadata": {"sample_sizes": [n]},
                               "points": [
                                   {"space": "MNI",
                                    "coordinates": [float(r.x), float(r.y), float(r.z)],
                                    "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                                   for r in rows.itertuples()],
                               "images": []}]}
        if sid in with_images and sid in images:
            gp, vp = images[sid]
            entry["analyses"][0]["images"] = [
                {"url": gp, "filename": os.path.basename(gp), "space": "MNI",
                 "value_type": "g"},
                {"url": vp, "filename": os.path.basename(vp), "space": "MNI",
                 "value_type": "g_var"}]
        studies.append(entry)
    return Studyset({"id": "pain", "name": "pain", "studies": studies}), len(trimmed)


print(f"{'peaks/study':>11s} {'images':>7s} {'peaks':>6s} {'covered':>8s} {'r_top':>7s} "
      f"{'mag ratio':>10s} {'scale':>7s}")
for cap in (3, 10, 20, 10**6):
    for n_img in (0, 1, 2, 5):
        studyset, n_peaks = build(cap, n_img)
        estimator = CBES(fwhm=10.0, null_method="none", threshold="study-min",
                         peak_bias="per-study",
                         peak_bias_scale="auto" if n_img else 1.0,
                         use_images=bool(n_img))
        try:
            est = estimator.fit(studyset)
        except Exception as exc:
            print(f"{('all' if cap > 10**5 else str(cap)):>11s} {n_img:7d} {n_peaks:6d} "
                  f"{'-':>8s} {'-':>7s} {type(exc).__name__:>10s} {'-':>7s}", flush=True)
            continue
        g_est = est.get_map("g", return_type="array").ravel()
        n = est.get_map("n_studies", return_type="array").ravel()
        scored = signal & (n > 0)
        scale = float(getattr(estimator, "_peak_bias_scale_", float("nan")))
        label = "all" if cap > 10**5 else str(cap)
        if scored.sum() < 50:
            print(f"{label:>11s} {n_img:7d} {n_peaks:6d} {(n>0).mean():8.1%} {'-':>7s} "
                  f"{'too few':>10s} {scale:7.3f}", flush=True)
            continue
        rho = stats.spearmanr(np.abs(g_est[scored]), np.abs(truth[scored])).statistic
        mag = float(np.median(np.abs(g_est[scored]) / np.maximum(np.abs(truth[scored]), 1e-6)))
        print(f"{label:>11s} {n_img:7d} {n_peaks:6d} {(n>0).mean():8.1%} {rho:7.3f} "
              f"{mag:10.3f} {scale:7.3f}", flush=True)
    print(flush=True)
