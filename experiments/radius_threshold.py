"""Test 5: is coverage_radius = 2 x FWHM defensible, and should threshold default to study-min?

coverage_radius decides which studies count as *silent* at a voxel, and silence is the entire
input to the selection model -- so this one rule of thumb sets how much correction happens. It
has never been varied. threshold decides how surprising that silence is and how far peaks are
discounted; the default "pooled-min" applies the most liberal study's threshold to everyone,
which should under-correct the strict ones.

Both scored against the image truth (inverse-variance mean of the 21 per-study Hedges' g maps),
on pattern (rank correlation over the truth's top quartile) and magnitude (paired median ratio
there). Run at 10 and 20 peaks per study as well as the full tables, because the realistic
regime is the one that matters and coverage differs enormously between them.

For the threshold half the studies are deliberately thresholded at *different* z values, which
is the case study-min exists for; with one common threshold the two rules cannot differ.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.studyset import Studyset
from nimare.transforms import ImagesToCoordinates, ImageTransformer

U = 3.2905
ss_full = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss_full.masker

# Image truth, from every study. Held fixed across every configuration below.
total = weight = None
for gp, vp in zip(ss_full.images["g"], ss_full.images["g_var"]):
    if gp is None or vp is None:
        continue
    g = masker.transform(str(gp)).ravel()
    v = masker.transform(str(vp)).ravel()
    ok = np.isfinite(g) & np.isfinite(v) & (v > 0)
    w = np.where(ok, 1.0 / np.maximum(v, 1e-6), 0.0)
    total = w * np.where(ok, g, 0.0) if total is None else total + w * np.where(ok, g, 0.0)
    weight = w if weight is None else weight + w
truth = np.divide(total, weight, out=np.zeros_like(total), where=weight > 0)
nz = truth != 0
signal = nz & (np.abs(truth) >= np.percentile(np.abs(truth[nz]), 75))

sizes = {r.study_id: float(np.mean(r.sample_sizes))
         for r in ss_full.metadata.itertuples() if r.sample_sizes}
es.CBES._load_image_studies = lambda self, dataset: {}  # coordinates only, always


def peaks_at(thresholds):
    """Peak tables with a per-study threshold, so the two threshold rules can differ.

    Thresholded by transforming the whole studyset once per distinct threshold and then taking
    each study's rows from the pass it was assigned. Rebuilding a one-study Studyset from
    ``to_dict()`` instead loses the images' absolute paths, so ImagesToCoordinates finds
    nothing.
    """
    out = {}
    for z_thr in sorted(set(thresholds.values())):
        got = ImagesToCoordinates(merge_strategy="demolish", z_threshold=z_thr,
                                  two_sided=True, remove_subpeaks=True).transform(ss_full)
        table = got.coordinates
        for sid, assigned in thresholds.items():
            if assigned == z_thr:
                rows = table[table["study_id"] == sid]
                if len(rows):
                    out[sid] = rows.copy()
    return out


def build(tables, cap):
    studies = []
    for sid, rows in tables.items():
        rows = rows.assign(_abs=rows["z_stat"].astype(float).abs()).nlargest(cap, "_abs")
        if not len(rows):
            continue
        n = sizes.get(sid, 20.0)
        studies.append({"id": sid, "name": sid, "metadata": {"sample_sizes": [n]},
                        "analyses": [{"id": f"{sid}-1", "name": "1",
                                      "metadata": {"sample_sizes": [n]},
                                      "points": [
                                          {"space": "MNI",
                                           "coordinates": [float(r.x), float(r.y), float(r.z)],
                                           "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                                          for r in rows.itertuples()],
                                      "images": []}]})
    return Studyset({"id": "p", "name": "p", "studies": studies})


def score(studyset, **kwargs):
    result = CBES(fwhm=10.0, mask=masker.mask_img, null_method="none",
                  peak_bias="per-study", **kwargs).fit(studyset)
    g = result.get_map("g", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    prevalence = result.get_map("prevalence", return_type="array").ravel()
    covered = k > 0
    both = covered & signal & np.isfinite(g)
    if both.sum() < 50:
        return dict(covered=covered.mean(), r_top=np.nan, mag=np.nan, prev=np.nan)
    ratios = np.abs(g[both]) / np.maximum(np.abs(truth[both]), 1e-9)
    return dict(covered=covered.mean(),
                r_top=stats.spearmanr(np.abs(g[both]), np.abs(truth[both])).statistic,
                mag=float(np.median(ratios)),
                prev=float(np.median(prevalence[covered])))


common = peaks_at({sid: U for sid in sizes})

print("=" * 78)
print("(a) coverage_radius, with one common reporting threshold")
print("=" * 78)
print(f"{'peaks':>6s} {'radius':>7s} {'covered':>8s} {'r_top':>7s} {'mag':>7s} {'prev':>7s}")
for cap in (10, 20, 10_000):
    studyset = build(common, cap)
    label = "all" if cap > 1000 else str(cap)
    for radius in (8.0, 12.0, 16.0, 20.0, 26.0, 34.0):
        s = score(studyset, coverage_radius=radius)
        print(f"{label:>6s} {radius:7.0f} {100 * s['covered']:7.1f}% {s['r_top']:7.3f} "
              f"{s['mag']:7.2f} {s['prev']:7.3f}", flush=True)

print()
print("=" * 78)
print("(b) threshold rule, with studies deliberately thresholded differently")
print("=" * 78)
rng = np.random.default_rng(0)
grid = np.round(np.arange(2.8, 4.21, 0.2), 2)
mixed_thresholds = {sid: float(z) for sid, z in
                    zip(sorted(sizes), rng.choice(grid, len(sizes)))}
mixed = peaks_at(mixed_thresholds)
print("per-study thresholds drawn from a 0.2-wide grid on [2.8, 4.2]; "
      f"true mean {np.mean(list(mixed_thresholds.values())):.2f}, "
      f"pooled minimum {min(mixed_thresholds.values()):.2f}")
print(f"{'peaks':>6s} {'rule':>14s} {'covered':>8s} {'r_top':>7s} {'mag':>7s} {'prev':>7s}")
for cap in (10, 20, 10_000):
    studyset = build(mixed, cap)
    label = "all" if cap > 1000 else str(cap)
    for rule_name, rule in (("pooled-min", "pooled-min"), ("study-min", "study-min")):
        s = score(studyset, threshold=rule)
        print(f"{label:>6s} {rule_name:>14s} {100 * s['covered']:7.1f}% {s['r_top']:7.3f} "
              f"{s['mag']:7.2f} {s['prev']:7.3f}", flush=True)
    # The ceiling: hand it the thresholds it would otherwise have to infer.
    s = score(studyset, threshold=float(np.mean(list(mixed_thresholds.values()))))
    print(f"{label:>6s} {'true mean':>14s} {100 * s['covered']:7.1f}% {s['r_top']:7.3f} "
          f"{s['mag']:7.2f} {s['prev']:7.3f}", flush=True)
