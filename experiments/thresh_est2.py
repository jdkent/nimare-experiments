"""Threshold recovery when studies report only a handful of peaks -- the real regime.

Papers report tables of 5-20 peaks, not 200. The smallest of m peaks above u sits above u by
an amount that grows as m shrinks, so that is exactly when the raw minimum fails and the
order-statistic correction should earn its keep.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import infer_threshold_from_minimum
from nimare.transforms import ImageTransformer, ImagesToCoordinates

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
rng = np.random.default_rng(0)

print(f"{'imposed u':>9s} {'peaks kept':>11s} {'studies':>8s} "
      f"{'raw min':>17s} {'corrected':>17s}")
for true_u in (2.3263, 3.0902, 3.2905, 4.2649):
    cs = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=true_u, two_sided=True, remove_subpeaks=True
    ).transform(ss)
    df = cs.coordinates
    peaks = {sid: np.abs(sub["z_stat"].astype(float).to_numpy()) for sid, sub in df.groupby("id")}
    for keep in (3, 5, 10, 20, None):
        raw, fixed, counts = [], [], []
        for z in peaks.values():
            # A paper reports its *largest* peaks, so keep the top m -- but the smallest of
            # those is still the minimum of m draws, which is what the correction models.
            m = len(z) if keep is None else min(keep, len(z))
            z = np.sort(z)[::-1][:m]
            raw.append(z.min())
            fixed.append(infer_threshold_from_minimum(z.min(), m))
            counts.append(m)
        label = "all" if keep is None else str(keep)
        print(f"{true_u:9.3f} {label:>11s} {len(raw):8d} "
              f"{np.mean(raw):9.3f} ({np.mean(raw)-true_u:+.3f}) "
              f"{np.nanmean(fixed):9.3f} ({np.nanmean(fixed)-true_u:+.3f})")
    print()
