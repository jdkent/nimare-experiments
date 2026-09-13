"""Threshold recovery under the filters papers actually apply.

Not arbitrary top-m truncation -- that is not how tables are built. What real reporting does
is (a) keep one local maximum per cluster, or a few separated by >= 8 mm, and (b) drop clusters
below an extent threshold. Both are correlated with peak height, so both could in principle
push a study's smallest reported statistic above its threshold. This measures whether they do.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import infer_threshold_from_minimum
from nimare.transforms import ImageTransformer, ImagesToCoordinates

ss = ImageTransformer(target="z").transform(load_pain())

# The extent threshold is the filter not yet varied; one peak per cluster is already known
# to recover the threshold to within 0.01. "All local maxima" is dropped: it is the slow case
# and the least like a published table.
CONVENTIONS = [
    ("one peak per cluster",    dict(remove_subpeaks=True,  cluster_threshold=None)),
    ("one peak, k>=10 voxels",  dict(remove_subpeaks=True,  cluster_threshold=10)),
    ("one peak, k>=20 voxels",  dict(remove_subpeaks=True,  cluster_threshold=20)),
    ("one peak, k>=50 voxels",  dict(remove_subpeaks=True,  cluster_threshold=50)),
]

print(f"{'imposed u':>10s} {'reporting convention':>22s} {'peaks/study':>12s} "
      f"{'raw min':>17s} {'corrected':>17s}")
for true_u in (2.3263, 3.0902, 3.2905, 4.2649):
    for label, kwargs in CONVENTIONS:
        kwargs = dict(kwargs)
        kwargs.setdefault("min_distance", 8.0)
        df = ImagesToCoordinates(
            merge_strategy="demolish", z_threshold=true_u, two_sided=True, **kwargs
        ).transform(ss).coordinates
        if not len(df):
            print(f"{true_u:10.3f} {label:>22s}   (nothing reported)")
            continue
        raw, fixed, counts = [], [], []
        for _, sub in df.groupby("id"):
            z = np.abs(sub["z_stat"].astype(float).to_numpy())
            raw.append(z.min())
            fixed.append(infer_threshold_from_minimum(z.min(), len(z)))
            counts.append(len(z))
        print(f"{true_u:10.3f} {label:>22s} {np.mean(counts):12.1f} "
              f"{np.mean(raw):9.3f} ({np.mean(raw)-true_u:+.3f}) "
              f"{np.nanmean(fixed):9.3f} ({np.nanmean(fixed)-true_u:+.3f})", flush=True)
    print(flush=True)
