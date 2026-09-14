"""The kernel's share, and what float32 buys, as the number of studies grows.

The pain collection is 21 studies and the kernel runs on 0.45M silent pairs per call there.
A realistic large CBMA is 50-200 studies, which is where the micro-benchmark's 1.43x was
measured (8M pairs). Studies are replicated with jittered coordinates to sweep the count
without changing anything else about the fit.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.studyset import Studyset
from nimare.transforms import ImagesToCoordinates, ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
coords = ss.coordinates
sizes = {r.study_id: float(np.mean(r.sample_sizes))
         for r in ss.metadata.itertuples() if r.sample_sizes}
mask = ss.masker.mask_img
rng = np.random.default_rng(0)

plain = es._censoring_terms
stats = {"calls": 0, "time": 0.0, "pairs": 0}


def timed(*a):
    t = time.perf_counter()
    out = plain(*a)
    stats["time"] += time.perf_counter() - t
    stats["calls"] += 1
    stats["pairs"] += a[0].size
    return out


es._censoring_terms = timed


def replicate(copies):
    studies = []
    for c in range(copies):
        jitter = rng.normal(0.0, 4.0, (len(coords), 3)) if c else np.zeros((len(coords), 3))
        offset = dict(zip(coords.index, jitter))
        for sid, rows in coords.groupby("study_id"):
            n = sizes.get(sid, 20.0)
            name = f"{sid}_r{c}"
            points = [
                {"space": "MNI",
                 "coordinates": [float(r.x + offset[i][0]), float(r.y + offset[i][1]),
                                 float(r.z + offset[i][2])],
                 "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                for i, r in zip(rows.index, rows.itertuples())
            ]
            studies.append({"id": name, "name": name, "metadata": {"sample_sizes": [n]},
                            "analyses": [{"id": f"{name}-1", "name": "1",
                                          "metadata": {"sample_sizes": [n]},
                                          "points": points, "images": []}]})
    return Studyset({"id": "big", "name": "big", "studies": studies})


for copies in (1, 3, 6):
    dataset = replicate(copies)
    n_studies = len(dataset.coordinates["study_id"].unique())
    row = []
    for dtype in (np.float64, np.float32):
        es._CENSORING_DTYPE = dtype
        best = None
        for _ in range(2):
            stats.update(calls=0, time=0.0, pairs=0)
            t = time.perf_counter()
            CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study").fit(dataset)
            total = time.perf_counter() - t
            if best is None or total < best[0]:
                best = (total, stats["time"], stats["pairs"] / max(stats["calls"], 1))
        row.append(best)
    (t64, k64, pairs), (t32, k32, _) = row
    print(f"{n_studies:4d} studies  {pairs/1e6:5.2f}M pairs/call   "
          f"kernel {k64/t64*100:4.1f}% of fit   "
          f"kernel {k64:6.2f}s -> {k32:6.2f}s ({k64/k32:.2f}x)   "
          f"fit {t64:6.2f}s -> {t32:6.2f}s ({t64/t32:.2f}x)", flush=True)
