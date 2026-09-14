import os, sys, time, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.transforms import ImagesToCoordinates, ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
es.CBES._load_image_studies = lambda self, dataset: {}

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
for dtype in (np.float64, np.float32, np.float64, np.float32):
    es._CENSORING_DTYPE = dtype
    stats.update(calls=0, time=0.0, pairs=0)
    t = time.perf_counter()
    CBES(fwhm=10.0, null_method="none", peak_bias="per-study").fit(ss)
    total = time.perf_counter() - t
    print(f"{np.dtype(dtype).name:8s} fit {total:6.2f}s   kernel {stats['time']:6.2f}s "
          f"({stats['time']/total*100:4.1f}%)   {stats['calls']} calls, "
          f"{stats['pairs']/max(stats['calls'],1)/1e6:.2f}M pairs/call", flush=True)
