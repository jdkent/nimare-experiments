"""Is single precision worth it in the censoring kernel?

The censoring terms are about 63% of a whole-brain fit, over the estimator's only array held
per silent ``(study, voxel)`` pair. Single precision halves the traffic. The question is what
it costs: the score is ``pdf(upper) - pdf(lower)``, which cancels badly near ``mu = 0``, and
the EM's stopping rule compares a likelihood gain against 1e-6 relative -- if the kernel's own
noise were larger than that, the rule would stop reading convergence and start reading
rounding.

Fitted ``g``, the number of censoring evaluations (a proxy for total EM iterations across
chunks), peak memory and wall time, with ``_CENSORING_DTYPE`` flipped between the two. Nothing
is scored against a truth: the double-precision fit *is* the reference. The collection must be
coordinates-only -- with an image for every study there are no silent pairs and the kernel is
never called at all.
"""
import os, sys, time, tracemalloc, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.transforms import ImagesToCoordinates, ImageTransformer

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                         remove_subpeaks=True).transform(ss)
es.CBES._load_image_studies = lambda self, dataset: {}

calls = [0]
plain = es._censoring_terms


def counted(*a):
    calls[0] += 1
    return plain(*a)


es._censoring_terms = counted


def run(dtype, trace=False):
    es._CENSORING_DTYPE = dtype
    calls[0] = 0
    # tracemalloc costs about 15% and dilutes the ratio, so timing and memory are separate runs.
    if trace:
        tracemalloc.start()
    t = time.perf_counter()
    result = CBES(fwhm=10.0, null_method="none", peak_bias="per-study").fit(ss)
    elapsed = time.perf_counter() - t
    peak = tracemalloc.get_traced_memory()[1] / 2**20 if trace else float("nan")
    if trace:
        tracemalloc.stop()
    g = result.get_map("g", return_type="array").ravel()
    print(f"{np.dtype(dtype).name:8s} {elapsed:6.1f}s  {calls[0]:4d} evaluations  "
          f"{peak:7.0f} MiB peak", flush=True)
    del result
    return elapsed, g


# Interleaved repeats, best of each. A single pair of fits is not a measurement: an earlier
# version read 1.40x from one pair and 0.90x from the next, which is machine noise, not dtype.
_, g64 = run(np.float64)
_, g32 = run(np.float32)
times = {np.float64: [], np.float32: []}
for _ in range(4):
    for dtype in (np.float64, np.float32):
        times[dtype].append(run(dtype)[0])
t64, t32 = min(times[np.float64]), min(times[np.float32])
print(f"\nfloat64 {['%.1f' % x for x in times[np.float64]]}  best {t64:.1f}s")
print(f"float32 {['%.1f' % x for x in times[np.float32]]}  best {t32:.1f}s")
mem64 = run(np.float64, trace=True)
mem32 = run(np.float32, trace=True)

both = np.isfinite(g64) & np.isfinite(g32) & (g64 != 0)
d = np.abs(g32[both] - g64[both])
rel = d / np.abs(g64[both])
print(f"\nspeedup         {t64 / t32:.2f}x")
print(f"voxels          {both.sum()}")
print(f"|dg|            median {np.median(d):.3e}  max {d.max():.3e}")
print(f"relative        median {np.median(rel):.3e}  max {rel.max():.3e}")
print(f"moved over 1%   {(rel > 0.01).mean() * 100:.3f}% of voxels")
print(f"mean g          {np.mean(g64[both]):.6f} -> {np.mean(g32[both]):.6f}")
