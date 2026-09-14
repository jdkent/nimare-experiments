import time
import numpy as np
from nimare.meta.cbma.effectsize import _censoring_terms

rng = np.random.default_rng(0)
n = 8_000_000
mu = np.abs(rng.normal(0.3, 0.2, n))
sigma = np.abs(rng.normal(0.25, 0.05, n)) + 0.05
cutoff = np.full(n, 0.55)
inv_sigma = 1.0 / sigma
args64 = (mu, cutoff * inv_sigma, 2.0 * cutoff * inv_sigma, inv_sigma, inv_sigma**2)
args32 = tuple(a.astype(np.float32) for a in args64)

def bench(args, repeats=5):
    _censoring_terms(*args)
    ts = []
    for _ in range(repeats):
        t = time.perf_counter()
        out = _censoring_terms(*args)
        ts.append(time.perf_counter() - t)
    return min(ts), out

t64, o64 = bench(args64)
t32, o32 = bench(args32)
print(f"float64 {t64*1e3:8.1f} ms   float32 {t32*1e3:8.1f} ms   speedup {t64/t32:.2f}x")
for k in o64:
    a, b = o64[k], o32[k].astype(np.float64)
    denom = np.maximum(np.abs(a), 1e-12)
    print(f"  {k:14s} max rel err {np.max(np.abs(a - b) / denom):.3e}  "
          f"median rel err {np.median(np.abs(a - b) / denom):.3e}")
print("bytes per pair:", sum(a.nbytes for a in args64) / n, "->", sum(a.nbytes for a in args32) / n)
