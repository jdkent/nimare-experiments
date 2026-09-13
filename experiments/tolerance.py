"""Every chunk runs to max_iter. Is that mu failing to converge, or pi dragging it there?

If it is pi, loosening its half of the test buys iterations at no cost to the maps -- but that
has to be shown against the maps, not assumed. Compares runtime and map agreement across
tolerances, with the 1e-5 fit as the reference.
"""
import time, warnings; warnings.simplefilter("ignore")
import numpy as np
import nimare.meta.cbma.effectsize as fx
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import get_template

mask = get_template(space="mni152_2mm", mask="brain")
ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0), (-40, -20, 50), (40, 20, -10)],
    effect_sizes=0.6, n_studies=40, sample_size=(20, 40), seed=1,
    n_noise_foci=6, noise_extent=60.0,
    threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
)

def run(tolerance, max_iter=25):
    fx._EM_TOLERANCE = tolerance
    est = CBES(fwhm=10.0, mask=mask, null_method="parametric", max_iter=max_iter,
               threshold="study-min-corrected", peak_bias="per-study")
    start = time.perf_counter()
    res = est.fit(ss)
    elapsed = time.perf_counter() - start
    return elapsed, {n: res.get_map(n, return_type="array").ravel()
                     for n in ("g", "prevalence", "z")}

base_time, base = run(1e-5)
print(f"{'tolerance':>10s} {'fit s':>7s} {'speedup':>8s} "
      f"{'max |dg|':>9s} {'max |dpi|':>10s} {'max |dz|':>9s}")
print(f"{1e-5:10.0e} {base_time:7.1f} {1.0:8.2f}x {0.0:9.5f} {0.0:10.5f} {0.0:9.5f}")
for tolerance in (1e-4, 1e-3, 1e-2):
    elapsed, maps = run(tolerance)
    cov = base["g"] != 0
    print(f"{tolerance:10.0e} {elapsed:7.1f} {base_time/elapsed:8.2f}x "
          f"{np.abs(maps['g'][cov]-base['g'][cov]).max():9.5f} "
          f"{np.abs(maps['prevalence'][cov]-base['prevalence'][cov]).max():10.5f} "
          f"{np.abs(maps['z'][cov]-base['z'][cov]).max():9.5f}", flush=True)
fx._EM_TOLERANCE = 1e-5
