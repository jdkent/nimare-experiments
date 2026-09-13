"""Re-validate the Pareto tail, averaging over repetitions rather than trusting one draw.

The failure mode that matters is anticonservatism: a corrected p reported smaller than the
truth. So the summary statistic is the *ratio* of reported p to true p, and what matters is
that it does not fall well below 1 -- being conservative is a cost, being anticonservative is
a bug.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from nimare.meta.cbma.effectsize import _gpd_tail_p

rng = np.random.default_rng(0)
N_SMALL, N_REPS = 500, 25

def report(name, draw):
    reference = np.sort(draw(40000))
    print(f"\n{name}")
    print(f"{'target p':>9s} {'stat':>8s} {'empirical':>12s} {'GPD median':>12s} "
          f"{'GPD/target':>11s} {'anticons.':>10s}")
    for target in (0.05, 0.01, 0.002, 0.0005):
        stat = float(np.quantile(reference, 1.0 - target))
        emp, gpd = [], []
        for _ in range(N_REPS):
            small = draw(N_SMALL)
            emp.append((1 + np.sum(small >= stat)) / (1 + N_SMALL))
            out = _gpd_tail_p(np.array([stat]), small)
            gpd.append(float(out[0]) if out is not None else np.nan)
        gpd = np.array(gpd)
        usable = np.isfinite(gpd)
        median = float(np.median(gpd[usable])) if usable.any() else np.nan
        # Fraction of repetitions reporting a p more than 2x smaller than the truth.
        bad = float(np.mean(gpd[usable] < target / 2.0)) if usable.any() else np.nan
        print(f"{target:9.4g} {stat:8.3f} {np.median(emp):12.4g} {median:12.4g} "
              f"{median / target:11.2f} {bad:10.2f}")

report("max of 2k correlated normals", lambda n: np.array(
    [np.abs(rng.normal(0, 1, 2000)).max() for _ in range(n)]))
report("gumbel", lambda n: rng.gumbel(3.0, 0.6, n))
report("heavy tail (t3)", lambda n: np.abs(rng.standard_t(3, n)) + 2.0)
