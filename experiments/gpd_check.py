"""Does the Pareto tail recover the p-values a much larger permutation set would give?

The claim is that a few hundred relocations plus an extreme-value model reproduce the tail of
thousands. So: build a large reference null, then ask whether the GPD fitted to a small subset
of it agrees with the reference's own empirical quantiles, out past where the subset alone
could resolve anything.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from nimare.meta.cbma.effectsize import _gpd_tail_p

rng = np.random.default_rng(0)

def report(name, draw):
    reference = np.sort(draw(20000))
    print(f"\n{name}: reference null of {len(reference):,} maxima")
    print(f"{'target p':>9s} {'reference stat':>15s} "
          f"{'empirical p (500)':>18s} {'GPD p (500)':>13s} {'GPD/target':>11s}")
    for target in (0.05, 0.01, 0.002, 0.0005):
        stat = float(np.quantile(reference, 1.0 - target))
        for n_small in (500,):
            small = draw(n_small)
            empirical = (1 + np.sum(small >= stat)) / (1 + n_small)
            gpd = _gpd_tail_p(np.array([stat]), small)
            gpd_p = float(gpd[0]) if gpd is not None else np.nan
            print(f"{target:9.4g} {stat:15.3f} {empirical:18.4g} {gpd_p:13.4g} "
                  f"{gpd_p / target:11.2f}")

# A maximum statistic over many correlated voxels is a Gumbel-ish extreme; test that and two
# harder shapes, since the method must not be quietly relying on one parent distribution.
report("max of 10k correlated normals", lambda n: np.array(
    [np.abs(rng.normal(0, 1, 2000)).max() for _ in range(n)]))
report("gumbel", lambda n: rng.gumbel(3.0, 0.6, n))
report("heavy tail (t3)", lambda n: np.abs(rng.standard_t(3, n)) + 2.0)
