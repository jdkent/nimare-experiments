"""The same bias in three lines of scalar arithmetic, with no brain and no estimator.

Hedges' variance is `1/n + g^2 / (2(n-1))`. Inverse-variance weighting with a variance computed
from the *observed* effect gives a study that drew high a larger variance and less weight, so the
pooled estimate is pulled toward zero. If the effect in the coverage bed is this, it must appear
here too, at the same magnitude, in a model with nothing in it but the weights.
"""
import numpy as np

rng = np.random.default_rng(0)
REPS = 200_000
for n_studies in (12, 24):
    for mu in (0.4, 0.8, 1.2):
        n = rng.integers(20, 41, size=(REPS, n_studies))
        g = mu + rng.standard_normal((REPS, n_studies)) / np.sqrt(n)
        v_hedges = 1.0 / n + g**2 / (2.0 * (n - 1))
        v_flat = 1.0 / n.astype(float)
        out = []
        for v in (v_hedges, v_flat):
            w = 1.0 / v
            out.append((w * g).sum(1) / w.sum(1))
        bias_h = out[0].mean() - mu
        bias_f = out[1].mean() - mu
        print(f"{n_studies:2d} studies, mu {mu:.1f}:  Hedges weights {bias_h:+.4f}  "
              f"({100*bias_h/mu:+.1f}%)   draw-independent {bias_f:+.4f}")
print("\nThe coverage bed measured -0.026 on mu = 0.800 with twelve studies, and -0.002 with a")
print("draw-independent variance. If those numbers appear below, the mechanism is the weights.")
