"""Is the relocation null itself wrong for selection_model="none", or is the observed map hot?

At span 20 the observed and relocated foci have the same spatial distribution and the same
coverage fraction (0.752 vs 0.784), so the relocation null ought to be an exact randomization
null and the rejection rate ought to be 0.050. It is 0.078. Something in the pair (observed
statistic, relocation null) is asymmetric and the coverage bookkeeping is not it.

The arbiter is a null built without any relocation at all: fit R *independent* global-null
datasets and pool their |z|. That is the true sampling distribution of the statistic. Compare
its quantiles to the relocation null's on the same statistic:

    relocation quantile < truth quantile  ->  the null is too tight, p is too small, and the
                                              relocation is the defect
    relocation quantile = truth quantile  ->  the null is right and the observed map is hot,
                                              which would mean the statistic itself is not
                                              exchangeable under relocation

zero-inflated is run alongside as the control that is known to come out near nominal.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/nonetruth"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
R = int(sys.argv[1]) if len(sys.argv) > 1 else 60


def build(seed, n_coord=15, span=20.0):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        pts = [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-span, span, 3)],
                "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
               for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


QS = [0.90, 0.95, 0.99, 0.999]
print(f"{R} independent global-null datasets give the truth; one of them also supplies a "
      f"100-relocation null.\n", flush=True)

for sel in ("none", "zero-inflated"):
    truth = []
    est = None
    for seed in range(R):
        e = CBES(fwhm=10.0, mask=MASK, null_method="montecarlo", selection_model=sel,
                 threshold="study-min", peak_bias="per-study", n_iters=100, seed=seed)
        res = e.fit(build(seed))
        truth.append(np.abs(res.get_map("z", return_type="array")))
        if seed == 0:
            est = e
    truth = np.concatenate(truth)
    hist = est.null_distributions_["histweights_corr-none_method-montecarlo"]
    edges = est.null_distributions_["histogram_bins"]
    centres = edges if len(edges) == len(hist) else 0.5 * (edges[1:] + edges[:-1])
    cdf = np.cumsum(hist) / hist.sum()
    print(f"  {sel}: |z| quantiles", flush=True)
    print(f"    {'q':>7s} {'truth':>8s} {'relocation':>11s} {'ratio':>7s}", flush=True)
    for q in QS:
        t = float(np.quantile(truth, q))
        r = float(np.interp(q, cdf, centres))
        print(f"    {q:7.3f} {t:8.3f} {r:11.3f} {r / t if t else np.nan:7.3f}", flush=True)
    # the tail probability the relocation null assigns to the truth's own 95th percentile
    t95 = float(np.quantile(truth, 0.95))
    p_reloc = 1.0 - float(np.interp(t95, centres, cdf))
    print(f"    relocation p at the truth's 0.95 point: {p_reloc:.4f} (should be 0.0500)\n",
          flush=True)
