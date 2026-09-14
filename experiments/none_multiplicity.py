"""Do observed and relocated foci give the same number of studies per covered voxel?

Leading hypothesis for the selection_model="none" inflation. The relocation null spreads foci
uniformly over the whole mask; the observed foci occupy whatever sub-volume the generator (or a
real literature) put them in. Concentration raises the number of studies whose kernel reaches a
given covered voxel, and with "none" the standard error is a pure function of that multiplicity
and the weights -- nothing pulls it back. More studies per covered voxel in the observed map
than in the null means a systematically smaller se on the observed side, hence larger |z|, hence
an inflated rejection rate with no bug anywhere in the null machinery.

"zero-inflated" would be partly protected: its censoring term reads silence over a 2 x FWHM
sphere and pushes the estimate toward zero exactly where few studies reported, which works
against the same gradient.

Reported per condition, restricted to voxels where the statistic exists:
    n_studies   mean studies contributing
    se          mean standard error
    |z|         mean absolute statistic
for the observed fit and for 30 relocations of the same foci.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/nonemult"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)


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


def profile(est, table):
    fit, z = est._statistic(table, est._sample_sizes_, est._thresholds_, est._image_studies_)
    c = fit["covered"]
    return (float(fit["n_studies"][c].mean()), float(np.median(fit["se"][c])),
            float(np.abs(z[c]).mean()), float(c.mean()))


S = int(sys.argv[1]) if len(sys.argv) > 1 else 8
print(f"{S} datasets, 30 relocations each; statistics over covered voxels only.\n", flush=True)
print(f"{'selection':>14s} {'span':>5s} {'source':>11s} {'n_studies':>10s} {'se':>8s} "
      f"{'mean |z|':>9s} {'covered':>8s}", flush=True)

for sel in ("none", "zero-inflated"):
    for span in (20.0, 12.0):
        obs, nul = [], []
        for seed in range(S):
            est = CBES(fwhm=10.0, mask=MASK, null_method="none", selection_model=sel,
                       threshold="study-min", peak_bias="per-study", seed=seed)
            est.fit(build(seed, span=span))
            obs.append(profile(est, est._focus_table_))
            ijk = est._in_mask_ijk()
            rng = np.random.default_rng(1000 + seed)
            for _ in range(30):
                t = est._focus_table_.copy()
                t[["i", "j", "k"]] = ijk[rng.integers(0, len(ijk), size=len(t))]
                nul.append(profile(est, t))
        for name, rows in (("observed", obs), ("relocated", nul)):
            n, se, z, c = np.mean(np.array(rows), axis=0)
            print(f"{sel:>14s} {span:5.0f} {name:>11s} {n:10.2f} {se:8.4f} {z:9.4f} "
                  f"{c:8.3f}", flush=True)
        print(flush=True)
