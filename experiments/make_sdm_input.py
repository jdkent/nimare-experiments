"""Prepare a small paired dataset: derived coordinates plus the images they came from.

Both methods get exactly the same coordinates, derived from the same images at a known
threshold, so the comparison is not confounded by curation. The images are the truth both are
trying to recover, which is what makes magnitude comparable rather than just pattern.

SDM reads a tab-separated ``sdm_table.txt`` naming the studies, their sample sizes and the
*t* threshold each was reported at, plus one ``<study>.spm_mni.txt`` per study listing
``x,y,z,t``. It infers a one-sample design from the absence of an ``n2`` column, so the
threshold and the peak statistics are both t on ``N - 1`` degrees of freedom -- not z, which
is what NiMARE carries. Studies whose images yield no suprathreshold peak get the
``.no_peaks.txt`` marker file rather than being dropped, so both methods see the same 21.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates, z_to_t

U = 3.2905  # one-tailed p < .001, the modal reporting threshold
OUT = "/tmp/claude-0/sdm_input"
os.makedirs(OUT, exist_ok=True)

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
sizes = {str(i): int(n) for i, n in zip(ss.ids, ss.sample_sizes())}
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates
coords["id"] = coords["id"].astype(str)

by_study = {str(k): v for k, v in coords.groupby("id")}
rows = []
for sid in sorted(sizes):
    n = sizes[sid]
    dof = n - 1
    # SDM's threshold column is t, on the same dof as the peaks it censors.
    t_thr = float(z_to_t(U, dof))
    sub = by_study.get(sid)
    safe = sid.replace(" ", "_").replace("/", "_")
    if sub is None or len(sub) == 0:
        open(f"{OUT}/{safe}.no_peaks.txt", "w").close()
        rows.append((safe, n, t_thr, 0))
        continue
    with open(f"{OUT}/{safe}.spm_mni.txt", "w") as handle:
        for r in sub.itertuples():
            t = float(z_to_t(float(r.z_stat), dof))
            handle.write(f"{int(round(r.x))},{int(round(r.y))},{int(round(r.z))},{t:.2f}\n")
    rows.append((safe, n, t_thr, len(sub)))

with open(f"{OUT}/sdm_table.txt", "w") as handle:
    handle.write("study\tn1\tt_thr\n")
    for safe, n, t_thr, _ in rows:
        handle.write(f"{safe}\t{n}\t{t_thr:.4f}\n")

print(f"{len(rows)} studies written to {OUT}")
print(f"{'study':>22s} {'N':>4s} {'t_thr':>7s} {'peaks':>6s}")
for safe, n, t_thr, k in rows:
    print(f"{safe[:22]:>22s} {n:4d} {t_thr:7.2f} {k:6d}")
print(f"total peaks {sum(r[3] for r in rows)}")
