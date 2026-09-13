"""How much does the calibrated scale depend on the assumed smoothness?

The regional censoring term needs the FWHM of the studies' statistic maps, which coordinates do
not carry. If the calibrated rho barely moves with it, the assumption is harmless. If it decides
whether a solution exists at all, then one unknown constant has been traded for another -- and
estimating that one from the same peak counts that set the calibration target would be circular,
manufacturing an interior solution that merely looks like success.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import _coverage_resels, _ec_peak, _rft_censoring_terms
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates
ids = sorted(set(str(i) for i in coords["id"]))

studies = []
for sid in ids:
    sub = coords[coords["id"].astype(str) == sid]
    meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
    studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [{
        "id": sid, "name": "1", "metadata": meta,
        "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                    "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                   for r in sub.itertuples()]}]})
cs = Studyset({"id": "pain", "name": "pain", "studies": studies}, target=None,
              mask=masker.mask_img)

print("predicted reporting rate at a fixed effect, by assumed smoothness")
print("(observed on this collection: 0.677)\n")
print(f"{'FWHM':>6s} {'resels':>9s} {'g=0.1':>7s} {'g=0.2':>7s} {'g=0.3':>7s} {'g=0.5':>7s}")
sqrt_n = np.array([np.sqrt(sizes[s]) for s in ids])
for fwhm in (6.0, 8.0, 10.0, 12.0, 14.0):
    resels = _coverage_resels(20.0, fwhm)
    peak = _ec_peak(resels)
    row = [f"{fwhm:6.1f}", f"{resels[2]:9.1f}"]
    for g in (0.1, 0.2, 0.3, 0.5):
        p = _rft_censoring_terms(np.full(len(ids), g), np.full(len(ids), U), sqrt_n, resels, peak)
        row.append(f"{(1.0 - p['prob']).mean():7.3f}")
    print(" ".join(row), flush=True)

print(f"\n{'FWHM':>6s} {'calibrated rho':>15s}  note")
for fwhm in (8.0, 10.0, 12.0, 14.0):
    est = CBES(fwhm=10.0, null_method="parametric", threshold="reporting_threshold",
               peak_bias="per-study", peak_bias_scale="rate-match", use_images=False,
               censoring="rft", smoothness_fwhm=fwhm, mask=masker)
    est.fit(cs)
    got = est._peak_bias_scale_
    note = "railed, fell back" if got == 1.0 else "interior solution"
    print(f"{fwhm:6.1f} {got:15.3f}  {note}", flush=True)
print("\nimage-calibrated reference on this collection: 0.33-0.48")
