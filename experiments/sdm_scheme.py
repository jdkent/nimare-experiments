"""Score CBES on SDM's own validation scheme, and show what that scheme cannot see.

The AES-SDM validation recreates a *single study's* effect size map from its peaks, compares it
to the map from the raw SPM, and reports the mean square error relative to the ES-SDM default
isotropic kernel (FWHM = 20 mm): "a relative MSE <100% would indicate an improvement".

Three recreations, per study, scored against that study's own g map:

  iso-20mm    the paper's baseline. Each peak contributes ``ES_i * K_i(v)`` with an
              un-normalized Gaussian, and peaks are combined by a K-weighted average, which is
              this script's reading of "combined using a weighted average".
  SDM-PSI     the real thing: the midpoint of the lower/upper bound maps sdm pp wrote.
  CBES        the kernel-weighted peak values scaled by the per-study rho the estimator fits.

A fourth row, ``predict zero``, is not part of their scheme and is the point of including it: a
relative metric cannot tell you whether *either* method beats saying nothing. MSE penalizes bias,
so it is a fair absolute reference in a way correlation would not be.
"""
import os, pickle, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import peak_stat_to_hedges_g
from nimare.meta.utils import get_ale_kernel
from nimare.transforms import ImageTransformer
from nimare.utils import mm2vox

U, ISO_FWHM = 3.2905, 20.0
PP = "/tmp/claude-0/sdm_input/pp"

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = {str(i): float(n) for i, n in zip(ss.ids, ss.sample_sizes())}
truth = {}
for sid, gp in zip(ss.images["id"].astype(str), ss.images["g"]):
    if gp is not None:
        g = masker.transform(str(gp)).ravel()
        truth[sid] = np.nan_to_num(g)

coords = pickle.load(open("/tmp/claude-0/cmp/pain_coords_ss.pkl", "rb"))
table = coords.coordinates
table["id"] = table["id"].astype(str)

# The per-study rho the estimator actually fits, read off a real fit.
fit = CBES(fwhm=None, mask=masker.mask_img, use_images=False, null_method="none",
           threshold=U, peak_bias="per-study").fit(coords)
rho = {str(k): float(v) for k, v in fit.estimator._peak_bias_.items()}
print(f"per-study rho: median {np.median(list(rho.values())):.3f}\n", flush=True)

shape = masker.mask_img.shape[:3]
mask_flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[mask_flat] = np.arange(mask_flat.size)
n_vox = mask_flat.size


def support(fwhm=None, sample_size=None):
    _, kern = (get_ale_kernel(masker.mask_img, fwhm=fwhm) if fwhm is not None
               else get_ale_kernel(masker.mask_img, sample_size=sample_size))
    kern = kern / kern.max()
    half = (np.array(kern.shape) - 1) // 2
    sel = kern >= 0.01
    return np.array(np.nonzero(sel)).T - half, kern[sel]


def recreate(sub, sid, offs, wts, scale):
    """K-weighted combination of peak values over the kernel support."""
    num = np.zeros(n_vox)
    den = np.zeros(n_vox)
    n = sizes[sid]
    z = sub["z_stat"].astype(float).to_numpy()
    g_k, _ = peak_stat_to_hedges_g(np.abs(z), np.full(len(z), n), stat_type="z")
    g_k = np.sign(z) * g_k * scale
    ijk = mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine)
    for centre, value in zip(ijk, g_k):
        vox = centre + offs
        inside = np.all((vox >= 0) & (vox < shape), axis=1)
        cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
        keep = cols >= 0
        cols, w = cols[keep], wts[inside][keep]
        np.add.at(num, cols, value * w * w)   # ES_i * K_i, combined with weight K_i
        np.add.at(den, cols, w)
    return np.divide(num, den, out=np.zeros(n_vox), where=den > 0)


iso_offs, iso_wts = support(fwhm=ISO_FWHM)
ORACLE_RHO = 0.220  # measured in weighting4.py from the same studies
rows = {k: [] for k in ("iso", "sdm", "cbes", "zero", "oracle", "iso_oracle")}
for sid, sub in table.groupby("id"):
    if sid not in truth:
        continue
    t = truth[sid]
    ok = np.isfinite(t)
    rows["zero"].append(np.mean(t[ok] ** 2))
    rows["iso"].append(np.mean((recreate(sub, sid, iso_offs, iso_wts, 1.0)[ok] - t[ok]) ** 2))

    c_offs, c_wts = support(sample_size=sizes[sid])
    rows["cbes"].append(
        np.mean((recreate(sub, sid, c_offs, c_wts, rho.get(sid, 1.0))[ok] - t[ok]) ** 2))
    # What the same recreation scores once the scale is right. ORACLE_RHO is measured from the
    # images (mean |truth| at the peaks over mean reported |g|), so this arm is not available
    # from coordinates alone -- it isolates how much of CBES's MSE is the scale rather than
    # the kernel.
    rows["oracle"].append(
        np.mean((recreate(sub, sid, c_offs, c_wts, ORACLE_RHO)[ok] - t[ok]) ** 2))
    rows["iso_oracle"].append(
        np.mean((recreate(sub, sid, iso_offs, iso_wts, ORACLE_RHO)[ok] - t[ok]) ** 2))

    lo, hi = f"{PP}/{sid}_lower.nii.gz", f"{PP}/{sid}_upper.nii.gz"
    if os.path.exists(lo) and os.path.exists(hi):
        mid = 0.5 * (nib.load(lo).get_fdata() + nib.load(hi).get_fdata())
        vec = masker.transform(
            resample_to_img(nib.Nifti1Image(mid, nib.load(lo).affine), masker.mask_img,
                            interpolation="continuous")).ravel()
        rows["sdm"].append(np.mean((np.nan_to_num(vec)[ok] - t[ok]) ** 2))

base = float(np.mean(rows["iso"]))
print(f"{'recreation':>26s} {'MSE':>9s} {'relative MSE':>14s}")
for key, label in (("iso", "ES-SDM iso 20mm (base)"), ("sdm", "SDM-PSI anisotropic"),
                   ("cbes", "CBES kernel + rho"), ("oracle", "CBES + oracle scale"),
                   ("iso_oracle", "iso 20mm + oracle scale"), ("zero", "predict zero")):
    if not rows[key]:
        continue
    m = float(np.mean(rows[key]))
    print(f"{label:>26s} {m:9.4f} {100 * m / base:13.0f}%", flush=True)
print("\n(<100% is what their scheme calls an improvement; 'predict zero' is not in it)")
