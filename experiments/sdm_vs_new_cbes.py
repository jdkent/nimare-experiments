"""SDM-PSI's coefficient against the redesigned CBES, on the same 21-study pain collection.

Uses the SDM coefficient already on disk from a prior run -- 50 imputations over the same 21
NIDM pain studies, from coordinates built at U = 3.2905 by `make_sdm_input.py` -- so no SDM
re-run is needed. That matters because SDM's coefficient *is* the imputation machinery: its
output lives at `analysis_MyMean/mi/coeff.nii.gz`, above 50 directories each holding a complete
set of imputed per-study g and g_var maps, and `coeff` is their Rubin's-rules pooling. Only PSI
(permutation of subject images) is inference-only. And `pp/` writes `<study>_lower.nii.gz` and
`<study>_upper.nii.gz` -- SDM derives exactly the effect-size bounds CBES's censored likelihood
uses, then imputes *within* them where CBES integrates *over* them. That is the whole difference
between the two methods, and it is what this compares.

The inputs are not identical, and cannot be:

  SDM        21 coordinate tables. It needs nothing else.
  CBES       19 coordinate tables + 2 g/g_var images, because the redesign requires at least
             one image: coordinates carry no magnitude in it.

That asymmetry is the real choice a user faces, not a flaw in the comparison, so it is reported
rather than engineered away.

**The reference has to exclude CBES's two images.** Scoring against a mean of all 21 would score
CBES partly against its own input. So the truth is the inverse-variance mean of the *other 19*
images. SDM is then mildly flattered -- 2 of its 21 coordinate tables came from studies in the
reference -- and CBES is not flattered at all, which is the conservative direction for the claim
being made.

Pattern and magnitude are scored separately because they are separate claims, and magnitude by a
paired median ratio over the reference's top quartile rather than a ratio of means: the
reference's whole-brain mean is near zero and dividing by it manufactures ratios that describe
the denominator.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.image import resample_to_img
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates, d_to_g, t_to_d

U = 3.2905
SDM_COEFF = "/tmp/claude-0/sdm_input/analysis_MyMean/mi/coeff.nii.gz"
N_IMAGES = int(os.environ.get("NIMAGES", 2))
WORKDIR = f"/tmp/claude-0/sdmcmp_{os.getpid()}"
os.makedirs(WORKDIR, exist_ok=True)


def g_and_var(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    g = d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)
    return g, 1.0 / n + g**2 / (2.0 * n)


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    quartile = t >= np.percentile(t, 75)
    ratio = float(np.median(a[quartile] / np.maximum(t[quartile], 1e-6)))
    return (float(stats.pearsonr(a, t)[0]), float(stats.spearmanr(a, t)[0]), float(auc),
            ratio, float(np.median(a[quartile])), float(np.median(t[quartile])))


if __name__ == "__main__":
    ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
    masker = ss.masker
    ids = [str(i) for i in ss.images["id"]]
    sizes = {str(i): int(np.atleast_1d(n).ravel()[0]) for i, n in zip(ss.ids, ss.sample_sizes())}

    zmaps, gmaps, vmaps = {}, {}, {}
    for sid, zp, gp, vp in zip(ids, ss.images["z"], ss.images["g"], ss.images["g_var"]):
        if zp is None or gp is None:
            continue
        zmaps[sid] = np.nan_to_num(masker.transform(str(zp)).ravel())
        gmaps[sid] = np.nan_to_num(masker.transform(str(gp)).ravel())
        vmaps[sid] = np.nan_to_num(masker.transform(str(vp)).ravel(), nan=np.inf)

    usable = sorted(zmaps)
    donors = usable[:N_IMAGES]
    held = [s for s in usable if s not in donors]
    num = den = 0.0
    for sid in held:
        w = 1.0 / np.maximum(vmaps[sid], 1e-6)
        num = num + w * gmaps[sid]
        den = den + w
    truth = np.abs(num / np.maximum(den, 1e-12))

    # Exactly SDM's coordinates: the same extraction at the same threshold.
    coords = ImagesToCoordinates(
        merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
    ).transform(ss).coordinates
    coords["id"] = coords["id"].astype(str)
    by_study = {str(k): v for k, v in coords.groupby("id")}

    def write(values, name):
        path = f"{WORKDIR}/{name}.nii.gz"
        nib.save(masker.inverse_transform(values), path)
        return path

    studies = []
    for sid in usable:
        meta = {"sample_sizes": [sizes.get(sid, 25)]}
        analysis = {"id": sid, "name": "1", "metadata": meta, "points": [], "images": []}
        if sid in donors:
            analysis["images"] = [
                {"url": write(gmaps[sid], f"{sid}_g"), "filename": "g",
                 "space": "MNI", "value_type": "g"},
                {"url": write(np.minimum(vmaps[sid], 1e6), f"{sid}_v"), "filename": "v",
                 "space": "MNI", "value_type": "g_var"},
            ]
        else:
            sub = by_study.get(sid)
            if sub is None or not len(sub):
                continue  # reported nothing: still on the roster, contributing silence
            analysis["points"] = [
                {"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                 "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                for r in sub.itertuples()
            ]
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [analysis]})

    print(f"NIDM pain: {len(usable)} studies with images; CBES gets {len(donors)} of them as "
          f"images and {len(studies) - len(donors)} as coordinate tables")
    print(f"reference: inverse-variance mean of the {len(held)} images CBES never saw\n")

    result = CBES(mask=masker, null_method="none").fit(
        Studyset({"id": "cmp", "name": "cmp", "studies": studies},
                 target=None, mask=masker.mask_img))
    have = set(result.maps)
    cbes_g = np.abs(result.get_map("g", return_type="array").ravel())
    cbes_m = (np.abs(result.get_map("g_marginal", return_type="array").ravel())
              if "g_marginal" in have else None)
    covered = result.get_map("n_studies", return_type="array").ravel() > 0

    sdm_img = resample_to_img(nib.load(SDM_COEFF), masker.mask_img, interpolation="continuous")
    sdm = np.abs(masker.transform(sdm_img).ravel())

    # Two donors pooled by inverse variance: what a user would do with the images alone.
    w = sum(1.0 / np.maximum(vmaps[s], 1e-6) for s in donors)
    images_only = np.abs(
        sum(gmaps[s] / np.maximum(vmaps[s], 1e-6) for s in donors) / np.maximum(w, 1e-12))

    use = covered & np.isfinite(truth) & np.isfinite(sdm) & (sdm != 0) & np.isfinite(cbes_g)
    print(f"scoring over {int(use.sum())} voxels both methods cover\n")
    print(f"{'estimate':>26} {'r':>7} {'rank r':>7} {'AUC':>6} {'mag ratio':>10} "
          f"{'|est| top':>10} {'|ref| top':>10}")
    arms = [("SDM-PSI coeff (21 tables)", sdm),
            ("images only (2 images)", images_only),
            ("CBES g (2 img + 19 tab)", cbes_g)]
    if cbes_m is not None:
        arms.append(("CBES g_marginal", cbes_m))
    for label, est in arms:
        r, rho, auc, ratio, top_est, top_ref = score(est, truth, use)
        print(f"{label:>26} {r:+7.3f} {rho:+7.3f} {auc:6.3f} {ratio:10.2f} "
              f"{top_est:10.3f} {top_ref:10.3f}")
