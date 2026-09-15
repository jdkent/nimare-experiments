"""CBES against SDM-PSI on a NeuroVault collection, which neither had been compared on.

SDM-PSI had been run here on HCP and on NIDM pain only. jdkent asked whether a NeuroVault-based
studyset had been covered, and it had not. This closes that: the cached "animal" collection --
11 movie-watching studies sharing an as-Animal contrast, each with a z map and derived g/g_var.

Two things make this the weakest of the three beds, and both are worth stating before the numbers
rather than after:

  * **11 studies**, so a half is five or six. That is a thin roster for a split-half design and
    the split-to-split spread will be wide.
  * **No published tables.** NeuroVault carries maps, not the papers' coordinate tables, so the
    tables are extracted from each study's own z map with `report_peaks`. The audit found that
    extraction recovers only 23% of a paper's published peaks within 8 mm, so these are *not*
    what a literature would have handed the estimator. Every arm sees the same extracted tables,
    so the comparison between estimators is fair; what it cannot speak to is behaviour on real
    transcribed coordinates.

The collection also does not record per-study sample sizes, so one value is assumed for every
study. That makes the z-to-g conversion uniform across the roster rather than accidentally
informative -- and it puts this bed in the *unidentified* regime for the prevalence, since only
sample-size spread separates pi from mu.

No cap on peaks: each study reports whatever survives its own correction.
"""
import os, sys, glob, subprocess, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma.effectsize import CBES
from nimare.studyset import Studyset
from nimare.transforms import z_to_t
from reporting import report_peaks
from sdm_params import PARAMS

SCHEME, FOCUS = "cluster", "max"
N_IMAGES = int(os.environ.get("NIMAGES", 1))
N_SPLITS = int(os.environ.get("NSPLITS", 6))
HOME = os.path.expanduser("~/.nimare")
ASSUMED_N = 20
ASSUMED_Z = 3.09
SDM_HOME = "/tmp/claude-0/sdm/SdmPsiGui-linux64-v6.23"
WORK = f"/tmp/claude-0/animalsdm_{os.getpid()}"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def load(pattern):
    out = {}
    for path in sorted(glob.glob(f"{HOME}/{pattern}")):
        study = os.path.basename(path).split("-")[1]
        img = resample_to_img(nib.load(path), mask_img, interpolation="continuous",
                              force_resample=True, copy_header=True)
        out[study] = np.nan_to_num(masker.transform(img).ravel())
    return out


def pooled(members, g_maps, var_maps):
    stack = np.array([g_maps[i] for i in members])
    weights = np.array([1.0 / np.maximum(var_maps[i], 1e-9) for i in members])
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def write_sdm_t_map(t_values, n, path):
    """A study's t map, headered so `sdm_parse` accepts it rather than dropping it silently."""
    template = nib.load(f"{SDM_HOME}/share/sdm_template.nii.gz")
    img = resample_to_img(masker.inverse_transform(t_values), template,
                          interpolation="continuous", force_resample=True, copy_header=True)
    out = nib.Nifti1Image(np.nan_to_num(np.asarray(img.dataobj, dtype=np.float32)),
                          template.affine)
    out.header.set_intent("t test", (float(n) - 1.0,), name="")
    out.header.set_sform(template.affine, code=4)
    out.header.set_qform(template.affine, code=2)
    nib.save(out, path)


def run_sdm(rows, workdir):
    os.makedirs(workdir, exist_ok=True)
    for name, n, t_thr, t_map, found in rows:
        if t_map is not None:
            write_sdm_t_map(t_map, n, f"{workdir}/{name}.nii.gz")
            continue
        if not found:
            open(f"{workdir}/{name}.no_peaks.txt", "w").close()
            continue
        with open(f"{workdir}/{name}.spm_mni.txt", "w") as handle:
            for xyz, value in found:
                t = float(z_to_t(np.array([abs(value)]), n - 1)[0])
                handle.write(f"{int(round(xyz[0]))},{int(round(xyz[1]))},"
                             f"{int(round(xyz[2]))},{t:.2f}\n")
    with open(f"{workdir}/sdm_table.txt", "w") as handle:
        handle.write("study\tn1\tt_thr\n")
        for name, n, t_thr, _, _ in rows:
            handle.write(f"{name}\t{n}\t{t_thr:.4f}\n")
    with open(f"{workdir}/sdmpsi_params.xml", "w") as handle:
        handle.write(PARAMS.format(n_studies=len(rows)))
    for command in ("pp", "mi"):
        subprocess.run([f"{SDM_HOME}/bin/linux64/sdm_parse", command], cwd=workdir,
                       capture_output=True, text=True, timeout=5400)
    coeff = f"{workdir}/analysis_MyMean/mi/coeff.nii.gz"
    if not os.path.exists(coeff):
        return None
    return np.abs(masker.transform(resample_to_img(
        nib.load(coeff), mask_img, interpolation="continuous",
        force_resample=True, copy_header=True)).ravel())


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
           if n_pos and n_neg else np.nan)
    quartile = t >= np.percentile(t, 75)
    return (float(stats.pearsonr(a, t)[0]), float(stats.spearmanr(a, t)[0]), float(auc),
            float(np.median(a[quartile] / np.maximum(t[quartile], 1e-6))))


if __name__ == "__main__":
    z_maps = load("study-*-animal_*_z.nii.gz")
    g_maps = load("study-*-animal_*[0-9]_g.nii.gz")
    var_maps = load("study-*-animal_*_g_var.nii.gz")
    ids = sorted(set(z_maps) & set(g_maps) & set(var_maps))
    print(f"NeuroVault animal: {len(ids)} studies, {N_IMAGES} image(s), {N_SPLITS} splits, "
          f"assumed N = {ASSUMED_N}, tables extracted with {SCHEME}/{FOCUS}\n", flush=True)

    rng = np.random.default_rng(0)
    rows = {}
    for split in range(N_SPLITS):
        order = list(rng.permutation(ids))
        half = len(ids) // 2
        work, held = order[:half], order[half:]
        truth = np.abs(pooled(held, g_maps, var_maps))
        images, tables = work[:N_IMAGES], work[N_IMAGES:]

        studies, sdm_rows = [], []
        for i in images:
            os.makedirs(f"{WORK}/cbes", exist_ok=True)
            gp, vp = f"{WORK}/cbes/g{i}.nii.gz", f"{WORK}/cbes/v{i}.nii.gz"
            nib.save(masker.inverse_transform(g_maps[i]), gp)
            nib.save(masker.inverse_transform(var_maps[i]), vp)
            meta = {"sample_sizes": [ASSUMED_N]}
            studies.append({"id": f"i{i}", "name": f"i{i}", "metadata": meta, "analyses": [
                {"id": f"i{i}", "name": "1", "metadata": meta, "points": [], "images": [
                    {"url": gp, "filename": "g", "value_type": "g", "space": "MNI"},
                    {"url": vp, "filename": "v", "value_type": "g_var", "space": "MNI"}]}]})
            t = np.sign(z_maps[i]) * np.abs(stats.t.isf(
                stats.norm.sf(np.abs(z_maps[i])), ASSUMED_N - 1))
            sdm_rows.append((f"i{i}", ASSUMED_N,
                             float(z_to_t(np.array([ASSUMED_Z]), ASSUMED_N - 1)[0]),
                             np.nan_to_num(t), []))

        for i in tables:
            peaks, height = report_peaks(z_maps[i], mask_bool, shape, zooms,
                                         scheme=SCHEME, focus=FOCUS)
            found = [(nib.affines.apply_affine(affine, np.asarray(ijk, float)), float(v))
                     for ijk, v in peaks]
            sdm_rows.append((f"c{i}", ASSUMED_N,
                             float(z_to_t(np.array([height]), ASSUMED_N - 1)[0]), None, found))
            if not found:
                continue
            meta = {"sample_sizes": [ASSUMED_N], "reporting_threshold": float(height)}
            studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
                {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                    {"space": "MNI", "coordinates": [float(v) for v in xyz],
                     "values": [{"kind": "Z", "value": float(value)}]}
                    for xyz, value in found]}]})

        if len(studies) < 2:
            print(f"  split {split + 1}: too few studies reported, skipped", flush=True)
            continue
        result = CBES(mask=masker, null_method="none",
                      threshold="reporting_threshold").fit(
            Studyset({"id": "a", "name": "a", "studies": studies}, target=None, mask=mask_img))
        cbes_g = np.abs(result.get_map("g", return_type="array").ravel())
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        only = np.abs(pooled(images, g_maps, var_maps))
        sdm = run_sdm(sdm_rows, f"{WORK}/sdm_{split}")

        use = (covered & np.isfinite(cbes_g) & (cbes_g > 0) & np.isfinite(truth)
               & np.isfinite(only))
        if sdm is not None:
            use = use & np.isfinite(sdm)
        arms = [("images only", only), ("CBES g", cbes_g)]
        if sdm is not None:
            arms.append(("SDM-PSI", sdm))
        for name, est_map in arms:
            rows.setdefault(name, []).append(score(est_map, truth, use))
        print(f"  split {split + 1}: {int(use.sum())} voxels, "
              f"sdm {'ok' if sdm is not None else 'FAILED'}", flush=True)

    print(f"\n{'estimate':>14} {'r':>7} {'rank r':>8} {'AUC':>7} {'magnitude':>10} {'n':>3}")
    for name in ("images only", "CBES g", "SDM-PSI"):
        if name not in rows:
            continue
        arr = np.array(rows[name])
        m = arr.mean(axis=0)
        print(f"{name:>14} {m[0]:+7.3f} {m[1]:+8.3f} {m[2]:7.3f} {m[3]:10.3f} "
              f"{len(arr):3d}")
