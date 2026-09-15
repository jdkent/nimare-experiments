"""CBES against SDM-PSI on the NIDM pain collection, at one image and at two.

SDM-PSI had only ever been run here on HCP, and only with two images. That leaves the question
jdkent actually asked unanswered on the dataset that matters most -- pain is the collection the
headline claim rests on, and it is the only one carrying coordinates a human transcribed from a
paper alongside the maps.

Both estimators get the same parity they got on HCP: the same studies supply images, the rest
supply tables, and the reference is the half of the collection neither one saw. ES-SDM was built
for exactly this mixture, so supplying it images is not a courtesy but its intended input.

Two things this bed can do that the HCP one cannot:

  * **the tables are the published ones.** `TABLES=published` uses the collection's own 267
    transcribed peaks; `extracted` re-derives them with `report_peaks`. The audit found those
    differ sharply -- the cluster scheme recovers 23% of published peaks within 8 mm -- so the
    comparison is run on the real tables by default.
  * **the image count can go to one**, which is the configuration under question.

Scored against an inverse-variance pool of the held-out half, which is a reference built from
studies that contributed no coordinate. Peaks are never capped: a published table is whatever the
paper printed, and an extracted one is whatever survives its own correction.
"""
import os, sys, subprocess, warnings; warnings.simplefilter("ignore")
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
from nimare.transforms import d_to_g, t_to_d, z_to_t
from load_pain import load_pain
from reporting import report_peaks
from sdm_params import PARAMS

TABLES = os.environ.get("TABLES", "published")
N_IMAGES = int(os.environ.get("NIMAGES", 1))
N_SPLITS = int(os.environ.get("NSPLITS", 6))
SCHEME, FOCUS = "cluster", "max"
SDM_HOME = "/tmp/claude-0/sdm/SdmPsiGui-linux64-v6.23"
ASSUMED_Z = 3.09

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
WORK = f"/tmp/claude-0/painsdm_{os.getpid()}"


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def g_and_var(z, n):
    g = to_g(z, n)
    return g, 1.0 / n + g**2 / (2.0 * n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g, var = g_and_var(maps[i], sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(var, 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def write_sdm_t_map(t_values, n, path):
    """A study's t map, headered so `sdm_parse` accepts it rather than silently dropping it."""
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
    """`rows` is [(name, n, t_thr, t_map_or_None, [(xyz, z), ...])]; returns |coeff| or None."""
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
        done = subprocess.run([f"{SDM_HOME}/bin/linux64/sdm_parse", command],
                              cwd=workdir, capture_output=True, text=True, timeout=5400)
        print(f"    sdm {command}: "
              f"{' '.join(done.stdout.strip().splitlines()[-1:])[:100]}", flush=True)
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
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    quartile = t >= np.percentile(t, 75)
    return (float(stats.pearsonr(a, t)[0]), float(stats.spearmanr(a, t)[0]), float(auc),
            float(np.median(a[quartile] / np.maximum(t[quartile], 1e-6))))


if __name__ == "__main__":
    ss = load_pain()
    maps, sizes, study_ids = [], [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
        study_ids.append(row.study_id)
    maps = np.array(maps)
    frame = ss.coordinates
    by_study = {sid: g[["x", "y", "z"]].astype(float).to_numpy()
                for sid, g in frame.groupby("study_id")}
    total = len(maps)
    print(f"NIDM pain: {total} studies, tables={TABLES}, {N_IMAGES} image(s), "
          f"{N_SPLITS} splits\n", flush=True)

    rng = np.random.default_rng(0)
    rows = {}
    for split in range(N_SPLITS):
        order = rng.permutation(total)
        half = total // 2
        work, held = order[:half], order[half:]
        truth = np.abs(pooled(held, maps, sizes))
        image_members, table_members = list(work[:N_IMAGES]), list(work[N_IMAGES:])

        studies, sdm_rows, donor_g, donor_v = [], [], [], []
        for i in image_members:
            g, var = g_and_var(maps[i], sizes[i])
            donor_g.append(g)
            donor_v.append(var)
            os.makedirs(f"{WORK}/cbes", exist_ok=True)
            gp, vp = f"{WORK}/cbes/g{i}.nii.gz", f"{WORK}/cbes/v{i}.nii.gz"
            nib.save(masker.inverse_transform(g), gp)
            nib.save(masker.inverse_transform(var), vp)
            meta = {"sample_sizes": [int(sizes[i])]}
            studies.append({"id": f"i{i}", "name": f"i{i}", "metadata": meta, "analyses": [
                {"id": f"i{i}", "name": "1", "metadata": meta, "points": [], "images": [
                    {"url": gp, "filename": "g", "value_type": "g", "space": "MNI"},
                    {"url": vp, "filename": "v", "value_type": "g_var", "space": "MNI"}]}]})
            t = np.sign(maps[i]) * np.abs(stats.t.isf(
                stats.norm.sf(np.abs(maps[i])), sizes[i] - 1))
            sdm_rows.append((f"i{i}", sizes[i],
                             float(z_to_t(np.array([ASSUMED_Z]), sizes[i] - 1)[0]),
                             np.nan_to_num(t), []))

        for i in table_members:
            n = int(sizes[i])
            if TABLES == "published":
                pts = by_study.get(study_ids[i])
                found = [(xyz, float(ASSUMED_Z)) for xyz in (pts if pts is not None else [])]
                height = ASSUMED_Z
            else:
                peaks, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                             scheme=SCHEME, focus=FOCUS)
                found = [(nib.affines.apply_affine(affine, np.asarray(ijk, float)), float(v))
                         for ijk, v in peaks]
            sdm_rows.append((f"c{i}", n, float(z_to_t(np.array([height]), n - 1)[0]),
                             None, found))
            if not found:
                continue
            meta = {"sample_sizes": [n]}
            studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
                {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                    {"space": "MNI", "coordinates": [float(v) for v in xyz],
                     "values": [{"kind": "Z", "value": float(value)}]}
                    for xyz, value in found]}]})

        est = CBES(mask=masker, null_method="none", threshold=ASSUMED_Z)
        result = est.fit(Studyset({"id": "p", "name": "p", "studies": studies},
                                  target=None, mask=mask_img))
        cbes_g = np.abs(result.get_map("g", return_type="array").ravel())
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        only = np.abs(pooled(image_members, maps, sizes))
        sdm = run_sdm(sdm_rows, f"{WORK}/sdm_{split}")

        use = covered & np.isfinite(cbes_g) & (cbes_g > 0) & np.isfinite(truth) & np.isfinite(only)
        if sdm is not None:
            use = use & np.isfinite(sdm)
        arms = [("images only", only), ("CBES g", cbes_g)]
        if sdm is not None:
            arms.append(("SDM-PSI", sdm))
        for name, est_map in arms:
            rows.setdefault(name, []).append(score(est_map, truth, use))
        print(f"  split {split + 1}: {int(use.sum())} voxels scored, "
              f"sdm {'ok' if sdm is not None else 'FAILED'}", flush=True)

    print(f"\n{'estimate':>14} {'r':>7} {'rank r':>8} {'AUC':>7} {'magnitude':>10}")
    for name in ("images only", "CBES g", "SDM-PSI"):
        if name not in rows:
            continue
        m = np.mean(np.array(rows[name]), axis=0)
        print(f"{name:>14} {m[0]:+7.3f} {m[1]:+8.3f} {m[2]:7.3f} {m[3]:10.3f}")
