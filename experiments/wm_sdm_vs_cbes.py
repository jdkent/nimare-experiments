"""CBES against SDM-PSI on a working-memory studyset drawn from NeuroStore.

Pain and HCP are the two beds the headline claim rests on, and both are narrow in the same way:
pain is one collection curated by hand, HCP is one task split into synthetic studies. Neither
answers whether the estimator works on a studyset assembled the way a meta-analyst would
actually assemble one -- by asking a literature index for the studies on a topic and taking
whatever maps they happen to have shared.

`fetch_working_memory.py` does that asking. This bed runs both estimators on the result.

Parity with the other beds: the same studies supply images, the rest supply coordinates, and the
reference is the half of the corpus neither estimator saw, pooled inverse-variance. Coordinates
are extracted with `report_peaks`, so each study's table is whatever survives its own correction
and a study with nothing surviving drops out. Each study hands CBES its own real height through
`reporting_threshold` rather than letting the estimator infer one, and hands SDM the same height
converted to t.

Unlike pain there are no transcribed tables here, so `extracted` is the only option -- which
makes this the bed that measures what the extraction scheme costs on real, uncurated data.

One property of a real studyset that neither pain nor HCP has: the sample sizes span 21 to 1369,
so a single study carries more inverse-variance weight than all the others together. When that
study is held out the reference is very nearly that study alone; when it is in the working half
the reference is a real pool of typical studies, but a much noisier one. Those are two different
questions, so every summary is also reported split by which happened. Averaging over them reports
a number that describes neither.
"""
import os, sys, subprocess, json, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma.effectsize import CBES, DEFAULT_REPORT_RADIUS_MM
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, z_to_t
from reporting import report_peaks
from sdm_params import PARAMS

CORPUS = "/tmp/claude-0/wm"
N_IMAGES = int(os.environ.get("NIMAGES", 1))
N_SPLITS = int(os.environ.get("NSPLITS", 6))
SCHEME = os.environ.get("SCHEME", "cluster")
FOCUS = "max"
RUN_SDM = os.environ.get("SDM", "1") == "1"
#: Radius over which a report asserts its lower bound; None is the named voxel alone.
REPORT_RADIUS = (None if os.environ.get("REPORTR", "default") in ("none", "None")
                 else (float(os.environ["REPORTR"]) if os.environ.get("REPORTR", "default")
                       not in ("default", "") else DEFAULT_REPORT_RADIUS_MM))

CURATE = os.environ.get("CURATE", "0") == "1"

# Map names that mark a map as something other than a group activation contrast. A searchlight
# decoding accuracy map, a white-matter FA association and a resting fALFF difference are all
# "related to working memory" in the index's sense, but they are not the same estimand as a task
# contrast, so a reference pooled across them is pooling different quantities. Excluding them is
# a judgement about what a map measures, not about how strong it is, so it leaves untouched the
# rule that the peak count is whatever survives a correction.
NON_ACTIVATION = ("searchlight", "decoding", "falff", "alff", "network", "microstructure",
                  "adjusted for", " fa ", "fa,")
SDM_HOME = "/tmp/claude-0/sdm/SdmPsiGui-linux64-v6.23"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
WORK = f"/tmp/claude-0/wmsdm_{os.getpid()}"


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


def sdm_name(prefix, index):
    """A study label SDM can look up: letters only, no digits.

    `sdm_parse` fails to find a study's `.no_peaks.txt` when the study name carries two or more
    digits -- `c16` and `c99` are reported missing while the file sits next to the table, but
    `c3`, `c7` and any digit-free name are found. Studies that report peaks are unaffected, so
    the bug only bites the silent studies, which is exactly the arm this bed needs. Encoding the
    index in letters sidesteps it, and is applied to every study in every arm so nothing about
    the comparison depends on it.
    """
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(ord("a") + rem) + letters
    return f"{prefix}{letters}"


def write_sdm_t_map(t_values, n, path):
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
    maps = np.load(f"{CORPUS}/wm_z.npy").astype(np.float64)
    meta = json.load(open(f"{CORPUS}/wm_meta.json"))
    if CURATE:
        keep = [i for i, r in enumerate(meta)
                if not any(term in r["name"].lower() for term in NON_ACTIVATION)]
        print(f"curated to {len(keep)} activation contrasts of {len(meta)}", flush=True)
        maps, meta = maps[keep], [meta[i] for i in keep]
    sizes = [int(r["n"]) for r in meta]
    total = len(maps)
    print(f"working memory: {total} studies, {N_IMAGES} image(s), {SCHEME}/{FOCUS} reporting, "
          f"{N_SPLITS} splits", flush=True)
    print(f"sample sizes {min(sizes)}-{max(sizes)}, median {sorted(sizes)[total // 2]}",
          flush=True)

    # How coherent is an uncurated studyset? A bed whose maps do not agree with each other
    # cannot separate two estimators, because the held-out reference is then mostly noise.
    corr = np.corrcoef(maps)
    off = corr[np.triu_indices(total, 1)]
    print(f"pairwise spatial r between study maps: median {np.median(off):+.3f}, "
          f"{np.mean(off > 0) * 100:.0f}% positive\n", flush=True)

    rng = np.random.default_rng(0)
    rows, strata, skipped = {}, {}, 0
    for split in range(N_SPLITS):
        order = rng.permutation(total)
        half = total // 2
        work, held = order[:half], order[half:]
        truth = np.abs(pooled(held, maps, sizes))
        image_members, table_members = list(work[:N_IMAGES]), list(work[N_IMAGES:])
        # Inverse-variance weight goes as the sample size, so one very large study can make the
        # reference a single study wearing a pool's clothes. Report that before believing a split.
        held_n = np.array([sizes[i] for i in held], dtype=float)
        share = held_n.max() / held_n.sum()

        studies, sdm_rows = [], []
        for i in image_members:
            g, var = g_and_var(maps[i], sizes[i])
            os.makedirs(f"{WORK}/cbes", exist_ok=True)
            label = sdm_name("i", int(i))
            gp, vp = f"{WORK}/cbes/g{i}.nii.gz", f"{WORK}/cbes/v{i}.nii.gz"
            nib.save(masker.inverse_transform(g), gp)
            nib.save(masker.inverse_transform(var), vp)
            # An image study declares no reporting threshold: it reported everything.
            m = {"sample_sizes": [int(sizes[i])], "reporting_threshold": float("nan")}
            studies.append({"id": label, "name": label, "metadata": m, "analyses": [
                {"id": label, "name": "1", "metadata": m, "points": [], "images": [
                    {"url": gp, "filename": "g", "value_type": "g", "space": "MNI"},
                    {"url": vp, "filename": "v", "value_type": "g_var", "space": "MNI"}]}]})
            t = np.sign(maps[i]) * np.abs(stats.t.isf(
                stats.norm.sf(np.abs(maps[i])), sizes[i] - 1))
            sdm_rows.append((label, sizes[i], 0.0, np.nan_to_num(t), []))

        reported = 0
        for i in table_members:
            n = int(sizes[i])
            peaks, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                         scheme=SCHEME, focus=FOCUS)
            found = [(nib.affines.apply_affine(affine, np.asarray(ijk, float)), float(v))
                     for ijk, v in peaks]
            label = sdm_name("c", int(i))
            sdm_rows.append((label, n, float(z_to_t(np.array([height]), n - 1)[0]),
                             None, found))
            if not found:
                continue
            reported += 1
            m = {"sample_sizes": [n], "reporting_threshold": float(height)}
            studies.append({"id": label, "name": label, "metadata": m, "analyses": [
                {"id": label, "name": "1", "metadata": m, "points": [
                    {"space": "MNI", "coordinates": [float(v) for v in xyz],
                     "values": [{"kind": "Z", "value": float(value)}]}
                    for xyz, value in found]}]})

        npeaks = sum(len(r[4]) for r in sdm_rows)
        print(f"  split {split + 1}: {reported}/{len(table_members)} tables reported, "
              f"{npeaks} peaks, {len(held)} studies held out "
              f"(largest holds {share * 100:.0f}% of the reference weight, "
              f"image donor n={max(sizes[i] for i in image_members)})", flush=True)
        if reported == 0:
            skipped += 1
            print("    no table reported anything, skipped", flush=True)
            continue

        est = CBES(mask=masker, null_method="none", threshold="reporting_threshold",
                   report_radius=REPORT_RADIUS)
        result = est.fit(Studyset({"id": "wm", "name": "wm", "studies": studies},
                                  target=None, mask=mask_img))
        cbes_g = np.abs(result.get_map("g", return_type="array").ravel())
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        # The estimator's own answer to "did the coordinates matter here": the share of the
        # likelihood weight at a voxel that came from the silence and peak terms rather than
        # from the image. If this is near zero, g is the donated image and nothing else.
        try:
            cshare = result.get_map("coordinate_share", return_type="array").ravel()
            strong = cbes_g >= np.nanpercentile(cbes_g[covered], 90)
            print(f"    coordinate share: median {np.nanmedian(cshare[covered]):.3f}, "
                  f"at the strongest decile {np.nanmedian(cshare[covered & strong]):.3f}",
                  flush=True)
        except Exception:
            pass
        only = np.abs(pooled(image_members, maps, sizes))
        sdm = run_sdm(sdm_rows, f"{WORK}/sdm_{split}") if RUN_SDM else None

        use = (covered & np.isfinite(cbes_g) & (cbes_g > 0)
               & np.isfinite(truth) & np.isfinite(only))
        if sdm is not None:
            use = use & np.isfinite(sdm)
        arms = [("images only", only), ("CBES g", cbes_g)]
        if sdm is not None:
            arms.append(("SDM-PSI", sdm))
        # The regime that matters is where the one very large study landed: held out, it makes
        # the reference almost a single study; in the working half, the reference is a genuine
        # pool of typical studies but a weaker one.
        regime = "reference dominated" if share > 0.5 else "reference pooled"
        for name, est_map in arms:
            rows.setdefault(name, []).append(score(est_map, truth, use))
            strata.setdefault((regime, name), []).append(score(est_map, truth, use))
        print(f"    {int(use.sum())} voxels scored, "
              f"sdm {'ok' if sdm is not None else ('FAILED' if RUN_SDM else 'skipped')}",
              flush=True)

    print(f"\n{'estimate':>14} {'r':>7} {'rank r':>8} {'AUC':>7} {'magnitude':>10} {'n':>4}")
    for name in ("images only", "CBES g", "SDM-PSI"):
        if name not in rows:
            continue
        arr = np.array(rows[name])
        m = np.mean(arr, axis=0)
        print(f"{name:>14} {m[0]:+7.3f} {m[1]:+8.3f} {m[2]:7.3f} {m[3]:10.3f} {len(arr):4d}")
    for regime in ("reference pooled", "reference dominated"):
        present = [k for k in strata if k[0] == regime]
        if not present:
            continue
        print(f"\n{regime}:")
        for name in ("images only", "CBES g", "SDM-PSI"):
            if (regime, name) not in strata:
                continue
            arr = np.array(strata[(regime, name)])
            m = np.mean(arr, axis=0)
            print(f"{name:>14} {m[0]:+7.3f} {m[1]:+8.3f} {m[2]:7.3f} {m[3]:10.3f} "
                  f"{len(arr):4d}")
    if skipped:
        print(f"\n{skipped} of {N_SPLITS} splits had no reporting table at all")
