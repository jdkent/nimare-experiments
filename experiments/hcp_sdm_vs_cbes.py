"""SDM-PSI against the redesigned CBES, refereed by subjects neither of them saw.

Every earlier version of this comparison was confounded. The reference was built from the same
maps the coordinates were extracted from, so both methods were scored partly against their own
input -- and asymmetrically, since SDM imputes effect mass at the peaks of the reference studies
while CBES reads those tables only as indicators. HCP fixes it, because it has per-subject maps:

  * `used` subjects are cut into `n_studies` synthetic studies of `n_per_study` each. Each study's
    own one-sample map is computed and its coordinates extracted the way a paper would produce
    them (multiplicity-corrected, whole clusters, one focus per cluster, 8 mm apart, no cap).
  * `held` subjects -- never used to make any study -- give the truth at every voxel.

The selection that produces the peaks is then statistically independent of the quantity being
compared against. Prevalence is 1 by construction (every study draws from one population), so
mu is the marginal and the pi/mu split is not under test here.

Three arms on the same studies:

  SDM-PSI       all `n_studies` coordinate tables. Its coefficient is the Rubin's-rules pooling
                over 50 imputations -- `analysis_MyMean/mi/coeff.nii.gz` sits above i00000..i00049,
                each a full set of imputed per-study maps. Only PSI is inference-only. And `pp/`
                writes each study's `_lower`/`_upper` effect-size bounds: SDM derives exactly the
                bounds CBES's censored likelihood uses, then imputes *within* them where CBES
                integrates *over* them.
  SDM, parity   the same `n_images` studies supplied as t maps instead of peak files -- ES-SDM
                takes a mixture by design -- so neither method is handed more studies-as-maps
                than the other.
  CBES          `n_images` studies contribute their g/g_var maps, the rest their coordinate
                tables. The redesign requires at least one image.
  images only   those same `n_images` maps pooled by inverse variance, i.e. what a user would do
                instead of running CBES at all.

With the images given to both, the only remaining asymmetry is the form each reads them in --
t maps for SDM, g/g_var for CBES -- which is each method's own native input for the same
underlying data.
"""
import os, sys, glob, subprocess, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z, z_to_t
from reporting import report_peaks
from sdm_params import PARAMS

CONTRAST = os.environ.get("CONTRAST", "MOTOR_LH")
SCHEME = os.environ.get("SCHEME", "cluster")
FOCUS = os.environ.get("FOCUS", "max")
N_PER_STUDY = int(os.environ.get("NPER", 30))
N_STUDIES = int(os.environ.get("NSTUDIES", 16))
N_IMAGES = int(os.environ.get("NIMAGES", 2))
N_SPLITS = int(os.environ.get("NSPLITS", 3))
SDM_HOME = "/tmp/claude-0/sdm/SdmPsiGui-linux64-v6.23"
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
ZOOMS = mask_img.header.get_zooms()[:3]




def hedges(mean, sd, n):
    correction = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
    return correction * mean / np.maximum(sd, 1e-9)


def write_sdm_t_map(t_values, n, path):
    """Write a study's t map where SDM will accept it as a supplied image.

    ES-SDM was built for exactly this mixture -- its title is "combines reported peak
    coordinates and statistical parametric maps" -- and `sdm_parse` looks for `<study>.nii.gz`
    beside the peak files. It is strict about the header, and refuses the study outright rather
    than warning: the NIfTI intent must say "t test" and carry the degrees of freedom, and the
    xform codes must be flagged MNI or Talairach (a fresh header defaults to 'aligned').
    """
    template = nib.load(f"{SDM_HOME}/share/sdm_template.nii.gz")
    img = resample_to_img(masker.inverse_transform(t_values), template,
                          interpolation="continuous", force_resample=True, copy_header=True)
    out = nib.Nifti1Image(np.nan_to_num(np.asarray(img.dataobj, dtype=np.float32)),
                          template.affine)
    out.header.set_intent("t test", (float(n) - 1.0,), name="")
    out.header.set_sform(template.affine, code=4)
    out.header.set_qform(template.affine, code=2)
    nib.save(out, path)


def run_sdm(tables, workdir, images=()):
    """Write SDM's inputs, run pp then mi, and return its coefficient in masker space.

    `tables` is [(name, n, height_z, t_map, [(ijk, z), ...]), ...]; `images` names the studies
    to supply as t maps instead of peak files, so SDM gets exactly the parity CBES gets. Its
    threshold column and peak statistics are t on n - 1 df -- it infers a one-sample design
    from the absent n2 column -- so both are mapped back from z, which is what NiMARE carries.
    """
    os.makedirs(workdir, exist_ok=True)
    images = set(images)
    rows = []
    for name, n, height, t_map, found in tables:
        dof = n - 1
        rows.append((name, n, float(z_to_t(np.array([height]), dof)[0])))
        if name in images:
            write_sdm_t_map(t_map, n, f"{workdir}/{name}.nii.gz")
            continue
        if not found:
            open(f"{workdir}/{name}.no_peaks.txt", "w").close()
            continue
        with open(f"{workdir}/{name}.spm_mni.txt", "w") as handle:
            for ijk, value in found:
                xyz = nib.affines.apply_affine(affine, np.asarray(ijk, dtype=float))
                t = float(z_to_t(np.array([abs(value)]), dof)[0]) * np.sign(value or 1.0)
                handle.write(f"{int(round(xyz[0]))},{int(round(xyz[1]))},"
                             f"{int(round(xyz[2]))},{t:.2f}\n")
    with open(f"{workdir}/sdm_table.txt", "w") as handle:
        handle.write("study\tn1\tt_thr\n")
        for name, n, t_thr in rows:
            handle.write(f"{name}\t{n}\t{t_thr:.4f}\n")
    with open(f"{workdir}/sdmpsi_params.xml", "w") as handle:
        handle.write(PARAMS.format(n_studies=len(rows)))

    for command in ("pp", "mi"):
        done = subprocess.run([f"{SDM_HOME}/bin/linux64/sdm_parse", command],
                              cwd=workdir, capture_output=True, text=True, timeout=5400)
        tail = "\n".join(done.stdout.strip().splitlines()[-2:])
        print(f"    sdm {command}: {tail[:120]}", flush=True)
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
            float(np.median(a[quartile] / np.maximum(t[quartile], 1e-6))),
            float(np.median(a[quartile])), float(np.median(t[quartile])))


if __name__ == "__main__":
    S = np.load(CACHE)
    print(f"{CONTRAST}: {S.shape[0]} subjects, {S.shape[1]} voxels at 4 mm")
    print(f"{N_PER_STUDY} subjects x {N_STUDIES} studies, {N_IMAGES} of them sharing images, "
          f"{SCHEME}/{FOCUS} reporting, {N_SPLITS} splits\n", flush=True)

    rng = np.random.default_rng(0)
    rows = {}
    for split in range(N_SPLITS):
        need = N_PER_STUDY * N_STUDIES
        order = rng.permutation(S.shape[0])
        used, held = order[:need], order[need:]
        truth = np.abs(hedges(S[held].mean(0), S[held].std(0, ddof=1), len(held)))

        tables, studies, donor_g, donor_v = [], [], [], []
        work = f"/tmp/claude-0/hcpsdm_{os.getpid()}_{split}"
        for k in range(N_STUDIES):
            block = S[used[k * N_PER_STUDY:(k + 1) * N_PER_STUDY]]
            mean, sd = block.mean(0), block.std(0, ddof=1)
            t = mean / np.maximum(sd / np.sqrt(N_PER_STUDY), 1e-9)
            z = np.nan_to_num(t_to_z(t, dof=N_PER_STUDY - 1))
            found, height = report_peaks(z, mask_bool, shape, ZOOMS, scheme=SCHEME, focus=FOCUS)
            tables.append((f"s{k:02d}", N_PER_STUDY, float(height), t, found))

            meta = {"sample_sizes": [N_PER_STUDY], "reporting_threshold": float(height)}
            analysis = {"id": f"s{k}", "name": "1", "metadata": meta,
                        "points": [], "images": []}
            if k < N_IMAGES:
                g = hedges(mean, sd, N_PER_STUDY)
                var = 1.0 / N_PER_STUDY + g**2 / (2.0 * N_PER_STUDY)
                donor_g.append(g)
                donor_v.append(var)
                # A directory of its own: `sdm_parse pp` scans its working directory for
                # `<study>.nii.gz` and would otherwise see CBES's donor maps beside its own.
                os.makedirs(f"{work}/cbes", exist_ok=True)
                gp = f"{work}/cbes/s{k}_g.nii.gz"
                vp = f"{work}/cbes/s{k}_v.nii.gz"
                nib.save(masker.inverse_transform(g), gp)
                nib.save(masker.inverse_transform(var), vp)
                analysis["images"] = [
                    {"url": gp, "filename": "g", "space": "MNI", "value_type": "g"},
                    {"url": vp, "filename": "v", "space": "MNI", "value_type": "g_var"}]
            elif found:
                analysis["points"] = [
                    {"space": "MNI",
                     "coordinates": [float(c) for c in nib.affines.apply_affine(
                         affine, np.asarray(ijk, dtype=float))],
                     "values": [{"kind": "Z", "value": float(value)}]}
                    for ijk, value in found]
            else:
                continue  # reported nothing at all, and has no image: off the roster
            studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                            "analyses": [analysis]})

        peaks = sum(len(f) for _, _, _, _, f in tables)
        print(f"  split {split + 1}: {peaks} peaks over {N_STUDIES} tables "
              f"({peaks / N_STUDIES:.1f} each), {len(held)} subjects held out", flush=True)

        collection = Studyset({"id": "hcp", "name": "hcp", "studies": studies},
                              target=None, mask=mask_img)
        # The plain-Tobit arm (prevalence pinned at 1) was measured here and was *worse* --
        # 0.60 against 0.63 -- so the option was reverted from the estimator rather than
        # shipped. See notes: the censoring term over-shrinks mu whatever the prevalence does.
        result = CBES(mask=masker, null_method="none",
                      threshold="reporting_threshold").fit(collection)
        cbes_g = np.abs(result.get_map("g", return_type="array").ravel())
        covered = result.get_map("n_studies", return_type="array").ravel() > 0
        cbes_m = (np.abs(result.get_map("g_marginal", return_type="array").ravel())
                  if "g_marginal" in set(result.maps) else None)
        # HCP's prevalence is 1 by construction: every synthetic study draws from the same
        # population. So the fitted prevalence is a direct check on whether the correction
        # invents absence that is not there -- which is the one way it can only do harm here.
        if "prevalence" in set(result.maps):
            pi = result.get_map("prevalence", return_type="array").ravel()
            strong = truth >= np.percentile(truth[covered], 90)
            print(f"    fitted prevalence (true 1.0): median {np.median(pi[covered]):.3f}, "
                  f"at the strongest decile {np.median(pi[covered & strong]):.3f}", flush=True)

        weights = sum(1.0 / np.maximum(v, 1e-9) for v in donor_v)
        images_only = np.abs(sum(g / np.maximum(v, 1e-9) for g, v in zip(donor_g, donor_v))
                             / np.maximum(weights, 1e-12))

        # SDM is ~15 min of pp plus the imputations per split, so the rest of the pipeline is
        # debuggable without it.
        # A coefficient already computed for this split can be reused: an SDM run is minutes
        # of pp plus the imputations, and the split is reproducible from the seed.
        reuse = os.environ.get("REUSE_SDM")
        if reuse and os.path.exists(reuse):
            sdm = np.abs(masker.transform(resample_to_img(
                nib.load(reuse), mask_img, interpolation="continuous",
                force_resample=True, copy_header=True)).ravel())
            print(f"    reusing {reuse}", flush=True)
        else:
            # Parity: SDM gets the same studies as maps that CBES does, as t maps rather than
            # g maps, which is each method's native input for the same information.
            sdm = None if os.environ.get("SKIP_SDM") else run_sdm(
                tables, work, images=[f"s{k:02d}" for k in range(N_IMAGES)])
        use = covered & np.isfinite(truth) & np.isfinite(cbes_g)
        if sdm is not None:
            use = use & np.isfinite(sdm) & (sdm != 0)
        arms = [("images only", images_only), ("CBES g", cbes_g)]
        if cbes_m is not None:
            arms.append(("CBES g_marginal", cbes_m))
        if sdm is not None:
            arms.insert(0, ("SDM-PSI coeff", sdm))
        else:
            print("    SDM produced no coefficient for this split", flush=True)
        for name, est in arms:
            rows.setdefault(name, []).append(score(est, truth, use))
        print(f"    scored over {int(use.sum())} voxels", flush=True)

    print(f"\n{'estimate':>18} {'r':>7} {'rank r':>7} {'AUC':>6} {'mag ratio':>10} "
          f"{'|est| top':>10} {'|ref| top':>10}")
    for name in ("SDM-PSI coeff", "images only", "CBES g", "CBES g_marginal"):
        if name not in rows:
            continue
        r, rho, auc, ratio, top_est, top_ref = np.nanmean(np.array(rows[name]), axis=0)
        print(f"{name:>18} {r:+7.3f} {rho:+7.3f} {auc:6.3f} {ratio:10.2f} "
              f"{top_est:10.3f} {top_ref:10.3f}")
