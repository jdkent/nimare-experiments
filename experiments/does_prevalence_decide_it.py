"""Where does the coordinate channel start earning its keep? Vary the prevalence and find out.

Two beds bracketed the regime and disagreed. On the 21-study pain collection -- different
paradigms, different populations, so a prevalence below 1 -- the silence correction cut rmse 23%.
On HCP synthetic studies drawn from one population, prevalence exactly 1, it made the magnitude
worse than doing nothing (0.63 against 0.85), and SDM-PSI given the same data lost too. The
docstring now admits a user cannot tell which regime a real collection is in. This builds the
regime as a dial.

A collection of `n_studies`, of which `round(prevalence * n_studies)` are **MOTOR_LH** studies
that genuinely have the effect, and the rest are **EMOTION_FACES** studies that do not -- real
subjects, real noise structure, real reported peaks, just in the wrong places. So they are silent
at the motor voxels for the right reason rather than by construction, which is what a
heterogeneous literature looks like.

Held-out MOTOR_LH subjects, used to make no coordinate, give the truth. That fixes both estimands
exactly, which is the point:

    mu(v)         the effect among studies that have one   = held-out MOTOR g
    pi(v) mu(v)   the effect averaged over all studies     = prevalence x held-out MOTOR g

`g` claims to estimate the first and `g_marginal` the second, and both claims are checkable here
rather than arguable. Image donors are drawn at random across *both* kinds, so `images only`
estimates the marginal and is the honest comparator for `g_marginal`.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_PER_STUDY = int(os.environ.get("NPER", 30))
N_STUDIES = int(os.environ.get("NSTUDIES", 16))
N_IMAGES = int(os.environ.get("NIMAGES", 2))
N_SPLITS = int(os.environ.get("NSPLITS", 4))
PREVALENCES = [float(x) for x in os.environ.get("PREV", "1.0,0.75,0.5,0.25").split(",")]

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
ZOOMS = mask_img.header.get_zooms()[:3]

EFFECT = np.load("/tmp/claude-0/hcp/MOTOR_LH_masked.npy")
NULL = np.load("/tmp/claude-0/hcp/EMOTION_FACES_masked.npy")


def hedges(block):
    n = block.shape[0]
    correction = 1.0 - 3.0 / (4.0 * (n - 1) - 1.0)
    return correction * block.mean(0) / np.maximum(block.std(0, ddof=1), 1e-9)


def build_study(block, name, as_image, work):
    """One synthetic study: its own one-sample map, reported the way a paper would report it."""
    n = block.shape[0]
    mean, sd = block.mean(0), block.std(0, ddof=1)
    z = np.nan_to_num(t_to_z(mean / np.maximum(sd / np.sqrt(n), 1e-9), dof=n - 1))
    found, height = report_peaks(z, mask_bool, shape, ZOOMS, scheme=SCHEME, focus=FOCUS)
    meta = {"sample_sizes": [n], "reporting_threshold": float(height)}
    analysis = {"id": name, "name": "1", "metadata": meta, "points": [], "images": []}
    if as_image:
        g = hedges(block)
        var = 1.0 / n + g**2 / (2.0 * n)
        os.makedirs(work, exist_ok=True)
        gp, vp = f"{work}/{name}_g.nii.gz", f"{work}/{name}_v.nii.gz"
        nib.save(masker.inverse_transform(g), gp)
        nib.save(masker.inverse_transform(var), vp)
        analysis["images"] = [
            {"url": gp, "filename": "g", "space": "MNI", "value_type": "g"},
            {"url": vp, "filename": "v", "space": "MNI", "value_type": "g_var"}]
        return analysis, meta, (g, var), len(found)
    if not found:
        return None, meta, None, 0
    analysis["points"] = [
        {"space": "MNI",
         "coordinates": [float(c) for c in nib.affines.apply_affine(
             affine, np.asarray(ijk, dtype=float))],
         "values": [{"kind": "Z", "value": float(value)}]} for ijk, value in found]
    return analysis, meta, None, len(found)


def ratio_at_top(est, truth, use):
    """Median ratio over the truth's top quartile, which is the only place a ratio means much."""
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 75)
    return float(np.median(a[top] / np.maximum(t[top], 1e-6)))


if __name__ == "__main__":
    print(f"MOTOR_LH studies carry the effect, EMOTION_FACES studies do not.")
    print(f"{N_PER_STUDY} subjects x {N_STUDIES} studies, {N_IMAGES} random donors, "
          f"{N_SPLITS} splits per prevalence\n")
    print(f"{'true pi':>8} {'reported pi':>11}   " +
          "  ".join(f"{n:>14}" for n in ("g / mu", "g_marg / pi*mu", "images / pi*mu")))

    for prevalence in PREVALENCES:
        n_effect = int(round(prevalence * N_STUDIES))
        rows = {"g": [], "marg": [], "img": [], "pi": []}
        for split in range(N_SPLITS):
            rng = np.random.default_rng(1000 + split)
            need = N_PER_STUDY * N_STUDIES
            eff_order = rng.permutation(EFFECT.shape[0])
            null_order = rng.permutation(NULL.shape[0])
            used_eff = eff_order[: N_PER_STUDY * n_effect]
            held = eff_order[N_PER_STUDY * n_effect:]
            if len(held) < 100:
                print(f"{prevalence:8.2f}  not enough held-out subjects, skipped")
                break
            mu_truth = np.abs(hedges(EFFECT[held]))
            marginal_truth = prevalence * mu_truth

            donors = set(rng.choice(N_STUDIES, size=N_IMAGES, replace=False).tolist())
            work = f"/tmp/claude-0/prev_{os.getpid()}_{split}"
            studies, donor_maps = [], []
            for k in range(N_STUDIES):
                if k < n_effect:
                    block = EFFECT[used_eff[k * N_PER_STUDY:(k + 1) * N_PER_STUDY]]
                else:
                    j = k - n_effect
                    block = NULL[null_order[j * N_PER_STUDY:(j + 1) * N_PER_STUDY]]
                analysis, meta, maps, _ = build_study(
                    block, f"s{k:02d}", k in donors, work)
                if analysis is None:
                    continue
                if maps is not None:
                    donor_maps.append(maps)
                studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                                "analyses": [analysis]})
            if not donor_maps:
                continue

            result = CBES(mask=masker, null_method="none",
                          threshold="reporting_threshold").fit(
                Studyset({"id": "p", "name": "p", "studies": studies},
                         target=None, mask=mask_img))
            covered = result.get_map("n_studies", return_type="array").ravel() > 0
            g = np.abs(result.get_map("g", return_type="array").ravel())
            have = set(result.maps)
            marg = (np.abs(result.get_map("g_marginal", return_type="array").ravel())
                    if "g_marginal" in have else g)
            pi = (result.get_map("prevalence", return_type="array").ravel()
                  if "prevalence" in have else np.ones_like(g))

            w = sum(1.0 / np.maximum(v, 1e-9) for _, v in donor_maps)
            images = np.abs(sum(gm / np.maximum(v, 1e-9) for gm, v in donor_maps)
                            / np.maximum(w, 1e-12))

            use = covered & np.isfinite(mu_truth) & (mu_truth > 0)
            strong = use & (mu_truth >= np.percentile(mu_truth[use], 90))
            # The pain bed reported rmse and bias over *all* voxels; this bed reports a ratio
            # at the truth's top quartile. Those are not the same statistic, and on a map that
            # is mostly near-zero truth read as |g| the first is dominated by the absolute-value
            # floor. Report both so the two beds can be compared at all.
            rows.setdefault("rmse_g", []).append(
                float(np.sqrt(np.mean((g[use] - marginal_truth[use]) ** 2))))
            rows.setdefault("rmse_img", []).append(
                float(np.sqrt(np.mean((images[use] - marginal_truth[use]) ** 2))))
            rows.setdefault("bias_g", []).append(float(np.mean(g[use] - marginal_truth[use])))
            rows.setdefault("bias_img", []).append(
                float(np.mean(images[use] - marginal_truth[use])))
            rows["g"].append(ratio_at_top(g, mu_truth, use))
            rows["marg"].append(ratio_at_top(marg, marginal_truth, use))
            rows["img"].append(ratio_at_top(images, marginal_truth, use))
            rows["pi"].append(float(np.median(pi[strong])))

        if not rows["g"]:
            continue
        print(f"{prevalence:8.2f} {np.mean(rows['pi']):11.3f}   " +
              "  ".join(f"{np.mean(rows[k]):14.2f}" for k in ("g", "marg", "img")) +
              f"   | whole-map vs pi*mu: rmse g {np.mean(rows['rmse_g']):.3f} "
              f"img {np.mean(rows['rmse_img']):.3f}, "
              f"bias g {np.mean(rows['bias_g']):+.3f} img {np.mean(rows['bias_img']):+.3f}")
