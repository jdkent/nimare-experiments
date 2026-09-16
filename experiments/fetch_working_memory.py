"""Assemble a working-memory image corpus from the NeuroStore base-study endpoint.

`fetch_corpus.py` built a general cognitive corpus and stopped at its target, so it resolved
NeuroVault collections for only the first slice of the 684 base studies NeuroStore lists as
carrying group images. Four of the nineteen studies whose title, abstract or keywords name a
working-memory paradigm fell inside that slice; the rest were never resolved, and the maps that
did come down were deleted after conversion. So a working-memory bed needs its own fetch.

Selection is on the study text (title, abstract, NeuroStore keywords) and then on the image name
within the collection, because a collection carries every contrast a paper reported and only some
of them are the working-memory one. Choosing which contrast a study contributes is not a peak cap:
it fixes what the study is about, not how much of it survives a correction.

Maps are put on the same 4 mm MNI brain mask the pain and HCP beds use, and stored as z so
`report_peaks` can threshold them the way a paper would.
"""
import json, os, subprocess, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.transforms import t_to_z

OUT = "/tmp/claude-0/wm"
INDEX = "/tmp/claude-0/cmp/neurostore_images.json"
os.makedirs(OUT, exist_ok=True)

# Terms that name a working-memory paradigm rather than merely mentioning memory.
STUDY_TERMS = ("working memory", "working-memory", "n-back", "nback", "2-back", "3-back",
               "delayed match", "delayed-match", "digit span", "sternberg",
               "change detection", "memory load", "verbal memory load")
# Terms that pick the working-memory contrast out of a collection's other maps.
IMAGE_TERMS = STUDY_TERMS + ("back >", "back minus", "load", "wm ", " wm", "maintenance",
                             "retention", "delay period", "high load", "encoding")
# Names that mark a map as something other than a group activation contrast.
IMAGE_REJECT = ("mask", "roi", "anatom", "gray matter", "grey matter", "vbm", "fa ",
                "structural", "seed", "resting", "rest ", "connectivity")


def normalise_map_type(raw):
    """NeuroVault reports `T map`/`Z map`, not the short codes the old fetcher assumed."""
    text = str(raw or "").strip().lower()
    if text in ("t", "z"):
        return text
    if "t map" in text or "t-map" in text or "tmap" in text:
        return "t"
    if "z map" in text or "z-map" in text or "zmap" in text:
        return "z"
    return None


def get(url, timeout=60):
    try:
        raw = subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                             capture_output=True, text=True).stdout
        return json.loads(raw)
    except Exception:
        return None


def study_text(rec):
    md = rec.get("metadata") or {}
    return " ".join([rec.get("name") or "", rec.get("description") or "",
                     " ".join(md.get("keywords") or [])]).lower()


def resolve_collection(rec):
    """NeuroStore base study -> NeuroVault collection id, via a version's source url."""
    for vid in (rec.get("versions") or []):
        if not isinstance(vid, str):
            continue
        study = get(f"https://neurostore.org/api/studies/{vid}")
        url = ((study or {}).get("metadata") or {}).get("url") or ""
        if "neurovault.org/collections/" in url:
            cid = url.rstrip("/").split("/")[-1]
            if cid.isdigit():
                return cid
    return None


def pick_image(cid):
    """The working-memory contrast a collection contributes, with its sample size."""
    usable = []
    page = get(f"https://neurovault.org/api/collections/{cid}/images/?format=json")
    while page:
        for img in (page.get("results") or []):
            kind = normalise_map_type(img.get("map_type"))
            n = img.get("number_of_subjects")
            name = str(img.get("name", ""))
            if kind is None or not n or int(n) < 8:
                continue
            if str(img.get("analysis_level", "group")).lower() not in ("group", "", "none"):
                continue
            if any(bad in name.lower() for bad in IMAGE_REJECT):
                continue
            hits = sum(term in name.lower() for term in IMAGE_TERMS)
            usable.append((hits, {"collection": cid, "image": img.get("id"),
                                  "url": img.get("file"), "map_type": kind, "n": int(n),
                                  "name": name[:70]}))
        nxt = page.get("next")
        page = get(nxt) if nxt else None
    if not usable:
        return None
    best = max(hit for hit, _ in usable)
    # Among equally on-topic maps, the largest sample is the most informative reference.
    return max((rec for hit, rec in usable if hit == best), key=lambda r: r["n"])


def to_z(path, map_type, n, grid, mask):
    img = nib.load(path)
    if img.ndim > 3:
        img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine)
    img = resample_to_img(img, grid, interpolation="continuous", force_resample=True,
                          copy_header=True)
    data = np.asarray(img.dataobj, dtype=np.float64)
    values = np.where(np.isfinite(data), data, 0.0)[mask]
    if not np.any(values):
        return None
    return values if map_type == "z" else t_to_z(values, n - 1)


if __name__ == "__main__":
    mask_img = load_mni152_brain_mask(resolution=4)
    mask = np.asarray(mask_img.get_fdata() > 0)
    NiftiMasker(mask_img).fit()
    print(f"grid: {int(mask.sum())} in-mask voxels at 4 mm", flush=True)

    base = json.load(open(INDEX))["results"]
    hits = [r for r in base
            if r.get("level") == "group"
            and (r.get("has_t_maps") or r.get("has_z_maps"))
            and any(term in study_text(r) for term in STUDY_TERMS)]
    print(f"{len(hits)} working-memory base studies with group t/z maps", flush=True)

    cached = f"{OUT}/wm_manifest.json"
    manifest = json.load(open(cached)) if os.path.exists(cached) else []
    done = {rec["study"] for rec in manifest}
    for rec in hits:
        if rec["id"] in done:
            continue
        cid = resolve_collection(rec)
        if cid is None:
            print(f"  {rec['id']}: no NeuroVault collection", flush=True)
            continue
        chosen = pick_image(cid)
        if chosen is None:
            print(f"  {rec['id']} col {cid}: no usable group t/z map", flush=True)
            continue
        chosen["study"] = rec["id"]
        chosen["study_name"] = (rec.get("name") or "")[:70]
        manifest.append(chosen)
        print(f"  {rec['id']} col {cid}: n={chosen['n']} {chosen['map_type']} "
              f"| {chosen['name']}", flush=True)
        json.dump(manifest, open(cached, "w"))

    print(f"\n{len(manifest)} collections resolved; converting", flush=True)
    vectors, meta = [], []
    for rec in manifest:
        path = f"{OUT}/img_{rec['image']}.nii.gz"
        if not os.path.exists(path) or os.path.getsize(path) < 2000:
            subprocess.run(["curl", "-sL", "--max-time", "180", "-o", path, rec["url"]],
                           capture_output=True)
        if not os.path.exists(path) or os.path.getsize(path) < 2000:
            print(f"  {rec['study']}: download failed", flush=True)
            continue
        try:
            z = to_z(path, rec["map_type"], rec["n"], mask_img, mask)
        except Exception as exc:
            print(f"  {rec['study']}: {type(exc).__name__} {exc}", flush=True)
            z = None
        if z is None or not np.isfinite(z).all():
            print(f"  {rec['study']}: unusable after conversion", flush=True)
            continue
        vectors.append(z.astype(np.float32))
        meta.append(rec)

    X = np.array(vectors)
    np.save(f"{OUT}/wm_z.npy", X)
    json.dump(meta, open(f"{OUT}/wm_meta.json", "w"))
    print(f"\n{len(X)} working-memory z maps on the 4 mm grid", flush=True)
    if len(X):
        ns = sorted(r["n"] for r in meta)
        print(f"sample sizes: median {ns[len(ns) // 2]}, range {ns[0]}-{ns[-1]}")
        peak = [float(np.abs(v).max()) for v in X]
        print(f"max |z| per map: median {np.median(peak):.2f}, "
              f"range {min(peak):.2f}-{max(peak):.2f}")
