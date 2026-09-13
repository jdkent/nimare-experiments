"""Build the manifest of usable reference maps.

Usable means: an unthresholded group-level T or Z map that declares a sample size. Thresholded
maps are excluded because the point of the corpus is to serve as ground truth, and a map that
has already been thresholded is exactly what we are trying to correct for. NeuroVault also
records smoothness_fwhm, which is worth keeping: it is the quantity the regional censoring term
needs and that coordinates cannot supply.
"""
import json, subprocess, collections

def get(url, timeout=45):
    try:
        return json.loads(subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                                         capture_output=True, text=True).stdout)
    except Exception:
        return None

cols = json.load(open("/tmp/claude-0/corpus/collections.json"))
manifest, per_collection = [], collections.Counter()
for study_id, cid, name in cols:
    page = get(f"https://neurovault.org/api/collections/{cid}/images/?format=json")
    for img in ((page or {}).get("results") or []):
        kind = str(img.get("map_type", "")).strip().lower()
        if kind not in ("t map", "z map"):
            continue
        if img.get("is_thresholded"):
            continue          # already thresholded: not ground truth
        level = str(img.get("analysis_level") or "group").lower()
        if level not in ("group", "none", ""):
            continue
        n = img.get("number_of_subjects")
        if not n or int(n) < 10:
            continue
        if not img.get("file"):
            continue
        manifest.append({
            "collection": cid, "image": img.get("id"), "url": img.get("file"),
            "map_type": "t" if kind.startswith("t") else "z", "n": int(n),
            "smoothness": img.get("smoothness_fwhm"),
            "paradigm": img.get("cognitive_paradigm_cogatlas"),
            "name": str(img.get("name", ""))[:50], "study": study_id,
        })
        per_collection[cid] += 1
        if per_collection[cid] >= 2:
            break             # cap per collection so no study dominates the corpus
    if len(manifest) >= 260:
        break

json.dump(manifest, open("/tmp/claude-0/corpus/manifest.json", "w"))
ns = [m["n"] for m in manifest]
sm = [m["smoothness"] for m in manifest if m["smoothness"]]
print(f"{len(manifest)} usable maps from {len(per_collection)} collections")
print(f"  t maps {sum(m['map_type']=='t' for m in manifest)}, "
      f"z maps {sum(m['map_type']=='z' for m in manifest)}")
print(f"  sample size: median {sorted(ns)[len(ns)//2]}, range {min(ns)}-{max(ns)}")
print(f"  smoothness declared for {len(sm)}/{len(manifest)}"
      + (f", median {sorted(sm)[len(sm)//2]:.1f} mm" if sm else ""))
para = collections.Counter(m["paradigm"] for m in manifest if m["paradigm"])
print(f"  distinct cognitive paradigms: {len(para)}; commonest {para.most_common(4)}")
