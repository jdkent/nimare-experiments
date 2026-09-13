"""Collect group-level t/z maps with sample sizes, for the similarity-vs-magnitude test.

NeuroStore says which base studies have images and links each to its NeuroVault collection;
NeuroVault holds the maps and, crucially, the sample size, which NeuroStore does not carry.
Only images that declare a sample size and a t or z type are usable, since Hedges' g needs
both.
"""
import json, os, subprocess, sys, time

TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else 120
OUT = "/tmp/claude-0/corpus"
os.makedirs(OUT, exist_ok=True)

def get(url, timeout=45):
    try:
        raw = subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                             capture_output=True, text=True).stdout
        return json.loads(raw)
    except Exception:
        return None

base = json.load(open("/tmp/claude-0/cmp/neurostore_images.json"))["results"]
candidates = [r for r in base
              if r.get("level") == "group" and (r.get("has_t_maps") or r.get("has_z_maps"))]
print(f"{len(candidates)} candidate base studies", flush=True)

# NeuroStore -> NeuroVault collection id, via a version's metadata.
collections, seen = [], set()
for rec in candidates:
    for vid in (rec.get("versions") or []):
        if not isinstance(vid, str):
            continue
        study = get(f"https://neurostore.org/api/studies/{vid}")
        url = ((study or {}).get("metadata") or {}).get("url") or ""
        if "neurovault.org/collections/" in url:
            cid = url.rstrip("/").split("/")[-1]
            if cid.isdigit() and cid not in seen:
                seen.add(cid)
                collections.append((rec["id"], cid, rec.get("name", "")[:60]))
            break
    if len(collections) >= TARGET * 2:
        break
print(f"{len(collections)} distinct NeuroVault collections resolved", flush=True)
json.dump(collections, open(f"{OUT}/collections.json", "w"))

manifest = []
for study_id, cid, name in collections:
    if len(manifest) >= TARGET:
        break
    page = get(f"https://neurovault.org/api/collections/{cid}/images/?format=json")
    for img in ((page or {}).get("results") or []):
        kind = str(img.get("map_type", "")).lower()
        n = img.get("number_of_subjects")
        if kind not in ("t", "z") or not n or int(n) < 8:
            continue
        if str(img.get("analysis_level", "group")).lower() not in ("group", "", "none"):
            continue
        manifest.append({
            "collection": cid, "image": img.get("id"), "url": img.get("file"),
            "map_type": kind, "n": int(n), "name": str(img.get("name", ""))[:60],
            "study": study_id, "study_name": name,
        })
        break   # one map per collection, so no single study dominates
    if len(manifest) % 20 == 0 and manifest:
        print(f"  {len(manifest)} usable maps so far", flush=True)

json.dump(manifest, open(f"{OUT}/manifest.json", "w"))
print(f"\n{len(manifest)} usable group t/z maps with sample sizes", flush=True)
if manifest:
    ns = [m["n"] for m in manifest]
    print(f"sample sizes: median {sorted(ns)[len(ns)//2]}, range {min(ns)}-{max(ns)}")
    print(f"map types: {sum(m['map_type']=='t' for m in manifest)} t, "
          f"{sum(m['map_type']=='z' for m in manifest)} z")
