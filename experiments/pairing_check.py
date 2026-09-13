VERSIONS = [
{
"id": "68gGCDXkSQns",
"name": "A neurocognitive investigation of the impact of so",
"versions": [
"7X2QiUEUmMsd",
"4gJvH557287z",
"73h3Urvj2vrL",
"fmBQRVitw5Rj"
]
},
{
"id": "4nwC6Yu36KmY",
"name": "Neural correlates of explicit social judgments on ",
"versions": [
"52uV83WXV2x2",
"5VAseLjzrvou",
"6peujgCVbEx9",
"ZosndnjwXiyP"
]
},
{
"id": "35cg7umQRuQ3",
"name": "The neural basis of conceptualizing the same actio",
"versions": [
"5JZjdGDBWxSQ",
"CwipaUyrDWis"
]
},
{
"id": "Ymurgs6Gj9Mp",
"name": "Neural reward related-reactions to monetary gains ",
"versions": [
"6gjQCQ7y6s3r",
"LrAWrQuxrNYw",
"XNjNov5qMsPf"
]
},
{
"id": "xoSUFosddghF",
"name": "Great expectations: neural computations underlying",
"versions": [
"7UFKY77CuPen",
"4hBW53426tCD"
]
},
{
"id": "XeY5DxcDsN8Y",
"name": "Comparison of gray matter volume between migraine ",
"versions": [
"4FK9gW4ZHfUe",
"tF2sgPoxkg2A",
"VBjotpHmYVJn",
"ZqJ8zSfu8Qmz"
]
},
{
"id": "WutohPkH76VS",
"name": "Amygdala functional connectivity in major depressi",
"versions": [
"iD8KEUe9M97k",
"5vwo6xXZCP74"
]
},
{
"id": "WrZpmZtCtFre",
"name": "Anterior insula coordinates hierarchical processin",
"versions": [
"4FasVqEFPupp",
"kTn6aAHPBk6M",
"6hkTXKeqEdkd",
"7N5TjHmj4LyP"
]
},
{
"id": "wiJD9ggf6F3V",
"name": "Functional MRI of emotional memory in adolescent d",
"versions": [
"3X2CqaffSFdN",
"An2CrszcrLnH"
]
},
{
"id": "V2fkJy8ePKyv",
"name": "Congruency and reactivation aid memory integration",
"versions": [
"5JxXu5vYCcfb",
"jTMy4AQYsJPs",
"bKDMD4z7PuN2"
]
},
{
"id": "UdCx5uAz2cop",
"name": "Neural Correlates of Attentional Flexibility durin",
"versions": [
"538k9i3pBNgS",
"4HKbxRVwrV8W",
"XHkqPfmMfpTL"
]
},
{
"id": "Ub78mJc9PobF",
"name": "Deactivation in the posterior mid-cingulate cortex",
"versions": [
"4deTV8phfZY4",
"7k9urUQE82k3",
"wBAjDigzhrtA"
]
}
]
import json, subprocess, collections
totals = collections.Counter()
for rec in VERSIONS:
    for vid in rec['versions']:
        if not isinstance(vid, str):
            continue
        try:
            raw = subprocess.run(
                ["curl", "-s", "--max-time", "45",
                 f"https://neurostore.org/api/studies/{vid}?nested=true"],
                capture_output=True, text=True).stdout
            s = json.loads(raw)
        except Exception:
            continue
        for a in (s.get('analyses') or []):
            imgs = a.get('images') or []
            pts = a.get('points') or []
            stat = [i for i in imgs
                    if str(i.get('value_type')).lower() in ('t', 'z', 't map', 'z map')]
            totals['analyses'] += 1
            if imgs: totals['with any image'] += 1
            if stat: totals['with t/z image'] += 1
            if pts: totals['with points'] += 1
            if stat and pts: totals['BOTH t/z image and points'] += 1
print(f"\nacross {len(VERSIONS)} base studies, every version, all analyses:")
for k in ('analyses','with any image','with t/z image','with points','BOTH t/z image and points'):
    print(f"  {k:32s} {totals[k]:5d}")
