"""A4: Tabelle eps0(b) = b ln b dT'_b(0)/db / T'_b(0) aus allen d2_rand/d2_fein-Dateien; Vergleich b->1 mit 1 - 1/ln(1/ln b),
Steigung in ln ln b fuer grosse b.  Aufruf: python a4_eps0_tabelle.py"""
import json, glob, mpmath as mp
mp.mp.dps = 40


def bwert(bs):
    return mp.mpf(eval(bs.replace("^", "**").replace("exp(", "mp.exp(").replace("/", "*mp.mpf(1)/")))


tab = {}
for f in sorted(glob.glob("cache/d2_*.json")):
    for bs, rows in json.load(open(f)).items():
        r0 = [r for r in rows if mp.mpf(r["x"]) == 0][0]
        b = bwert(bs); e0 = b * mp.log(b) * mp.mpf(r0["D1"]) / mp.mpf(r0["T1"])
        tab.setdefault(bs, []).append((f.split("_")[1] + "_" + f.split("_")[2] + "_" + f.split("_")[3], e0, mp.mpf(r0["T1"])))
rows = sorted(tab.items(), key=lambda kv: bwert(kv[0]))
for bs, v in rows:
    b = bwert(bs); lb = mp.log(b)
    sp = max(e for _, e, _ in v) - min(e for _, e, _ in v) if len(v) > 1 else None
    extra = ""
    if b < mp.mpf("1.06"):
        pred = 1 - 1 / mp.log(1 / lb)
        extra = f"  1-1/ln(1/ln b)={mp.nstr(pred,10)}  Abw/ln b={mp.nstr((v[0][1]-pred)/lb,6)}"
    print(f"b={bs:>22}  eps0={mp.nstr(v[0][1],15)}  T'(0)={mp.nstr(v[0][2],15)}  Einstellungen {len(v)}, Spanne {mp.nstr(sp,3) if sp is not None else '-'}{extra}")
big = [(bwert(bs), v[0][1]) for bs, v in rows if bwert(bs) >= 10]
for (b1, e1), (b2, e2) in zip(big, big[1:]):
    print(f"Steigung d eps0 / d lnln b zwischen {mp.nstr(b1,4)} und {mp.nstr(b2,4)}: {mp.nstr((e2-e1)/(mp.log(mp.log(b2))-mp.log(mp.log(b1))),6)}")
