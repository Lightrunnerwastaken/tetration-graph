"""Unabhängige Nachprüfung der Gegenbeispiele mit höherer Stellenzahl und feinerer Quadratur.
Aufruf: python s5_pruefen.py DPS K  (liest daten/gegenbeispiel_*_40.json)"""
import json
import sys
import time

import mpmath as mp

import e4_lib as E

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
E.init(dps)
for fall in ("a", "b", "c", "d", "ab", "alle"):
    f = E.HIER / "daten" / f"gegenbeispiel_{fall}_40.json"
    if not f.exists():
        continue
    d = json.loads(f.read_text(encoding="utf-8"))
    stoer = ("trig", [(kk, mp.mpf(a), mp.mpf(b)) for kk, a, b in d["p"]])
    t = time.time()
    r = E.kenngroessen(stoer, k, moden=8, momente=4)
    r["b2"] = (r["b"] - E.L ** 2 * r["a"]) / E.L ** 4
    mn = min(1 + E.p_wert(stoer, mp.mpf(j) / 1000)[1] for j in range(1000))
    amp = max(abs(E.p_wert(stoer, mp.mpf(j) / 1000)[0]) for j in range(1000))
    print(f"=== Fall {fall} (Bedingungen {d['bed']}), dps {dps}, k {k}, {time.time()-t:.0f}s; max|p| = {mp.nstr(amp, 4)}, min(1+p') = {mp.nstr(mn, 4)}")
    print("    p:", ", ".join(f"k={kk}: sin {mp.nstr(mp.mpf(a), 8)}, cos−1 {mp.nstr(mp.mpf(b), 8)}" for kk, a, b in d["p"]))
    for c in ("a", "b", "b2", "c", "d", "neg1", "neg2", "neg3", "neg4", "V3", "V4"):
        mark = "  <- Bedingung" if c in d["bed"] else ""
        print(f"    {c:5s} |.| = {mp.nstr(abs(r[c]), 4)}{mark}")
    N = mp.fsum(abs(r[f'neg{m}']) ** 2 for m in range(1, 9))
    print(f"    N_8 = {mp.nstr(N, 6)}")
