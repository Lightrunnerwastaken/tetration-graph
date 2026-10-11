"""zert_job — EIN Übergang N → N+1 für die parallele Verlängerung des Zweig-Zertifikats (Strang X, Auftrag Leitung).
  1. hochgenaue Mitte von p_N per Newton (nicht streng): Start aus Ps Box (N = 193), aus den Zweig-Daten (polynome_*.json, N ≤ 300)
     oder aus T_N f_K (zertifizierte Kneser-Koeffizienten, N > 300);
  2. Startbox B_N = Mitte ± EPS·ρ^{−k} (EPS = 1e-120) → daten/boxen/B_N.json;
  3. Ps Krawczyk-Übergang (zert_kern.uebergang) mit NC Knoten; scheitert er, NC + 3 (bis zu zweimal) → erg/K_{N+1}.json mit Protokoll.
Aufruf: python zert_job.py N [NC]
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import arb, ctx

import zert_kern as Z

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
DATEN = HIER / "daten"
(HIER / "boxen").mkdir(exist_ok=True)
(HIER / "erg").mkdir(exist_ok=True)
N = int(sys.argv[1])
M = N + 1
NC0 = int(sys.argv[2]) if len(sys.argv) > 2 else (7 if M <= 205 else 7 + math.ceil((M - 205) / 30))
EPS_EXP = 150
t0 = time.time()
_alt = HIER / "erg" / f"K_{M}.json"
if _alt.exists() and (HIER / "boxen" / f"B_{N}.json").exists() and json.loads(_alt.read_text(encoding="utf-8")).get("ok"):
    print(f"ENDE N = {N} → {M}: schon vorhanden (ok)", flush=True)
    sys.exit(0)
ctx.prec = 2 * M + 384


def startwerte():
    if N == 193:
        z = json.loads((DATEN / "p24_zustand.json").read_text(encoding="utf-8"))
        assert z["N"] == 193
        return [Z.dekodiere(x) for x in z["box"]], "P-Box p24"
    if N <= 194:
        return json.loads((DATEN / "polynome_bis_194.json").read_text())[str(N)], "polynome_bis_194"
    if N <= 300:
        return json.loads((DATEN / "polynome_195_bis_300.json").read_text())[str(N)], "polynome_195_bis_300"
    a = json.loads((DATEN / "kneser_halbexp_koeff_zert_640.json").read_text())["a"]
    return [e[0] for e in a[: N + 1]], "T_N f_K (zert. Koeff.)"


start, quelle = startwerte()
mitte, newton = Z.newton_pN(start, N)
konv = newton[-1] < 1e-150     # nur Information: die Kette (zert_kette.py) entscheidet, ob p_N in B_N liegt
# Radius der Startbox: 2^{−(N+150)}·ρ^{−k} — 1000-mal grösser als der Radius 2^{−(N+160)}·ρ^{−k} der Ausgabebox des Vorgängers
# (Ps Verengung: r1 = 100|Y F| + 2^{−prec/2}·w, prec = 2N + 320), damit K_N ⊂ B_N gelten kann; klein genug für die Röhre.
eps = arb(2) ** (-(N + 150))
B = [arb(mitte[k].mid(), eps * arb(Z.RHO) ** (-k)) for k in range(N + 1)]
(HIER / "boxen" / f"B_{N}.json").write_text(json.dumps({"N": N, "quelle": quelle, "eps_bits": N + 150, "newton": newton,
                                                         "konvergiert": bool(konv), "box": [Z.kodiere(t) for t in B]}),
                                             encoding="utf-8")
print(json.dumps({"N": N, "quelle": quelle, "newton": [f"{x:.1e}" for x in newton], "konv": bool(konv), "t_newton": round(time.time() - t0, 1)}),
      flush=True)
versuche = []
ok, K = False, None
for NC in (NC0, NC0 + 3, NC0 + 6):
    ok, K, proto = Z.uebergang(B, N, NC)
    versuche.append(proto)
    print(json.dumps(proto, ensure_ascii=False), flush=True)
    if ok:
        break
aus = {"M": M, "ok": bool(ok), "newton_konvergiert": bool(konv), "versuche": versuche, "t_gesamt": round(time.time() - t0, 1)}
if ok:
    aus["box"] = [Z.kodiere(t) for t in K]
(HIER / "erg" / f"K_{M}.json").write_text(json.dumps(aus), encoding="utf-8")
print(f"ENDE N = {N} → {M}: ok = {aus['ok']}, NC = {versuche[-1].get('NC')}, t = {aus['t_gesamt']} s", flush=True)
