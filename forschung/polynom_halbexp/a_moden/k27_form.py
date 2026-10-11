"""Strang A: Formfunktion der Eckschicht R_N(v) = E_N(L(1 − v/N^α)) / E_N(L) aus gemeinsam/eckschicht_profil.json;
Konvergenz in N für α = 1 und α = 0.88 (max über das Gitter |v| ≤ 4)."""
import json, sys
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from lin_kern import GEM
d = json.load(open(GEM / "eckschicht_profil.json", encoding="utf-8"))["daten"]
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))
for al in ("1.0", "0.88"):
    Ns = sorted(int(k.split("|")[1]) for k in d if k.startswith(al + "|"))
    R = {}
    for N in Ns:
        t = d[f"{al}|{N}"]
        e0 = cx(t["0|0"]["E"])
        R[N] = {k: cx(v["E"]) / e0 for k, v in t.items() if float(k.split("|")[0]) <= 4}
    ref = R[Ns[-1]]
    print(f"α = {al}: max_v |R_N(v) − R_1000(v)| (|v| ≤ 4):", ", ".join(
        f"N={N}: {mp.nstr(max(abs(R[N][k] - ref[k]) for k in ref if k in R[N]), 3)}" for N in Ns[:-1]))
    print("   R_1000(v) für φ = 0, ±45°, ±90°:")
    for rho in ("0.5", "1", "2", "4"):
        print("    ρ=" + rho + ": " + ", ".join(f"φ={ph}: {mp.nstr(ref[f'{rho}|{ph}'], 4)}" for ph in (-90, -45, 0, 45, 90)))
