"""X14 — Konstanten für die Schwanzabschätzung (Satz X3): M_R = max_{|x|=R}|f_K(x)| für R = 1.65 und R = 2, über den ganzen Kreis.

f_K ausserhalb D über f_K(x) = exp(f_K(Log x)) (B Lemma 8.4(ii)), rekursiv bis Log^{(j)} x ∈ D (dort Kneser-Reihe, 4000 Koeffizienten);
Hauptzweig (Schlitze von L, L̄ radial; Werte auf dem Schlitz sind ein Ufer). Ausserdem |f_K(L(1+u))|/ρ entlang des Schlitzes (Konstante c₁).
"""
import json
import sys
from pathlib import Path

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 30
A = [mp.mpf(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"][:2500]]
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(L)


def fK(z, tiefe=0):
    if abs(z) < 0.985 * rho:
        s = mp.mpc(0)
        for c in reversed(A):
            s = s * z + c
        return s
    if tiefe > 6:
        raise RuntimeError("zu tief")
    return mp.exp(fK(mp.log(z), tiefe + 1))


for R in (mp.mpf("1.65"), mp.mpf(2)):
    werte = []
    for k in range(0, 721):
        ph = mp.pi * k / 720
        z = R * mp.expj(ph)
        try:
            werte.append((abs(fK(z)), float(ph)))
        except RuntimeError:
            werte.append((mp.nan, float(ph)))
    M, phm = max(werte, key=lambda t: t[0])
    ausw = [werte[k] for k in (0, 90, 180, 240, 300, 360, 480, 600, 720)]
    print(f"R = {R}: max|f_K| = {mp.nstr(M, 8)} bei φ = {phm:.3f};  Stichproben |f_K| bei φ = "
          + ", ".join(f"{p:.2f}:{mp.nstr(v, 5)}" for v, p in ausw))
    for q in (2, 3, 4):
        print(f"     Basis (M_R/ρ)(ρ/R)^(1+q), q = {q}: {mp.nstr((M / rho) * (rho / R) ** (1 + q), 5)}")
print("Schlitz: |f_K(L(1+u))|/ρ und (|f/ρ| − 1)/u:")
for u in ("0.01", "0.05", "0.1", "0.15", "0.2"):
    u = mp.mpf(u)
    z = L * (1 + u) * mp.expj(mp.mpf("1e-12"))
    v = abs(fK(z)) / rho
    print(f"   u = {u}: {mp.nstr(v, 8)}, c₁(u) = {mp.nstr((v - 1) / u, 5)}")
