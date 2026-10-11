"""W45 — Explizite (gemessene) Konstanten für Satz W1'' (Modell-Linse x = 0.7), N = 150…1000, für Xs KONSTANTEN.md.
H*(N) := max_{γ*} ln(Σ|Y_i|²)^{1/2} (γ* = Rand der Modell-Linse ∪ h(A-Bogen)); M_T(N) := max_{∂K_T} ln(Σ|T_i|²)^{1/2} (Modell-Linse ∪ A-Bogen ∪ L).
ε_eff(N, d) = ε_Γ(d) − M_T(N)/((N+2)d) mit ε_Γ ≥ 1 − M_ℬ/d (w41: M_ℬ = 0.01590).
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
R36, R41 = {}, {}
for f in ("w36_pruefung_N150_300.json", "w36_pruefung_N400_1000.json"):
    R36.update(json.load(open(f)))
for f in ("w41_hstern_N150_300.json", "w41_hstern_N400_1000.json"):
    R41.update(json.load(open(f)))
NS = sorted(int(k) for k in R41)
MB = 0.01590
rows = []
for N in NS:
    a, b = R36[str(N)], R41[str(N)]
    hs = max([v for k, v in b.items() if k.startswith("Y_x0.7")] + [v for k, v in a.items() if k.startswith("V3_")])
    mt = max([v for k, v in b.items() if k.startswith("T_x0.7")] + [v for k, v in a.items() if k.startswith("T_A")] + [a["T_L"]])
    rows.append((N, hs, mt))
print("  N   H*(N)=max ln|Y| auf γ*   M_T(N)=max ln|T| auf ∂K_T   ε_eff(d=0.05)   ε_eff(d=0.1)")
for N, hs, mt in rows:
    e5 = (1 - MB / 0.05) - mt / ((N + 2) * 0.05)
    e10 = (1 - MB / 0.1) - mt / ((N + 2) * 0.1)
    print(f"{N:5d}   {hs:8.3f}                   {mt:8.3f}                 {e5:6.3f}          {e10:6.3f}")
Ns = np.array([r[0] for r in rows], float)
for nm, col in (("H*", 1), ("M_T", 2)):
    y = np.array([r[col] for r in rows])
    k, c = np.polyfit(np.log(Ns), y, 1)
    off = np.max(y - (c + k * np.log(Ns)))
    print(f"{nm}: Fit ln-Wert ≈ {c:+.3f} + {k:.3f}·ln N  ⇒  Schranke e^{{{c+off:.3f}}}·N^{{{k:.3f}}} (Hülle über alle N)")
