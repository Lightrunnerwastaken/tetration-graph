"""Strang P, Prüfung 15b: Θ'_N im K'-Rahmen für grosse N, KONSERVATIV:
Γin-Daten = Cs gemessene R_N auf dem D-Bogen (c_literatur/c27_cache.json, nur gelesen; t = 0 … 1, τ = 103.38°·t) bis N = 1000,
auf f_K(Γout) λ := 0 gesetzt (wahr ist λ stark negativ, p14/p15 bis N = 400) — das überschätzt Θ'.
Zusätzlich: dieselbe Rechnung im W̄-Rahmen wäre Cs c29 (dort mit dem A-Bogen)."""
import json
import sys
from pathlib import Path

import numpy as np

import kprime as K
from p15_theta import omega, gprime_tau, T0

sys.stdout.reconfigure(encoding="utf-8")
C27 = json.loads((Path(__file__).resolve().parent.parent / "c_literatur" / "c27_cache.json").read_text())
rho, L = K.rho, K.L
ts = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0]
slit = [(u, K.tau(L * (1 + u))) for u in np.logspace(-5, np.log10(0.2043), 25)]
s_in = {t: K.s_Gin(rho * np.exp(1j * (np.angle(L) + np.radians(103.38 * t)))) for t in ts[1:]}
pos = np.array([[s, 0.0] for s in np.logspace(-9, np.log10(1.07), 30)])
for key in ["500_1000", "600_1140", "700_1280", "800_1420", "900_1560", "1000_1700"]:
    N = int(key.split("_")[0])
    R = C27[key]["R"]
    RL = R["D0.0"]
    neg = np.array([[s_in[t], np.log(R[f"D{t}"] / RL) / N] for t in ts[1:]])
    # zusätzliche Stützstelle sehr nahe L: λ → 0
    neg = np.vstack([[[-1e-9, 0.0]], neg])
    rs = [omega(t, neg, pos) / gprime_tau(t) for (u, t) in slit]
    out = []
    for dlev in [0.05, 0.1]:
        kk = np.exp(dlev)
        c = (kk ** 2 * T0 - np.conj(T0)) / (kk ** 2 - 1)
        r = kk * abs(T0 - np.conj(T0)) / (kk ** 2 - 1)
        vals = [omega(c + r * np.exp(1j * ph), neg, pos) / dlev for ph in np.linspace(0, 2 * np.pi, 181)[:-1]]
        out.append(max(vals))
    print(f"N = {N}: λ(t=0.005/0.05/0.2/1) = " + "/".join(f"{np.log(R[f'D{t}'] / RL) / N:+.5f}" for t in [0.005, 0.05, 0.2, 1.0]) +
          f";  Schlitz max ω'/g' = {max(rs):+.4f} (bei u = {slit[int(np.argmax(rs))][0]:.1e});  Γ'_0.05: {out[0]:+.4f}, Γ'_0.1: {out[1]:+.4f}"
          f"  ⇒ Θ'_N ≤ {max(max(rs), *out):+.4f}")
