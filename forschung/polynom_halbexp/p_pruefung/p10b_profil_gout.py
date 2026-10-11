"""Strang P, Prüfung 10b: Randprofil der adjungierten Eckmode beidseitig von L (Γin: τ > 0, Γout: τ < 0), kubisch (Daten W)."""
import sys
from pathlib import Path
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
WD = Path(__file__).resolve().parent.parent / "w_grenzschicht"
L = complex(0.31813150520476413, 1.3372357014306895); rho = abs(L); thL = np.angle(L)
for N in [400, 1000]:
    d = np.load(WD / f"w01_N{N}.npz"); v = d["Vt"][0]; s1 = d["S"][0]
    def V(x):
        w = rho / x; s = 0j
        for c in v[::-1]:
            s = s * w + c
        return s
    out = []
    for t in [-120, -60, -30, -16, -8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, 16, 30, 60, 103]:
        out.append(f"{t:+g}°:{(np.log(abs(V(rho*np.exp(1j*(thL+np.deg2rad(t)))))) + np.log(s1))/N:+.4f}")
    print(f"N = {N}: log(σ₁|V₁(ρe^(i(argL+τ)))|)/N  " + "  ".join(out))
