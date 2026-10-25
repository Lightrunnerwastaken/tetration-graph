"""E6: tau = z + theta(z) auf Im z = y0 nahe der Achse; Kurve ueber [-3,4] (tau(z+1)=tau(z)+1) auf Doppelpunkte.
Gleichwertig: G(q) = q exp(2 pi i theta) schlicht auf |q| < exp(-2 pi y0)."""
import json, numpy as np, lib21 as S
erg = {}
for y0 in [0.01, 0.005, 0.003]:
    x = np.concatenate([np.linspace(0, 0.02, 2000, endpoint=False), np.linspace(0.02, 0.98, 3000, endpoint=False), np.linspace(0.98, 1, 2001)])
    g = S.tau(x + 1j * y0)
    kurve = np.concatenate([g[:-1] + k for k in range(-3, 4)])
    dp = S.doppelpunkte(kurve, False, 1)
    erg[str(y0)] = {"doppelpunkte": str(dp), "Re_tau_minus_x": [float((g - x).real.min()), float((g - x).real.max())],
                    "Im_tau": [float(g.imag.min()), float(g.imag.max())], "q_radius": float(np.exp(-2*np.pi*y0))}
    print(y0, erg[str(y0)], flush=True)
json.dump(erg, open("daten/e6_tau_nah.json", "w", encoding="utf-8"), indent=1)
