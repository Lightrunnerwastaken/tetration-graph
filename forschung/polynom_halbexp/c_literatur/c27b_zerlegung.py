"""Strang C, Runde 4, c27b: Zerlegung der Randwerte R_N(t) aus c27 in Eckteil und Spitzenteil.
Modell je Randpunkt t: R_N(t) = a·R_N(L) + b·e^{rN} (a, b, r frei; Fit N = 400…1000, log-Residuen).
Frage: geht die Spitzenrate r(t) zur Ecke hin gegen 0 (WKB, (E) hält für alle N) oder bleibt sie fest (|s|^p e^{cN}, (E) bricht)?
"""
import re, sys
import numpy as np
from scipy.optimize import least_squares
sys.stdout.reconfigure(encoding="utf-8")
txt = open("c27_rand_gross.log", encoding="utf-8").read().split("N=")
daten = {}
for blk in txt[1:]:
    m = re.match(r"(\d+) prec=(\d+)", blk)
    if not m:
        continue
    N, prec = int(m.group(1)), int(m.group(2))
    R = dict((k, float(v)) for k, v in re.findall(r"([DAS][\d.]+):([\d.e+-]+)", blk.split("Profil")[0]))
    daten.setdefault(N, R)
Ns = np.array(sorted(n for n in daten if n >= 400), dtype=float)
RL = np.array([daten[int(n)]["D0.0"] for n in Ns])
print("N:", Ns.astype(int).tolist())
for name in ["D0.005", "D0.01", "D0.02", "D0.05", "D0.1", "D0.2", "D0.35", "D0.5", "D0.7", "D1.0", "A0.02", "S0.3", "S1.0", "S4.0"]:
    y = np.array([daten[int(n)][name] for n in Ns])
    def res(p):
        a, lb, r = p
        return np.log(a * RL + np.exp(lb + r * (Ns - 1000))) - np.log(y)
    best = None
    for r0 in (0.002, 0.006, 0.0118):
        for a0 in (0.5, 0.9):
            try:
                sol = least_squares(res, [a0, np.log(max(y[-1] * 0.05, 1e-6)), r0], bounds=([0, -60, -0.05], [5, 20, 0.05]))
            except Exception:
                continue
            if best is None or sol.cost < best.cost:
                best = sol
    a, lb, r = best.x
    anteil = np.exp(lb) / y[-1]
    print(f"{name}: a = {a:.4f}, Spitzenrate r = {r:.4f}, Spitzenanteil bei N=1000 = {anteil:.3f}, max|log-Residuum| = {np.max(np.abs(best.fun)):.1e}")
