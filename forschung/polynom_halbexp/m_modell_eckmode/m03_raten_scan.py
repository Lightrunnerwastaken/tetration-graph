"""M03 — Eckmoden-Rate im normierten radialen Quadrat-Modell f(z) = z + ε(z − e^{iθ})(z − e^{−iθ}) (Einheitskreis).
σ₁, σ₂ der Inversen (ℓ², letzte B Spalten) für mehrere N; Fit ln σ₁ = A + cN − a ln N.
Aufruf: python m03_raten_scan.py "eps1,eps2,..." "theta_grad1,..." "N1,N2,..." [B]
Ausgabe: Zeilen + JSON-Anhang an m03_raten.jsonl
"""
import json
import sys
import time

import numpy as np

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
EPS = [float(x) for x in sys.argv[1].split(",")]
THS = [float(x) for x in sys.argv[2].split(",")]
NS = [int(x) for x in sys.argv[3].split(",")]
B = int(sys.argv[4]) if len(sys.argv) > 4 else 150


def modell(eps, th_grad):
    th = np.deg2rad(th_grad)
    return [eps, 1 - 2 * eps * np.cos(th), eps]


for th in THS:
    for eps in EPS:
        coef = modell(eps, th)
        mu = 1 + 2j * eps * np.sin(np.deg2rad(th))
        t0 = time.time()
        res = []
        for N in NS:
            S, U, Vt, cols = K.ecksingulaer(coef, N, rho=1.0, B=B, k=3)
            res.append((N, float(S[0]), float(S[1]), float(S[2])))
            print(f"θ={th}° ε={eps} |μ|={abs(mu):.4f} N={N}: σ = {S[0]:.5g}, {S[1]:.5g}, {S[2]:.5g} ({time.time()-t0:.0f}s)", flush=True)
        Ns = np.array([r[0] for r in res], float)
        ls = np.log([r[1] for r in res])
        fit = None
        if len(Ns) >= 3:
            A = np.vstack([np.ones_like(Ns), Ns, -np.log(Ns)]).T
            fit = np.linalg.lstsq(A, ls, rcond=None)[0].tolist()
            lok = np.diff(ls) / np.diff(Ns)
            print(f"   Fit: c = {fit[1]:.5f}, a = {fit[2]:.3f}; lokale Raten {', '.join(f'{x:.5f}' for x in lok)}; "
                  f"c/ln|μ| = {fit[1]/np.log(abs(mu)):.4f}", flush=True)
        with open("m03_raten.jsonl", "a", encoding="utf-8") as fh:
            fh.write(json.dumps({"eps": eps, "theta": th, "mu": [mu.real, mu.imag], "B": B, "res": res, "fit": fit}) + "\n")
