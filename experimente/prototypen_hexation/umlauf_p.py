"""Monodromie der Kneser-Abel-Funktion α_K um den anziehenden Fixpunkt p der Pentation (Basis 2).
 α_K über ihex (Pentation/Umkehrung bis in den Rechenbogen, je Punkt unabhängig) auf Kreisen um p; Sprünge = Zweigwechsel.
 Dazu α_K(P(w)) − α_K(w) − 1 als Probe (muss 0 sein, wo beide auf demselben Blatt liegen).
 Fork nur gelesen. Aufruf: umlauf_p.py r1,r2,...
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

radien = [float(r) for r in sys.argv[1].split(",")]
mp.mp.dps = 40
t0 = time.time()
P_FIX = complex(1.7937303432391165676, 1.5736302435118265996)
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)


def ihex_liste(ws, grenze=3):
    ex = [f"iferr(if(1, kv=alarm({grenze}, ihex({w.real:.25f}+I*{w.imag:.25f})); if(type(kv)==\"t_ERROR\", 10^9, kv)),E,10^9)" for w in ws]
    return np.array([complex(v) for v in gp.eval_batch("2", ex)])


try:
    k = gp.eval_batch("2", ["hexinit(sexp(1),0,1)", "hexerr", "knc", "kngood"])
    print(f"Stellen {mp.nstr(k[1], 4)}, Bogen um {mp.nstr(k[2], 6)} mit Radius {mp.nstr(k[3], 5)} [{time.time() - t0:.0f}s]", flush=True)
    for r in radien:
        phi = np.linspace(0, 2 * np.pi, 73)
        ws = P_FIX + r * np.exp(1j * phi)
        A = ihex_liste(ws)
        ok = np.abs(A) < 1e8
        spr = [j for j in range(len(A) - 1) if ok[j] and ok[j + 1] and abs(A[j + 1] - A[j]) > 0.5]
        print(f"r = {r}: {ok.sum()}/{len(A)} Punkte berechnet [{time.time() - t0:.0f}s]", flush=True)
        for j in range(0, len(A), 6):
            print(f"   φ = {phi[j] * 180 / np.pi:5.0f}°: α_K = {A[j]:.6f}" if ok[j] else f"   φ = {phi[j] * 180 / np.pi:5.0f}°: —", flush=True)
        if spr:
            print("   Sprünge:", [(round(phi[j] * 180 / np.pi), f"{A[j]:.4f} → {A[j + 1]:.4f}", f"Δ = {A[j + 1] - A[j]:.5f}") for j in spr], flush=True)
        # Probe Abel-Gleichung an einigen Punkten: α_K(P(w)) − α_K(w) − 1
        idx = [j for j in range(0, len(A), 12) if ok[j]]
        Pw = np.array([complex(v) for v in gp.eval_batch("2", [f"iferr(pent({ws[j].real:.25f}+I*{ws[j].imag:.25f}),E,10^9)" for j in idx])])
        APw = ihex_liste(Pw)
        print("   Probe α_K(P(w)) − α_K(w) − 1:", [f"{(APw[i] - A[j] - 1):.2e}" if abs(APw[i]) < 1e8 else "—" for i, j in enumerate(idx)], flush=True)
finally:
    gp.close()
print(f"fertig [{time.time() - t0:.0f}s]")
