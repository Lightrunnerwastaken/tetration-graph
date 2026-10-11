"""Sprungvektor der q-Abel-Koordinate W = α_q auf einem Kreis um den anziehenden Fixpunkt p (Fork: knalphap je Punkt, eigene Stetigkeitskorrektur
 nur um Vielfache von Π). Vergleich des verbleibenden Sprungs mit Π_p = 2πi/ln P'(p) und mit der Abel-Gleichung W(P(w)) = W(w) + 1.
 Fork nur gelesen. Aufruf: umlauf_w2.py r
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import numpy as np
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

r = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
mp.mp.dps = 40
t0 = time.time()
P_FIX = complex(1.7937303432391165676, 1.5736302435118265996)
MU_P = complex(-0.02631375492, 0.6081415615)
PI_P = 2j * np.pi / np.log(MU_P)
PI_ = 2j * np.pi / np.log(complex(1.368932747339784991147632, 2.488144667623474107366636))
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=30)
try:
    gp.eval_batch("2", ["hexinit(sexp(1),0,1)"])
    phi = np.linspace(0, 2 * np.pi, 721)
    ws = P_FIX + r * np.exp(1j * phi)
    roh = np.array([complex(v) for v in gp.eval_batch("2", [f"iferr(knalphap({w.real:.20f}+I*{w.imag:.20f}),E,10^9)" for w in ws])])
    Pw = np.array([complex(v) for v in gp.eval_batch("2", [f"pent({w.real:.20f}+I*{w.imag:.20f})" for w in ws[::60]])])
    rohP = np.array([complex(v) for v in gp.eval_batch("2", [f"iferr(knalphap({w.real:.20f}+I*{w.imag:.20f}),E,10^9)" for w in Pw])])
finally:
    gp.close()
ok = np.abs(roh) < 1e8
print(f"r = {r}: {ok.sum()}/{len(ws)} Punkte; Π_p = {PI_P:.6f}, |Π_p| = {abs(PI_P):.4f}; Π = {PI_:.6f} [{time.time() - t0:.0f}s]")
W = roh.copy()
for j in range(1, len(W)):
    if not (ok[j] and ok[j - 1]):
        continue
    best = min(range(-4, 5), key=lambda k: abs(W[j] + k * PI_ - W[j - 1]))
    W[j:] = W[j:] + best * PI_ if False else W[j:]
    W[j] = W[j] + best * PI_
spr = [(j, W[j] - W[j - 1]) for j in range(1, len(W)) if ok[j] and ok[j - 1] and abs(W[j] - W[j - 1]) > 0.3]
for j, d in spr:
    # Zerlegung d = a·Π_p + b·Π (reell), plus Rest
    M = np.array([[PI_P.real, PI_.real], [PI_P.imag, PI_.imag]])
    a, b = np.linalg.solve(M, [d.real, d.imag])
    print(f"  Sprung bei φ = {phi[j] * 180 / np.pi:6.1f}°: ΔW = {d:.6f}  = {a:.4f}·Π_p + {b:.4f}·Π")
print("  Abel-Gleichung W(P(w)) − W(w) − 1 (bis auf kΠ) an 12 Punkten:",
      [f"{min((rohP[i] - roh[60 * i] - 1 + k * PI_ for k in range(-4, 5)), key=abs):.1e}" if abs(rohP[i]) < 1e8 and ok[60 * i] else "—" for i in range(len(Pw))])
