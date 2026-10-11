"""W12 — Rechte Eckmode u₁ (Polynom in t = x/ρ): (1/N)ln|U(t)| gegen g_{Ĉ∖W̄}(t) (unskaliert ρt), aussen und innen.
u₁ hochgenau: u = J̃^{-1} y / ... hier einfach float-SVD-Vektor (für die n-te-Wurzel-Struktur ausreichend ausserhalb).
Aufruf: python w12_u_green.py N
"""
import sys
import numpy as np
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(W.HIER + f"/w01_N{N}.npz")
u = d["U"][:, 0]; coef = d["coef"]
rho = W.RHO; L = W.L; a0 = coef[0]; R0 = abs(L - a0)
exec(open(W.HIER + "/w11_green_vergleich.py", encoding="utf-8").read().split("def green_W")[1].join(["def green_W", ""]) if False else "")
def green_W(zeta):
    m = (zeta - L) / (zeta - np.conj(L))
    th1 = np.angle((-rho - L) / (-rho - np.conj(L))); th2 = np.angle((a0 + R0 - L) / (a0 + R0 - np.conj(L)))
    def inside(lo, hi, t): return ((t - lo) % (2*np.pi)) < ((hi - lo) % (2*np.pi))
    if inside(th2, th1, 0.0): lo, wid = th2, (th1 - th2) % (2*np.pi)
    else: lo, wid = th1, (th2 - th1) % (2*np.pi)
    ang = (np.angle(m) - lo) % (2*np.pi)
    if ang > wid: return -1.0
    tau = (abs(m) ** (np.pi/wid)) * np.exp(1j*ang*np.pi/wid); tau0 = np.exp(1j*((0.0 - lo) % (2*np.pi))*np.pi/wid)
    return float(np.log(abs((tau - np.conj(tau0)) / (tau - tau0))))
phis = [0, 30, 60, 70, 76.6, 80, 90, 120, 150, 180]
for r in [1.0, 1.05, 1.1, 1.2, 1.4, 1.7, 2.0]:
    vals = [np.log(abs(W.U_auswerten(u, rho*r*np.exp(1j*np.deg2rad(p)))))/N for p in phis]
    gs = [green_W(rho*r*np.exp(1j*np.deg2rad(p))) for p in phis]
    print(f"r={r:4.2f} u: " + " ".join(f"{x:7.4f}" for x in vals))
    print(f"        g: " + " ".join(f"{x:7.4f}" for x in gs))
