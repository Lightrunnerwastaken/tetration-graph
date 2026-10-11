"""W21 — Eckanteil T(η) = ⟨u₁, y_η⟩ (kopffreie Pol-Antwort) auf einem Polargitter um L (unskaliert), genau (arb).
T(η) = ⟨Z, T_N J∞k_η⟩ − ⟨Du, T_N k_η⟩,  Jᵀ Z = D u₁.  Ausgabe (1/N)ln|T| als Tabelle; Datei w21_N{N}.npz.
Aufruf: python w21_eckanteil_gitter.py N
"""
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
n = N + 1
ctx.prec = int(1.5 * N) + 400
ctx.cap = N + 10            # WICHTIG: sonst schneidet acb_series bei 10 Termen ab
coef, _ = W.kubik_koeff()
rho = W.RHO
L = W.L
d = np.load(W.HIER + f"/w01_N{N}.npz")
u1 = d["U"][:, 0]
t0 = time.time()
J = W.J_mat(W.kubik_a(N, coef), N)
rr = arb(rho)
Dk = [rr ** k for k in range(n)]
Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)]), algorithm="approx")
fpol = acb_poly([acb(float(c)) for c in coef])
fp = fpol.derivative()
mser = acb_poly([0])
Pw = acb_poly([1])
for c in fp.coeffs():
    mser += c * Pw
    Pw = (Pw * fpol).truncate(n)
mser = mser.truncate(n)
uD = [arb(float(u1[k])) * Dk[k] for k in range(n)]
Zl = [Z[j, 0] for j in range(n)]
print(f"N={N}: gelöst ({time.time()-t0:.0f}s)", flush=True)


def T(eta):
    et = acb(eta.real, eta.imag)
    s1 = list(acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs())
    kz = [et ** (-(j + 1)) for j in range(n)]
    s2 = list((mser * acb_poly(kz)).truncate(n).coeffs())
    s1 += [acb(0)] * (n - len(s1))
    s2 += [acb(0)] * (n - len(s2))
    s = acb(0)
    for j in range(n):
        s += Zl[j] * (s1[j] + s2[j]) - uD[j] * kz[j]
    m = abs(s).mid()
    return float(m.log()) / N if m > 0 else -1e9


radii = [0.003, 0.01, 0.03, 0.06, 0.1, 0.2]
# Winkel relativ zur Aussennormalen arg L (Richtung von L weg): ψ ∈ (−90°, 90°) ist ausserhalb von D
psis = np.arange(-180, 181, 15)
argL = np.angle(L)
tab = np.full((len(radii), len(psis)), np.nan)
for i, r in enumerate(radii):
    for j, p in enumerate(psis):
        eta = L + r * np.exp(1j * (argL + np.deg2rad(p)))
        tab[i, j] = T(eta)
    print(f"r={r}: " + " ".join(f"{p}:{v:+.4f}" for p, v in zip(psis, tab[i])), flush=True)
np.savez(W.HIER + f"/w21_N{N}.npz", tab=tab, radii=radii, psis=psis)
print(f"fertig ({time.time()-t0:.0f}s)")
