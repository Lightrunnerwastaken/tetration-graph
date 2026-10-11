"""W47 — Explizite D6-Buchhaltung für Satz W1'' (kubisch, Modell-Linse x = 0.7): max_{∂K_T}|T_i| ≤ C(f)·sup_{Ext γ*}|Y_i| + max(1/|ζ|)·√(N+1),
C(f) = max_{ζ∈∂K_T}[|m(ζ)/ζ| + Σ_j 1/|x_j f'(x_j)| + Σ_k (2ρ)^k |[Q_ζ]_k|], Q_ζ = (m(x) − m(ζ))/(ζ − x), m = f'∘f (Cauchy für y_k auf |x| = 2ρ).
"""
import sys

import numpy as np

import w_kern as W
import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
src = open(W.HIER + "/w41_hstern_modell.py", encoding="utf-8").read()
exec(src.split("XS = [0.7, 0.9]")[0].split("import w_kern as W")[1].replace("import w24_rand_genau as M", ""))
coef, _ = W.kubik_koeff()
F = np.polynomial.Polynomial(coef)
dF = F.deriv()
m = dF(F)
pts = []
for t in np.linspace(0.05, 180 - np.degrees(argL), 300):
    ph = argL + np.deg2rad(t)
    pts.append(rho * np.exp(s_lens(ph, 0.7) + 1j * ph))
for dd in np.linspace(0.002, 1.9, 300):
    pts.append(M.a_bogen(dd))
R2 = 2 * rho
Cmax, parts, head = 0.0, None, 0.0
for z in pts:
    Q = (m - m(z)) // np.polynomial.Polynomial([z, -1])
    r = (F - z).roots()
    t1 = abs(m(z) / z)
    t2 = sum(1 / abs(x * dF(x)) for x in r)
    t3 = sum(R2 ** k * abs(c) for k, c in enumerate(Q.coef))
    if t1 + t2 + t3 > Cmax:
        Cmax, parts = t1 + t2 + t3, (t1, t2, t3, z)
    head = max(head, 1 / abs(z))
print(f"C(f) = {Cmax:.3f} (|m/ζ| {parts[0]:.3f}, Urbilder {parts[1]:.3f}, niedrige Koeff. {parts[2]:.3f} bei ζ = {parts[3]:.3f}); Kopf max 1/|ζ| = {head:.3f}")
for N in (150, 193, 300, 600, 1000):
    CH = np.exp(-4.06) * N ** 0.883
    b = Cmax * CH + head * np.sqrt(N + 1)
    print(f"N = {N}: (H*)-Weg ⇒ max_{{∂K_T}}|T| ≤ {b:.3f} (ln {np.log(b):.2f})")
