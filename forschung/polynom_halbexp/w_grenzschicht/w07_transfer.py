"""W07 — Adjungierte Eckmode genau (arb) und ihr Transferbild 𝒱 = 𝓛V (Laurent-Polynom um ã₀), kubisches Modell.

Umskaliert (t = x/ρ): F̃(t) = f(ρt)/ρ, ã₀ = F̃(0). y := J̃^{-T} u₁ (= σ₁ v₁ + Rest), V(t) = Σ_k y_k t^{−k},
c_m := ℓ_y((F̃ − ã₀)^m) = Σ_k y_k [t^k](F̃ − ã₀)^m,  𝒱(s) = Σ_m c_m (s − ã₀)^{−m−1}  (= Σ_{F̃(z)=s} V(z)/(z F̃'(z))).
Ausgabe: Profile |V| auf |t| = 1, |𝒱| auf dem Kreis C₀ (Mitte ã₀, Radius R̃₀ = |L̃ − ã₀|), Index-/Phasenraumprofile von c,
Kontrolle 𝒱 über Urbilder. Speichert w07_N{N}.npz (y, c als float + log-Beträge).
Aufruf: python w07_transfer.py N
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, acb, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
n = N + 1
ctx.prec = int(1.2 * N) + 400
coef, _ = W.kubik_koeff()
rho = W.RHO
d = np.load(W.HIER + f"/w01_N{N}.npz")
u1 = d["U"][:, 0]
S1 = d["S"][0]
t0 = time.time()
a = W.kubik_a(N, coef)
J = W.J_mat(a, N)
# umskaliert: J̃ = D J D^{-1}, D = diag ρ^k; J̃^T y = u  ⇔  J^T (D^{-1} y) = D^{-1}... sorgfältig:
# J̃ = D J D⁻¹ ⇒ J̃ᵀ = D⁻¹ Jᵀ D. J̃ᵀ y = u ⇔ Jᵀ (D y) = D u.
rr = arb(rho)
Dk = [rr ** k for k in range(n)]
rhs = arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)])
Z = J.transpose().solve(rhs, algorithm="approx")      # Z = D y
y = [Z[k, 0] / Dk[k] for k in range(n)]
print(f"N={N}: Lösen {time.time()-t0:.0f}s; |y|_2 ≈ {np.sqrt(sum(float(v.mid())**2 for v in y)):.4f} (σ1 = {S1:.4f})", flush=True)
# F̃ in t: Koeffizienten c_j ρ^{j-1}
Ft = [arb(float(coef[j])) * rr ** (j - 1) for j in range(4)]
a0t = Ft[0]
G = arb_poly([arb(0)] + Ft[1:])          # F̃ − ã₀
c = []
P = arb_poly([1])
for m in range(n):
    s = arb(0)
    co = P.coeffs()
    for k in range(min(len(co), n)):
        s += y[k] * co[k]
    c.append(s)
    P = (P * G).truncate(n)
print(f"c_m berechnet ({time.time()-t0:.0f}s)", flush=True)
yf = np.array([float(v.mid()) for v in y])
cf = np.array([float(v.mid()) for v in c])
lc = np.array([float(abs(v).mid().log()) if abs(v).mid() > 0 else -1e9 for v in c])
ly = np.array([float(abs(v).mid().log()) if abs(v).mid() > 0 else -1e9 for v in y])
a0 = float(a0t.mid())
Lt = W.L / rho
R0t = abs(Lt - a0)
np.savez(W.HIER + f"/w07_N{N}.npz", y=yf, c=cf, lc=lc, ly=ly, a0=a0, R0=R0t)
print(f"ã₀ = {a0:.5f}, R̃₀ = {R0t:.5f}")
for x in (0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.99, 1.0):
    k = int(x * N)
    print(f"  k/N={x:.2f}: (1/N)ln|y_k| = {ly[k]/N:+.4f}   (1/N)ln|c_k| = {lc[k]/N:+.4f}   (1/N)ln(|c_k|/R̃₀^k) = {(lc[k] - k*np.log(R0t))/N:+.4f}")


def Vf(t):
    s = 0
    for v in yf[::-1]:
        s = s / 1 * 1
    w = 1 / t
    s = 0j
    for v in yf[::-1]:
        s = s * w + v
    return s


def Vcal(sv):
    w = 1 / (sv - a0)
    s = 0j
    for v in cf[::-1]:
        s = s * w + v
    return s * w


argL = np.angle(Lt)
print("\n|V| auf |t|=1 (log10):", " ".join(f"{p:.0f}:{np.log10(abs(Vf(np.exp(1j*np.deg2rad(p))))):.2f}" for p in (0, 40, 70, 76.6, 80, 90, 120, 150, 180)))
# Kreis C0: Winkel ψ um ã₀; L̃ bei ψ_L = arg(L̃ − ã₀)
psiL = np.degrees(np.angle(Lt - a0))
print(f"ψ_L = arg(L̃ − ã₀) = {psiL:.2f}°")
psis = [0, 20, 40, 60, psiL - 10, psiL - 2, psiL, psiL + 2, psiL + 10, 120, 150, 180]
print("|𝒱| auf C₀ (log10):", " ".join(f"{p:.1f}:{np.log10(abs(Vcal(a0 + R0t*np.exp(1j*np.deg2rad(p))))):.2f}" for p in psis))
# Kontrolle über Urbilder
Fp = np.polynomial.Polynomial([float(v.mid()) for v in Ft])
for sv in (a0 + R0t * np.exp(1j * np.deg2rad(psiL + 10)), 1.2 * np.exp(1j * 2.0), -1.1):
    roots = (Fp - sv).roots()
    dirk = sum(Vf(z) / (z * Fp.deriv()(z)) for z in roots)
    print(f"  Kontrolle s={sv:.3f}: 𝒱 = {Vcal(sv):.4e}, Σ V(z)/(zF̃'(z)) = {dirk:.4e}")
