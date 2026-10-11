"""W19 — Leitordnungs-Vorhersage des Stegprofils λ(φ) und der Rate c (Vermutung W2):
Λ_h := abgeschlossene Linse zwischen A-Bogen (∂D(ã₀,R̃₀) ausserhalb D) und h(A-Bogen) (Urbild unter dem Hauptzweig, innerhalb D);
Φ_B := harmonisch auf Ĉ∖Λ_h, = 0 auf ∂Λ_h, Φ_B ~ −log|t| bei 0 und bei ∞ (Dipol-Green-Funktion).
Vorhersage: λ(φ) = Φ_B(e^{iφ}) auf Γin, c = Φ_B(−1); Stegfrequenz x(φ) = −∂_rΦ_B.
Ansatz: Φ_B = −log|t| + h₀, h₀ harmonisch beschränkt auf Ĉ∖Λ_h mit h₀ = log|t| auf ∂Λ_h (Ladungen in Λ_h, Σq = 0).
Aufruf: python w19_dipol.py [mu_abs mu_arg]
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
mu = float(sys.argv[1]) * np.exp(1j * float(sys.argv[2])) if len(sys.argv) > 2 else np.sqrt(W.L)
coef, Lf = W.kubik_koeff(mu=mu)
rho = W.RHO
Ft = np.polynomial.Polynomial([coef[j] * rho ** (j - 1) for j in range(4)])
Lt = Lf / rho
a0 = Ft(0).real
R0 = abs(Lt - a0)
psiL = np.angle(Lt - a0)
K = 1600
psis = np.linspace(psiL, -psiL, K + 1)
A = a0 + R0 * np.exp(1j * psis)                       # von L̃ nach L̃̄ (ausserhalb D)
hs = []
prev = Lt
for s in A:
    r = (Ft - s).roots()
    z = r[np.argmin(abs(r - prev))]
    hs.append(z)
    prev = z
hs = np.array(hs)                                     # von L̃ nach L̃̄ (innerhalb D)
# Rand der Linse: A (L̃→L̃̄) dann h(A) rückwärts (L̃̄→L̃)
cur = np.concatenate([A, hs[::-1][1:-1]])
t_ = np.gradient(np.concatenate([cur[-1:], cur, cur[:1]]))[1:-1]
# Innennormale der Linse: Linse liegt rechts von A (A läuft im Uhrzeigersinn um ã₀?) — bestimme per Testpunkt
mid = (A[K // 2] + hs[K // 2]) / 2                    # Punkt in der Linse (auf der reellen Achse)
nr = 1j * t_ / abs(t_)
# Vorzeichen so, dass Ladungen in die Linse zeigen: teste am Scheitel
i0 = K // 2
if abs(cur[i0] + nr[i0] * 0.01 - mid) > abs(cur[i0] - nr[i0] * 0.01 - mid):
    nr = -nr
ds = abs(t_)
src = cur + nr * np.minimum(5 * ds, 0.3 * abs(cur - Lt) * abs(cur - np.conj(Lt)) + 1e-12)
Mat = np.log(np.abs(cur[:, None] - src[None, :]))
Mat = np.vstack([np.hstack([Mat, np.ones((len(cur), 1))]), np.hstack([np.ones((1, len(src))), np.zeros((1, 1))])])
rhs = np.concatenate([np.log(np.abs(cur)), [0.0]])
sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=None)
q, C = sol[:-1], sol[-1]
res = (Mat @ sol - rhs)[:-1]
print(f"μ = {mu:.4f}: ã₀ = {a0:.4f}, R̃₀ = {R0:.4f}, A-Scheitel {A[K//2].real:.4f}, h(A)-Scheitel {hs[K//2].real:.4f}")
print(f"Kollokationsrest max {np.abs(res).max():.2e} (an den Ecken grösser), Σq = {q.sum():.2e}")


def h0(t):
    return np.sum(q * np.log(np.abs(t - src))) + C


def PhiB(t):
    return -np.log(abs(t)) + h0(t)


argL = np.degrees(np.angle(Lt))
print(f"Vorhersage c = Φ_B(−1) = {PhiB(-1+0j):.5f}")
print(" τ(°)   λ_pred(τ) = Φ_B(e^{i(argL+τ)})   x_pred = −∂_rΦ_B")
for tau in (0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 103.38):
    t = np.exp(1j * np.deg2rad(argL + tau))
    hh = 1e-5
    xr = -(PhiB(t * (1 + hh)) - PhiB(t * (1 - hh))) / (2 * hh)
    print(f" {tau:6.2f}   {PhiB(t):.5f}     {xr:.4f}")
print("Γout-Seite (τ<0):", " ".join(f"{tau}:{PhiB(np.exp(1j*np.deg2rad(argL+tau))):.4f}" for tau in (-1, -4, -10, -30)))
print("\nΦ_B auf Gitter (Vergleich mit w11, N=600 Messwerten):")
phis_g = [0, 30, 60, 76.6, 90, 120, 180]
for r in (0.3, 0.5, 0.7, 0.9):
    print(f" r={r}: " + " ".join(f"{PhiB(r*np.exp(1j*np.deg2rad(p))):.4f}" for p in phis_g))
print(" Robin-Werte: h0(0) =", round(h0(0j+1e-12), 4), " h0(∞) ≈", round(h0(1e6+0j), 4))
