"""W18 — Leitordnung der adjungierten Eckmode: Green-Funktion g_{Ω₀}(t, 0) des Gebiets Ω₀, berandet von Γin (Einheitskreisbogen,
umskaliert) und h(A-Bogen) (Urbild des äusseren Bogens von C₀ unter dem Hauptzweig h = F̃⁻¹). Vorhersage: Stegfrequenz x(φ) = −∂_r g_{Ω₀}
auf Γin; Vergleich mit Messung (w09) und mit ∂_n g_{Ĉ∖W̄}. Ladungssimulation (Kollokation), kubisches Modell.
Aufruf: python w18_green_omega0.py [mu_abs mu_arg]
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
if len(sys.argv) > 2:
    mu = float(sys.argv[1]) * np.exp(1j * float(sys.argv[2]))
else:
    mu = np.sqrt(W.L)
coef, Lf = W.kubik_koeff(mu=mu)
rho = W.RHO
Ft = np.polynomial.Polynomial([coef[j] * rho ** (j - 1) for j in range(4)])
Lt = Lf / rho
a0 = Ft(0).real
R0 = abs(Lt - a0)
psiL = np.angle(Lt - a0)
# A-Bogen: ψ von psiL (bei L̃) über 0 nach −psiL (bei L̃̄), Punkte ausserhalb des Einheitskreises
psis = np.linspace(psiL, -psiL, 801)
A = a0 + R0 * np.exp(1j * psis)
# h längs des Bogens stetig verfolgen, Start bei L̃
hs = []
prev = Lt
for s in A:
    r = (Ft - s).roots()
    z = r[np.argmin(abs(r - prev))]
    hs.append(z)
    prev = z
hs = np.array(hs)
print(f"μ = {mu:.4f}; ã₀ = {a0:.4f}, R̃₀ = {R0:.4f}; h(A-Scheitel) = {hs[400]:.4f}; |h| min = {abs(hs).min():.4f}")
# Rand von Ω₀: Γin (φ von arg L̃ bis 2π − arg L̃) und h(A-Bogen) (von L̃̄ zurück nach L̃)
argL = np.angle(Lt)
phis = np.linspace(argL, 2 * np.pi - argL, 1201)
G1 = np.exp(1j * phis)
G2 = hs[::-1]          # von L̃̄ nach L̃
bd = np.concatenate([G1, G2[1:-1]])
# Ladungen: nach aussen verschoben (Γin: radial nach aussen; h(A): Richtung von 0 weg?) — Normale numerisch
def normale(c):
    t = np.gradient(c)
    nrm = -1j * t / abs(t)       # rechts der Laufrichtung
    return nrm
# Orientierung: Γin gegen Uhrzeigersinn von L nach L̄ (Gebiet links) ⇒ äussere Normale = rechts
cur = np.concatenate([G1, G2[1:-1]])
nr = normale(np.concatenate([cur, cur[:1]]))[:-1]
# Ladungsabstand proportional zum lokalen Punktabstand
ds = abs(np.gradient(cur))
src = cur + nr * 6 * ds
# Kollokation: g(t) = −log|t| + Σ q_j log|t − s_j| + C = 0 auf dem Rand
Mat = np.log(np.abs(cur[:, None] - src[None, :]))
Mat = np.hstack([Mat, np.ones((len(cur), 1))])
rhs = np.log(np.abs(cur))
sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=None)
q, C = sol[:-1], sol[-1]
res = Mat @ sol - rhs
print(f"Kollokationsrest max {np.abs(res).max():.2e}, Σq = {q.sum():.4f} (soll 1), Robin C = {C:.4f}")


def g(t):
    return -np.log(abs(t)) + np.sum(q * np.log(np.abs(t - src))) + C


def x_rad(phi_deg, h=1e-4):
    t = np.exp(1j * np.deg2rad(phi_deg))
    return (g(t * (1 - h)) - g(t * (1 - 2 * h))) / h * 0 + (g(t * (1 - h)) - 0.0) / h   # −∂_r g bei r = 1 (g(1) = 0)


print("φ    x_pred(Ω₀)")
for p in (80, 85, 90, 100, 120, 150, 180):
    print(f"{p:4d}  {x_rad(p):.4f}")
print("Kontrolle innen: g(0.5·e^{iπ}) =", round(g(-0.5 + 0j), 4), " g(0.5) =", round(g(0.5 + 0j), 4), " g(0.3) =", round(g(0.3 + 0j), 4))
