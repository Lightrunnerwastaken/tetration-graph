"""Strang P, Prüfung 16: Zutaten von „(a′)' ⇒ (E_R')“ im K'-Rahmen.
(1) Univalenz/Bi-Lipschitz von f_K auf ∂D: min_{x≠y∈∂D} |f(x) − f(y)|/|x − y| (Gitter) — mit Argumentprinzip ⇒ f_K univalent auf D̄.
(2) Eingaben w_η = [J_∞ k_η]_{≤N} für η ∈ ∂K' = Γin ∪ f_K(Γout): |w_η|_{ℓ²_ρ} wächst höchstens polynomial (gegen A-Bogen: 1.21^N).
(3) (E_L2)-Test (kubisch, Ws Daten): σ₁|V₁| auf |x| = ρ(1+δ) fern von L, L̄ — exponentiell klein?
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = [mp.mpf(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
mp.mp.dps = 30
L = complex(mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337)))
rho = abs(L)
argL = np.angle(L)


def fK(z):
    s = mp.mpc(0)
    for c in reversed(A):
        s = s * z + c
    return complex(s)


M = 720
th = 2 * np.pi * np.arange(M) / M
xs = rho * np.exp(1j * th)
fx = np.array([fK(mp.mpc(x)) for x in xs])
D = np.abs(fx[:, None] - fx[None, :]) / np.where(np.eye(M) > 0, 1, np.abs(xs[:, None] - xs[None, :]))
np.fill_diagonal(D, np.inf)
i, j = np.unravel_index(np.argmin(D), D.shape)
print(f"(1) min |f(x)−f(y)|/|x−y| auf ∂D (720 Punkte) = {D[i, j]:.4f} bei arg x = {np.degrees(th[i]):.1f}°, arg y = {np.degrees(th[j]):.1f}°")
# Windungszahl um einen Innenpunkt (f(0)) = 1?
wn = np.sum(np.angle(np.roll(fx, -1) - A[0].__float__()) - np.angle(fx - float(A[0])))
dw = np.angle((np.roll(fx, -1) - float(A[0])) / (fx - float(A[0])))
print(f"    Windungszahl von f(∂D) um f(0): {np.sum(dw) / (2 * np.pi):.6f}")

# (2) Eingaben auf ∂K'
a_np = [complex(c) for c in A[:900]]


def inv_ser(eta, N):
    b = [1 / (eta - a_np[0])]
    for k in range(1, N + 1):
        s = sum(a_np[jj] * b[k - jj] for jj in range(1, k + 1))
        b.append(s / (eta - a_np[0]))
    return np.array(b)


print("(2) |[1/(η − f)]_{≤N}|_{ℓ²_ρ} auf ∂K' (N = 100/200/400):")
for name, eta in [("Γin τ=0.5°", rho * np.exp(1j * (argL + np.radians(0.5)))), ("Γin τ=30°", rho * np.exp(1j * (argL + np.radians(30)))),
                  ("Γin −ρ", -rho), ("f(Γout) 0.5° vor L", fK(mp.mpf(rho) * mp.expj(argL - np.radians(0.5)))),
                  ("f(Γout) 30° vor L", fK(mp.mpf(rho) * mp.expj(argL - np.radians(30)))), ("f(Γout) bei ρ", fK(mp.mpf(rho))),
                  ("A-Bogen Scheitel (zum Vergleich)", float(A[0]) + abs(L - float(A[0])))]:
    b = inv_ser(eta, 400)
    w = np.array([abs(x) * rho ** k for k, x in enumerate(b)])
    print(f"   {name:34s}: " + " / ".join(f"{np.sqrt(np.sum(w[:N + 1] ** 2)):.3e}" for N in [100, 200, 400]))

# (3) (E_L2) kubisch
WD = HIER.parent / "w_grenzschicht"
for N in [400, 1000]:
    d = np.load(WD / f"w01_N{N}.npz")
    v, s1 = d["Vt"][0], d["S"][0]

    def V(x):
        w = rho / x
        s = 0j
        for c in v[::-1]:
            s = s * w + c
        return s
    for delta in [0.0, 0.02, 0.05, 0.1]:
        best = -np.inf
        for t in np.linspace(0, 2 * np.pi, 1441)[:-1]:
            if min(abs(np.angle(np.exp(1j * (t - argL)))), abs(np.angle(np.exp(1j * (t + argL))))) < np.radians(15):
                continue
            best = max(best, (np.log(abs(V(rho * (1 + delta) * np.exp(1j * t)))) + np.log(s1)) / N)
        print(f"(3) kubisch N = {N}, δ = {delta}: max log(σ₁|V₁|)/N auf |x| = ρ(1+δ), |arg x ∓ arg L| ≥ 15°: {best:+.5f}")
