"""D20 — Toeplitz-Modell der Eckmode: endliche Abschnitte T_N(a) mit stückweise stetigem Symbol
a(φ) = m(x) = f'(f(x)) auf Γin (|f(x)| < |L|), a = Λ·e^{iψ} (konstant, gross) auf Γout, x = |L|e^{iφ}.
Frage: entsteht bei konstantem Λ eine wachsende Eckmode (σ1(T_N(a)^{-1}) → ∞)?  Und mit „wachsendem“ Γout-Symbol
a_j(φ) = Λ·exp(j·κ·dist(φ, Γin)) (nicht-Toeplitz, imitiert die Komposition |f(x)/x|^j)?
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
L = 1.3745570107437075
a_coef = [float(x) for x in json.load(open("../gemeinsam/polynome_bis_194.json"))["194"]]
pp = np.polynomial.Polynomial(a_coef)
dp = pp.deriv()
M = 1 << 14
phi = 2 * np.pi * np.arange(M) / M
x = L * np.exp(1j * phi)
fx = pp(x)
m = dp(fx)                       # f'(f(x)) = e^x / f'(x)
inn = np.abs(fx) < L
phiL = np.angle(complex(0.31813150520476413, 1.3372357014306895))


def toeplitz_inv_svd(sym, N):
    c = np.fft.fft(sym) / M          # c[d] = Fourier-Koeff. bei e^{-i d φ}?  Konvention: a(φ) = Σ a_d e^{i d φ}
    ad = lambda d: c[d % M]            # a_d = (1/2π)∫ a e^{-idφ}
    T = np.array([[ad(i - j) for j in range(N + 1)] for i in range(N + 1)])
    S = np.linalg.svd(np.linalg.inv(T), compute_uv=False)
    return S[:3]


for Lam, psi in ((10.0, 0.0), (100.0, 0.0), (1e4, 0.0), (100.0, np.pi / 2), (100.0, np.pi), (100.0, -np.pi / 2)):
    sym = np.where(inn, m, Lam * np.exp(1j * psi))
    out = [f"N={N}: " + ", ".join(f"{s:.3f}" for s in toeplitz_inv_svd(sym, N)) for N in (50, 100, 200, 400)]
    print(f"Λ={Lam:g}, ψ={psi:.2f}:  " + " | ".join(out), flush=True)

# „wachsendes“ Γout-Symbol: a_{ij} = Fourier-Koeffizient (i−j) von m·1in + Λ e^{j κ dist}·1out  (spaltenabhängig)
dist = np.minimum(np.abs(((phi - phiL + np.pi) % (2 * np.pi)) - np.pi), np.abs(((phi + phiL + np.pi) % (2 * np.pi)) - np.pi))
for kappa in (0.73, 0.3):
    out = []
    for N in (50, 100, 200, 300):
        T = np.zeros((N + 1, N + 1), dtype=complex)
        for j in range(N + 1):
            sym = np.where(inn, m, np.exp(np.minimum(j * kappa * dist, 600.0)))
            c = np.fft.fft(sym) / M
            for i in range(N + 1):
                T[i, j] = c[(i - j) % M]
        S = np.linalg.svd(np.linalg.inv(T), compute_uv=False)
        out.append(f"N={N}: {S[0]:.3f}, {S[1]:.3f}, {S[2]:.3f}")
    print(f"wachsend κ={kappa}: " + " | ".join(out), flush=True)
