"""X21 — Explizite Konstanten für X4a/X5/X6 (KONSTANTEN.md), Kneser, ℓ²_ρ / L²(dφ/2π) auf |x| = ρ.
 (1) Randwerte auf ∂D: max/min |f'|, max/min |m| (m = eˣ/f'), max |f|.
 (2) M2 := sup_{x∈∂D} |f''(f(x))| (für |M_P − m| ≤ ½M2·R); f'' ausserhalb D über f(w) = exp(f(log w)).
 (3) Bernstein-Konstante K_N := max_{δ∈P_N} ‖δ'∘f‖/(N‖δ∘f‖) in L²(∂D, dφ/2π) (Vandermonde mit Arnoldi auf z_j = f(ρe^{iφ_j})).
 (4) Schwanz τ_N aus Lemma X5 (Zeilen > 5N): Kreisteil, Schlitzteil, m-Teil — explizit.
Aufruf: python x21_konstanten.py
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = np.array([float(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]])
K = len(A)
L = complex(mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337)))
rho = abs(L)
k = np.arange(K)
A1 = A[1:] * k[1:]            # f'
A2 = A1[1:] * k[1:-1]         # f''


def reihe(c, z):
    return np.polynomial.polynomial.polyval(z, c)


def f_in(z):
    return reihe(A, z), reihe(A1, z), reihe(A2, z)


def f_all(w, tiefe=0):
    """(f, f', f'') bei w; innen Reihe, aussen f(w) = exp(f(log w))."""
    if abs(w) <= rho * 1.0000001:
        return f_in(w)
    if tiefe > 8:
        raise RuntimeError
    lw = np.log(w)
    g0, g1, g2 = f_all(lw, tiefe + 1)
    G1 = g1 / w
    G2 = (g2 - g1) / w ** 2
    F = np.exp(g0)
    return F, F * G1, F * (G1 ** 2 + G2)


Mpts = 4096
phi = 2 * np.pi * (np.arange(Mpts) + 0.5) / Mpts
x = rho * np.exp(1j * phi)
F0, F1, F2 = f_in(x)
m = np.exp(x) / F1
print(f"(1) auf |x| = ρ: max|f| = {np.abs(F0).max():.4f}, max|f'| = {np.abs(F1).max():.4f}, min|f'| = {np.abs(F1).min():.4f}, "
      f"max|m| = {np.abs(m).max():.4f}, min|m| = {np.abs(m).min():.4f}; L²-Norm ‖m‖ = {np.sqrt(np.mean(np.abs(m) ** 2)):.4f}")
f2f = np.array([f_all(w)[2] for w in F0])
print(f"(2) M2 = sup|f''(f(x))| = {np.abs(f2f).max():.4f} (bei φ = {np.degrees(phi[np.argmax(np.abs(f2f))]):.1f}°); Mittel {np.mean(np.abs(f2f)):.4f}")

# (3) Arnoldi auf z = f(x), Gewicht 1/M
z = F0
for N in (() if len(sys.argv) > 1 else (25, 50, 100, 150, 193)):
    Q = np.zeros((Mpts, N + 1), complex)
    D = np.zeros((Mpts, N + 1), complex)
    Q[:, 0] = 1.0
    for j in range(N):
        v = z * Q[:, j]
        dv = Q[:, j] + z * D[:, j]
        for _ in range(2):
            for i in range(j + 1):
                h = np.vdot(Q[:, i], v) / Mpts
                v = v - h * Q[:, i]
                dv = dv - h * D[:, i]
        hn = np.sqrt(np.mean(np.abs(v) ** 2))
        Q[:, j + 1] = v / hn
        D[:, j + 1] = dv / hn
    s = np.linalg.svd(D / np.sqrt(Mpts), compute_uv=False)[0]
    print(f"(3) N = {N}: max ‖δ'∘f‖/‖δ∘f‖ = {s:.2f} = {s / N:.4f}·N", flush=True)

# (4) Schwanz Lemma X5
R = 1.649
MR = 2.7251556
beta = (MR / rho) * (rho / R) ** 5
sig, CJ, c1 = 5.447, 0.95, 1.2
mk = None
# Taylor-Koeffizienten von m (float, aus eˣ/f'), ℓ²_ρ-skaliert: über FFT auf |x| = r < ρ ist instabil — direkte Reihe:
# m-Koeffizienten ℓ²_ρ-skaliert per FFT auf |x| = ρ (m stetig auf ∂D, Koeffizienten fallen ~k^{-5.4}); f' mit skalierten Koeffizienten
MF = 1 << 15
w = np.exp(2j * np.pi * np.arange(MF) / MF)
_rs = mp.mpf(abs(mp.findroot(lambda t: mp.exp(t) - t, mp.mpc(0.318, 1.337))))
_as = json.loads((HIER.parent / 'gemeinsam' / 'kneser_halbexp_koeff_4000.json').read_text())['a']
b1 = np.array([float(kk * mp.mpf(_as[kk]) * _rs ** kk) for kk in range(1, K)])   # ρ^k-skaliert: f'(ρw) = Σ b1[k-1] w^{k-1}·(1/ρ)
b1 = b1 / rho
fp = np.polynomial.polynomial.polyval(w, b1)
mvals = np.exp(rho * w) / fp
mr = np.abs(np.fft.fft(mvals) / MF)[: MF // 2]
for N in (50, 100, 150, 193, 300):
    kreis = np.sqrt(N + 1) * beta ** N / np.sqrt(1 - (rho / R) ** 2)
    ks = np.arange(5 * N + 1, 200 * N)
    a = (ks + 1) / 1.2 - (N - 1) * c1
    Fk = (1 / np.pi) * N * np.sqrt(N + 1) / rho * CJ * float(mp.gamma(sig + 1)) * a ** (-sig - 1) * 1.01
    schlitz = np.sqrt(np.sum(Fk ** 2))
    mteil = mr[4 * N:].sum()
    print(f"(4) N = {N}: τ_N ≤ Kreis {kreis:.2e} + Schlitz {schlitz:.2e} + m-Teil {mteil:.2e} = {kreis + schlitz + mteil:.2e}", flush=True)
