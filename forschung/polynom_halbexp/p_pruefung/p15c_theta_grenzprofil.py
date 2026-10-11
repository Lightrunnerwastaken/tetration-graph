"""Strang P, Prüfung 15c: Θ' im K'-Rahmen für ein Grenzprofil auf Γin (N → ∞-Extrapolation), λ = 0 auf f_K(Γout) (konservativ):
λ(t) = min(A t^p, λ_sp), t = τ/103.38°, mit A, p, λ_sp aus D (t^{1.44…1.52}, Spitzenrate 0.0118) — und Varianten."""
import sys
import numpy as np
import kprime as K
import p15_theta as P
sys.stdout.reconfigure(encoding="utf-8")
rho, L = K.rho, K.L
tg = np.concatenate([np.logspace(-7, 0, 120)])
s_t = np.array([K.s_Gin(rho * np.exp(1j * (np.angle(L) + np.radians(103.38 * t)))) for t in tg])
pos = np.array([[s, 0.0] for s in np.logspace(-9, np.log10(1.07), 30)])
slit = [(u, K.tau(L * (1 + u))) for u in np.logspace(-6, np.log10(0.2043), 25)]
for (A, p, lsp, name) in [(0.263, 1.5, 0.0118, "D-Profil (p=1.5)"), (0.263, 1.44, 0.0118, "p = 1.44"),
                          (0.263, 1.3, 0.0118, "p = 1.30"), (0.263, 1.2, 0.0118, "p = 1.20 (< 1/α')"),
                          (0.526, 1.5, 0.0236, "doppelt"), (0.04, 1.0, 0.0118, "linear 0.04·t (Ds alte Lesart)")]:
    lam = np.minimum(A * tg ** p, lsp)
    neg = np.column_stack([s_t, lam])
    neg = np.vstack([[[-1e-12, 0.0]], neg])
    rs = [P.omega(t, neg, pos) / P.gprime_tau(t) for (u, t) in slit]
    out = []
    for d in [0.05, 0.1]:
        kk = np.exp(d); T0 = P.T0
        c = (kk ** 2 * T0 - np.conj(T0)) / (kk ** 2 - 1); r = kk * abs(T0 - np.conj(T0)) / (kk ** 2 - 1)
        out.append(max(P.omega(c + r * np.exp(1j * ph), neg, pos) / d for ph in np.linspace(0, 2 * np.pi, 121)[:-1]))
    print(f"{name:32s}: Schlitz max ω'/g' = {max(rs):.4f} (u = {slit[int(np.argmax(rs))][0]:.0e}), bei u = 1e-6: {rs[0]:.4f};"
          f"  Γ'_0.05: {out[0]:.4f}, Γ'_0.1: {out[1]:.4f}")
