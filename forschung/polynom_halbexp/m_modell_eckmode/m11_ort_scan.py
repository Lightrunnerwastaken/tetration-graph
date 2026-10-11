"""M11 — Eckanregung als Funktion des Singularitätsorts (Prüfstein für (E_L)'): normiertes radiales Quadrat-Modell, Exp-Lage.
Eingabe u_τ = (1 − z/P)^σ + (1 − z/P̄)^σ, P = e^{i(θ+τ)} (τ = 0: Ecke; τ > 0: D-Bogen). Π_N u = J_N⁻¹T_N[u∘f + m u].
Gemessen: Fehler sup_{|z|=1}|Π_N u − u|, Eckamplituden a_i = ⟨u_i, Π_N u⟩ (u_i = linke Eck-Singulärvektoren aus m07),
Vorhersage aus Ds Eckformel D6 (nur Pol-Term, Schlitz radial nach aussen):
  a_i ≈ σ_i · (sin πσ/π) · Γ(σ+1) · [m(P) V_i(P)/P · P + konj.-Term] · N_eff^{−σ−1},  N_eff = Schwerpunkt der v_i-Koeffizienten.
Aufruf: python m11_ort_scan.py sigma "N1,N2,..." "tau1,tau2,..."
"""
import math
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, ctx

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
sig = float(sys.argv[1])
NS = [int(x) for x in sys.argv[2].split(",")]
RADIAL = sys.argv[3].startswith("r")
TAUS = [float(x) for x in sys.argv[3].lstrip("r").split(",")]
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)
coef = [eps, c0, eps]


def pot_reihe(A, s, n):
    g = [A[0] ** s]
    for nn in range(1, n):
        acc = acb(0)
        for k in range(1, min(len(A) - 1, nn) + 1):
            acc += (s * k - (nn - k)) * A[k] * g[nn - k]
        g.append(acc / (nn * A[0]))
    return g


def eingabe(P, N):
    Pa = acb(P.real, P.imag)
    s = arb(sig)
    u = pot_reihe([acb(1), -1 / Pa], s, N + 1)
    A = [1 - acb(coef[0]) / Pa, -acb(coef[1]) / Pa, -acb(coef[2]) / Pa]
    uf = pot_reihe(A, s, N + 1)
    m = [acb(c0 + 2 * eps * coef[0]), acb(2 * eps * coef[1]), acb(2 * eps * coef[2])]
    rhs = [uf[k] + sum((m[i] * u[k - i] for i in range(3) if k - i >= 0), acb(0)) for k in range(N + 1)]
    return [2 * x.real for x in rhs]


phis = np.linspace(0, 2 * np.pi, 4001)
zs = np.exp(1j * phis)
pref = abs(math.sin(math.pi * sig)) / math.pi * math.gamma(sig + 1)
for N in NS:
    t0 = time.time()
    d = np.load(f"m07_N{N}_e{eps}_t{thg}.npz")
    S, U, Vt = d["S"], d["U"], d["Vt"]
    k = np.arange(N + 1)
    ctx.prec = int(0.9 * N) + 300
    J = K.jmat(coef, N, reell=True)
    PUNKTE = [(np.exp(1j * th) * (1 + t)) if RADIAL else np.exp(1j * (th + np.deg2rad(t))) for t in TAUS]
    rhs_cols = [eingabe(P_, N) for P_ in PUNKTE]
    B = arb_mat(N + 1, len(TAUS))
    for c, col in enumerate(rhs_cols):
        for i in range(N + 1):
            B[i, c] = col[i]
    X = J.solve(B, algorithm="approx")
    print(f"N={N} ({time.time()-t0:.0f}s): σ₁ = {S[0]:.4g}, σ₂ = {S[1]:.4g}", flush=True)
    for c, t in enumerate(TAUS):
        P = PUNKTE[c]
        p = np.array([float(X[i, c].mid()) for i in range(N + 1)])
        uex = (1 - zs / P) ** sig + (1 - zs / np.conj(P)) ** sig
        err = np.abs(np.polyval(p[::-1], zs) - uex)
        a = [float(U[:, i] @ p) for i in range(2)]
        vorh = []
        from scipy.special import gammaln
        kk = k[k > sig + 1.5]
        lb = gammaln(sig + 1) + gammaln(kk - sig - 1) - gammaln(kk)       # log B(σ+1, k−σ−1)
        mP = c0 + 2 * eps * K.polywert(coef, P)
        for i in range(2):
            vv = Vt[i][k > sig + 1.5]
            I1 = np.sum(vv * P ** (-kk) * np.exp(lb))                       # ∫ t^σ V_i(P(1+t)) dt (Pol-Term)
            I2 = np.sum(vv * np.conj(P) ** (-kk) * np.exp(lb))
            vorh.append(S[i] * abs(math.sin(math.pi * sig)) / math.pi * (abs(mP * I1) + abs(np.conj(mP) * I2)))
        print(f"   {'δ' if RADIAL else 'τ'}={t:7.4f}: Fehler {err.max():.3e} (bei {np.rad2deg(phis[np.argmax(err)]):.1f}°) | a₁ = {a[0]:+.3e}, a₂ = {a[1]:+.3e} | "
              f"Vorhersage |a₁| ≲ {vorh[0]:.3e}, |a₂| ≲ {vorh[1]:.3e}", flush=True)
