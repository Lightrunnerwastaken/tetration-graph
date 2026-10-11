"""M09 — Vorhersage-Test: Eingabe mit schwacher Singularität (1 − z/P)^σ an einem Randpunkt P.
  P in der Ecke (P = e^{±iθ}, reelle Kombination)  gegen  P auf dem D-Bogen (P = −1).
Erwartung (Bericht §5): Ecke → Konvergenz; D-Bogen → Eckanteil σ₁·V₁(P)·N^{−σ−1} wächst ⇒ Divergenz für grosse N.
Lineares Problem im normierten radialen Quadrat-Modell: Π_N u = J_N⁻¹ T_N[u∘f + m u], Fehler sup_{|z|=1}|Π_N u − u|.
Aufruf: python m09_singularitaet_ort.py eps theta_grad sigma "N1,N2,..."
"""
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, ctx

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
eps = float(sys.argv[1])
thg = float(sys.argv[2])
sig = float(sys.argv[3])
NS = [int(x) for x in sys.argv[4].split(",")]
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)
coef = [eps, c0, eps]


def pot_reihe(A, s, n):
    """Taylor-Koeffizienten von A(z)^s (A Polynom-Liste acb, A[0] ≠ 0) bis Grad n−1: n g_n A0 = Σ_k (s k − (n−k)) A_k g_{n−k}."""
    g = [A[0] ** s]
    for nn in range(1, n):
        acc = acb(0)
        for k in range(1, min(len(A) - 1, nn) + 1):
            acc += (s * k - (nn - k)) * A[k] * g[nn - k]
        g.append(acc / (nn * A[0]))
    return g


def eingabe(P, N, reell_kombi):
    """Koeffizienten (bis Grad N) von u und von T_N[u∘f + m u] für u = (1 − z/P)^σ (+ konj. falls reell_kombi)."""
    Pa = acb(P.real, P.imag)
    s = arb(sig)
    u = pot_reihe([acb(1), -1 / Pa], s, N + 1)
    # u∘f: A(z) = 1 − f(z)/P
    A = [1 - acb(coef[0]) / Pa, -acb(coef[1]) / Pa, -acb(coef[2]) / Pa]
    uf = pot_reihe(A, s, N + 1)
    m = [acb(c0 + 2 * eps * coef[0]), acb(2 * eps * coef[1]), acb(2 * eps * coef[2])]
    mu_ = [sum((m[i] * u[k - i] for i in range(3) if k - i >= 0), acb(0)) for k in range(N + 1)]
    rhs = [uf[k] + mu_[k] for k in range(N + 1)]
    if reell_kombi:
        u = [2 * x.real for x in u]
        rhs = [2 * x.real for x in rhs]
    else:
        u = [x.real for x in u]
        rhs = [x.real for x in rhs]
    return u, rhs


phis = np.linspace(0, 2 * np.pi, 4001)
zs = np.exp(1j * phis)
R0 = abs(np.exp(1j * th) - eps)
fall = [("Ecke e^{±iθ}", np.exp(1j * th), True), ("D-Bogen −1", complex(-1, 0), False),
        (f"A-Bogen a0+R0={eps + R0:.4f}", complex(eps + R0, 0), False)]
for N in NS:
    t0 = time.time()
    ctx.prec = int(0.9 * N) + 300
    J = K.jmat(coef, N, reell=True)
    for name, P, rk in fall:
        u, rhs = eingabe(P, N, rk)
        X = J.solve(arb_mat(N + 1, 1, rhs), algorithm="approx")
        p = np.array([float(X[k, 0].mid()) for k in range(N + 1)])
        # u exakt auf dem Kreis: (1 − z/P)^σ (+ konj)
        if rk:
            uex = (1 - zs / P) ** sig + (1 - zs / np.conj(P)) ** sig
        else:
            uex = (1 - zs / P) ** sig
        err = np.abs(np.polyval(p[::-1], zs) - uex)
        j = int(np.argmax(err))
        # auch im Inneren r = 0.9
        zi = 0.9 * zs
        ui = ((1 - zi / P) ** sig + (1 - zi / np.conj(P)) ** sig) if rk else (1 - zi / P) ** sig
        erri = np.max(np.abs(np.polyval(p[::-1], zi) - ui))
        print(f"N={N} {name}: sup|Π_N u − u| auf |z|=1: {err.max():.3e} (bei {np.rad2deg(phis[j]):.1f}°), auf |z|=0.9: {erri:.3e}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
