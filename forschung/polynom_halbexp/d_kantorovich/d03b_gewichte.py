"""D03b — Auswertung der affin-invarianten Kantorovich-Grösse h = α·ω für Gewichtsfamilien und optimale Gewichte.

log h(u) mit u_k = log w_k ist konvex (log-sum-exp minus linear, Maximum davon, plus log-sum-exp) -> globales Minimum.
Aufruf: python d03b_gewichte.py pad 10,20,...  (liest d03_daten/<modus>_<N>.npz)
"""
import json
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

sys.stdout.reconfigure(encoding="utf-8")
MODUS = sys.argv[1] if len(sys.argv) > 1 else "pad"
NS = [int(x) for x in sys.argv[2].split(",")]
absL = 1.3745570107437075
OPT = len(sys.argv) <= 3 or sys.argv[3] != "ohne_opt"


def lade(N):
    D = np.load(f"d03_daten/{MODUS}_{N}.npz")
    return D["Z"], D["d"], D["paare"]


def h_von(u, A, ld, I, Jx):
    la = logsumexp(ld + u)
    s = logsumexp(A + u[:, None], axis=0) - u[I] - u[Jx]
    c = int(np.argmax(s))
    return la, s[c], c


def optimiere(A, ld, I, Jx, u0, tau=0.05):
    n = len(u0)

    def fg(u):
        # glatt: log α + τ·LSE(s/τ)
        la_vec = ld + u
        la = logsumexp(la_vec)
        ga = np.exp(la_vec - la)
        M = A + u[:, None]
        lse = logsumexp(M, axis=0)
        s = lse - u[I] - u[Jx]
        sm = tau * logsumexp(s / tau)
        pc = np.exp(s / tau - sm / tau)              # Gewichte der Paare
        W = np.exp(M - lse[None, :])                 # d lse_c / d u_k
        g = ga + W @ pc
        np.subtract.at(g, I, pc)
        np.subtract.at(g, Jx, pc)
        return la + sm, g

    r = minimize(fg, u0, jac=True, method="L-BFGS-B", options={"maxiter": 3000})
    return r.x


aus = {}
for N in NS:
    Z, d, paare = lade(N)
    with np.errstate(divide="ignore"):
        A = np.log(np.abs(Z))
        ld = np.log(np.abs(d))
    I, Jx = paare[:, 0], paare[:, 1]
    k = np.arange(N + 1)
    zeile = {}
    fam = {f"L,g={g}": k * np.log(absL) + g * np.log1p(k) for g in (-2, 0, 1, 2, 3, 4)}
    fam.update({f"rho={r}": k * np.log(r) for r in (1.0, 1.2, 1.3, 1.45)})
    for name, u in fam.items():
        la, lo, c = h_von(u, A, ld, I, Jx)
        zeile[name] = (float(np.exp(la)), float(np.exp(lo)), float(np.exp(la + lo)), [int(I[c]), int(Jx[c])])
    if OPT:
        u = fam["L,g=0"].copy()
        for tau in (0.3, 0.1, 0.03, 0.01):
            u = optimiere(A, ld, I, Jx, u, tau)
        la, lo, c = h_von(u, A, ld, I, Jx)
        zeile["optimal"] = (float(np.exp(la)), float(np.exp(lo)), float(np.exp(la + lo)), [int(I[c]), int(Jx[c])])
        u = u - u[0]
        zeile["opt_u"] = [float(x) for x in u]
    aus[N] = zeile
    print(f"N={N}")
    for name, z in zeile.items():
        if name != "opt_u":
            print(f"   {name:9s} α={z[0]:.3e} ω={z[1]:.3e} h={z[2]:.3e} bei {z[3]}")
    if OPT:
        uu = np.array(zeile["opt_u"])
        # effektiver Radius der optimalen Gewichte: Steigung von u
        print("   opt. log-Gewicht − k·log|L| bei k=0,N/4,N/2,3N/4,N-1,N:",
              [round(float(uu[m] - m * np.log(absL)), 2) for m in (0, N // 4, N // 2, 3 * N // 4, N - 1, N)])
json.dump({str(k): v for k, v in aus.items()}, open(f"d03b_{MODUS}.json", "w"))
