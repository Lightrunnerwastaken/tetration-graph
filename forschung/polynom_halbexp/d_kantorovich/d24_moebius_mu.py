"""D24 — Ursache der Eckmode: komplexes Möbius-Modell mit Fixpunkten L, L̄ und Multiplikator μ an L.
f_μ = M⁻¹(μ·M(x)), M(x) = (x − L)/(x − L̄);  f_μ(L) = L, f_μ'(L) = μ, f_μ'(L̄) = 1/μ.
μ = e^{iθ}: elliptisch (wie Strang Es Möbius-exp, dort ohne Eckmode); |μ| > 1: an L abstossend wie bei exp (μ = √L).
Linearisierung J_N an T_N f_μ (komplex, acb), letzte BL Spalten der Inversen, |L|-umskaliert: σ1, σ2 gegen N.
"""
import sys
import time

import numpy as np
from flint import acb, acb_mat, acb_poly, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
MUS = sys.argv[2].split(";")
BL0 = 100
Lr = 1.3745570107437075


def J_acb(a, N):
    n = N + 1
    f = acb_poly(a)
    P = [acb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(n))
    fpf = acb_poly([0])
    for m in range(1, n):
        fpf += m * a[m] * P[m - 1]
    J = acb_mat(n, n)
    for k in range(n):
        for j in range(n):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return J


t0 = time.time()
for mus in MUS:
    for N in NS:
        n = N + 1
        BL = min(BL0, n)
        ctx.prec = int(0.9 * N) + 300
        L = acb(arb("0.31813150520476413531265425158766451720351761387139986692237"),
                arb("1.33723570143068940890116214319371061253950213846051241687"))
        Lb = L.conjugate()
        mu = eval(mus, {"L": L, "acb": acb, "arb": arb, "sqrt": lambda z: z.sqrt(), "expi": lambda t: acb(0, t).exp()})
        al, be, ga, de = L - mu * Lb, -(L * Lb) * (1 - mu), 1 - mu, mu * L - Lb
        # f = (al x + be)/(ga x + de) = (al x + be) (1/de) Σ (−ga/de)^k x^k
        q = -ga / de
        c, pw = [], acb(1) / de
        for k in range(n):
            c.append(be * pw + (al * (pw / q) if k > 0 else 0))
            pw = pw * q
        pol = -de / ga
        J = J_acb(c, N)
        E = acb_mat(n, BL)
        for cc in range(BL):
            E[N - BL + 1 + cc, cc] = 1
        X = J.solve(E, algorithm="approx")
        k = np.arange(n)
        cols = np.arange(N - BL + 1, n)
        G = np.array([[complex(float(X[i, cc].real.mid()), float(X[i, cc].imag.mid())) for cc in range(BL)]
                      for i in range(n)]) * (Lr ** k)[:, None] / (Lr ** cols)[None, :]
        S = np.linalg.svd(G, compute_uv=False)
        # Kontrolle f(L) = L
        fl = sum((c[i] * L ** i for i in range(min(n, 400))), acb(0))
        print(f"μ={mus:22s} N={N}: |μ|={float(abs(mu).mid()):.4f} arg={float(mu.imag.mid() and np.angle(complex(float(mu.real.mid()), float(mu.imag.mid())))):+.3f} "
              f"Pol |{complex(float(pol.real.mid()), float(pol.imag.mid())):.3f}| f(0)={complex(float(c[0].real.mid()), float(c[0].imag.mid())):.3f} | "
              f"σ1={S[0]:.4f} σ2={S[1]:.4f} σ3={S[2]:.4f}  ({time.time()-t0:.0f}s)", flush=True)
