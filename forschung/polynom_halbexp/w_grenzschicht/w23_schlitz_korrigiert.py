"""W23 — Nachmessung von Ds Schlitz-/Eckwerten (d29) MIT ctx.cap (acb_series schnitt sonst bei 10 Termen ab).
Eckanteil |⟨u1, J⁻¹P_N J∞k_η⟩| (mit Kopf, wie D) und kopffrei |⟨u1, y_η⟩| (y_η = Π_N k_η − T_N k_η), kubisch.
Aufruf: python w23_schlitz_korrigiert.py N1,N2,...
"""
import sys, time, json
import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff(); rho = W.RHO; L = W.L; B = 100
res = {}
for N in [int(x) for x in sys.argv[1].split(",")]:
    n = N + 1; ctx.prec = int(0.9 * N) + 300; ctx.cap = n + 10
    t0 = time.time()
    a = W.kubik_a(N, coef); P = W.potenzen(a, N); J = W.J_mat(a, N, P)
    x0 = acb_poly([acb(c) for c in a])
    mser = acb_poly([0])
    for j in range(1, n): mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
    mser = mser.truncate(n)
    pole = [("L", L)] + [(f"u={u}", L + u * np.exp(1j * np.deg2rad(87.15))) for u in (0.001, 0.003, 0.01, 0.03, 0.1)] \
        + [(f"tau={t}", rho * np.exp(1j * (np.angle(L) + np.deg2rad(t)))) for t in (2, 8, 24)]
    E = arb_mat(n, B + 2 * len(pole)); gls = []
    for c in range(B): E[N - B + 1 + c, c] = 1
    for i, (nm, z) in enumerate(pole):
        zz = acb(z.real, z.imag)
        s1 = list(acb_series((acb_poly([zz]) - x0).coeffs(), prec=n).inv().coeffs()); s1 += [acb(0)] * (n - len(s1))
        gl = [zz ** (-(k + 1)) for k in range(n)]
        s2 = list((mser * acb_poly(gl)).truncate(n).coeffs()); s2 += [acb(0)] * (n - len(s2))
        for k in range(n):
            w = s1[k] + s2[k]; E[k, B + 2 * i] = w.real; E[k, B + 2 * i + 1] = w.imag
        gls.append(np.array([complex(float(g.real.mid()), float(g.imag.mid())) for g in gl]))
    X = J.solve(E, algorithm="approx")
    k = np.arange(n); sc = rho ** k; cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * sc[:, None] / (rho ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    z = {"s1": float(S[0])}
    for i, (nm, zc) in enumerate(pole):
        resp = (np.array([float(X[r, B + 2 * i].mid()) for r in range(n)]) + 1j * np.array([float(X[r, B + 2 * i + 1].mid()) for r in range(n)])) * sc
        kz = gls[i] * sc
        z[nm] = {"mitKopf": float(abs(U[:, 0] @ resp)), "kopffrei": float(abs(U[:, 0] @ (resp - kz)))}
    res[N] = z
    json.dump(res, open(W.HIER + "/w23_schlitz.json", "w"), indent=1)
    print(f"N={N} ({time.time()-t0:.0f}s) σ1={S[0]:.3f} | " + " ".join(f"{nm}: {z[nm]['mitKopf']:.3g}/{z[nm]['kopffrei']:.3g}" for nm, _ in pole), flush=True)
