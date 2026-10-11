"""Prüft, ob Ds d28/d29 (acb_series ohne ctx.cap) vom 10-Term-Abschnitt betroffen sind: |⟨u1, J⁻¹P_N J∞k_η⟩| mit und ohne cap."""
import sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx
import w_kern as W
N = 400; n = N + 1; ctx.prec = int(0.9 * N) + 300
coef, _ = W.kubik_koeff(); rho = W.RHO; L = W.L
d = np.load(W.HIER + f"/w01_N{N}.npz"); U = d["U"]
a = W.kubik_a(N, coef); P = W.potenzen(a, N); J = W.J_mat(a, N, P)
x0 = acb_poly([acb(c) for c in a])
mser = acb_poly([0])
for j in range(1, n): mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
mser = mser.truncate(n)
pole = [("u=0.003", L + 0.003*np.exp(1j*np.deg2rad(87.15))), ("u=0.03", L + 0.03*np.exp(1j*np.deg2rad(87.15))), ("τ=4°", rho*np.exp(1j*(np.angle(L)+np.deg2rad(4)))), ("L", L)]
for cap in (10, n + 10):
    ctx.cap = cap
    E = arb_mat(n, 2*len(pole))
    for i,(nm, z) in enumerate(pole):
        zz = acb(z.real, z.imag)
        s1 = acb_series((acb_poly([zz]) - x0).coeffs(), prec=n).inv().coeffs(); s1 = list(s1) + [acb(0)]*(n-len(s1))
        gl = [zz**(-(k+1)) for k in range(n)]
        s2 = (mser * acb_poly(gl)).truncate(n).coeffs(); s2 = list(s2) + [acb(0)]*(n-len(s2))
        for k in range(n):
            w = s1[k] + s2[k]; E[k, 2*i] = w.real; E[k, 2*i+1] = w.imag
    X = J.solve(E, algorithm="approx")
    sc = rho ** np.arange(n)
    out = []
    for i,(nm, z) in enumerate(pole):
        resp = (np.array([float(X[r,2*i].mid()) for r in range(n)]) + 1j*np.array([float(X[r,2*i+1].mid()) for r in range(n)])) * sc
        out.append(f"{nm}: {abs(U[:,0] @ resp):.4g}")
    print(f"cap={cap}: " + "  ".join(out))
