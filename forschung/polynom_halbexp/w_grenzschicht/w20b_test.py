import sys, numpy as np
from flint import acb, acb_series, arb, arb_mat, acb_poly, acb_mat, ctx
import w_kern as W
N = 300; n = N + 1; ctx.prec = int(1.5*N) + 400
coef, _ = W.kubik_koeff(); rho = W.RHO; L = W.L
d = np.load(W.HIER + f"/w01_N{N}.npz"); u1 = d["U"][:, 0]; S = d["S"]
J = W.J_mat(W.kubik_a(N, coef), N)
fpol = acb_poly([acb(float(c)) for c in coef])
fN1 = acb_poly([1])
for _ in range(N+1): fN1 = (fN1 * fpol).truncate(n)
for eta in [L, L*1.05, -rho*1.0]:
    et = acb(eta.real, eta.imag)
    inv = acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv()
    ser = (fN1 * acb_poly(inv.coeffs())).truncate(n).coeffs(); ser = ser + [acb(0)]*(n-len(ser))
    g = [ser[j] / et**(N+1) for j in range(n)]
    # direkte Lösung, Real- und Imaginärteil getrennt
    Er = arb_mat(n, 2)
    for j in range(n): Er[j,0] = g[j].real; Er[j,1] = g[j].imag
    X = J.solve(Er, algorithm="approx")
    resp = np.array([complex(float(X[j,0].mid()), float(X[j,1].mid())) * rho**j for j in range(n)])
    gn = np.array([complex(float(g[j].real.mid()), float(g[j].imag.mid())) * rho**j for j in range(n)])
    print(f"eta={eta:.3f}: |g|_rescaled = {np.linalg.norm(gn):.3e}, |y_eta| (rescaled l2) = {np.linalg.norm(resp):.3e}, l1_rho = {np.sum(np.abs(resp)):.3e}, Eckanteil <u1,resp> = {abs(u1@resp):.3e}")
