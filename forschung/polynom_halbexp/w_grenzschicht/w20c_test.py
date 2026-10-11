import numpy as np
from flint import acb, acb_series, arb, arb_mat, acb_poly, ctx
import w_kern as W
N = 300; n = N + 1; ctx.prec = int(1.5*N) + 400
coef, _ = W.kubik_koeff(); rho = W.RHO; L = W.L
a = W.kubik_a(N, coef)
J = W.J_mat(a, N)
fpol = acb_poly([acb(float(c)) for c in coef])
fp = fpol.derivative()
mser = acb_poly([0])
P = acb_poly([1]); fpcoeffs = fp.coeffs()
for i, c in enumerate(fpcoeffs):
    mser += c * P; P = (P * fpol).truncate(n)
mser = mser.truncate(n)
et = acb(L.real, L.imag)
s1 = acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs()
kz = [et ** (-(j + 1)) for j in range(n)]
s2 = (mser * acb_poly(kz)).truncate(n).coeffs()
s1 = list(s1) + [acb(0)]*(n-len(s1)); s2 = list(s2) + [acb(0)]*(n-len(s2)); TJk = [s1[j] + s2[j] for j in range(n)]
# J T_N k
Jk = [sum((J[j, i] * kz[i] for i in range(n)), acb(0)) for j in range(n)]
g2 = [TJk[j] - Jk[j] for j in range(n)]
fN1 = acb_poly([1])
for _ in range(N + 1): fN1 = (fN1 * fpol).truncate(n)
inv = acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv()
ser = (fN1 * acb_poly(inv.coeffs())).truncate(n).coeffs()
g1 = [ser[j] / et ** (N + 1) for j in range(n)]
for j in (0, 1, 10, 100, 200, 290, 300):
    print(j, complex(float(g1[j].real.mid()), float(g1[j].imag.mid())) * rho**j, complex(float(g2[j].real.mid()), float(g2[j].imag.mid())) * rho**j)
# Direkt: Π_N k = J^{-1} T_N J∞ k ; y = Π_N k − T_N k
E = arb_mat(n, 4)
for j in range(n):
    E[j,0] = TJk[j].real; E[j,1] = TJk[j].imag; E[j,2] = g1[j].real; E[j,3] = g1[j].imag
X = J.solve(E, algorithm="approx")
Pi = np.array([complex(float(X[j,0].mid()), float(X[j,1].mid())) for j in range(n)])
Yg = np.array([complex(float(X[j,2].mid()), float(X[j,3].mid())) for j in range(n)])
kzf = np.array([complex(float(kz[j].real.mid()), float(kz[j].imag.mid())) for j in range(n)])
sc = rho ** np.arange(n)
print("|Pi_N k_L| =", np.linalg.norm(Pi*sc), " |y_L| = |Pi - T_N k| =", np.linalg.norm((Pi-kzf)*sc), " |J^{-1} g1| =", np.linalg.norm(Yg*sc))
print("Prüfe J X = E Rest:", max(abs(float((J*X - E)[j, c].mid())) for j in range(0, n, 37) for c in range(4)))
