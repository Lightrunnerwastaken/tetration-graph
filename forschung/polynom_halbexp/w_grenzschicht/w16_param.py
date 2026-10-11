"""W16 — Kubisches Modell mit anderem Multiplikator μ = |μ| e^{iβ} (Fixpunkt bei arg L): volle Inverse, SVD, Speichern.
Aufruf: python w16_param.py TAG MU_ABS MU_ARG N1,N2,...   → w16_{TAG}_N{N}.npz  (S, U[:,:4], Vt[:4], coef)
"""
import sys, time
import numpy as np
from flint import ctx
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
tag, mabs, marg = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
mu = mabs * np.exp(1j * marg)
coef, Lf = W.kubik_koeff(mu=mu)
F = np.polynomial.Polynomial(coef)
print(f"{tag}: μ = {mu:.4f}, coef = {coef}, f(L)−L = {abs(F(Lf)-Lf):.1e}, f'(L) = {F.deriv()(Lf):.4f}", flush=True)
t0 = time.time()
for N in [int(x) for x in sys.argv[4].split(",")]:
    ctx.prec = int(0.9 * N) + 300
    G = W.inverse_umskaliert(W.J_mat(W.kubik_a(N, coef), N), N)
    U, S, Vt = np.linalg.svd(G)
    np.savez(W.HIER + f"/w16_{tag}_N{N}.npz", S=S, U=U[:, :4], Vt=Vt[:4], coef=np.array(coef), mu=mu)
    print(f"N={N} ({time.time()-t0:.0f}s): σ = {S[0]:.5g} {S[1]:.5g} {S[2]:.5g}", flush=True)
