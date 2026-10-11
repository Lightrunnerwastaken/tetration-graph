"""M13c — Akkretivität mit Poisson-Fortsetzung auf ganz L²(∂D) (Frequenzen −N..N)?"""
import sys, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
M = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg); c0 = 1 - 2*eps*np.cos(th)
f = lambda z: eps + c0*z + eps*z*z
fp = lambda z: c0 + 2*eps*z
h = lambda w: (-c0 + np.sqrt(c0*c0 - 4*eps*(eps - w) + 0j))/(2*eps)
def bogen(a, b, n):
    t = (np.arange(n) + 0.5)/n
    return a + (b-a)*t, (b-a)/n
nin = int(M*(1 - th/np.pi)); nout = M - nin
phin, win = bogen(th, 2*np.pi - th, nin); phout, wout = bogen(-th, th, nout)
zin, zout = np.exp(1j*phin), np.exp(1j*phout)
gin, min_ = f(zin), fp(f(zin)); gout, fpo = h(zout), fp(zout)
for N in NS:
    k = np.arange(-N, N+1)
    def harm(p, kk):  # Poisson-Fortsetzung von e^{ikψ} in p: p^k (k≥0), conj(p)^{|k|} (k<0)
        return np.where(kk[None, :] >= 0, p[:, None]**np.maximum(kk, 0)[None, :], np.conj(p)[:, None]**np.maximum(-kk, 0)[None, :])
    Ain = zin[:, None]**k + harm(gin, k)/min_[:, None]
    Aout = zout[:, None]**k + fpo[:, None]*harm(gout, k)
    Ein, Eout = zin[:, None]**k, zout[:, None]**k
    B = (Ein.conj().T @ (Ain*(win/(2*np.pi)))) + (Eout.conj().T @ (Aout*(wout/(2*np.pi))))
    H = 0.5*(B + B.conj().T)
    lam = np.linalg.eigvalsh(H)
    # analytischer Teil zum Vergleich
    idx = k >= 0
    lamH = np.linalg.eigvalsh(H[np.ix_(idx, idx)])
    print(f"N={N}: L²: λ_min(Re) = {lam[0]:.6f}; nur H² (k≥0): {lamH[0]:.6f}", flush=True)
