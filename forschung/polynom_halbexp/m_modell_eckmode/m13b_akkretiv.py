"""M13b — Ist 𝒯 akkretiv? λ_min(Re ⟨𝒯x, x⟩) auf P_N; dazu σ_min(𝒯|P_N) bis N = 800 (float, Gitter)."""
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
    t0 = time.time()
    k = np.arange(N+1)
    Ain = zin[:, None]**k + (gin[:, None]**k)/min_[:, None]
    Aout = zout[:, None]**k + fpo[:, None]*gout[:, None]**k
    Ein, Eout = zin[:, None]**k, zout[:, None]**k
    # B_{lk} = <𝒯 z^k, z^l> = (1/2π)∫ (𝒯z^k) conj(z^l) dφ
    B = (Ein.conj().T @ (Ain*(win/(2*np.pi)))) + (Eout.conj().T @ (Aout*(wout/(2*np.pi))))
    H = 0.5*(B + B.conj().T)
    lam = np.linalg.eigvalsh(H)
    A = np.vstack([Ain*np.sqrt(win/(2*np.pi)), Aout*np.sqrt(wout/(2*np.pi))])
    S = np.linalg.svd(A, compute_uv=False)
    print(f"N={N}: λ_min(Re P𝒯P) = {lam[0]:.6f}, λ_max = {lam[-1]:.4f}; σ_min(𝒯|P_N) = {S[-1]:.6f} ({time.time()-t0:.0f}s)", flush=True)
