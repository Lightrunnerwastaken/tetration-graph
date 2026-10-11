"""M13d — Blockstruktur der Form Re<𝒯x,x> (H²): c_lo(M), c_hi(M) (Frequenzen M+1..M'), Kopplung δ(M)."""
import sys, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
MS = [int(x) for x in sys.argv[1].split(",")]
MP = int(sys.argv[2]); G = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg); c0 = 1 - 2*eps*np.cos(th)
f = lambda z: eps + c0*z + eps*z*z
fp = lambda z: c0 + 2*eps*z
h = lambda w: (-c0 + np.sqrt(c0*c0 - 4*eps*(eps - w) + 0j))/(2*eps)
def bogen(a, b, n):
    t = (np.arange(n) + 0.5)/n
    return a + (b-a)*t, (b-a)/n
nin = int(G*(1 - th/np.pi)); nout = G - nin
phin, win = bogen(th, 2*np.pi - th, nin); phout, wout = bogen(-th, th, nout)
zin, zout = np.exp(1j*phin), np.exp(1j*phout)
gin, min_ = f(zin), fp(f(zin)); gout, fpo = h(zout), fp(zout)
k = np.arange(MP+1)
Kin = (gin[:, None]**k)/min_[:, None]; Kout = fpo[:, None]*gout[:, None]**k
Ein, Eout = zin[:, None]**k, zout[:, None]**k
# Matrix von P_{0..MP} 𝒦 auf H²-Frequenzen: Kmat[l,k] = <𝒦 z^k, z^l>
Kmat = (Ein.conj().T @ (Kin*(win/(2*np.pi)))) + (Eout.conj().T @ (Kout*(wout/(2*np.pi))))
B = np.eye(MP+1) + Kmat
H = 0.5*(B + B.conj().T)
print(f"M'={MP}: λ_min(Re P𝒯P) gesamt = {np.linalg.eigvalsh(H)[0]:.6f}", flush=True)
for M in MS:
    lo = slice(0, M+1); hi = slice(M+1, MP+1)
    clo = np.linalg.eigvalsh(H[lo, lo])[0]
    chi = np.linalg.eigvalsh(H[hi, hi])[0]
    d1 = np.linalg.norm(Kmat[hi, lo], 2); d2 = np.linalg.norm(Kmat[lo, hi], 2)
    # Kopplung in der Form: |Re<Ky,z>+Re<Kz,y>| <= ||(K_hl + K_lh^*)/2|| * 2 ||y|| ||z||
    dsym = np.linalg.norm(0.5*(Kmat[hi, lo] + Kmat[lo, hi].conj().T), 2)
    print(f"  M={M}: c_lo = {clo:.6f}, c_hi = {chi:.6f}, ‖P_hi𝒦P_lo‖ = {d1:.4f}, ‖P_lo𝒦P_hi‖ = {d2:.4f}, "
          f"symm. Kopplung ‖(K_hl+K_lh*)/2‖ = {dsym:.4f}", flush=True)
