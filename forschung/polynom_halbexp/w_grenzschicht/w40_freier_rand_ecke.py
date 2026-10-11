"""W40 — Freies Randwertproblem (§15) MIT der lokalen Eckstruktur bei L, L̄ (Deckverschiebung 2πi − log μ, κ₁ = 2πi/(2πi − log μ)).
Gegenüber w33b: (a) Eckfunktionen Re(c_n Z^{nκ₁}) mit Z = e^{iγ}(t − L̃)/(t − P), Schnitt = Kreisbogen L̃ → P IN der Bucht
(jede Potenz ζ^{nκ₁} erfüllt H = H∘h samt Normalableitung auf Γout lokal exakt); (b) Eckwinkel der Bucht frei (B nicht mehr an h(A-Bogen) geheftet).
Unskaliert auf ρ = 1. H = −log|t| + Σq log|t − s| + C + Eckteil;  H = 0 auf B;  Cauchy-Daten H = H∘h auf Γout;  Σq = 0.
Aufruf: python w40_freier_rand_ecke.py [K] [n_ecke] [psi_cut]
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 3
PSI_CUT = float(sys.argv[3]) if len(sys.argv) > 3 else -100.0
QCUT = int(sys.argv[4]) if len(sys.argv) > 4 else 12      # Ladungen erst ab Index QCUT (Abstand zur Ecke)
WD = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3     # Gewicht der Ableitungsbedingung
NG = int(sys.argv[6]) if len(sys.argv) > 6 else 241       # Kollokationspunkte auf Γout
QF = [float(v) for v in sys.argv[7].split(",")] if len(sys.argv) > 7 else [0.2, 0.5, 0.8]   # Ladungslagen (Anteil Buchtbreite)
QS = int(sys.argv[8]) if len(sys.argv) > 8 else 2         # Schrittweite der Ladungen
RC = float(sys.argv[9]) if len(sys.argv) > 9 else 1e-13   # rcond
BS = int(sys.argv[10]) if len(sys.argv) > 10 else 2         # Schrittweite Kollokation auf B
DMU = complex(*[float(v) for v in sys.argv[11].split(',')]) if len(sys.argv) > 11 else 0j   # Störung von μ = f'(L)
coef, Lf = W.kubik_koeff(mu=np.sqrt(W.L + 0j) + DMU)
rho = W.RHO
EPS = float(sys.argv[12]) if len(sys.argv) > 12 else 0.0   # Störung ε·(x² − 2Re L·x + |L|²)² (verschwindet doppelt in L, L̄)
if EPS != 0.0:
    qL = np.polynomial.Polynomial([abs(Lf) ** 2, -2 * Lf.real, 1.0])
    coef = list((np.polynomial.Polynomial(coef) + EPS * qL ** 2).coef)
Ft = np.polynomial.Polynomial([coef[j] * rho ** (j - 1) for j in range(len(coef))])
Lt = Lf / rho
phiL = np.angle(Lt)
a0 = Ft(0).real
R0 = abs(Lt - a0)
mu = Ft.deriv()(Lt)
K1 = 2j * np.pi / (2j * np.pi - np.log(mu))


def h(t):
    r = (Ft - t).roots()
    return r[np.argmin(np.abs(r))]


hv = np.vectorize(h)
# Start-Buchtrand: h(A-Bogen)
psis = np.linspace(np.angle(Lt - a0), -np.angle(Lt - a0), 2001)
B0 = hv(a0 + R0 * np.exp(1j * psis))
o = np.argsort(np.angle(B0))
uu = np.linspace(-1, 1, 401)
phi_grid = phiL * np.tanh(3.0 * uu) / np.tanh(3.0)
r0 = np.interp(phi_grid, np.angle(B0)[o], np.abs(B0)[o])
r0[0] = r0[-1] = 1.0
# Γout-Kollokation (zur Ecke gradiert)
ug = np.linspace(-1, 1, NG)[3:-3]
Gphi = phiL * np.tanh(3.0 * ug) / np.tanh(3.0)
Gc = np.exp(1j * Gphi)
dlt = 1e-4
hG, hGp, hGm = hv(Gc), hv(Gc * (1 + dlt)), hv(Gc * (1 - dlt))
# Eckfunktionen: Schnitt von L̃ in Richtung ψ_cut (relativ zu arg L) als Kreisbogen nach P (in der Bucht)
P = 0.5 * (1 + np.interp(0, phi_grid, r0)) + 0j
gam = np.pi + np.angle(Lt - P) - (phiL + np.deg2rad(PSI_CUT))


def Zl(t):
    return np.exp(1j * gam) * (t - Lt) / (t - P)


def Zu(t):            # Spiegelpartner bei L̄
    return np.conj(Zl(np.conj(t)))


def ecke_cols(t):
    """Spalten (Re c_n, Im c_n) → Re(c Z^{nκ}) + Re(c̄ Z̄^{nκ̄})."""
    cols = []
    z1, z2 = Zl(t), Zu(t)
    for n in range(1, NC + 1):
        w1 = np.exp(n * K1 * np.log(z1))
        w2 = np.exp(n * np.conj(K1) * np.log(z2))
        cols += [w1.real + w2.real, -w1.imag + w2.imag]
    if not cols:
        return np.zeros((len(np.atleast_1d(t)), 0))
    return np.array(cols).T


def rB_of(a):
    x = phi_grid / phiL
    w = np.exp(sum(a[k] * np.cos(k * np.pi * x / 2) ** 2 for k in range(len(a))))   # am Rand ≠ 1 ⇒ Eckwinkel frei
    return 1 - w * (1 - r0)


def baue(a, ret_all=False):
    rB = rB_of(a)
    Bc = (rB * np.exp(1j * phi_grid))[3:-3:BS]
    srcs = [((rB + f * (1 - rB)) * np.exp(1j * phi_grid))[QCUT:-QCUT:QS] for f in QF]
    src = np.concatenate(srcs)

    def Lm(p):
        return np.log(np.abs(p[:, None] - src[None, :]))

    def cols(p):
        return np.hstack([Lm(p), np.ones((len(p), 1)), ecke_cols(p)])
    MB = cols(Bc)
    bB = np.log(np.abs(Bc))
    MG = cols(Gc) - cols(hG)
    bG = np.log(np.abs(Gc)) - np.log(np.abs(hG))
    MD = (cols(Gc * (1 + dlt)) - cols(Gc * (1 - dlt))) / (2 * dlt) - (cols(hGp) - cols(hGm)) / (2 * dlt)
    bD = (np.log(np.abs(Gc * (1 + dlt))) - np.log(np.abs(Gc * (1 - dlt)))) / (2 * dlt) - (np.log(np.abs(hGp)) - np.log(np.abs(hGm))) / (2 * dlt)
    nc = MB.shape[1]
    MS = np.zeros((1, nc))
    MS[0, :len(src)] = 1
    Mat = np.vstack([MB, MG, WD * MD, 1e3 * MS])
    rhs = np.concatenate([bB, bG, WD * bD, [0.0]])
    sc = np.linalg.norm(Mat, axis=0) + 1e-30
    sol, *_ = np.linalg.lstsq(Mat / sc, rhs, rcond=RC)
    sol = sol / sc
    res = Mat @ sol - rhs
    if ret_all:
        return sol, src, res, rB
    return res


def H_fun(sol, src):
    def H(t):
        t = np.atleast_1d(np.asarray(t, complex))
        M = np.hstack([np.log(np.abs(t[:, None] - src[None, :])), np.ones((len(t), 1)), ecke_cols(t)])
        return -np.log(np.abs(t)) + M @ sol
    return H


a_init = np.zeros(K)
r_init = baue(a_init)
fit = least_squares(baue, a_init, diff_step=1e-3, max_nfev=300)
sol, src, res, rB = baue(fit.x, True)
H = H_fun(sol, src)
if len(sys.argv) > 13 and sys.argv[13] == "jac":
    Jm = fit.jac
    sv, Vh = np.linalg.svd(Jm, full_matrices=False)[1:]
    print("Jacobi-Singulärwerte (B-Parameter):", np.array2string(sv, precision=3))
    v = Vh[-1]
    print("schwächste Richtung:", np.array2string(v, precision=3))
    for eps_ in (0.05, 0.1, 0.2):
        for sgn in (+1, -1):
            sol2, src2, res2, rB2 = baue(fit.x + sgn * eps_ * v, True)
            H2 = H_fun(sol2, src2)
            print(f"  Schritt {sgn*eps_:+.2f}: Rest {np.linalg.norm(res2):.3e} (Basis {np.linalg.norm(res):.3e}), c = {H2(-1.0)[0]:.5f}, Bucht(0°) = {np.interp(0, phi_grid, rB2):.4f}")
nq = len(src)
cc = sol[nq + 1:]
nB, nG = len(phi_grid[3:-3:2]), len(Gc)
rb_, rg_, rd_ = res[:nB], res[nB:nB + nG], res[nB + nG:nB + 2 * nG]
print(f"Rest: B {np.linalg.norm(rb_):.2e} (max bei φ = {np.degrees(phi_grid[3:-3:2][np.argmax(abs(rb_))]):.1f}°), Γout-Wert {np.linalg.norm(rg_):.2e} "
      f"(max bei {np.degrees(Gphi[np.argmax(abs(rg_))]):.1f}°), Γout-Ableitung {np.linalg.norm(rd_):.2e} (max bei {np.degrees(Gphi[np.argmax(abs(rd_))]):.1f}°); "
      f"erste Ladung {np.degrees(phi_grid[QCUT]):.2f}° (L bei {np.degrees(phiL):.2f}°)")
Bm = (rB * np.exp(1j * phi_grid))[4:-4:2]
Gm = np.exp(1j * 0.5 * (Gphi[1:] + Gphi[:-1]))
hGm_ = hv(Gm)
dl = 1e-5
dH = lambda z: (H(z * (1 + dl)) - H(z * (1 - dl))) / (2 * dl)
dHh = lambda z: np.array([(H(h(zz * (1 + dl)))[0] - H(h(zz * (1 - dl)))[0]) / (2 * dl) for zz in z])
print(f"Validierung zwischen Kollokationspunkten: max|H| auf B {np.abs(H(Bm)).max():.2e}, max|H − H∘h| auf Γout {np.abs(H(Gm) - H(hGm_)).max():.2e}, "
      f"max|∂_r(H − H∘h)| {np.abs(dH(Gm) - dHh(Gm)).max():.2e}; max|q| = {np.abs(sol[:len(src)]).max():.1f}, H(0.80) = {H(0.80)[0]:+.4f}")
print(f"K = {K}, Eckterme n ≤ {NC}, ψ_cut = {PSI_CUT}°: Rest {np.linalg.norm(r_init):.3e} → {np.linalg.norm(res):.3e}; B-Parameter {np.round(fit.x, 4)}")
# Schnitt in der Bucht? Punkte t(u) = (L̃ + u e^{−iγ}P)/(1 + u e^{−iγ}), u ∈ (0, ∞)
uu_ = np.geomspace(1e-4, 1e4, 400)
tc = (Lt + uu_ * np.exp(-1j * gam) * P) / (1 + uu_ * np.exp(-1j * gam))
phc = np.angle(tc)
rBc = np.interp(np.abs(phc), phi_grid[phi_grid >= 0], rB[phi_grid >= 0])
inn = (np.abs(tc) < 1) & (np.abs(tc) > rBc) & (np.abs(phc) <= phiL)
print(f"Schnitt L̃ → P in der Bucht: {inn.all()} (Anteil {inn.mean():.3f}; P = {P.real:.4f})")
np.savez(W.HIER + f"/w40_modell_K{K}_n{NC}_c{int(abs(PSI_CUT))}.npz", sol=sol, src=src, rB=rB, phi_grid=phi_grid, gam=gam, P=P, K1=K1, NC=NC)
print(f"Eckkoeffizienten c_n (ρ = 1): " + ", ".join(f"{cc[2*k]:+.4f}{cc[2*k+1]:+.4f}i" for k in range(NC)))
print(f"Buchtrand(0°) = {np.interp(0, phi_grid, rB):.4f} (gemessen 0.815);  c = H(−1) = {H(-1.0)[0]:+.5f} (gemessen 0.0120)")
print("Profil Γin: " + " ".join(f"{t}°:{H(np.exp(1j*(phiL+np.deg2rad(t))))[0]:+.4f}" for t in (1, 2, 4, 8, 15, 30, 60, 103.38)))
print("   gemessen: 1°:+0.0003 2°:+0.0009 4°:+0.0024 8°:+0.0049 15°:+0.0067 30°:+0.0091 60°:+0.0108 103°:+0.0116")
hh = 1e-5
print("Stegfrequenz: " + " ".join(f"{p}°:{-(H((1+hh)*np.exp(1j*np.deg2rad(p)))[0]-H((1-hh)*np.exp(1j*np.deg2rad(p)))[0])/(2*hh):.4f}" for p in (90, 120, 180)))
print("H(1.05), H(1.2), H(1.3) (gemessen 0.30/0.117/0.010 bei N=600): " + " ".join(f"{H(x)[0]:.4f}" for x in (1.05, 1.2, 1.3)))
# lokale Knoten um L (Radien wie w38, unskaliert 0.03/0.06/0.12)
for rr in (0.005, 0.01, 0.03, 0.06, 0.12):
    ps = np.arange(-180, 180.1, 0.5)
    tt = Lt + (rr / rho) * np.exp(1j * (phiL + np.deg2rad(ps)))
    inb = np.array([(abs(t) < 1) and (abs(t) > np.interp(abs(np.angle(t)), phi_grid[phi_grid >= 0], rB[phi_grid >= 0])) and abs(np.angle(t)) < phiL for t in tt])
    Hv = H(tt)
    Hv[inb] = np.nan
    zs = [ps[i] + 0.5 * Hv[i] / (Hv[i] - Hv[i + 1]) for i in range(len(ps) - 1)
          if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1]) and Hv[i] * Hv[i + 1] < 0]
    print(f"   r = {rr}: Nullstellen ψ = " + ", ".join(f"{z:.1f}" for z in zs) + f";  H(ψ=0) = {H(Lt + (rr/rho)*np.exp(1j*phiL))[0]:+.5f}, "
          f"H(180) = {H(Lt - (rr/rho)*np.exp(1j*phiL))[0]:+.5f}  (gemessen w38: siehe w39)")
