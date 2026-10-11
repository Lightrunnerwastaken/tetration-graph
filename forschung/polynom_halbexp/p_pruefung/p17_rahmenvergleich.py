"""Strang P, Prüfung 17: Eck-Bedingung im W̄- und im K'-Rahmen für das gemessene Γin-Profil (W 00:24: λ − λ(0) = 0.0014/0.0024/
0.0032/0.0054/0.0089/0.0121 bei τ = 2/3/4/8/24/103°, lokaler Exponent ≈ 1.27 ≈ 1/α'), λ = 0 auf dem anderen Bogen (A-Bogen bzw. f_K(Γout)).
Grösse: sup_Schlitz ω_ζ(λ)/g(ζ) bis u → 0 (Ecken-Grenzwert). Erwartung: W̄ (Eckenexponent 1/α = 1.1326) endlich; K' (1/α' = 1.2704) logarithmisch.
"""
import sys

import numpy as np

import kprime as K

sys.stdout.reconfigure(encoding="utf-8")
L, rho = K.L, K.rho
a0 = 0.49856328794111443
R0 = abs(L - a0)
Lb = np.conj(L)
beta = np.angle((a0 + R0 - L) / (a0 + R0 - Lb))
om = np.angle((-rho - L) / (-rho - Lb)) - beta


def tauW(z):
    w = (z - L) / (z - Lb) * np.exp(-1j * beta)
    ang = np.angle(w)
    if ang < -1e-12:
        ang += 2 * np.pi
    return abs(w) ** (np.pi / om) * np.exp(1j * ang * np.pi / om)


T0W = tauW(1e14)
T0K = K.tau0
taus_deg = np.array([0.0, 2, 3, 4, 8, 24, 103.38])
lam_W = np.array([0.0, 0.0014, 0.0024, 0.0032, 0.0054, 0.0089, 0.0121])


def lam_of_tau(tdeg, p_lokal=1.27):
    """Interpolation; unter 2° Potenzgesetz mit Exponent p_lokal."""
    tdeg = np.asarray(tdeg)
    out = np.interp(tdeg, taus_deg[1:], lam_W[1:])
    small = tdeg < 2
    out[small] = 0.0014 * (tdeg[small] / 2) ** p_lokal
    return out


def setup(rahmen, p_lokal):
    tg = np.concatenate([np.logspace(-9, np.log10(103.38), 400)])
    if rahmen == "W":
        s = np.array([tauW(rho * np.exp(1j * (np.angle(L) + np.radians(t)))).real for t in tg])
    else:
        s = np.array([K.s_Gin(rho * np.exp(1j * (np.angle(L) + np.radians(t)))) for t in tg])
    return np.abs(s), lam_of_tau(tg, p_lokal)


U = np.linspace(np.log(1e-16), np.log(1e14), 80001)
SP = np.exp(U)
DS = SP * (U[1] - U[0])


def omega_neg(tz, sabs, lam, T2):
    a = np.where(SP > np.sqrt(T2), T2 / SP, SP)
    la = np.interp(np.log(a), np.log(sabs), lam, left=0.0)
    P = (tz.imag / np.pi) / ((-SP - tz.real) ** 2 + tz.imag ** 2)
    return np.sum(la * P * DS)


for p_lokal in [1.27, 1.40, 1.13]:
    for rahmen in ["W", "K"]:
        sabs, lam = setup(rahmen, p_lokal)
        T0 = T0W if rahmen == "W" else T0K
        T2 = abs(T0) ** 2
        res = []
        for u in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]:
            z = L * (1 + u)
            tz = tauW(z) if rahmen == "W" else K.tau(z)
            gz = np.log(abs((tz - np.conj(T0)) / (tz - T0)))
            res.append(omega_neg(tz, sabs, lam, T2) / gz)
        print(f"lokaler Exponent {p_lokal:.2f}, Rahmen {'W̄ ' if rahmen == 'W' else 'K′'}: sup-Kandidaten ω(λ)/g am Schlitz "
              f"u = 1e-2/1e-4/1e-6/1e-8/1e-10: " + " / ".join(f"{r:.4f}" for r in res))
