"""O39 — Das feste Grenzobjekt auf der Anpassungsskala z = (t − L)/(ρ·N^{−1/Re κ₁}): Ψ_Y(z) = Y_1(t)/Y_1(L), Ψ_T(z) = T_1(t)/T_1(L).

Je N: (a) Werte von Y_1, T_1, 𝒰 (Erzeugende der Vorwärts-Eckmode u_1: 𝒰(η) = Σ u_k ρ^k η^{−k−1}) auf Kreisen |z| = A (5°-Gitter, A bis 8);
(b) Taylor-Koeffizienten von Ψ_Y, Ψ_T, 𝒰 in z um 0 (FFT auf |z| = 1, 64 Punkte);
(c) nur kubisch: Probe der exakten Eckformel D6 (W §17): R(ζ) := T(ζ) − Y(ζ)m(ζ)/ζ − Y(h(ζ))/(h f'(h)) + 𝒰(ζ) = ferne Urbilder + niedrige Koeff.
Speichert o39_<modell>_<N>.npz.  Aufruf: python o39_zskala.py N1,N2,... [modell]
"""
import os
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(v) for v in sys.argv[1].split(",")]
name = sys.argv[2] if len(sys.argv) > 2 else "kubisch"
L = K.Lc
RHO = K.RHO
argL = np.angle(L)
REK = 1.1182
mu = np.sqrt(L)
AS = [0.5, 1, 2, 3, 4, 6, 8]
PS = np.deg2rad(np.arange(-180, 180, 5))
NF = 64

for N in NS:
    t0 = time.time()
    n = N + 1
    J, M, komplex = K.matrizen(name, N, int(1.3 * N) + 300)
    I = arb_mat(n, n)
    for i in range(n):
        I[i, i] = arb(1)
    Jinv = K.to_np(J.solve(I, algorithm="approx"), komplex)
    U, s, Vt = np.linalg.svd(Jinv)
    u1 = U[:, 0]
    ctx.prec = 200 + N // 4
    ctx.cap = n + 10                  # Reihenlänge (sonst kappt flint bei 10 Termen!)
    fk, mk, _ = K.modell(name, N)
    fpol = acb_poly([acb(c) for c in fk])
    mpol = acb_poly([acb(c) for c in mk])
    rp = RHO ** np.arange(n)
    kk = np.arange(n)
    skala = RHO * N ** (-1 / REK)

    def werte(t):
        ev = (RHO / t) ** kk
        g = Jinv @ ev
        et = acb(t.real, t.imag)
        s1 = list(acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs())
        s1 += [acb(0)] * (n - len(s1))
        kz = acb_poly([et ** (-(j + 1)) for j in range(n)])
        s2 = list((mpol * kz).coeffs())[:n]
        s2 += [acb(0)] * (n - len(s2))
        c = np.array([complex(float(s1[k].real.mid()) + float(s2[k].real.mid()),
                              float(s1[k].imag.mid()) + float(s2[k].imag.mid())) for k in range(n)]) * rp
        d = t ** (-(kk + 1.0)) * rp
        gt = Jinv @ c - d
        return u1 @ g, u1 @ gt, (u1 @ ev) / t

    yL, tL, uL = werte(L)
    kreis = {}
    for A_ in AS:
        kreis[A_] = np.array([werte(L + A_ * skala * np.exp(1j * (argL + p))) for p in PS])
    zf = np.exp(2j * np.pi * np.arange(NF) / NF)
    fw = np.array([werte(L + skala * np.exp(1j * argL) * z) for z in zf])    # z gemessen in Einheiten skala·e^{i arg L}
    tay = np.fft.fft(fw, axis=0) / NF                                       # Taylor-Koeff. in z (Radius 1)
    d6 = None
    if name == "kubisch":
        fn = np.polynomial.Polynomial([float(c.mid()) for c in fk[:4]])
        fpn = fn.deriv()
        mn = lambda x: fpn(fn(x))
        d6 = []
        for A_ in [1, 2, 4]:
            for p in PS[::6]:
                zeta = L + A_ * skala * np.exp(1j * (argL + p))
                x = L + (zeta - L) / mu
                for _ in range(40):
                    x = x - (fn(x) - zeta) / fpn(x)
                yz, tz, uz = werte(zeta)
                yh, _, _ = werte(x)
                R = tz - yz * mn(zeta) / zeta - yh / (x * fpn(x)) + uz
                d6.append([A_, p, abs(R), abs(yz), abs(tz), abs(uz), abs(yh)])
        d6 = np.array(d6)
    np.savez(f"o39_{name}_{N}.npz", N=N, s=s[:8], yL=yL, tL=tL, uL=uL, AS=np.array(AS), PS=PS,
             kreis=np.array([kreis[a] for a in AS]), tay=tay, d6=d6 if d6 is not None else np.zeros(0))
    print(f"{name} N={N} ({time.time() - t0:.0f}s): ln|Y(L)| {np.log(abs(yL)):.3f}  ln|T(L)| {np.log(abs(tL)):.3f}  ln|𝒰(L)| {np.log(abs(uL)):.3f}  "
          f"(ln √N = {0.5 * np.log(N):.3f})", flush=True)
    ty = tay[:, 0] / yL
    tt = tay[:, 1] / tL
    print("   Taylor Ψ_Y: " + " ".join(f"{ty[j].real:+.3f}{ty[j].imag:+.3f}i" for j in range(6)), flush=True)
    print("   Taylor Ψ_T: " + " ".join(f"{tt[j].real:+.3f}{tt[j].imag:+.3f}i" for j in range(6)), flush=True)
    if d6 is not None:
        print(f"   D6-Rest |R|/|Y|: max {np.max(d6[:, 2] / d6[:, 3]):.2e}, |R| max {np.max(d6[:, 2]):.3e}; |𝒰|/|Y| auf A=1: "
              f"{np.median(d6[d6[:, 0] == 1, 5] / d6[d6[:, 0] == 1, 3]):.3f}", flush=True)
