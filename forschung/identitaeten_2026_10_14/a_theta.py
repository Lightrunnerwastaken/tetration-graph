"""Etappe 1, Teil A3: Fourier-Koeffizienten von Kneisers θ zur Basis e.

  sexp_e(z) = F(z + θ(z)),  F(z) = σ⁻¹(L^z) := σ⁻¹(e^{Lz}),  σ Schröder-Funktion von exp bei L (σ(L) = 0, σ'(L) = 1),
  θ(z) = a_0 + Σ_{n≥1} a_n e^{2πinz}  im oberen Halbraum (a_0 nur modulo 2πi/L festgelegt).
  θ(z) = log σ(sexp z)/L − z,  σ(w) = L^m s(u_m), u_{m+1} = log(1 + u_m/L), u_0 = w − L, s = Schröder-Reihe bei 0.
Abtastung auf Im z = y (Motor erlaubt |Im z| ≤ 0.2), zwei Höhen y und zwei Knotenzahlen N als Abgleich:
a_n = e^{2πny}·(1/N) Σ_j θ(t_j + iy) e^{−2πin t_j} muss von y unabhängig sein.
Aufruf: python a_theta.py <dps>
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

from motor import gpzahl, werte

DPS = int(sys.argv[1])
mp.mp.dps = DPS + 30
HIER = Path(__file__).resolve().parent
L = mp.conj(-mp.lambertw(-1))   # Fixpunkt e^L = L mit Im L > 0 (−W₀(−1) hat Im < 0)
t0 = time.time()


def schroeder_reihe(K):
    """s(L(e^u − 1)) = L s(u), s(u) = u + Σ_{k≥2} s_k u^k (Koeffizienten numerisch)."""
    # g(u) = L(e^u − 1) als Reihe; Potenzen von g
    g = [mp.mpf(0)] + [L / mp.factorial(k) for k in range(1, K + 1)]

    def mul(a, b):
        c = [mp.mpc(0)] * (K + 1)
        for i in range(K + 1):
            if a[i] != 0:
                for j in range(K + 1 - i):
                    c[i + j] += a[i] * b[j]
        return c
    gp = [None, g]
    for k in range(2, K + 1):
        gp.append(mul(gp[-1], g))
    s = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
    for m in range(2, K + 1):
        # [u^m]: Σ_k s_k gp[k][m] = L s_m ;  gp[m][m] = L^m
        rest = mp.fsum(s[k] * gp[k][m] for k in range(1, m))
        s[m] = rest / (L - L ** m)
    return s


S = schroeder_reihe(30)
MSTEP = 130


def theta_roh(w, z):
    for _ in range(60):                        # erst Hauptzweig-Logarithmus (stetig in w), bis nahe L
        w = mp.log(w)
    u = w - L
    assert abs(u) < mp.mpf("1e-4"), mp.nstr(u, 5)
    for _ in range(MSTEP - 60):
        u = mp.log(1 + u / L)
    su = mp.polyval(S[::-1], u)
    return MSTEP + mp.log(su) / L - z          # log L = L (Hauptzweig)


def sexp(zs):
    return werte("exp(1)", [f"sexp({gpzahl(z, DPS + 10)})" for z in zs], DPS)


def koeffizienten(y, N, nmax=64):
    ts = [mp.mpf(j) / N - mp.mpf(1) / 2 for j in range(N)]      # t ∈ [−½, ½)
    zs = [t + 1j * y for t in ts]
    w = sexp(zs)
    th = [theta_roh(wi, zi) for wi, zi in zip(w, zs)]
    per = 2j * mp.pi / L
    # stetig machen (Zweig von log)
    for j in range(1, N):
        best = min(range(-3, 4), key=lambda m: abs(th[j] + m * per - th[j - 1]))
        th[j] += best * per
    best = min(range(-3, 4), key=lambda m: abs(th[0] + m * per - th[-1]))
    sprung = abs(th[0] - th[-1])
    print("   max Schritt", mp.nstr(max(abs(th[j] - th[j - 1]) for j in range(1, N)), 3), "Randsprung", mp.nstr(sprung, 3), "Zweig am Rand", best, flush=True)
    a = []
    for n in range(nmax + 1):
        s = mp.fsum(th[j] * mp.exp(-2j * mp.pi * n * zs[j]) for j in range(N)) / N
        a.append(s)
    return a, th


ergebnis = {"L": mp.nstr(L, DPS)}
res = {}
for (y, N) in ((mp.mpf("0.2"), 384), (mp.mpf("0.2"), 512), (mp.mpf("0.18"), 512)):
    a, th = koeffizienten(y, N)
    res[(str(y), N)] = a
    print(f"[{time.time()-t0:.0f}s] y={y} N={N}", flush=True)
    for n in range(6):
        print(f"   a_{n} = {mp.nstr(a[n], 30)}", flush=True)
(a1, a2, a3) = res.values()
for n in range(65):
    ergebnis[f"a_{n}"] = [mp.nstr(mp.re(a2[n]), DPS + 5), mp.nstr(mp.im(a2[n]), DPS + 5)]
    d_N = abs(a1[n] - a2[n])
    d_y = abs(a3[n] - a2[n]) if n > 0 else min(abs(a3[0] - a2[0] + m * 2j * mp.pi / L) for m in range(-3, 4))
    ergebnis[f"abgleich_a_{n}"] = {"N": mp.nstr(d_N, 3), "y": mp.nstr(d_y, 3)}
    print(f"a_{n}: |ΔN| = {mp.nstr(d_N, 3)}  |Δy| = {mp.nstr(d_y, 3)}", flush=True)
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / f"theta_{DPS}.json").write_text(json.dumps(ergebnis, indent=1), encoding="utf-8")
