"""Strang V — v59_eval: adjungierte Gewichte w_n(t) = |Λ_n(t)/Λ_n(L)| auf der ganzen Kontur γ (Ks K-Form-Linse über Γin, Γout-Bogen),
aus den exakten Vektoren ℓ_n von v59_elim.py. t_j = k(e^{iψ_j}), k = ρ·w·e^{F(w)} (Theodorsen wie v12), j < MB.
Auswertung Λ_n(t) = Σ ℓ_k (ρ/t)^k in Ballarithmetik (Radius = Genauigkeitsnachweis).
Ausgabe v59_w_<MODELL>_N<N>_S<S>.npz (psi, t, lw_n) und Ratenprofil (max über die Linse ausserhalb ±δ um L, L̄; Steg; Γout).
Aufruf: python v59_eval.py MODELL N S [MB] [tau1] [p]
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, acb, acb_poly, ctx, fmpz

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
MODELL, N, S_ = sys.argv[1], int(sys.argv[2]), sys.argv[3]
MB = int(sys.argv[4]) if len(sys.argv) > 4 else 16384
T1_ = sys.argv[5] if len(sys.argv) > 5 else "12"
P_ = sys.argv[6] if len(sys.argv) > 6 else "1.3"
SUF = ("" if (T1_, P_) == ("12", "1.3") else f"_t{T1_}_p{P_}") + ("" if MB == 16384 else f"_MB{MB}")
import os
NUR = [int(x) for x in os.environ["V59_NUR"].split(",")] if os.environ.get("V59_NUR") else None
_argv = sys.argv
sys.argv = ["v12", "kubisch", "10", S_, T1_, P_, "K"]
sys.path.insert(0, str(HIER))
import v12_konform_S as V12  # noqa: E402
sys.argv = _argv


def lese(me):
    return arb(fmpz(int(me[0]))) * arb(2) ** int(me[1])


def main():
    t0 = time.time()
    D = json.loads((HIER / f"v59_ell_{MODELL}_N{N}.json").read_text(encoding="utf-8"))
    P = D["P"]
    ctx.prec = P + 64
    F, it, dd = V12.theodorsen()
    res, kmin, kmax, mono, dL = V12.pruefe_abbildung(F)
    m_ = len(F)
    big = 4 * m_
    coef = np.zeros(big, dtype=complex)
    coef[:m_] = F
    vals = (np.fft.ifft(coef) * big)[::big // MB]
    psi = 2 * np.pi * np.arange(MB) / MB
    u = np.exp(-vals - 1j * psi)
    t = V12.rho_f * np.exp(vals + 1j * psi)
    thL = V12.thL
    ang = np.mod(np.angle(t), 2 * np.pi)
    gin = (ang > thL) & (ang < 2 * np.pi - thL)
    dist = np.degrees(np.minimum(np.abs(np.angle(t / V12.Lf)), np.abs(np.angle(t / np.conj(V12.Lf)))))
    print(f"[{MODELL}, N = {N}, Linse S = {S_}, MB = {MB}] Theodorsen {it} It., Randfehler {res:.1e}, |k|/ρ ∈ [{kmin:.4f}, {kmax:.4f}], "
          f"F₀ = {F[0].real:.5f}, Γin-Punkte {gin.sum()}", flush=True)
    out = dict(psi=psi, t=t, F=F[:2 * N + 8])
    ns = sorted(int(n) for n in D["ell"] if NUR is None or int(n) in NUR)
    for n in ns:
        ctx.prec = P + 64
        ell = [lese(me) for me in D["ell"][str(n)]]
        LamL = acb(lese(D["LamL"][str(n)][0]), lese(D["LamL"][str(n)][1]))
        bits = max(200, int(math.log2(max(D["maxell"][str(n)], 1.0) / float(abs(LamL).mid()))) + 300 + n // 2)   # + n/2: acb-Horner bläht Radien um ≤ √2 je Schritt
        ctx.prec = min(P, bits)
        poly = acb_poly([acb(x) for x in ell])
        aL = abs(LamL)
        lw = np.empty(MB)
        relmax = 0.0
        for j in range(MB):
            v = poly(acb(float(u[j].real), float(u[j].imag)))
            av = abs(v)
            if av.mid() > 0:
                relmax = max(relmax, float((v.rad() / av).mid()))
                lw[j] = float((av / aL).log().mid())
            else:
                relmax = float("inf")
                lw[j] = -1e300
        out[f"lw_{n}"] = lw
        jl = np.where(gin & (dist > 3))[0]
        jm = jl[np.argmax(lw[jl])]
        print(f"  n = {n} ({ctx.prec} Bit, rel. Radius ≤ {relmax:.1e}): ln w auf Γin (±3° um L frei) max {lw[jm]:+.3f} bei arg t = {math.degrees(ang[jm]):.1f}°, "
              f"Steg {lw[np.argmin(np.abs(ang - np.pi))]:+.3f}, Γout max {lw[~gin & (dist > 3)].max():+.3f} ({time.time() - t0:.0f}s)", flush=True)
    np.savez(HIER / f"v59_w_{MODELL}_N{N}_S{S_}{SUF}.npz", **out)
    print("Ratenprofil (ln w_{n2} − ln w_{n1})/(n2 − n1):", flush=True)
    for n1, n2 in [(600, 1000), (1000, 1500), (1500, 2000), (1600, 2000), (1800, 2000)]:
        if f"lw_{n1}" in out and f"lw_{n2}" in out:
            r = (out[f"lw_{n2}"] - out[f"lw_{n1}"]) / (n2 - n1)
            zeile = []
            for dmin in (1, 3, 10, 30):
                jl = np.where(gin & (dist > dmin))[0]
                jm = jl[np.argmax(r[jl])]
                zeile.append(f"±{dmin}° frei: max {r[jm]:+.5f} bei {math.degrees(ang[jm]):.1f}°")
            js = np.argmin(np.abs(ang - np.pi))
            jo = np.where(~gin & (dist > 3))[0]
            print(f"  {n1}→{n2}: Linse " + " | ".join(zeile) + f" | Steg {r[js]:+.5f} | Γout max {r[jo].max():+.5f}", flush=True)
    print(f"fertig ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
