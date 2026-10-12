"""Strang V — v59_gewichtet: gewichtetes CF-Problem (Auftrag Leitung 13:00)
  CF_w(N) = min_{P̃ ≡ p_N mod x^{N+2}} sup_γ |Λ_N(t)/Λ_N(L)|·|P̃(t) − f_K(t)|.
Mit k: 𝔻 → Int γ (v12) und der äusseren Funktion Φ auf 𝔻, |Φ| = w_N∘k auf ∂𝔻 (Φ = e^h, Re h = ln w_N, Im h = konjugierte Funktion; Φ und 1/Φ
analytisch ⇒ Φ·w^{N+2}H^∞ = w^{N+2}H^∞, also bleibt die Nebenklasse erhalten):
  CF_w(N) = ‖T_{N+2}((Φ·c) mod w^{N+2})‖₂ (Carathéodory–Fejér/Schur), c_m = [w^m](p_N − f_K)∘k aus v59c.
Eingaben: v59c_S<S>_N<N>.npy, v59_w_<MODELL>_N2001_S<S>.npz (ln w_n auf MB Randpunkten). Prüfungen: Re h = ln w, Aliasing von Φ,
Halbierung von MB (jeder zweite Punkt), Gewicht ≡ 1 ⇒ CF wie v57.
Aufruf: python v59_gewichtet.py MODELL S1,S2,… [N1,N2,…]
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
MODELL = sys.argv[1]
SL = sys.argv[2].split(",")              # S oder S_t<τ₁>_p<p>
NS = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [200, 400, 600, 800, 1000, 1500, 2000]
ORDNER = HIER / "cfw"


def toepnorm(c):
    M = len(c)
    ij = np.subtract.outer(np.arange(M), np.arange(M))
    T = np.where(ij >= 0, c[np.maximum(ij, 0)], 0)
    return float(np.linalg.norm(T, 2))


def faltung(a, b):
    """Faltung; mit V59G_EXAKT=1 exakt aus den Double-Werten (flint, 300 Bit), sonst numpy."""
    if not os.environ.get("V59G_EXAKT"):
        return np.convolve(a, b)
    from flint import acb, acb_poly, ctx
    ctx.prec = 300
    pa = acb_poly([acb(complex(x).real, complex(x).imag) for x in a])
    pb = acb_poly([acb(complex(x).real, complex(x).imag) for x in b])
    pc = pa * pb
    return np.array([complex(float(pc[i].real.mid()), float(pc[i].imag.mid())) for i in range(pc.degree() + 1)])


def aussen(lw):
    MB = len(lw)
    a = np.fft.fft(lw) / MB
    h = np.zeros(MB, dtype=complex)
    h[0] = a[0]
    h[1:MB // 2] = 2 * a[1:MB // 2]
    h[MB // 2] = a[MB // 2]
    hv = np.fft.ifft(h) * MB
    fehler_re = float(np.max(np.abs(hv.real - lw)))
    Phv = np.exp(hv)
    Phk = np.fft.fft(Phv) / MB
    alias = float(np.max(np.abs(Phk[MB // 4: 3 * MB // 4])) / np.max(np.abs(Phk)))
    return Phk, fehler_re, alias


def main():
    erg = {}
    for S in SL:
        W = np.load(ORDNER / f"v59_w_{MODELL}_N2001_S{S}{os.environ.get('V59G_SUF', '')}.npz")
        for N in NS:
            q = ORDNER / f"v59c_S{S}_N{N}.npy"
            if not q.exists() or f"lw_{N}" not in W.files:
                continue
            c = np.load(q)
            M = N + 2
            lw = W[f"lw_{N}"]
            Phk, fre, alias = aussen(lw)
            g = faltung(Phk[:M], c)[:M]
            cfw = toepnorm(g)
            Phk2, _, _ = aussen(lw[::2])
            cfw2 = toepnorm(faltung(Phk2[:M], c)[:M])
            cf = toepnorm(c)
            ga = np.abs(g)
            erg[f"{S}|{N}"] = dict(S=S, N=N, CF=cf, CFw=cfw, CFw_halbe=cfw2, re_fehler=fre, alias=alias,
                                   lw_max=float(lw.max()), max_g=float(ga.max()), m_max=int(np.argmax(ga)))
            print(f"S = {S}, N = {N}: CF = {cf:.3e}, CF_w = {cfw:.3e} (MB/2: {cfw2:.3e}; Re h − ln w ≤ {fre:.0e}, Aliasing {alias:.0e}); "
                  f"max ln w = {lw.max():+.2f}; max|(Φc)_m| = {ga.max():.2e} bei m = {int(np.argmax(ga))}", flush=True)
        ks = sorted(int(k.split("|")[1]) for k in erg if k.split("|")[0] == S)
        for n1, n2 in zip(ks[:-1], ks[1:]):
            e1, e2 = erg[f"{S}|{n1}"], erg[f"{S}|{n2}"]
            print(f"  S = {S}, {n1}→{n2}: Rate CF {math.log(e2['CF'] / e1['CF']) / (n2 - n1):+.5f}, Rate CF_w {math.log(e2['CFw'] / e1['CFw']) / (n2 - n1):+.5f}, "
                  f"Potenz CF_w {math.log(e2['CFw'] / e1['CFw']) / math.log(n2 / n1):+.2f}", flush=True)
    json.dump(erg, open(ORDNER / f"v59_cfw_{MODELL}{os.environ.get('V59G_SUF', '')}{'_exakt' if os.environ.get('V59G_EXAKT') else ''}_{'_'.join(SL)}.json", "w", encoding="utf-8"), indent=0)


if __name__ == "__main__":
    main()
