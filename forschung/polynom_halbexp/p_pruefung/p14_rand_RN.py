"""Strang P, Prüfung 14: Pol-Antwort R_N(η) = |y_η|_{ℓ¹_ρ} (y_η = A_f⁻¹B_N[k_η]_{>N} = Π_N k_η − [k_η]_{≤N}) auf dem Rand ∂K' =
Γin ∪ f_K(Γout) (obere Hälfte; untere per Konjugation), für die einteilige Bedingung im K'-Rahmen (Θ'_N, Prüfung 15).
A_f = J_N(f_K) aus den 4000 Kneser-Koeffizienten (≈ 50 Stellen), python-flint, Präzision 320 Bit.
Aufruf: python p14_rand_RN.py 100,200,300   → p14_RN_N{N}.json
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, acb_mat, acb_poly, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
ctx.prec = 320
mp.mp.dps = 60
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
A_mp = [mp.mpf(s) for s in A_str]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(Lm)
argL = mp.arg(Lm)


def fK_mp(z):
    s = mp.mpc(0)
    for c in reversed(A_mp):
        s = s * z + c
    return s


def to_acb(z):
    return acb(arb(mp.nstr(mp.re(z), 55)), arb(mp.nstr(mp.im(z), 55)))


def randpunkte():
    """(Name, Parameter, η) — Γin: τ hinter L (Grad), Γout: t vor L (Grad); dicht nahe L."""
    taus = sorted(set(np.round(np.concatenate([np.logspace(-3, np.log10(103.38), 60), [103.38]]), 6)))
    ts = sorted(set(np.round(np.concatenate([np.logspace(-3, np.log10(76.618), 50), [76.618]]), 6)))
    pts = [("L", 0.0, Lm)]
    for t in taus:
        pts.append(("in", float(t), rho * mp.expj(argL + mp.radians(t))))
    for t in ts:
        pts.append(("out", float(t), fK_mp(rho * mp.expj(argL - mp.radians(t)))))
    return pts


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n, P):
    """Newton: 1/a mod x^n (P = arb_poly oder acb_poly)."""
    g = P([1 / a[0]])
    m = 1
    while m < n:
        m = min(2 * m, n)
        g = mullow(g, P([2]) - mullow(a.truncate(m), g, m), m)
    return g


Ns = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [100]
PTS = randpunkte()
for N in Ns:
    t0 = time.time()
    ctx.prec = max(320, 2 * N + 300)
    f = arb_poly([arb(s) for s in A_str[:N + 1]])
    fp = f.derivative()
    ex = arb_poly([1 / arb.fac_ui(kk) for kk in range(N + 1)])
    m = mullow(ex, inv_series(fp, N + 1, arb_poly), N + 1)               # m = eˣ/f' = f'∘f
    J = arb_mat(N + 1, N + 1)
    pw = arb_poly([1])
    for j in range(N + 1):
        for kk in range(N + 1):
            J[kk, j] = pw[kk] + (m[kk - j] if kk >= j else 0)
        pw = mullow(pw, f, N + 1)
    Ji = acb_mat(J.inv())
    rho_a = arb(mp.nstr(rho, 55))
    wts = [rho_a ** kk for kk in range(N + 1)]
    print(f"N = {N}: Matrix + Inverse {time.time() - t0:.0f} s", flush=True)
    out = []
    fa = acb_poly([acb(c) for c in f.coeffs()])
    ma = acb_poly([acb(c) for c in m.coeffs()])
    for name, par, eta in PTS:
        e = to_acb(eta)
        inv_f = inv_series(acb_poly([e]) - fa, N + 1, acb_poly)
        kz = acb_poly([1 / e ** (kk + 1) for kk in range(N + 1)])
        w = inv_f + mullow(ma, kz, N + 1)
        wv = acb_mat([[w[kk]] for kk in range(N + 1)])
        y = Ji * wv
        R = arb(0)
        for kk in range(N + 1):
            R += abs(y[kk, 0] - kz[kk]) * wts[kk]
        out.append({"rand": name, "par": par, "eta": [float(mp.re(eta)), float(mp.im(eta))],
                    "R": float(R.mid()), "rad": float(R.rad())})
    (HIER / f"p14_RN_N{N}.json").write_text(json.dumps({"N": N, "punkte": out}), encoding="utf-8")
    RL = out[0]["R"]
    print(f"   R_N(L) = {RL:.4f}; Γin τ = 0.01/1/10/103°: " +
          ", ".join(f"{o['R']:.3f}" for o in out if o["rand"] == "in" and any(abs(o["par"] - v) / v < 0.05 for v in [0.01, 1, 10, 103.38])) +
          "; Γout t = 0.01/1/10/76.6°: " +
          ", ".join(f"{o['R']:.3e}" for o in out if o["rand"] == "out" and any(abs(o["par"] - v) / v < 0.05 for v in [0.01, 1, 10, 76.618])) +
          f"  (max Radius {max(o['rad'] for o in out):.1e}; {time.time() - t0:.0f} s)", flush=True)
