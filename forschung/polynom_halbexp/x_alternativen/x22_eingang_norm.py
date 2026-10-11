"""X22 — Eingangsnorm für den eckfreien Teil (Satz P-II / X1): |w_η|_{ℓ²_ρ}, w_η = T_M[1/(η − f_K) + m·k_η], auf Ps ∂K'-Gitter
(L, Γin-Punkte, f_K(Γout)-Punkte), M = 100 … 800. P-II beweist |w_η| ≤ C·M^{3/2}; Frage: gilt |w_η| ≤ C_w·√M (Faktorisierung
1/(η − f(x)) = 1/((ξ − x)·Q_ξ(x)), ξ = h(η), |ξ| ≥ ρ)? Ausgabe je M: max über η von |w_η|/√(M+1) und von |T_M[1/(η−f)]|, Ort des Maximums.
Aufruf: python x22_eingang_norm.py M1 M2 ...
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, acb_poly, arb, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 60
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
A_mp = [mp.mpf(s) for s in A_str]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho_m = abs(Lm)
argL = mp.arg(Lm)


def fK_mp(z):
    s = mp.mpc(0)
    for c in reversed(A_mp):
        s = s * z + c
    return s


def to_acb(z):
    return acb(arb(mp.nstr(mp.re(z), 55)), arb(mp.nstr(mp.im(z), 55)))


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n, P):
    g = P([1 / a[0]])
    k = 1
    while k < n:
        k = min(2 * k, n)
        g = mullow(g, P([2]) - mullow(a.truncate(k), g, k), k)
    return g


taus = list(np.logspace(-3, np.log10(103.38), 24)) + [103.38]
ts = list(np.logspace(-3, np.log10(76.618), 24)) + [76.618]
PTS = [("L", 0.0, Lm)] + [("in", t, rho_m * mp.expj(argL + mp.radians(t))) for t in taus] + \
      [("out", t, fK_mp(rho_m * mp.expj(argL - mp.radians(t)))) for t in ts]
for M in (int(a) for a in sys.argv[1:]):
    t0 = time.time()
    ctx.prec = int(1.8 * M) + 300
    F = acb_poly([acb(arb(s)) for s in A_str[: M + 1]])
    Fr = arb_poly([arb(s) for s in A_str[: M + 1]])
    ex = arb_poly([1 / arb.fac_ui(k) for k in range(M + 1)])
    m = acb_poly([acb(c) for c in [mullow(ex, inv_series(Fr.derivative(), M + 1, arb_poly), M + 1)[k] for k in range(M + 1)]])
    rho = arb(mp.nstr(rho_m, 55))
    rk = [arb(1)]
    for k in range(1, M + 1):
        rk.append(rk[-1] * rho)
    best_w, best_q = (0, None), (0, None)
    for name, par, eta in PTS:
        e = to_acb(eta)
        q = inv_series(acb_poly([e]) - F, M + 1, acb_poly)
        kz = acb_poly([1 / e ** (k + 1) for k in range(M + 1)])
        w = q + mullow(m, kz, M + 1)
        nq = float(np.sqrt(sum(abs(complex((q[k] * rk[k]).mid())) ** 2 for k in range(M + 1))))
        nw = float(np.sqrt(sum(abs(complex((w[k] * rk[k]).mid())) ** 2 for k in range(M + 1))))
        if nw / np.sqrt(M + 1) > best_w[0]:
            best_w = (nw / np.sqrt(M + 1), (name, par))
        if nq > best_q[0]:
            best_q = (nq, (name, par))
    print(f"M = {M}: max_η |w_η|/√(M+1) = {best_w[0]:.4f} bei {best_w[1]};  max_η |T_M[1/(η−f)]|_ℓ² = {best_q[0]:.4f} bei {best_q[1]}  "
          f"({time.time() - t0:.0f}s)", flush=True)
