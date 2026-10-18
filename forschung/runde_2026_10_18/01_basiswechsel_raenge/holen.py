"""Holt Kneser-Werte T_b(x) an Tschebyschow-Knoten von [-1, 0] aus dem Motor (eigener Zwischenspeicher).

Aufruf: python holen.py <basis> [N] [dps] [fn]     z. B.  python holen.py 2 120 60 hexx   (fn = sexp | hexx)
Ergebnis: cache/cheb_<basis>_<N>_<dps>.json  (Knoten x_j = -1/2 + cos(pi (j+1/2)/N)/2, Werte T_b(x_j))
"""
import json
import sys
import time

import mpmath as mp

import motor

basis = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 120
dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
fn = sys.argv[4] if len(sys.argv) > 4 else "sexp"
mp.mp.dps = dps + 20
xs = [mp.mpf(-1) / 2 + mp.cos(mp.pi * (j + mp.mpf(1) / 2) / N) / 2 for j in range(N)]
t = time.time()
def ausdruck(arg):
    return f"hyper({fn[5:]}, {arg})" if fn.startswith("hyper") else f"{fn}({arg})"
ws = motor.werte(basis, [ausdruck(motor.gpzahl(x, dps + 15)) for x in xs], dps)
# Kontrollen: T(0) = 1, T(-1) = 0 direkt aus dem Motor
k = motor.werte(basis, [ausdruck("0"), ausdruck("-1"), ausdruck("1/2")], dps)
print(f"{basis}: {N} Werte in {time.time()-t:.1f}s; T(0)-1 = {mp.nstr(k[0]-1, 3)}, T(-1) = {mp.nstr(k[1], 3)}, max|Im| = "
      f"{mp.nstr(max(abs(mp.im(w)) for w in ws), 3)}")
tag = basis.replace("(", "").replace(")", "").replace("/", "_")
out = motor.CDIR / (f"cheb_{tag}_{N}_{dps}.json" if fn == "sexp" else f"{fn}_{tag}_{N}_{dps}.json")
out.write_text(json.dumps({"basis": basis, "N": N, "dps": dps, "x": [mp.nstr(x, dps + 15) for x in xs],
                           "T": [mp.nstr(mp.re(w), dps + 10) for w in ws],
                           "T_half": mp.nstr(mp.re(k[2]), dps + 10)}), encoding="utf-8")
motor.schliessen()
