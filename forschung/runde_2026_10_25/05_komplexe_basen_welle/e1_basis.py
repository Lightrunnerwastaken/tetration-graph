"""E1: Grunddaten einer komplexen Basis: Fixpunkte L± (Newton ab sexp(±4i)), Multiplikatoren, Funktionalgleichung, A = sexp'(0).
Aufruf: python e1_basis.py "<gp-ausdruck der basis>" <dps>
Ausgabe: daten/basis_<tag>_<dps>.json
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
from kern import Fix, sexp, schliessen, _tag, HIER

basis, DPS = sys.argv[1], int(sys.argv[2])
mp.mp.dps = DPS + 20
t0 = time.time()
# Basis als Zahl: aus dem Motor selbst (sexp(1) = b)
pts = [mp.mpc(0), mp.mpc(-1), mp.mpc(1), mp.mpc("0.5"), mp.mpc("-0.5"), mp.mpc(0, 4), mp.mpc(0, -4), mp.mpc("0.3", "0.7"), mp.mpc("1.3", "0.7"),
       mp.mpc("0.3", "-0.7"), mp.mpc("1.3", "-0.7")]
v = sexp(basis, pts, DPS)
b = v[2]
lb = mp.log(b)
print(f"[{time.time()-t0:.0f}s] b = {mp.nstr(b, 20)}, ln b = {mp.nstr(lb, 20)}")
print("sexp(0) =", mp.nstr(v[0], 12), " sexp(-1) =", mp.nstr(v[1], 12))
print("FG: |sexp(1.3+0.7i) - b^sexp(0.3+0.7i)| =", mp.nstr(abs(v[8] - mp.exp(lb * v[7])), 5),
      " unten:", mp.nstr(abs(v[10] - mp.exp(lb * v[9])), 5))
fp = Fix(b, v[5])
fm = Fix(b, v[6])
for nm, f, w in (("L+", fp, v[5]), ("L-", fm, v[6])):
    print(f"{nm} = {mp.nstr(f.L, 25)}  |sexp(±4i) - L| = {mp.nstr(abs(w - f.L), 3)}  λ = {mp.nstr(f.lam, 15)}  |λ| = {mp.nstr(abs(f.lam), 12)}  k = {f.k}")
# A = sexp'(0), dazu sexp''(0)/2 per Cauchy
K = 32
rho = mp.mpf("0.1")
cp = [rho * mp.expjpi(2 * mp.mpf(j) / K) for j in range(K)]
cv = sexp(basis, cp, DPS)
co = [mp.fsum(cv[j] * mp.expjpi(-2 * mp.mpf(j * n) / K) for j in range(K)) / K / rho ** n for n in range(4)]
print("Taylor sexp bei 0:", [mp.nstr(c, 20) for c in co])
out = {"basis": basis, "dps": DPS, "b": [mp.nstr(mp.re(b), DPS), mp.nstr(mp.im(b), DPS)],
       "Lp": [mp.nstr(mp.re(fp.L), DPS), mp.nstr(mp.im(fp.L), DPS)], "Lm": [mp.nstr(mp.re(fm.L), DPS), mp.nstr(mp.im(fm.L), DPS)],
       "lam_p": mp.nstr(fp.lam, 20), "lam_m": mp.nstr(fm.lam, 20), "kp": fp.k, "km": fm.k,
       "taylor0": [[mp.nstr(mp.re(c), DPS), mp.nstr(mp.im(c), DPS)] for c in co]}
(HIER / "daten").mkdir(exist_ok=True)
fn = HIER / "daten" / f"basis_{_tag(basis)}_{DPS}.json"
fn.write_text(json.dumps(out, indent=0), encoding="utf-8")
print(f"[{time.time()-t0:.0f}s] geschrieben {fn.name}")
schliessen()
