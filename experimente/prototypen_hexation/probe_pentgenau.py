"""Wie genau sind die Eingangsdaten des Fork-Verfahrens? sexp-Motor (sexpinit), Pentation (pentinit-Schätzung) und pent/ipent an
 komplexen Punkten auf dem Bogen um Re p — Vergleich dps D gegen dps D + 15. Aufruf: probe_pentgenau.py basis D Re_p Im_p
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, D = sys.argv[1], int(sys.argv[2])
c, r = mp.mpf(sys.argv[3]), mp.mpf(sys.argv[4])
mp.mp.dps = D + 30
pts = [c + mp.mpf("0.65") * r * mp.expjpi(mp.mpf(k) / 8) for k in range(1, 8)] + [c + mp.mpf("0.65") * r * mp.expjpi(mp.mpf(1) / 2) * f for f in (mp.mpf("0.5"), mp.mpf("1.3"))]
ausdr = ["sexpinit(sexp(1))[1]", "pentinit(sexp(1),1)"]
for z in pts:
    zs = f"({mp.nstr(mp.re(z), D + 20)}+I*{mp.nstr(mp.im(z), D + 20)})"
    ausdr += [f"iferr(pent({zs}),E,0)", f"iferr(ipent({zs}),E,0)"]
erg = {}
for d in (D, D + 15):
    t0 = time.time()
    gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=d)
    try:
        erg[d] = gp.eval_batch(basis, ausdr)
    finally:
        gp.close()
    print(f"dps {d}: {time.time() - t0:.1f} s; sexp-Motor {mp.nstr(erg[d][0], 6)} Stellen, Pentation {mp.nstr(erg[d][1], 6)} Stellen", flush=True)
a, b = erg[D], erg[D + 15]
for i, z in enumerate(pts):
    dp = abs(a[2 + 2 * i] - b[2 + 2 * i]) / abs(b[2 + 2 * i]) if b[2 + 2 * i] != 0 else mp.nan
    di = abs(a[3 + 2 * i] - b[3 + 2 * i]) / abs(b[3 + 2 * i]) if b[3 + 2 * i] != 0 else mp.nan
    print(f"  z = {mp.nstr(z, 8):28s} pent: rel. Abw. {mp.nstr(dp, 3):10s} ipent: rel. Abw. {mp.nstr(di, 3)}")
