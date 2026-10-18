"""A6: Wo kippt d_n vom Sattel-Regime (positiv, Sattel auf der positiven Achse) ins Regime der Singularität bei -3 (Vorzeichenwechsel)?
Sattel: S_n = g(R) R^{-n}/sqrt(2 pi b);  -3-Anteil: H_n ~ (-1/3)^n /(n log n) * (1+(g0-gamma)/log n)  (g0 = log sexp'(-1)).
"""
import json, mpmath as mp
mp.mp.dps = 60
C = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))][:150]
def sexp(x):
    k = int(mp.floor(x + mp.mpf(0.5))); v = mp.polyval(C[::-1], x - k)
    for _ in range(k): v = mp.exp(v)
    return v
def lg(t): return mp.log(sexp(mp.exp(t)) - mp.log(mp.exp(t) + 2))
g0 = mp.log(C[1])
t = mp.log(2.5)
for n in [140, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000]:
    t = mp.findroot(lambda s: mp.diff(lg, s) - n, t)
    b = mp.diff(lg, t, 2); R = mp.exp(t)
    logS = lg(t) - n * t - mp.log(2 * mp.pi * b) / 2
    logH = -n * mp.log(3) - mp.log(n * mp.log(n)) + mp.log(1 + (g0 - mp.euler) / mp.log(n))
    print(f"n={n:5d} R={mp.nstr(R,7)}  log10 S_n={mp.nstr(logS/mp.log(10),7)}  log10 |H_n|={mp.nstr(logH/mp.log(10),7)}")
