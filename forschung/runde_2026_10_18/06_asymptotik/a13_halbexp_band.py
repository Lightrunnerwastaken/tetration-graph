"""A13: Halbe e-Funktion f = sexp(slog x + 1/2) im Band der linearen Höhe:
|slog_lin(f(x)) - slog_lin(x) - 1/2| <= A = max P - min P (A11). Stichprobe x in (0, 1e6]."""
import json, mpmath as mp
mp.mp.dps = 40
C = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))][:180]
S = lambda u: mp.polyval(C[::-1], u)
def sexp(x):
    k = int(mp.floor(x)) + 1; v = S(x - k)      # x-k in [-1,0)
    for _ in range(k): v = mp.exp(v)
    return v
def slog(y):
    n = 0
    while y > 1: y = mp.log(y); n += 1
    return n + mp.findroot(lambda u: S(u) - y, y - 1)
def slog_lin(y):
    n = 0
    while y > 1: y = mp.log(y); n += 1
    return n - 1 + y
A = mp.mpf("0.0170358261009277477145")
worst = 0
for i in range(1, 400):
    x = mp.mpf(10) ** (mp.mpf(i) / 66) * mp.mpf("0.05")
    fx = sexp(slog(x) + mp.mpf(1) / 2)
    dev = slog_lin(fx) - slog_lin(x) - mp.mpf(1) / 2
    worst = max(worst, abs(dev))
print("max |Abweichung| =", mp.nstr(worst, 10), " Schranke A =", mp.nstr(A, 10))
