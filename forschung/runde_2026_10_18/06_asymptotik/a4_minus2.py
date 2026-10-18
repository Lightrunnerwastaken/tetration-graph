"""A4: Vollständige Entwicklung bei -2: sexp(-2+e) = log e + sum_{n>=0} g_n e^n, g = log(sexp(-1+e)/e), Radius genau 1.
Taylor um -1 aus zwei Einstellungen (r=0.6,N=400,dps60 / r=0.8,N=300,dps45).
"""
import json, mpmath as mp
from motor import sexp
mp.mp.dps = 70
def lade(f): return [mp.mpf(s) for s in json.load(open(f))]
TA = lade("daten/taylor_-1_0.6_400_60.json"); TB = lade("daten/taylor_-1_0.8_300_45.json")
def gko(T, K):
    a = [T[k + 1] for k in range(K + 1)]            # sexp(-1+e)/e
    a0 = a[0]; q = [x / a0 for x in a]; q[0] = mp.mpf(0)
    # log(1+Q): l' = Q'/(1+Q)
    l = [mp.log(a0)] + [mp.mpf(0)] * K
    for m in range(1, K + 1):
        l[m] = q[m] - mp.fsum(j * l[j] * q[m - j] for j in range(1, m)) / m
    return l
K = 140
gA, gB = gko(TA, K), gko(TB, K)
print("g_0 = log sexp'(-1) =", mp.nstr(gA[0], 30), "  sexp'(-1) =", mp.nstr(TA[1], 30))
for n in [1, 2, 3, 4, 5, 6, 8, 10, 20, 40, 60, 80, 100, 120, 140]:
    print(f"{n:4d} g_n={mp.nstr(gA[n],22):>28}  |A-B|={mp.nstr(abs(gA[n]-gB[n]),3):>9}  (-1)^n g_n n log n = {mp.nstr((-1)**n*gA[n]*n*mp.log(n) if n>1 else 0,8)}")
json.dump([mp.nstr(x, 40) for x in gA], open("daten/g_minus2.json", "w"))
# Gegenprobe gegen Motor an Punkten links von -2 (oberes Ufer bzw. komplex)
mp.mp.dps = 50
for e in [mp.mpf("0.5"), mp.mpc(-0.5, 0.3), mp.mpc(-0.8, 0.05), mp.mpc(0.2, -0.7)]:
    z = -2 + e
    ser = mp.log(e) + mp.polyval(gA[:K][::-1], e)
    mot = sexp([z], 45)[0]
    print("z =", mp.nstr(z, 6), " Reihe-Motor =", mp.nstr(abs(ser - mot), 3), " Motor =", mp.nstr(mot, 15))
