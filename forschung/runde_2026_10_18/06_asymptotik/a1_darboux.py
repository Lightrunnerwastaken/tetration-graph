"""A1: Taylor-Koeffizienten c_n um 0 gegen die Darboux-Formel (-1)^{n+1}/(n 2^n) aus der Singularität bei -2.
d_n = c_n - [z^n] log(z+2). Satz: d_n = O((3-eps)^{-n}); Leitform d_n ~ (-1/3)^n /(n log n) (Singularität log log bei -3).
"""
import json
import mpmath as mp
mp.mp.dps = 80
A = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1_512_60.json"))]
B = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))]
def d(c, n):
    return c[n] - (mp.log(2) if n == 0 else (-1) ** (n + 1) / (n * mp.mpf(2) ** n))
print(" n   c_n(A)                    rel.Abw. Darboux    d_n*(-3)^n*n*log n   |dA-dB|/|d|")
for n in [1,2,3,5,10,15,20,30,40,50,60,70,80,90,100,110,120,140]:
    if n >= len(B): break
    dA, dB = d(A, n), d(B, n)
    dar = (-1) ** (n + 1) / (n * mp.mpf(2) ** n)
    print(f"{n:3d} {mp.nstr(A[n],18):>26} {mp.nstr((A[n]-dar)/dar,6):>14} {mp.nstr(dB*(-3)**n*n*mp.log(n),10):>16} {mp.nstr(abs(dA-dB)/abs(dB),3):>10}")
