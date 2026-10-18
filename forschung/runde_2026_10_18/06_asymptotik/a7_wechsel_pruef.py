"""A7: d_n aus zwei Einstellungen; Vergleich mit Sattel S_n und -3-Anteil H_n; Vorzeichenwechsel."""
import json, mpmath as mp
mp.mp.dps = 60
P = [mp.mpf(s) for s in json.load(open("daten/gkoeff_1.9_1400_140.json"))]
Q = [mp.mpf(s) for s in json.load(open("daten/gkoeff_1.8_1300_150.json"))]
C = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))][:150]
def sexp(x):
    k = int(mp.floor(x + mp.mpf(0.5))); v = mp.polyval(C[::-1], x - k)
    for _ in range(k): v = mp.exp(v)
    return v
def lg(t): return mp.log(sexp(mp.exp(t)) - mp.log(mp.exp(t) + 2))
g0 = mp.log(C[1]); t = mp.log(2.5)
print("n   d_n(P)              rel|P-Q|     S_n/d_n    (d_n-S_n)/H_n")
for n in [100,200,300,400,450,480,500,520,540,560,600,640]:
    t = mp.findroot(lambda s: mp.diff(lg, s) - n, t)
    b = mp.diff(lg, t, 2)
    S = mp.exp(lg(t) - n * t) / mp.sqrt(2 * mp.pi * b)
    H = (-mp.mpf(1)/3)**n / (n*mp.log(n)) * (1 + (g0-mp.euler)/mp.log(n))
    print(n, mp.nstr(P[n],12), mp.nstr(abs(P[n]-Q[n])/abs(P[n]),3), mp.nstr(S/P[n],8), mp.nstr((P[n]-S)/H,6))
s=''.join('+' if P[n]>0 else '-' for n in range(400,650))
print(s)
print("erster Index ab dem d_n alterniert:", next(n for n in range(300,650) if all(mp.sign(P[m])==(-1)**m*mp.sign(P[640]) for m in range(n,650))))
