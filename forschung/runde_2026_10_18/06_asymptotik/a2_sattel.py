"""A2: Sattelpunkt-(Hayman-)Näherung für d_n = c_n - [z^n]log(z+2).
d_n ≈ g(R) R^{-n} / sqrt(2 pi b(R)),  a(R) = R g'(R)/g(R) = n,  b = R a'(R);  g = sexp - log(z+2) reell auf R>0.
sexp(x) für x in [0,4] über die Taylorreihe um 0 (|x-2|..., hier: sexp(x)=exp^k(sexp(x-k)), x-k in [-0.5,0.5]).
"""
import json, mpmath as mp
mp.mp.dps = 60
C = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))][:150]
def sexp0(t):
    return mp.polyval(C[::-1], t)
def sexp(x):
    k = int(mp.floor(x + mp.mpf(0.5)))
    v = sexp0(x - k)
    for _ in range(k): v = mp.exp(v)
    return v
def g(x): return sexp(x) - mp.log(x + 2)
def lg(t):   # log g(e^t)
    return mp.log(g(mp.exp(t)))
B = [mp.mpf(s) for s in json.load(open("daten/taylor_0_1.5_400_50.json"))]
def d(n): return B[n] - (-1)**(n+1)/(n*mp.mpf(2)**n)
for n in [20, 40, 60, 80, 100, 120, 140]:
    t = mp.findroot(lambda t: mp.diff(lg, t) - n, mp.log(2.5))
    b = mp.diff(lg, t, 2)
    R = mp.exp(t)
    S = mp.exp(lg(t) - n * t) / mp.sqrt(2 * mp.pi * b)
    print(n, "R=", mp.nstr(R, 8), " Sattel/d_n =", mp.nstr(S / d(n), 8))
