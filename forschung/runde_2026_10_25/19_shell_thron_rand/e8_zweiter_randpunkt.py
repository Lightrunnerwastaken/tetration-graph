"""E8: Gegenprobe am zweiten Randpunkt alpha = sqrt2 - 1 (bzw. beliebiges alpha als Argument).
(i) konforme Radius r aus Koeffizienten-Hülle, (ii) |phi(-lambda)| aus Birkhoff-Mittel,
(iii) Motor-Lösung (Kneser-Art) auf dem Rand eps=0 und knapp aussen eps=0.002: kappa_K,
(iv) reguläre Koenigs-Konstante knapp innen eps=-0.002.
Aufruf: python e8_zweiter_randpunkt.py [sqrt2|gold]
"""
import sys
import numpy as np
import mpmath as mp
from e1_siegel import koeff_np, radius, birkhoff
from motor import werte, schliessen, gpzahl

name = sys.argv[1] if len(sys.argv) > 1 else "sqrt2"
mp.mp.dps = 45
alpha = mp.sqrt(2) - 1 if name == "sqrt2" else (3 - mp.sqrt(5)) / 2
lam0 = mp.exp(2j * mp.pi * alpha)
S = 1.7
c = koeff_np(complex(lam0), 6000, S)
est, rfit, pw = radius(c, 20)
print(f"alpha={mp.nstr(alpha,12)}  r(Fenster) = {[round(S*e,6) for e in est]}  r(Hülle) = {S*rfit:.6f}  Potenz {pw:.3f}")
w0, mx = birkhoff(lam0, -lam0, 8000, 25)
print(f"Birkhoff N=8000: phi(-lam) = {mp.nstr(w0,10)}  |.| = {mp.nstr(abs(w0),8)}   lam*phi = {mp.nstr(lam0*w0,8)}")


def poincare(t, NC=40):
    cc = [mp.mpc(0)] * (NC + 1); e = [mp.mpc(0)] * (NC + 1)
    cc[1] = e[0] = e[1] = mp.mpc(1)
    for n in range(2, NC + 1):
        R = mp.fsum(k * cc[k] * e[n - k] for k in range(1, n)) / n
        cc[n] = R / (t ** (n - 1) - 1)
        e[n] = cc[n] + R
    return (lambda w: mp.polyval(cc[::-1], w)), (lambda w: mp.polyval([k * cc[k] for k in range(NC, 0, -1)], w))


for eps in [mp.mpf(0), mp.mpf("0.002")]:
    t = (1 + eps) * lam0
    a = t * mp.exp(-t); b = mp.exp(a); L = mp.exp(t)
    z = mp.mpc(0.5, 3)
    T3, Tl = werte(gpzahl(b, 40), [f"sexp({gpzahl(z, 35)})", "sexp(0.5-10*I)"], 30)
    schliessen()
    v, dv = poincare(t)
    target = a * (T3 - L); w = target
    for _ in range(40):
        w = w - (v(w) - target) / dv(w)
    kap = w / t ** z
    print(f"Motor eps={mp.nstr(eps,3)}: kappa_K = {mp.nstr(kap,10)}  |kappa_K| = {mp.nstr(abs(kap),10)}  T(0.5-10i) = {mp.nstr(Tl,8)}")

eps = mp.mpf("-0.002")
t = (1 + eps) * lam0
a = t * mp.exp(-t)
v, dv = poincare(t)
p = a - t; k = 0
while abs(p) > 1e-6:
    p = t * (mp.exp(p) - 1); k += 1
w = p
for _ in range(30):
    w = w - (v(w) - p) / dv(w)
kap = w / t ** k
print(f"regulär eps=-0.002: kappa_reg = {mp.nstr(kap,10)}  |kappa_reg| = {mp.nstr(abs(kap),10)}")
