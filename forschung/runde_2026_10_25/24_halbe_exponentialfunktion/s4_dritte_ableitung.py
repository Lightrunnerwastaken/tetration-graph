"""S4: f_K''' > 0 auf ganz R (numerisch, mit strenger Rueckfuehrung auf zwei kompakte Pruefungen, Satz K5).
  (i)  x <= 0:  f(x) = g(e^x), g = log f.  f'''(x) = u*h3(u), u = e^x, h3 = g' + 3u g'' + u^2 g'''.  Pruefe h3 > 0 auf [0, 1].
  (ii) [0, e^{x2}] mit f'(x2) = 2: direkte Gitterpruefung von f'''.
  (iii) [x2, oo): Fortpflanzung  e^{3x} f'''(e^x) = e^f [f''' + 3 f''(f'-1) + f'(f'-1)(f'-2)]  (alle Glieder >= 0 fuer f' >= 2).
Dazu: Nullstellen von f^(4) (Vorzeichenwechsel) und Ableitungen bei x* fuer den Bericht.
"""
import sys
import mpmath as mp
from s2_ableitungen import taylor0, abl

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30

# (i)
mn, wo = None, None
for i in range(0, 2001):
    u = mp.mpf(i) / 2000
    T = taylor0(u)
    f0, f1, f2, f3 = T[0], T[1], 2 * T[2], 6 * T[3]
    g1 = f1 / f0
    g2 = f2 / f0 - g1 ** 2
    g3 = f3 / f0 - 3 * f1 * f2 / f0 ** 2 + 2 * g1 ** 3
    h3 = g1 + 3 * u * g2 + u ** 2 * g3
    if mn is None or h3 < mn:
        mn, wo = h3, u
print("(i)  min h3 auf [0,1] =", mp.nstr(mn, 10), "bei u =", mp.nstr(wo, 5), "  (h3(0) = a1/a0 =", mp.nstr(taylor0(0)[1] / taylor0(0)[0], 10), ")")

# (ii)
x2 = mp.findroot(lambda x: abl(x)[1] - 2, mp.mpf("1.2"))
E = mp.exp(x2)
print("(ii) x2 (f' = 2) =", mp.nstr(x2, 15), "  e^{x2} =", mp.nstr(E, 10))
mn, wo = None, None
N = 1500
for i in range(N + 1):
    x = E * i / N
    v = abl(x)[3]
    if mn is None or v < mn:
        mn, wo = v, x
print("     min f''' auf [0, e^{x2}] (Gitter", N + 1, "Punkte) =", mp.nstr(mn, 10), "bei x =", mp.nstr(wo, 6))

# f^(4): Vorzeichenwechsel
for s in ["-1.05", "0.65"]:
    z = mp.findroot(lambda x: abl(x)[4], mp.mpf(s))
    print("f^(4) = 0 bei x =", mp.nstr(z, 15), "  f(x) =", mp.nstr(abl(z)[0], 12))
xs = mp.findroot(lambda x: abl(x)[1] - 1, mp.mpf("0.24"))
print("bei x* =", mp.nstr(xs, 20), ": f, f', f'', f''' =", [mp.nstr(v, 12) for v in abl(xs)[:4]])
