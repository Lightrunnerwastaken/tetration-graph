"""E4: reguläre (Koenigs-)Tetration knapp innerhalb des Rands am goldenen Randpunkt.

t = (1-eps) e^{2 pi i alpha} (anziehend). kappa_reg = Koenigs-Koordinate von v(1) = a - t:
kappa = lim t^{-k} g^k(a - t),  g(v) = t(e^v - 1).  Dann a(T_reg(z) - L) = v_t(kappa t^z).
Gegenprobe: Motor (rechnet im Inneren regulär) bei z = 0.5+3i.
Aufruf: python e4_regulaer_innen.py eps [motor]
"""
import sys
import mpmath as mp

eps = mp.mpf(sys.argv[1])
mp.mp.dps = 45
alpha = (3 - mp.sqrt(5)) / 2
t = (1 - eps) * mp.exp(2j * mp.pi * alpha)
a = t * mp.exp(-t)
L = mp.exp(t)
NC = 40
c = [mp.mpc(0)] * (NC + 1); e = [mp.mpc(0)] * (NC + 1)
c[1] = e[0] = e[1] = mp.mpc(1)
for n in range(2, NC + 1):
    R = mp.fsum(k * c[k] * e[n - k] for k in range(1, n)) / n
    c[n] = R / (t ** (n - 1) - 1)
    e[n] = c[n] + R
v = lambda w: mp.polyval(c[::-1], w)
dv = lambda w: mp.polyval([k * c[k] for k in range(NC, 0, -1)], w)


def phi(p):
    k = 0
    while abs(p) > 1e-6:
        p = t * (mp.exp(p) - 1)
        k += 1
        if k > 10 ** 6:
            raise RuntimeError
    w = p
    for _ in range(30):
        w = w - (v(w) - p) / dv(w)
    return w / t ** k, k


kap, k = phi(a - t)
print(f"eps={eps}  |t|={mp.nstr(abs(t), 8)}  Iterationen {k}  kappa_reg = {mp.nstr(kap, 12)}  |kappa_reg| = {mp.nstr(abs(kap), 12)}")
z = mp.mpc(0.5, 3)
Tz = L + v(kap * t ** z) / a
print("  T_reg(0.5+3i) =", mp.nstr(Tz, 15))
if len(sys.argv) > 2:
    from motor import werte, schliessen, gpzahl
    b = mp.exp(a)
    Tm = werte(gpzahl(b, 40), [f"sexp({gpzahl(z, 35)})"], 30)[0]
    schliessen()
    print("  Motor        =", mp.nstr(Tm, 15), "  Abw.", mp.nstr(abs(Tm - Tz), 3))
