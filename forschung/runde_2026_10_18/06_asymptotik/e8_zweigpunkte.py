"""E8: Lösungen von sexp(w) = 2 pi i k (k != 0) im Hauptblatt -> logarithmische Verzweigungspunkte w-1 auf Nebenblättern.
Newton mit Motor; Ableitung zentral mit h = 1e-12 (dps 40, Fehler ~h^2).
Startwerte: Gitter in Re in [-1.5,3], Im in [0.2,3.5].
"""
import mpmath as mp, itertools, json
from motor import sexp
DPS = 40
mp.mp.dps = DPS + 10
def F(w): return sexp([w], DPS)[0]
def newton(w, c, it=40):
    for _ in range(it):
        h = mp.mpf("1e-12")
        a, b, f0 = sexp([w + h, w - h, w], DPS)
        d = (a - b) / (2 * h)
        st = (f0 - c) / d
        if abs(st) > 0.5: st = st / abs(st) * 0.5
        w = w - st
        if abs(st) < mp.mpf(10) ** (-30): return w, abs(f0 - c)
    return None, None
gef = []
for c in [2j * mp.pi, 4j * mp.pi]:
    for x, y in itertools.product([-1.5, -0.5, 0.5, 1.5, 2.5], [0.3, 1.0, 2.0, 3.0]):
        try:
            w, res = newton(mp.mpc(x, y), c)
        except Exception as e:
            continue
        if w is None or abs(mp.re(w)) > 20: continue
        if not any(abs(w - u[1]) < 1e-15 for u in gef if u[0] == c):
            gef.append((c, w))
            print("sexp(w) =", mp.nstr(c, 6), " w =", mp.nstr(w, 25), " Zweigpunkt w-1 =", mp.nstr(w - 1, 12), flush=True)
json.dump([[mp.nstr(mp.im(c), 20), mp.nstr(mp.re(w), 35), mp.nstr(mp.im(w), 35)] for c, w in gef], open("daten/zweigpunkte_2pik.json", "w"))
