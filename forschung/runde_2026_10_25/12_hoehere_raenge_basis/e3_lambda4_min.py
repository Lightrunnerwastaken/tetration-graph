"""Multiplikator lambda_4(b) = T_b'(p_4(b)) am negativen Fixpunkt der Tetration (er steuert die Pentation:
H_5(b,x) - p_4 ~ c * lambda_4^x fuer x -> -inf). Raster und Minimum per Goldenem Schnitt.
Aufruf: python e3_lambda4_min.py <dps>"""
import sys, json, time
import mpmath as mp
import motor

dps = int(sys.argv[1]); mp.mp.dps = dps + 10


def lam4(b):
    bs = motor.gpzahl(b, dps + 8)
    T = lambda x: mp.re(motor.werte(bs, ["hyper(4,%s)" % motor.gpzahl(x, dps + 8)], dps)[0])
    p = mp.findroot(lambda x: T(x) - x, (mp.mpf("-1.5"), mp.mpf("-1.8")), solver="secant", tol=mp.mpf(10) ** (-2 * dps + 8), verify=False)
    h = mp.mpf(10) ** (-dps // 3)
    l = (T(p + h) - T(p - h)) / (2 * h)
    motor.schliessen()
    return p, l


t0 = time.time()
raster = {}
for k in range(0, 9):
    b = mp.mpf("1.45") + mp.mpf("0.05") * k
    p, l = lam4(b)
    raster[mp.nstr(b, 6)] = [mp.nstr(p, 18), mp.nstr(l, 16)]
    print(mp.nstr(b, 6), mp.nstr(p, 15), mp.nstr(l, 15), f"{time.time()-t0:.0f}s", flush=True)
# Goldener Schnitt
g = (mp.sqrt(5) - 1) / 2
a, c = mp.mpf("1.5"), mp.mpf("1.75")
x1, x2 = c - g * (c - a), a + g * (c - a)
f1, f2 = lam4(x1)[1], lam4(x2)[1]
for it in range(28):
    if f1 < f2:
        c, x2, f2 = x2, x1, f1; x1 = c - g * (c - a); f1 = lam4(x1)[1]
    else:
        a, x1, f1 = x1, x2, f2; x2 = a + g * (c - a); f2 = lam4(x2)[1]
bmin = (a + c) / 2
pm, lm = lam4(bmin)
print("Minimum: b* =", mp.nstr(bmin, 12), "+-", mp.nstr(c - a, 2), " lambda_4 =", mp.nstr(lm, 16), " p_4 =", mp.nstr(pm, 15), f"{time.time()-t0:.0f}s")
(motor.HIER / "daten" / f"lambda4_min_{dps}.json").write_text(json.dumps(
    {"raster": raster, "bstern": mp.nstr(bmin, 14), "breite": mp.nstr(c - a, 3), "lambda4_min": mp.nstr(lm, 18), "p4": mp.nstr(pm, 18)}, indent=1), encoding="utf-8")
