"""Reelle Fixpunkt-Landkarte der Raenge 4-6 als Funktion der Basis.
p_4(b): negativer Fixpunkt der Tetration (Baustelle der Pentation), lambda_4 = T'(p_4).
Pentation: Minimum von H_5(x)-x auf [-12,0] (kein negativer Fixpunkt erwartet).
Hexation: Singularitaet sigma_6 = ihyper(6,p_4) - 1 (dort H_6 -> -inf), negativer Fixpunkt p_6, lambda_6,
Kontrolle der ganzzahligen Identitaet H_6(-3) = -2.
Aufruf: python e1_struktur.py <dps> <basis> [<basis> ...]"""
import sys, json, time
import mpmath as mp
import motor

dps = int(sys.argv[1]); BASEN = sys.argv[2:]
mp.mp.dps = dps + 10
OUT = motor.HIER / "daten" / "struktur.json"
alles = json.loads(OUT.read_text(encoding="utf-8")) if OUT.is_file() else {}


def H(b, r, x):
    return mp.re(motor.werte(b, ["hyper(%d,%s)" % (r, motor.gpzahl(x, dps + 8))], dps)[0])


def iH(b, r, y):
    return mp.re(motor.werte(b, ["ihyper(%d,%s)" % (r, motor.gpzahl(y, dps + 8))], dps)[0])


def fix(b, r, x0, x1):
    p = mp.findroot(lambda x: H(b, r, x) - x, (mp.mpf(x0), mp.mpf(x1)), solver="secant", tol=mp.mpf(10) ** (-2 * dps + 4))
    h = mp.mpf(10) ** (-dps // 3)
    lam = (H(b, r, p + h) - H(b, r, p - h)) / (2 * h)
    return p, lam


t0 = time.time()
for b in BASEN:
    e = {}
    p4, l4 = fix(b, 4, "-1.5", "-1.8")
    e["p4"], e["lambda4"] = mp.nstr(p4, dps - 4), mp.nstr(l4, 12)
    g = [mp.mpf(k) / 2 for k in range(-24, 1)]
    d5 = [H(b, 5, x) - x for x in g]
    e["min_H5_minus_x_auf_[-12,0]"] = mp.nstr(min(d5), 10)
    e["H5(-12)"] = mp.nstr(H(b, 5, -12), 15)
    print(b, "p4", e["p4"], "l4", e["lambda4"], "min(H5-x)", e["min_H5_minus_x_auf_[-12,0]"], f"{time.time()-t0:.0f}s", flush=True)
    try:
        y6 = mp.findroot(lambda y: H(b, 6, y) - p4, (mp.mpf(-3), mp.mpf(-2)), solver="illinois", tol=mp.mpf(10) ** (-2 * dps + 4), verify=False)
        s6 = y6 - 1
        e["sigma6"] = mp.nstr(s6, dps - 6)
        e["H6(-3)+2"] = mp.nstr(H(b, 6, -3) + 2, 3)
        e["H6(sigma6+1)-p4"] = mp.nstr(H(b, 6, s6 + 1) - p4, 3)
        p6 = mp.findroot(lambda x: H(b, 6, x) - x, (s6 + mp.mpf(10) ** -8, mp.mpf(-3)), solver="illinois", tol=mp.mpf(10) ** (-2 * dps + 4), verify=False)
        h = mp.mpf(10) ** (-dps // 3)
        l6 = (H(b, 6, p6 + h) - H(b, 6, p6 - h)) / (2 * h)
        e["p6"], e["lambda6"] = mp.nstr(p6, dps - 6), mp.nstr(l6, 12)
        e["sigma6<p6<-3"] = bool(s6 < p6 < -3)
        print(b, "sigma6", e["sigma6"], "p6", e["p6"], "l6", e["lambda6"], "H6(-3)+2", e["H6(-3)+2"], f"{time.time()-t0:.0f}s", flush=True)
    except Exception as ex:
        e["hex_fehler"] = str(ex)[:300]
        print(b, "Hexation-Fehler", str(ex)[:300], flush=True)
    alles[b + "|" + str(dps)] = e
    OUT.write_text(json.dumps(alles, indent=1), encoding="utf-8")
    motor.schliessen()
