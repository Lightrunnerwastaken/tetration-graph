"""Scharfe Ungleichungen für Knesers sexp (Basis e oder 2) mit Ballarithmetik.

Baut auf w_reell.py (Kopie aus forschung/identitaeten_2026_10_14/luecken/wendepunkt) auf: gleiche Annahme H
(Motorwerte an den Kreisknoten auf EPS genau, |sexp| <= C auf den Kreisscheiben). Alles andere streng.

Geprüft wird:
  (W)  Wendepunkt x_w (Neuberechnung, w_reell.zertifiziere).
  (M)  m = T'(x_w) = min T' auf (-2, oo);  d = T'(0) = max T' auf [-1, 0].
  (TH) theta(x) = slog_lin(T(x)) - x: Extrema auf [-1, 0] (Nullstellen von T' - 1), scharfe Schranken.
  (E)  Schwelle y* von E(y) = (T(y) + T(y+1))/2 - T(y+1/2) auf [x_w - 1, x_w].
  (P)  TPID 7: T'(-1+y) < T'(-1-y) für y in [x_w+1, 1)  =>  T(y)T(-y) streng fallend auf [0, 1).
  (LK) (log T')'' > 0 auf [x_w, x_w + 1], Minimum mu.
Aufruf: python u_zert.py e   (bzw. 2)
"""
import json
import sys
import time
from fractions import Fraction as F

from flint import arb

import w_reell
from w_reell import KONF, Sexp, intervall, arb_bruch, _frac_unten, _frac_oben, HIER

HALB = F(1, 2)


def stuecke(A, B, n_pro_einheit, schnitte_extra=()):
    """Zerlegt [A, B] in Stücke, die keine Vielfachen von 1/2 (und keine Extra-Schnitte) überqueren."""
    s = {A, B}
    k = int((A - 1) * 2) - 1
    while F(k, 2) <= B + 1:
        if A < F(k, 2) < B:
            s.add(F(k, 2))
        k += 1
    for e in schnitte_extra:
        if A < e < B:
            s.add(e)
    s = sorted(s)
    out = []
    for lo, hi in zip(s[:-1], s[1:]):
        n = max(1, int((hi - lo) * n_pro_einheit))
        out += [(lo + (hi - lo) * i / n, lo + (hi - lo) * (i + 1) / n) for i in range(n)]
    return out


def halbteilen(lo, hi, extra):
    """Mittelpunkt, ggf. so verschoben, dass keine halbganze Stelle (auch nach Verschiebung um extra) überquert wird."""
    return (lo + hi) / 2


class Werkzeug:
    def __init__(self, basis):
        self.basis = basis
        self.cfg = KONF[basis]
        self.S = Sexp(self.cfg)
        self.lam = self.cfg["lnb"]()

    def W(self, lo, hi=None):
        hi = lo if hi is None else hi
        return self.S.werte(lo, hi)

    def sicher(self, f, lo, hi, tiefe=0, maxt=14):
        """f(lo, hi) -> Ball; gibt Liste (lo, hi, vorzeichen) mit adaptiver Halbierung zurück."""
        try:
            b = f(lo, hi)
            vz = 1 if b > 0 else (-1 if b < 0 else 0)
        except RuntimeError:
            vz = 0
        if vz == 0 and tiefe < maxt:
            m = (lo + hi) / 2
            return self.sicher(f, lo, m, tiefe + 1, maxt) + self.sicher(f, m, hi, tiefe + 1, maxt)
        return [(lo, hi, vz)]

    def newton(self, f, df, lo, hi, it=80):
        """Intervall-Newton: f, df: (lo, hi) -> Ball. Gibt (lo, hi) mit genau einer Nullstelle zurück."""
        d0 = df(lo, hi)
        if not (d0 > 0 or d0 < 0):
            raise RuntimeError("Newton: Ableitung nicht vorzeichenfest")
        for _ in range(it):
            m = (lo + hi) / 2
            N = arb_bruch(m) - f(m, m) / df(lo, hi)
            nl, nh = max(lo, _frac_unten(N)), min(hi, _frac_oben(N))
            if nl > nh:
                raise RuntimeError("Newton: leer")
            fertig = (nh - nl) >= (hi - lo) * F(99, 100)
            lo, hi = nl, nh
            if fertig:
                break
        return lo, hi


def vz_verlauf(erg):
    """Fasst Liste (lo, hi, vz) zu Blöcken gleichen Vorzeichens zusammen."""
    bl = []
    for lo, hi, v in erg:
        if bl and bl[-1][2] == v and bl[-1][1] == lo:
            bl[-1] = (bl[-1][0], hi, v)
        else:
            bl.append((lo, hi, v))
    return [(float(a), float(b), v) for a, b, v in bl]


def main(basis):
    t0 = time.time()
    out = {"basis": basis}
    print(f"=== Basis {basis} ===", flush=True)
    wz = Werkzeug(basis)
    W = wz.W
    lam = wz.lam
    # (W) Wendepunkt neu
    zw = w_reell.zertifiziere(basis)
    Jlo, Jhi = F(zw["J"][0]), F(zw["J"][1])
    xw_lo, xw_hi = wz.newton(lambda a, b: W(a, b)[2], lambda a, b: W(a, b)[3], Jlo, Jhi)
    xw = intervall(xw_lo, xw_hi)
    out["x_w"] = xw.str(30)
    print("x_w =", xw.str(30), flush=True)

    # (M)
    m = W(xw_lo, xw_hi)[1]
    d = W(F(0))[1]
    Tw = W(xw_lo, xw_hi)[0]
    out["m=T'(x_w)"] = m.str(30)
    out["d=T'(0)"] = d.str(30)
    out["T(x_w)"] = Tw.str(30)
    out["T(x_w+1)"] = (lam * Tw).exp().str(30)
    print("m = T'(x_w) =", m.str(30), "  d = T'(0) =", d.str(30), flush=True)

    # (TH) theta = T(x) - x - 1 auf [-1, 0]; kritische Punkte T' = 1
    gitter = [F(-1) + F(k, 200) for k in range(201)]
    g = [W(x)[1] - 1 for x in gitter]
    krit = []
    for i in range(200):
        if (g[i] > 0 and g[i + 1] < 0) or (g[i] < 0 and g[i + 1] > 0):
            a, b = gitter[i], gitter[i + 1]
            if a < HALB * -1 < b:
                raise RuntimeError("Gitter überquert -1/2")
            lo, hi = wz.newton(lambda u, v: W(u, v)[1] - 1, lambda u, v: W(u, v)[2], a, b)
            th = W(lo, hi)[0] - intervall(lo, hi) - 1
            krit.append(dict(x=intervall(lo, hi).str(30), theta=th.str(30)))
            print("  T'=1 bei", intervall(lo, hi).str(25), " theta =", th.str(25), flush=True)
    out["theta_krit"] = krit
    # Nullstellen von theta im Inneren von (-1, 0): Schnittpunkte Kneser / lineare Tetration
    th_g = [W(x)[0] - arb_bruch(x) - 1 for x in gitter]
    nullen = []
    for i in range(1, 199):
        if (th_g[i] > 0 and th_g[i + 1] < 0) or (th_g[i] < 0 and th_g[i + 1] > 0):
            a, b = gitter[i], gitter[i + 1]
            lo, hi = wz.newton(lambda u, v: W(u, v)[0] - intervall(u, v) - 1, lambda u, v: W(u, v)[1] - 1, a, b)
            nullen.append(intervall(lo, hi).str(30))
            print("  theta = 0 bei", intervall(lo, hi).str(25), flush=True)
    out["theta_null_innen"] = nullen
    # Vergleich mit der C^2-kubischen Polynomtetration (nur Basis e): p(x) = 1 + a1 x + a1^2/4 x^2 + a1^2/6 x^3
    if basis == "e":
        a1 = 6 - 2 * arb(6).sqrt()
        p = lambda x: 1 + a1 * x + a1 ** 2 / 4 * x ** 2 + a1 ** 2 / 6 * x ** 3
        dp = lambda x: a1 + a1 ** 2 / 2 * x + a1 ** 2 / 2 * x ** 2
        def pinv(w):
            x = arb(w.mid()) - 1
            for _ in range(60):
                x = x - (p(x) - w) / dp(x)
                x = arb(x.mid())
            return x
        thc = []
        for k in range(1, 400):
            x = F(-1) + F(k, 400)
            thc.append((float((pinv(W(x)[0]) - arb_bruch(x)).mid()), float(x)))
        out["theta_kubisch_max"] = max(thc)
        out["theta_kubisch_min"] = min(thc)
        out["a1_kubisch"] = a1.str(20)
        print("  kubisch C^2: theta in", min(thc), max(thc), flush=True)
    # Produktwerte T(y)T(-y) (Kontrolle)
    out["T(y)T(-y)"] = {str(F(k, 10)): (W(F(k, 10))[0] * W(-F(k, 10))[0]).str(15) for k in range(1, 10)}
    # Schwellen y_t* von E_t(y) = (1-t)T(y) + tT(y+1) - T(y+t) (Punktwerte, Bisektion auf [x_w-1, x_w])
    def Et(t, y):
        ta = arb_bruch(t)
        return float(((1 - ta) * W(y)[0] + ta * W(y + 1)[0] - W(y + t)[0]).mid())
    ys = {}
    for k in range(1, 10):
        t = F(k, 10)
        lo, hi = xw_lo - 1, xw_hi
        for _ in range(45):
            mm = (lo + hi) / 2
            mm = F(round(float(mm) * 10 ** 13), 10 ** 13)
            if Et(t, mm) < 0:
                lo = mm
            else:
                hi = mm
        ys[str(t)] = float(lo)
    out["y_t*"] = ys
    print("  y_t*:", ys, flush=True)
    # Vollständigkeit: T' - 1 hat keine weiteren Nullstellen (T' streng monoton auf [-1, x_w] und [x_w, 0])
    if not (xw_hi < 0 and xw_lo > -1):
        print("  Achtung: x_w nicht in (-1, 0)")
    out["theta_vollstaendig"] = "T' streng fallend auf [-1,x_w], streng wachsend auf [x_w,0] (Zertifikat W) => hoechstens je eine Nullstelle von T'-1"

    # (E) Schwelle y* von E(y) = (T(y)+T(y+1))/2 - T(y+1/2) auf [x_w-1, x_w]
    def E(a, b):
        return (W(a, b)[0] + W(a + 1, b + 1)[0]) / 2 - W(a + HALB, b + HALB)[0]

    def dE(a, b):
        return (W(a, b)[1] + W(a + 1, b + 1)[1]) / 2 - W(a + HALB, b + HALB)[1]

    A, B = xw_lo - 1, xw_hi
    erg = []
    for lo, hi in stuecke(A, B, 400):
        erg += wz.sicher(E, lo, hi)
    bl = vz_verlauf(erg)
    print("  E-Vorzeichenbloecke:", bl, flush=True)
    # Nullstelle eingrenzen: genau ein Vorzeichenwechsel -1 -> (0 ...) -> +1, auf der Klammer dE > 0
    res_E = dict(bloecke=bl)
    i0 = max(i for i, e in enumerate(erg) if e[2] == -1)
    i1 = min(i for i, e in enumerate(erg) if e[2] == 1)
    if all(e[2] == -1 for e in erg[:i0 + 1]) and all(e[2] == 1 for e in erg[i1:]) and all(e[2] == 0 for e in erg[i0 + 1:i1]):
        a = erg[i0][0]; b = erg[i1][1]
        if dE(a, b) > 0:
            lo, hi = wz.newton(E, dE, a, b)
            res_E["y*"] = intervall(lo, hi).str(30)
            res_E["x*=y*+1"] = (intervall(lo, hi) + 1).str(30)
            res_E["dE>0 auf Klammer"] = True
            print("  y* =", intervall(lo, hi).str(30), flush=True)
    out["E_halb"] = res_E

    # (P) TPID 7: T'(-1+y) < T'(-1-y) auf [x_w + 1, 1 - delta]; Rest analytisch
    delta = F(1, 100)
    def P(a, b):   # obere Schranke von T'(-1+y) - T'(-1-y) auf [a, b]
        return W(-1 + a, -1 + b)[1] - W(-1 - b, -1 - a)[1]
    erg = []
    for lo, hi in stuecke(xw_lo + 1, 1 - delta, 200):
        erg += wz.sicher(P, lo, hi, maxt=8)
    out["P_bloecke"] = vz_verlauf(erg)
    ok_P = all(e[2] == -1 for e in erg)
    rest = W(-2 + delta)[1] - W(F(0))[1]
    out["P_ok"] = bool(ok_P)
    out["P_rest: T'(-2+delta) - T'(0)"] = rest.str(10)
    print("  TPID7-Ableitungsbedingung:", ok_P, "Rest", rest.str(10), flush=True)
    # Kontrollwerte des Produkts
    out["T(1/2)T(-1/2)"] = (W(HALB)[0] * W(-HALB)[0]).str(25)

    # (LK) (log T')'' = (T''' T' - T''^2)/T'^2 > 0 auf [x_w, x_w + 1]
    def LK(a, b):
        w = W(a, b)
        return (w[3] * w[1] - w[2] ** 2) / w[1] ** 2
    erg = []
    lows = []
    for lo, hi in stuecke(xw_lo, xw_hi + 1, 500):
        r = wz.sicher(LK, lo, hi, maxt=6)
        erg += r
        for a, b, v in r:
            lows.append((LK(a, b).lower(), a, b))
    ok_LK = all(e[2] == 1 for e in erg)
    unten = min(lows, key=lambda t: t[0])
    # Punktwert-Minimum auf feinem Gitter um die Stelle
    a0 = unten[1]
    pts = [a0 - F(1, 100) + F(k, 10000) for k in range(201)]
    pv = [(LK(x, x), x) for x in pts if not any(x == F(j, 2) for j in range(-4, 4))]
    pmin = min(pv, key=lambda t: t[0].mid())
    # Verfeinern: Stücke mit unterer Schranke < pmin - tol halbieren
    tol = arb("1e-7")
    ziel = pmin[0].lower() - tol
    offen = [(a, b) for (lw, a, b) in lows if lw < ziel]
    besser = min(lw for (lw, a, b) in lows if not lw < ziel) if any(not lw < ziel for (lw, a, b) in lows) else None
    runden = 0
    while offen and runden < 30:
        runden += 1
        neu = []
        for a, b in offen:
            mm = (a + b) / 2
            for u, v in ((a, mm), (mm, b)):
                lw = LK(u, v).lower()
                if lw < ziel:
                    neu.append((u, v))
                else:
                    besser = lw if besser is None or lw < besser else besser
        offen = neu
    if not offen:
        unten = (besser, None, None)
    out["LK_verfeinert_runden"] = runden
    out["LK_ok"] = bool(ok_LK)
    out["mu_untere_Schranke"] = unten[0].str(12)
    out["mu_obere_Schranke(Punktwert)"] = pmin[0].str(15)
    out["mu_Stelle"] = float(pmin[1])
    print("  (log T')''>0:", ok_LK, " mu in [", unten[0].str(10), ",", pmin[0].str(12), "] bei", float(pmin[1]), flush=True)

    out["sekunden"] = round(time.time() - t0, 1)
    (HIER / f"ergebnis_{basis}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print("fertig in", out["sekunden"], "s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "e")
