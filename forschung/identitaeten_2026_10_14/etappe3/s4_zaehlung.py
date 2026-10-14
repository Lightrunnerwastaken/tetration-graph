"""Argumentprinzip: Nullstellen von sexp' und sexp'' im Rechteck [x0, x0+1] × [−Y, Y].

Wegen der Spiegelsymmetrie ist die Windungszahl = Δarg / π längs des oberen Halbwegs
x0+1 → x0+1+iY → x0+iY → x0.  Werte: Taylor-Reihen nahe der Achse, sonst sexp = F(z+θ(z)).
Aufruf: python s4_zaehlung.py
"""
import json

import mpmath as mp

from lib3 import HIER, Basis, ThetaReihe, lade_theta

mp.mp.dps = 30
ZEN = {}
for f in ("taylor_130_256.json", "taylor_links_80_256.json"):
    d = json.loads((HIER / "daten" / f).read_text(encoding="utf-8"))
    for k, v in d["zentren"].items():
        ZEN[mp.mpf(k)] = [mp.mpf(x) for x in v[:110]]
B = Basis(50)
L = B.L
TR = ThetaReihe(lade_theta(), B)


def Psi012(u):
    n = 0
    while abs(u) > mp.mpf("0.02"):
        u /= L; n += 1
    P = B.P
    w = mp.polyval(P[::-1], u)
    d1 = mp.polyval([k * P[k] for k in range(len(P) - 1, 0, -1)], u)
    d2 = mp.polyval([k * (k - 1) * P[k] for k in range(len(P) - 1, 1, -1)], u)
    for _ in range(n):
        ew = mp.exp(w)
        d2 = ew * (d1 ** 2 + d2) / L ** 2
        d1 = ew * d1 / L
        w = ew
    return w, d1, d2


def werte(z):
    """(sexp, sexp', sexp'') für Re z in [−1.9, 0.5], Im z >= 0."""
    c = min(ZEN, key=lambda c: abs(z - c))
    if mp.im(z) <= mp.mpf("0.09") and abs(z - c) <= mp.mpf("0.1"):
        co = ZEN[c]; u = z - c
        s0 = mp.polyval(co[::-1], u)
        s1 = mp.polyval([k * co[k] for k in range(len(co) - 1, 0, -1)], u)
        s2 = mp.polyval([k * (k - 1) * co[k] for k in range(len(co) - 1, 1, -1)], u)
        return s0, s1, s2
    th, t1, t2 = TR.theta(z, 2)
    s = z + th
    u = mp.exp(L * s)
    w, p1, p2 = Psi012(u)
    F1 = p1 * L * u
    F2 = p2 * (L * u) ** 2 + p1 * L ** 2 * u
    return w, F1 * (1 + t1), F2 * (1 + t1) ** 2 + F1 * t2


def pfad(x0, Y):
    """Liste von Funktionen-Auswertern entlang des oberen Halbwegs (Parameter σ in [0,3])."""
    def punkt(sig):
        if sig <= 1:                  # rechte Seite, aufwärts: x0+1+i·sig·Y
            return "r", mp.mpc(x0 + 1, sig * Y)
        if sig <= 2:                  # oben, nach links
            return "o", mp.mpc(x0 + 1 - (sig - 1), Y)
        return "l", mp.mpc(x0, (3 - sig) * Y)
    return punkt


def funktionen(x0, Y):
    punkt = pfad(x0, Y)

    def f(sig):
        art, z = punkt(sig)
        if art == "r":                # über die Funktionalgleichung aus der linken Linie
            s0, s1, s2 = werte(z - 1)
            e = mp.exp(s0)
            d1 = e * s1
            d2 = e * (s1 ** 2 + s2)
            return d1, d2
        if art == "o":
            # oben: Re z kann > 0.5 sein -> θ-Reihe direkt
            s0, s1, s2 = werte(z) if mp.re(z) <= 0.5 else werte(z - 1)
            if mp.re(z) > 0.5:
                e = mp.exp(s0)
                return e * s1, e * (s1 ** 2 + s2)
            return s1, s2
        s0, s1, s2 = werte(z)
        return s1, s2
    return f


def windung(x0, Y=mp.mpf(4), n0=600):
    f = funktionen(x0, Y)
    sig = [mp.mpf(3) * j / n0 for j in range(n0 + 1)]
    vals = [f(s) for s in sig]
    tot1 = mp.mpf(0); tot2 = mp.mpf(0)
    i = 0
    maxsprung = 0
    while i < len(sig) - 1:
        a1 = mp.arg(vals[i + 1][0] / vals[i][0]); a2 = mp.arg(vals[i + 1][1] / vals[i][1])
        if max(abs(a1), abs(a2)) > 0.3 and sig[i + 1] - sig[i] > mp.mpf(1e-9):
            m = (sig[i] + sig[i + 1]) / 2
            sig.insert(i + 1, m); vals.insert(i + 1, f(m))
            continue
        maxsprung = max(maxsprung, abs(a1), abs(a2))
        tot1 += a1; tot2 += a2; i += 1
    return tot1 / mp.pi, tot2 / mp.pi, len(sig), maxsprung


if __name__ == "__main__":
    print("x0      N(sexp')   N(sexp'')   Punkte  max Δarg")
    erg = {}
    for x0 in ("-1.8", "-1.6", "-1.4", "-1.2", "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4"):
        x = mp.mpf(x0)
        n1, n2, npt, ms = windung(x)
        erg[x0] = [mp.nstr(n1, 8), mp.nstr(n2, 8)]
        print(f"{x0:6s} {mp.nstr(n1, 6):>10s} {mp.nstr(n2, 6):>10s} {npt:8d} {mp.nstr(ms, 3)}", flush=True)
    (HIER / "daten" / "zaehlung.json").write_text(json.dumps(erg, indent=1), encoding="utf-8")
