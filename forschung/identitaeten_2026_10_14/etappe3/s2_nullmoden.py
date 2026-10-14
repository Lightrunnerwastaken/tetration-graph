"""Nullmoden-Identitäten auf der reellen Achse (Basis e).

θ(x) = h(sexp x) − x auf (0,1) (Randwerte von oben, stetiger Zweig), θ'(x) = sexp'(x) h'(sexp x) − 1.
Geprüft wird, dass der Periodenmittelwert multiplikativ ist (θ ist Randwert einer in der oberen Halbebene
holomorphen 1-periodischen Funktion mit Grenzwert a_0):
  M[θ²] = M[θ]²,  M[θ³] = M[θ]³,  M[e^{cθ}] = e^{c M[θ]},  M[e^{Lθ}]·M[e^{−Lθ}] = 1,
  M[log(1+θ')] = 0,  M[(1+θ')^s] = 1 (Re s < 1),  M[θ e^{+2πimx}] = 0,  M[θ e^{−2πimx}] = a_m.
Quadratur: tanh-sinh auf (0,1) mit Schrittweiten 2^{-k}.
Aufruf: python s2_nullmoden.py DPS  (Taylor-Daten aus daten/taylor_130_256.json)
"""
import json
import sys
import time

import mpmath as mp

from lib3 import HIER, Basis, lade_theta

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 80
mp.mp.dps = dps
T = json.loads((HIER / "daten" / "taylor_130_256.json").read_text(encoding="utf-8"))
ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in T["zentren"].items()}
B = Basis(70)
L = B.L
PER = 2j * mp.pi / L


def sexp_reihe(x):
    """(sexp(x) − 1 bzw. sexp(x), sexp'(x)) für x in [−0.5, 0.5] über die nächste Taylor-Reihe."""
    c = min(ZEN, key=lambda c: abs(x - c))
    co = ZEN[c]
    u = x - c
    if c == 0:
        # sexp − 1 ohne Auslöschung
        s1 = mp.polyval(co[:0:-1], u) * u
    else:
        s1 = mp.polyval(co[::-1], u) - 1
    d = mp.polyval([k * co[k] for k in range(len(co) - 1, 0, -1)], u)
    return s1, d


def theta_und_ableitung(x, d_links, d_rechts):
    """x in (0,1); d_links = x, d_rechts = 1 − x (genau)."""
    if d_links <= mp.mpf("0.5"):
        xs = d_links
    else:
        xs = -d_rechts            # x − 1
    s1, ds = sexp_reihe(xs)
    l1 = mp.log1p(s1)              # log sexp(xs)
    hv, hd = B.h_mit_ableitung(l1)  # h(sexp) = 1 + h(log sexp)
    th = 1 + hv - xs
    w = 1 + s1
    thp = ds * hd / w - 1          # h'(w) = h'(log w)/w
    return th, thp


def knoten(k, tmax=5):
    h = mp.mpf(2) ** (-k)
    n = int(tmax / h)
    pts = []
    for j in range(-n, n + 1):
        t = j * h
        a = mp.pi * mp.sinh(t)
        # x = 1/(1+e^{-a}), 1-x = 1/(1+e^{a})
        dl = 1 / (1 + mp.exp(-a))
        dr = 1 / (1 + mp.exp(a))
        wgt = h * mp.pi * mp.cosh(t) * dl * dr
        pts.append((dl, dr, wgt))
    return pts


def entwirren(werte, per):
    """Stetigen Zweig wählen: Werte sind bis auf Vielfache von per bestimmt; Knoten nach x sortiert."""
    out = [werte[0]]
    for v in werte[1:]:
        k = mp.nint(mp.re((out[-1] - v) / per)) if False else None
        # per ist komplex: kleinste |v + m per − vorher| über m
        best = min(range(-3, 4), key=lambda m: abs(v + m * per - out[-1]))
        m0 = best
        # grosse Sprünge: grob über Projektion
        r = (out[-1] - v) / per
        mr = int(mp.nint(mp.re(r)))
        cand = [mr + m for m in range(-2, 3)] + [m0]
        best = min(cand, key=lambda m: abs(v + m * per - out[-1]))
        out.append(v + best * per)
    return out


def lauf(k):
    pts = knoten(k)
    th, thp = [], []
    for dl, dr, w in pts:
        a, b = theta_und_ableitung(dl, dl, dr)
        th.append(a); thp.append(b)
    # Zweig: von der Mitte aus stetig in beide Richtungen
    mid = len(pts) // 2
    rechts = entwirren(th[mid:], PER)
    links = entwirren(th[mid::-1], PER)[::-1]
    th = links[:-1] + rechts
    lg = [mp.log(1 + v) for v in thp]
    lr = entwirren(lg[mid:], 2j * mp.pi)
    ll = entwirren(lg[mid::-1], 2j * mp.pi)[::-1]
    lg = ll[:-1] + lr
    W = [p[2] for p in pts]
    X = [p[0] for p in pts]

    def M(f):
        return mp.fsum(w * f(i) for i, w in enumerate(W))
    res = {
        "M_theta": M(lambda i: th[i]),
        "M_theta2": M(lambda i: th[i] ** 2),
        "M_theta3": M(lambda i: th[i] ** 3),
        "M_eLth": M(lambda i: mp.exp(L * th[i])),
        "M_e-Lth": M(lambda i: mp.exp(-L * th[i])),
        "M_e1th": M(lambda i: mp.exp(th[i])),
        "M_eith": M(lambda i: mp.exp(1j * th[i])),
        "M_log1thp": M(lambda i: lg[i]),
        "M_inv1thp": M(lambda i: mp.exp(-lg[i])),
        "M_sqrt1thp": M(lambda i: mp.exp(lg[i] / 2)),
        "M_pow_-2+i": M(lambda i: mp.exp((-2 + 1j) * lg[i])),
    }
    for m in (1, 2, 3):
        res[f"M_theta_e+{m}"] = M(lambda i: th[i] * mp.expjpi(2 * m * X[i]))
        res[f"M_theta_e-{m}"] = M(lambda i: th[i] * mp.expjpi(-2 * m * X[i]))
        res[f"M_inv1thp_e+{m}"] = M(lambda i: mp.exp(-lg[i]) * mp.expjpi(2 * m * X[i]))
    return res, len(pts)


if __name__ == "__main__":
    a = lade_theta()
    ergebnisse = {}
    for k in (5, 6, 7):
        t = time.time()
        r, n = lauf(k)
        ergebnisse[k] = r
        print(f"--- k={k}, {n} Knoten, {time.time()-t:.1f}s", flush=True)
    r = ergebnisse[7]; r6 = ergebnisse[6]
    mt = r["M_theta"]
    def zeile(name, wert, soll, alt):
        print(f"{name:28s} {mp.nstr(wert, 25):>60s}  |Abw|={mp.nstr(abs(wert-soll),3):>9s}  |k6-k7|={mp.nstr(abs(wert-alt),3)}")
    print("M[θ] =", mp.nstr(mt, 40), "  a_0 (etappe2) =", mp.nstr(a[0], 25), " Differenz/(2πi/L) =", mp.nstr((mt - a[0]) / PER, 12))
    zeile("M[θ²] vs M[θ]²", r["M_theta2"], mt ** 2, r6["M_theta2"])
    zeile("M[θ³] vs M[θ]³", r["M_theta3"], mt ** 3, r6["M_theta3"])
    zeile("M[e^θ] vs e^{M θ}", r["M_e1th"], mp.exp(mt), r6["M_e1th"])
    zeile("M[e^{iθ}] vs e^{i M θ}", r["M_eith"], mp.exp(1j * mt), r6["M_eith"])
    zeile("M[e^{Lθ}]·M[e^{-Lθ}] vs 1", r["M_eLth"] * r["M_e-Lth"], 1, r6["M_eLth"] * r6["M_e-Lth"])
    zeile("M[e^{Lθ}] vs e^{L a0}", r["M_eLth"], mp.exp(L * a[0]), r6["M_eLth"])
    zeile("M[log(1+θ')] vs 0", r["M_log1thp"], 0, r6["M_log1thp"])
    zeile("M[1/(1+θ')] vs 1", r["M_inv1thp"], 1, r6["M_inv1thp"])
    zeile("M[(1+θ')^{1/2}] vs 1", r["M_sqrt1thp"], 1, r6["M_sqrt1thp"])
    zeile("M[(1+θ')^{-2+i}] vs 1", r["M_pow_-2+i"], 1, r6["M_pow_-2+i"])
    for m in (1, 2, 3):
        zeile(f"M[θ e^(+2πi{m}x)] vs 0", r[f"M_theta_e+{m}"], 0, r6[f"M_theta_e+{m}"])
        zeile(f"M[θ e^(-2πi{m}x)] vs a_{m}", r[f"M_theta_e-{m}"], a[m], r6[f"M_theta_e-{m}"])
        zeile(f"M[e^(+2πi{m}x)/(1+θ')] vs 0", r[f"M_inv1thp_e+{m}"], 0, r6[f"M_inv1thp_e+{m}"])
    out = {str(k): {n: [mp.nstr(mp.re(v), dps), mp.nstr(mp.im(v), dps)] for n, v in rr.items()} for k, rr in ergebnisse.items()}
    (HIER / "daten" / f"nullmoden_{dps}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
