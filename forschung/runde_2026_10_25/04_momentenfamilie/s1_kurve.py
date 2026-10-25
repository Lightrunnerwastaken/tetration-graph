"""Runde 25.10.2026, Gebiet 04: Momentenfamilie.

Prueft numerisch die Bausteine des Jordan-Kurven-Arguments:
  (1) Die Kurve Γ_T = θ_T([0,1]) ist doppelpunktfrei (Kneser, Gegenbeispiele, Störungsfamilie).
  (2) Cauchy-Transformierte C_T(ζ) = M[1/(ζ − θ_T)] gegen 1/(ζ − c) links der Kurve (Ω) und innen (D).
  (3) Orientierung (umlaufene Fläche) und Lage von c.
Benutzt nur Lesen von etappe4/e4_lib.py (Kneser auf der reellen Achse aus Taylor-Daten, kein Motor).

Aufruf: python s1_kurve.py <dps> <k>
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E4 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe4"
sys.path.insert(0, str(E4))
import e4_lib as E  # noqa: E402

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
E.init(dps)


def lade(fall):
    d = json.loads((E4 / "daten" / f"gegenbeispiel_{fall}_40.json").read_text(encoding="utf-8"))
    return ("trig", [(int(m), mp.mpf(a), mp.mpf(b)) for m, a, b in d["p"]])


def segs_schneiden(P):
    """Anzahl echter Schnitte zwischen nicht benachbarten Strecken des offenen Polygonzugs P (Liste komplexer Zahlen).
    O(n^2) mit Grobfilter; n ~ 600 ist ok."""
    n = len(P)
    pts = [complex(z) for z in P]
    bb = []
    for i in range(n - 1):
        a, b = pts[i], pts[i + 1]
        bb.append((min(a.real, b.real), max(a.real, b.real), min(a.imag, b.imag), max(a.imag, b.imag)))

    def kreuz(o, a, b):
        return (a - o).real * (b - o).imag - (a - o).imag * (b - o).real

    treffer = []
    for i in range(n - 1):
        for j in range(i + 2, n - 1):
            if bb[i][1] < bb[j][0] or bb[j][1] < bb[i][0] or bb[i][3] < bb[j][2] or bb[j][3] < bb[i][2]:
                continue
            a, b, c, d = pts[i], pts[i + 1], pts[j], pts[j + 1]
            d1, d2 = kreuz(c, d, a), kreuz(c, d, b)
            d3, d4 = kreuz(a, b, c), kreuz(a, b, d)
            if d1 * d2 < 0 and d3 * d4 < 0:
                treffer.append((i, j))
    return treffer


def min_abstand(P, X, luecke=0.02):
    """kleinster Abstand |θ(x)−θ(y)| über Knotenpaare mit Kreisabstand von x,y >= luecke (Mass der Doppelpunktfreiheit)."""
    pts = [complex(z) for z in P]
    xs = [float(x) for x in X]
    best = (1e99, None, None)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            dx = abs(xs[i] - xs[j]); dx = min(dx, 1 - dx)
            if dx < luecke:
                continue
            dd = abs(pts[i] - pts[j])
            if dd < best[0]:
                best = (dd, xs[i], xs[j])
    return best


def flaeche(P):
    pts = [complex(z) for z in P]
    s = 0.0
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        s += (a.conjugate() * b).imag
    return s / 2


def analyse(name, stoer):
    X, W, th, _ = E.auswerten(stoer, k)
    M = lambda vals: mp.fsum(w * v for w, v in zip(W, vals))  # noqa: E731
    c = M(th)
    sch = segs_schneiden(th)
    ma = min_abstand(th, X)
    A = flaeche(th)
    re = [float(mp.re(t)) for t in th]; im = [float(mp.im(t)) for t in th]
    print(f"== {name}: c = {mp.nstr(c, 12)}  Knoten {len(th)}")
    print(f"   Re θ_T in [{min(re):.4f}, {max(re):.4f}], Im in [{min(im):.4f}, {max(im):.4f}]")
    print(f"   Endwerte θ_T(x→0+) = {mp.nstr(th[0], 8)},  θ_T(x→1−) = {mp.nstr(th[-1], 8)}")
    print(f"   Polygonschnitte: {len(sch)}  (Beispiele {sch[:3]});  min Abstand (|x−y|≥0.02): {ma[0]:.3e} bei x={ma[1]:.4f}, y={ma[2]:.4f}")
    print(f"   umlaufene Fläche (Polygon, offen geschlossen) = {A:.6e}")
    # Cauchy-Transformierte
    for zeta in [mp.mpc(-5, 0), mp.mpc(20, 5), mp.mpc(25, 10), mp.mpc(26, -3), mp.mpc(26, 5.4), mp.mpc(27, 5.0)]:
        C = M([1 / (zeta - t) for t in th])
        print(f"   ζ = {mp.nstr(zeta, 4):>14}:  C − 1/(ζ−c) = {mp.nstr(C - 1 / (zeta - c), 5)}")
    return X, th, c


if __name__ == "__main__":
    faelle = [("Kneser", ("trig", [])), ("Gegenbeispiel a (V2=0)", lade("a")), ("Gegenbeispiel c", lade("c")),
              ("Moebius a=0.3", ("moeb", 0.3))]
    out = {}
    for name, st in faelle:
        X, th, c = analyse(name, st)
        out[name] = {"x": [float(x) for x in X], "re": [float(mp.re(t)) for t in th], "im": [float(mp.im(t)) for t in th],
                     "c": [float(mp.re(c)), float(mp.im(c))]}
    (HIER / "daten").mkdir(exist_ok=True)
    (HIER / "daten" / f"kurven_{dps}_{k}.json").write_text(json.dumps(out), encoding="utf-8")
