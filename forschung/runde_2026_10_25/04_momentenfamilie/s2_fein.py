"""Feine Abtastung der Kurve Γ_T = θ_T((0,1)) und Doppelpunkt-Suche.

Gleichmässiges Gitter in x plus geometrische Verdichtung an den Enden. θ_T(x) = 1 + h(log sexp(x+p)) − x (wie e4_lib).
Aufruf: python s2_fein.py <dps> <N> <fall> [eps]
  fall: kneser | a | c | d | moeb<a> | sin<eps> (p = eps·sin 2πx) | sin2<eps> (p = eps·sin 4πx)
Ausgabe: Schnittpaare des Polygonzugs (mit x-Werten), kleinster Abstand zwischen fernen Kurvenstücken, Fläche, Lage von c.
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E4 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe4"
sys.path.insert(0, str(E4))
import e4_lib as E  # noqa: E402

dps = int(sys.argv[1]); N = int(sys.argv[2]); fall = sys.argv[3]
E.init(dps)


def stoerung(fall):
    if fall == "kneser":
        return ("trig", [])
    if fall in ("a", "b", "c", "d"):
        d = json.loads((E4 / "daten" / f"gegenbeispiel_{fall}_40.json").read_text(encoding="utf-8"))
        return ("trig", [(int(m), mp.mpf(a), mp.mpf(b)) for m, a, b in d["p"]])
    if fall.startswith("moeb"):
        return ("moeb", mp.mpf(fall[4:]))
    if fall.startswith("sin2"):
        return ("trig", [(2, mp.mpf(fall[4:]), 0)])
    if fall.startswith("sin"):
        return ("trig", [(1, mp.mpf(fall[3:]), 0)])
    raise ValueError(fall)


def gitter(N):
    xs = [mp.mpf(j) / N for j in range(1, N)]
    rand = [mp.mpf(10) ** (-e / 4) for e in range(int(4 * mp.log10(N)) + 1, 4 * 12)]
    xs = sorted(set(rand + xs + [1 - r for r in rand]))
    return xs


def theta_T(stoer, xs):
    vals = []
    for x in xs:
        xs_ = x if x <= mp.mpf("0.5") else x - 1
        p, dp = E.p_wert(stoer, xs_)
        s1, ds = E.sexp_reihe(xs_ + p)
        hv, hd = E.B.h_mit_ableitung(mp.log1p(s1))
        vals.append(1 + hv - xs_)
    # stetiger Zweig, ausgehend von der Mitte
    mid = min(range(len(xs)), key=lambda i: abs(xs[i] - mp.mpf("0.5")))
    links = E.entwirren(vals[mid::-1], E.PER)[::-1]
    rechts = E.entwirren(vals[mid:], E.PER)
    return links[:-1] + rechts


def kreuz(o, a, b):
    return (a - o).real * (b - o).imag - (a - o).imag * (b - o).real


def schnitte(pts):
    n = len(pts); out = []
    bb = [(min(pts[i].real, pts[i + 1].real), max(pts[i].real, pts[i + 1].real),
           min(pts[i].imag, pts[i + 1].imag), max(pts[i].imag, pts[i + 1].imag)) for i in range(n - 1)]
    for i in range(n - 1):
        for j in range(i + 2, n - 1):
            if bb[i][1] < bb[j][0] or bb[j][1] < bb[i][0] or bb[i][3] < bb[j][2] or bb[j][3] < bb[i][2]:
                continue
            a, b, c, d = pts[i], pts[i + 1], pts[j], pts[j + 1]
            if kreuz(c, d, a) * kreuz(c, d, b) < 0 and kreuz(a, b, c) * kreuz(a, b, d) < 0:
                out.append((i, j))
    return out


def windung(pts, z):
    """Umlaufzahl des (über ∞ geschlossenen) Polygons um z: Summe der Winkeländerungen."""
    s = 0.0
    import cmath
    for i in range(len(pts) - 1):
        s += cmath.phase((pts[i + 1] - z) / (pts[i] - z))
    # Schluss: vom letzten zum ersten Punkt (beide nahe der Spitze bei ∞) über die kurze Verbindung
    s += cmath.phase((pts[0] - z) / (pts[-1] - z))
    return s / (2 * cmath.pi)


if __name__ == "__main__":
    st = stoerung(fall)
    xs = gitter(N)
    th = theta_T(st, xs)
    pts = [complex(t) for t in th]
    X = [float(x) for x in xs]
    sch = schnitte(pts)
    A = sum((pts[i].conjugate() * pts[i + 1]).imag for i in range(len(pts) - 1)) / 2
    A += (pts[-1].conjugate() * pts[0]).imag / 2
    print(f"{fall}: Punkte {len(pts)}, Schnitte {len(sch)}, Fläche {A:.6e}")
    for i, j in sch[:10]:
        print(f"   Schnitt: x ∈ [{X[i]:.5f},{X[i+1]:.5f}] mit x ∈ [{X[j]:.5f},{X[j+1]:.5f}],  Punkt ≈ {pts[i]:.6f}")
    # Mittelwert c mit Trapez auf dem gleichmässigen Teil ist ungenau; nur Lage-Information
    (HIER / "daten").mkdir(exist_ok=True)
    (HIER / "daten" / f"fein_{fall}_{dps}_{N}.json").write_text(
        json.dumps({"x": X, "re": [p.real for p in pts], "im": [p.imag for p in pts], "schnitte": sch}), encoding="utf-8")
    for z in [complex(26.193289, 5.425333), complex(26.25, 5.1), complex(26.3, 4.9), complex(26.15, 5.3)]:
        print(f"   Umlaufzahl um {z}: {windung(pts, z):+.4f}")
