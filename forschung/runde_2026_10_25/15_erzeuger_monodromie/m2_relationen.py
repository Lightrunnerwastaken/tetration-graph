"""m2: Relationen in der Monodromiegruppe G direkt geometrisch pruefen (ohne Kodierung).

Ein Wort in a, b1, b2 wird als zusammengesetzter Weg (erst links, dann rechts) von mehreren Blaettern aus hochgehoben.
Kehrt die Hebung auf JEDEM Probeblatt zurueck, wirkt das Wort (auf diesen Blaettern) trivial.
Vorhersage aus Satz M2: [b1, a^n b1 a^-n] = 1 (n != 0), [b2, b1 b2 b1^-1] = 1, [b2, a b1 a^-1] = 1;
dagegen [a, b1] != 1, [b1, b2] != 1  -> G ist nicht frei.
Aufruf: python m2_relationen.py   (Gleitkomma, nutzt m1 Einstellung A)
"""
import sys
sys.argv = [sys.argv[0], "A"]
import m1_monodromie as m

W = {"a": m.WEGE["a"], "b1": m.WEGE["b1"], "b2": m.WEGE["b2"]}


def weg(wort):
    """wort: Liste von (name, +-1)."""
    pkt = [m.L]
    for name, e in wort:
        w = W[name] if e > 0 else list(reversed(W[name]))
        pkt += w[1:]
    return pkt


def inv(wort):
    return [(n, -e) for n, e in reversed(wort)]


def komm(x, y):
    return x + y + inv(x) + inv(y)


A = [("a", 1)]; B1 = [("b1", 1)]; B2 = [("b2", 1)]
konj = lambda g, x: g + x + inv(g)
tests = {
    "[b1, a b1 a^-1]": komm(B1, konj(A, B1)),
    "[b1, a^-2 b1 a^2]": komm(B1, konj(A[:] * 0 + [("a", -1), ("a", -1)], B1)),
    "[b2, b1 b2 b1^-1]": komm(B2, konj(B1, B2)),
    "[b2, a b1 a^-1]": komm(B2, konj(A, B1)),
    "[a, b1]": komm(A, B1),
    "[b1, b2]": komm(B1, B2),
    "[a, b2]": komm(A, B2),
}
blaetter = [(), (1,), (-1,), (0, 1), (0, -1), (0, 0, 1), (1, 1), (2, -1, 3)]
us = {c: m.dekodiere(c) for c in blaetter}
for name, wort in tests.items():
    sprung = []
    for c, u in us.items():
        u2 = m.verfolge(u, weg(wort))
        sprung.append(abs(u2 - u))
    triv = max(sprung) < 1e-8
    print(f"{name:22s} max Sprung {max(sprung):.2e}  min Sprung {min(sprung):.2e}  -> {'trivial auf allen Probeblaettern' if triv else 'NICHT trivial'}")
    sys.stdout.flush()
