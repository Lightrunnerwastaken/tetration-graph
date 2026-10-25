"""S2: Vorzeichen der Ableitungen f_K^(n) auf ganz R (numerisch, mpmath).
Taylor-Reihe von f_K an beliebigem reellem x0 aus den Koeffizienten um 0 (kneser_halbexp_koeff_1000.json, Fehler <= 4e-47):
  |x0| <= 1.0 : Umentwicklung der Reihe um 0
  x0 > 1.0    : f(x0+t) = exp(f(log(x0+t)))   (rekursiv)
  x0 < -1.0   : f(x0+t) = log(f(exp(x0+t)))
Ausgabe: fuer n = 1..NMAX das Minimum von f^(n) auf einem Gitter in [-12, 40] (log-Gitter rechts) und die Vorzeichenwechsel.
Gegenprobe: an x0 in [0.6, 1.0] beide Wege (Umentwicklung und exp∘f∘log) vergleichen.
"""
import json
import sys
from pathlib import Path
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
QUELLE = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_halbexp_koeff_1000.json")
A = [mp.mpf(s) for s in json.loads(QUELLE.read_text(encoding="utf-8"))["a"]][:330]
NMAX = 10
ORD = NMAX + 1


def mul(p, q):
    r = [mp.mpf(0)] * ORD
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j in range(ORD - i):
            r[i + j] += a * q[j]
    return r


def compose(F, G):
    """F: Taylor-Koeffizienten an G[0]; G: Reihe in t. Ergebnis F(G(t))."""
    H = [mp.mpf(0)] * ORD
    d = [mp.mpf(0)] + list(G[1:])
    pot = [mp.mpf(1)] + [mp.mpf(0)] * (ORD - 1)
    for k in range(ORD):
        for j in range(ORD):
            H[j] += F[k] * pot[j]
        pot = mul(pot, d)
    return H


def exp_reihe(c):
    e = mp.exp(c)
    return [e / mp.factorial(k) for k in range(ORD)]


def log_reihe(c):
    return [mp.log(c)] + [(-1) ** (k + 1) / (k * c ** k) for k in range(1, ORD)]


def taylor0(x0):
    """Taylor-Verschiebung per wiederholtem Horner-Schema (synthetische Division), O(K * ORD)."""
    x0 = mp.mpf(x0)
    b = list(A)
    out = []
    for j in range(ORD):
        n = len(b)
        q = [mp.mpf(0)] * n
        acc = mp.mpf(0)
        for k in range(n - 1, -1, -1):
            acc = acc * x0 + b[k]
            q[k] = acc
        out.append(q[0])
        b = q[1:]
    return out


def taylor(x0, weg=None):
    x0 = mp.mpf(x0)
    if weg == "direkt" or (weg is None and abs(x0) <= 1):
        return taylor0(x0)
    if x0 > 0:
        lg = log_reihe(x0)
        F = taylor(lg[0])
        return compose(exp_reihe(compose(F, lg)[0]), compose(F, lg))
    ex = [mp.exp(x0) / mp.factorial(k) for k in range(ORD)]
    F = compose(taylor(ex[0]), ex)
    return compose(log_reihe(F[0]), F)


def abl(x0):
    T = taylor(x0)
    return [T[n] * mp.factorial(n) for n in range(ORD)]


if __name__ == "__main__":
    # Gegenprobe
    for x0 in ["0.6", "0.8", "1.0"]:
        a1 = abl(x0)
        a2 = [c * mp.factorial(n) for n, c in enumerate(taylor(x0, weg="log"))] if False else None
        x = mp.mpf(x0)
        lg = log_reihe(x)
        F = taylor0(lg[0])
        G = compose(F, lg)
        T2 = compose(exp_reihe(G[0]), G)
        d = max(abs(a1[n] - T2[n] * mp.factorial(n)) / (1 + abs(a1[n])) for n in range(ORD))
        print(f"Gegenprobe x0={x0}: max rel. Abweichung der Ableitungen 0..{NMAX} = {mp.nstr(d, 3)}")
    gitter = [mp.mpf(i) / 10 for i in range(-120, 41)] + [mp.mpf(4) * mp.mpf(10) ** (mp.mpf(i) / 20) for i in range(1, 21)]
    werte = [(x, abl(x)) for x in gitter]
    for n in range(1, NMAX + 1):
        vz = [(x, w[n]) for x, w in werte]
        neg = [x for x, v in vz if v < 0]
        wechsel = [mp.nstr((vz[i][0] + vz[i + 1][0]) / 2, 4) for i in range(len(vz) - 1) if (vz[i][1] < 0) != (vz[i + 1][1] < 0)]
        mn = min(vz, key=lambda p: p[1])
        print(f"n={n:2d}: min f^(n) auf Gitter = {mp.nstr(mn[1], 6)} bei x={mp.nstr(mn[0], 5)};  negativ an {len(neg)} Gitterpunkten;"
              f" Wechsel bei {wechsel[:8]}")
