"""Etappe 1, Teil A2: parabolische Konstanten ohne Kneser-Motor.

Fatou-Koordinate α(f(x)) = α(x) + 1 am parabolischen Punkt 0 für f = x·eˣ (Produkt-Tetration) und f = eˣ − 1,
asymptotisch α(x) = A/x + B·ln x + Σ c_k x^k (exakte Brüche), dann Abstieg mit f⁻¹ bis x ≈ 1/M und zurück.
Konstanten: α_f(1), prodtet_e(½) = f^{½}(1) für x·eˣ, Écalle-Halbiterierte von eˣ − 1 bei 1, f^{½}(1).
Zwei Abstiegstiefen M und zwei Reihenlängen K zum Abgleich. Ausgabe daten/fatou.json.
"""
import json
from fractions import Fraction as Fr
from math import factorial
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
KMAX = 70


def ser_inv(h, n):          # 1/h, h[0] = 1
    r = [Fr(0)] * n
    r[0] = Fr(1)
    for m in range(1, n):
        r[m] = -sum(h[j] * r[m - j] for j in range(1, m + 1) if j < len(h))
    return r


def ser_mul(a, b, n):
    c = [Fr(0)] * n
    for i in range(n):
        if a[i]:
            for j in range(n - i):
                c[i + j] += a[i] * b[j]
    return c


def ser_log(h, n):          # ln h, h[0] = 1:  (ln h)' = h'/h
    hi = ser_inv(h, n)
    hp = [Fr(k + 1) * h[k + 1] for k in range(n - 1)] + [Fr(0)]
    q = ser_mul(hp, hi, n)
    return [Fr(0)] + [q[k - 1] / k for k in range(1, n)]


def fatou_reihe(h, K):
    """h = f(x)/x als Koeffizientenliste; Rückgabe A, B, c[1..K] mit α = A/x + B ln x + Σ c_k x^k."""
    n = K + 3
    hi = ser_inv(h, n + 1)
    lh = ser_log(h, n)
    hk = [None, [Fr(1)] + [Fr(0)] * (n - 1)]
    pw = [Fr(1)] + [Fr(0)] * (n - 1)
    hpow = [pw]
    for k in range(1, K + 2):
        pw = ser_mul(pw, h, n)
        hpow.append(pw)
    # E_m = A·hi[m+1] + B·lh[m] + Σ_{k=1}^{m−1} c_k·hpow[k][m−k] − [m == 0]
    A = Fr(1) / hi[1] * 1      # m = 0: A·hi[1] − 1 = 0
    A = 1 / hi[1]
    B = -A * hi[2] / lh[1]    # m = 1
    c = [None]
    for m in range(2, K + 2):
        s = A * hi[m + 1] + B * lh[m] + sum(c[k] * hpow[k][m - k] for k in range(1, m - 1))
        kk = m - 1
        c.append(-s / hpow[kk][1])
    return A, B, c


def h_xexp(n):
    return [Fr(1, factorial(k)) for k in range(n)]


def h_expm1(n):
    return [Fr(1, factorial(k + 1)) for k in range(n)]


ergebnis = {}
for name, hfun, f, finv in (
        ("xexp", h_xexp, lambda x: x * mp.exp(x), lambda y: mp.lambertw(y).real),
        ("expm1", h_expm1, lambda x: mp.expm1(x), lambda y: mp.log1p(y))):
    A, B, c = fatou_reihe(hfun(KMAX + 6), KMAX)
    ergebnis[name + ".A"] = str(A)
    ergebnis[name + ".B"] = str(B)
    ergebnis[name + ".c"] = [str(v) for v in c[1:]]
    print(name, "A =", A, "B =", B, "c1..6 =", [str(v) for v in c[1:7]], flush=True)
    werte = {}
    for (M, K, dps) in ((4000, 60, 190), (8000, 70, 200)):
        mp.mp.dps = dps
        cm = [mp.mpf(v.numerator) / v.denominator if v is not None else None for v in c]

        def alpha(x):
            return A.numerator / mp.mpf(A.denominator) / x + B.numerator / mp.mpf(B.denominator) * mp.log(x) + mp.fsum(cm[k] * x ** k for k in range(1, K + 1))

        x = mp.mpf(1)
        for _ in range(M):
            x = finv(x)
        a1 = alpha(x) + M
        # letzter Summand als Fehlermass der asymptotischen Reihe
        letzter = abs(cm[K] * x ** K)
        ziel = alpha(x) + mp.mpf(1) / 2
        t = mp.findroot(lambda s: alpha(s) - ziel, x + x ** 2 / (2 * abs(A.numerator / mp.mpf(A.denominator))))
        for _ in range(M):
            t = f(t)
        werte[(M, K)] = (a1, t, letzter)
        print(f"  M={M} K={K}: alpha(1) = {mp.nstr(a1, 50)}  f^(1/2)(1) = {mp.nstr(t, 50)}  letzter Term {mp.nstr(letzter, 3)}", flush=True)
    (a_1, t_1, _), (a_2, t_2, _) = werte.values()
    mp.mp.dps = 200
    ergebnis[name + ".alpha(1)"] = mp.nstr(a_2, 180)
    ergebnis[name + ".halb(1)"] = mp.nstr(t_2, 180)
    ergebnis[name + ".abweichung_alpha"] = mp.nstr(abs(a_1 - a_2), 5)
    ergebnis[name + ".abweichung_halb"] = mp.nstr(abs(t_1 - t_2), 5)
    print(name, "Abgleich:", ergebnis[name + ".abweichung_alpha"], ergebnis[name + ".abweichung_halb"], flush=True)

(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / "fatou.json").write_text(json.dumps(ergebnis, indent=1), encoding="utf-8")
