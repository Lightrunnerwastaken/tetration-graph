"""X25c — wie x25, aber N und Antworten wählbar (Aufruf: python x25c_bogen_N.py N); Bogenmodell für (Fam-Eck): Eckantwort δσ₁/σ₁ als lineare Funktion GEGLÄTTETER Fixpunktdaten.
Lineare Fixpunktdaten einer Störung η (erste Ordnung): δL = η(L)/(1 − μ), δr = Re(δL·L̄/ρ), δμ = η'(L) + f''(L)·δL.
Bogenmodell: η(L), η'(L) ersetzt durch die harmonische (Poisson-)Glättung η(L(1−ε)), η'(L(1−ε)) — für ein Polynom = Mittel über einen Bogen
der Länge ~ερ um L. Anpassung von (a, c, d) in  Antwort = a·δr + c·Re δμ + d·Im δμ  an den glatten Formen (1, x, q, x², x³, x⁴; N = 400),
Vorhersage für q², q·x^{N/2}, p_N und die Zufallsstörungen st/sv (Seeds 1–3, wie x18) — für ε = 0 … 0.2.
Gemessene Antworten (x18-Logs, N = 400, in Einheiten R).
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N = int(sys.argv[1])
mp.mp.dps = 50
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
PD = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"][str(N)]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(Lm)
mu = mp.sqrt(Lm)
f2L = sum(k * (k - 1) * mp.mpf(A_str[k]) * Lm ** (k - 2) for k in range(2, len(A_str)))   # f''(L) (Reihe auf dem Rand)
rho_f = float(rho)


def sup_rand(c):
    z = rho_f * np.exp(1j * np.linspace(0, 2 * np.pi, 4001)[:-1])
    return np.max(np.abs(np.polynomial.polynomial.polyval(z, np.array([float(v) for v in c]))))


TNf = [mp.mpf(s) for s in A_str[: N + 1]]
eta_p = [mp.mpf(a) - b for a, b in zip(PD, TNf)]
R = sup_rand(eta_p)
LRe = float(mp.re(Lm))
q = [rho_f ** 2, -2 * LRe, 1.0]


def normiert(g):
    g = np.concatenate([np.array(g, float), np.zeros(N + 1 - len(g))])
    return g * (R / sup_rand(g))


def zufall(kind, seed):
    rng = np.random.default_rng(seed)
    lo, hi = {"t": (0, N // 2), "v": (0, N)}[kind]
    g = np.zeros(N + 1)
    g[lo:hi + 1] = rng.standard_normal(hi - lo + 1) / rho_f ** np.arange(lo, hi + 1)
    return g * (R / sup_rand(g))


ANTW = {
    200: {"1": 32.0, "x": 4.26, "q": 1.96, "x²": -25.3, "x³": -20.1, "x⁴": 7.74, "p_N": -0.127,
          "st1": 2.64, "st2": 2.72, "st3": 2.66, "sv1": 1.55, "sv2": 1.78, "sv3": 1.72},
    400: {"1": 131.0, "x": 15.9, "q": 7.27, "x²": -106.0, "x³": -80.0, "x⁴": 39.0, "p_N": -0.338,
          "st1": 9.18, "st2": 3.15, "st3": 7.17, "sv1": 6.23, "sv2": 1.78, "sv3": 4.80},
    600: {"1": 230.0, "x": 27.5, "q": 12.1, "x²": -188.0, "x³": -140.0, "x⁴": 74.2, "p_N": -0.337, "st1": 16.1, "st2": 1.07, "st3": 10.6, "sv2": 0.77},
}[N]
FORM = {"1": [1.0], "x": [0, 1.0], "q": q, "x²": [0, 0, 1.0], "x³": [0, 0, 0, 1.0], "x⁴": [0, 0, 0, 0, 1.0]}
formen = {}
for k, v in ANTW.items():
    if k in FORM:
        formen[k] = (normiert(FORM[k]), v)
    elif k == "p_N":
        formen[k] = (None, v)
    else:
        formen[k] = (zufall(k[1], int(k[2])), v)
glatt = [k for k in ["1", "x", "q", "x²", "x³", "x⁴"] if k in ANTW]


def daten(c, eps):
    z = Lm * (1 - eps)
    v = mp.mpc(0)
    dv = mp.mpc(0)
    for a in reversed(c):
        dv = dv * z + v
        v = v * z + a
    v, dv = v / R, dv / R
    dL = v / (1 - mu)
    dr = mp.re(dL * mp.conj(Lm) / rho)
    dmu = dv + f2L * dL
    return [float(dr), float(mp.re(dmu)), float(mp.im(dmu))]


for eps in (0.0, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15):
    X, y = {}, {}
    for name, (c, resp) in formen.items():
        coef = eta_p if name == "p_N" else [mp.mpf(float(t)) for t in c]
        X[name] = daten(coef, eps)
        y[name] = resp
    A = np.array([X[n] for n in glatt])
    b = np.array([y[n] for n in glatt])
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    fehler_glatt = np.max(np.abs(A @ sol - b)) / np.max(np.abs(b))
    vor = {n: float(np.dot(sol, X[n])) for n in formen if n not in glatt}
    print(f"ε = {eps}: a, c, d = {sol[0]:.1f}, {sol[1]:.1f}, {sol[2]:.1f}; Anpassfehler glatt {fehler_glatt:.3f}", flush=True)
    res = np.array([vor[n] - y[n] for n in vor])
    print(f"    RMS Rest {np.sqrt(np.mean(res ** 2)):.2f}; Vorhersage/gemessen: " + "; ".join(f"{n}: {vor[n]:+.2f}/{y[n]:+.2f}" for n in vor), flush=True)
