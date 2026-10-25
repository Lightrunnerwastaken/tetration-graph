"""F3: Woher die Divergenz der Sprungreihe kommt.
H(t) = sexp(t) - sexp(t-1) - 1 = e^{g} - g - 1, g(t) = sexp(t-1); g hat bei t = -1 einen Logarithmus (sexp-Singularitaet bei -2).
Zerlegung H = -log(1+t) + R(t).  r_k = [t^k] R = J_k/k! - (-1)^k/k.  Wurzeltest fuer den Radius von R.
Singulaerer Anteil der Fourier-Koeffizienten exakt: sigma_n = -i int_0^oo -log(1-iy) e^{-2 pi n y} dy = i E1(2 pi i n)/(2 pi n).
Rest a_n - sigma_n gegen sum_k k! r_k/(2 pi i n)^{k+1} (beste Kuerzung).  Dazu: Konstantensuche fuer Extremwerte (Beinahe-Treffer).
"""
import json
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
mp.mp.dps = 32
C0 = [mp.mpf(s) for s in json.load(open(HIER / "daten" / "taylor_0_1_512_60.json"))][:200]
Cm = [mp.mpf(s) for s in json.load(open(HIER / "daten" / "taylor_-1_0.6_400_60.json"))][:120]
K = 80
h = [C0[k] - Cm[k] for k in range(K)]  # [t^k] (sexp(t) - sexp(t-1)), k>=1; k=0: 1
r = [mp.mpf(0)] + [h[k] - (-1) ** k / mp.mpf(k) for k in range(1, K)]
print("r_k = J_k/k! - (-1)^k/k  und |r_k|^(-1/k) (Radius-Schaetzung):")
for k in list(range(1, 12)) + list(range(12, K, 6)):
    print(f"  k={k:2d}  r_k = {mp.nstr(r[k], 14):>22}   |r_k|^(-1/k) = {mp.nstr(abs(r[k]) ** (-mp.mpf(1) / k), 8)}")
# Fehlergrenze der Koeffizienten: Taylor(-1) auf Kreis 0.6, 60 Stellen -> ~1e-60/0.6^k
print("  Fehlerschaetzung Taylor(-1) bei k=79: ~", mp.nstr(mp.mpf(10) ** -58 / mp.mpf("0.6") ** 79, 3))

# direkte Fourier-Koeffizienten (wie f1)
P = lambda u: u + 1 - mp.polyval(C0[::-1], u)
_gl = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)
KN, GW = [], []
for j in range(8):
    lo, hi = mp.mpf(-1) + mp.mpf(j) / 8, mp.mpf(-1) + mp.mpf(j + 1) / 8
    for x, w in _gl:
        KN.append((lo + hi) / 2 + (hi - lo) / 2 * x); GW.append(w * (hi - lo) / 2)
PK = [P(u) for u in KN]
a = lambda n: mp.fsum(w * p * mp.expj(-2 * mp.pi * n * u) for u, w, p in zip(KN, GW, PK))
print("\nn | a_n (direkt) | sigma_n = i E1(2 pi i n)/(2 pi n) | Rest - beste Kuerzung (k<=kmax) | e^{-2 pi n * 2}")
for n in [1, 2, 3, 4, 5, 6, 8]:
    an = a(n)
    s = 2 * mp.pi * n
    sig = 1j * mp.e1(1j * s) / s
    rest = an - sig
    best = None
    teil = mp.mpc(0)
    for k in range(1, K):
        teil += mp.factorial(k) * r[k] / (1j * s) ** (k + 1)
        e = abs(rest - teil)
        if best is None or e < best[1]:
            best = (k, e)
    print(f"  n={n}: a_n = {mp.nstr(an, 16)}  sigma_n = {mp.nstr(sig, 12)}  |Rest-Reihe| = {mp.nstr(best[1], 3)} (k<={best[0]})  e^(-4 pi n) = {mp.nstr(mp.exp(-2*s), 3)}")
    if n == 1:
        print(f"         a_1 voll: {mp.nstr(an.real, 28)} + {mp.nstr(an.imag, 28)} i")
    if n == 2:
        print(f"         a_2 voll: {mp.nstr(an.real, 28)} + {mp.nstr(an.imag, 28)} i")

# Beinahe-Treffer-Pruefung
print("\nKonstantensuche:")
S1 = lambda u: mp.polyval([k * C0[k] for k in range(1, len(C0))][::-1], u)
S2 = lambda u: mp.polyval([k * (k - 1) * C0[k] for k in range(2, len(C0))][::-1], u)
xw = mp.findroot(S2, -0.52); m = S1(xw)
print("  m = min sexp' =", mp.nstr(m, 25), "  22/23 =", mp.nstr(mp.mpf(22) / 23, 25), "  Differenz", mp.nstr(m - mp.mpf(22) / 23, 5))
ex = [mp.findroot(lambda u: S1(u) - 1, u0) for u0 in (-0.8, -0.2)]
werte = {"u_min": ex[0], "u_max": ex[1], "P_min": P(ex[0]), "P_max": P(ex[1]), "x_w": xw, "m": m, "c": C0[1],
         "a_0": mp.mpf(1) / 2 - mp.quad(lambda u: mp.polyval(C0[::-1], u), [-1, 0]), "u_min+u_max": ex[0] + ex[1], "P_max+P_min": P(ex[0]) + P(ex[1])}
for name, v in werte.items():
    with mp.workdps(20):
        idf = mp.identify(v, ["pi", "e", "log(2)"], tol=mp.mpf(10) ** -15)
    print(f"  {name:12s} = {mp.nstr(v, 22):>26}  identify(15 St.): {idf}")
