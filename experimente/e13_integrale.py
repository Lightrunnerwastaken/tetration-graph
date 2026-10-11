"""exp.integrale — Integralidentitäten an Kneser (Basis e, Motor) und regulär (√2, exakt).

  Young:        ∫_0^1 T + ∫_1^a A = a                      (A = slog)
  unten:        ∫_{-1}^0 T = −∫_0^1 A
  andere Form:  ∫_0^1 T = 1 − ln a · ∫_0^1 a^u A(u) du
  Produkt:      ∫_{-1}^0 T'(u) T(u+1) du = (a − 1)/ln a      (T' unabhängig über Cauchy-Integrale)
Gauss–Legendre mit 40 Knoten; Gegencheck mit 30 Knoten. PSLQ-Suche nach geschlossener Form für ∫_0^1 T_e.
"""
import mpmath as mp

from gemeinsam import Lauf, cache_speichern, cauchy_ableitungen, gauss_legendre, kneser, plt
from regulaer import Regulaer

mp.mp.dps = 40
lauf = Lauf("exp.integrale", __file__)
grenze = mp.mpf(10) ** -28


def pruefe_basis(name, T, A, a, n=40):
    lna = mp.log(a)
    I_T01 = gauss_legendre(T, mp.mpf(0), mp.mpf(1), n)
    I_T01b = gauss_legendre(T, mp.mpf(0), mp.mpf(1), n - 10)
    I_A1a = gauss_legendre(A, mp.mpf(1), a, n)
    I_Tm10 = gauss_legendre(T, mp.mpf(-1), mp.mpf(0), n)
    I_A01 = gauss_legendre(A, mp.mpf(0), mp.mpf(1), n)
    I_aA = gauss_legendre(lambda us: [mp.power(a, u) * v for u, v in zip(us, A(us))], mp.mpf(0), mp.mpf(1), n)
    lauf.pruefe(f"{name}: Quadratur stabil (40 vs 30 Knoten)", abs(I_T01 - I_T01b) < grenze, mp.nstr(abs(I_T01 - I_T01b), 3))
    lauf.pruefe(f"{name}: ∫_0^1 T + ∫_1^a A = a", abs(I_T01 + I_A1a - a) < grenze, mp.nstr(abs(I_T01 + I_A1a - a), 3))
    lauf.pruefe(f"{name}: ∫_-1^0 T = −∫_0^1 A", abs(I_Tm10 + I_A01) < grenze, mp.nstr(abs(I_Tm10 + I_A01), 3))
    lauf.pruefe(f"{name}: ∫_0^1 T = 1 − ln a·∫_0^1 a^u A(u) du", abs(I_T01 - (1 - lna * I_aA)) < grenze, mp.nstr(abs(I_T01 - (1 - lna * I_aA)), 3))
    return {"int_T_0_1": I_T01, "int_A_1_a": I_A1a, "int_T_m1_0": I_Tm10, "int_A_0_1": I_A01, "int_auA": I_aA}


# ------------------------------------------------------------------ Kneser e
Te = lambda zs: kneser("sexp", "exp(1)", zs, dps=40)
Ae = lambda ys: kneser("slog", "exp(1)", ys, dps=40)
werte_e = pruefe_basis("Kneser e", Te, Ae, mp.e)
# Produkt-Integral mit unabhängig bestimmter Ableitung
x, w = __import__("gemeinsam")._gl_knoten(30)
knoten = [(xi - 1) / 2 for xi in x]                      # auf [−1, 0]
abl = [cauchy_ableitungen(Te, u, min(mp.mpf("0.3"), (u + 2) / 3), 1, 64)[1] for u in knoten]
T_hoch = Te([u + 1 for u in knoten])
I_prod = mp.fsum(wi * d * mp.re(t) for wi, d, t in zip(w, abl, T_hoch)) / 2
lauf.pruefe("Kneser e: ∫_-1^0 T'(u)T(u+1) du = (e − 1)/ln e", abs(I_prod - (mp.e - 1)) < mp.mpf(10) ** -25, mp.nstr(abs(I_prod - (mp.e - 1)), 3))

# ------------------------------------------------------------------ Regulär √2
R = Regulaer(mp.sqrt(2), N=70, tol=mp.mpf("1e-3"))
Tr = lambda zs: [R.tet(z) for z in zs]
Ar = lambda ys: [R.slog(y) if y > 0 else R.slog(mp.power(R.b, y)) - 1 for y in ys]
werte_r = pruefe_basis("regulär √2", Tr, Ar, R.b, n=30)

# ------------------------------------------------------------------ Geschlossene Form? (Gegenbeispielsuche)
# Regel: 9 Zahlen × log10(maxcoeff) muss deutlich unter der Genauigkeit (~33 Stellen) bleiben.
I = werte_e["int_T_0_1"]
kand = lambda: [I, 1, mp.e, mp.e ** 2, mp.pi, mp.euler, mp.log(2), mp.re(-mp.lambertw(-1, -1)), mp.im(-mp.lambertw(-1, -1))]
mp.mp.dps = 30
rel_gross = mp.pslq(kand(), maxcoeff=10 ** 5, maxsteps=10 ** 6)       # bewusst zu grob: zeigt die Scheinrelation
mp.mp.dps = 38
rel = mp.pslq(kand(), tol=mp.mpf(10) ** -31, maxcoeff=10 ** 3, maxsteps=10 ** 6)
mp.mp.dps = 40
def rest(r):
    return None if r is None else abs(mp.fsum(c * v for c, v in zip(r, kand())))
lauf.pruefe("PSLQ zu grob (Koeff. ≤ 1e5 bei 30 Stellen) liefert nur eine Scheinrelation (Rest gross bei 40 Stellen)",
            rel_gross is None or rest(rel_gross) > mp.mpf(10) ** -28, mp.nstr(rest(rel_gross), 3) if rel_gross else "keine")
lauf.pruefe("PSLQ (38 Stellen, Toleranz 1e-31, Koeff. ≤ 1e3): keine echte Relation von ∫_0^1 T_e zu 1, e, e², π, γ, ln 2, Re L, Im L",
            rel is None or rest(rel) > mp.mpf(10) ** -28, "keine" if rel is None else f"{rel}, Rest {mp.nstr(rest(rel), 3)}")

lauf.daten = {"kneser_e": {k: mp.nstr(v, 30) for k, v in werte_e.items()}, "regulaer_sqrt2": {k: mp.nstr(v, 30) for k, v in werte_r.items()},
              "produkt_integral_e": mp.nstr(I_prod, 30), "pslq_1e3": rel, "pslq_1e5_schein": rel_gross}

# ------------------------------------------------------------------ Bild: Flächenbild der Young-Identität
xs = [mp.mpf(-1.2) + i * mp.mpf("0.02") for i in range(121)]
ys = Te(xs)
fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.plot([float(x) for x in xs], [float(mp.re(y)) for y in ys], color="#22408f", lw=2, label="y = T_e(x)")
xf = [float(x) for x in xs if 0 <= x <= 1]
yf = [float(mp.re(y)) for x, y in zip(xs, ys) if 0 <= x <= 1]
ax.fill_between(xf, 0, yf, color="#7fa8ff", alpha=0.45, label=f"∫₀¹ T_e = {float(werte_e['int_T_0_1']):.6f}")
ax.fill_betweenx(yf, 0, xf, color="#f2b84b", alpha=0.45, label=f"∫₁ᵉ slog_e = {float(werte_e['int_A_1_a']):.6f}")
ax.plot([0, 1, 1], [mp.e, mp.e, 0], color="#555", lw=1, ls="--")
ax.set_xlim(-1.2, 1.2); ax.set_ylim(0, 3.2)
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_title("Blau + Gelb = Rechteck 1 × e: ∫₀¹T + ∫₁ᵉ slog = e")
ax.legend(fontsize=8, loc="upper left")
lauf.bild(fig)
cache_speichern()
lauf.zusammenfassung = (f"Kneser e: ∫₀¹T = {mp.nstr(werte_e['int_T_0_1'], 25)}, ∫₁ᵉ slog = {mp.nstr(werte_e['int_A_1_a'], 20)}, "
                        f"∫₋₁⁰T = {mp.nstr(werte_e['int_T_m1_0'], 20)}; alle Identitäten auf ≤ 1e-28 (Produktform 1e-25). "
                        f"Regulär √2: ∫₀¹T = {mp.nstr(werte_r['int_T_0_1'], 20)}. PSLQ ohne Relation.")
lauf.schreibe()
