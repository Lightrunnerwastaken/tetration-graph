"""exp.stufe5 — η₅ und e₅: wo die Pentation P_b(x) = b↑↑↑x die Gerade y = x berührt (TPID 10/11).

 Pentation wie in exp.pentation (regulär am negativen Fixpunkt der Kneser-Tetration T_b, Hilfsmodul stufen.py).
 Für jede Basis: Taylorreihe von P_b um x0 aus Werten auf einem Kreis (komplexe Pentation), Stelle x* mit P_b'(x*) = 1,
 g(b) = P_b(x*) − x*; η₅ = Nullstelle von g (Sekanten in b), e₅ = x*(η₅).
 TPID 10 (dantheman163, Tetration Forum tid 162, 31.10.2010): „Eulers“ e, e₄, e₅, … → 4 und „Etas“ η, η₄, η₅, … → 2.
 TPID 11 (tommy1729, 01.12.2010): η_n² < e_{n−1}.
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from stufen import p_abl, p_nullstelle, p_verschieben, p_wert, pentation

mp.mp.dps = 50
lauf = Lauf("exp.stufe5", __file__)
E4, ETA4 = mp.mpf("3.0885322718067176544821807826411"), mp.mpf("1.6353244967152763993453446183062")
ETA = mp.e ** (1 / mp.e)
R, NK = mp.mpf("0.25"), 64


def stufe(b, x0, L0, M=8):
    """Pentation zur Basis b; Minimum von P_b(x) − x nahe x0 aus der Taylorreihe."""
    L, lam, psi, c, kette, P = pentation(b, L0, mp.mpf("0.1"), M=M)
    zs = [x0 + R * mp.expjpi(2 * mp.mpf(j) / NK) for j in range(NK)]
    w = P(zs)
    a = [mp.re(mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / NK) for j in range(NK)) / NK) / R ** m for m in range(26)]
    g1 = p_abl(a)
    g1[0] -= 1
    t = p_nullstelle(g1)
    return {"x": x0 + t, "g": p_wert(a, t) - (x0 + t), "a": p_verschieben(a, t), "L": L, "lam": lam, "P": P}


# ------------------------------------------------------------------ η₅ durch Sekanten in b
b0, b1 = mp.mpf("1.72"), mp.mpf("1.73")
s0 = stufe(b0, mp.mpf("3.0"), mp.mpf("-1.67"))
s1 = stufe(b1, s0["x"], s0["L"])
for _ in range(14):
    b2 = b1 - s1["g"] * (b1 - b0) / (s1["g"] - s0["g"])
    b0, s0, b1 = b1, s1, b2
    s1 = stufe(b1, s0["x"], s0["L"])
    print(f"  b = {mp.nstr(b1, 30)}  g = {mp.nstr(s1['g'], 3)}  x* = {mp.nstr(s1['x'], 20)}")
    if abs(b1 - b0) < mp.mpf(10) ** -28:
        break
eta5, e5, a5 = b1, s1["x"], s1["a"]
lauf.pruefe("Berührung: P(e₅) = e₅, P'(e₅) = 1, P''(e₅) > 0", abs(a5[0] - e5) < mp.mpf(10) ** -24 and abs(a5[1] - 1) < mp.mpf(10) ** -24 and a5[2] > 0,
            f"η₅ = {mp.nstr(eta5, 25)}, e₅ = {mp.nstr(e5, 25)}, P''/2 = {mp.nstr(a5[2], 10)}")
tief = stufe(eta5, e5, s1["L"], M=10)
lauf.pruefe("unabhängig von der Rechentiefe (M = 8 gegen 10): Minimum bleibt bei 0 und e₅", abs(tief["g"]) < mp.mpf(10) ** -22 and abs(tief["x"] - e5) < mp.mpf(10) ** -18,
            f"g = {mp.nstr(tief['g'], 3)}, Δx = {mp.nstr(tief['x'] - e5, 3)}")
ganz = s1["P"]([mp.mpf(0), mp.mpf(1), mp.mpf(2)])
lauf.pruefe("Kontrolle ganzer Höhen bei η₅: P(0) = 1, P(1) = η₅, P(2) = η₅↑↑η₅", abs(ganz[0] - 1) + abs(ganz[1] - eta5) < mp.mpf(10) ** -24)

# Fixpunktzahl der Pentation: P_b(x) − x ab x = −6 in Schritten von 1/8, bis P über x davonwächst (P > 6);
# Stapel von je 4 Punkten, damit P am Stapelende höchstens T_b(≈6) ist (T_b(7.5) läuft im Motor über)
def schnitte(b, L0):
    s = stufe(b, mp.mpf("3.0"), L0)
    x, vor, n, stellen = mp.mpf(-6), None, 0, []
    while True:
        xs = [x + mp.mpf(k) / 8 for k in range(4)]
        for xi, yi in zip(xs, s["P"](xs)):
            d = yi - xi
            if vor is not None and d * vor < 0:
                n += 1
                stellen.append(mp.nstr(xi - mp.mpf(1) / 16, 3))
            vor = d
        if yi > 6 and d > 0:
            return n, stellen
        x = xs[-1] + mp.mpf(1) / 8


z_unter, z_ueber = schnitte(mp.mpf("1.70"), mp.mpf("-1.66")), schnitte(mp.mpf("1.74"), mp.mpf("-1.67"))
lauf.pruefe("Fixpunkte der Pentation: zwei für b = 1.70 < η₅, keiner für b = 1.74 > η₅ (P > x für alle x ≤ 0)",
            z_unter[0] == 2 and z_ueber[0] == 0, {"1.70": f"{z_unter[0]} bei {z_unter[1]}", "1.74": z_ueber[0]})

# TPID 10 / 11
folge_e = [mp.e, E4, e5]
folge_eta = [ETA, ETA4, eta5]
lauf.pruefe("TPID 11: η_n² < e_{n−1} (η₄² < e, η₅² < e₄)", ETA4 ** 2 < mp.e and eta5 ** 2 < E4,
            f"η₄² = {mp.nstr(ETA4 ** 2, 8)} < {mp.nstr(mp.e, 8)};  η₅² = {mp.nstr(eta5 ** 2, 8)} < {mp.nstr(E4, 8)}")
lauf.pruefe("Etas wachsen: η < η₄ < η₅ < 2", folge_eta[0] < folge_eta[1] < folge_eta[2] < 2, [mp.nstr(v, 10) for v in folge_eta])
lauf.pruefe("Eulers wachsen: e < e₄ < e₅ < 4", folge_e[0] < folge_e[1] < folge_e[2] < 4, [mp.nstr(v, 10) for v in folge_e])

lauf.daten = {"eta5": mp.nstr(eta5, 30), "e5": mp.nstr(e5, 30), "P2_halbe": mp.nstr(a5[2], 15), "L_eta5": mp.nstr(s1["L"], 25),
              "lam_eta5": mp.nstr(s1["lam"], 20), "eulers": [mp.nstr(v, 20) for v in folge_e], "etas": [mp.nstr(v, 20) for v in folge_eta],
              "schnitte": {"1.70": [z_unter[0], z_unter[1]], "1.74": z_ueber[0]}}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
gx = [mp.mpf("-3") + mp.mpf(i) / 20 for i in range(0, 141)]
for b, farbe, name, L0 in ((mp.mpf("1.70"), "#2f8f83", "b = 1.70", mp.mpf("-1.66")), (eta5, "#b5523b", "b = η₅", s1["L"]), (mp.mpf("1.74"), "#7a49c2", "b = 1.74", mp.mpf("-1.67"))):
    s = stufe(b, e5, L0)
    ax[0].plot([float(x) for x in gx], [float(y) for y in s["P"](gx)], color=farbe, label=name)
ax[0].plot([-3, 4], [-3, 4], ":", color="#8a93a3", label="y = x")
ax[0].plot([float(e5)], [float(e5)], "o", color="#b5523b", ms=4)
ax[0].set_ylim(-2.5, 5); ax[0].set_xlim(-3, 4); ax[0].set_xlabel("x"); ax[0].set_title("Pentation b↑↑↑x berührt y = x bei e₅"); ax[0].legend(fontsize=8)
n = [3, 4, 5]
ax[1].plot(n, [float(v) for v in folge_e], "o-", color="#b5523b", label="„Eulers“ e, e₄, e₅")
ax[1].plot(n, [float(v) for v in folge_eta], "o-", color="#3d6fb6", label="„Etas“ η, η₄, η₅")
ax[1].axhline(4, ls=":", color="#b5523b"); ax[1].axhline(2, ls=":", color="#3d6fb6")
ax[1].set_xticks(n); ax[1].set_xticklabels(["Stufe 3\n(Potenz)", "Stufe 4\n(Tetration)", "Stufe 5\n(Pentation)"])
ax[1].set_ylim(1, 4.3); ax[1].set_title("TPID 10: gegen 4 und 2?"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"η₅ = {mp.nstr(eta5, 20)}, e₅ = {mp.nstr(e5, 20)} (Pentation regulär am negativen Fixpunkt von T_b); "
                        f"Etas {', '.join(mp.nstr(v, 6) for v in folge_eta)}; Eulers {', '.join(mp.nstr(v, 6) for v in folge_e)}.")
lauf.schreibe()
