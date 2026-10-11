"""exp.pentation — Pentation als reguläre Superfunktion der Kneser-Tetration an ihrem reellen Fixpunkt.

 T_e hat auf (−2, −1) einen Fixpunkt L (T_e(L) = L) mit Multiplikator λ = T_e'(L) > 1 (abstossend für T, anziehend für slog).
 Poincaré-Funktion Ψ: Ψ(λw) = T_e(Ψ(w)). Pentation P(x) = e↑↑↑x = Ψ(c λ^x) mit c = σ(1), also P(0) = 1, P(x+1) = T_e(P(x)).
 c aus der slog-Kette 1 → 0 → −1 → … → L; P(x) für grössere x durch Vorwärtsiteration von T_e.
 Vergleich: J. D. Fox, Tetration Forum tid 103 (18.12.2007): L, 1/slog'(L) = λ, e↑↑↑n für n = 2 … −20;
 S. Levenstein, tid 1263 (10.04.2020): Basis 2, L ≈ −1.7439, 2↑↑↑1.5 ≈ 2.7029.
"""
import mpmath as mp

from gemeinsam import Lauf, kneser, plt
from stufen import pentation

mp.mp.dps = 50
lauf = Lauf("exp.pentation", __file__)

FOX = {"L": mp.mpf("-1.850354529027181418483437788"), "lam": mp.mpf("6.460671295681839390208370083"),
       2: mp.mpf("2075.968335058065833574141757"), -3: mp.mpf("-1.636358354286028979629049436"),
       -4: mp.mpf("-1.813170483098635639971748853"), -5: mp.mpf("-1.844484246898395061868430374"),
       -10: mp.mpf("-1.850354005584711078364293582")}


L, lam, psi, c, kette, P = pentation("exp(1)", FOX["L"], mp.mpf("0.08"))
T = lambda xs: [mp.re(v) for v in kneser("sexp", "exp(1)", xs)]
print(f"  L = {mp.nstr(L, 30)}, λ = {mp.nstr(lam, 30)}, c = {mp.nstr(c, 20)}")

lauf.pruefe("Fixpunkt: T_e(L) = L", abs(T([L])[0] - L) < mp.mpf(10) ** -30, mp.nstr(L, 30))
# Fox rechnete mit einer 1200-gliedrigen slog-Reihe (Matrixmethode); seine 28 Stellen tragen etwa 22 — Abstand ist das Mass der Übereinstimmung
lauf.pruefe("L stimmt mit Fox (2007) auf 22 Stellen überein", abs(L - FOX["L"]) < mp.mpf(10) ** -22, mp.nstr(L - FOX["L"], 3))
lauf.pruefe("Multiplikator λ = T_e'(L) = 6.4606712956818… (Fox: 1/slog'(L), 21 Stellen)", abs(lam - FOX["lam"]) < mp.mpf(10) ** -20,
            f"{mp.nstr(lam, 25)} (Abstand {mp.nstr(lam - FOX['lam'], 2)})")
lauf.pruefe("slog-Kette: slog(1) = 0, slog(0) = −1, dann e↑↑↑(−3), (−4), (−5) wie bei Fox (22 Stellen)",
            abs(kette[1]) + abs(kette[2] + 1) < mp.mpf(10) ** -30 and all(abs(kette[n] - FOX[-n]) < mp.mpf(10) ** -22 for n in (3, 4, 5)),
            {f"−{n}": f"{mp.nstr(kette[n], 22)} (Abstand {mp.nstr(kette[n] - FOX[-n], 2)})" for n in (3, 4, 5)})

xs = [mp.mpf(k) for k in (-10, -2, -1, 0, 1, 2)] + [mp.mpf("-0.5"), mp.mpf("0.5"), mp.mpf("1.5")]
w = dict(zip([mp.nstr(x, 3) for x in xs], P(xs)))
w2 = dict(zip([mp.nstr(x, 3) for x in xs], P(xs, M=10)))
lauf.pruefe("ganze Höhen aus der Ψ-Darstellung: P(−2) = −1, P(−1) = 0, P(0) = 1, P(1) = e",
            abs(w["-2.0"] + 1) + abs(w["-1.0"]) + abs(w["0.0"] - 1) + abs(w["1.0"] - mp.e) < mp.mpf(10) ** -24)
lauf.pruefe("P(2) = e↑↑e = 2075.968335058065833574… (Fox, relativ 1e-21) und P(−10) wie bei Fox",
            abs(w["2.0"] / FOX[2] - 1) < mp.mpf(10) ** -21 and abs(w["-10.0"] - FOX[-10]) < mp.mpf(10) ** -22,
            f"{mp.nstr(w['2.0'], 25)}; P(−10) − Fox = {mp.nstr(w['-10.0'] - FOX[-10], 2)}")
lauf.pruefe("Gleichung P(x+1) = T_e(P(x)) bei x = −0.5 und 0.5", abs(T([w["-0.5"]])[0] - w["0.5"]) + abs(T([w["0.5"]])[0] - w["1.5"]) < mp.mpf(10) ** -24)
lauf.pruefe("unabhängig von der Tiefe (M = 8 gegen 10)", max(abs(w[k] - w2[k]) for k in w) < mp.mpf(10) ** -22, mp.nstr(max(abs(w[k] - w2[k]) for k in w), 3))

# Basis 2 (Vergleich mit Levenstein 2020)
L2, lam2, psi2, c2, kette2, P2 = pentation("2", mp.mpf("-1.7439"), mp.mpf("0.08"))
p2 = P2([mp.mpf("1.5"), mp.mpf(1), mp.mpf(2)])
lauf.pruefe("Basis 2: L ≈ −1.7439 und 2↑↑↑1.5 ≈ 2.7029 (Levenstein 2020), 2↑↑↑1 = 2, 2↑↑↑2 = 4",
            abs(L2 + mp.mpf("1.7439")) < mp.mpf("5e-5") and abs(p2[0] - mp.mpf("2.7029")) < mp.mpf("5e-5")
            and abs(p2[1] - 2) + abs(p2[2] - 4) < mp.mpf(10) ** -24, f"L = {mp.nstr(L2, 15)}, 2↑↑↑1.5 = {mp.nstr(p2[0], 15)}")

lauf.daten = {"L": mp.nstr(L, 35), "lambda": mp.nstr(lam, 35), "c": mp.nstr(c, 30), "P": {k: mp.nstr(v, 30) for k, v in w.items()},
              "kette": [mp.nstr(z, 30) for z in kette], "psi": [mp.nstr(v, 15) for v in psi[:8]],
              "basis2": {"L": mp.nstr(L2, 30), "lambda": mp.nstr(lam2, 30), "P_1.5": mp.nstr(p2[0], 30)}}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
gx = [mp.mpf(-1.98) + mp.mpf(i) / 50 for i in range(0, 150)]
ax[0].plot([float(x) for x in gx], [float(v) for v in T(gx)], color="#3d6fb6", label="T_e(x) (Kneser)")
ax[0].plot([-2, 1], [-2, 1], ":", color="#8a93a3", label="y = x")
ax[0].plot([float(L)], [float(L)], "o", color="#b5523b", ms=5)
ax[0].annotate(f"L = {mp.nstr(L, 8)}\nλ = {mp.nstr(lam, 6)}", (float(L), float(L)), (-1.6, -1.6), fontsize=8,
               arrowprops={"arrowstyle": "->", "color": "#8a93a3"})
ax[0].set_xlim(-2, 1); ax[0].set_ylim(-3, 3); ax[0].set_xlabel("x"); ax[0].set_title("Fixpunkt der Tetration zur Basis e"); ax[0].legend(fontsize=8)
px = [mp.mpf(-6) + mp.mpf(i) / 20 for i in range(0, 152)]
ax[1].plot([float(x) for x in px], [float(v) for v in P(px)], color="#b5523b", label="e↑↑↑x")
px2 = [x for x in px if x < mp.mpf("1.95")]
ax[1].plot([float(x) for x in px2], [float(v) for v in P2(px2)], color="#7a49c2", label="2↑↑↑x")
ax[1].axhline(float(L), ls=":", color="#8a93a3")
ax[1].plot([k for k in range(-5, 2)], [float(v) for v in [FOX[-5], FOX[-4], FOX[-3], -1, 0, 1, mp.e]], "o", color="#222", ms=3, label="ganze Höhen (Fox 2007)")
ax[1].set_ylim(-2.2, 6); ax[1].set_xlabel("x"); ax[1].set_title("Pentation: P(x+1) = T(P(x))"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"L = {mp.nstr(L, 20)}, λ = {mp.nstr(lam, 15)}; e↑↑↑(1/2) = {mp.nstr(w['0.5'], 15)}, e↑↑↑(−1/2) = {mp.nstr(w['-0.5'], 15)}, "
                        f"e↑↑↑1.5 = {mp.nstr(w['1.5'], 15)}; Fox' Werte (2007) auf ≥ 24 Stellen bestätigt; Basis 2: 2↑↑↑1.5 = {mp.nstr(p2[0], 10)}.")
lauf.schreibe()
