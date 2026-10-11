"""exp.hexation_kneser — Kneser-artige Hexation über η₅ (Basis 2) und die Berührbasis η₆ der Hexation (TPID 10/11).

 Rechenwerk: knhexinit im gp-tetration-Fork (exp-079, 09.10.2026): reelle Taylorreihe der Abel-Funktion um Re q plus θ-Koeffizienten
 als ein Ausgleichsproblem, zweimal auf verschiedenen Bögen gelöst (Fehler-Vektor = erreichte Stellen). q ist der abstossende komplexe
 Fixpunkt der Pentation (kn.satz.kein_kneser_anziehend).
 η₆: grösste Basis, für die H_b einen reellen Fixpunkt hat. Für jede Basis Newton auf H_b'(x) = 1 (zentrale Differenzen),
 m(b) = H_b(x*) − x*, Sekanten in b auf m = 0; e₆ = x*(η₆).
"""
import sys

import mpmath as mp

from gemeinsam import GP_SRC, Lauf, plt

sys.path.insert(0, str(GP_SRC))
from fatou_backend import FatouGP, find_default_gp_exe  # noqa: E402

DPS = 30
mp.mp.dps = DPS + 15
lauf = Lauf("exp.hexation_kneser", __file__)
ETA = mp.e ** (1 / mp.e)
ETA4, E4 = mp.mpf("1.6353244967152763993453446183062"), mp.mpf("3.0885322718067176544821807826411")
ETA5, E5 = mp.mpf("1.737355925151687573101390"), mp.mpf("3.301407278523367493377305")
ETA6_KONTROLLE = mp.mpf("1.79912123612881268672356700702897599")      # Prototyp eta6.py bei dps 45 (Motor 43.7 Stellen)
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)
INIT = "hexinit(sexp(1),0,1)"


def s(x, n=DPS + 5):
    return mp.nstr(x, n)


def hexation(b):
    """hexinit zur Basis b; Rückgabe (H(xs) als Funktion, hexreg, Stellen). Überläufe per iferr → 10^100."""
    bs = s(b)
    kopf = gp.eval_batch(bs, [INIT, "hexreg", "hexerr"])
    H = lambda xs: [w.real for w in gp.eval_batch(bs, [f"iferr(hex({s(x)}),E,10^100)" for x in xs])]
    return H, int(kopf[1].real), kopf[2].real


def minimum(H, x0=None):
    """x* mit H'(x*) = 1 (Start: Gitterminimum von H(x) − x auf [2, 6]); Rückgabe m, x*, H', H''."""
    if x0 is None:
        gitter = [mp.mpf(2) + mp.mpf(k) / 8 for k in range(33)]
        g = [h - x for h, x in zip(H(gitter), gitter)]
        x0 = gitter[min(range(len(g)), key=lambda j: g[j])]
    x, h = x0, mp.mpf(10) ** -(DPS // 2 - 2)
    for _ in range(40):
        hm, h0, hp = H([x - h, x, x + h])
        d1, d2 = (hp - hm) / (2 * h), (hp - 2 * h0 + hm) / h ** 2
        dx = max(min((d1 - 1) / d2, mp.mpf("0.25")), mp.mpf("-0.25"))
        x -= dx
        if abs(dx) < mp.mpf(10) ** -(DPS // 2 - 4):
            break
    return H([x])[0] - x, x, d1, d2


try:
    # ------------------------------------------------------------------ A: Basis 2
    H2, reg2, stellen2 = hexation(2)
    lauf.pruefe("Basis 2 nimmt den Kneser-Weg (hexreg = 2) mit mindestens 30 Stellen (Fehler-Vektor zweier Bögen)", reg2 == 2 and stellen2 >= 30,
                f"{mp.nstr(stellen2, 4)} Stellen")
    werte2 = H2([mp.mpf(x) / 2 for x in range(-6, 5)])          # −3, −2.5, …, 2
    ganz = {x: werte2[2 * x + 6] for x in (-3, -2, -1, 0, 1, 2)}
    soll = {-3: -2, -2: -1, -1: 0, 0: 1, 1: 2, 2: 4}
    lauf.pruefe("ganze Höhen: H(−3) = −2, H(−2) = −1, H(−1) = 0, H(0) = 1, H(1) = 2, H(2) = 2↑↑↑2 = 4",
                max(abs(ganz[x] - soll[x]) for x in soll) < mp.mpf(10) ** -25, {str(x): mp.nstr(ganz[x], 12) for x in soll})
    p = gp.eval_batch("2", ["kp", "klam", "kncc[1]", "hex(1.5)", "hex(0.5)"])
    h15, h05 = p[3].real, p[4].real
    lauf.pruefe("2↑↑↑↑1.5 stimmt mit Levensteins 2.6729 (Tetration Forum tid 1263, 2020) auf 4 Stellen", abs(h15 - mp.mpf("2.6729")) < mp.mpf("5e-5"),
                f"2↑↑↑↑1.5 = {mp.nstr(h15, 30)}")
    lauf.pruefe("Fixpunkt q der Pentation abstossend, |P'(q)| > 1 (Kneser braucht das)", abs(p[1]) > 1,
                f"q = {mp.nstr(p[0], 15)}, P'(q) = {mp.nstr(p[1], 12)}, |P'(q)| = {mp.nstr(abs(p[1]), 8)}")
    halbe = H2([mp.mpf(x) / 4 for x in range(-11, 9)])
    lauf.pruefe("H steigt auf (−2.75, 2]", all(a < b for a, b in zip(halbe, halbe[1:])))
    periode = 2j * mp.pi / mp.log(p[1])
    K = -p[2] - periode                                          # ζ(z) − z für grosse Im z (Normierung von F_q wie in prototypen_hexation)
    zeta_c = mp.mpc("0.40312930", "-1.70978835")                 # Ecke der reellen Kurve (prototypen_hexation/hex_ecke.py)
    z0 = zeta_c - 2 + periode - K

    # ------------------------------------------------------------------ B: η₆ und e₆
    rand = {}
    for b in (mp.mpf("1.78"), mp.mpf("1.80")):
        Hb, reg, st = hexation(b)
        rand[b] = minimum(Hb) + (reg, st)
    lauf.pruefe("Einschluss: H_1.78 hat einen reellen Fixpunkt (min(H − x) < 0), H_1.80 keinen (> 0); beide Kneser",
                rand[mp.mpf("1.78")][0] < 0 < rand[mp.mpf("1.80")][0] and all(v[4] == 2 for v in rand.values()),
                {mp.nstr(b, 3): f"min = {mp.nstr(v[0], 8)} bei x = {mp.nstr(v[1], 8)}" for b, v in rand.items()})
    b0, b1 = mp.mpf("1.795"), mp.mpf("1.80")
    m0 = minimum(hexation(b0)[0])[0]
    m1, x1 = rand[b1][0], rand[b1][1]
    verlauf = []
    for _ in range(12):
        b2 = b1 - m1 * (b1 - b0) / (m1 - m0)
        Hb, reg, st = hexation(b2)
        m2, x2, d1, d2 = minimum(Hb, x1)
        verlauf.append([mp.nstr(b2, DPS), mp.nstr(m2, 3)])
        print(f"  b = {mp.nstr(b2, DPS)}  m = {mp.nstr(m2, 3)}  x* = {mp.nstr(x2, 20)}", flush=True)
        b0, m0, b1, m1, x1 = b1, m1, b2, m2, x2
        if abs(b1 - b0) < mp.mpf(10) ** -(DPS - 6):
            break
    eta6, e6 = b1, x1
    lauf.pruefe("Berührung: H(e₆) = e₆, H'(e₆) = 1, H''(e₆) > 0 bei b = η₆",
                abs(m1) < mp.mpf(10) ** -25 and abs(d1 - 1) < mp.mpf(10) ** -10 and d2 > 0,
                f"η₆ = {mp.nstr(eta6, 26)}, e₆ = {mp.nstr(e6, 18)}, H(e₆) − e₆ = {mp.nstr(m1, 3)}, H''/2 = {mp.nstr(d2 / 2, 8)}")
    lauf.pruefe("η₆ stimmt mit dem Kontrollwert bei dps 45 auf 25 Stellen", abs(eta6 - ETA6_KONTROLLE) < mp.mpf(10) ** -25,
                f"Δ = {mp.nstr(eta6 - ETA6_KONTROLLE, 3)}")
    etas, eulers = [ETA, ETA4, ETA5, eta6], [mp.e, E4, E5, e6]
    lauf.pruefe("TPID 11: η₆² < e₅", eta6 ** 2 < E5, f"η₆² = {mp.nstr(eta6 ** 2, 10)} < e₅ = {mp.nstr(E5, 10)}")
    lauf.pruefe("Etas wachsen: η < η₄ < η₅ < η₆ < 2", all(a < b for a, b in zip(etas, etas[1:] + [2])), [mp.nstr(v, 10) for v in etas])
    lauf.pruefe("Eulers wachsen: e < e₄ < e₅ < e₆ < 4", all(a < b for a, b in zip(eulers, eulers[1:] + [4])), [mp.nstr(v, 10) for v in eulers])
    exp_eta = [mp.log((2 - etas[i]) / (2 - etas[i + 1])) / mp.log(mp.mpf(i + 4) / (i + 3)) for i in range(3)]
    exp_e = [mp.log((4 - eulers[i]) / (4 - eulers[i + 1])) / mp.log(mp.mpf(i + 4) / (i + 3)) for i in range(3)]

    # ------------------------------------------------------------------ Bild
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    gx = [mp.mpf(2) + mp.mpf(k) / 20 for k in range(0, 61)]
    for b, farbe, name in ((mp.mpf("1.78"), "#2f8f83", "b = 1.78"), (eta6, "#b5523b", "b = η₆"), (mp.mpf("1.80"), "#7a49c2", "b = 1.80")):
        Hb = hexation(b)[0]
        ys = Hb(gx)
        paare = [(float(x), float(y - x)) for x, y in zip(gx, ys) if y < 50]
        ax[0].plot([a for a, _ in paare], [c for _, c in paare], color=farbe, label=name)
    ax[0].axhline(0, ls=":", color="#8a93a3")
    ax[0].plot([float(e6)], [0], "o", color="#b5523b", ms=4)
    ax[0].set_ylim(-0.9, 0.6); ax[0].set_xlim(2, 5); ax[0].set_xlabel("x"); ax[0].set_ylabel("H_b(x) − x")
    ax[0].set_title("Hexation b↑↑↑↑x berührt y = x bei e₆"); ax[0].legend(fontsize=8)
    n = [3, 4, 5, 6]
    ax[1].plot(n, [float(4 - v) for v in eulers], "o-", color="#b5523b", label="4 − e_n")
    ax[1].plot(n, [float(2 - v) for v in etas], "o-", color="#3d6fb6", label="2 − η_n")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xticks(n); ax[1].set_xticklabels(["3", "4", "5", "6"]); ax[1].minorticks_off()
    ax[1].set_xlabel("Stufe n (3 = Potenz, 6 = Hexation)"); ax[1].set_title("TPID 10: Abstände zu 4 und 2 (doppelt-log)"); ax[1].legend(fontsize=8)
    lauf.bild(fig)
finally:
    gp.close()

lauf.daten = {"hex2": {"stellen": mp.nstr(stellen2, 5), "H_1.5": mp.nstr(h15, 30), "H_0.5": mp.nstr(h05, 30), "q": mp.nstr(p[0], 25),
                       "P_strich_q": mp.nstr(p[1], 25), "theta_c0": mp.nstr(p[2], 20), "periode": mp.nstr(periode, 15),
                       "verzweigung_z0": mp.nstr(z0, 10)},
              "eta6": mp.nstr(eta6, 28), "e6": mp.nstr(e6, 18), "secanten": verlauf,
              "etas": [mp.nstr(v, 20) for v in etas], "eulers": [mp.nstr(v, 20) for v in eulers],
              "abstand_exponenten_eta": [mp.nstr(v, 5) for v in exp_eta], "abstand_exponenten_e": [mp.nstr(v, 5) for v in exp_e]}
lauf.zusammenfassung = (f"2↑↑↑↑1.5 = {mp.nstr(h15, 20)} (Kneser an q, {mp.nstr(stellen2, 3)} Stellen); η₆ = {mp.nstr(eta6, 20)}, e₆ = {mp.nstr(e6, 15)}; "
                        f"Etas {', '.join(mp.nstr(v, 6) for v in etas)}; Eulers {', '.join(mp.nstr(v, 6) for v in eulers)}.")
lauf.schreibe()
