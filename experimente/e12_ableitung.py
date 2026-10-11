"""exp.ableitung — Ableitungsidentitäten an Kneser (Basis e) und regulär (√2); Wendepunkte über viele Basen.

Ableitungen von Kneser über Cauchy-Integrale (Kreise um reelle Punkte, Motorwerte, 48 Knoten).
Geprüft: Kettenregel, Produktformel (n = 2, 3), Rückwärtsform, T'(n) aus T'(0), Formel für T'',
(ln T)''(x+1) = ln a · T''(x). Wendepunkte: Vorzeichenwechsel von T'' auf (−2, x_max).
"""
import mpmath as mp

from gemeinsam import Lauf, cache_speichern, cauchy_ableitungen, kneser, plt
from regulaer import Regulaer

mp.mp.dps = 40
lauf = Lauf("exp.ableitung", __file__)


def kneser_T(basis):
    return lambda zs: kneser("sexp", basis, zs, dps=40)


def abl(basis, x, m=2):
    x = mp.mpf(x)
    r = min(mp.mpf("0.3"), (x + 2) / 3)
    return cauchy_ableitungen(kneser_T(basis), mp.mpf(x), r, m, 64)


# ------------------------------------------------------------------ Identitäten an Kneser e
c = mp.mpf(1)                     # ln e
D = {x: abl("exp(1)", x, 2) for x in ["-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2"]}
D = {mp.mpf(k): v for k, v in D.items()}
T = lambda x: D[mp.mpf(x)][0]
T1 = lambda x: D[mp.mpf(x)][1]
T2 = lambda x: D[mp.mpf(x)][2]
grenze = mp.mpf(10) ** -22
res = max(abs(T1(x + 1) - c * T(x + 1) * T1(x)) / abs(T1(x + 1)) for x in map(mp.mpf, ["-1.5", "-1", "-0.5", "0", "0.5", "1"]))
lauf.pruefe("Kettenregel T'(x+1) = ln a·T(x+1)·T'(x) (Kneser e, 6 Stellen x, relativ)", res < grenze, mp.nstr(res, 3))
x = mp.mpf("-1")
res = max(abs(T1(x + n) - c ** n * mp.fprod(T(x + k) for k in range(1, n + 1)) * T1(x)) for n in (2, 3))
lauf.pruefe("Produktformel T'(x+n) = (ln a)^n ∏ T(x+k) · T'(x), n = 2, 3", res < grenze * 1e6, mp.nstr(res, 3))
res = abs(T1(mp.mpf(-1)) - T1(mp.mpf(0)) / c)
lauf.pruefe("T'(−1) = T'(0)/ln a", res < grenze, mp.nstr(res, 3))
res = abs(T1(mp.mpf(2)) - c ** 2 * mp.e * mp.e ** mp.e * T1(mp.mpf(0)))
lauf.pruefe("T'(2) = (ln a)² · a · a^a · T'(0)", res < grenze * 1e3, mp.nstr(res, 3))
falsch = max(abs(T2(x + 1) - c * T(x + 1) * (c * T(x + 1) * T1(x) ** 2 + T2(x))) for x in map(mp.mpf, ["-1", "-0.5", "0", "0.5"]))
lauf.pruefe("erste eigene Fassung mit Faktor T(x+1)² ist FALSCH (Gegenprobe)", falsch > 1, mp.nstr(falsch, 3))
res = max(abs(T2(x + 1) - c * T(x + 1) * (c * T1(x) ** 2 + T2(x))) / abs(T2(x + 1)) for x in map(mp.mpf, ["-1", "-0.5", "0", "0.5"]))
lauf.pruefe("T''(x+1) = ln a·T(x+1)·[ln a·T'(x)² + T''(x)] (relativ)", res < grenze, mp.nstr(res, 3))
lnT2 = lambda x: (T2(x) * T(x) - T1(x) ** 2) / T(x) ** 2
res = max(abs(lnT2(x + 1) - c * T2(x)) for x in map(mp.mpf, ["-0.5", "0", "0.5"]))
lauf.pruefe("(ln T)''(x+1) = ln a · T''(x)", res < grenze * 1e3, mp.nstr(res, 3))
taylor = cauchy_ableitungen(kneser_T("exp(1)"), mp.mpf(0), mp.mpf("0.5"), 8, 64)
koef = [t / mp.factorial(k) for k, t in enumerate(taylor)]
lauf.pruefe("Taylor bei 0: Koeffizient 0 ist T(0) = 1", abs(koef[0] - 1) < grenze, mp.nstr(koef[0], 12))

# ------------------------------------------------------------------ Regulär √2: dieselben Formeln, exakt
R = Regulaer(mp.sqrt(2), N=70, tol=mp.mpf("1e-3"))
lb = R.lb


def reg_abl(x):
    y = R.tet(x)
    a1 = mp.diff(R.slog, y)
    a2 = mp.diff(R.slog, y, 2)
    return y, 1 / a1, -a2 / a1 ** 3


res = 0
for x in map(mp.mpf, ["-1.3", "-0.4", "0.2", "0.9"]):
    y0, d0, _ = reg_abl(x)
    y1, d1, _ = reg_abl(x + 1)
    res = max(res, abs(d1 - lb * y1 * d0))
lauf.pruefe("Kettenregel an der regulären √2-Tetration", res < mp.mpf(10) ** -20, mp.nstr(res, 3))


# ------------------------------------------------------------------ Wendepunkte
def wendepunkte(f2, x0, x1, schritt, T0=None):
    xs, w, alt, out = [], [], None, []
    x = mp.mpf(x0)
    while x <= x1:
        if T0 is not None and abs(T0(x + mp.mpf("0.6"))) > 1e5:      # Motor nicht in astronomische Höhen treiben
            break
        v = f2(x)
        xs.append(x); w.append(v)
        if alt is not None and (alt[1] < 0) != (v < 0):
            lo, hi = alt[0], x
            for _ in range(40):
                m = (lo + hi) / 2
                if (f2(m) < 0) == (alt[1] < 0): lo = m
                else: hi = m
            out.append((lo + hi) / 2)
        alt = (x, v)
        x += schritt
    return out, xs, w


basen = [("1.5", "1.5", 18), ("2", "2", 4), ("e", "exp(1)", 3), ("3", "3", 2.5), ("5", "5", 2.2), ("10", "10", 2), ("100", "100", 1.6)]
tabelle, kurven = [], {}
for name, ausdruck, xmax in basen:
    f2 = lambda x, a=ausdruck: abl(a, x, 2)[2]
    T0 = lambda x, a=ausdruck: mp.re(kneser("sexp", a, [mp.mpf(x)], dps=40)[0])
    wp, xs, w = wendepunkte(f2, "-1.94", xmax, mp.mpf("0.06") if xmax < 5 else mp.mpf("0.2"), T0)
    xmax = float(xs[-1]) if xs else xmax
    tabelle.append({"basis": name, "wendepunkte": [mp.nstr(p, 10) for p in wp], "bereich": [-1.94, xmax]})
    kurven[name] = (xs, w)
    lauf.pruefe(f"Kneser b={name}: genau ein Wendepunkt auf (−1.94, {xmax})", len(wp) == 1, [mp.nstr(p, 8) for p in wp])
    cache_speichern()

reg_tab = []
for b in [mp.mpf("1.2"), mp.mpf("1.3"), mp.sqrt(2), mp.mpf("1.44")]:
    Rb = Regulaer(b, N=50, tol=mp.mpf("1e-3"))

    def f2(x, Rb=Rb):
        y = Rb.tet(x)
        a1 = mp.diff(Rb.slog, y); a2 = mp.diff(Rb.slog, y, 2)
        return -a2 / a1 ** 3
    mp.mp.dps = 30
    wp, _, _ = wendepunkte(f2, "-1.9", 8, mp.mpf("0.25"))
    mp.mp.dps = 40
    reg_tab.append({"basis": mp.nstr(b, 6), "wendepunkte": [mp.nstr(p, 8) for p in wp]})
    lauf.pruefe(f"regulär b={mp.nstr(b, 5)}: kein Wendepunkt (überall konkav) auf (−1.9, 8)", not wp, [mp.nstr(p, 6) for p in wp])

lauf.daten = {"T_e_strich_0": mp.nstr(D[mp.mpf(0)][1], 25), "T_e_2strich_0": mp.nstr(D[mp.mpf(0)][2], 25),
              "taylor_e_bei_0": [mp.nstr(k, 18) for k in koef], "wendepunkte_kneser": tabelle, "regulaer": reg_tab}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
xs = [mp.mpf(-1.9) + i * mp.mpf("0.05") for i in range(76)]
vals = [abl("exp(1)", x, 2) for x in xs]
ax[0].plot([float(x) for x in xs], [float(v[0]) for v in vals], label="T_e")
ax[0].plot([float(x) for x in xs], [float(v[1]) for v in vals], label="T_e′")
ax[0].plot([float(x) for x in xs], [float(v[2]) for v in vals], label="T_e″")
wp_e = [mp.mpf(p) for p in tabelle[2]["wendepunkte"]]
for p in wp_e:
    ax[0].axvline(float(p), ls="--", color="#888")
    ax[0].annotate(f"Wendepunkt {float(p):.4f}", (float(p), 3.2), textcoords="offset points", xytext=(4, 0), fontsize=8)
ax[0].set_ylim(-4, 8)
ax[0].axhline(0, color="#bbb", lw=0.8)
ax[0].set_xlabel("x")
ax[0].set_title("Kneser-Tetration zur Basis e und ihre Ableitungen")
ax[0].legend(fontsize=8)
bx = [z["basis"] for z in tabelle]
by = [float(mp.mpf(z["wendepunkte"][0])) if z["wendepunkte"] else float("nan") for z in tabelle]
ax[1].plot(range(len(bx)), by, "o-")
ax[1].set_xticks(range(len(bx)), bx)
ax[1].set_xlabel("Basis b")
ax[1].set_ylabel("Lage des Wendepunkts x_w")
ax[1].set_title("Genau ein Wendepunkt je Basis (Kneser)")
lauf.bild(fig)
cache_speichern()
lauf.zusammenfassung = (f"Alle Ableitungsidentitäten an Kneser e auf ≤ 1e-20 bestätigt. T_e'(0) = {mp.nstr(D[mp.mpf(0)][1], 16)}, "
                        f"T_e''(0) = {mp.nstr(D[mp.mpf(0)][2], 14)}. Wendepunkte (Kneser): "
                        + ", ".join(f"b={z['basis']}: {z['wendepunkte']}" for z in tabelle)
                        + ". Regulär (1<b<η): " + ", ".join(f"b={z['basis']}: {z['wendepunkte'] or 'keiner'}" for z in reg_tab) + ".")
lauf.schreibe()
