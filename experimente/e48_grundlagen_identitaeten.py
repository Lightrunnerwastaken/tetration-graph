"""exp.grundlagen_identitaeten — Identitäten, die jede Tetration erfüllt, an Kneser-Werten nachgerechnet.

 1. Tauschregel: T(x+1)^{T(y)} = T(y+1)^{T(x)} = a^{T(x)T(y)}; daraus a = T(x+1)^{1/T(x)} und T(x+2) = T(x+1)^{T(x+1)/T(x)}.
 2. Taylor-Leiter: aus den Taylor-Koeffizienten t_n von T um x folgen die um x+1 (Miller-Rekursion für exp einer Reihe):
      e_0 = T(x+1),  n·e_n = c·Σ_{k=1}^n k·t_k·e_{n−k}   (c = ln a)   ⇔   T^{(n)}(x+1) = T(x+1)·B_n(cT'(x), …, cT^{(n)}(x)).
 3. Stufenleiter: P = Pentation (P(x+1) = T(P(x)), P(0) = 1, reguläre Lösung am Fixpunkt von T, stufen.pentation):
      P(−1) = 0, P(−2) = −1, und ln ln a = ln T'(0) − ln T'(−1) = ln P'(0) − 2 ln P'(−1) + ln P'(−2).
 4. Frei und fest: T̃ = T∘(id + p), p 1-periodisch, p(0) = 0 — Verhältnisse zwischen Stellen im Abstand 1 bleiben, alles andere ändert sich.
 5. Exponentialität: χ(f) = lim (slog f(x) − slog x) zählt Exponentiale: x² → 0, x^{ln x} → 0, 2^x → 1, x^x → 1.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, kneser, cauchy_ableitungen
from stufen import pentation, taylor

mp.mp.dps = 45
lauf = Lauf("exp.grundlagen_identitaeten", __file__)
E = "exp(1)"
FOX_L = mp.mpf("-1.850354529027181418483437788")
daten = {}


def T(basis, xs):
    return [mp.re(v) for v in kneser("sexp", basis, xs)]


# ------------------------------------------------------------------ 1. Tauschregel
x, y = mp.mpf("0.3"), mp.mpf("-0.6")
Tx, Ty, Tx1, Ty1, Tx2 = T("2", [x, y, x + 1, y + 1, x + 2])
tausch = [Tx1 ** Ty, Ty1 ** Tx, mp.mpf(2) ** (Tx * Ty)]
ganz = [16 ** 2 == 4 ** 4 == 2 ** 8, 65536 ** 2 == 4 ** 16 == 2 ** 32, 65536 ** 4 == 16 ** 16 == 2 ** 64]   # (m,n) = (3,2), (4,2), (4,3)
lauf.pruefe("Tauschregel T(x+1)^T(y) = T(y+1)^T(x) = a^(T(x)T(y)) (Basis 2, x = 0.3, y = −0.6) und ganz: 16² = 4⁴, 65536² = 4¹⁶",
            max(abs(tausch[0] - tausch[2]), abs(tausch[1] - tausch[2])) < mp.mpf(10) ** -35 and all(ganz), mp.nstr(tausch[2], 25))
lauf.pruefe("Basis aus zwei Nachbarwerten: T(x+1)^(1/T(x)) = 2, und T(x+2) = T(x+1)^(T(x+1)/T(x)) ohne a",
            abs(Tx1 ** (1 / Tx) - 2) < mp.mpf(10) ** -35 and abs(Tx1 ** (Tx1 / Tx) - Tx2) < mp.mpf(10) ** -35,
            {"T(1.3)^(1/T(0.3))": mp.nstr(Tx1 ** (1 / Tx), 25), "T(2.3)": mp.nstr(Tx2, 20)})

# ------------------------------------------------------------------ 2. Taylor-Leiter
K = 16
c2 = mp.log(2)
t0 = taylor(2, mp.mpf(0), mp.mpf("0.5"), N=96, K=K)
t1 = taylor(2, mp.mpf(1), mp.mpf("0.5"), N=96, K=K)
e = [mp.exp(c2 * t0[0])]
for n in range(1, K + 1):
    e.append(c2 * mp.fsum(k * t0[k] * e[n - k] for k in range(1, n + 1)) / n)
rel = max(abs(e[n] - t1[n]) / abs(t1[n]) for n in range(K + 1))
lauf.pruefe(f"Taylor-Leiter: die {K + 1} Koeffizienten von T₂ um 1 folgen aus denen um 0 (Miller-Rekursion = Bell-Polynome)",
            rel < mp.mpf(10) ** -20, f"grösster relativer Abstand {mp.nstr(rel, 3)}")
bell2 = t1[2] * 2 - t1[0] * (c2 * 2 * t0[2] + (c2 * t0[1]) ** 2)       # T''(1) = T(1)(cT''(0) + c²T'(0)²), mit T'' = 2·t₂
lauf.pruefe("zweite Stufe ausgeschrieben: T''(x+1) = T(x+1)·(c·T''(x) + c²·T'(x)²)", abs(bell2) < mp.mpf(10) ** -25, mp.nstr(bell2, 3))
daten["taylor"] = {"um_0": [mp.nstr(v, 15) for v in t0[:8]], "um_1": [mp.nstr(v, 15) for v in t1[:8]], "rel_abstand": mp.nstr(rel, 3)}

# ------------------------------------------------------------------ 3. Stufenleiter
leiter = {}
for basis, x0, a in ((E, FOX_L, mp.e), ("2", mp.mpf("-1.7439"), mp.mpf(2))):
    L, lam, psi, cc, kette, P = pentation(basis, x0, mp.mpf("0.08"))
    w = P([mp.mpf(-2), mp.mpf(-1), mp.mpf(0), mp.mpf(1), mp.mpf(2)])
    dP = [cauchy_ableitungen(P, mp.mpf(-j), mp.mpf("0.1"), m_max=1, N=48)[1] for j in (0, 1, 2)]
    dT = [cauchy_ableitungen(lambda zs: kneser("sexp", basis, zs), mp.mpf(-j), mp.mpf("0.3"), m_max=1, N=128)[1] for j in (0, 1)]
    leiter[basis] = {"P(-2..2)": [mp.nstr(v, 20) for v in w], "P'(0,-1,-2)": [mp.nstr(v, 20) for v in dP], "T'(0,-1)": [mp.nstr(v, 20) for v in dT],
                     "lnln_a": mp.nstr(mp.log(mp.log(a)), 20),
                     "erste_Differenz_T": mp.nstr(mp.log(dT[0]) - mp.log(dT[1]), 20),
                     "zweite_Differenz_P": mp.nstr(mp.log(dP[0]) - 2 * mp.log(dP[1]) + mp.log(dP[2]), 20)}
    T_a = T(basis, [a])[0]
    lauf.pruefe(f"Basis {basis}: Werte der Pentation P(−2) = −1, P(−1) = 0, P(0) = 1, P(1) = a, P(2) = T(a)",
                abs(w[0] + 1) + abs(w[1]) + abs(w[2] - 1) + abs(w[3] - a) + abs(w[4] / T_a - 1) < mp.mpf(10) ** -24, leiter[basis]["P(-2..2)"])
    lauf.pruefe(f"Basis {basis}: ln ln a = 1. Differenz von ln T′ = 2. Differenz von ln P′ (an 0, −1, −2), P′(0)/P′(−1) = T′(0), P′(−1)/P′(−2) = T′(−1)",
                abs(mp.log(dP[0]) - 2 * mp.log(dP[1]) + mp.log(dP[2]) - mp.log(mp.log(a))) < mp.mpf(10) ** -25
                and abs(dP[0] / dP[1] - dT[0]) + abs(dP[1] / dP[2] - dT[1]) < mp.mpf(10) ** -25,
                {k: leiter[basis][k] for k in ("lnln_a", "erste_Differenz_T", "zweite_Differenz_P")})
daten["stufenleiter"] = leiter

# ------------------------------------------------------------------ 4. Frei und fest
eps = mp.mpf("0.05")


def p(z):
    return eps * (mp.sin(2 * mp.pi * z) + (1 - mp.cos(2 * mp.pi * z)) / 2)


def T_gestoert(zs):
    return kneser("sexp", E, [z + p(z) for z in zs])


def jet(f, x0):
    return cauchy_ableitungen(f, x0, mp.mpf("0.2"), m_max=2, N=96)


stellen = [mp.mpf(-1), mp.mpf(0), mp.mpf("0.5"), mp.mpf(1)]
J = {float(s): jet(lambda zs: kneser("sexp", E, zs), s) for s in stellen}
G = {float(s): jet(T_gestoert, s) for s in stellen}
fest = {"T'(0)/T'(-1)": (J[0][1] / J[-1][1], G[0][1] / G[-1][1]),
        "(T''(0)-T'(0)^2)/T''(-1)": ((J[0][2] - J[0][1] ** 2) / J[-1][2], (G[0][2] - G[0][1] ** 2) / G[-1][2]),
        "T'(1)/T'(0)": (J[1][1] / J[0][1], G[1][1] / G[0][1])}
frei = {"T(1/2)": (J[0.5][0], G[0.5][0]), "T'(1/2)/T'(0)": (J[0.5][1] / J[0][1], G[0.5][1] / G[0][1]), "T'(0)": (J[0][1], G[0][1])}
lauf.pruefe("gestörte Tetration T(x + p(x)) (p 1-periodisch, p(0) = 0): Verhältnisse über Abstand 1 bleiben gleich (1, 1, e)",
            all(abs(u - v) < mp.mpf(10) ** -25 for u, v in fest.values()) and abs(fest["T'(1)/T'(0)"][0] - mp.e) < mp.mpf(10) ** -25,
            {k: [mp.nstr(u, 15), mp.nstr(v, 15)] for k, (u, v) in fest.items()})
lauf.pruefe("…aber Wert und Steigung bei ½ und die Steigung bei 0 ändern sich",
            all(abs(u - v) > mp.mpf("1e-3") for u, v in frei.values()), {k: [mp.nstr(u, 12), mp.nstr(v, 12)] for k, (u, v) in frei.items()})
daten["frei_und_fest"] = {"fest": {k: [mp.nstr(u, 20), mp.nstr(v, 20)] for k, (u, v) in fest.items()},
                          "frei": {k: [mp.nstr(u, 15), mp.nstr(v, 15)] for k, (u, v) in frei.items()}}


# ------------------------------------------------------------------ 5. Exponentialität
def slog_e(yv):
    n = 0
    while yv > 1:
        yv, n = mp.log(yv), n + 1
    return n + mp.re(kneser("slog", E, [yv])[0])


# x = T(k+2), Y = ln ln x = T(k), Z = ln Y = T(k−1):  slog f(x) − slog x = slog(ln ln ln f(x)) − slog Z
def klein(t, eZ):
    return mp.log1p(t) if eZ < 10 ** 4 else mp.mpf(0)          # e^{−T(4)} liegt unter jeder Rechengenauigkeit


LL2 = mp.log(mp.log(2))
FUNKTIONEN = {"x²": (lambda Z: Z + mp.log1p(mp.log(2) * mp.exp(-Z)), 0), "x^(ln x)": (lambda Z: Z + mp.log(2), 0),
              "2^x": (lambda Z: mp.exp(Z) + klein(LL2 * mp.exp(-mp.exp(Z)) if Z < 10 else 0, mp.exp(Z)), 1),
              "x^x": (lambda Z: mp.exp(Z) + klein(mp.exp(Z - mp.exp(Z)) if Z < 10 else 0, mp.exp(Z)), 1)}
Zs = [mp.mpf(0)]
for _ in range(4):
    Zs.append(mp.exp(Zs[-1]))                                   # Z = T(−1..3), also x = T(2..6)
chi = {name: [slog_e(g(Z)) - slog_e(Z) for Z in Zs] for name, (g, _) in FUNKTIONEN.items()}
lauf.pruefe("Exponentialität: slog f(x) − slog x bei x = e↑↑2 … e↑↑6 strebt gegen 0 (x², x^ln x) bzw. 1 (2^x, x^x)",
            all(abs(v[-1] - ziel) < abs(v[0] - ziel) and abs(v[-1] - ziel) < mp.mpf("1e-5") for (name, v), (_, ziel) in
                zip(chi.items(), FUNKTIONEN.values())),
            {name: [mp.nstr(z, 6) for z in v] for name, v in chi.items()})
daten["exponentialitaet"] = {name: [mp.nstr(z, 15) for z in v] for name, v in chi.items()}

lauf.daten = daten
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ns = list(range(K + 1))
ax[0].semilogy(ns, [abs(float(v)) for v in t0], "o-", label="Koeffizienten um 0 (gemessen)")
ax[0].semilogy(ns, [abs(float(v)) for v in t1], "s", ms=8, mfc="none", label="um 1 (gemessen)")
ax[0].semilogy(ns, [abs(float(v)) for v in e], "x", color="k", label="um 1 (aus denen um 0 gerechnet)")
ax[0].set_xlabel("n"); ax[0].set_title("Taylor-Leiter, Basis 2: |t_n|"); ax[0].legend(fontsize=8)
for (name, v), (_, ziel) in zip(chi.items(), FUNKTIONEN.values()):
    ax[1].semilogy(range(2, 7), [max(abs(float(z - ziel)), 1e-40) for z in v], "o-", label=f"{name} (Ziel {ziel})")
ax[1].set_xlabel("k  (x = e↑↑k)"); ax[1].set_title("|slog f(x) − slog x − χ|"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = ("Tauschregel, Taylor-Leiter (Bell-Polynome), Stufenleiter ln ln a = m-te Differenz von ln H′, Frei-und-fest unter "
                        "T∘(id+p) und Exponentialität über slog an Kneser-Werten bestätigt.")
lauf.schreibe()
