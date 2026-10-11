"""exp.polynom_schranke — die Schranke für die Polynom-Wurzeln p_N: streng für endlich viele N, und warum der einfache Beweis für alle N scheitert.

 1. Krawczyk-Test (arb, Intervallrechnung) für N = 1 … 120: In der Box x̃ ± 2^−250·|L|^−k um die Näherung liegt genau eine Lösung von
    p∘p = exp + O(x^{N+1}); dazu die strenge Schranke max_{|x|≤|L|} |p_N(x)| ≤ Σ_k |a_k| |L|^k.
 2. Newton–Kantorovich um Knesers abgeschnittene Reihe q_N, Normen X_ρ (ρ = 1.3) und Y_s (s = 0.5):
    η_N = |J⁻¹F(q_N)|_X sagt den wahren Abstand |p_N − q_N|_X voraus und fällt; β_N = |J⁻¹|_{Y→X} wächst schneller, h_N = β γ η ≫ 1/2.
 3. Spektrum der gewichteten Linearisierung: eine Leiter kleiner Singulärwerte; der kleinste Singulärvektor sitzt bei arg x ≈ ±arg L.
"""
import mpmath as mp
from flint import arb, arb_mat, arb_poly, ctx

from gemeinsam import Lauf, kneser, plt
from polynom_wurzel import loesungen

mp.mp.dps = 180
ctx.prec = 600
lauf = Lauf("exp.polynom_schranke", __file__)
NMAX = 120
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
L_OBEN = arb("1.37455701074370750324")            # ≥ |L| = 1.374557010743707503228…


def F_J(a, N):
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    g, fpf = arb_poly([0]), arb_poly([0])
    for m in range(N + 1):
        g += a[m] * P[m]
        if m:
            fpf += m * a[m] * P[m - 1]
    F, J, fak = arb_mat(N + 1, 1), arb_mat(N + 1, N + 1), arb(1)
    for k in range(N + 1):
        fak = fak * k if k else fak
        F[k, 0] = g[k] - 1 / fak
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def mitte(M):
    return arb_mat([[arb(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


loes = loesungen(NMAX, bits=600)

# 1. Krawczyk
ergebnis, schranke = {}, {}
for N in range(1, NMAX + 1):
    xm = [arb(mp.nstr(c, 170)) for c in loes[N]]
    xm = [arb(c.mid()) for c in xm]
    r = [arb(2) ** (-250) * L_OBEN ** (-k) for k in range(N + 1)]
    X = [arb(xm[k].mid(), r[k].upper()) for k in range(N + 1)]
    Fm, Jm = F_J(xm, N)
    _, JX = F_J(X, N)
    C = mitte(mitte(Jm).inv())
    d = arb_mat([[X[k] - xm[k]] for k in range(N + 1)])
    eins = arb_mat([[1 if i == j else 0 for j in range(N + 1)] for i in range(N + 1)])
    K = arb_mat([[xm[k]] for k in range(N + 1)]) - C * Fm + (eins - C * JX) * d
    ergebnis[N] = all(abs(K[k, 0] - xm[k]).upper() < r[k].upper() for k in range(N + 1))
    schranke[N] = sum((abs(X[k]) * L_OBEN ** k for k in range(N + 1)), arb(0)).upper()
lauf.pruefe("Krawczyk für N = 1 … 120: in der Box um jede Näherung liegt genau eine Lösung (streng, Intervallrechnung)",
            all(ergebnis.values()), {"fehlgeschlagen": [N for N, ok in ergebnis.items() if not ok]})
B = max(schranke.values())
lauf.pruefe("streng: max_{|x|≤|L|} |p_N(x)| ≤ Σ|a_k||L|^k ≤ 2.2405 für alle N ≤ 120 (gemessenes Maximum 2.2327 bei x = |L|)",
            B < arb("2.2405"), {N: schranke[N].str(10, radius=False) for N in (1, 2, 3, 5, 10, 40, 120)})

# 2. Kantorovich um Knesers abgeschnittene Reihe
r_c, M, K_ = mp.mpf("0.9"), 256, 80
zs = [r_c * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
fz = kneser("sexp", "exp(1)", [w + mp.mpf("0.5") for w in kneser("slog", "exp(1)", zs)])
kn = [mp.re(mp.fsum(fz[j] * mp.expjpi(-2 * k * mp.mpf(j) / M) for j in range(M)) / M) / r_c ** k for k in range(K_ + 1)]
s, rho = mp.mpf("0.5"), mp.mpf("1.3")
zahl = lambda x: mp.mpf(x.mid().str(60, radius=False))
kant = {}
for N in (10, 20, 30, 40, 50, 60):
    F, J = F_J([arb(mp.nstr(c, 60)) for c in kn[:N + 1]], N)
    Ji = J.inv()
    Jinv = [[zahl(Ji[i, j]) for j in range(N + 1)] for i in range(N + 1)]
    Fv = [zahl(F[k, 0]) for k in range(N + 1)]
    beta = max(mp.fsum(abs(Jinv[k][j]) * rho ** k for k in range(N + 1)) / s ** j for j in range(N + 1))
    eta = mp.fsum(abs(mp.fsum(Jinv[k][j] * Fv[j] for j in range(N + 1))) * rho ** k for k in range(N + 1))
    wahr = mp.fsum(abs(loes[N][k] - kn[k]) * rho ** k for k in range(N + 1))
    R = mp.fsum(abs(c) * s ** k for k, c in enumerate(kn[:N + 1]))
    gamma = 2 * max(k * R ** (k - 1) / rho ** k for k in range(1, N + 1)) + mp.fsum(k * (k - 1) * abs(c) * R ** max(k - 2, 0) for k, c in enumerate(kn[:N + 1]))
    kant[N] = {"res": mp.fsum(abs(x) * s ** k for k, x in enumerate(Fv)), "beta": beta, "eta": eta, "wahr": wahr, "h": beta * gamma * eta}
lauf.pruefe("Newton-Schritt von Knesers Reihe aus sagt den wahren Abstand voraus: η_N = |p_N − q_N|_X auf 1 % (N = 10 … 60), und η fällt",
            all(abs(v["eta"] / v["wahr"] - 1) < 0.01 for v in kant.values()) and kant[60]["eta"] < kant[10]["eta"] / 1000,
            {N: mp.nstr(v["eta"], 3) for N, v in kant.items()})
wachs = (kant[60]["beta"] / kant[20]["beta"]) ** (mp.mpf(1) / 40)
lauf.pruefe("…aber |J⁻¹| wächst (Faktor ≈ 2.6 pro Grad) und h_N = βγη liegt ab N = 20 weit über 1/2 — Kantorovich mit Normen scheitert",
            all(kant[N]["h"] > mp.mpf("0.5") for N in (20, 30, 40, 50, 60)) and kant[60]["h"] > kant[20]["h"],
            {"Wachstum β pro Grad": mp.nstr(wachs, 3), "h": {N: mp.nstr(v["h"], 3) for N, v in kant.items()}})

# 3. Spektrum der gewichteten Linearisierung bei N = 30
N = 30
_, J = F_J([arb(mp.nstr(c, 60)) for c in kn[:N + 1]], N)
mp.mp.dps = 80
W = mp.matrix(N + 1, N + 1)
for k in range(N + 1):
    for j in range(N + 1):
        W[k, j] = zahl(J[k, j]) * s ** k / rho ** j
U, S, V = mp.svd_r(W)
sv = sorted(S[i] for i in range(N + 1))
i0 = min(range(N + 1), key=lambda t: S[t])
v = [V[i0, j] / rho ** j for j in range(N + 1)]
phis = [2 * mp.pi * t / 360 for t in range(360)]
betr = [abs(mp.polyval(v[::-1], rho * mp.expj(p))) for p in phis]
spitze = float(phis[max(range(180), key=lambda t: betr[t])]) * 180 / float(mp.pi)
lauf.pruefe("Spektrum (N = 30): Leiter kleiner Singulärwerte (je Faktor ≈ 3.4); der kleinste Singulärvektor ist auf |x| = 1.3 bei arg x ≈ arg L = 76.6° am grössten",
            sv[1] / sv[0] > 2 and sv[4] / sv[3] > 2 and abs(spitze - 76.6) < 10,
            {"kleinste": [mp.nstr(x, 3) for x in sv[:6]], "Spitze bei Grad": round(spitze, 1)})

lauf.daten = {"krawczyk_ok_bis": max(N for N, ok in ergebnis.items() if ok), "schranke": {N: schranke[N].str(12, radius=False) for N in schranke},
              "kantorovich": {N: {k: mp.nstr(x, 4) for k, x in v.items()} for N, v in kant.items()},
              "singulaerwerte_N30": [mp.nstr(x, 4) for x in sv[:10]], "spitze_grad": spitze}
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ns = sorted(kant)
for key, lab in (("res", "|F(q_N)|_Y"), ("beta", "β_N = |J⁻¹|"), ("eta", "η_N = |J⁻¹F|"), ("h", "h_N = βγη")):
    ax[0].semilogy(ns, [float(kant[N][key]) for N in ns], "o-", label=lab)
ax[0].axhline(0.5, color="k", ls=":", lw=0.8)
ax[0].set_xlabel("Grad N"); ax[0].set_title("Newton–Kantorovich um Knesers Reihe"); ax[0].legend(fontsize=8)
ax[1].plot([float(p) * 180 / float(mp.pi) for p in phis], [float(b) for b in betr])
ax[1].axvline(76.6, color="k", ls=":", lw=0.8); ax[1].axvline(360 - 76.6, color="k", ls=":", lw=0.8)
ax[1].set_xlabel("arg x (Grad), |x| = 1.3"); ax[1].set_title("kleinster Singulärvektor als Polynom (N = 30)")
lauf.bild(fig)
lauf.zusammenfassung = (f"Streng (Krawczyk) für N ≤ {NMAX}: eindeutige Lösung und |p_N| ≤ 2.2405 auf |x| ≤ |L|. Für alle N scheitert Kantorovich: "
                        f"|J⁻¹| wächst um ≈ {mp.nstr(wachs, 3)} pro Grad, der Rest steht fast senkrecht auf den fast-singulären Richtungen bei L und L̄.")
lauf.schreibe()
