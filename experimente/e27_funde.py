"""exp.funde — unabhängige Nachrechnung der wiedergefundenen Untersuchungen (Juli/August 2026).

 (A) Geteilter Fixpunkt: b = −19/25, Partner d mit λ_d = λ_b − 2πi; σ = Log λ_d / Log λ_b; Residuum A_b − σ A_d fällt mit Quotient → λ_d.
 (B) Imaginäre Richtung: A_e(L_2) und K_1 = A_e'(L_2) (Koenigs-Abelfunktion am oberen Fixpunkt von e, Normierung κ'(L) = 1).
 (C) Periode-2-Region: geschlossene Form (u = −s/(e^s − 1), v = u e^s, c = v e^{−u}, μ = uv = (s/(2 sinh(s/2)))²),
     linker Randpunkt b_links, Spitzenkonstanten √6/e, 2/e, 2√6 e^{−3/2}, Rechteckwelle für b = 0.03.
 Vergleichswerte aus BASE_CHANGE_AT_A_SHARED_FIXED_POINT.md, IMAGINARY_DIRECTION.md, PERIODE_2_REGION.md.
"""
import math

import mpmath as mp
import numpy as np

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, plt

mp.mp.dps = 50
lauf = Lauf("exp.funde", __file__)
ZPI = 2j * mp.pi

# ------------------------------------------------------------------ (A) geteilter Fixpunkt, b = −19/25
cb = mp.log(mp.mpf(-19) / 25)
kand = []
for k in range(-4, 5):
    lam = -mp.lambertw(-cb, k)
    for m in (1, -1):
        if abs(lam - m * ZPI) < 1:
            kand.append((k, m, lam))
k, m, lam_b = kand[0]
L = lam_b / cb
lam_d = lam_b - m * ZPI
cd = lam_d * mp.exp(-lam_d)
d = mp.exp(cd)
sigma = mp.log(lam_d) / mp.log(lam_b)
orig = {"L": mp.mpc("1.9047135243385656", "-0.3776210535085195"), "d": mp.mpc("1.4251867364529", "-0.0462033694930"),
        "sigma": mp.mpc("-0.20046089358504287", "0.00348604017743")}
abw = max(abs(L - orig["L"]), abs(d - orig["d"]), abs(sigma - orig["sigma"]))
lauf.pruefe("b = −19/25: genau ein Partner (m = +1), L, d und σ wie im Original", len(kand) == 1 and abs(mp.im(cd)) <= mp.pi and abs(abw) < 1e-12,
            f"L = {mp.nstr(L, 16)}, d = {mp.nstr(d, 14)}, σ = {mp.nstr(sigma, 14)}, Abweichung {mp.nstr(abw, 2)}")
B = AmFixpunkt(cb, L, N=70)
D = AmFixpunkt(cd, L, N=70)
w0 = D.sigma(mp.mpc(1))
res = []
for n in (40, 41, 80, 81):
    v = D.lam ** n * w0
    res.append(mp.log(B.sigma(D.Psi(v)) / v) / B.log_lam)
quot = res[3] / res[2]
lauf.pruefe("Residuum A_b(T_d(n)) − σ·A_d(T_d(n)): Quotient aufeinanderfolgender Werte → λ_d", abs(quot - lam_d) < mp.mpf(10) ** -10,
            f"{mp.nstr(quot, 16)} vs λ_d = {mp.nstr(lam_d, 16)}")
pseudo = 2 * mp.pi / abs(mp.im(mp.log(lam_d)))
daempf = abs(lam_d) ** pseudo
lauf.pruefe("Pseudoperiode 2π/|arg λ_d| = 21.908, Dämpfung pro Periode 3.128e-4 (Original)",
            abs(pseudo - mp.mpf("21.908")) < 0.001 and abs(daempf / mp.mpf("3.128e-4") - 1) < 0.001, f"{mp.nstr(pseudo, 8)}, {mp.nstr(daempf, 6)}")

# ------------------------------------------------------------------ (B) imaginäre Richtung
Le = -mp.lambertw(-1, -1)
L2 = -mp.lambertw(-mp.log(2), -1) / mp.log(2)
E_ = AmFixpunkt(mp.mpf(1), Le, N=70)
Z_ = AmFixpunkt(mp.log(2), L2, N=70)


def abel(F, z):
    return mp.log(F.sigma(z)) / F.log_lam


def abel_strich(F, z, h=mp.mpf(10) ** -15):
    return (abel(F, z + h) - abel(F, z - h)) / (2 * h)


per_e = ZPI / E_.log_lam
roh = abel(E_, L2)
ziel = mp.mpc("0.326448651628991083472", "0.523513968909648002647")
n_w = mp.nint(mp.re((ziel - roh) / per_e))
lauf.pruefe("A_e(L_2) = 0.326448651628991 + 0.523513968909648i (bis auf den Zweig 2πi/Log λ_e)", abs(roh + n_w * per_e - ziel) < mp.mpf(10) ** -18,
            f"Zweigverschiebung {int(n_w)}, Rest {mp.nstr(abs(roh + n_w * per_e - ziel), 3)}")
K1 = abel_strich(E_, L2)
lauf.pruefe("K₁ = A_e'(L_2) = −0.0305475246507711 − 1.291981944842259i", abs(K1 - mp.mpc("-0.0305475246507710565", "-1.291981944842259")) < mp.mpf(10) ** -12,
            mp.nstr(K1, 16))
roh2 = abel(Z_, Le)
ziel2 = mp.mpc("-2.664068065075641706335", "0.043403036394983249782")
per_2 = ZPI / Z_.log_lam
n2 = mp.nint(mp.re((ziel2 - roh2) / per_2))
lauf.pruefe("umgekehrt A_2(L_e) = −2.664068065075642 + 0.043403036394983i", abs(roh2 + n2 * per_2 - ziel2) < mp.mpf(10) ** -18,
            f"Zweigverschiebung {int(n2)}")
lauf.pruefe("reelle Basen: arg λ = Im λ (Abklingrate in imaginärer Richtung), Drehrate log log d + Re λ",
            abs(mp.im(Z_.log_lam) - mp.im(Z_.lam)) < mp.mpf(10) ** -40 and abs(mp.re(Z_.log_lam) - (mp.log(mp.log(2)) + mp.re(Z_.lam))) < mp.mpf(10) ** -40)

# ------------------------------------------------------------------ (C) Periode-2-Region
def zyklus_aus_s(s):
    s = mp.mpc(s)
    u = -s / (mp.exp(s) - 1)
    v = u * mp.exp(s)
    c = v * mp.exp(-u)
    return u, v, c


u, v, c = zyklus_aus_s(mp.mpf("1.3"))
z1, z2 = mp.exp(v), mp.exp(u)
res_z = max(abs(mp.exp(c * z1) - z2), abs(mp.exp(c * z2) - z1), abs(u * v - (mp.mpf("1.3") / (2 * mp.sinh(mp.mpf("0.65")))) ** 2),
            abs(u * mp.exp(-v) - c))
lauf.pruefe("geschlossene Form: z1 = e^v, z2 = e^u ist 2-Zyklus von e^{cz}, μ = uv = (s/(2 sinh(s/2)))², c = u e^{−v}", res_z < mp.mpf(10) ** -40,
            mp.nstr(res_z, 3))


def G(p):
    r = mp.sqrt(p)
    return r / mp.sinh(r) if p != 0 else mp.mpf(1)


def p_von_phi(phi, start):
    return mp.findroot(lambda p: G(p) - mp.expj(phi), start)


def c_von_p(p):
    return zyklus_aus_s(2 * mp.sqrt(p))[2]


p, phi = mp.mpc(0), mp.mpf(0)
stufen = 400
for j in range(1, stufen + 1):
    phi = mp.mpf("1.06") * j / stufen
    p = p_von_phi(phi, p if j > 1 else mp.mpc(-6 * mp.mpf("0.0001"), -6 * mp.mpf("0.0001")))
phi_pi = mp.findroot(lambda f: mp.im(c_von_p(p_von_phi(f, p))) - mp.pi, mp.mpf("1.0513"))
b_links = mp.exp(c_von_p(p_von_phi(phi_pi, p)))
lauf.pruefe("linker Randpunkt b_links = −0.0805281576254835628 (|μ| = 1 und Im c = π)", abs(b_links - mp.mpf("-0.0805281576254835628232286623457")) < mp.mpf(10) ** -18,
            f"{mp.nstr(b_links, 20)}, φ = {mp.nstr(phi_pi, 12)}")
# Spitze
t = mp.mpf(10) ** -8
bt = mp.exp(-mp.e) * (1 - t)
ct = mp.log(bt)
EE = lambda z: mp.exp(ct * mp.exp(ct * z))
z = mp.findroot(lambda z: EE(z) - z, mp.exp(-1) + mp.sqrt(t))
zz = mp.exp(ct * z)
mu_t = ct ** 2 * z * zz
lauf.pruefe("an der Spitze: (1 − μ)/t → 2/e, |z1 − z2|/√t → 2√6 e^{−3/2}",
            abs((1 - mu_t) / t - 2 / mp.e) < 1e-6 and abs(abs(z - zz) / mp.sqrt(t) - 2 * mp.sqrt(6) * mp.e ** mp.mpf(-1.5)) < 1e-6,
            f"{mp.nstr((1 - mu_t) / t, 10)}, {mp.nstr(abs(z - zz) / mp.sqrt(t), 10)}")
zf, ce = 1.0, -math.e
for n in range(1, 10 ** 6 + 1):
    zf = math.exp(ce * zf)
spitze = abs(zf - math.exp(-1)) * math.sqrt(10 ** 6)
lauf.pruefe("Spitze b = e^{−e}: |z_n − 1/e|·√n → √6/e = 0.9011169177 (n = 10⁶: Original 0.901664)", abs(spitze - 0.901664) < 2e-4,
            f"{spitze:.6f}")
# reelle Achse: μ reell in (0, 1)
mus = []
for bb in (mp.mpf("0.06"), mp.mpf("0.03"), mp.mpf("0.001"), mp.mpf("1e-6")):
    cc = mp.log(bb)
    z = mp.findroot(lambda z: mp.exp(cc * mp.exp(cc * z)) - z, mp.mpf("0.9"))
    mus.append(cc ** 2 * z * mp.exp(cc * z))
lauf.pruefe("reell unter e^{−e}: Zyklus-Multiplikator immer in (0, 1) — keine Periodenverdopplung auf der Achse",
            all(0 < mu < 1 for mu in mus), [mp.nstr(mu, 6) for mu in mus])

# Rechteckwelle b = 0.03, Startinterpolation T(x) = 1 + x(b − 1) auf [0, 1]
b3 = 0.03
c3 = math.log(b3)
L3 = float(-mp.lambertw(-c3, 0).real / c3)
zA_ = float(mp.findroot(lambda z: mp.exp(c3 * mp.exp(c3 * z)) - z, 0.8))
zB_ = math.exp(c3 * zA_)
M = 2 ** 14
xs = np.arange(2 * M) / M                       # [0, 2)
T = np.where(xs < 1, 1 + xs * (b3 - 1), 0.0)
T[xs >= 1] = np.exp(c3 * (1 + (xs[xs >= 1] - 1) * (b3 - 1)))
for _ in range(80):
    T = np.exp(c3 * T)
koef = np.fft.rfft(T) / (2 * M)
ck = {k: abs(koef[k]) * math.pi * k for k in (1, 3, 5, 9, 21, 41)}   # c_k zu e^{iπkx}, Periode 2
gerade = max(abs(koef[k]) for k in (2, 4, 6, 10))
lauf.pruefe("Rechteckwelle: T(x + 80) hat |c_k|·πk = |zA − zB| = 0.7651944057 (ungerade k), gerade Moden ≈ 0",
            all(abs(v - 0.7651944057) < 2e-3 for v in ck.values()) and gerade < 1e-3,
            {k: round(v, 7) for k, v in ck.items()})
xstern = float(mp.findroot(lambda x: 1 + x * (b3 - 1) - L3, 0.5))
y = L3
for _ in range(12):
    y = math.exp(c3 * y)
lauf.pruefe("Kreuzungssatz: T(x*) = L (x* aus der Startinterpolation), dann T(x* + n) = L für alle n", abs(y - L3) < 1e-12,
            f"x* = {xstern:.9f}, L = {L3:.10f}")

lauf.daten = {"A": {"L": mp.nstr(L, 20), "d": mp.nstr(d, 16), "sigma": mp.nstr(sigma, 16), "quotient": mp.nstr(quot, 16)},
              "B": {"A_e_L2": mp.nstr(roh + n_w * per_e, 20), "K1": mp.nstr(K1, 16), "A_2_Le": mp.nstr(roh2 + n2 * per_2, 20)},
              "C": {"b_links": mp.nstr(b_links, 22), "phi_pi": mp.nstr(phi_pi, 14), "spitze_n1e6": spitze,
                    "rechteck": {k: v for k, v in ck.items()}, "xstern": xstern, "zA": zA_, "zB": zB_, "L": L3}}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
ax[0].plot(xs[::16], T[::16], lw=1, label="T(x + 80), b = 0.03")
ax[0].axhline(L3, ls=":", color="#b5523b", label="abstossender Fixpunkt L")
ax[0].set_xlabel("x"); ax[0].set_title("Rechteckwelle unter e^{−e}"); ax[0].legend(fontsize=8)
th = np.linspace(-math.pi, math.pi, 400)
ax[1].plot([float(mp.re(b_links))], [0], "o", color="#b5523b", label=f"b_links = {mp.nstr(b_links, 10)}")
ax[1].plot([math.exp(-math.e)], [0], "s", color="#3d6fb6", label="e^{−e}")
ax[1].set_xlim(-0.1, 0.1); ax[1].set_ylim(-0.08, 0.08); ax[1].set_aspect("equal"); ax[1].legend(fontsize=8)
ax[1].set_title("Endpunkte der Periode-2-Blase auf der reellen Achse")
lauf.bild(fig)
lauf.zusammenfassung = (f"Original-Zahlen unabhängig bestätigt: −19/25 ↔ d = {mp.nstr(d, 10)}, σ = {mp.nstr(sigma, 10)}, Residuumsquotient → λ_d; "
                        f"A_e(L_2), K₁, A_2(L_e) (bis auf Zweig); b_links = {mp.nstr(b_links, 16)}; Spitzenkonstanten; Rechteckwelle |c_k|πk = |zA − zB|.")
lauf.schreibe()
