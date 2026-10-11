"""exp.abel_bruchsumme — Bruchsumme des Superlogarithmus über η (Basis e, Kneser).

 A = slog_e. Auf [0, 1] als Tschebyschow-Interpolant aus Kneser-Werten, sonst A(y) = n + A(ln^n y).
 Kompensierte Bruchsumme (Müller–Schleicher, Summand mit A(N+x) − A(N) → 0):
   S(x) = lim_N [ x A(N+1) + Σ_{k=1}^N (A(k) − A(k+x)) ] = Σ_{k≥1} [ x ΔA(k) − (A(k+x) − A(k)) ].
 Behauptung: Glieder = (x − x²)/2 · A''(k) + …, |A''(y)| ≤ C/y² ⇒ absolut konvergent.
 Ableitungen: S''(x) = −Σ_{k≥1} A''(k+x);  Konstante γ_A := lim_N [Σ_{k≤N} A'(k) − A(N)] (Gegenstück zu γ für log).
"""
import mpmath as mp

from gemeinsam import Lauf, kneser, plt

mp.mp.dps = 30
lauf = Lauf("exp.abel_bruchsumme", __file__)
NK = 44
knoten = [mp.mpf(1) / 2 + mp.cos(mp.pi * (j + mp.mpf(1) / 2) / NK) / 2 for j in range(NK)]
werte = [mp.re(w) for w in kneser("slog", "exp(1)", knoten, dps=40)]
koef = []
for m in range(NK):
    koef.append(2 * mp.fsum(werte[j] * mp.cos(mp.pi * m * (j + mp.mpf(1) / 2) / NK) for j in range(NK)) / NK)
koef[0] /= 2


def A_basis(y):
    t = 2 * y - 1
    b1 = b2 = mp.mpf(0)
    for c in reversed(koef[1:]):
        b1, b2 = 2 * t * b1 - b2 + c, b1
    return t * b1 - b2 + koef[0]


def A(y):
    y, n = mp.mpf(y), 0
    while y > 1:
        y = mp.log(y)
        n += 1
    return n + A_basis(y)


# Ableitung der Tschebyschow-Reihe (d_{m-1} = d_{m+1} + 2m c_m), dann Kettenregel über die Logarithmen.
# Fassung 1 nahm mp.diff(A, y): bei y = 1 wechselt die Darstellung, der winzige Schritt machte A'(1) = 33067.
dkoef = [mp.mpf(0)] * (NK + 1)
for m in range(NK - 1, 0, -1):
    dkoef[m - 1] = dkoef[m + 1] + 2 * m * koef[m]
dkoef[0] /= 2


def A_basis_strich(y):
    t = 2 * y - 1
    b1 = b2 = mp.mpf(0)
    for c in reversed(dkoef[1:NK - 1]):
        b1, b2 = 2 * t * b1 - b2 + c, b1
    return 2 * (t * b1 - b2 + dkoef[0])


def A1(y):
    y, fak = mp.mpf(y), mp.mpf(1)
    while y > 1:
        fak /= y
        y = mp.log(y)
    return fak * A_basis_strich(y)


def A2(y):
    return mp.diff(A1, y, h=mp.mpf(10) ** -8 * (1 + abs(y)))


test = [mp.mpf("0.37"), mp.mpf("2.5"), mp.mpf("40")]
ref = [mp.re(w) for w in kneser("slog", "exp(1)", test, dps=40)]
fehler = max(abs(A(y) - r) for y, r in zip(test, ref))
lauf.pruefe("Interpolant + Rückführung A(y) = n + A(lnⁿ y) trifft Kneser-slog_e an freien Punkten (0.37, 2.5, 40)", fehler < mp.mpf(10) ** -20,
            mp.nstr(fehler, 3))
lauf.pruefe("Normierung A(1) = 0, A(e) = 1, A(0) = −1", max(abs(A(1)), abs(A(mp.e) - 1), abs(A_basis(0) + 1)) < mp.mpf(10) ** -20)
lauf.pruefe("abgeleitete Abel-Gleichung e·A'(e) = A'(1) über beide Darstellungen; A'(1) = 1/T'(0)", abs(A1(mp.e) * mp.e - A1(mp.mpf(1))) < mp.mpf(10) ** -20,
            mp.nstr(A1(mp.mpf(1)), 15))
schranke = {k: abs(A2(mp.mpf(k))) * k ** 2 for k in (10, 100, 1000, 10000)}
lauf.pruefe("Abschätzung |A''(y)| ≤ C/y² (y² |A''(y)| beschränkt und fallend)", all(schranke[k] < 1 for k in schranke) and schranke[10000] < schranke[10],
            {k: mp.nstr(v, 4) for k, v in schranke.items()})


def S(x, N=4000):
    x = mp.mpf(x)
    s = mp.fsum(x * (A(k + 1) - A(k)) - (A(k + x) - A(k)) for k in range(1, N + 1))
    m = N + mp.mpf("0.5")                                         # Schwanz: Glieder (x − x²)/2·A2 + (x − x³)/6·A3 + ...
    return s - (x - x * x) / 2 * A1(m) - (x - x ** 3) / 6 * A2(m)


xh = mp.mpf("0.5")
Sa, Sb = S(xh, 1000), S(xh, 4000)
lauf.pruefe("kompensierte Reihe konvergiert (N = 1000 vs 4000 mit Schwanzkorrektur)", abs(Sa - Sb) < mp.mpf(10) ** -9,
            f"S(1/2) = {mp.nstr(Sb, 15)}, Unterschied {mp.nstr(abs(Sa - Sb), 3)}")
S15 = S(mp.mpf("1.5"))
lauf.pruefe("Differenzengleichung S(3/2) − S(1/2) = A(3/2)", abs(S15 - Sb - A(mp.mpf("1.5"))) < mp.mpf(10) ** -9, mp.nstr(abs(S15 - Sb - A(mp.mpf("1.5"))), 3))
lauf.pruefe("ganze Stellen: S(1) = A(1) = 0, S(2) = A(1) + A(2)", abs(S(1)) < mp.mpf(10) ** -9 and abs(S(2) - A(2)) < mp.mpf(10) ** -9,
            f"S(2) = {mp.nstr(S(2), 12)}")
hd = mp.mpf("0.01")
S2_fd = (S(xh + hd, 2000) - 2 * S(xh, 2000) + S(xh - hd, 2000)) / hd ** 2
S2_reihe = -mp.fsum(A2(k + xh) for k in range(1, 2001)) - (-A1(2001 + xh - mp.mpf("0.5")))
lauf.pruefe("S''(x) = −Σ_k A''(k+x) (x = 1/2)", abs(S2_fd - S2_reihe) < mp.mpf(10) ** -4, f"{mp.nstr(S2_fd, 8)} vs {mp.nstr(S2_reihe, 8)}")
gam = {N: mp.fsum(A1(k) for k in range(1, N + 1)) - A(N) for N in (1000, 4000)}
gamma_A = gam[4000] - A1(mp.mpf(4000)) / 2                        # Euler–Maclaurin-Korrektur
lauf.pruefe("Konstante γ_A = lim [Σ_{k≤N} A'(k) − A(N)] existiert", abs(gam[1000] - A1(mp.mpf(1000)) / 2 - gamma_A) < mp.mpf(10) ** -8,
            f"γ_A = {mp.nstr(gamma_A, 12)}")
lauf.daten = {"S_halb": mp.nstr(Sb, 15), "S_2": mp.nstr(S(2), 15), "gamma_A": mp.nstr(gamma_A, 12),
              "y2_A2": {k: mp.nstr(v, 6) for k, v in schranke.items()}, "interpolationsfehler": mp.nstr(fehler, 3)}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
xs = [mp.mpf(i) / 10 for i in range(0, 31)]
ax[0].plot([float(x) for x in xs], [float(S(x, 600)) for x in xs], label="S(x) (kompensierte Bruchsumme)")
ax[0].plot([1, 2, 3], [0, float(A(2)), float(A(2) + A(3))], "o", color="#b5523b", label="Teilsummen A(1)+…+A(n)")
ax[0].set_xlabel("x"); ax[0].set_title("Bruchsumme von slog_e"); ax[0].legend(fontsize=8)
ks = [10 ** (mp.mpf(j) / 4) for j in range(4, 17)]
ax[1].loglog([float(k) for k in ks], [float(abs(A2(k)) * k ** 2) for k in ks], "o-", ms=3, label="y²·|A''(y)|")
ax[1].set_xlabel("y"); ax[1].set_title("A'' fällt schneller als 1/y²"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Kompensierte Bruchsumme von slog_e konvergiert (A'' = O(1/y²)); S(1/2) = {mp.nstr(Sb, 12)}, γ_A = {mp.nstr(gamma_A, 10)}; "
                        "Differenzengleichung und S'' = −Σ A'' bestätigt.")
lauf.schreibe()
