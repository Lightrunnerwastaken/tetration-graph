"""exp.carleman — intuitiver slog (abgeschnittenes Carleman-System) gegen Kneser-Referenzwerte.

slog_n(x) = −1 + Σ_{k=1}^n a_k x^k mit Σ_k a_k k^j/j! − a_j = δ_{j0} (j = 0..n−1).
Prüfpunkte: Kneser-Werte T_e(x) aus gp-tetration (80 Stellen, |T_e(x)| < |L_e| ≈ 1.3746).
Gegenbeispielkandidat wäre ein Fehler, der über mehrere Verdopplungen von n stagniert.
"""
import mpmath as mp

from gemeinsam import Lauf, plt, referenzwerte

lauf = Lauf("exp.carleman", __file__)
ref = referenzwerte()
W = ref["werte"]
punkte = []  # (y, erwartetes slog(y))
for x in ("-0.4375", "-0.3125", "-0.1875", "-0.0625"):
    punkte.append((W[f"sexp|e|{x}|80"], x))
punkte.append(("1.25", W["slog|e|1.25|80"]))

NS = [4, 8, 12, 16, 24, 32, 48, 64, 96, 128]
fehler = {}
for n in NS:
    mp.mp.dps = 80 + 5 * n
    M = mp.matrix(n, n)
    for j in range(n):
        fj = mp.factorial(j)
        for k in range(1, n + 1):
            M[j, k - 1] = mp.mpf(k) ** j / fj - (1 if (j >= 1 and j == k) else 0)
    rhs = mp.matrix(n, 1)
    rhs[0] = 1
    a = mp.lu_solve(M, rhs)
    err = []
    for y, s in punkte:
        y, s = mp.mpf(y), mp.mpf(s)
        val = -1 + mp.fsum(a[k - 1] * y ** k for k in range(1, n + 1))
        err.append(abs(val - s))
    fehler[n] = [mp.nstr(e, 4) for e in err]
    print(n, fehler[n])

mp.mp.dps = 30
inner = {n: max(mp.mpf(e) for e in fehler[n][:4]) for n in NS}
rand = {n: mp.mpf(fehler[n][4]) for n in NS}
lauf.pruefe("Fehler (innere Punkte) sinkt von n=16 bis n=128 um Faktor > 100", inner[128] < inner[16] / 100,
            f"{mp.nstr(inner[16], 3)} → {mp.nstr(inner[128], 3)}")
stagn = [n for a_, n in zip(NS, NS[1:]) if n >= 16 and inner[n] > inner[a_] * 0.9]
lauf.pruefe("keine Stagnation (Maximalfehler sinkt bei jedem Schritt ab n=16 um ≥ 10 %)", not stagn, stagn)
stellen = -mp.log10(inner[128])


def gerade(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    m = sxy / sxx
    res = sum((y - my - m * (x - mx)) ** 2 for x, y in zip(xs, ys))
    return m, res


ab = [n for n in NS if n >= 16]
p_alg, r_alg = gerade([mp.log(n) for n in ab], [mp.log(inner[n]) for n in ab])
p_geo, r_geo = gerade([mp.mpf(n) for n in ab], [mp.log(inner[n]) for n in ab])
lauf.pruefe("algebraischer Abfall n^p passt besser als geometrischer e^{−cn}", r_alg < r_geo,
            f"p = {mp.nstr(p_alg, 3)}; Residuen alg {mp.nstr(r_alg, 3)} vs geo {mp.nstr(r_geo, 3)}")
lauf.daten = {"n": NS, "fehler": fehler, "punkte": [p[1] if isinstance(p[1], str) else "slog(1.25)" for p in punkte],
              "exponent_p": float(p_alg), "stellen_n128": float(stellen),
              "vorab_erwartet": "≥ 10 Stellen bei n = 128 und Faktor 10^4 ab n = 16 — beides nicht eingetreten",
              "referenz_sha256": ref["sha256"][:16], "referenz_herkunft": ref["herkunft"]}

fig, ax = plt.subplots(figsize=(6.5, 4))
for i, (y, s) in enumerate(punkte):
    lab = f"slog(T_e({s}))" if i < 4 else "slog(1.25) (nahe |L|)"
    ax.semilogy(NS, [float(mp.mpf(fehler[n][i])) for n in NS], "o-", ms=4, label=lab)
ax.set_xlabel("Systemgrösse n")
ax.set_ylabel("|slog_n − slog_Kneser|")
ax.set_title("Intuitiver slog gegen Kneser (Basis e)")
ax.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Innere Punkte: Fehler {mp.nstr(inner[16], 2)} (n=16) → {mp.nstr(inner[128], 2)} (n=128), "
                        f"≈ {mp.nstr(stellen, 3)} Stellen; Abfall algebraisch ≈ n^{mp.nstr(p_alg, 3)}, nicht geometrisch. "
                        f"Punkt 1.25 nahe |L_e|: {mp.nstr(rand[128], 2)}. Konvergenz gestützt, aber langsam; "
                        "vorab erwartete 10 Stellen nicht erreicht.")
lauf.schreibe()
