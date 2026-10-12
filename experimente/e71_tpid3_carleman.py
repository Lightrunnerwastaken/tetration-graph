"""exp.tpid3 — TPID 3 (bo198214, 29.06.2008): Eigenwerte abgeschnittener Carleman-Matrizen von f(x) = (x+s)^p − s gehen gegen {p^k}.

 f hat den Fixpunkt x₀ = 1 − s mit Multiplikator f'(x₀) = p. Carleman-Matrix um 0: M[j][k] = Koeffizient von x^k in f(x)^j (j, k < N).
 Gegenstück zu TPID 1 (exp.vermutungen, Funktion b^x). Gemessen: grösster Abstand der ersten vier Eigenwerte (nach Nähe sortiert) zu 1, p, p², p³.
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 60
lauf = Lauf("exp.tpid3", __file__)


def reihe(s, p, N):
    """Taylorkoeffizienten von f(x) = (x+s)^p − s um 0 bis x^(N−1)."""
    a = [s ** p * mp.binomial(p, n) / s ** n for n in range(N)]
    a[0] -= s
    return a


def mal(a, b, N):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N)]


def carleman(s, p, N):
    f = reihe(s, p, N)
    zeilen, pot = [], [mp.mpf(1)] + [mp.mpf(0)] * (N - 1)
    for _ in range(N):
        zeilen.append(pot)
        pot = mal(pot, f, N)
    return mp.matrix(zeilen)


def abstand(s, p, N, m=4):
    ev = mp.eig(carleman(s, p, N), left=False, right=False)
    ziele = [mp.mpf(p) ** k for k in range(m)]
    frei, worst = list(ev), mp.mpf(0)
    for z in ziele:                      # gierige Zuordnung: nächster freier Eigenwert
        j = min(range(len(frei)), key=lambda i: abs(frei[i] - z))
        worst = max(worst, abs(frei[j] - z) / max(1, abs(z)))
        frei.pop(j)
    return worst


faelle = {"s=3/4, p=2/3 (anziehend, x0 im Kreis)": (mp.mpf("0.75"), mp.mpf(2) / 3),
          "s=3/4, p=3/2 (abstossend, x0 im Kreis)": (mp.mpf("0.75"), mp.mpf("1.5")),
          "s=3/10, p=2 (Polynom, x0 = 0.7)": (mp.mpf("0.3"), mp.mpf(2)),
          "s=1/2, p=1/2 (x0 auf dem Rand)": (mp.mpf("0.5"), mp.mpf("0.5")),
          "s=3/10, p=1/2 (x0 ausserhalb)": (mp.mpf("0.3"), mp.mpf("0.5"))}
Ns = [8, 12, 16, 24, 32]
werte = {}
for name, (s, p) in faelle.items():
    werte[name] = [abstand(s, p, N) for N in Ns]
    print(f"  {name}: " + ", ".join(f"N={N}: {mp.nstr(w, 3)}" for N, w in zip(Ns, werte[name])))

namen = list(faelle)
for name in namen[:2]:
    w = werte[name]
    lauf.pruefe(f"TPID 3, {name}: Abstand zu 1, p, p², p³ fällt und liegt bei N = 32 unter 10⁻³",
                all(w[i + 1] < w[i] for i in range(len(w) - 1)) and w[-1] < mp.mpf(10) ** -3, [mp.nstr(x, 3) for x in w])
ev = mp.eig(carleman(*faelle[namen[2]], 24), left=False, right=False)
gross = max(abs(e) for e in ev)
lauf.pruefe(f"Gegenbeispiel-Kandidat {namen[2]}: kein Eigenwert nahe p = 2, grösster Betrag ≈ 1 (Spektrum des superanziehenden Fixpunkts −s)",
            abs(gross - 1) < mp.mpf("0.01") and min(abs(e - 2) for e in ev) > mp.mpf("0.4"), f"max |λ| = {mp.nstr(gross, 6)}")
for name in namen[3:]:
    w = werte[name]
    lauf.pruefe(f"Grenzfall {name}: bis N = 32 keine Annäherung (Abstand bleibt > 0.1) — Fixpunkt 1 − s nicht im Konvergenzkreis |x| < s",
                min(w) > mp.mpf("0.1"), [mp.nstr(x, 3) for x in w])

lauf.daten = {"N": Ns, "abstand": werte}
lauf.zusammenfassung = "TPID 3: " + "; ".join(f"{k}: {mp.nstr(v[-1], 2)} bei N = {Ns[-1]}" for k, v in werte.items())
fig, ax = plt.subplots(figsize=(5, 3.4))
for name, w in werte.items():
    ax.semilogy(Ns, [float(x) for x in w], "o-", label=name)
ax.set_xlabel("Matrixgrösse N"); ax.set_ylabel("Abstand zu 1, p, p², p³"); ax.legend(fontsize=7); ax.set_title("TPID 3")
lauf.bild(fig)
lauf.schreibe()
