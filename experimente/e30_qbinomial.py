"""exp.qbinomial — q-Binomial-Reihe der regulären Tetration (Reshetnikov, MathOverflow 259278) und Möbius-Iteration.

 t(z) = Σ_n Σ_{k≤n} (−1)^{n−k} q^{C(n−k,2)} [z; n]_q [n; k]_q ^k a,   q = −W(−ln a) = λ (Multiplikator am unteren Fixpunkt).
 [z; n]_q = Π_{j<n} (1 − q^{z−j})/(1 − q^{j+1}).  Teilsummen = Interpolation in Potenzen von q^z durch die Punkte 0..N.
 Die Glieder sind von der Grösse q^{−n²/2}: man braucht etwa 0.08 n² Stellen Rechengenauigkeit (sonst täuscht Auslöschung Divergenz vor).
 Vergleich: reguläre Tetration Ψ(λ^z σ(1)) am unteren Fixpunkt (experimente/fixpunkt.py).
 Möbius: f(z) = (3z+2)/(z+2), Fixpunkte 2 (λ = 1/4) und −1 (λ = 4); reguläre Iterierte an beiden Fixpunkten stimmen überein.
"""
import mpmath as mp

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, plt
from regulaer import schroeder_reihe

lauf = Lauf("exp.qbinomial", __file__)
mp.mp.dps = 60


def referenz(a, zs, dps=60):
    with mp.workdps(dps):
        u = mp.log(a)
        L = mp.re(-mp.lambertw(-u) / u)
        F = AmFixpunkt(u, L, N=90)
        s1 = F.sigma(mp.mpf(1))
        return {z: F.Psi(F.lam ** z * s1) for z in zs}, F.lam


def qreihe(a, zs, N, dps):
    """Teilsummen der q-Binomial-Reihe für n = 0..N an den Stellen zs."""
    with mp.workdps(dps):
        a = mp.mpf(a)
        q = mp.re(-mp.lambertw(-mp.log(a)))
        tet = [mp.mpf(1)]
        for _ in range(N):
            tet.append(a ** tet[-1])
        # [n; k]_q über die Pascal-Regel [n; k] = [n−1; k−1] + q^k [n−1; k]
        qb = [[mp.mpf(1)]]
        for n in range(1, N + 1):
            zeile = [mp.mpf(1)] * (n + 1)
            for k in range(1, n):
                zeile[k] = qb[n - 1][k - 1] + q ** k * qb[n - 1][k]
            qb.append(zeile)
        teil = {z: [] for z in zs}
        for z in zs:
            s, zn = mp.mpf(0), mp.mpf(1)                     # zn = [z; n]_q
            for n in range(N + 1):
                if n > 0:
                    zn *= (1 - q ** (z - (n - 1))) / (1 - q ** n)
                d = mp.fsum((-1) ** (n - k) * q ** ((n - k) * (n - k - 1) // 2) * qb[n][k] * tet[k] for k in range(n + 1))
                s += zn * d
                teil[z].append(+s)
        return teil, q


ergebnisse = {}
for a, zs in ((mp.sqrt(2), [mp.mpf(1) / 2, mp.mpf("2.5"), mp.mpf("-0.7"), mp.mpc("0.5", "0.5")]), (mp.mpf("1.3"), [mp.mpf(1) / 2, mp.mpf("1.7")])):
    ref, lam = referenz(a, zs)
    N = 80 if a > mp.mpf("1.4") else 60                      # kleineres q braucht mehr Stellen je Glied
    stellen = int(N * N / 2 * mp.log10(1 / mp.re(-mp.lambertw(-mp.log(a))))) + 80   # Glieder ~ q^(−n²/2)
    teil, q = qreihe(a, zs, N, dps=stellen)
    fehler = {z: [abs(teil[z][n] - ref[z]) for n in range(N + 1)] for z in zs}
    fehler = {z: f + [f[-1]] * (80 - N) for z, f in fehler.items()}
    ergebnisse[mp.nstr(a, 6)] = (zs, ref, fehler, lam, q)
    print(f"  a = {mp.nstr(a, 6)}: q = {mp.nstr(q, 12)}, λ = {mp.nstr(lam, 12)}")
    for z in zs:
        print("    z =", mp.nstr(z, 4), " Fehler n=10/20/40/60/80:", [mp.nstr(fehler[z][n], 2) for n in (10, 20, 40, 60, 80)])

zs2, ref2, fe2, lam2, q2 = ergebnisse[mp.nstr(mp.sqrt(2), 6)]
lauf.pruefe("q = −W(−ln a) ist der Multiplikator λ am unteren Fixpunkt (a = √2: ln 2)", abs(q2 - lam2) < mp.mpf(10) ** -40 and abs(q2 - mp.log(2)) < mp.mpf(10) ** -40)
ganz, _ = qreihe(mp.sqrt(2), [mp.mpf(2), mp.mpf(3), mp.mpf(5)], 8, dps=80)
tet_ganz = {2: mp.sqrt(2) ** mp.sqrt(2), 3: mp.sqrt(2) ** (mp.sqrt(2) ** mp.sqrt(2))}
tet_ganz[5] = mp.sqrt(2) ** (mp.sqrt(2) ** tet_ganz[3])
lauf.pruefe("Teilsummen treffen die ganzen Höhen exakt (Interpolation: ab n > z bleibt die Summe stehen)",
            all(abs(ganz[mp.mpf(z)][n] - tet_ganz[z]) < mp.mpf(10) ** -50 for z in (2, 3, 5) for n in range(z, 9)))
lauf.pruefe("a = √2: Reihe gegen reguläre Tetration bei z = 1/2, 2.5, 0.5 + 0.5i (N = 80: unter 1e-30) und −0.7 (langsamer: unter 1e-15)",
            all(fe2[z][80] < mp.mpf(10) ** (-15 if mp.re(z) < 0 else -30) for z in zs2), {mp.nstr(z, 3): mp.nstr(fe2[z][80], 2) for z in zs2})
zs3, ref3, fe3, _, _ = ergebnisse[mp.nstr(mp.mpf("1.3"), 6)]
lauf.pruefe("a = 1.3: ebenso bei z = 1/2 und 1.7 (Fehler bei N = 60 unter 1e-30)", all(fe3[z][60] < mp.mpf(10) ** -30 for z in zs3),
            {mp.nstr(z, 3): mp.nstr(fe3[z][60], 2) for z in zs3})
# Rate: Fehler fällt geometrisch in N?
z = mp.mpf(1) / 2
raten = [(fe2[z][n + 10] / fe2[z][n]) ** (mp.mpf(1) / 10) for n in (30, 50, 65)]
lauf.pruefe("Fehler fällt geometrisch (Faktor je Glied stabil)", max(raten) / min(raten) < 1.5, [mp.nstr(r, 4) for r in raten])
# Faktor je Glied = |q^z| / q^{−2}: G(w) = T(z) mit w = q^z hat die nächste Singularität bei w = q^{−2} (z = −2)
rate_soll, rate_ist = {}, {}
for z, n in ((mp.mpf(1) / 2, 65), (mp.mpf("2.5"), 30), (mp.mpf("-0.7"), 65)):
    rate_ist[mp.nstr(z, 3)] = (fe2[z][n + 10] / fe2[z][n]) ** (mp.mpf(1) / 10)
    rate_soll[mp.nstr(z, 3)] = q2 ** (z + 2)
lauf.pruefe("Faktor je Glied ≈ q^{Re z + 2} (Abstand von w = q^z zur Singularität w = q^{−2}; auf 5 %)",
            all(abs(rate_ist[k] / rate_soll[k] - 1) < mp.mpf("0.05") for k in rate_ist),
            {k: f"{mp.nstr(rate_ist[k], 4)} gegen {mp.nstr(rate_soll[k], 4)}" for k in rate_ist})

# Reshetnikov-Vermutung 3 (T′ vollständig monoton) ⇔ alle Koenigs-Koeffizienten c_k = p_k σ(1)^k < 0
vorzeichen = {}
for a in (mp.sqrt(2), mp.mpf("1.3"), mp.mpf("1.1"), mp.mpf("1.44")):
    u = mp.log(a)
    Fk = AmFixpunkt(u, mp.re(-mp.lambertw(-u) / u), N=60)
    s1k = Fk.sigma(mp.mpf(1))
    ck = [mp.re(Fk.p[k] * s1k ** k) for k in range(1, 61)]
    vorzeichen[mp.nstr(a, 4)] = (all(c < 0 for c in ck), mp.nstr(ck[0], 4), mp.nstr(ck[59], 3))
lauf.pruefe("Koenigs-Koeffizienten c_1 … c_60 alle negativ (a = √2, 1.3, 1.1, 1.44) ⇒ T′ ist bis zu dieser Ordnung eine positive Summe fallender Exponentiale",
            all(v[0] for v in vorzeichen.values()), {k: f"c₁ = {v[1]}, c₆₀ = {v[2]}" for k, v in vorzeichen.items()})

# ------------------------------------------------------------------ Möbius
L1, L2 = mp.mpf(2), mp.mpf(-1)


def mf(z):
    return (3 * z + 2) / (z + 2)


lam1, lam2 = mp.mpf(1) / 4, mp.mpf(4)
NS = 40
g1 = [mp.mpf(0)] + [mp.taylor(lambda w: mf(L1 + w) - L1, 0, NS)[k] for k in range(1, NS + 1)]
s1 = schroeder_reihe(g1, lam1, NS)
geschlossen = mp.taylor(lambda w: (L1 - L2) * w / (L1 + w - L2), 0, NS)
lauf.pruefe("Möbius: allgemeine Schröder-Reihe am Fixpunkt 2 = geschlossene Form (L₁ − L₂)(z − L₁)/(z − L₂)",
            max(abs(s1[k] - geschlossen[k]) for k in range(NS + 1)) < mp.mpf(10) ** -30)


def halb(z, L, Lo, lam):
    s = (L - Lo) * (z - L) / (z - Lo)
    w = mp.sqrt(lam) * s
    return (L * (L - Lo) - Lo * w) / ((L - Lo) - w)                     # Umkehrung von s


testz = [mp.mpf("0.7"), mp.mpf("5"), mp.mpc("1", "1")]
diff = max(abs(halb(z, L1, L2, lam1) - halb(z, L2, L1, lam2)) for z in testz)
zweimal = max(abs(halb(halb(z, L1, L2, lam1), L1, L2, lam1) - mf(z)) for z in testz)
lauf.pruefe("Möbius: reguläre halbe Iterierte am Fixpunkt 2 = am Fixpunkt −1", diff < mp.mpf(10) ** -35, mp.nstr(diff, 3))
lauf.pruefe("Möbius: zweimal angewandt gibt f", zweimal < mp.mpf(10) ** -35, mp.nstr(zweimal, 3))

lauf.daten = {k: {"q": mp.nstr(v[4], 20), "fehler_N": {mp.nstr(z, 4): {n: mp.nstr(v[2][z][n], 3) for n in (10, 20, 40, 60, 80)} for z in v[0]}}
              for k, v in ergebnisse.items()}
lauf.daten["raten_sqrt2_z_halb"] = [mp.nstr(r, 6) for r in raten]
lauf.daten["rate_gegen_q_hoch_z_plus_2"] = {k: [mp.nstr(rate_ist[k], 6), mp.nstr(rate_soll[k], 6)] for k in rate_ist}
lauf.daten["koenigs_vorzeichen"] = {k: {"alle_negativ": v[0], "c1": v[1], "c60": v[2]} for k, v in vorzeichen.items()}
lauf.daten["T_sqrt2_halb"] = mp.nstr(mp.re(ref2[mp.mpf(1) / 2]), 40)

fig, ax = plt.subplots(figsize=(6.4, 4))
for z, farbe in zip(zs2, ("#b5523b", "#3d6fb6", "#2f8f83", "#7a49c2")):
    ax.semilogy(range(81), [float(max(fe2[z][n], mp.mpf(10) ** -60)) for n in range(81)], color=farbe, lw=1.2, label=f"√2, z = {mp.nstr(z, 3)}")
ax.set_xlabel("Glieder N"); ax.set_ylabel("Abstand zur regulären Tetration"); ax.set_title("q-Binomial-Reihe (Reshetnikov)"); ax.legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"q-Binomial-Reihe = reguläre Tetration (√2 und 1.3, reell und komplex), Fehler bei N = 80 unter 1e-30, geometrisch mit Faktor ≈ "
                        f"{mp.nstr(raten[-1], 3)} je Glied (= q^(z+2), Singularität bei z = −2); nötig sind etwa 0.08 N² Stellen. "
                        f"Koenigs-Koeffizienten c_1…c_60 negativ für vier Basen. Möbius: reguläre Iterierte an beiden Fixpunkten gleich.")
lauf.schreibe()
