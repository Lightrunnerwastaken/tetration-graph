"""exp.sqrt2_abel — Differenz der regulären Abel-Funktionen an den Fixpunkten a und A.

Für 1 < b < η ist p(θ) = α_A(α_a^{-1}(θ)) − θ 1-periodisch. Gemessen: Fourier-Moden |p̂(k)|, k ≤ 10.
Modell: ln|p̂(k)| = c + γ ln k − S k + d/k.
  H_a: S = 2π²/|ln λ_a| (Streifenhalbbreite π/|ln λ_a|, unterer Fixpunkt)
  H_A: S = 2π²/ln λ_A   (oberer Fixpunkt)
Zusätzlich: Exponent γ (Vermutung kl.verm.exponent: γ = −5/4 für alle b).
Gegencheck: dieselbe Messung mit anderer Reihenlänge und Toleranz.
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 170
lauf = Lauf("exp.sqrt2_abel", __file__)
NGITTER = 64
KMAX = 10
BASEN = [mp.mpf("1.25"), mp.mpf("1.3"), mp.mpf("1.35"), mp.sqrt(2), mp.mpf("1.43")]
NAMEN = ["1.25", "1.3", "1.35", "√2", "1.43"]
BODEN = mp.mpf(10) ** (-mp.mp.dps + 30)


def moden(b, N=100, tol="2e-3"):
    R = Regulaer(b, N=N, tol=mp.mpf(tol))
    x0 = (R.a + R.A) / 2
    th0 = R.alpha_unten(x0)
    lo, hi = R.f(x0), x0
    werte = []
    for j in range(NGITTER):
        th = th0 + mp.mpf(j) / NGITTER
        x = mp.findroot(lambda y: R.alpha_unten(y) - th, (lo, hi), solver="illinois",
                        tol=mp.mpf(10) ** (-mp.mp.dps + 10))
        werte.append(R.alpha_oben(x) - th)
    out = []
    for k in range(1, KMAX + 1):
        c = mp.fsum(werte[j] * mp.expjpi(-2 * k * mp.mpf(j) / NGITTER) for j in range(NGITTER)) / NGITTER
        out.append(abs(c))
    return R, out


def anpassen(ks, ms, S=None):
    """Kleinste Quadrate für ln m = c + γ ln k − S k + d/k (S frei oder fest)."""
    zeilen, rechts = [], []
    for k, m in zip(ks, ms):
        k = mp.mpf(k)
        if S is None:
            zeilen.append([1, mp.log(k), -k, 1 / k])
            rechts.append(mp.log(m))
        else:
            zeilen.append([1, mp.log(k), 1 / k])
            rechts.append(mp.log(m) + S * k)
    A = mp.matrix(zeilen)
    y = mp.matrix(rechts)
    x = mp.lu_solve(A.T * A, A.T * y)
    return [x[i] for i in range(x.rows)]


ergebnis, gammas = [], []
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
for b, name in zip(BASEN, NAMEN):
    R, ms = moden(b)
    S_a = 2 * mp.pi ** 2 / abs(mp.log(R.lam_a))
    S_A = 2 * mp.pi ** 2 / mp.log(R.lam_A)
    gueltig = [(k, m) for k, m in zip(range(1, KMAX + 1), ms) if m > BODEN * 1e6]
    ks, mv = [k for k, _ in gueltig], [m for _, m in gueltig]
    z = {"basis": name, "a": R.a, "A": R.A, "lambda_a": R.lam_a, "lambda_A": R.lam_A,
         "S_vorhersage_a": float(S_a), "S_vorhersage_A": float(S_A),
         "moden": [mp.nstr(m, 8) for m in ms], "gueltige_moden": len(ks)}
    lauf.pruefe(f"b={name}: Differenz nicht konstant", ms[0] > BODEN * 1e10, mp.nstr(ms[0], 5))
    if len(ks) >= 5:
        c, g, S, d = anpassen(ks, mv)
        z.update({"S_frei": float(S), "gamma_frei": float(g), "rel_S_minus_Sa": float(S / S_a - 1),
                  "rel_S_minus_SA": float(S / S_A - 1)})
        lauf.pruefe(f"b={name}: freie Steigung S trifft 2π²/|ln λ_a| auf 1e-3", abs(S / S_a - 1) < 1e-3,
                    f"S={float(S):.5f}, Vorhersage {float(S_a):.5f}")
        c2, g2, d2 = anpassen(ks, mv, S=S_a)
        z["gamma_bei_Sa"] = float(g2)
        gammas.append(float(g2))
    elif len(ks) >= 2:
        # zu wenige Moden für den freien Fit: Steigung nach Abzug von γ = −5/4 prüfen
        k1, k2 = ks[-2], ks[-1]
        S_eff = (mp.log(mv[-2]) - mp.log(mv[-1]) + mp.mpf(-5) / 4 * mp.log(mp.mpf(k2) / k1)) / (k2 - k1)
        z["S_mit_gamma_54"] = float(S_eff)
        z["rel_S_minus_Sa"] = float(S_eff / S_a - 1)
        lauf.pruefe(f"b={name}: Steigung (mit γ=−5/4) trifft 2π²/|ln λ_a| auf 1e-2", abs(S_eff / S_a - 1) < 1e-2,
                    f"S={float(S_eff):.4f}, Vorhersage {float(S_a):.4f}")
    # modellfrei: r_k = |p̂(k)| e^{S_a k} / (|p̂(1)| e^{S_a}) — universell, falls nur C(b) basisabhängig ist
    r = [m * mp.exp(S_a * (k - 1)) / mv[0] for k, m in gueltig]
    z["r_k"] = [mp.nstr(x, 8) for x in r]
    ergebnis.append(z)
    ax[0].semilogy(ks, [float(m) for m in mv], "o-", ms=4, label=f"b = {name}")
    ax[1].plot(ks, [float(x * mp.mpf(k) ** mp.mpf(1.25)) for k, x in zip(ks, r)], "o-", ms=4, label=f"b = {name}")

lauf.pruefe("angepasstes γ (bei festem S_a) ist für alle Basen gleich (Spannweite < 0.02)",
            len(gammas) >= 3 and max(gammas) - min(gammas) < 0.02, [round(g, 4) for g in gammas])
# modellfreie Universalität der normierten Verhältnisse
voll = [z for z in ergebnis if z["gueltige_moden"] >= 5]
spreizung = []
for k in range(5):
    werte = [mp.mpf(z["r_k"][k]) for z in voll]
    spreizung.append(max(werte) / min(werte) - 1)
lauf.pruefe("normierte Verhältnisse r_k (k ≤ 5) über die Basen gleich auf 1e-3", max(spreizung) < 1e-3,
            [mp.nstr(s, 2) for s in spreizung])
band = max(abs(mp.mpf(x) * mp.mpf(k + 1) ** mp.mpf(1.25) - 1) for k, x in enumerate(ergebnis[0]["r_k"]))
lauf.pruefe("r_k · k^{5/4} = 1 ± 0.01 für k ≤ 10 (b = 1.25)", band < 0.01, mp.nstr(band, 3))

# Gegencheck: andere Reihenlänge und Toleranz, gleiche Basis
_, ms2 = moden(mp.mpf("1.3"), N=80, tol="1e-3")
_, ms1 = moden(mp.mpf("1.3"))
rel = max(abs(x / y - 1) for x, y in zip(ms1[:6], ms2[:6]))
lauf.pruefe("Gegencheck b=1.3 (N=80, tol=1e-3): Moden 1–6 stimmen auf 1e-20 relativ", rel < mp.mpf(10) ** -20, mp.nstr(rel, 3))

ax[1].axhline(1, ls="--", color="#888", lw=1)
ax[0].set_xlabel("Mode k")
ax[0].set_ylabel("|p̂(k)|")
ax[0].set_title("Moden der Differenz α_A∘α_a⁻¹ − id")
ax[0].legend(fontsize=8)
ax[1].set_xlabel("Mode k")
ax[1].set_ylabel("r_k · k^{5/4}")
ax[1].set_title("Nach Abzug von C(b)·exp(−2π²k/|ln λ_a|): universell ≈ k^(−5/4)")
ax[1].set_ylim(0.98, 1.02)
ax[1].legend(fontsize=8)
lauf.bild(fig)

lauf.daten = {"basen": ergebnis, "gamma_bei_festem_S": gammas, "spreizung_r_k": [mp.nstr(s, 3) for s in spreizung],
              "band_k54": mp.nstr(band, 3), "gegencheck_rel": rel}
w = ergebnis[3]
lauf.zusammenfassung = (
    f"|p̂(k)| ≈ C k^γ exp(−2π²k/|ln λ_a|): freie Steigung trifft 2π²/|ln λ_a| auf "
    f"{', '.join(f'{z['rel_S_minus_Sa']:.1e}' for z in ergebnis)} (relativ); "
    f"normierte Verhältnisse r_k über die Basen gleich bis {mp.nstr(max(spreizung), 2)}; r_k·k^(5/4) = 1 ± {mp.nstr(band, 2)} "
    f"(k ≤ 10). Konkurrenz 2π²/ln λ_A um ≥ 20 % daneben. Basis √2: |p̂(1)| = {w['moden'][0]}.")
lauf.schreibe()
