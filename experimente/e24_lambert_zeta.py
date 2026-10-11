"""exp.lambert_zeta — Übergang Lambert → Digamma für b → η, und Stetigkeit der Tetration bei η.

 Möbius-Modell M mit Fixpunkten a < A und Multiplikator λ bei a: (M(y) − a)/(M(y) − A) = λ (y − a)/(y − A).
 Bruchsumme exakt als q-Digamma (q = λ):  S^M(x) = a x + (A − a)[ψ_λ(x+1+s) − ψ_λ(1+s)]/ln λ,  λ^s = (a − 1)/(A − 1),
 ψ_q(z) = −ln(1 − q) + ln q Σ_{n≥0} q^{n+z}/(1 − q^{n+z}) → ψ(z) für q → 1.
 Grenzwert λ → 1: (A − a)/ln λ → −2e, s → Z(0) = 2e/(e − 1):  S^M → e x − 2e[ψ(x+1+Z0) − ψ(1+Z0)].
 Echte Tetration: Basen b = exp(λ e^{−λ}) (Rezeptur, a = e^λ), δ = 1 − λ → 0; T_b(1/2) → T_η(1/2)? S_b(1/2) → S_η(1/2)?
 Von oben: Kneser T_b(1/2) für b = η + 10^{−k}.
 Eingabe: T_η(1/2), S_η(1/2) aus experimente/ergebnisse/exp.eta_digamma.json.
"""
import json
from pathlib import Path

import mpmath as mp

from gemeinsam import Lauf, kneser, plt
from regulaer import Regulaer

mp.mp.dps = 30
lauf = Lauf("exp.lambert_zeta", __file__)
E = mp.e
eta_daten = json.loads((Path(__file__).parent / "ergebnisse" / "exp.eta_digamma.json").read_text(encoding="utf-8"))["daten"]
T_eta_h = mp.mpf(eta_daten["T_eta_halb"])
S_eta_h = mp.mpf(eta_daten["S_eta_halb"])
Z0 = 2 * E / (E - 1)


def psi_q(q, z, tol=mp.mpf(10) ** -28):
    s, n = mp.mpf(0), 0
    while True:
        t = q ** (n + z) / (1 - q ** (n + z))
        s += t
        n += 1
        if abs(t) < tol * (1 - q) and n > 10:
            break
    return -mp.log(1 - q) + mp.log(q) * s


def basis(lam):
    lam = mp.mpf(lam)
    ub = lam * mp.exp(-lam)
    a = mp.exp(lam)
    A = -mp.lambertw(-ub, -1).real / ub
    return mp.exp(ub), a, A


def moebius_summe(lam, x, K=200000):
    b, a, A = basis(lam)
    r0 = (a - 1) / (A - 1)
    y = lambda t: (a - A * r0 * lam ** t) / (1 - r0 * lam ** t)
    direkt = a * x + mp.fsum(y(k) - y(k + x) for k in range(1, K))
    s = mp.log(r0) / mp.log(lam)
    formel = a * x + (A - a) * (psi_q(lam, x + 1 + s) - psi_q(lam, 1 + s)) / mp.log(lam)
    return direkt, formel, s, (A - a) / mp.log(lam)


xh = mp.mpf("0.5")
d, fo, s, quot = moebius_summe(mp.mpf("0.9"), xh, K=3000)
lauf.pruefe("Möbius-Modell: Bruchsumme = q-Digamma-Formel mit q = λ (λ = 0.9)", abs(d - fo) < mp.mpf(10) ** -25, mp.nstr(abs(d - fo), 3))
abw_psi = [abs(psi_q(mp.mpf(1) - dd, mp.mpf("2.5")) - mp.digamma(mp.mpf("2.5"))) for dd in (mp.mpf("0.01"), mp.mpf("0.001"))]
lauf.pruefe("ψ_q(z) → ψ(z) für q → 1 (linear in 1 − q)", abw_psi[1] < abw_psi[0] / 8, [mp.nstr(v, 3) for v in abw_psi])

grenz_moeb = E * xh - 2 * E * (mp.digamma(xh + 1 + Z0) - mp.digamma(1 + Z0))
deltas = [mp.mpf(v) for v in ("0.2", "0.1", "0.05", "0.02", "0.01", "0.005")]
moeb, s_werte, q_werte = {}, {}, {}
for dd in deltas:
    lam = 1 - dd
    _, fo, s, quot = moebius_summe(lam, xh, K=2)
    moeb[dd], s_werte[dd], q_werte[dd] = fo, s, quot
lauf.pruefe("Grenzübergang im Modell: (A − a)/ln λ → −2e, s → 2e/(e − 1), S^M(1/2) → e/2 − 2e[ψ(3/2+Z0) − ψ(1+Z0)]",
            abs(q_werte[deltas[-1]] + 2 * E) < 0.02 and abs(s_werte[deltas[-1]] - Z0) < 0.02 and abs(moeb[deltas[-1]] - grenz_moeb) < 0.01,
            f"(A−a)/lnλ = {mp.nstr(q_werte[deltas[-1]], 6)}, s = {mp.nstr(s_werte[deltas[-1]], 6)}, S^M = {mp.nstr(moeb[deltas[-1]], 8)} vs {mp.nstr(grenz_moeb, 8)}")

# ------------------------------------------------------------------ echte reguläre Tetration für b ↑ η
echt = {}
for dd in deltas:
    lam = 1 - dd
    b, a, A = basis(lam)
    R = Regulaer(b, N=90, tol=(A - a) / 12)
    Th = R.tet(xh)
    y0, yh, summe = mp.mpf(1), Th, a * xh
    for k in range(1, 400000):
        y0, yh = mp.power(b, y0), mp.power(b, yh)
        t = y0 - yh
        summe += t
        if abs(t) < mp.mpf(10) ** -24 and k > 10 / dd:
            break
    echt[dd] = (Th, summe)
    print("  δ =", dd, "T_b(1/2) =", mp.nstr(Th, 15), "S_b(1/2) =", mp.nstr(summe, 15))
dT = {dd: abs(echt[dd][0] - T_eta_h) for dd in deltas}
dS = {dd: abs(echt[dd][1] - S_eta_h) for dd in deltas}
ordnung_T = mp.log(dT[deltas[-2]] / dT[deltas[-1]]) / mp.log(deltas[-2] / deltas[-1])
# Fassung 1 erwartete Ordnung 1 in δ; gemessen 2.0 — richtig, denn η − b ≈ η δ²/(2e): linear in η − b.
lauf.pruefe("von unten: reguläres T_b(1/2) → T_η(1/2) für b ↑ η, linear in η − b (Ordnung 2 in δ = 1 − λ)",
            all(dT[deltas[i + 1]] < dT[deltas[i]] for i in range(len(deltas) - 1)) and abs(ordnung_T - 2) < 0.1,
            f"Abstand {mp.nstr(dT[deltas[-1]], 3)} bei δ = 0.005, Ordnung {mp.nstr(ordnung_T, 4)}")
ETA = E ** (1 / E)
steig_u = {dd: (T_eta_h - echt[dd][0]) / (ETA - basis(1 - dd)[0]) for dd in deltas}
s_unten = 2 * steig_u[deltas[-1]] - steig_u[deltas[-2]]          # linear in δ extrapoliert
ordnung_S = mp.log(dS[deltas[-2]] / dS[deltas[-1]]) / mp.log(deltas[-2] / deltas[-1])
lauf.pruefe("Lambert-Bruchsumme → Digamma-Bruchsumme: S_b(1/2) → S_η(1/2)", all(dS[deltas[i + 1]] < dS[deltas[i]] for i in range(len(deltas) - 1)),
            f"Abstand {mp.nstr(dS[deltas[-1]], 3)} bei δ = 0.005, Ordnung {mp.nstr(ordnung_S, 4)}")
rest_echt = {dd: echt[dd][1] - moeb[dd] for dd in deltas}
ziel = S_eta_h - grenz_moeb
lauf.pruefe("Rest echt − Möbius konvergiert gegen den η-Rest (der Teil, den das Modell nicht trifft: (1/3) log k und C_T)",
            abs(rest_echt[deltas[-1]] - ziel) < abs(rest_echt[deltas[0]] - ziel) / 10, f"{mp.nstr(rest_echt[deltas[-1]], 8)} → {mp.nstr(ziel, 8)}")

# ------------------------------------------------------------------ von oben: Kneser für b ↓ η
oben = {}
for k in (1, 2, 3):
    try:
        w = kneser("sexp", f"exp(exp(-1))+10^-{k}", [xh], dps=30)[0]
        oben[k] = mp.re(w)
    except Exception as fehler:                                  # Motor kann nahe η versagen
        print("  Kneser b = η + 1e-%d:" % k, type(fehler).__name__, str(fehler)[:120])
if len(oben) >= 2:
    ks = sorted(oben)
    abst = [abs(oben[k] - T_eta_h) for k in ks]
    lauf.pruefe("von oben: Kneser-T_b(1/2) → T_η(1/2) für b ↓ η", all(abst[i + 1] < abst[i] for i in range(len(abst) - 1)),
                {f"η+1e-{k}": mp.nstr(oben[k] - T_eta_h, 4) for k in ks})
    steig_o = {k: (oben[k] - T_eta_h) / mp.mpf(10) ** -k for k in ks}
    s_oben = (10 * steig_o[ks[-1]] - steig_o[ks[-2]]) / 9          # linear in b − η extrapoliert
    lauf.pruefe("Steigung dT_b(1/2)/db bei η von unten und von oben gleich (numerisch C¹ in der Basis)", abs(s_unten - s_oben) < 0.01,
                f"unten {mp.nstr(s_unten, 6)}, oben {mp.nstr(s_oben, 6)}")

lauf.daten = {"T_b_halb": {mp.nstr(dd, 3): mp.nstr(echt[dd][0], 20) for dd in deltas},
              "S_b_halb": {mp.nstr(dd, 3): mp.nstr(echt[dd][1], 20) for dd in deltas},
              "moebius": {mp.nstr(dd, 3): mp.nstr(moeb[dd], 15) for dd in deltas}, "moebius_grenze": mp.nstr(grenz_moeb, 20),
              "ordnung_T": mp.nstr(ordnung_T, 6), "ordnung_S": mp.nstr(ordnung_S, 6),
              "kneser_oben": {f"1e-{k}": mp.nstr(v, 20) for k, v in oben.items()},
              "steigung_unten": mp.nstr(s_unten, 8), "steigung_oben": mp.nstr(s_oben, 8) if len(oben) >= 2 else None}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
ax[0].loglog([float(dd) for dd in deltas], [float(dT[dd]) for dd in deltas], "o-", label="|T_b(1/2) − T_η(1/2)| (regulär, von unten)")
ax[0].loglog([float(dd) for dd in deltas], [float(dS[dd]) for dd in deltas], "s-", label="|S_b(1/2) − S_η(1/2)|")
if oben:
    ax[0].loglog([float(E ** (-1 / E) * 0 + mp.mpf(10) ** -k) for k in sorted(oben)], [float(abs(oben[k] - T_eta_h)) for k in sorted(oben)], "^--",
                 label="|T_b(1/2) − T_η(1/2)| (Kneser, b − η)")
ax[0].set_xlabel("δ = 1 − λ  bzw.  b − η"); ax[0].set_title("Stetigkeit bei η"); ax[0].legend(fontsize=7)
qs = [1 - mp.mpf(10) ** (-mp.mpf(k) / 4) for k in range(4, 17)]
ax[1].semilogx([float(1 - q) for q in qs], [float(psi_q(q, mp.mpf("2.5"))) for q in qs], "o-", ms=3, label="ψ_q(2.5)")
ax[1].axhline(float(mp.digamma(mp.mpf("2.5"))), ls="--", color="#b5523b", label="ψ(2.5)")
ax[1].set_xlabel("1 − q"); ax[1].set_title("q-Digamma → Digamma"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"Möbius-Modell: Lambert-Bruchsumme = q-Digamma (q = λ), Grenzwert ψ. Echte reguläre Tetration: T_b(1/2) → T_η(1/2) "
                        f"(Ordnung {mp.nstr(ordnung_T, 3)} in 1 − λ), S_b(1/2) → S_η(1/2). Kneser von oben: "
                        + (", ".join(f"η+1e-{k}: {mp.nstr(oben[k] - T_eta_h, 3)}" for k in sorted(oben)) if oben else "nicht berechenbar") + ".")
lauf.schreibe()
