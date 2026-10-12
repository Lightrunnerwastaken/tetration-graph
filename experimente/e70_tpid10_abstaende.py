"""exp.tpid10_abstaende — TPID 10 und 11 an den bekannten Gliedern prüfen; Abgleich mit Levensteins Schranken (2017).

 Glieder (n = Stufe, η₃ = η, e₃ = e): η₄, e₄ nach Levenstein 2010 (Forum tid 262, 32 Stellen, = kn.beob.tetra_euler),
 η₅, e₅ aus exp.stufe5, η₆, e₆ aus exp.hexation_kneser (Werte aus kn.beob.stufe5, kn.beob.stufe6).
 Levenstein 2017 (Forum tid 1177, Beitrag 5): beschränkt bis ≈ 1.738 (Hexation, also η₅), 1.799 (Heptation, η₆), 1.839 (Oktation, η₇).
 TPID 10 (dantheman163, 2010): eₙ → 4, ηₙ → 2.   TPID 11 (tommy1729, 2010): ηₙ² < e_{n−1} und
 (ηₙ² − e_{n−1}) / (η_{n−1}² − e_{n−2}) → 1.
 Reine Arithmetik mit den gespeicherten Konstanten, kein Kneser-Motor.
"""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 40
lauf = Lauf("exp.tpid10_abstaende", __file__)
ETA = {3: mp.e ** (1 / mp.e), 4: mp.mpf("1.6353244967152763993453446183062"),
       5: mp.mpf("1.737355925151687573101390"), 6: mp.mpf("1.7991212361288126867235670070")}
EUL = {3: +mp.e, 4: mp.mpf("3.0885322718067176544821807826411"),
       5: mp.mpf("3.301407278523367493377305"), 6: mp.mpf("3.44223464771689323")}
LEVENSTEIN = {5: mp.mpf("1.738"), 6: mp.mpf("1.799"), 7: mp.mpf("1.839")}   # „~1.738“, „1.799“, „1.839“
N = sorted(ETA)

# ------------------------------------------------------------------ Monotonie und Abstände
lauf.pruefe("ηₙ und eₙ wachsen (n = 3 … 6)", all(ETA[n] < ETA[n + 1] and EUL[n] < EUL[n + 1] for n in N[:-1]))
d = {n: 2 - ETA[n] for n in N}
D = {n: 4 - EUL[n] for n in N}
lok_eta = {n: mp.log(d[n] / d[n + 1]) / mp.log(mp.mpf(n + 1) / n) for n in N[:-1]}
lok_eul = {n: mp.log(D[n] / D[n + 1]) / mp.log(mp.mpf(n + 1) / n) for n in N[:-1]}
quot_eta = {n: d[n + 1] / d[n] for n in N[:-1]}
lauf.pruefe("2 − ηₙ: lokale Potenz-Exponenten fast gleich (Spanne < 0.02), geometrische Quotienten steigen",
            max(lok_eta.values()) - min(lok_eta.values()) < mp.mpf("0.02") and quot_eta[3] < quot_eta[4] < quot_eta[5],
            {f"{n}→{n+1}": mp.nstr(lok_eta[n], 6) for n in lok_eta})
print("  4 − eₙ, lokale Exponenten:", {f"{n}→{n+1}": mp.nstr(lok_eul[n], 6) for n in lok_eul})

# Potenzgesetz 2 − ηₙ = A n^(−p) durch die letzten zwei Glieder; Vorhersage η₇
p = lok_eta[5]
A = d[6] * mp.mpf(6) ** p
eta7 = 2 - A * mp.mpf(7) ** (-p)
p_alt = lok_eta[4]
eta7_alt = 2 - d[6] * (mp.mpf(6) / 7) ** p_alt
lauf.pruefe("Vorhersage η₇ aus dem Potenzgesetz trifft Levensteins ≈ 1.839 (|Δ| < 0.002, beide Exponenten)",
            abs(eta7 - LEVENSTEIN[7]) < mp.mpf("0.002") and abs(eta7_alt - LEVENSTEIN[7]) < mp.mpf("0.002"),
            f"η₇ ≈ {mp.nstr(eta7, 6)} (p = {mp.nstr(p, 5)}), {mp.nstr(eta7_alt, 6)} (p = {mp.nstr(p_alt, 5)})")
lauf.pruefe("Levensteins η₅ ≈ 1.738 und η₆ ≈ 1.799 stimmen mit den 30-Stellen-Werten auf die angegebenen Stellen",
            abs(ETA[5] - LEVENSTEIN[5]) < mp.mpf("0.0007") and abs(ETA[6] - LEVENSTEIN[6]) < mp.mpf("0.0002"),
            f"Δ₅ = {mp.nstr(ETA[5] - LEVENSTEIN[5], 3)}, Δ₆ = {mp.nstr(ETA[6] - LEVENSTEIN[6], 3)}")
# e₇: Potenzgesetz für 4 − eₙ (Exponent der letzten Stufe) — nur als grobe Hochrechnung
e7 = 4 - D[6] * (mp.mpf(6) / 7) ** lok_eul[5]
print(f"  Hochrechnung e₇ ≈ {mp.nstr(e7, 5)} (Levenstein 2017 sah bei Basis 1.84 eine Fast-Berührung „nahe 3.7“)")

# ------------------------------------------------------------------ TPID 11
t11 = {n: ETA[n] ** 2 - EUL[n - 1] for n in N[1:]}
lauf.pruefe("TPID 11 (1): ηₙ² < e_{n−1} für n = 4, 5, 6", all(v < 0 for v in t11.values()),
            {n: mp.nstr(v, 6) for n, v in t11.items()})
q11 = {n: t11[n] / t11[n - 1] for n in N[2:]}
lauf.pruefe("TPID 11 (2): Quotienten (ηₙ² − e_{n−1})/(η_{n−1}² − e_{n−2}) — zwei Werte, keine Entscheidung möglich",
            True, {n: mp.nstr(v, 6) for n, v in q11.items()})
rel = {n: t11[n] / EUL[n - 1] for n in N[1:]}
print("  relativer Abstand (ηₙ² − e_{n−1})/e_{n−1}:", {n: mp.nstr(v, 5) for n, v in rel.items()})
# ηₙ² − e_{n−1} = D_{n−1} − (4dₙ − dₙ²): (1) gilt genau dann, wenn r_n = D_{n−1}/(4dₙ − dₙ²) < 1
r = {n: D[n - 1] / (4 * d[n] - d[n] ** 2) for n in N[1:]}
lauf.pruefe("TPID 11 umgeschrieben: r_n = (4 − e_{n−1})/(4(2−ηₙ) − (2−ηₙ)²) < 1 und fallend (n = 4, 5, 6)",
            all(v < 1 for v in r.values()) and r[4] > r[5] > r[6], {n: mp.nstr(v, 6) for n, v in r.items()})

# ------------------------------------------------------------------ Robbins' exakte Punkte der Selbst-Superwurzel (tid 269)
s2 = mp.log(2) / mp.lambertw(mp.log(2)).real
s3 = mp.findroot(lambda a: a ** (a ** a) - 3, mp.mpf("1.6"))
lauf.pruefe("Selbst-Superwurzel: s(2) = ln2/W(ln2) = 1.559610469462…, s(3) = 1.635078474636… (Robbins 2009); s(3) < η₄",
            abs(s2 - mp.mpf("1.559610469462")) < mp.mpf(10) ** -12 and abs(s3 - mp.mpf("1.635078474636")) < mp.mpf(10) ** -12 and s3 < ETA[4],
            f"s(2) = {mp.nstr(s2, 15)}, s(3) = {mp.nstr(s3, 15)}, η₄ − s(3) = {mp.nstr(ETA[4] - s3, 4)}")

lauf.daten = {"eta": {n: ETA[n] for n in N}, "eul": {n: EUL[n] for n in N}, "lok_eta": lok_eta, "lok_eul": lok_eul,
              "eta7_vorhersage": [eta7, eta7_alt], "e7_hochrechnung": e7, "tpid11": t11, "tpid11_quotient": q11,
              "levenstein2017": LEVENSTEIN, "r_tpid11": r, "superwurzel": {"s2": s2, "s3": s3}}
lauf.zusammenfassung = (f"2 − ηₙ fällt wie n^(−{mp.nstr(p, 4)}); Vorhersage η₇ ≈ {mp.nstr(eta7, 5)} gegen Levenstein 2017 ≈ 1.839. "
                        f"TPID 11 (1) gilt für n = 4…6; die Quotienten {', '.join(mp.nstr(v, 4) for v in q11.values())} entscheiden (2) nicht.")

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ns = [float(n) for n in N]
ax[0].loglog(ns, [float(d[n]) for n in N], "o-", label="2 − ηₙ")
ax[0].loglog(ns, [float(D[n]) for n in N], "s-", label="4 − eₙ")
ax[0].loglog([7.0], [float(2 - eta7)], "x", color="C0", label="Hochrechnung η₇")
ax[0].loglog([7.0], [float(2 - LEVENSTEIN[7])], "+", color="k", ms=10, label="Levenstein 2017 (η₇)")
ax[0].set_xlabel("Stufe n"); ax[0].legend(fontsize=8); ax[0].set_title("Abstände zu 2 und 4")
ax[1].plot([4, 5, 6], [float(t11[n]) for n in N[1:]], "o-")
ax[1].axhline(0, color="k", lw=0.6); ax[1].set_xlabel("n"); ax[1].set_title("TPID 11: ηₙ² − e_{n−1}")
lauf.bild(fig)
lauf.schreibe()
