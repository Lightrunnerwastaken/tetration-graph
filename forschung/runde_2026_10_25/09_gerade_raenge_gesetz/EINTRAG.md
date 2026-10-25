# Eintrag Runde 25.10.2026 — 09 Gesetz für die θ-Singularität bei geraden Rängen

Dateien: `daten/aussagen/r25_09_kn.yaml`; Theorie `daten/theorie/r25_09_kn.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler.

## Angelegte ids
- kn.lem.umkehr_minus_unendlich | Kern-Lemma: H_s⁻¹(V) = z_s + Φ_s(λ_{s−2}^V), Φ_s′(0) = c_{s−1}/A_s (G2, G3) — **bewiesen** (unter den genannten Voraussetzungen)
- kn.satz.theta_gerade_raenge_struktur | θ_r = −z + Θ_r((z − x_r*)^{α_r} g_r), α_r = ln λ_{r−4}/ln λ_{r−2}, Exponenten kα_r + m, Koeffizientenreihe (G4, G5) — **bedingt** (Zweigwahl)
- kn.beob.theta10_singularitaet | Rang 10: Ort x₁₀* = −6.6936050663, kp₁₀, μ, λ₈, A₁₀ (G1, G7)
- kn.beob.theta_gerade_raenge_exponenten | Fit α₈ = 0.7302, α₁₀ = 0.8360; Folge α_r; roher Exponent taugt nicht (G6, G7, G9)
- kn.verm.alpha_r_gegen_1 | Vermutung α_r → 1, z_{r−2} − p_{r−2} → 0 (G8)

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| G3 (bewiesen) | Beweis gelesen | richtig als lokale Folgerung: λ^{H_s} = σ(H_s(z+1) − p)/c ist bei z_s analytisch mit einfacher Nullstelle; V nur modulo 2πi/ln λ bestimmt, das legt den Zweig fest. Voraussetzung „Funktionalgleichung auf dem Zweig“ in die Aussage geschrieben → **bewiesen** als Lemma mit expliziten Voraussetzungen |
| G2 (numerisch) | `p2_lemma.py 20 rang9` neu (Motor ohne Zwischenspeicher, 78 s) | `daten/lemma_20.json` byte-gleich mit der Fassung des Berichts; Quotienten 2.00013049 / 5.17179703 / 9.49303563, Konstanten 0.779424 / 29.3046705 / 11330.7097 (Ziele 0.779373 / 29.3046709 / 11330.7170) |
| G4 (hergeleitet mit Lücke) | Herleitung gelesen | richtig bis auf die genannte Zweigwahl; zusätzlich gehört die Analytizität von h bei ψ_r in die Lücke (Bericht begründet sie mit h′ ≠ 0, das genügt nicht) → **bedingt**, Beweis als `skizze` mit drei Lücken. Für r = 6 stimmen θ* und K mit kn.satz.theta6_singularitaet überein (Formel nachgerechnet: (A₆/λ₄^{pcent})^α = A₆^α b^{−pcent}) |
| G1/G6 (numerisch, Rang 10) | **zweite Einstellung**, die im Bericht fehlte: `p4_theta_r.py 20 0.025 64 0.03 10` (Linie y = 0.025, 3 min), dann dasselbe Fit-Modell wie `p5_fit.py`; ausserdem `p5_fit.py` erneut | aₙ beider Linien gleich auf ≤ 4.8·10⁻⁶ relativ (n ≤ 31); arg(a₃₁/a₃₀) = −1.925615 wie im Bericht; α₁₀ = 0.8360 auf beiden Linien, c₁ = 0.00214 − 0.0268i, Rest mit Vorhersage-α 7.9·10⁻⁶ (Bericht 7.8·10⁻⁶); p5_fit.py: alle Zahlen gleich. kp₁₀, μ, p₈, λ₈, x₁₀* dabei neu gerechnet und gleich |
| Motor-Grenze | Versuche y = 0.05, y = 0.04 (64 Punkte) und y = 0.03 mit 128 Punkten | scheitern an einzelnen Punkten x ≈ 0.27–0.30 (Rang-10-Auswertung gibt Fehler) — als Hinweis in kn.beob.theta10_singularitaet vermerkt |
| G8 (Vermutung) | Differenzen aus den Tabellenwerten nachgerechnet | x_r* − p_{r−2} und z_{r−2} − p_{r−2} stimmen |

Zwischenspeicher der Gegenprobe nur im Berichtsordner: neue Datei `daten/theta10_20_0.025_64.json`, neue Einträge in `cache/2_20.json`; `daten/lemma_20.json` wurde gleich wieder geschrieben.

## Änderungen am Bestand (Vorschläge)
1. **kn.verm.theta_gerade_raenge**, Feld `status`: `vermutung` → `ueberholt` (ersetzt durch kn.satz.theta_gerade_raenge_struktur, Kante `ersetzt` ist angelegt). Feld `nachgeprueft` ergänzen: „25.10.2026: Rang 8 α = 0.7302 (Fit mit Korrekturgliedern), Rang 10 Ort x₁₀* = −6.6936050663 und α = 0.8360 (zwei Linien).“ Feld `gegenbeispielsuche`: „Rang 10 nicht gerechnet“ → „Rang 10 gerechnet (kn.beob.theta10_singularitaet), Ort und Exponent getroffen“. Der dort genannte lokale Exponent 1.713 (Ziel 1.730) ist durch das Korrekturglied n^{−α} erklärt (kn.beob.theta_gerade_raenge_exponenten, Kante `korrigiert`). Alternativ, falls die Vermutung stehen bleiben soll: `status` → `bedingt` mit dem `bedingungen_text` des neuen Satzes.
2. **kn.satz.theta6_singularitaet**, Feld `pruefung.f4_uebertragbar` („Auf gerade Ränge allgemein, siehe kn.verm.theta_gerade_raenge“) → „… hergeleitet in kn.satz.theta_gerade_raenge_struktur“.
3. **kn.satz.theta6_koeffizienten**, Feld `pruefung.f4_uebertragbar`: Verweis auf kn.satz.theta_gerade_raenge_struktur (ii) statt kn.verm.theta_gerade_raenge.

## Vorschlag Stränge
- In den bestehenden Strang `kn.theta_raenge` („Knesers Welle auf allen Rängen“): alle fünf neuen ids. Neue Zeile `verlauf` (gemeinsam mit Bericht 08):
  „Auf geraden Rängen ist die Welle dort singulär, wo die Kneser-Lösung den Fixpunkt zwei Stufen tiefer trifft, mit Exponent ln λ_{r−4}/ln λ_{r−2} (hergeleitet bis auf die Zweigwahl, geprüft bis Rang 10); an einem anziehenden Fixpunkt (Rang 5) endet sie an einer natürlichen Grenze, dem Beckenrand.“

## Vorschlag Ideenkarte
- Unteridee `stufen.welle` („Knesers Welle auf höheren Rängen“): `kurz` unverändert passend; eventuell ergänzen „… — ein Kern-Lemma erklärt alle geraden Ränge zugleich.“
