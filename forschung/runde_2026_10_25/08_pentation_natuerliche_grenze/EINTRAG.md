# Eintrag Runde 25.10.2026 — 08 Natürliche Grenze der θ-Welle der Pentation

Dateien: `daten/aussagen/r25_08_kn.yaml`; Theorie `daten/theorie/r25_08_kn.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler.

## Angelegte ids
- kn.satz.theta_anziehend_natuerliche_grenze | Durchschiebe-Satz: an einem anziehenden Fixpunkt ist ∂Ω eine natürliche Grenze von θ; limsup|aₙ|^{1/n} = e^{2πy_c} (N1, N3, N10) — **bewiesen**
- kn.satz.theta5_natuerliche_grenze | Anwendung auf die Kneser-artige Pentation zur Basis 2 (N2, N3) — **bedingt** (K₅ nur numerisch)
- kn.beob.theta5_dach | höchster Randpunkt y_max ≈ 0.2645 bei x ≈ 0.031; Rate und Phase der aₙ treffen ihn; 96 % L̄₄-Becken unter dem Dach, kein Fraktalbeleg ≥ 10⁻³; Entkommen-Kriterium korrigiert (N4, N5, N8, N9)
- kn.beob.theta5_unter_konvergenzkreis | Konvergenzkreis keine natürliche Grenze; Lückensätze greifen nicht (N6, N7)
- kn.frage.theta5_rand | offen: Nullstellen von σ∘H₅ in Ω, ∂Ω ∩ ℍ ≠ ∅ streng, Hausdorff-Dimension, Dachspitze, e6-Hexation (Offene Fragen 1–5)

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| N1 (bewiesen) | Lemma und Satz Schritt für Schritt gelesen | richtig: E_n → 0 gleichmässig auf U′, Identitätssatz auf int U′ (zusammenhängend, U′ + n ⊂ ℍ), Becken rückwärts invariant; im Satz ist γ(t*) ∈ ∂Ω, weil {t : γ(t) ∉ Ω} abgeschlossen ist. Keine Lücke gefunden → **bewiesen** |
| N3 (bewiesen) | gelesen | richtig, sofern θ 1-periodisch ist (für die Fourier-Reihe nötig) — so in die Aussage geschrieben; {Im z > y_max} ⊂ Ω folgt aus der Definition von y_max |
| N2 (hergeleitet mit Lücke) | `s1_rand.py 100 0.01` erneut (29 s) | Ausgabe identisch mit `cache/s1_out.txt`: Funktionalgleichung max 4.08·10⁻¹⁵, Schnitt gemieden, Klassen {2288, 179, 33, 0}, Spitzen 0.26350 (x = 0.02350) und 0.26450 (x = 0.03100) → als **bedingt** eingetragen |
| N4/N5 (numerisch) | `s2_koeff.py 0.30 0.33 512 0.6 0.12` erneut | Tabelle des Berichts genau reproduziert (|a₁₀| = 6.534·10⁵, n = 18: Rate 0.2769, Phase 0.0021; n = 24: 0.2670, 0.0283). Anmerkung: Die beiden Linien weichen bei n = 32 schon um 0.18 relativ ab; verlässlich nur bis n ≈ 24 (Bericht: „bis n ≈ 30“) — so eingetragen |
| N6 (numerisch) | derselbe Lauf | Strecke [0, 200, 0, 0], Mittelwertfehler 6.81·10⁻⁸ / 1.70·10⁻⁷, |a₆₀q⁶⁰| = 1.79·10²² — wie im Bericht |
| N8 (numerisch grob) | nicht neu gerechnet | übernommen als „grob, ein Fenster“ |
| N9 | Bahn von L₃ im Lauf von s1 sichtbar | läuft nach L₄ (1.4804 + 1.6677i nach 12 Schritten) |
| N10 | gelesen | Strukturteil folgt aus N1 (in f4 des Satzes); die Behauptung über die e6-Hexation ist nicht gerechnet → nur als offene Frage (5) in kn.frage.theta5_rand |

## Änderungen am Bestand (Vorschläge)
1. **kn.verm.theta5_natuerliche_grenze**, Feld `status`: `vermutung` → `ueberholt` (ersetzt durch kn.satz.theta5_natuerliche_grenze, Kante `ersetzt` ist angelegt); `versionen` ergänzen: „25.10.2026: bewiesen unter den numerischen Annahmen über K₅ (Durchschiebe-Satz kn.satz.theta_anziehend_natuerliche_grenze); die Lücke ‚Ψ nicht ganz‘ entfällt, weil nur die lokale Umkehrung von σ nahe L₄ gebraucht wird.“ Alternativ, falls die Vermutung stehen bleiben soll: `status` → `bedingt` mit demselben `bedingungen_text` wie der neue Satz.
2. **kn.beob.theta5_wachsende_koeffizienten**, Feld `aussage`: „Die θ-Reihe konvergiert nur für Im z > y_c ≈ 0.27“ → „… für Im z > y_c = y_max ≈ 0.2645 (kn.beob.theta5_dach)“; „darunter liegen fraktal durchmischt Punkte im Becken von L̄₄ und entkommende Bahnen, höchster Punkt ausserhalb bei y ≈ 0.26“ → „darunter liegen vor allem Punkte im Becken von L̄₄ (96 % unter dem Dach) und einige entkommende Bahnen; einen fraktalen Rand zeigt die Kastenzählung auf Skalen ≥ 10⁻³ nicht; höchster Punkt ausserhalb bei y ≈ 0.2645, x ≈ 0.031“. Begründung: N4, N8 und N9 (das Kriterium |w| > 60 in s4_becken5.py zählt Punkte mit grossem |Im w|, die nach L₄ laufen, falsch als entkommend).
   Feld `pruefung.f5_neue_frage`: erste Frage beantwortet (ja, bedingt: kn.satz.theta5_natuerliche_grenze).
3. **kn.satz.kein_kneser_anziehend**, Feld `kanten`: `{typ: verwandt, ziel: kn.satz.theta_anziehend_natuerliche_grenze}` (Gegenstück: was die Welle an einem anziehenden Fixpunkt stattdessen tut).

## Vorschlag Stränge
- In den bestehenden Strang `kn.theta_raenge` („Knesers Welle auf allen Rängen“): alle fünf neuen ids. Neue Zeile `verlauf` (gemeinsam mit Bericht 09):
  „Auf geraden Rängen ist die Welle dort singulär, wo die Kneser-Lösung den Fixpunkt zwei Stufen tiefer trifft, mit Exponent ln λ_{r−4}/ln λ_{r−2} (hergeleitet bis auf die Zweigwahl, geprüft bis Rang 10); an einem anziehenden Fixpunkt (Rang 5) endet sie an einer natürlichen Grenze, dem Beckenrand.“
- Falls der Strang zu gross wird: eigener Strang `kn.theta_anziehend` „Wellen an anziehenden Fixpunkten“ mit kn.beob.theta5_wachsende_koeffizienten, kn.verm.theta5_natuerliche_grenze und den fünf neuen ids; verlauf: „An einem anziehenden Fixpunkt lebt die Welle nur auf der Becken-Komponente Ω; deren Rand ist eine natürliche Grenze, der Konvergenzkreis nicht (25.10.2026).“

## Vorschlag Ideenkarte
- Unteridee `stufen.welle` („Knesers Welle auf höheren Rängen“), `kurz` ergänzen: „… ; an einem anziehenden Fixpunkt endet die Welle an einer natürlichen Grenze, dem Rand des Beckens.“
