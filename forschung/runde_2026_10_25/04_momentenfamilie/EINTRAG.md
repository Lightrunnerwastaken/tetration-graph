# Eintrag Bericht 04 (Momentenfamilie) — 25.10.2026

Dateien: `daten/aussagen/r25_04_br.yaml`, `daten/theorie/r25_04_br.yaml`. `python baue.py --nur-pruefen`: 0 Fehler (nur erwartete Hinweise „kein Kapitel/Strang/Idee“).

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| br.satz.momente_cauchy | Momentenfamilie ⇔ Cauchy-Transformierte = 1/(ζ − c) aussen (Befund C) | bedingt |
| br.satz.sehnenkriterium | Bogen-Bild, Doppelpunkt-Kriterium, sup\|1/φ′ − 1\| < m* ⇒ doppelpunktfrei (B, D, S) | bedingt |
| br.beob.kneser_kurve_sehnen | Kneser-Kurve doppelpunktfrei, ≤ 2 Schnitte je Höhe, m* = 0.3010 (K, m*) | numerisch |
| br.satz.momente_jordan | Jordan-Fall: Momentenfamilie ⇒ Kneser; Gegenbeispiel-Bauarten scheitern (J, R) | skizze |
| br.kor.momente_umgebung | 0.7686 < φ′ < 1.4306 + Momentenfamilie ⇒ Kneser (N) | skizze |
| br.beob.gegenbeispiel_doppelpunkt | Etappe-4-Gegenbeispiele a, c haben je einen Doppelpunkt (X) | numerisch |
| br.frage.momente_doppelpunkt | Offener Fall mit Doppelpunkt (O, Offene Fragen 1–2) | offen |

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| C | Beweis gelesen (Laplace/Fubini, Umkehrung über Eindeutigkeit der Laplace-Transformation) | richtig. Die Integrierbarkeit von e^{λ\|θ_T\|} folgt aus der bewiesenen Majorante von kn.satz.theta_mittelwerte (über θ_T = θ∘φ + p, min φ′ > 0); „Re θ_T nach unten beschränkt“ braucht die Lokalform von kn.satz.theta_singularitaet, die dort an eine Gleitkomma-Zweigprüfung gebunden ist → **bedingt** statt „bewiesen“. |
| B, D, S | Beweise gelesen | richtig; das Kriterium mit m* als Infimum ist vollständig, der Zahlenwert numerisch, die Spitze qualitativ → **bedingt**. Umrechnung 1/(1 ± 0.3010) = 0.7686 / 1.4306 nachgerechnet. |
| K, m* | `s2_fein.py 16 600 kneser` (neue Abtastung, 671 Punkte) und `s4_sehnen.py fein_kneser_16_600.json 4000` bzw. `8000` | 0 Schnitte, Fläche 0.076155, Umlaufzahl +1 um c; m* = 0.3010 bei Höhe 5.45318, t = 0.2414, s = 0.74122 — wie im Bericht. |
| J | Beweis kritisch gelesen | Gerüst stimmt (F = qΦ(G)G′, F(0) = 1, F dτ = 2π dx); Lücken (a)–(d) des Berichts bestätigt; zusätzlich nicht eigens gezeigt: c ∈ D_T (C gibt nur c ∉ Ω_T). → **skizze**. |
| N | folgt aus S und J | **skizze** (erbt beide Lücken). |
| X | `s2_fein.py 16 600 a` und `… c`; max p′ aus den Fourier-Daten (`gegenbeispiel_*_40.json`) | je genau ein Schnitt: a bei x ∈ [0.3717, 0.3733] / [0.6133, 0.6150], c bei [0.3400, 0.3417] / [0.6367, 0.6383]; max p′ = 0.5985 (a), 0.6662 (c). Kleine Abweichung: max\|p\| ist 0.058 (a) bzw. 0.069 (c), Bericht schreibt „\|p\| ≤ 0.065“ für a — eingetragen als \|p\| < 0.07. |
| O, R | gelesen | in br.frage.momente_doppelpunkt bzw. br.satz.momente_jordan eingearbeitet, nicht einzeln eingetragen. |

Neue Dateien der Gegenprobe (nur im Berichtsordner): `daten/fein_kneser_16_600.json`, `daten/fein_a_16_600.json`, `daten/fein_c_16_600.json`.

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **br.frage.momentenfamilie**, Feld `aussage`, anhängen: „**Teilantwort (Runde 25.10.2026):** Ist die Wellenkurve θ_T doppelpunktfrei, so erzwingt die Momentenfamilie Kneser (br.satz.momente_jordan, Skizze), insbesondere für alle T mit 0.7686 < φ′ < 1.4306 (br.kor.momente_umgebung) — ohne Streifenbedingung. Offen bleibt nur der Fall mit Doppelpunkt (br.frage.momente_doppelpunkt).“ Status bleibt `offen`. Feld `kanten` ergänzen: `{typ: verwandt, ziel: br.satz.momente_jordan}`, `{typ: verwandt, ziel: br.frage.momente_doppelpunkt}`. Begründung: Bericht 04, Befunde J, N, O.
2. **br.satz.reelle_kennzeichnung**, `pruefung.f5_neue_frage`: „Legt die Momentenfamilie Kneser fest (br.frage.momentenfamilie)?“ → „… (br.frage.momentenfamilie; im Jordan-Fall ja, br.satz.momente_jordan)“.

## Stränge

- Vorschlag: neuer Strang **`br.momente`** (Gebiet br), Name „Die Momentenfamilie“, Mitglieder: br.frage.momentenfamilie (aus br.reell herübernehmen), br.satz.momente_cauchy, br.satz.sehnenkriterium, br.beob.kneser_kurve_sehnen, br.satz.momente_jordan, br.kor.momente_umgebung, br.beob.gegenbeispiel_doppelpunkt, br.frage.momente_doppelpunkt.
  `verlauf: Die Momente sind eine Cauchy-Transformierte; ist die Wellenkurve doppelpunktfrei, erzwingen sie Kneser (Skizze), ausdrücklich für 0.77 < φ′ < 1.43; die alten Doppelgänger haben einen Doppelpunkt, und dieser Fall ist offen.`
- Alternative: alles in **`br.reell`** aufnehmen und dort `verlauf` um „die Momentenfamilie erzwingt Kneser im Jordan-Fall (Skizze)“ ergänzen.

## Ideenkarte

Unteridee **`fortsetzung.reell`** („Reelle Kennzeichnung?“) — `straenge` um `br.momente` ergänzen (bzw. nichts zu tun, falls alles in `br.reell` geht).
