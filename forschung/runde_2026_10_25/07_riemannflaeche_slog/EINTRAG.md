# Eintrag Runde 25.10.2026 — 07 Riemannsche Fläche von Knesers slog

Dateien: `daten/aussagen/r25_07_kn.yaml`, `daten/aussagen/r25_07_it.yaml`; Theorie `daten/theorie/r25_07_kn.yaml`, `daten/theorie/r25_07_it.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler (die zwei Warnungen gehören zu fremden bw-Einträgen).

## Angelegte ids
- kn.beob.slog_taylor_radius_L | Taylorreihe von slog um 0: Leitglied −2Re(1/(nL^{n+1})), Radius |L|, Rest n^{−Re γ} (R1)
- kn.satz.slog_periodische_halbebene | slog = slog(eʷ) − 1 holomorph und 2πi-periodisch auf Re w < Re L; Log-Gitter L + 2πiℤ, L̄ + 2πiℤ (R2) — **bedingt**
- kn.beob.sexp_trifft_L | Kneser-Blatt an L + 2πik regulär, s*ₖ auf 28 Stellen; sexp(s*₁ + 1) = L (R3, R4)
- kn.kor.slog_monodromie_L | Umlauf um L: u ↦ u + γ, M_L(s) = (u+γ) + ϑ(u+γ); Naht des Motors (R5, R6, R14)
- kn.satz.sexp_asymptotischer_wert_null | 0 asymptotischer Wert von sexp; Log-Punkt der slog-Fläche über 0 (R8)
- kn.satz.slog_verzweigung_dicht | singuläre Punkte voll exp-invariant, enthalten G(L) ∪ G(L̄) ∪ G(0), dicht in ℂ (R9, R10, R11)
- kn.verm.slog_verzweigung_genau | Vermutung: genau die drei Bahnen, nur logarithmische Punkte (R12, R13)
- it.kor.poincare_exp_ableitung_produkt | Ψ′(ζ) = exp Σ(Ψ(ζ/Lᵏ) − L), Ψ ohne kritische Punkte (R7) — Neuheit `folklore`

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| R1 (numerisch) | `s1_radius.py` erneut | alle Zahlen gleich (Quotient 1.00000000002 bei n = 250, 300; Steigungen −5.4 … −4.5) |
| R2 (Bericht: bewiesen) | Beweis gelesen | richtig, aber nicht unbedingt: braucht Holomorphie von A auf \|w\| < \|L\| (nur numerisch, kn.beob.slog_scheibe) und für „genau“ die Lage der Randsingularitäten (R1, numerisch) → als **bedingt** eingetragen |
| R3/R4 (numerisch) | `s4_schnitt_und_sstern.py 35` (neue Motorwerte, dps 35) und **unabhängig ohne Motor**: Taylorreihe von sexp um 0 aus `runde_2026_10_18/06_asymptotik/daten/taylor_0_1.5_400_50.json` an s*ₖ − 2, dann drei exp-Schritte | s*₁, s*₂, s*₃, s*₋₁ gleich auf allen 28 Stellen; \|sexp(s*ₖ) − (L+2πik)\| ≤ 1.2·10⁻³⁰, \|sexp(s*ₖ+1) − L\| ≤ 1.6·10⁻³⁰ (k = 1, 2, −1). Zusatz: Dass sexp den Wert L annimmt, ist über kn.satz.wertevorrat sogar bewiesen |
| R5 (bewiesen + numerisch) | Beweis gelesen; `s5b_monodromie_genau.py 35` | Beweis richtig (Im γ > 0, also bleiben Umläufe gegen den Uhrzeiger im Holomorphiegebiet von ϑ; für k < 0 nahe genug bei L). Als **Korollar** von kn.kor.slog_bei_L eingetragen. Abweichung 2.8·10⁻¹⁵ / 1.9·10⁻¹⁵ / 2.7·10⁻¹¹ / 6.2·10⁻⁸ (y = 1.8 / 2 / 2.5 / 3), wie im Bericht |
| R7 (bewiesen) | Rechnung nachvollzogen | richtig; die Nullstellenfreiheit stand schon in kn.satz.theta_singularitaet (i), neu ist nur die Produktformel (allgemein bekannt für Poincaré-Funktionen) → `folklore` |
| R8 (bewiesen, „Lücke nur …“) | Beweis gelesen | richtig; die genannte Lücke ist geschlossen: der Hauptzweig von slog setzt sich über slog(v) = slog(log v) + 1 auf Scheiben mit nach unten beschränktem Radius um alle v ≥ 0 fort, und log ζ liegt für Re ζ gross, \|Im ζ\| < 3π/2 nahe der Achse → **bewiesen** |
| R9–R11 (bewiesen) | Beweis gelesen | richtig. Eine Stelle verbessert: „alle Urbilder“ führt der Bericht über den Weg 0 → 2πik auf der imaginären Achse, das braucht A holomorph auf dem Einheitskreis (nur numerisch). Ersetzt durch eine k-fache kleine Schleife um 0, wo der Hauptzweig sicher holomorph ist (sexp biholomorph nahe [−1, 0]). L ∈ B über den Weg sexp(it), t → ∞ (nur sexp′ ≠ 0 und sexp(it) → L nötig) → **bewiesen** |
| R10 | Literatur | J(exp) = ℂ (misiurewicz1981) und Dichte der Rückwärtsbahnen in J (Standardsatz, mit fatou1926 zitiert, Stelle nicht nachgeschlagen) |
| R12, R13 | — | Vermutung; keine Gegenbeispielsuche im Bericht, so vermerkt |

Nicht eingetragen: R6 und R14 als eigene Aussagen (stehen als Nachprüfung bzw. Grenze der ϑ-Reihe in kn.kor.slog_monodromie_L).

## Änderungen am Bestand (Vorschläge)
1. **kn.satz.wertevorrat**, Feld `beweis.text`: anhängen „Für c = L, L̄ greift log^[k] c → L nicht (σ(L) = 0); stattdessen c′ = L + 2πi: sexp(w) = L + 2πi ⇒ sexp(w+1) = L (Beispiel kn.beob.sexp_trifft_L).“ Begründung: Der Beweis setzt σ(log^[k] c) ≠ 0 voraus; für c = L ist das falsch, die Aussage bleibt aber richtig.
2. **kn.satz.nebenblatt_zweigpunkte**, Feld `pruefung.f5_neue_frage`: „Liegen die Verzweigungspunkte dicht in ℂ?“ → ergänzen „Für slog ja (kn.satz.slog_verzweigung_dicht); für sexp offen, folgt nicht direkt.“
3. **kn.beob.slog_scheibe**, Feld `kanten`: Kante `{typ: verwandt, ziel: kn.satz.slog_periodische_halbebene}` (die Holomorphie auf D ist dort die Bedingung).
4. **kn.kor.slog_bei_L**, Feld `pruefung.f5_neue_frage` („Wie sieht slog auf den Nebenblättern nahe L aus?“): beantwortet durch kn.kor.slog_monodromie_L (Zweige sₖ) und kn.beob.sexp_trifft_L (reguläre Blätter über L).

## Vorschlag Stränge
- Neuer Strang `kn.slog_flaeche` „Die Riemannsche Fläche von Knesers slog“ (Gebiet kn), Mitglieder: kn.beob.slog_taylor_radius_L, kn.satz.slog_periodische_halbebene, kn.beob.sexp_trifft_L, kn.kor.slog_monodromie_L, kn.satz.sexp_asymptotischer_wert_null, kn.satz.slog_verzweigung_dicht, kn.verm.slog_verzweigung_genau.
  verlauf: „Von der Scheibe |w| < |L| über das 2πi-periodische Blatt mit seinem Log-Gitter und die exakte Monodromie um L zu den drei Grundpunkten L, L̄, 0 und einer dichten, exp-invarianten Verzweigungsmenge (25.10.2026).“
  (Falls Bericht 06, Riemannsche Fläche von sexp, einen eigenen Strang vorschlägt: beide zu einem Strang „Riemannsche Flächen von sexp und slog“ zusammenlegen.)
- it.kor.poincare_exp_ableitung_produkt in den bestehenden Strang `it.hyperbolisch` (oder `it.erzeuger`, neben it.satz.regulaerer_erzeuger_fortsetzung).

## Vorschlag Ideenkarte
- Den neuen Strang unter der bestehenden Unteridee `fortsetzung.kneser` („Kneser über η“) anhängen; `kurz` ergänzen um „… und die Riemannsche Fläche von slog: periodisches Blatt, Monodromie um L, dichte Verzweigung.“

## Fehlende Quellen (nicht angelegt)
Bergweiler, „Iteration of meromorphic functions“, Bull. AMS 29 (1993) — für die Dichte der Rückwärtsbahnen in J. (Für J(b^z) = ℂ bei b > η könnte das vorhandene `devaney1984` zitiert werden; nicht geprüft.)
