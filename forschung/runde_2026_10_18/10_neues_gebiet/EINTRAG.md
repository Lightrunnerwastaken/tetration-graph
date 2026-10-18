# Eintrag Bericht 10 (Erzeuger und Zufallstürme) in den Graphen

Dateien: `daten/aussagen/r18_10_it.yaml`, `r18_10_kl.yaml`, `r18_10_fp.yaml`; Theorie `daten/theorie/r18_10_{it,kl,fp}.yaml`.
Prüfung `baue.py --nur-pruefen`: 0 Fehler.

## Angelegte ids
- it.satz.kein_ganzer_erzeuger | kein ganzer/meromorpher Erzeuger von b^z (bewiesen, folklore)
- it.lem.erzeuger_starr | holomorpher Erzeuger am Fixpunkt = c·J_p (bewiesen, klassisch)
- it.satz.erzeuger_resonanz | gemeinsamer Erzeuger an zwei Fixpunkten nur bei λ_p = λ_q^m; über η nie, Möbius m = −1 (bewiesen)
- it.satz.kneser_erzeuger | J_K = J_L(1+θ′∘slog_K), (w−L)^(1+2πi/L), Basis e: C⁵ nicht C⁶ (bedingt: a₁ ≠ 0)
- it.satz.regulaerer_erzeuger_fortsetzung | J_L fortsetzbar über ℂ∖{0,1,e,…}; 0, 1 verzweigen (bewiesen + numerisch)
- kl.satz.resonanzbasen | λ_aλ_A < 1, Resonanzbasen b₂ = 1.04480…; sonst reguläre Iterationen verschieden (bewiesen, Werte numerisch)
- kl.beob.resonanz_abweichung | an b₂, b₃ trotzdem verschieden (numerisch)
- fp.satz.zufallsturm_unter_eins | Basen in (0,1): zwei Häufungspunkte, 0-1-Gesetz, robuste Zweierschwingung (bewiesen)
- fp.satz.zufallsturm_divergenz | Basen ≥ 1 mit Masse über η: Divergenz f.s. (bewiesen)
- fp.beob.zufallsturm_verschmelzung | U(0,β)-Turm konvergiert auch unter e^(−e) (numerisch)
- fp.verm.zufallsturm_konvergenz | Konvergenz f.s. für jedes β > 0 (Vermutung)
- fp.verm.zufallsturm_slog_schwelle | slog-Schwelle κ > 1 (Vermutung, heuristisch)

Nicht eigens eingetragen: A6 (keine Resonanz zwischen Fixpunkten von exp, numerisch) — steht in f4 von it.satz.erzeuger_resonanz. Ehrlichkeitsvermerk: die Eindeutigkeit von b_k je k ist nur numerisch (im Bericht als bewiesen nur „Existenz“).

## Vorschlag Stränge
- Neuer Strang `it.erzeuger` (Gebiet it), Name „Der Erzeuger — Tetration als Strömung“, verlauf: „Die Julia-Gleichung ist linear: kein ganzer Erzeuger, Starrheit am Fixpunkt, zwei Fixpunkte nur bei Resonanz; Knesers Feld ist an L genau C⁵.“ Mitglieder: die fünf it-ids.
- kl.satz.resonanzbasen, kl.beob.resonanz_abweichung → bestehender Strang `kl.zwei_abel`.
- fp.satz.zufallsturm_unter_eins, fp.beob.zufallsturm_verschmelzung, fp.verm.zufallsturm_konvergenz → `fp.periode2`; fp.satz.zufallsturm_divergenz, fp.verm.zufallsturm_slog_schwelle → `fp.turm`. Alternativ ein neuer Strang `fp.zufall` „Türme mit zufälligen Basen“, verlauf: „Zufall zerstört Eulers Zweierschwingung; über 1 entscheidet die Häufigkeit von Basen nahe 1, vermutlich mit slog-Schwelle.“

## Vorschlag Ideenkarte
- Teil A: unter `dynamik.koenigs` oder neue Unteridee `dynamik.erzeuger` („Der Erzeuger J: wo eine gebrochene Iteration singulär sein muss“).
- Teil B: `turm.konvergenz` und `turm.periode2`; oder neue Unteridee `turm.zufall` („Unendliche Türme mit zufälligen Basen“).
