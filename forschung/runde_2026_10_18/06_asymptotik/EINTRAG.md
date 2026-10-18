# Eintrag Runde 18.10.2026 — 06 Asymptotik

Dateien: `daten/aussagen/r18_06_kn.yaml`, `r18_06_gr.yaml`; Theorie `daten/theorie/r18_06_{kn,gr}.yaml`.
Die Taylor-Aussagen (kn.satz.minus2_hebbar, kn.satz.taylor_darboux, kn.beob.taylor_drei_stufen) zitieren auch Bericht 05 (T1–T4), der sie unabhängig fand.

## Angelegte ids
- kn.satz.minus2_hebbar | sexp − log(z+2) holomorph auf ℂ∖(−∞,−3]; einzige Nullstelle −1 (B1; 05: T1)
- kn.kor.minus2_entwicklung | Entwicklung bei −2, Radius 1, gₙ ≈ (−1)ⁿ/(n log n) (B2, B3)
- kn.satz.taylor_darboux | cₙ = (−1)^{n+1}/(n2ⁿ) + O(ρ⁻ⁿ); c₂₂ erste negative Zahl; beantwortet gr.beob.kneser_ableitung f5 (B4, B5; 05: T1, T4)
- kn.beob.taylor_drei_stufen | Rest dₙ: Hayman-Sattel bis n = 505, dann Singularität bei −3 (B6; 05: T2, T3)
- kn.satz.wertevorrat | jeder Wert ≠ 0 unendlich oft, 0 nur bei −1 (B7)
- kn.satz.nebenblatt_zweigpunkte | unendlich viele log-Verzweigungspunkte auf Nebenblättern (B8, Skizze)
- kn.kor.doppelreihe_oben | sexp oben als Doppelreihe, reine Koeffizienten pₖκᵏ (B9)
- kn.kor.slog_bei_L | slog nahe L (B10)
- kn.kor.halbexp_band | halbe e-Funktion: Höhe ½ ± 0.0170358 im linearen log* (B12)
- gr.satz.slog_logstern_welle | slog = linearer log* + Welle in [−0.0078653, 0.0091705], kein Grenzwert (B11)

## Vorschlag Stränge
- Neuer Strang `kn.singularitaeten` „Singularitäten und Taylor-Koeffizienten von Knesers sexp“: kn.satz.singulaer_minus2, kn.satz.minus2_hebbar, kn.kor.minus2_entwicklung, kn.satz.taylor_darboux, kn.beob.taylor_drei_stufen, kn.satz.wertevorrat, kn.satz.nebenblatt_zweigpunkte, kn.satz.ableitung_nullstellenfrei, gr.beob.kneser_ableitung.
  verlauf: „Vom Logarithmus bei −2 über den hebbaren Schnitt bis −3 zu den Taylor-Vorzeichen (c₂₂ < 0) und den Nebenblättern (18.10.2026).“
- kn.kor.halbexp_band in den bestehenden Strang `kn.halbe_exp`; kn.kor.doppelreihe_oben und kn.kor.slog_bei_L zum Strang um kn.satz.binet_kappa / kn.def.theta_welle (Struktur oben); gr.satz.slog_logstern_welle zu den gr-Strängen um gr.satz.logstern.

## Vorschlag Ideenkarte
- Unteridee „Analytische Struktur von Knesers sexp“ (falls nicht vorhanden neu), kurz: „Schnitte, Nebenblätter und Taylor-Asymptotik; offen: dichte Verzweigungspunkte, zertifizierte Schranke |g| ≤ 180 auf |z| = 2.5.“

## Fehlende Quellen (nicht angelegt)
Flajolet–Odlyzko 1990 (Singularity analysis), Flajolet–Sedgewick (Analytic Combinatorics, VI.2), Hayman 1956 (admissible functions), Denjoy–Wolff.
