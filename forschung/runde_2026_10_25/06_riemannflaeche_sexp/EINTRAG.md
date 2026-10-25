# Eintrag Bericht 06 (Riemannsche Fläche von sexp) — 25.10.2026

Dateien: `daten/aussagen/r25_06_kn.yaml`, `daten/theorie/r25_06_kn.yaml`. `python baue.py --nur-pruefen`: 0 Fehler (nur erwartete Hinweise).

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| kn.satz.riemannflaeche_logueberlagerung | Fläche = volle Log-Überlagerung der eigenen Verschiebung; Zweigwerte = alle Log-Ketten von sexp(z+K) (R1) | bewiesen |
| kn.satz.verzweigungsorte_charakterisierung | z verzweigt ⇔ sexp(z+M) = sexp(M−n), n ≥ 2; Typ log^[n−1]; Nullstellen einfach (R2) | bewiesen |
| kn.satz.nebenblatt_nullstellen_rechts | Re sexp(ζ+N+k) ≥ sexp(N+k−1); feste Tiefe nirgends dicht (R3, R4) | bewiesen |
| kn.satz.verzweigung_dicht | Verzweigungsorte dicht in ℂ, abzählbar; je Zweig lokal endlich (R5, R6) | bewiesen |
| kn.beob.verzweigung_gitter | Reihen mit Abstand γ = 2πi/L, Fehlerverhältnis e^{2πiγ}; keine echte Gitterperiode (R7, R10) | numerisch |
| kn.beob.verzweigung_nah | Beispiele: Tiefe 1 bei 0.300141 − 0.006244i usw. (R8, R9) | numerisch |
| kn.frage.verzweigung_abstand | Abstand fester Tiefe, Dimension, Streifen −2 < Re ≤ −1, andere Basen (Offene Fragen 1–4) | offen |

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| R1 | Beweis gelesen | richtig: Fortsetzbarkeit von sexp ⇔ Fortsetzbarkeit von S = sexp(·+1) ohne Nullstelle; Überlagerung zusammenhängend wegen der einfachen Nullstelle −2 von S; Induktion gibt alle Ketten. **bewiesen**. |
| R2 | Beweis gelesen | richtig (Eindeutigkeit von exp zwingt die Kette nach einer 0 auf 1, e, e^e, …; einfache Nullstellen über sexp′ ≠ 0; Singularitäten nur, wo ein Zwischenwert 0 wird). **bewiesen**. |
| R3, R4 | Rückwärtsinduktion nachgerechnet; Normalitätsargument gelesen | richtig; R4 stützt sich auf „Fatou-Menge von exp leer“ (Misiurewicz 1981, Quelle `misiurewicz1981`, ungeprüft). **bewiesen**. |
| R5, R6 | gelesen; bei R6 geprüft, dass auch Fortsetzungen innerhalb der Scheibe nur in F_{N′}⁻¹({E_1,…,E_{N′}}) singulär werden | richtig; R5 braucht den Satz über dichte Rückwärtsbahnen (Bergweiler 1993, nicht in der Quellenliste). **bewiesen**. |
| R7 | `s1_familie.py 45` (neue Motorwerte bei dps 45) | alle w_j (k = 1, 2, −1) auf den 22 gedruckten Stellen wie bei dps 40 und 50; δ-Werte gleich; δ₁/δ₀ von Hand nachgerechnet: −0.0012284 + 0.0004250i (k = 1), −0.0012326 + 0.0004181i (k = 2); e^{2πiγ} = −0.0012262 + 0.0004246i — bestätigt. Status bleibt **numerisch** (Restglied nicht abgeschätzt). |
| R8, R9 | `s2_dicht.py 45` | alle drei Punkte auf ≥ 15 Stellen wie bei dps 40 (z. B. 0.300141218682681473 − 0.00624399328462789499i), Kettenprobe ≤ 4·10⁻¹³⁴; Bedingungen nachvollzogen: 0.3 → sexp(z+4) = sexp(2), 1.3+0.1i → sexp(z+3) = sexp(0), −0.5+0.25i → sexp(z+4) = sexp(0). |
| R10 | gelesen | trivial richtig; in kn.beob.verzweigung_gitter eingearbeitet. |

Neue Dateien der Gegenprobe (nur im Berichtsordner): `daten/familie_45.json`, `daten/dicht_45.json`, zusätzliche Einträge in `daten/motor_cache.json` (diesmal ohne Dateisperre).

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **kn.satz.nebenblatt_zweigpunkte** (r18_06_kn.yaml):
   - `status`: skizze → **bewiesen**; `beweis.art`: skizze → vollstaendig; `beweis.luecken` streichen („Dass jeder Turm … nur für n ≤ 2 ausgeführt“). Begründung: Die Erreichbarkeit aller Türme folgt aus kn.satz.riemannflaeche_logueberlagerung, die Vollständigkeit der Liste aus kn.satz.verzweigungsorte_charakterisierung (Bericht 06, §1).
   - `beweis.text` anhängen: „Vollständigkeit und Erreichbarkeit für alle n: kn.satz.verzweigungsorte_charakterisierung.“
   - `pruefung.f5_neue_frage`: „Liegen die Verzweigungspunkte dicht in ℂ?“ → „Dicht: ja (kn.satz.verzweigung_dicht). Offen: Gibt es Lösungen von sexp(w) ∈ 2πiℤ∖{0} im Streifen −2 < Re w ≤ −1 (kn.frage.verzweigung_abstand)?“
   - `bewertung.richtigkeit`: „Skizze; Beispielpunkte auf 50 Stellen“ → „vollständig (mit kn.satz.verzweigungsorte_charakterisierung); Beispielpunkte auf 50 Stellen“.
   - `kanten` ergänzen: `{typ: verwandt, ziel: kn.satz.verzweigung_dicht}`.
2. **kn.satz.wertevorrat**, `pruefung.f5_neue_frage`: ergänzen „; die ganze Fläche: kn.satz.riemannflaeche_logueberlagerung“.
3. Mit Bericht 07 (slog, parallel eingetragen in `r25_07_kn.yaml`): eine Kante `kn.satz.slog_verzweigung_dicht` —analog_zu→ `kn.satz.verzweigung_dicht` wäre passend (nicht von mir gesetzt, da fremde Datei).

## Stränge

Bestehender Strang **`kn.singularitaeten`** („Schnitte und Taylor-Koeffizienten von sexp“): Mitglieder um die sieben ids oben ergänzen; `verlauf` neu:
`Der Logarithmus bei −2 ist hebbar bis −3; daraus folgen die Taylor-Vorzeichen (c₂₂ < 0), der Wertevorrat und die ganze Riemannsche Fläche: eine volle Log-Überlagerung, deren Verzweigungsorte genau die Punkte sind, deren Turm auf den einer negativen ganzen Zahl fällt — dicht in ℂ, aber abzählbar.`
(Falls der Strang zu lang wird: eigener Strang `kn.riemannflaeche`, Name „Die Riemannsche Fläche von sexp“, gleiche Zeile ab „die ganze Riemannsche Fläche“.)

## Ideenkarte

Unteridee **`fortsetzung.kneser`** (enthält `kn.singularitaeten`) — nichts Neues nötig; bei eigenem Strang `kn.riemannflaeche` dort in `straenge` aufnehmen.

## Quellen (Vorschlag, nicht angelegt)

`bergweiler1993`: W. Bergweiler, *Iteration of meromorphic functions*, Bull. Amer. Math. Soc. 29 (1993), 151–188 — Rückwärtsbahnen nicht ausgenommener Punkte liegen dicht in der Julia-Menge (im Bericht „Thm. 4“, nicht geprüft). Benutzt in kn.satz.verzweigung_dicht (bisher nur im Text genannt).
