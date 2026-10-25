# Eintrag Bericht 01 (Halbebenen-Eindeutigkeit Q*) — 25.10.2026

Dateien: `daten/aussagen/r25_01_br.yaml`, `daten/theorie/r25_01_br.yaml`. `baue.py --nur-pruefen`: 0 Fehler (nur die erwarteten Hinweise).

## Angelegte Aussagen

| id | Titel kurz |
|---|---|
| br.satz.halbebene_liouville_scharf | Q* gilt, wenn p_T oben eindeutig ist und o(e^{2π Im z}) wächst; Schwelle scharf (bewiesen) |
| br.satz.schroeder_kugel_genuegt | Bleibt T in einem Halbstreifen nahe L, so ist T Kneser; keine Zusatzwindungen um L (bewiesen, über Paulsen–Cowgill) |
| br.satz.enden_klassifikation | Mögliche Enden bei i∞: Exponent log L ± 2πik bzw. log L̄ + 2πik; nur L, k = 0 kommt vor (skizze, 3 Lücken) |
| br.satz.verbiegung_schwarz | Verbiegung mit σ(U) ⊂ U ist trivial (bewiesen) |
| br.lem.zungen_bis_unendlich | Komponenten von {Im σ < 0} reichen bis i∞, σ am Rand beidseitig unbeschränkt (bewiesen) |
| br.beob.hauptblatt_singularitaet_trig | 48 von 48 trigonometrischen Verbiegungen: Log-Singularität auf dem Hauptblatt in H (numerisch) |
| br.verm.verbiegung_hauptblatt | Vermutung: das gilt für jedes nicht konstante trigonometrische Polynom |
| br.kor.nebenblatt_minus2_regulaer | Nach Kreuzen von (−4, −3) ist −2 regulär, Zählen allein beweist nichts (bewiesen) |
| br.frage.vertauschtes_kneser | Gibt es eine Tetration in Hol(H) mit Grenzwert L̄ oben? (offen) |

Befund 5 ist als Folgerung in br.satz.schroeder_kugel_genuegt eingearbeitet, Befund 2 in br.satz.halbebene_liouville_scharf, Befund 6 in br.satz.enden_klassifikation, Befund 12 (negativ) in der Frage und der Klassifikation.

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1 Satz A | Beweis Zeile für Zeile nachgerechnet (Koeffizientenformel, \|Re f\| ≤ 2A − Re f, Grenzübergang r → 0, Spiegelung) | richtig → bewiesen. Bemerkung: die Holomorphie von T auf H wird gar nicht benutzt. |
| B2 Schärfe | direkte Rechnung | richtig; „scharf“ nur als Aussage über periodische Funktionen (das Beispiel selbst liegt nicht in Hol(H), Befund 9) — so eingetragen |
| B3 Satz A′ | Beweis nachgelesen (Schröder-Gleichung auf B, Periodizität von E, \|E\| ≤ C\|q\|^{−0.2128}, Index-Argument links, c ≤ −1 nötig) | richtig → bewiesen (stützt sich auf Paulsen–Cowgill) |
| B4/B6 Enden | nachgelesen | Kern richtig; zusätzlich zur Lücke des Berichts eine zweite gefunden: beim Polfall ist „T hat keinen Grenzwert“ nicht bewiesen (σ(T) → ∞ längs einer Spirale, Ψ könnte dort einen asymptotischen Wert haben) → skizze mit 3 Lücken |
| B5 | Folgerung aus A′ geprüft (Kneser hat E → κ ≠ 0, also k = 0) | richtig |
| B7 Satz B | nachgelesen (Windungszahl 1, hebbare Null, Randbetrag 1, Max-/Min-Prinzip) | richtig → bewiesen |
| B8 Lemma C | nachgelesen | richtig; (b) braucht σ′ > 0 auf ℝ (damit die an ℝ anliegende Komponente einen Streifen enthält) — als Voraussetzung ergänzt |
| Folgerung „Gegenbeispiel muss auf Nebenblätter ausweichen“ | geprüft | folgt aus B + C **nicht** zwingend (Verschiebung um ganze Zahlen ändert σ-Werte); nur „σ verlässt U“ als Satz eingetragen, der Rest als Vermutung/Frage |
| B9, B10 | `hauptblatt.py 120 7` erneut gelaufen (52 s): die 5 Treffer z* und „43 mit Hauptblatt-Singularität, 0 ohne“ identisch. Zusätzlich eigener Code für ε = 0.1 (mpmath findroot dps 40, 20000 Wegpunkte, Bisektion der Vorzeichenwechsel von Im σ) | z* = −1.3444755008890480019 + 0.43879214577272809932i, \|σ+2\| = 3·10⁻⁴³; Im σ wechselt längs des Wegs nie das Vorzeichen (Bild ganz oben), Abstand zum Schnitt ≥ 0.019 → bestätigt |
| B11 | lokale Rechnung mit kn.satz.minus2_hebbar nachvollzogen | richtig → bewiesen (Korollar) |

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **br.frage.halbebene_eindeutigkeit**, Feld `aussage`:
   - Satz „Ansatz: θ_T = h∘T − z ist dann auf der oberen Halbebene bis auf harmlose Log-Punkte (T = L) periodisch holomorph“ →
     „… bis auf Log-Punkte (T = L) und Stellen, an denen T singuläre Werte von σ (vermutlich 0, 1, e, eᵉ, …) umläuft, periodisch holomorph (br.satz.enden_klassifikation)“.
     Grund: Befund 4 — die Mehrdeutigkeit kommt nicht nur von T = L.
   - Satz „nicht ausgeschlossen ist, dass alle auf Nebenblättern liegen“ → „für 48 trigonometrische Verbiegungen liegt ein Verzweigungspunkt nachweislich auf dem Hauptblatt
     (br.beob.hauptblatt_singularitaet_trig, br.verm.verbiegung_hauptblatt)“.
   - Anhängen: „**Stand 25.10.2026:** Die fehlende Bedingung bei i∞ ist genau ‚p_T eindeutig und o(e^{2π Im z})‘ (br.satz.halbebene_liouville_scharf); statt des Grenzwerts genügt ‚nahe L
     in einem Halbstreifen‘ (br.satz.schroeder_kugel_genuegt); jede echte Verbiegung verlässt mit σ die obere Halbebene (br.satz.verbiegung_schwarz); die möglichen Enden sind klassifiziert,
     natürlichster Kandidat für ein Gegenbeispiel ist ein vertauschtes Kneser (br.frage.vertauschtes_kneser).“ Status bleibt `offen`.
   - `kanten`: `{typ: verwandt, ziel: br.satz.halbebene_liouville_scharf}`, `{typ: verwandt, ziel: br.frage.vertauschtes_kneser}` ergänzen.
2. **br.beob.verbiegungen_singularitaeten**, `pruefung.f5_neue_frage`: „Gibt es eine Verbiegung, deren Verzweigungspunkte alle auf Nebenblättern liegen?“ → anhängen:
   „Für trigonometrische Polynome nein (48 Fälle, br.beob.hauptblatt_singularitaet_trig); allgemein br.verm.verbiegung_hauptblatt.“

## Stränge

Neuer Strang **br.halbebene** „Q*: erzwingt die Halbebene schon Kneser?“, `verlauf`: „Was bei i∞ fehlt, ist genau benannt (p_T eindeutig und langsamer als e^{2π Im z});
echte Verbiegungen müssen die obere Halbebene verlassen, ihre Enden sind klassifiziert — offen bleibt vor allem ein ‚vertauschtes Kneser‘.“
Mitglieder: die neun ids oben. (Alternative: in br.reell aufnehmen, wo br.frage.halbebene_eindeutigkeit steht.)

## Ideenkarte

Unteridee **fortsetzung.reell** („Reelle Kennzeichnung?“), `straenge` um br.halbebene ergänzen (Q* ist genau die offene Stelle der lokalen Kennzeichnung).
Alternative: fortsetzung.eindeutigkeit (Q* ist Paulsen–Cowgill ohne Grenzwertbedingung).
