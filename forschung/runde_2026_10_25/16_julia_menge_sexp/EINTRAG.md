# Eintrag Bericht 16 (Julia-Menge von eᶻ unter sexp) in den Graphen

Stand 25.10.2026 17:55. Dateien: `daten/aussagen/r25_16_kn.yaml`, `daten/aussagen/r25_16_fp.yaml`; Theorie: `daten/theorie/r25_16_kn.yaml`, `daten/theorie/r25_16_fp.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler in diesen Dateien (nur die erwarteten Hinweise „kein Kapitel / keine Idee / kein Strang“).

## Angelegte ids
- kn.satz.flucht_urbild_dicht | E = sexp⁻¹(I(exp)) 1-invariant, dicht, Fläche null, enthält ℝ (B1, B2) — bewiesen
- kn.satz.flucht_gitter | Fluchtmuster oben Λ-periodisch, Λ = ℤ + (2πi/L)ℤ; N = 4.2034 Im z − Re z + periodisch, N(z+ω) = N(z) + 4 (B3–B5) — bewiesen
- kn.beob.flucht_gitter_numerik | Anteil N(z+ω) = N(z)+4 je Höhenband, 1.0 ab Im z ≥ 3; Steigung 4.244 (B4, B5) — numerisch
- kn.beob.streifen_linse | Grundstreifen → Linse in |Im w| < 1.756; senkrechte Linien Bögen L̄→L, keine Kreise, Spiralen bei L (B6, B11) — numerisch
- kn.beob.reelle_urbilder | Im F > 0 auf dem oberen Streifen (≥ 0.178 Im z); erste nichtreelle Urbilder von ℝ bei Re z ≈ 2.52 (B7, B8) — numerisch
- fp.satz.kein_haar_bei_L | kein Haar von eᶻ hat L oder L̄ im Abschluss (B9) — skizze
- kn.kor.haar_urbilder_beschraenkt | Haar-Urbilder in senkrechten Streifen nach oben beschränkt (B10) — skizze
- kn.frage.flucht_torus | Dimension, Torus-Bild, Linse, {T = ∞}, andere Basen (Offene Fragen 1–5) — offen

Nicht eingetragen: Nebenbefund „θ an jeder ganzen Zahl unbeschränkt“ — steht schon in kn.satz.theta_singularitaet (ii).

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1, B2 | Beweis gelesen. (i) aus vollständiger Invarianz von I; Dichte über J(exp) = ℂ (Rückwärtsbahn von 1 dicht, F offen); Mass 0 über fp.satz.messbare_dynamik (Rezensionen) und lokale Biholomorphie | bestanden → bewiesen, Neuheit folklore; Mass-Aussage ausdrücklich an die Literatur gebunden |
| B3, B4/B5 (Formeln) | Beweis Zeile für Zeile: log L = L (Hauptzweig, Im L = 1.337 < π), F_L(u+1) = exp F_L(u), F_L(u+P) = F_L(u), Re(LP) = 0, ω = P − 4; P = 4.446951 + 1.057940i stimmt mit γ aus kn.satz.wertevorrat | bestanden → bewiesen; die Kneser-Nähe als „ausser nahe den Niveaugrenzen“ formuliert |
| B4 (Zahlen) | `e3_gegenprobe.py` neu gerechnet (R = 20, Gitter 150×900) | 0.8006 / 0.9463 / 0.9975 / 0.9978 / 0.9999 / 1.0 — gleich wie im Bericht |
| B6 | nicht neu gerechnet (e1 überschreibt Bilder); Wert 1.7558 aus `cache/e1_zahlen.json` übernommen | als numerisch eingetragen |
| B7, B8 | `e2_reell.py` neu gerechnet | min Im F/Im z = 0.17831; erstes Im F = π bei 1.51625 + 0.63456i — gleich |
| B11 | Kreistest in `e3_gegenprobe.py` neu gerechnet; Spiralen aus kn.kor.doppelreihe_oben nachgeprüft (Drehung Re L, Faktor e^{−Im L} = 0.263) | 2.88 / 3.77 / 4.78 / 5.97 % — gleich; „bewiesen (asymptotisch)“ nur für die Spirale, Rest numerisch |
| B9 | Beweis gelesen: (a) Strahl ohne Kontakt zu Im = (2k+1)π liegt in einem Streifen, L im Innern von S₀ erzwingt s = 0^∞; (b) schliesslich-0-Adressen landen in ℝ. Schlüssig, aber nur unter den Struktursätzen von Schleicher–Zimmer für den Fall „singulärer Wert flüchtet auf ℝ“ | als skizze eingetragen (Lücke benannt), Neuheit folklore |
| B10 | folgt aus B9 und F → L gleichmässig in senkrechten Streifen | skizze (erbt die Lücke) |

## Änderungen am Bestand
Keine nötig. (Optional: kn.satz.kneser_eigenschaften, Feld `aussage` — „für jedes feste x > −2“ ergänzen um „gleichmässig in jedem Streifen a ≤ Re z ≤ b“; Begründung: folgt aus F = F_L∘h mit θ → a₀ und wird in kn.kor.haar_urbilder_beschraenkt benutzt.)

## Quellen, die im Bestand fehlen (nicht angelegt)
- L. Rempe, The escaping set of the exponential, Ergodic Theory Dynam. Systems 30 (2010) — I(exp) zusammenhängend (in fp.satz.kein_haar_bei_L im Text genannt). Angaben aus dem Gedächtnis.
- R. L. Devaney, Knaster-like continua and complex dynamics, ETDS 13 (1993); L. Rempe, A landing theorem for periodic rays of exponential maps, Proc. AMS 134 (2006) — für die Neuheitsfrage von B9. Angaben aus dem Gedächtnis.

## Vorschlag Stränge
- Neuer Strang `kn.flucht` (gebiet kn), Name „Flucht und Haare unter sexp“, Mitglieder: kn.satz.flucht_urbild_dicht, kn.satz.flucht_gitter, kn.beob.flucht_gitter_numerik,
  kn.beob.streifen_linse, kn.beob.reelle_urbilder, kn.kor.haar_urbilder_beschraenkt, kn.frage.flucht_torus.
  `verlauf`: „sexp zieht die Fluchtmenge von eᶻ zurück: dicht, Fläche null, oben doppelt periodisch mit dem Gitter ℤ + (2πi/L)ℤ und vier Schritten mehr pro schrägem Schritt ω; der Grundstreifen landet in einer schmalen Linse, und weil kein Haar L erreicht, steigen die Haar-Urbilder nicht ins Unendliche.“
- fp.satz.kein_haar_bei_L in den bestehenden Strang `fp.fixpunkte` (verlauf unverändert).

## Vorschlag Ideenkarte
- Neue Unteridee `dynamik.flucht` unter `dynamik`, Name „Fluchtmenge und Haare“, `kurz`: „Welche Stellen die Tetration in die Flucht schickt — ein schräges Gitter in Kneser-Koordinaten — und warum kein Haar von eᶻ den Fixpunkt L erreicht.“, straenge [kn.flucht].
