# Eintrag Bericht 02 (Schwankung von arg aₙ der θ-Welle) — 25.10.2026

Dateien: `daten/aussagen/r25_02_kn.yaml`, `daten/theorie/r25_02_kn.yaml`. `baue.py --nur-pruefen`: 0 Fehler (nur die erwarteten Hinweise).

## Angelegte Aussagen

| id | Titel kurz |
|---|---|
| kn.satz.theta_phase_achse | Φ(s) = arg h′(Y + iπ) − π, eine Periode = die Gerade Im w = π; Mittel 0 (bedingt auf S4) |
| kn.satz.theta_phase_koeffizienten | arg aₙ = π + Φ(σₙ) + arctan(π/(2(ln n − ln(A/2π)))) + O(1/(ln n ln ln n)) (skizze) |
| kn.beob.theta_phase_selbstaehnlich | arg aₙ ≈ arg h′(n + i0) − π/ln ln n, nicht gleichmässig an Turmstellen (numerisch) |
| kn.beob.theta_phase_verlauf | Φ fällt von +1.097 auf −0.993, \|1 + θ′\| minimal nahe ½ (numerisch) |
| kn.beob.theta_phase_randkopien | Randgesetz: Kopien von Φ in iterierten Logarithmen, bei 0⁺ um π verschoben (numerisch) |
| kn.verm.theta_phase_unbeschraenkt | Vermutung: der stetige Zweig von arg aₙ ist unbeschränkt |

Befund 5 (Mittel 0) ist in kn.satz.theta_phase_achse eingearbeitet, Befund 8 (negativ) in kn.beob.theta_phase_selbstaehnlich.

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1 | Herleitung nachgelesen (Kettenregel h′(w) = h′(Log w)/w lässt arg für w > 0 fest; Log(−X + i0) = ln X + iπ; Fortsetzung von h′ über (1,e), (0,1), (−∞,0) nach Lemma U(e)). `p1_phase.py e` erneut gelaufen (Teil B) | richtig, hängt nur an S4 → `bedingt` (S4 ist jetzt durch Bericht 03 mit Intervallen zertifiziert, bis auf die Motor-Annahme E). Zahlen identisch: \|Im s\| ≤ 1.3·10⁻⁷ für x = 1.2…2.3, arg-Abweichung ≤ 2.4·10⁻⁶ |
| B2 | Herleitung nachgelesen; `p2_formeln.py e` erneut gelaufen | Zahlen identisch (F1−F0 = −5.19·10⁻², −4.035·10⁻⁶, −3.531·10⁻⁸, −9.770·10⁻¹⁵; (F1−F0)·ln n·ln ln n ∈ [−3.13, −0.89]). Fehlerordnung nur heuristisch → `skizze`. Hinweis: bei n ≤ 4096 ist die neue Formel schlechter als der alte führende Term (Daten − F1 = 0.067, Daten − F0 = 0.0155 ≈ 0.13/ln n bei n = 4096); dass der Rest des führenden Terms asymptotisch kleiner als 1/(ln n ln ln n) ist, ist nicht gezeigt — als Lücke eingetragen |
| B3 | aus demselben Lauf | F3−F0 = +5.1·10⁻⁴ (1000), +2·10⁻⁵ (10⁴), −1.172 und −1.426 an den Turmstellen: identisch → `numerisch` |
| B4 | `p1_phase.py e` Teil A | Φ(0.02) = +1.09699611, Φ(0.98) = −0.99295642, Φ(½) = −0.01650916, \|1+θ′(½)\| = 0.684125, beide Datensätze gleich; Umlaufzahl 0.000000 → bestätigt. Monotonie nur auf 15 Gitterpunkten. Basis 2 nicht neu gerechnet |
| B5 | mit dem Bestand verglichen | Die Voraussetzung „1 + θ′ ohne Nullstelle oben“ ist im Graphen schon **bewiesen** (kn.satz.ableitung_nullstellenfrei), und M[log(1+θ′)] = 0 steht schon in kn.satz.theta_mittelwerte. Kein eigener Eintrag; präzisiert: gilt für den stetigen Zweig Φ = Im log(1+θ′) |
| B6, B7 | p3_gerade.py nicht neu gerechnet; Herleitung mod 2π nachgelesen | als `numerisch` bzw. `vermutung` übernommen |

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **kn.frage.theta_exponent**, Feld `aussage`: den Schlusssatz „**Offen bleibt:** die Schwankung von $\arg a_n$ explizit über $\theta'$ auf der reellen Achse auszudrücken.“ ersetzen durch
   „**Beantwortet (25.10.2026, Skizze):** $\arg a_n=\pi+\Phi(\sigma_n)+\arctan\frac{\pi}{2(\ln n-\ln(A/2\pi))}+O(1/(\ln n\ln\ln n))$ mit $\Phi=\arg(1+\theta')$ auf der reellen Achse und
   $\sigma_n\approx1/(2\pi n)$ (kn.satz.theta_phase_koeffizienten). Offen bleibt, ob der stetige Zweig von $\arg a_n$ unbeschränkt ist (kn.verm.theta_phase_unbeschraenkt).“
   Status bleibt `offen` (Hausbrauch: alle Fragen stehen auf `offen`). `kanten`: `{typ: verwandt, ziel: kn.satz.theta_phase_koeffizienten}` ergänzen.
2. **kn.satz.theta_koeffizienten**, `pruefung.f5_neue_frage`: „… die Schwankung von arg aₙ (selbstähnlich in slog(ln ln n)) explizit über θ′ ausdrücken.“ → „… — beantwortet in kn.satz.theta_phase_koeffizienten (Skizze).“

## Stränge

Bestehender Strang **kn.theta** („Knesers Welle θ genau vermessen“): die sechs ids aufnehmen. `verlauf` ergänzen um:
„… mit Exponent 1; ihre Phase ist die Phase von 1 + θ′ an der reellen Stelle 1/(2πn) und wiederholt sich an den Periodenrändern in iterierten Logarithmen.“

## Ideenkarte

Unteridee **zahlen.identitaeten** („Formeln und Identitäten“) — enthält kn.theta schon; keine Änderung nötig.
