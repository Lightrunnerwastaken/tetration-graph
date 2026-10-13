# Literaturkarte — alle Gebiete so genau kartieren wie die eigene Forschung (ab 12.10.2026)

Auftrag des Autors: Die Bibliothek enthält viel eigene Forschung, aber die Arbeiten anderer nur dünn. Jedes Gebiet soll so genau kartiert
werden wie unsere eigene Arbeit: wer woran gearbeitet hat, mit welcher Methode und welchem Ergebnis, und wo mehrere Leute dieselbe Frage
verschieden angegangen sind.

Stand vorher: 321 Aussagen, davon 158 nur mit eigenen Quellen, 80 nur mit fremden, 77 gemischt; 75 Literatur- und 55 Forumsquellen.

## Stränge und Zuständigkeit (jede Datei hat genau einen Besitzer)

| Strang | Gebiete | Eigene Dateien (nur dieser Strang ändert sie) | Experiment-Nummern |
|---|---|---|---|
| L1 Iteration, Funktionalgleichungen, Lösungsklassen, Fixpunkte (Agent) | it, fg, kl, fp | daten/aussagen/01_, 02_, 03_, 04_, 09_, 10_*.yaml; daten/theorie/it, fg, kl, fp.yaml; daten/straenge/it, fg, kl, fp.yaml | e60–e69 |
| L2 Kneser, slog, Hyperoperationen, Basiswechsel (Agent) | kn, bw | daten/aussagen/05_, 06_*.yaml; daten/theorie/kn, bw.yaml; daten/straenge/kn, bw.yaml; **einziger mit Kneser-Motor** | e70–e79 |
| L3 Grundlagen, Arithmetik, Brücken, Gesamtverzeichnis (Leitung) | gr, ar, br | daten/aussagen/00_, 07_, 08_*.yaml; daten/theorie/gr, ar, br.yaml; daten/straenge/gr, ar, br.yaml; GESAMTVERZEICHNIS.md | e80–e89 |

Neue Aussagen dürfen auch in neuen Dateien `daten/aussagen/lit_<gebiet>_*.yaml` stehen (Kopf `gebiet: <gebiet>` und `aussagen:`), neue
Theorie-Blöcke in `daten/theorie/lit_<gebiet>.yaml`, neue Stränge in `daten/straenge/lit_<gebiet>.yaml`.

**Gemeinsame Dateien:** `daten/quellen.yaml`, `daten/objekte.yaml`, `daten/experimente.yaml`, `daten/lehrbuch.yaml`, `daten/gebiete.yaml` ändert
nur die Leitung. Neue Quellen schreibt jeder Strang in `daten/quellen/lit_L<n>.yaml`, neue Objekte in `daten/objekte/lit_L<n>.yaml`, neue
Experimente in `daten/experimente/lit_L<n>.yaml` (der Lader führt alles zusammen; ein Schlüssel in zwei Dateien ist ein Fehler — vorher prüfen,
ob die Quelle schon existiert, und dann den vorhandenen Schlüssel benutzen). Änderungswünsche an fremden Dateien: Eintrag im BRETT.

## Arbeitsdateien

- `L<n>_.../KARTE_<gebiet>.md` — die Literaturkarte (siehe Auftrag), `L<n>_.../BERICHT.md` — Stand, Entscheidungen, Lücken.
- `BRETT.md` — gemeinsames Brett (Einträge oben unter der Kopfzeile).
- `GESAMTVERZEICHNIS.md` — Gesamtliste aller Tetrations-Literatur (L5), mit Zuordnung zu den Strängen.
