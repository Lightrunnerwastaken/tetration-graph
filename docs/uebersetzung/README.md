# Englische Fassung — Arbeitsanleitung

Stand: Werkzeuge und Prüfregeln stehen (12.10.2026); übersetzt ist noch nichts (der Pilot in `pilot/` ist nicht eingelesen).

## Ablauf je Gebiet und Feldgruppe

```
python -m wissen.uebersetzung stand                                    # Fortschritt je Gebiet
python -m wissen.uebersetzung export --gebiet fp --felder kern --max 20 -o stapel.json
#   → Übersetzer (Sonnet) bekommt GLOSSAR.md + stapel.json und liefert dieselbe Struktur mit "en" statt "de"
python -m wissen.uebersetzung einlesen uebersetzt.json --modell sonnet          # prüft, stempelt, schreibt daten/uebersetzung/fp.yaml
python -m wissen.uebersetzung einlesen uebersetzt.json --modell sonnet --trocken  # nur prüfen
python baue.py --nur-pruefen                                           # volle Prüfung inkl. Übersetzungsregeln
```

Feldgruppen (`--felder`): `kern` (titel, aussage, bereich, einfach, formen[].notiz) · `beweis` (beweis.text, gegenbeispielsuche, naechste_frage) ·
`intern` (pruefung, bewertung, voraussetzungen) · `theorie` (worum, genau, warum, beispiel, einordnung; `bilder` bleiben). `--nur-fehlend` lässt
Fertiges aus. Ein veralteter Eintrag wird nur mit `--ersetzen` überschrieben.

Format des Übersetzer-Stapels (Eingabe): `[{"id", "art", "status", "de": {…}}]`; Ausgabe: `[{"id", "en": {…}}]` mit **derselben Struktur** wie `de`
(bei `theorie`: `de_theorie` → `en_theorie`). Formeln, Zahlen, `[[Verweise]]`, Kennungen bleiben Zeichen für Zeichen gleich.

## Datei `daten/uebersetzung/<gebiet>.yaml` (maschinell geführt)

Pro Aussage: `id`, `quelle_hash` (Hash der übersetzten deutschen Texte), optional `theorie_hash`, `herkunft` (Modell, Datum), `geprueft`
(`nein` | `zweitpruefung` | `mensch`), `en`, `en_theorie`. Einzelne Texte dürfen von Hand verbessert werden; danach `einlesen --ersetzen` mit dem
Gesamteintrag oder den Hash mitführen (sonst meldet die Prüfung „veraltet“).

## Prüfregeln (`wissen/uebersetzung.py`, Tests `tests/test_uebersetzung.py`)

| Regel | Stufe | Bedeutung |
|---|---|---|
| `uebersetzung` | FEHLER | falsche Gestalt: unbekanntes Feld, Listenlänge ≠ Deutsch, `tex` übersetzt, leerer Text, Eintrag für unbekannte/doppelte Aussage, `en` zugleich in Aussagedatei |
| `uebersetzung-invariante` | FEHLER | Formeln (ausser Prosa in `\text{…}`), `[[Verweise]]` oder Kennungen verändert |
| `uebersetzung-veraltet` | WARNUNG | deutscher Text seit der Übersetzung geändert, oder kein Hash |
| `uebersetzung-zahlen` | WARNUNG | Zahlen ausserhalb von Formeln verändert |
| `uebersetzung-wortwahl` | WARNUNG | Signalwort (nicht/nur/genau/jede/vermutlich/bedingt/Vermutung/widerlegt) ohne englisches Gegenstück; britische Schreibweise; Deutsch im englischen Text |
| `uebersetzung-absaetze` | WARNUNG | Absatzzahl verändert |
| `uebersetzung-fortschritt` | HINWEIS | je Gebiet: übersetzt / nicht aktuell / gegengelesen / Theorie-Blöcke |

Die Regeln fangen formale Fehler, **keine Sinnfehler** (Pilot: Haiku übersetzte „Messbare Folge“ als „testable sequence“). Darum zusätzlich:
unabhängige Zweitübersetzung (Haiku) und Abgleich der Abweichungen, Stichproben, alle Einstiegstexte und Kapitel 0 von einem Menschen.

## Pilot (12.10.2026)

15 Aussagen und 2 Theorie-Blöcke, von Sonnet, Haiku und mir unabhängig übersetzt (`pilot/`). Alle drei bestehen die formalen Prüfungen.
Sonnet: kein Sinnfehler, „Blätter“ richtig als *petals* (meine Referenz hatte *sheets*). Haiku: ein Sinnfehler (*Folge*), vier Stilmängel
(britisch *neighbouring*, ein eingefügtes *so*, „below <“, *lower* statt *deeper*). Empfehlung: Sonnet für den Hauptteil, Haiku nur für kurze
nicht-mathematische Texte und als Zweitübersetzer. Das Glossar wurde aus dem Pilot ergänzt (time-1, petals, level, domain, consequence/sequence …).
