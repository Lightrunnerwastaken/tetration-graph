# Regeln für das Eintragen der Runde 25.10.2026 in den Graphen

Repo `D:\Users\janis\Documents\Matura\tetration-graph`, Python `C:/Users/janis/miniconda3/python.exe`.

## Aufgabe
Die Befunde eines Berichts `forschung/runde_2026_10_25/<nr>_*/BERICHT.md` als Aussagen in den Graphen übertragen — nur die substanziellen (bewiesen,
hergeleitet, gut belegte Beobachtungen und klare Vermutungen/Fragen), meist 5–9 je Bericht. Negative Ergebnisse, die etwas klären, dürfen rein (als `gescheitert` oder in eine Aussage eingearbeitet).

## Format (genau wie bestehende Dateien; Vorbilder: `daten/aussagen/id_kn_1.yaml`, `id_br_2.yaml`, `id_gr_1.yaml`, `lit_kn_3.yaml`)
- Eine Datei pro Gebiet und Bericht: `daten/aussagen/r25_<nr>_<gebiet>.yaml` mit Kopf `gebiet: <gebiet>` und `aussagen:`. Gebiete: gr, fg, it, fp, kl, kn, bw, ar, br (Bedeutung in `daten/gebiete.yaml`).
- Pflichtfelder und erlaubte Werte: `wissen/schema.py` (ARTEN, EBENEN, STATUS, NEUHEIT, NUTZEN, KANTEN, PRUEFFRAGEN). Bei Arten satz/lemma/korollar/konstruktion/gegenbeispiel/korrektur/vermutung/beobachtung die fünf Prüffragen `pruefung` ausfüllen.
- Status ehrlich: bewiesen / bedingt (+ `bedingungen_text`) / skizze (+ `beweis.luecken`) / numerisch / vermutung / offen.
- Neuheit: `neu_ungeprueft` für Neues, `folklore`/`bekannt`/`klassisch` sonst; `vermutlich_neu` nur mit dokumentierter Recherche.
- Ids: `<gebiet>.<art-kürzel>.<name>` (Kürzel wie im Bestand: def, satz, lem, kor, konstr, verm, beob, gegen, gescheitert, korr, frage). Vor Vergabe mit grep prüfen, dass die id frei ist.
- `einfach:` in einfacher Sprache (Schulniveau), `aussage:` genau, Formeln in $…$ (KaTeX; kein `\operatorname`, stattdessen `\mathrm`).
- Quellen: `{ref: diese_arbeit, stelle: 'forschung/runde_2026_10_25/<nr>_*/BERICHT.md, B<k>; <skript>.py'}`; Literatur/Forum nur mit bestehenden ids aus `daten/quellen.yaml`/`daten/quellen/*.yaml` (keine neuen Quellen anlegen — falls nötig, im Vorschlag nennen).
- Kanten zu bestehenden Aussagen (grep nach passenden ids), z. B. stuetzt, verallgemeinert, analog_zu, verwandt, motiviert, spezialfall_von. Nur existierende Ziele.
- Zu jeder Aussage einen Theorie-Block in `daten/theorie/r25_<nr>_<gebiet>.yaml` (Felder worum, genau, warum, beispiel, bilder: [], einordnung mit [[id]]-Links; Vorbild `daten/theorie/id_kn_2.yaml`). Kurz halten.

## Technik
- YAML mit Backslashes NUR mit dem Write-/Edit-Werkzeug schreiben, nie über Bash-Heredocs (Backslashes gehen verloren). In einfach gequoteten YAML-Strings ein Apostroph als '' schreiben (also $f''$ für f′ in '…').
- Danach `python baue.py --nur-pruefen` (dauert ~2 min; mit langem timeout) und alle FEHLER in den eigenen Dateien beheben. Warnungen zu fehlender Übersetzung und Hinweise „steht in keinem Strang/Kapitel“ sind erwartet und bleiben.
- Andere Dateien NICHT ändern (keine straenge, ideen.yaml, lehrbuch.yaml, uebersetzung*, keine fremden Aussagen). Kein commit/push. Prozesse nie nach Namen beenden.

## Gegenprobe vor dem Eintragen (Pflicht)
- Die Berichte stammen von Agenten und sind nicht geprüft. Für jeden Befund, der als **bewiesen** eingetragen werden soll: Beweis kritisch lesen und
  Lücken suchen; im Zweifel als `skizze` (mit `beweis.luecken`) oder `bedingt` eintragen. Für die 1–2 wichtigsten numerischen Befunde je Bericht ein Skript
  aus dem Berichtsordner erneut laufen lassen (höchstens 5 Minuten, nur ein Prozess; Zwischenspeicher nur im Berichtsordner) und die Zahl vergleichen.
- Ergebnis der Gegenprobe in EINTRAG.md festhalten (je Befund: geprüft wie, Ergebnis). Befunde, die die Gegenprobe nicht bestehen, nicht eintragen oder herabstufen.
- **Starte keine weiteren Agenten.**

## Änderungen an bestehenden Aussagen
- Bestehende Aussagen NICHT selbst ändern. Wenn ein Befund eine bestehende Aussage verbessert, widerlegt oder deren Status ändert (z. B. „bedingt“ → „bewiesen“,
  Vermutung widerlegt, offene Frage beantwortet), dies in EINTRAG.md unter „Änderungen am Bestand“ genau vorschlagen (id, Feld, alter → neuer Inhalt, Begründung).

## Abgabe
- `forschung/runde_2026_10_25/<nr>_*/EINTRAG.md`: Liste der angelegten ids (je Zeile: id | Titel kurz), Vorschlag für Stränge (neuer Strang oder bestehender, mit Name und einer Zeile `verlauf`) und für die Ideenkarte (welche Unteridee in `daten/ideen.yaml`, oder eine neue mit Name und `kurz`).
- Abschlussnachricht ≤ 8 Zeilen.
