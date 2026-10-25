# Forschungsrunde 25.10.2026 — gemeinsame Regeln für alle 70 Agenten

Repo: `D:\Users\janis\Documents\Matura\tetration-graph`. Python: `C:/Users/janis/miniconda3/python.exe` (mpmath, sympy, python-flint, numpy, scipy).

## Ziel und Zeit
- So viel echte Forschung wie möglich im eigenen Gebiet: neue Identitäten und Darstellungen, versteckte Strukturen, neue Muster, Vermutungen mit Belegen,
  Beweise wo möglich. Ein neues, bedeutendes Teilgebiet zu entdecken ist ausdrücklich erwünscht.
- **Zeitbudget etwa 30 Minuten.** Beim Start mit `date` die Uhrzeit notieren, zwischendurch prüfen, nach etwa 30 Minuten abschliessen und den Bericht schreiben.
  Lieber eine Sache sauber als fünf halb.
- **Starte keine weiteren Agenten.** Arbeite allein.

## Vor dem Start (kurz, höchstens 5 Minuten)
- Was der Graph schon weiss: grep in `daten/aussagen/*.yaml` und `daten/ideen.yaml` nach den Stichworten deines Gebiets. Bauen statt wiederholen.
- Vorarbeit: `forschung/runde_2026_10_18/*/BERICHT.md` (10 Berichte der letzten Runde), `forschung/identitaeten_2026_10_14/` (θ-Welle, Mittelwertformeln,
  reelle Kennzeichnung), `forschung/eta_hyperoperationen/` (höhere Ränge, halber Rang). Nur lesen, was dein Gebiet betrifft.

## Rechnen
- Kneser-Werte: gp-Motor (fatou.gp-Fork in `D:\Users\janis\Documents\Tetration\gp-tetration\src\fatou_backend\vendor\fatou_fork.gp`, nur lesen; dort gibt es
  auch `hyper(r, x)` für höhere Ränge). Vorlage für den Zugang: `forschung/runde_2026_10_18/03_basis_ebene/motor.py` **in den eigenen Ordner kopieren** —
  sie legt den Zwischenspeicher im eigenen Ordner an. Nie in fremde Zwischenspeicher schreiben; `experimente/gemeinsam.py` nicht benutzen (schreibt in einen geteilten Speicher).
- **Laptop-Schonung:** bis zu 16 Agenten laufen gleichzeitig. Jeder Rechenlauf höchstens 5 Minuten (timeout), nur ein Rechenprozess gleichzeitig,
  kleine Läufe bevorzugen. Nichts im Hintergrund vergessen. Achtung: `timeout` in Git-Bash beendet Windows-Prozesse nicht zuverlässig — eigene hängende Prozesse über die Prozessnummer beenden.
- **Prozesse:** nie nach Namen beenden (andere Sitzungen rechnen mit gp.exe). Nur eigene über die Prozessnummer (`taskkill /F /T /PID <pid>`).
- Keine Dateien herunterladen (Webseiten lesen ist erlaubt: OEIS, arXiv-Abstracts, Tetration Forum).

## Grenzen
- Nicht ändern: `daten/`, `experimente/`, `wissen/`, `ansicht/`, gp-tetration. Kein git commit, kein push.
- Nur im eigenen Ordner `forschung/runde_2026_10_25/<nr>_<name>/` schreiben. Grosse Zwischenspeicher in `<ordner>/cache/`.
- YAML/Markdown mit Backslashes nur über das Write-Werkzeug schreiben, nicht über Bash-Heredocs.

## Strenge
- Jede Aussage einordnen: **bewiesen** / **hergeleitet mit Lücke** (Lücke benennen) / **numerisch** (Genauigkeit + Gegenprobe mit zweiter Einstellung) /
  **Vermutung** (Belege + Gegenbeispielsuche) / **negativ** (was nicht geht, und warum).
- Beinahe-Treffer bei Konstanten erst glauben, wenn sie 15 weitere Stellen halten. Ehrlich über Fehlschläge berichten.

## Abgabe
`forschung/runde_2026_10_25/<nr>_<name>/BERICHT.md` auf Deutsch: oben 3–6 Sätze in einfacher Sprache, dann Tabelle aller Befunde
(Aussage | Status | Beleg/Genauigkeit | Neuheit-Vermutung), dann Beweise/Details, dann „Vorschläge für Grapheinträge“ (Titel + Art + Status, je eine Zeile)
und offene Fragen. Skripte so ablegen, dass jede Zahl reproduzierbar ist.
