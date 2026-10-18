# Forschungsrunde 18.10.2026 — gemeinsame Regeln für alle 10 Agenten

Repo: `D:\Users\janis\Documents\Matura\tetration-graph`. Python: `C:/Users/janis/miniconda3/python.exe` (mpmath, sympy, python-flint, numpy, scipy).

## Ziel
So viel echte Forschung wie möglich im eigenen Gebiet: neue Identitäten und Darstellungen, versteckte Strukturen, neue Muster, Vermutungen (mit Belegen),
Beweise (wo möglich). Auch erlaubt und erwünscht: ein neues, bedeutendes Teilgebiet entdecken und anfangen.

## Vor dem Start
- `README.md` und `docs/REPOSITORY_GUIDE.md` kurz lesen.
- Prüfen, was der Graph schon enthält (grep in `daten/aussagen/*.yaml`, `daten/ideen.yaml`), damit nichts doppelt gemacht wird. Bauen statt wiederholen.
- Relevante Vorarbeit: `forschung/identitaeten_2026_10_14/` (Etappen 1–4 + luecken: θ-Welle, Mittelwertformeln, reelle Kennzeichnung), `forschung/eta_hyperoperationen/` (höhere Ränge, halber Rang).

## Rechnen
- Kneser-Werte: gp-Motor (fatou.gp-Fork). Zugänge: `experimente/gemeinsam.py` (`kneser(...)`) **nur lesend importieren**, oder `forschung/identitaeten_2026_10_14/etappe3/lib3.py` / `k3.py` als Vorlage **in den eigenen Ordner kopieren** und den Zwischenspeicher auf den eigenen Ordner umstellen. Nie in fremde Zwischenspeicher schreiben (zehn Agenten laufen gleichzeitig; parallele Schreibzugriffe zerstören Dateien). Für `gemeinsam.kneser` gilt: es schreibt in `experimente/daten/kneser_cache.json` — daher nicht benutzen, ausser man kopiert die Funktion mit eigenem Speicherpfad.
- Aktueller Motor mit Fix für dps ≥ 60 bei |Im z| > 1/4: `D:\Users\janis\Documents\Tetration\gp-tetration\src\fatou_backend\vendor\fatou_fork.gp` (nur lesen). Höhere Ränge: dort gibt es `hyper(r, x)` (Kneser für gerade Ränge); nur lesen.
- **Laptop-Schonung:** zehn Agenten teilen sich einen Laptop. Jeder Rechenlauf höchstens 10 Minuten (timeout), höchstens ein schwerer Prozess gleichzeitig pro Agent, lieber kleine Läufe. Keine Läufe im Hintergrund vergessen.
- **Prozesse:** nie nach Namen beenden (andere Sitzungen rechnen mit gp.exe). Nur eigene Prozesse über die Prozessnummer (taskkill /F /T /PID <pid>).
- Keine Dateien herunterladen (Webseiten lesen ist erlaubt, z. B. OEIS, arXiv-Abstracts, Tetration Forum).

## Grenzen
- Nicht ändern: `daten/`, `experimente/`, `wissen/`, `ansicht/`, gp-tetration. Kein git commit, kein push.
- Nur im eigenen Ordner `forschung/runde_2026_10_18/<nr>_<name>/` schreiben.

## Strenge
- Jede Aussage klar einordnen: **bewiesen** / **hergeleitet mit Lücke** (Lücke benennen) / **numerisch** (Genauigkeit + Gegenprobe) / **Vermutung** (Belege, Gegenbeispielsuche).
- Numerische Behauptungen mit zwei unabhängigen Einstellungen prüfen. Beinahe-Treffer bei Konstanten erst glauben, wenn sie 15 weitere Stellen halten.
- Ehrlich über Fehlschläge berichten; negative Ergebnisse sind Ergebnisse.

## Abgabe
`forschung/runde_2026_10_18/<nr>_<name>/BERICHT.md` auf Deutsch: oben eine kurze Zusammenfassung in einfacher Sprache, dann eine Tabelle aller Befunde
(Aussage | Status | Beleg/Genauigkeit | Neuheit-Vermutung), dann Beweise/Details, dann „Vorschläge für Grapheinträge“ (Titel + Art + Status, je eine Zeile) und offene Fragen.
Skripte so ablegen, dass jede Zahl reproduzierbar ist. Abschlussnachricht: höchstens 15 Zeilen.
