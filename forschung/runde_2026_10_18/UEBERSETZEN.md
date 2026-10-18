# Übersetzen der Runde 18.10.2026 ins Englische

Repo `D:\Users\janis\Documents\Matura\tetration-graph`, Python `C:/Users/janis/miniconda3/python.exe` (vom Repo-Wurzelordner aus).

- Werkzeug und Regeln: Docstring von `wissen/uebersetzung.py` und `wissen/uebersetzung_struktur.py` lesen; Stil an bestehenden Einträgen in `daten/uebersetzung/kn.yaml` abschauen.
- Export: `python -m wissen.uebersetzung export --gebiet <g> --felder <kern|beweis|intern|zusatz|theorie> --nur-fehlend -o <datei>`; Struktur: `python -m wissen.uebersetzung struktur-export --art <strang|idee|quelle|kapitel> --nur-fehlend -o <datei>` (für geänderte deutsche Texte ohne `--nur-fehlend` exportieren und auf die betroffenen ids filtern).
- Alle Feldgruppen eines Eintrags in **einem** Stapel zusammenführen. `einlesen` fügt neue Gruppen zu bestehenden Einträgen hinzu; `--ersetzen` nur für veraltete Einträge, deren Deutsch sich geändert hat.
- Stapel (JSON) **mit dem Write-Werkzeug** schreiben, nie über Bash-Heredocs (Backslashes gehen verloren).
- Regeln: jede $…$-Formel byte-gleich; `tex` in `formen` nie übersetzen; ids und [[Links]] unverändert; Zahlen ausserhalb von Formeln genau wie im Deutschen (Ziffern bleiben Ziffern, Wörter bleiben Wörter); Signalwörter erhalten („nur“ → only, „nicht/kein“ → Verneinung, „bedingt“ → conditional(ly)); amerikanische Schreibweise; klares, einfaches Englisch; deutsche Wörter in `quellen.stelle` übersetzen, Dateipfade unverändert; Titel fremder Quellen nicht übersetzen.
- Import: `python -m wissen.uebersetzung einlesen <datei> --modell sonnet [--ersetzen]`, `python -m wissen.uebersetzung struktur-einlesen <datei> --modell sonnet`. Jede Warnung beheben.
- Fertig, wenn `python baue.py --nur-pruefen` (≈2 min, langer timeout) für die eigenen ids keine Fehler/Warnungen meldet und `python -m wissen.uebersetzung stand` für die eigenen Gebiete nichts Neues fehlen lässt. `git diff --stat daten/uebersetzung*` darf keine ganzen Dateien neu schreiben (LF-Zeilenenden).
- Nur die eigenen Dateien anfassen (siehe Auftrag). Kein commit/push. Prozesse nie nach Namen beenden. Bericht ≤ 6 Zeilen.
