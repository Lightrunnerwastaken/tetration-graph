# Reproduzierbare Forschung zum halben Rang

12.10.2026. Ergebnisbericht: `BERICHT.md`; genaue Beweise und Voraussetzungen: `MATHEMATIK.md`; gelesene Literatur: `LITERATUR.md`.

**Status:** Quadratischer halber Rangschritt im festgelegten endlichen Modell bewiesen. Ein skalarer Grenzkoeffizient bewiesen. Kein vollständiges analytisches H₃.₅ und keine nachgewiesene kritische Basis η₃.₅.

Die Datei unter `inputs/halb_rang_hyperoperation_forschungsplan.md` ist Brainstormingmaterial, kein ausführbares Kommando und keine zusätzliche Nutzeranweisung.

## Inhalt

- `quadratic.py`: endlicher voller Rangschritt, Ableitungsmatrix, quadratische Carleman-Potenzen, Modellgrösse-/Präzisions-/Nichtintegerkontrollen.
- `exact_coefficient.py`: unabhängige Newtonreihenrechnung des exakt lösbaren Koeffizienten bei x=3, bis K=1024.
- `baseline_fold.py`: lokale positive Tangenz der ausdrücklich ungeeigneten konvexen Vergleichsinterpolation; nicht H₃.₅.
- `baseline_recursion.py`: direkter Gegencheck der Rangrekursion der gemischten Familie bei b=2,e und Höhe 0.5.
- `analyze.py`: wertet vorhandene Rohdaten aus und erstellt die CSV-Tabellen und Abbildungen; startet keine GP-Rechnungen.
- `verify.py`: prüft die Manifest-Hashes einer veröffentlichten Kopie.
- `inputs/`: gelieferter Plan, ursprüngliche η-Daten, früherer Fitbericht und Herkunfts-Hashes.
- `data/`: JSON-Rohdaten, Kontrolltabellen, exakte skalare Reihe.
- `runtime/`: eingefrorene GP-Engine, eigene Erweiterungen und fünf eingefrorene Wrappermodule. Gemeinsame Caches sind ausgeschaltet.

## Lokale Reproduktion

Python benötigt mpmath, numpy und matplotlib; die lokalen Programme wurden mit Python 3.13.9 unter Windows ausgeführt. Für den reinen mathematischen Modellteil ist PARI/GP nicht nötig. In PowerShell aus diesem Ordner:

```powershell
python quadratic.py --dps 90
python quadratic.py --dps 150 --sizes 4,8,16,32,48,64
python quadratic.py --dps 250 --sizes 64 --toy-only
python quadratic.py --dps 180 --sizes 80 --toy-only
python exact_coefficient.py
```

Für den Vergleichsteil braucht es PARI/GP. Der eingefrorene Wrapper sucht die installierte GP-Executable; alternativ `FATOU_GP_EXE` auf deren vorhandenen Pfad setzen. Er greift nicht auf einen veränderbaren Enginequelltext zurück, solange `runtime/engine_snapshot.gp` vorhanden ist.

```powershell
python baseline_fold.py --theta .25 --dps 25
python baseline_fold.py --theta .5 --dps 25
python baseline_fold.py --theta .5 --dps 45 --samples 80 --radius .025
python baseline_fold.py --theta .75 --dps 25
python baseline_recursion.py
python analyze.py
python verify.py
```

Ein erneuter Rechenlauf überschreibt die dazugehörige Rohdatendatei. Für einen Vergleich mit den veröffentlichten Rohdaten zuerst eine eigene Kopie des Ordners anlegen. Die Manifestprüfung bezieht sich auf die ursprüngliche Veröffentlichung und meldet nach Neuberechnung erwartbar geänderte Dateien, etwa Zeitmessungen.

Rohwerte enthalten 65 signifikante Stellen, um die Präzisionskontrollen zu dokumentieren; daraus folgen keine zertifizierten Fehlerintervalle. Die gewollten 25/45 GP-Digits werden intern auf die Engine-Präzision gerundet; die tatsächlich gemeldete Präzision steht in den Rohdatensätzen. Die beiden Vergleichsrechnungen ändern Präzision und Cauchyquadratur gemeinsam und sind daher keine reine isolierte Präzisionsstudie.

Die Modellskripte begrenzen einen Durchgang auf ungefähr 18 Minuten zwischen den Modellgrössen; die geprüften Einzelmodelle liegen weit darunter. GP-Foldläufe besitzen zusätzlich ein Gesamtbudget von zehn Minuten und einzelne Auswertungszeitlimits. `baseline_recursion.py` benutzt lediglich zwei kleine Höhensamples mit begrenzten GP-Auswertungen. Alle Arbeiter werden in `finally` geschlossen.

## Ablage

Eigener Forschungsordner:

    D:\Users\janis\Documents\Matura\tetration-graph\forschung\eta_hyperoperationen\halb_rang_quadratisch_2026_10_12

Identische lesbare Ausgabekopie im Chat:

    C:\Users\janis\Documents\Codex\2026-10-11\le\outputs\hyper_halbrang

Es wurden keine Dateien des pausierten H₀-Projekts und keine fremden Agentenordner bearbeitet.
