# η höherer Hyperoperationen

Eigenständige Forschung ab 12.10.2026. Das Projekt `polynom_halbexp` und seine gemeinsamen Dateien wurden hierfür nicht bearbeitet. Dieser Ordner untersucht eine ausgewählte analytische Hierarchie; er behauptet keinen allgemeinen Kneser-Satz für alle Ränge.

- **BERICHT.md:** Ergebnisse, η-Tabelle, Genauigkeit, Konstruktionen und offene Fragen.
- **MATHEMATIK.md:** Definitionen und vollständig bewiesene bedingte Aussagen, einschließlich Auswahlabhängigkeit.
- **ZWISCHENSTAND_1.md / ZWISCHENSTAND_2.md:** unveränderte Zwischenstände während der Arbeit.
- **analysis.json / eta_folge.png:** numerische Folgenanalyse und Darstellung.
- **eta7_*.json / eta8_*.json:** einzelne Berührsuchen einschließlich gespeicherter Teilläufe.
- **probe_hyper_*.json:** Konstruktionen und auch fehlgeschlagene Versuche; `status` und `error` lesen.
- **real_decation.json / decation_extrapolation.json:** drei reelle Koenigs-Näherungen der Decation und ihr Präzisionsvergleich.
- **archiv/:** historische GP-Codekombinationen und der anfängliche Hilfsprototyp; diese Dateien sind unveränderte Arbeitsstände, keine empfohlenen Programme.

## Rechenwerk und Reproduktion

Python mit `mpmath`, `numpy`, `scipy`, `matplotlib` sowie PARI/GP wird benötigt. Die bereits lokal vorhandene Python-Anbindung ist unter `_backend/fatou_backend` eingefroren. Die GP-Engine ist der unveränderte `engine_snapshot.gp`; ihr SHA-256 steht in `provenance.json`. Die eigenen Erweiterungen stehen in `hyper.gp` und den generierten `even8*.gp`. Die produktiven und gemeinsamen Upstream-Caches bleiben ausgeschaltet. Beim Start entsteht eine **eigene** inhaltlich gehashte GP-Kombinationsdatei im Forschungsordner.

Aus diesem Ordner, beispielsweise mit `C:\Users\janis\miniconda3\python.exe -B` als Python-Aufruf:

```text
python -B probe_lower.py
python -B probe_hyper.py --base 2 --rank 7 --dps 45 --timeout 180
python -B critical.py --rank 7 --b0 1.8397586607 --b1 1.8397586608 --dps 45 --samples 64 --radius .025 --x0 3.5434897855 --steps 4 --output eta7_neu.json

python -B build_even.py --pilot --target 10
python -B probe_hyper.py --base 2 --rank 8 --dps 25 --timeout 600 --tag _neu

python -B build_even.py --pilot --target 22 --output even8_control.gp
python -B probe_hyper.py --base 2 --rank 9 --dps 45 --even even8_control.gp --timeout 600 --tag _neu
python -B real_attractor.py
python -B analyze.py
python -B summarize.py
```

`critical.py --resume DATEI` setzt ausschließlich mit den letzten zwei gespeicherten Basis-/Minimapunkten fort; bereits berechnete Stützpunkte werden nicht wiederholt. `--eval-seconds` begrenzt einen GP-Aufruf, `--max-seconds` eine Sitzung einschließlich Reserve für Initialisierung. Die gewählten Läufe bleiben unter ungefähr 20 Minuten. Ein Zeitabbruch ist ein Rechenbefund, kein Nichtexistenzresultat.

Die schnelle Bogenplanung mit `--pilot` ist eine explizite Änderung der von Hexation übernommenen Heuristik: wenige Mittelpunkte, bis 10 Sekunden pro Kandidat statt einer Sekunde, bis 10 Sekunden gemessene Planzeit statt 400 Millisekunden. Das Ziel `--target` ist eine Planungsgröße. Die erreichten Stellen werden separat aus zwei Bögen geschätzt.

`hinit` schaltet bei einer explorativen positiven Fixpunktsuche ab Höhe 4 auf die vorhandene komplexe Auswahl um. Im tatsächlich untersuchten Parameterbereich liegen die gemessenen anziehenden positiven Fixpunkte darunter. Diese Abbruchregel ist **kein globaler Beweis**, dass oberhalb keine weiteren Fixpunkte existieren. Für komplexe geradzahlige Stufen unterstützt die aktuelle Erweiterung nur H₈. Reguläre reelle Stufen können darüber hinaus untersucht werden.

## Grenzen der Kontrollen

Normalisierung und Funktionalgleichung werden teilweise durch die Fortsetzungskette selbst erzwungen. Kleine Restfehler dafür sind Konsistenzprüfungen, keine unabhängigen Genauigkeitsbeweise. Maßgeblich sind veränderte Präzision, Cauchy-Radius und Punktzahl, zwei Bogenlösungen sowie zusätzliche lokale Fortsetzungsschritte. Bei der komplexen H₈ hat das Argument `extra` gegenwärtig keine Wirkung; sein Nullunterschied wird **nicht** als unabhängige Kontrolle gezählt.

Alle höheren η-Werte sind numerische lokale Berührungen der gewählten Fortsetzungen. Für ihre Deutung als globale einzige Beschränktheitsschwellen gelten die in MATHEMATIK.md genannten Voraussetzungen. Die Vermutungen ηᵣ→2 und eᵣ→4 bleiben offen.

## In den Wissensgraphen übernommen (13.10.2026)

Datei `daten/aussagen/kn_hoehere_raenge.yaml` (Theorie: `daten/theorie/kn_hoehere_raenge.yaml`), Stränge `kn.etas`, `kn.hoehere_stufen` und neu `kn.halbrang`:
`kn.beob.eta7_eta8`, `kn.satz.schwellen_monoton` (M1, M2), `kn.satz.basis2_hoehe3` (M2a), `kn.satz.kritische_rate` (M3), `kn.satz.eta_auswahl` (M4),
`kn.satz.eta_kriterium` und `kn.beob.kritische_profile` (rang_gamma), `kn.beob.zwischenwerte_basis2`, `kn.beob.decation_reell`,
`kn.satz.halbrang_linear` (rang_gamma §5), `kn.satz.halbrang_quadratisch` (halb_rang_analytisch §6), `kn.gescheitert.polynom_rekonstruktion`,
`kn.gescheitert.mischkurve` (halb_rang_quadratisch), `kn.frage.halber_rang`.
Nicht übernommen: Fit-Formeln als eigene Aussage, die Gamma-Signaturen am Berührpunkt (keine unabhängige Evidenz, rang_gamma §4), die Vermutung
zur Randstruktur der Newtonkoeffizienten (halb_rang_analytisch §8, Lücke benannt), die experimentelle Zwei-Bogen-Auswahl als Satz.
