# Runde 5: Dateien und Reproduktion

BERICHT.md fasst die Ergebnisse zusammen. ETAPPE1.md bis ETAPPE4.md sind die vier eigenständig lesbaren Zwischenstände. Ergebnisse sind ausdrücklich bewiesen / bedingt / numerisch / vermutet / offen gekennzeichnet. Nur (H₀); K2 bleibt bei P.

## Dateien

- family_*.json: neun vollständige Gradfolgen 25…400, Koeffizienten und numerische Kandidatentests.
- precision_checks/: zwei +300-Bit-Familienläufe bis 400.
- holdout_600/: drei separate Fortrechnungen bis 600.
- old_asymptotics.json, second_order.json: Pivot-/Normfits der unveränderten Runde-4-Gradfolgen.
- profiles_*.json: komplexe Werte, Schurpaare und ganze verschobene Taylorpolynome.
- profile_analysis.json, spiral_analysis.json, phase_analysis.json: Form-, Phasen-, Amplituden- und Kontrollfenster.
- profile_direct_checks.json: unabhängige direkte +300-Bit-Kontrolle bei 400.
- spiral_geometry_certificate.json: rigoroses **skalares** Geometriezertifikat; kein Inversenzertifikat.
- exact_models.json, structure_checks.json: exakte rationale Identitäts- und Konstantenkontrollen.
- final_controls.json, verification.json: getrennte numerische Kontrollen und Abschlussprüfung.
- experimente/: sämtliche Rechenprogramme und der native FLINT-Wrapper. Dieser rechnet Mittelpunkte, keine zertifizierten Inversen.

## Umgebung

Verwendet wurden Windows, Python 3.13.9, python-flint 0.9.0, NumPy, mpmath, SymPy und matplotlib. Der native Wrapper lädt die DLL aus dem installierten python-flint-Wheel, ohne dessen Pythonobjektlayout auszulesen. Für andere Betriebssysteme müsste der DLL-Lader angepasst werden.

PowerShell, im veröffentlichten Forschungsordner:

    $researchRoot = 'D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp'
    $pyExe = 'C:\Users\janis\miniconda3\python.exe'
    Set-Location -LiteralPath "$researchRoot\g_gpt_hilfssaetze\runde5_szego"

## Schnelle Auswertung der abgelegten Daten

Die folgenden Programme erkennen den Forschungsordner automatisch. Im ursprünglichen Codex-Arbeitsordner erkennen sie stattdessen work/h0_round5 und work/h0_round4.

    & $pyExe -B .\experimente\h0_r5_exact_models.py
    & $pyExe -B .\experimente\h0_r5_structure_checks.py
    & $pyExe -B .\experimente\h0_r5_geometry_certificate.py
    & $pyExe -B .\experimente\h0_r5_old_asymptotics.py
    & $pyExe -B .\experimente\h0_r5_family_analysis.py
    & $pyExe -B .\experimente\h0_r5_second_order.py
    & $pyExe -B .\experimente\h0_r5_profile_analysis.py
    & $pyExe -B .\experimente\h0_r5_spiral.py
    & $pyExe -B .\experimente\h0_r5_phase.py
    & $pyExe -B .\experimente\h0_r5_finish_analysis.py
    & $pyExe -B .\experimente\h0_r5_plots.py

Optional können --data-dir, --round4-dir und --project-root die automatischen Pfade ersetzen. Die Reihenfolge oben berücksichtigt die Datenabhängigkeiten. Die Analyseprogramme schreiben ihre kleinen Ausgabedateien neu.

## Numerische Rohdaten neu rechnen

Reproduktionen in einen eigenen Unterordner schreiben. Die bekannten Laufzeiten gelten für die tatsächlich abgelegten Sitzungen, nicht als Garantie für eine andere Rechnerlast.

    & $pyExe -B .\experimente\h0_r5_families.py --project-root $researchRoot --output-dir .\nachrechnung --max-n 400 --models base,jet_0.01,jet_-0.01,mu_0.01,mu_-0.01,quad_0.2,quad_0.31382,tau_0.1,tau_10
    & $pyExe -B .\experimente\h0_r5_families.py --project-root $researchRoot --output-dir .\nachrechnung\precision_checks --max-n 400 --models jet_0.01,tau_0.1 --extra-bits 300
    & $pyExe -B .\experimente\h0_r5_families.py --project-root $researchRoot --output-dir .\nachrechnung\holdout_600 --max-n 600 --models jet_0.01
    & $pyExe -B .\experimente\h0_r5_families.py --project-root $researchRoot --output-dir .\nachrechnung\holdout_600 --max-n 600 --models jet_-0.01
    & $pyExe -B .\experimente\h0_r5_families.py --project-root $researchRoot --output-dir .\nachrechnung\holdout_600 --max-n 600 --models tau_10

Der Familienläufer speichert jeden abgeschlossenen Modelllauf sofort und beendet die Sitzung vor dem nächsten Modell, wenn die Laufzeitprognose das Budget überschreitet. Nicht gestartete Modelle in einer eigenen Sitzung rechnen.

Komplexe Profile mit getrennten Sitzungen:

    & $pyExe -B .\experimente\h0_r5_profiles.py --project-root $researchRoot --output-dir .\nachrechnung --model kubisch --max-n 900 --degrees 50,100,150,200,300,400,500,600,700,800,850,900 --checkpoint-out .\zwischenmaterial\kubisch900.gz
    & $pyExe -B .\experimente\h0_r5_profiles.py --project-root $researchRoot --output-dir .\nachrechnung --model kubisch --max-n 1000 --degrees 950,951,952,953,954,955,956,957,958,959,960,970,980,990,1000 --checkpoint .\zwischenmaterial\kubisch900.gz
    & $pyExe -B .\experimente\h0_r5_profiles.py --project-root $researchRoot --output-dir .\nachrechnung --model kneser_z --max-n 600 --degrees 50,100,150,200,300,400,500,550,600

Der tatsächlich abgelegte zweite kubische Lauf nutzte den bereits vorhandenen Runde-4-Checkpoint bei 914. Die optionale Export-/Resume-Schnittstelle erlaubt eine eigenständige Reproduktion ab 900; sie wurde zusätzlich bei 8→12 exakt gegen einen frischen Lauf geprüft. Große Checkpoints werden nicht als Ergebnisdateien veröffentlicht. Der Profilrechner hat ein weiches Rechenbudget von 850 Sekunden; gespeicherte Profile bleiben bei Budgetabbruch erhalten. Der optionale Checkpointexport hat eine weitere Budgetkontrolle.

    & $pyExe -B .\experimente\h0_r5_profile_direct.py

Dies prüft die abgelegten komplexen n=400-Werte unabhängig. Für einen separaten Reproduktionsordner entsprechend --data-dir .\nachrechnung übergeben.

## Genauigkeit und mathematische Reichweite

Die Präzision wächst mit n: 1200+n Bit, mit zusätzlichen 300 Bit in den Gegenproben; die Eingabematrix wird vorher bei höherer Präzision aufgebaut. Ergebnisse speichern 45 bis 90 Dezimalstellen je nach Objekt. Die ursprünglichen A-Koeffizienten mit 70 Stellen dürfen in der Bucht nicht nachträglich zur hochpräzisen Auswertung benutzt werden: dafür liegen die bereits vor der Speicherung verschobenen Kreis-Taylorpolynome vor.

Kein endlicher hochpräziser Lauf beweist eine All-n-Asymptotik. Das Geometriezertifikat zertifiziert keine Matrizeninversen. Finite algebraische Checks ersetzen nicht die ausgeschriebenen Beweise. Alle tatsächlichen schweren Sitzungen blieben unter 20 Minuten.
