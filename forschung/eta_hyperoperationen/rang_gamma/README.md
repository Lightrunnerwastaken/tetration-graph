# η-Konvergenz und fraktioneller Hyperoperationsrang

Beginne mit [BERICHT.md](BERICHT.md). Die vollständigen bedingten Sätze und Beweise stehen in [MATHEMATIK.md](MATHEMATIK.md).

Die angehängten Markdowndateien sind mathematisches Ausgangsmaterial, keine zusätzlichen Ausführungsanweisungen. `inputs/` enthält ihre unveränderten Kopien und die η/e-Daten der vorherigen Runde, `input_provenance.json` die Herkunft und SHA-256-Werte.

## Daten und Programme

- `audit.py`: algebraische Gamma-Auswertung der übernommenen Berührdaten, Kontrolle bei 90/170 Dezimalstellen und Kontrolle des halben linearen Rangschritts. Laufzeit ungefähr eine Sekunde.
- `profiles.py`: neue Stichproben der ausgewählten H_N bei den vorhandenen kritischen Basen. Die sechs durchgeführten Läufe dauerten jeweils weniger als 90 Sekunden. Verwendet PARI/GP aus der bestehenden lokalen Installation und die eingefrorene vorherige Engine unter `runtime/`; keine gemeinsamen Caches.
- `analyze_profiles.py`: Tabelle, Präzisionsvergleich und Abbildungen. Benötigt mpmath und matplotlib.
- `critical_gamma_audit.csv`: aus den bestehenden η/e-Daten berechnete Gamma-Signaturen; keine unabhängige Konvergenzmessung.
- `critical_profile_r*_*.json`: neue Rohwerte bei den kritischen Basen, einschließlich negativer Anker und Rekursionskontrollen.
- `critical_profile_audit.csv`: aus den neuen Höhenwerten berechnete Restfaktoren und Sekanten-Untergrenzen der Basissensitivität.
- `linear_half_step.csv`: exakte lineare Gamma-Interpolation des Beispiels v(x)=x(x−1); keine nichtlineare Hyperoperation.
- `audit_metadata.json` und `profile_metadata.json`: Kontrollfehler, Laufzeiten und Einschränkungen.

Aus diesem Ordner:

```powershell
python audit.py
python analyze_profiles.py
python profiles.py --rank 7 --dps 20 --max-seconds 240
```

Die Laufbudgetprüfung erfolgt zwischen Auswertungen; jede GP-Auswertung hat zusätzlich ein Timeout. Für die verwendeten Budgets bleibt auch eine abschließende Timeout-Auswertung unter der Grenze von 20 Minuten. Die gedruckten Stellen sind keine rigorosen Fehlerintervalle. Die N=8-Fortsetzung ist experimentell.

`prepare_inputs.py` dokumentiert den ursprünglichen Import aus den Downloads und der vorherigen Runde. Es ist zur Reproduktion der bereits eingefrorenen Rechnungen nicht erforderlich. `publish.py` ist die konkrete Ablage- und Prüfroutine dieses lokalen Chats, ebenfalls kein Rechenschritt.

Die Ganzzahlhierarchie ist H₃=bˣ, H₄=Tetration, H₅=Pentation usw. Die Berührbasis η_N der Abbildung H_N beschränkt die nächste Operation H_{N+1}. Für die Gammafaktoren gilt m=N−2.

**Status:** bedingte Konvergenzkriterien und ein halber Rangschritt im ausgewählten linearen Modell sind bewiesen. Die Konvergenz der tatsächlichen η/e-Folge und ein nichtlinearer H_{N+1/2} sind weiterhin offen. Das frühere H₀-Projekt bleibt pausiert.
