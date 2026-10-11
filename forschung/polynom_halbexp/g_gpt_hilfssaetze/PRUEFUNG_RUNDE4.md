# Prüfungen der Runde 4

Stand 2026-10-11T17:19:33+02:00.

- Lückenlose Gradranges und JSON-Syntax geprüft: 951 kubische, 551 Kneser-Grade.
- Schur-Rekursion mit C-Matrixfenstern, O(n²) pro Grad; Ausgangsmatrix einmal je Sitzung.
- Fortschreitende Präzision 1200+n Bit; Input bei 1300+max_n Bit.
- Direkte Python-FLINT-Lösungen mit 1500+n Bit, unabhängig vom nativen Bridge und alten Schur-Updates.
- Exakte rationale Paar-/Quellenrekursion für n=0…8; exakte Spaltenidentität für j=0…6.
- Native Bridge gegen rationale Inversen bei 384 Bit geprüft, Fehler <1e-100.
- 21 feste Punkte einschließlich L und 72 Kreiswinkel pro Grad; volle finite Taylormajorante jedes Kreises.
- Kein vollständiges γ*-Kurvenzertifikat und kein All-n-Beweis. Kneser-Mittelpunktradien nicht propagiert.
- Vorherige Berichte, Rohdaten und Einzelgrad-Intervallzertifikate unverändert archiviert.

Laufzeiten jeder endgültigen Scan-Sitzung: {"kubisch": [931.2176702022552, 335.6048803329468], "kneser_z": [428.48647594451904]} Sekunden.
Laufzeiten der direkten Prüfsitzungen: {"kubisch": [630.6668745, 701.6421437263489], "kneser_z": [117.05535292625427]} Sekunden.
Die zwei beendeten Entwicklungsproben waren kürzer als fünf Minuten; ein Vergleichszwischenstand ist archiviert.
Keine Prozesse werden nach Abschluss dieses Durchgangs zurückgelassen. Unabhängige P-Prüfung bleibt offen.
