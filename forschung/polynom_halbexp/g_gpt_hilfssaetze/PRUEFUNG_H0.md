# Prüfungen der zweiten Runde

Zeitpunkt: 2026-10-11T11:03:16+02:00.

- Exakte Symbolprüfung der Viererrekursion und der Rang-6-Korrektur erfolgreich.
- 512-Bit-Intervallprüfung des großen kubischen Kreises erfolgreich: Rouché-Marge > 95.15, exaktes L aus Lambert-W.
- Kubisch: N = 200, 300, 400, 600, 800, 1000, 720…2000 Bit, SVD-freie numerische Abschnittslösungen.
- Kneser: N = 200, 300, 400, 600 mit zertifizierten Mittelpunkten; dieselben Kreis- und festen γ*-Punkte geprüft.
- Größte relative Abweichung der KB-Werte vom vorhandenen O38-Raster < 2e-11.
- Gemeinsames Cauchy-Integral: 2400 Bit, 2048/4096 Punkte, N = 200, 600, 1000; Residuenform und Gitterverdoppelung verglichen.
- Das Prüfbild wurde aus den gespeicherten JSON-Daten erzeugt und visuell kontrolliert.
- Syntaxprüfung aller 9 Python-Dateien des Pakets erfolgreich.
- Keine K2-Läufe; keine Änderungen an fremden Forschungsprogrammen.

Keine Intervallzertifizierung der Matrixinversen und kein Beweis von (H₀) für alle N.
Die unabhängige mathematische P-Prüfung steht aus. Der genaue Status der einzelnen Aussagen steht im Bericht.

Bibliotheksversionen: {"python-flint": "0.9.0", "numpy": "2.2.6", "sympy": "1.14.0", "matplotlib": "3.10.9", "mpmath": "1.3.0"}.
