# Analytischer quadratischer Halbrangjet

Einstieg: BERICHT.md. Vollständige Beweisfassung: MATHEMATIK.md. Quellenstatus: LITERATUR.md. Die beiden ZWISCHENSTAND-Dateien dokumentieren die Etappen.

Python 3.13.9, mpmath und für die Grafiken numpy/matplotlib. Es ist kein GP-Worker erforderlich. Die Dateien in inputs sind eingefrorener Kontext aus der vorigen Runde, keine neue Beweisfassung.

Reproduktion aus diesem Ordner:

```powershell
python intrinsic.py --K 128 --dps 300
python newton_height.py data/intrinsic_K128_300.json --dps 400
python newton_height.py data/intrinsic_K128_300.json --dps 500
python intrinsic.py --K 32 --dps 220
python rank_newton.py --N 128 --dps 400
python rank_newton.py --N 32 --dps 220
python direct_coefficients.py --K 256 --dps 400
python direct_coefficients.py --K 512 --dps 500
python evaluate.py data/direct_coefficients_K512_500.json --dps 500
python evaluate.py data/direct_coefficients_K512_500.json --dps 600
python abel_plana.py --dps 50
python abel_plana.py --dps 75
python generating.py --dps 40
python generating.py --dps 65
python verify.py
python analyze.py
```

Die längste gemessene Rechnung dauert rund 4.5 Minuten. Die Schleifen haben eine 18-Minuten-Budgetprüfung. Größere Parameter sind nicht durch diese Messung freigegeben. Partialsummengenauigkeit und Abschneiderest werden ausdrücklich getrennt.

`manifest.json` verzeichnet SHA256 und Dateigrößen des publizierten Pakets. Die Veröffentlichung kopiert ausschließlich diesen eigenen neuen Ordner; der pausierte H₀-Strang wird nicht geändert.
