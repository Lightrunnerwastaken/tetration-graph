# Halbrang H₃.₅(2,x) an Basis 2

Einstieg: [BERICHT.md](BERICHT.md). Herleitungen: [MATHEMATIK.md](MATHEMATIK.md). Zwischenstände: ZWISCHENSTAND_01/02. Quellen: [LITERATUR.md](LITERATUR.md).

Python: `C:/Users/janis/miniconda3/python.exe` (3.13, python-flint 0.9, mpmath, numpy, matplotlib). Aus diesem Ordner:

```powershell
python kneser_werte.py                                   # Kneser 2↑↑x an 12 Punkten (ein GP-Lauf, ~5 s)
python rang1_test.py --J 40 --prec 1000                  # Rang 1, 40 Ordnungen (~2 min)
python rang1_test.py --J 52 --prec 1300                  # Rang 1, 52 Ordnungen (~12 min)
python raenge.py --N 40 --J 6  --prec 2400 --name N40J6
python raenge.py --N 24 --J 10 --prec 2400 --name N24J10
python raenge.py --N 40 --J 12 --prec 3200 --name N40J12  # ~7 min
python raenge.py --N 40 --J 8 --R 7 --prec 3000 --name N40J8R7
python pruefe_ganz.py N24J10                             # Superfunktionsgleichung an ganzen Höhen
python wachstum.py N40J12                                # Rangwachstum der Ordnungen
python auswertung.py N40J8R7                             # Halbrang-Jet, Rang 3/2, Rangrekursion
python summen.py                                         # Summen bei eps=1, Vergleich mit Kneser
python leiter.py                                         # Rangleiter H3..H5, Rekursion bei eps=1
python gevrey.py rang1_J52.json                          # Divergenz, Borel-Padé
python bilder.py
```

Dateien: `exppoly.py` (exakter Rangschritt auf Exponentialpolynomen, Ballarithmetik), `halbschritt.py` (Newtonreihe im Rang, Levin), `borel.py`. Rohdaten in `data/`, Abbildung in `figures/`.
