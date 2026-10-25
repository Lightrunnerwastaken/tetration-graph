"""m3: Knesers theta an der ganzen Zahl: z + theta(z) -> unendlich fuer z -> -1 von oben (Satz M5).

theta(s) = sum a_n e^{2 pi i n s} (Koeffizienten aus identitaeten_2026_10_14/etappe2, nur lesen, n <= 4096 verlaesslich).
Wegen F_L(-1 + i y + theta(i y)) = sexp_K(-1 + i y) -> 0 und F_L != 0 muss |theta(i y)| -> unendlich.
Gegenprobe: zwei Reihenlaengen (2048 / 4096 Glieder).
"""
import sys
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "runde_2026_10_18" / "10_neues_gebiet"))
import lib10 as lib

mp.mp.dps = 30
B = lib.Basis(K=40)
for nmax in (2048, 4096):
    T = lib.ThetaReihe(lib.lade_theta(nmax=nmax), B)
    print(f"Reihe bis n = {nmax}")
    for y in ["0.3", "0.1", "0.03", "0.01", "0.006", "0.004"]:
        y = mp.mpf(y)
        th = T.theta(1j * y)
        w = -1 + 1j * y + th
        Fw = B.F(w)
        print(f"  y={mp.nstr(y,3):6s} theta(iy) = {mp.nstr(th, 12):34s} |w| = {mp.nstr(abs(w), 8):10s} |F_L(w)| = {mp.nstr(abs(Fw), 6):10s}"
              f" |F_L(w)|/y = {mp.nstr(abs(Fw)/y, 8)}")
    sys.stdout.flush()
