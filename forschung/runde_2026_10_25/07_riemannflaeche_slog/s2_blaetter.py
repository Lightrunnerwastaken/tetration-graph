"""s2: Zwei Blaetter von slog im Vergleich.
 h(w) = A(e^w) - 1  (periodisches Halbebenen-Blatt, Re w < Re L, aus der Taylorreihe um 0)
 k(w) = slog des Motors (Knesers Hauptblatt mit den Schnitten des Motors)
Erwartung: wo beide Blaetter verschieden sind, unterscheiden sie sich um die Monodromie bei L:
 Phi(k) - Phi(h) = +-gamma, Phi = id + theta, theta = O(e^{-2 pi Im s}).
Aufruf: python s2_blaetter.py <dps>"""
import json
import sys
import mpmath as mp
import motor

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = 50
D = json.load(open(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_slog_koeff_300.json"))
A = [mp.mpf(a) for a in D["A"]]
L = mp.mpc("0.3181315052047641353126542515876645172035", "1.33723570143068940890116214319371061254")
g = 2j * mp.pi / L


def Aser(z):
    return mp.polyval(A[::-1], z)


pkte = [mp.mpc(-1, 2), mp.mpc(-0.5, 5), mp.mpc(0.2, 1.0), mp.mpc(-2, -3), mp.mpc(0.1, 1.6), mp.mpc(-3, 7.9), mp.mpc(0, 3)]
aus = [f"slog({motor.gpzahl(w, 30)})" for w in pkte]
res = motor.werte("exp(1)", aus, dps)
for w, k in zip(pkte, res):
    h = Aser(mp.exp(w)) - 1
    d = k - h
    print(f"w={mp.nstr(w,4):>14}  Motor k={mp.nstr(k,14):>36}  Halbebene h={mp.nstr(h,14):>36}")
    print(f"      k-h={mp.nstr(d,14)}   (k-h)-gamma={mp.nstr(d-g,6)}  (k-h)+conj(gamma)={mp.nstr(d+mp.conj(g),6)}  Im h={mp.nstr(h.imag,4)}")
motor.schliessen()
