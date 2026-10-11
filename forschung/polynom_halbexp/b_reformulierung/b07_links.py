"""b07: Konvergenz auf der linken Seite von D (z = -1.3, -1.35, Winkel) gegen Knesers Taylorreihe (1000 Koeff.)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
a0 = K[0]
def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s
for rr, phi in ((1.3, 180), (1.35, 180), (1.3, 135), (1.3, 120), (1.36, 150), (1.3, 0), (1.35, 0), (1.3, 45), (1.33, 70)):
    z = rr * mp.expjpi(mp.mpf(phi) / 180)
    ref = horner(K, z)
    lem = abs(z) * abs(z - a0) / (abs(L) * abs(L - a0))
    row = [mp.nstr(abs(horner(P[N], z) - ref), 3) for N in (20, 40, 80, 120, 160, 194)]
    print(f"z = {rr}·e^(i{phi}°): Lemniskate (τ=1/2) Wert {mp.nstr(lem,4)}; |p_N - f_K| N=20,40,80,120,160,194:", row)
