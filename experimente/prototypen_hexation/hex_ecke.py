"""Kritischer Punkt ζ_c von F_q (Ecke der Kurve) genau; Vergleich w_c = F(ζ_c) mit P(L_T)."""
import sys, time
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
import mpmath as mp
from hexq import HexQ
H = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
PL = H.P([H.L])[0]
print("L_T =", mp.nstr(H.L, 25), "  P(L_T) =", mp.nstr(PL, 25), flush=True)
z = mp.mpc("0.403068905978191", "-1.70987834222861")
# F' und F'' aus einem kleinen Kreis (Cauchy), robuster als Differenzen
for it in range(8):
    r, N = mp.mpf("1e-3"), 16
    pts = [z + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
    try:
        f = H.F(pts)
    except Exception as e:
        print("Fehler", str(e)[:60]); break
    c = [mp.fsum(f[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N / r ** m for m in range(5)]
    dz = -c[1] / (2 * c[2])
    z += dz
    print(f"{it}: ζ_c = {mp.nstr(z, 22)}  F(ζ) ≈ {mp.nstr(c[0], 22)}  |F'| = {mp.nstr(abs(c[1]), 3)}  F''/2 = {mp.nstr(c[2], 8)}  |Schritt| {mp.nstr(abs(dz), 3)}", flush=True)
    if abs(dz) < mp.mpf(10) ** -24:
        break
wc = c[0] - c[1] ** 2 / (4 * c[2])
print("w_c =", mp.nstr(wc, 22), "  w_c − P(L_T) =", mp.nstr(wc - PL, 3))
