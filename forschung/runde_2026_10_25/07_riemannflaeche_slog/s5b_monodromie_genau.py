"""s5b: wie s5, aber links/rechts nur 1e-15 neben der Senkrechten Re w = Re L (Motor dps 30 und 40)."""
import sys
import mpmath as mp
import motor
import s5_monodromie as s5

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
L = s5.L
d = mp.mpf("1e-15")
ys = [mp.mpf("1.8"), mp.mpf(2), mp.mpf("2.5"), mp.mpf(3)]
aus = []
for y in ys:
    aus += [f"slog({motor.gpzahl(mp.mpc(L.real - d, y), 40)})", f"slog({motor.gpzahl(mp.mpc(L.real + d, y), 40)})"]
res = motor.werte("exp(1)", aus, dps)
for i, y in enumerate(ys):
    sl, sr = res[2 * i], res[2 * i + 1]
    u = s5.u_von_s(sl)
    pred = (u + s5.g) + s5.th(u + s5.g)
    print(f"dps={dps} y={mp.nstr(y,3)}: Motor rechts {mp.nstr(sr,16)}  Vorhersage M_L(links) {mp.nstr(pred,16)}  Abw {mp.nstr(abs(pred-sr),3)}  Im u={mp.nstr(u.imag,4)}")
motor.schliessen()
