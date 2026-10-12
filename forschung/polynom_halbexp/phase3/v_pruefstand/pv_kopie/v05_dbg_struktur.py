"""V-Diagnose: warum liefert p32 im Strukturgebiet nan / riesige sup|∂⁷f|? Eine S-Zelle bei θ ≈ 0.5 einzeln auswerten.
Aufruf: python v05_dbg_struktur.py [theta] [h]"""
import sys
import time

from flint import acb, arb, ctx

ARGS = sys.argv[1:]
sys.argv = [sys.argv[0], "test"]
import p32_rand_ableitungen as P32

sys.stdout.reconfigure(encoding="utf-8")
th = float(ARGS[0]) if len(ARGS) > 0 else 0.5
h = float(ARGS[1]) if len(ARGS) > 1 else 0.004
t0 = time.time()
RD = P32.Rand()
K = RD.K
rho = RD.R.absL
I = acb(0, 1)
w = acb("0.31")
for tt in [0.05, 0.15, 0.25, 0.35, 0.45, th]:
    w = K.slog_num(rho * (I * arb(tt)).exp(), w)
print("slog-Start", w.str(10), f"({time.time()-t0:.0f}s)")
hb = h / 2
xm = rho * (I * arb(th)).exp()
xrad = rho * arb(hb) * arb("1.0001")
zc, zr = K.slog_krawczyk(xm, w, xrad)
print("Zellscheibe z:", zc.str(10), "Radius", zr.str(3))
zm, zmr = K.slog_krawczyk(xm, zc, arb(0))
print("Mitte z:", zm.str(10), "Radius", zmr.str(3))
for name, Z in [("Punkt", P32.disc(zm, zmr)), ("Zelle", P32.disc(zc, zr))]:
    t1 = time.time()
    try:
        tpiI = 2 * acb.pi() * acb(0, 1)
        q0 = (tpiI * Z).exp()
        th0 = RD.theta_jets(q0, 0)[0]
        w0 = (RD.L * (Z + th0)).exp()
        print(f"  [{name}] q0 = {q0.str(5)}, |w0| ≈ {abs(w0).str(5)}")
        Gw, R0 = RD.G_werte(q0, 7)
        print(f"  [{name}] G_k:", [g.str(3) for g in Gw])
        fak = [1]
        for k in range(1, 8):
            fak.append(fak[-1] * k)
        Hk = [RD.sqrtL * w0 * R0] + [RD.sqrtL * w0 ** (1 - k) * Gw[k] / fak[k] for k in range(1, 8)]
        print(f"  [{name}] H_k:", [x.str(3) for x in Hk])
        Pw = RD.P_jet(w0, 7)
        print(f"  [{name}] P_jet(w0):", [x.str(3) for x in Pw])
        PH = RD.P_jet(Hk[0], 7)
        print(f"  [{name}] P_jet(H0):", [x.str(3) for x in PH])
        aus, _, _ = RD._aus_H(w0, Hk, 7)
        print(f"  [{name}] ∂_θ^j f:", [x.str(4) for x in aus], f"({time.time()-t1:.0f}s)")
    except Exception as e:
        print(f"  [{name}] FEHLER {type(e).__name__}: {e}")
