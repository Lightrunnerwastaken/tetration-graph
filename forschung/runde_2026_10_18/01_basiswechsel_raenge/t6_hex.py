"""Motor-Hexation (Kneser über eta5): Geschwindigkeit, Genauigkeit, Konsistenz mit lib-Pentation."""
import time
import mpmath as mp
import lib, motor
mp.mp.dps = 100
for b in ["2", "exp(1)"]:
    T = lib.Tet(b)
    t = time.time()
    h = motor.werte(b, ["hexx(0.5)", "hexx(-0.5)", "hexx(0.25)", "hexx(1.25)", "hexx(1)", "hexx(0)"], 60)
    print(b, f"{time.time()-t:.1f}s", [mp.nstr(x, 25) for x in h])
    # Gleichung hex(x+1) = P(hex(x))
    print("  hex(1.25) - P(hex(0.25)) =", mp.nstr(h[3] - T.P(mp.re(h[2])), 3))
    t = time.time()
    ih = motor.werte(b, [f"ihexx({motor.gpzahl(mp.re(h[0]), 60)})"], 60)
    print("  ihex(hex(0.5)) - 0.5 =", mp.nstr(ih[0] - mp.mpf(0.5), 3), f"{time.time()-t:.1f}s")
motor.schliessen()
