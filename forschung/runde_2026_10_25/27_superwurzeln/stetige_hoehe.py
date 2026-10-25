# Superwurzel zu gebrochener Hoehe (Kneser): gesucht b > eta mit sexp_b(s) = Y.
# Sekantenverfahren in b; jede Basis ist eine eigene Kneser-Sitzung (motor.py, eigener Zwischenspeicher).
# Aufruf: python stetige_hoehe.py Y s b0 b1 [dps]
import sys, time
import mpmath as mp
import motor

Y = mp.mpf(sys.argv[1]); s = sys.argv[2]
b0 = mp.mpf(sys.argv[3]); b1 = mp.mpf(sys.argv[4])
dps = int(sys.argv[5]) if len(sys.argv) > 5 else 20
mp.mp.dps = dps + 5

def f(b):
    t = time.time()
    bs = mp.nstr(b, dps + 3)
    w = motor.werte(f"({bs})", [f"sexp({s})"], dps)[0]
    print(f"  b={mp.nstr(b,dps)}  sexp_b({s})={mp.nstr(w,dps)}  ({time.time()-t:.1f}s)", flush=True)
    return mp.re(w) - Y

try:
    f0, f1 = f(b0), f(b1)
    for it in range(8):
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
        b0, f0, b1 = b1, f1, b2
        f1 = f(b1)
        if abs(f1) < mp.mpf(10) ** (-dps + 3):
            break
    print(f"ERGEBNIS Y={Y} s={s}: b = {mp.nstr(b1, dps-3)}  Rest {mp.nstr(f1,3)}")
finally:
    motor.schliessen()
