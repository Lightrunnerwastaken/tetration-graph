import sys, mpmath as mp
import motor
b = sys.argv[1] if len(sys.argv) > 1 else "1/2"
try:
    r = motor.werte(b, ["sexp(0)", "sexp(1)", "sexp(0.5)", "sexp(0.25)", "sexp(1.5)", "sexp(0.5+0.5*I)"], 20)
    for x in r: print(mp.nstr(x, 18))
finally:
    motor.schliessen()
