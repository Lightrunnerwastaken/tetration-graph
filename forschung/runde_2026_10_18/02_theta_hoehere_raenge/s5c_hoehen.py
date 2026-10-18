"""Hilfsrechnung: ipent(L3) (Punkt, an dem die regulaere Pentation L3 annimmt) und Rang-7-Konstanten (p6, lam6). ihyper(7,kp) ist entartet (kp ist Grenzwert von hex), nicht verwenden."""
import time, motor, mpmath as mp
t=time.time()
L3 = "(0.82467854614207440+1.5674321238496478*I)"
r = motor.roh("2", ["pentx(0.5)", f"ipentx({L3})", f"ipentx(conj({L3}))"], 30)
print("ipent(L3) =", r[1], " ipent(conj L3) =", r[2], time.time()-t, flush=True)
r = motor.roh("2", ["hexx(0.5)", "kp", "hyperdigits(7)", "hyp[7]", "hylam[7]"], 30)
print("kp", r[1], "Rang7 Stellen", r[2], "p6", r[3], "lam6", r[4], time.time()-t, flush=True)
kp = r[1]
r = motor.roh("2", [f"iferr(ihyper(7, {motor.gpzahl(kp,30)}),E,1000)", f"iferr(ihyper(7, conj({motor.gpzahl(kp,30)})),E,1000)", "hyper(7,0.5)"], 30)
print("ihyper(7,kp) =", r[0], r[1], " H7(0.5)=", r[2], time.time()-t, flush=True)
