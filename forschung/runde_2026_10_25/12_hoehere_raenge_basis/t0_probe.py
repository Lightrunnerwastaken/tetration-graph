"""Probe: Laufzeit und Werte von H_4, H_5, H_6 bei einer Basis, auch bei negativen Hoehen.
Aufruf: python t0_probe.py <basis> <dps>"""
import sys, time
import mpmath as mp
import motor

b = sys.argv[1]; dps = int(sys.argv[2])
aus = ["hyper(4,0.5)", "hyper(5,0.5)", "hyper(5,-2.5)", "hyper(5,-6)", "hyper(6,0.5)", "hyper(6,-2.5)", "hyper(6,-3.2)"]
t = time.time()
for a in aus:
    try:
        v = motor.werte(b, [a], dps)[0]
        print(a, mp.nstr(v, 20), f"{time.time()-t:.1f}s", flush=True)
    except Exception as e:
        print(a, "FEHLER", str(e)[:200], flush=True)
motor.schliessen()
