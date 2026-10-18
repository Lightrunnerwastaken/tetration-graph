"""Gegenprobe der Nullstellen von S_K in H (Einstellung B) und Wert von 1+theta' dort (z0+1)."""
import numpy as np
from zoo import *
from s1_loesen import ABB
for nm, z0 in [('q1', 2.179932220454+2.597286162123j), ('q4', 1.043592393040+1.349856055413j)]:
    a = ABB[nm][0]()
    for e in 'AB':
        s = Loesung.laden(a, f'{nm}_{e}'); z = z0
        for _ in range(30):
            v = s.Sweit(np.array([z]))[0]; dv = (s.Sweit(np.array([z+1e-6]))[0]-s.Sweit(np.array([z-1e-6]))[0])/2e-6; z -= v/dv
        print(nm, e, 'Nullstelle S:', f'{z:.13f}', ' |S| =', abs(s.Sweit(np.array([z]))[0]))
