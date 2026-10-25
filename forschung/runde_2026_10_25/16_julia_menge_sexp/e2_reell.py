"""Pruefung: Im F(z) > 0 auf dem oberen Grundstreifen (0<=Re z<1, 0<Im z<=7.5); erste Stelle, wo |Im F(z+1)| = pi erreicht wird
(dann ist F(z+2) reell negativ): Lage der ersten nichtreellen Urbilder der reellen Achse."""
import numpy as np
from kneser_feld import feld
x=(np.arange(400)+0.5)/400; y=np.concatenate([np.geomspace(1e-4,0.05,200),np.linspace(0.05,7.5,600)])
Z=x[None,:]+1j*y[:,None]; F=feld(Z)
print("min Im F auf oberem Streifen:", F.imag.min(), " min Im F/Im z:", (F.imag/Z.imag).min())
F1=np.exp(F)
print("max |Im F(z+1)|:", np.abs(F1.imag).max())
m=np.abs(F1.imag)>=np.pi
i,j=np.nonzero(m)
k=np.argmin(x[j]); print("kleinstes Re z mit Im F(z+1)>=pi: z+1 =", 1+x[j][k], "+i", y[i][k])
k=np.argmax(y[i]); print("groesstes Im z mit Im F(z+1)>=pi: z+1 =", 1+x[j][k], "+i", y[i][k])
