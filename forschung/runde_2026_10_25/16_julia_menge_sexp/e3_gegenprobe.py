"""Gegenprobe zu Befund 4 mit anderer Schwelle (R=20) und anderem Gitter (150x900), plus Kreistest der Bilder senkrechter Linien."""
import numpy as np, e1_eintrittsindex as e1
from kneser_feld import feld
e1.RG = 20.0
x=(np.arange(150)+0.37)/150; y=np.linspace(0.01,7.0,900); Z=x[None,:]+1j*y[:,None]
N,_=e1.indizes(feld(Z)); om=2j*np.pi/e1.L-4
zs=Z+om; k=np.floor(zs.real); N2,_=e1.indizes(feld(zs-k)); N2=N2-k.astype(int)
g=(N<e1.NMAX-10)&(N2<e1.NMAX)
for lo in [0,0.5,1,1.5,2,3,4,5]:
    m=(y>=lo)&(y<lo+0.5); print("R=20 Band",lo,(N2==N+4)[m][g[m]].mean())
yy=np.linspace(-7,7,4001)
for x0 in [0,0.25,0.5,0.75]:
    c=feld(x0+1j*yy); A=np.c_[2*c.real,2*c.imag,np.ones_like(yy)]; b=c.real**2+c.imag**2
    s=np.linalg.lstsq(A,b,rcond=None)[0]; m0=s[0]+1j*s[1]; r=np.sqrt(s[2]+abs(m0)**2)
    d=np.abs(np.abs(c-m0)-r); print(f"x0={x0}: Kreis Mitte {m0:.4f} Radius {r:.4f}, max Abw {d.max():.4f} ({d.max()/r:.2%})")
