# Test: ist dPhi(eta) an H_4 (b=sqrt2) eine Kombination aus eta(p4), eta'(p4), eta''(p4)?  Daten aus exp.py grad
import numpy as np
p4 = 1.54921732984299
g = {('s',1):-6.12515,('s',2):8.06010,('s',3):-2.8438,('c',1):-8.04229,('c',2):-11.25798,('c',3):9.910}
def eta(a,n,x,d):
    w=2*np.pi*n
    if a=='s': return [np.sin(w*x), w*np.cos(w*x), -w*w*np.sin(w*x)][d]
    return [1-np.cos(w*x), w*np.sin(w*x), w*w*np.cos(w*x)][d]
keys=list(g); y=np.array([g[k] for k in keys])
for D in (1,2,3):
    A=np.array([[eta(a,n,p4,d) for d in range(D)] for a,n in keys])
    c,*_=np.linalg.lstsq(A,y,rcond=None)
    print('Ableitungen bis',D-1,'Koeff',np.round(c,4),'Rest',np.round(A@c-y,3))
