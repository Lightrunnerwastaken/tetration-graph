import numpy as np, time, sys
from zoo import *
a = abb_exp(1.0)
print('L', a.L, 'lam', a.lam)
for (T,P,n) in [(30,1.0,12),(40,1.0,16),(40,0.5,14)]:
    t0=time.time()
    s = Loesung(a, 0.0, 1.0, T=T, P=P, n=n)
    s.loese(maxit=60)
    v = s.S(np.array([0.5, 0.25+0j, 0.5+0.3j]))
    print(T,P,n, len(s.t), 'S(0.5)=%.17f'%v[0].real, v[0].imag, 'err', v[0].real-1.6463542337511945809, time.time()-t0)
