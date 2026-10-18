import mpmath as mp, rk_lib as R, sys
R.init(60); N=100
def rel(a,b): return max(abs(x-y)/abs(y) for x,y in zip(a,b) if y!=0)
import time; t=time.time()
a=R.shift(R._Z[mp.mpf(0)],mp.mpf('0.1'),N); b=R.shift(R._Z[mp.mpf('0.2')],mp.mpf('-0.1'),N); print('0.1', mp.nstr(rel(a,b),3))
a=R.exp_s(R.kneser(-0.3,N),N); b=R.shift(R._Z[mp.mpf('0.4')],mp.mpf('0.3'),N); print('0.7', mp.nstr(rel(a,b),3))
a=R.log_s(R.kneser(0.4,N),N); b=R.kneser(-0.6,N) ; c=R.shift(R._Z[mp.mpf('-0.4')],mp.mpf('-0.2'),N); print('-0.6', mp.nstr(rel(a,c),3))
for k in [1,2,3,10,50,99]: print(k, mp.nstr(a[k],12), mp.nstr(c[k],12))
print(time.time()-t)
