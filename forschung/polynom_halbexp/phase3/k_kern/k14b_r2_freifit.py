"""k14b — freie Exponenten-Anpassung der Profilkorrekturen S_n = S_∞ + Σ C_i n^{−α_i} (Gs Profile, n ≥ 200) und der Amplitude B_n."""
import numpy as np
from scipy.optimize import minimize
D=np.load(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'k14_profile_G_kubisch.npz')); n=D['n'].astype(float); S=D['S']; B=D['B']
L=0.31813150520476413+1.3372357014306895j; mu=np.sqrt(L)
kappa=2j*np.pi/(2j*np.pi-np.log(mu)); a=1/kappa
print('a =',a, ' 2a =',2*a)
K=12
sel=n>=200
nn=n[sel]; SS=S[sel][:,1:K]   # skip constant coefficient (=1)
# weights: equalize coefficient scales
w=1/np.maximum(np.abs(SS[-1]),1e-12)
def resid(alphas, ret=False):
    cols=[np.ones_like(nn,dtype=complex)]+[nn**(-al) for al in alphas]
    A=np.array(cols).T
    X,res,rk,sv=np.linalg.lstsq(A,SS*w,rcond=None)
    R=A@X-SS*w
    if ret: return X,R
    return np.sum(np.abs(R)**2)
def fit(m, starts):
    best=None
    for st in starts:
        x0=np.concatenate([np.real(st),np.imag(st)])
        f=lambda x: resid(x[:m]+1j*x[m:])
        r=minimize(f,x0,method='Nelder-Mead',options={'maxiter':40000,'xatol':1e-10,'fatol':1e-30})
        if best is None or r.fun<best.fun: best=r
    al=best.x[:m]+1j*best.x[m:]
    return al,best.fun
for m,starts in [(1,[[a],[1.0]]),(2,[[a,1.0],[0.8,1.2],[1.0,2.0]]),(3,[[a,1.0,2*a],[0.8,1.1,1.8],[0.9,1.0,1.5]])]:
    al,f=fit(m,[np.array(s,complex) for s in starts])
    order=np.argsort(al.real)
    print('m=%d: Exponenten %s, Rest %.3e'%(m, ', '.join('%.4f%+.4fi'%(x.real,x.imag) for x in al[order]), f))
# fixed G-exponents for comparison
for alist,name in [([a,1.0,2*a],'G: a,1,2a'),([a,1.0],'a,1'),([1.0,2.0],'1,2'),([a,2*a],'a,2a')]:
    print('%-10s Rest %.3e'%(name,resid(np.array(alist,complex))))
# amplitude B_n: real; fit B_n = B_inf + sum Re(c_i n^{-al_i}) with real exponents
nb=n[n>=200]; Bb=B[n>=200]
def residB(al):
    A=np.array([np.ones_like(nb)]+[nb**(-x) for x in al]).T
    X=np.linalg.lstsq(A,Bb,rcond=None)[0]; return np.sum((A@X-Bb)**2),X
for m,st in [(1,[0.9]),(2,[0.9,1.8]),(2,[0.5,1.0]),(3,[0.9,1.0,1.8])]:
    r=minimize(lambda x: residB(x)[0], st, method='Nelder-Mead',options={'maxiter':20000,'xatol':1e-10,'fatol':1e-30})
    f,X=residB(r.x)
    print('B: m=%d Exponenten %s, B_inf=%.5f, Rest %.2e'%(m, np.round(np.sort(r.x),4), X[0], f))
