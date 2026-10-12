"""Quadratic rank flow for a declared finite interpolation model.

This is a jet in perturbation amplitude, NOT a full hyperoperation.
Nodes are 0,...,K; d(0)=0. Interpolation uses the degree-K polynomial.
All arithmetic uses mpmath. No upstream inputs or caches are modified.
"""
import argparse, csv, json, math, time
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).absolute().parent

def prefix(c, exclusive=False):
    out=[]; t=mp.mpf(0)
    for z in c:
        if exclusive: out.append(t)
        t += z
        if not exclusive: out.append(t)
    return out

class RankJet:
    def __init__(self,K):
        self.K=K
        weights=[mp.mpf((-1)**j*math.comb(K,j)) for j in range(K+1)]
        self.D=[]
        for i in range(1,K+1):
            row=[weights[j]/weights[i]/(i-j) if j!=i else mp.mpf(0) for j in range(K+1)]
            row[i]=-sum(row)
            self.D.append(row[1:])

    def linear(self,c,s):
        s=mp.mpf(s); kernel=[mp.mpf(1)]
        for j in range(1,self.K): kernel.append(kernel[-1]*(s+j-1)/j)
        return [mp.fsum(kernel[i-j]*c[j] for j in range(i+1)) for i in range(self.K)]

    def B(self,q):
        # B(q)=L diag(W q D^T), W=L-I. W is essential: w(j-1).
        acc=[mp.mpf(0)]*self.K; t=[]
        for i in range(self.K):
            t.append(mp.fsum(acc[j]*self.D[i][j] for j in range(self.K)))
            acc=[acc[j]+q[i][j] for j in range(self.K)]
        return prefix(t)

    def quadratic(self,c,s):
        """Upper block of M**s applied to (0,c tensor c).

        M=[[L,B],[0,L tensor L]]. Its nilpotence index is <=3K-1.
        Prefix sums apply M-I without allocating a K^3 tensor.
        """
        s=mp.mpf(s); K=self.K
        q=[[c[i]*c[j] for j in range(K)] for i in range(K)]
        if s==0: return [mp.mpf(0)]*K
        if s==1: return self.B(q)
        y=[mp.mpf(0)]*K; out=[mp.mpf(0)]*K; choose=mp.mpf(1)
        for p in range(1,3*K-1):
            b=self.B(q)
            wy=prefix(y,exclusive=True)
            y=[wy[i]+b[i] for i in range(K)]
            rowprefix=[prefix(row) for row in q]
            qnext=[[mp.mpf(0)]*K for _ in range(K)]
            for j in range(K):
                t=mp.mpf(0)
                for i in range(K):
                    t+=rowprefix[i][j]
                    qnext[i][j]=t-q[i][j]
            q=qnext
            choose*= (s-p+1)/p
            for i in range(K): out[i]+=choose*y[i]
        return out

    def step(self,c,s):
        a=self.linear(c,s); b=self.quadratic(c,s)
        return [a[i]+b[i] for i in range(self.K)]

    def interpolate(self,c,x):
        # Barycentric degree-K interpolation with prescribed value d(0)=0.
        for j in range(self.K+1):
            if x==j: return mp.mpf(0) if j==0 else c[j-1]
        den=mp.mpf(0); num=mp.mpf(0)
        for j in range(self.K+1):
            term=mp.mpf((-1)**j*math.comb(self.K,j))/(x-j)
            den+=term
            if j: num+=term*c[j-1]
        return num/den

    def exact_full(self,c):
        x=mp.mpf(1); out=[]
        for k in range(1,self.K+1):
            x=x+1+self.interpolate(c,x)
            out.append(x-k-1)
        return out

def maxabs(v): return max(map(abs,v))
def fmt(x): return mp.nstr(x,65)

def run(K,dps,toy_only=False):
    mp.mp.dps=dps; start=time.monotonic(); model=RankJet(K)
    v=[mp.mpf(k*(k-1)) for k in range(1,K+1)]
    a=model.linear(v,mp.mpf('.5')); q=model.quadratic(v,mp.mpf('.5'))
    qfull=model.quadratic(v,1)
    rhs=[u+w for u,w in zip(model.linear(q,'.5'),model.quadratic(a,'.5'))]
    z=[mp.mpf(k*(k-1)*(k-2)*(k+1)*(8*k+1))/60 for k in range(1,K+1)]
    result={'K':K,'dps':dps,'half_linear':[fmt(x) for x in a],
            'half_quadratic':[fmt(x) for x in q],
            'full_quadratic':[fmt(x) for x in qfull],
            'cohomology_abs_error':fmt(maxabs([u-w for u,w in zip(qfull,rhs)])),
            'cohomology_relative_error':fmt(maxabs([u-w for u,w in zip(qfull,rhs)])/(1+maxabs(qfull))),
            'polynomial_full_error':fmt(maxabs([u-w for u,w in zip(qfull,z)])),
            'epsilon_checks':[],'multimode_checks':[]}
    result['noninteger_values']=[]
    for x in map(mp.mpf,['.5','1.5','2.5']):
        exact=2*mp.rf(x,mp.mpf('2.5'))/mp.gamma(mp.mpf('3.5'))-2*mp.rf(x,mp.mpf('1.5'))/mp.gamma(mp.mpf('2.5'))
        result['noninteger_values'].append({'x':fmt(x),'linear_interpolation':fmt(model.interpolate(a,x)),
                                           'linear_gamma':fmt(exact),'linear_error':fmt(model.interpolate(a,x)-exact),
                                           'quadratic_interpolation':fmt(model.interpolate(q,x))})
    # Different non-polynomial node profiles test the whole bilinear map.
    inputs=[] if toy_only else [('gamma_1.3',[mp.rf(k,mp.mpf('1.3'))/mp.gamma(mp.mpf('2.3')) for k in range(1,K+1)]),
                    ('mixed',[mp.sin(k)/7+mp.mpf(k)/11 for k in range(1,K+1)])]
    for label,c in inputs:
        h=model.linear(c,'.5'); b=model.quadratic(c,'.5'); f=model.quadratic(c,1)
        rhs=[u+w for u,w in zip(model.linear(b,'.5'),model.quadratic(h,'.5'))]
        result['multimode_checks'].append({'input':label,'absolute_error':fmt(maxabs([u-w for u,w in zip(f,rhs)]))})
    # Full composition includes cubic and quartic terms omitted in the jet.
    # Only small K is used, because high-degree extrapolation is not controlled.
    if K in (4,6,8):
        for eps in map(mp.mpf,['0.0001','0.00005','0.000025','0.0000125']):
            c=[eps*x for x in v]; halfhalf=model.step(model.step(c,'.5'),'.5')
            full=model.exact_full(c); err=maxabs([u-w for u,w in zip(halfhalf,full)])
            result['epsilon_checks'].append({'epsilon':fmt(eps),'error':fmt(err),'error_over_epsilon_cubed':fmt(err/eps**3)})
    result['seconds']=time.monotonic()-start
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dps',type=int,default=90)
    p.add_argument('--sizes',default='3,4,6,8,10,12,16,20,24,32')
    p.add_argument('--toy-only',action='store_true')
    a=p.parse_args(); results=[]; start=time.monotonic()
    for K in map(int,a.sizes.split(',')):
        if time.monotonic()-start>1080: raise TimeoutError('18-minute total budget reached')
        r=run(K,a.dps,a.toy_only);results.append(r)
        (HERE/f'data/quadratic_{a.dps}{"_toy" if a.toy_only else ""}.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
        print(K,'Q(2)=',r['half_quadratic'][1], 'relerr=',r['cohomology_relative_error'],'seconds=',round(r['seconds'],2),flush=True)
