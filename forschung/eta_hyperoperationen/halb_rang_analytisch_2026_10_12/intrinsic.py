"""Exact rational rank jets at integer heights, without height interpolation.

For v=x(x-1), q_n(k)=U_k(n)*H_(n+2)+V_k(n). Polynomial recurrences
are solved exactly over Q. Continuation to s=1/2 gives rational + log(2).
This is NOT yet an analytic height function or a full hyperoperation.
"""
import argparse,csv,json,math,time
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).absolute().parent
ZERO=[F(0)]
def trim(p):
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def add(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,v in enumerate(p):r[i]+=v
    for i,v in enumerate(q):r[i]+=v
    return trim(r)
def scale(p,a):return trim([v*a for v in p])
def mul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q):r[i+j]+=v*w
    return trim(r)
def value(p,x):
    r=F(0)
    for c in reversed(p):r=r*x+c
    return r
def shift(p,a):
    r=[F(0)]*len(p)
    for n,c in enumerate(p):
        for j in range(n+1):r[j]+=c*math.comb(n,j)*a**(n-j)
    return trim(r)
def divide_linear(p,a):
    # Divide exactly by n+a.
    r=[F(0)]*(len(p)-1);carry=p[-1];r[-1]=carry
    for j in range(len(p)-2,0,-1):carry=p[j]-a*carry;r[j-1]=carry
    if p[0]-a*carry:raise ValueError('nonzero polynomial division remainder')
    return trim(r)
def rising(a,count):
    r=[F(1)]
    for i in range(count):r=mul(r,[F(a+i),F(1)])
    return r
def bernoulli(N):
    b=[F(1)]
    for n in range(1,N+1):b.append(-sum(F(math.comb(n+1,k))*b[k] for k in range(n))/F(n+1))
    return b
def antidifference(p,b):
    # P(n)=sum_(j=0)^(n-1) p(j), hence P(0)=0 and Delta P=p.
    r=[F(0)]*(len(p)+1)
    for m,c in enumerate(p):
        for j in range(m+1):r[m+1-j]+=c*F(math.comb(m+1,j),m+1)*b[j]
    return trim(r)
def harmonic(k):return sum((F(1,j) for j in range(1,k+1)),F(0))
def asmp(q):return mp.mpf(q.numerator)/q.denominator
def fractionstr(q):return f'{q.numerator}/{q.denominator}'
def compute(K,dps):
    mp.mp.dps=dps;start=time.monotonic();b=bernoulli(2*K+2)
    U={1:ZERO[:],2:ZERO[:]};V={1:ZERO[:],2:ZERO[:]}
    sumU=ZERO[:];sumV=ZERO[:];R=ZERO[:];S=ZERO[:];records=[]
    for k in range(3,K+1):
        P=scale(mul(rising(4,k-3),rising(3,k-2)),F(4*(k-1)*(k-2),math.factorial(k-1)*math.factorial(k-2)))
        R=add(R,P)
        correction=scale(P,-harmonic(k-2))
        for h in range(3,k+1):correction=add(correction,divide_linear(P,h))
        S=add(S,correction)
        uk=antidifference(add(sumU,R),b)
        uk[0]-=value(uk,F(-2))
        quotient=divide_linear(shift(uk,F(1)),3)
        vk=antidifference(add(add(sumV,S),scale(quotient,-1)),b)
        vk[0]=-value(uk,F(0))*harmonic(2)
        U[k]=uk;V[k]=vk
        # Verify polynomial identities, not just sampled values.
        assert add(shift(uk,F(1)),scale(uk,-1))==add(sumU,R)
        assert add(shift(vk,F(1)),scale(vk,-1))==add(add(sumV,S),scale(quotient,-1))
        assert value(uk,F(0))*harmonic(2)+value(vk,F(0))==0
        sumU=add(sumU,uk);sumV=add(sumV,vk)
        uh=value(uk,F(1,2));vh=value(vk,F(1,2));r=F(46,15)*uh+vh;c=-2*uh
        q=asmp(r)+asmp(c)*mp.log(2)
        records.append({'height':k,'degree_U':len(uk)-1,'degree_V':len(vk)-1,
           'half_rational':fractionstr(r),'half_log2_coefficient':fractionstr(c),
           'half_value':mp.nstr(q,dps-20),
           'U_coefficients':[fractionstr(x) for x in uk],
           'V_coefficients':[fractionstr(x) for x in vk]})
        print(k,mp.nstr(q,18),round(time.monotonic()-start,2),flush=True)
        if time.monotonic()-start>1080:raise TimeoutError('18 minute budget')
        if k%4==0:
            (HERE/f'data/intrinsic_K{K}_{dps}.json').write_text(json.dumps({'dps':dps,'K':K,'seconds':time.monotonic()-start,'records':records},indent=2),encoding='utf-8')
    out={'dps':dps,'K':K,'seconds':time.monotonic()-start,'records':records}
    (HERE/f'data/intrinsic_K{K}_{dps}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--K',type=int,default=32);p.add_argument('--dps',type=int,default=200)
    a=p.parse_args();compute(a.K,a.dps)
