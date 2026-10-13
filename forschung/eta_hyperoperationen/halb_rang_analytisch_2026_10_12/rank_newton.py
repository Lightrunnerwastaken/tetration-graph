"""Reverse the order of continuation: true polynomial integer-rank jets."""
import argparse,json,math,time
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
from intrinsic import add,scale,mul,value,shift,rising,bernoulli,antidifference,asmp
HERE=Path(__file__).absolute().parent
def derivative(p):return [F(j)*p[j] for j in range(1,len(p))] or [F(0)]
def summation(p,b):
    q=antidifference(p,b);q=shift(q,F(1));q[0]-=value(antidifference(p,b),F(1));return q
def run(N,dps):
    mp.mp.dps=dps;start=time.monotonic();b=bernoulli(2*N+8);q=[F(0)];rows=[]
    points=[F(1,2),F(3,2),F(5,2),F(3)]
    vals={x:[] for x in points}
    for n in range(N+1):
        if time.monotonic()-start>1080:raise TimeoutError('18 minute budget')
        for x in points:vals[x].append(value(q,x))
        rows.append({'n':n,'degree':len(q)-1,'values':{str(x):str(vals[x][-1]) for x in points}})
        print(n,mp.nstr(asmp(vals[F(1,2)][-1]),18),flush=True)
        a=scale(mul([F(-1),F(1)],rising(0,n+1)),F(2,math.factorial(n+2)))
        w=summation(a,b);source=mul(add(w,scale(a,-1)),derivative(a))
        q=add(summation(q,b),summation(source,b))
    samples=[]
    for x,data in vals.items():
        coeff=[]
        while data:
            coeff.append(data[0]);data=[d-c for c,d in zip(data,data[1:])]
        total=mp.mpf(0);bj=mp.mpf(1)
        for j,c in enumerate(coeff):
            if j:bj*= (mp.mpf('.5')-j+1)/j
            total+=bj*asmp(c)
            if j in [4,8,12,16,24,32,48,64,96,128]:
                samples.append({'x':str(x),'cutoff':j,'value':mp.nstr(total,60),'last_term_abs':mp.nstr(abs(bj*asmp(c)),30)})
                print('partial',x,j,mp.nstr(total,24),mp.nstr(abs(bj*asmp(c)),10),flush=True)
    out={'N':N,'dps':dps,'seconds':time.monotonic()-start,'integer_rank':rows,'samples':samples}
    (HERE/f'data/rank_newton_N{N}_{dps}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=32);p.add_argument('--dps',type=int,default=200);a=p.parse_args();run(a.N,a.dps)
