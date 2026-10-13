"""Exact Legendre generating coefficients, without the growing U,V arrays."""
import argparse,json,math,time
from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
from intrinsic import bernoulli,add,mul,scale,asmp
HERE=Path(__file__).absolute().parent
def run(K,dps):
    start=time.monotonic();mp.mp.dps=dps;b=bernoulli(2*K+2)
    moments=[F(1,2)]+[(F(1,2**l)-2)*b[l+1]/(l+1)+F(1,2)**l+F(3,2)**l-1 for l in range(1,2*K+1)]
    h=[F(1)];p=[F(0),F(1)];W=[None,F(7,8)];H=F(0);dr=[F(0)];dc=[F(0)]
    for j in range(1,K-1):
        h.append(h[-1]*(F(5,2)-j+1)*(F(5,2)+j)/F(j*j));H+=F(1,j)
        w=sum((c*moments[l] for l,c in enumerate(p)),F(0))
        if j>=len(W):W.append(w)
        else:assert W[j]==w
        dr.append(h[j]*(F(46,15)-H)-(F(3) if j==1 else 0)-w)
        dc.append(-2*h[j])
        p=scale(mul(p,[F(-j*(j+1)),F(1),F(1)]),F(1,(j+1)**2))
        if time.monotonic()-start>1080:raise TimeoutError('18 minute budget')
    bc=[F(1)]
    for j in range(1,K-1):bc.append(bc[-1]*(F(1,2)-j+1)/j)
    rows=[{'j':j,'rational':'0','log2_coefficient':'0','value':'0'} for j in range(2)]
    for j in range(2,K+1):
        r=2*sum((bc[k]*dr[j-2-k] for k in range(j-1)),F(0))
        c=2*sum((bc[k]*dc[j-2-k] for k in range(j-1)),F(0))
        val=asmp(r)+asmp(c)*mp.log(2)
        rows.append({'j':j,'rational':str(r),'log2_coefficient':str(c),'value':mp.nstr(val,dps-20)})
    out={'K':K,'dps':dps,'seconds':time.monotonic()-start,'coefficients':rows}
    (HERE/f'data/direct_coefficients_K{K}_{dps}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print('K',K,'seconds',round(out['seconds'],3),'last',mp.nstr(val,24),flush=True)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--K',type=int,default=256);p.add_argument('--dps',type=int,default=400);a=p.parse_args();run(a.K,a.dps)
