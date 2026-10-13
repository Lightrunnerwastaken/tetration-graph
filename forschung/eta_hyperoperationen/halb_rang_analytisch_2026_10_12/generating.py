"""Legendre/Jacobi generating functions for the intrinsic quadratic jet."""
import json,time,argparse
from pathlib import Path
import mpmath as mp
from abel_plana import ap_sum
HERE=Path(__file__).absolute().parent
def P(n,t):return mp.legenp(n,0,t,type=3)
def J(n,t):return P(n,t)*mp.log((t+1)/(t-1))/2-mp.legenq(n,0,t,type=3)
def A_U(s,w):return 2*w*w*mp.power(1+w,s)*(P(s+2,1+2*w)-1)
def A_q(s,w):
    t=1+2*w
    T=mp.ceil((mp.mp.dps+20)*mp.log(10)/(2*mp.pi))
    W=ap_sum(lambda m:P(m,t)/(m+1),s+1,T=T)-t/2
    return 2*w*w*mp.power(1+w,s)*(J(s+2,t)-3*t/2-W)
def run(dps):
    start=time.monotonic();mp.mp.dps=dps
    data=json.loads((HERE/'data/height_newton_K128_500.json').read_text(encoding='utf-8'))
    rows=[]
    for w in [mp.mpf('.1'),mp.mpf('.3')]:
        aq=mp.fsum(mp.mpf(r['value'])*w**r['j'] for r in data['coefficients'])
        au=mp.fsum(-mp.mpf(r['log2_coefficient'].split('/')[0])/mp.mpf(r['log2_coefficient'].split('/')[1] if '/' in r['log2_coefficient'] else 1)*w**r['j']/2 for r in data['coefficients'])
        gu=A_U(mp.mpf('.5'),w);gq=A_q(mp.mpf('.5'),w)
        t=1+2*w;T=mp.ceil((dps+20)*mp.log(10)/(2*mp.pi));rmax=t+mp.sqrt(t*t-1)
        tail=2*(1+rmax**mp.mpf('1.5'))*mp.exp(-2*mp.pi*T)/(2*mp.pi*(1-mp.exp(-2*mp.pi*T)))
        tail*=2*w*w*mp.sqrt(1+w)
        rows.append({'w':str(w),'T':str(T),'AP_tail_bound':mp.nstr(tail,15),'U_expected':mp.nstr(au,dps-5),'U_formula':mp.nstr(gu,dps-5),
            'U_error':mp.nstr(abs(au-gu),15),'q_expected':mp.nstr(aq,dps-5),
            'q_formula_re':mp.nstr(gq.real,dps-5),'q_formula_im':mp.nstr(gq.imag,dps-5),
            'q_error':mp.nstr(abs(aq-gq),15)})
        print(w,mp.nstr(abs(au-gu),8),mp.nstr(abs(aq-gq),8),round(time.monotonic()-start,2),flush=True)
    out={'dps':dps,'seconds':time.monotonic()-start,'records':rows}
    (HERE/f'data/generating_{dps}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dps',type=int,default=40);a=p.parse_args();run(a.dps)
