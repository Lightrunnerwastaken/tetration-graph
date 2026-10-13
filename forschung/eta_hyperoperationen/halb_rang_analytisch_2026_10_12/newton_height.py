"""Exact height differences, followed by high precision Newton evaluation."""
import csv,json,math,time
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).absolute().parent
def asmp(q):return mp.mpf(q.numerator)/q.denominator
def run(source,dps=300):
    start=time.monotonic();mp.mp.dps=dps
    raw=json.loads(source.read_text(encoding='utf-8'));K=raw['records'][-1]['height']
    rs=[F(0)]*(K+1);cs=[F(0)]*(K+1)
    for r in raw['records']:
        k=r['height'];rs[k]=F(r['half_rational']);cs[k]=F(r['half_log2_coefficient'])
    ar=[];ac=[]
    while rs:
        ar.append(rs[0]);ac.append(cs[0])
        rs=[b-a for a,b in zip(rs,rs[1:])];cs=[b-a for a,b in zip(cs,cs[1:])]
    coeff=[asmp(r)+asmp(c)*mp.log(2) for r,c in zip(ar,ac)]
    out={'source':source.name,'K':K,'dps':dps,'coefficients':[
        {'j':j,'rational':str(r),'log2_coefficient':str(c),'value':mp.nstr(coeff[j],dps-20)}
        for j,(r,c) in enumerate(zip(ar,ac))],'samples':[]}
    points=[mp.mpf('.5'),mp.mpf('1.5'),mp.mpf('2.5'),mp.mpc('.5','.25'),mp.mpc('1.5','.5')]
    cuts=[n for n in [8,16,24,32,48,64,80,96,128] if n<=K]
    for x in points:
        total=mp.mpc(0);bj=mp.mpc(1)
        for j in range(K+1):
            if j:bj*= (x-j+1)/j
            total+=bj*coeff[j]
            if j in cuts:
                r={'x_re':mp.nstr(x.real),'x_im':mp.nstr(x.imag),'cutoff':j,
                  'q_re':mp.nstr(total.real,70),'q_im':mp.nstr(total.imag,70),
                  'last_term_abs':mp.nstr(abs(bj*coeff[j]),70)}
                out['samples'].append(r);print(x,j,mp.nstr(total,24),mp.nstr(abs(bj*coeff[j]),8),flush=True)
    out['seconds']=time.monotonic()-start
    (HERE/f'data/height_newton_K{K}_{dps}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    with (HERE/f'data/height_newton_K{K}_{dps}.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(out['samples'][0]));w.writeheader();w.writerows(out['samples'])
    return out
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--dps',type=int,default=300)
    a=p.parse_args();run(a.source,a.dps)
