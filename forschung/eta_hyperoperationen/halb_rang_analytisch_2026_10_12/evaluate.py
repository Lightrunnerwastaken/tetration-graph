"""Evaluate the convergent quadratic Newton function; partial sums, not intervals."""
import json,csv,time,argparse
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).absolute().parent
def mpf(q):q=F(q);return mp.mpf(q.numerator)/q.denominator
def run(source,dps):
    mp.mp.dps=dps;start=time.monotonic();raw=json.loads(source.read_text(encoding='utf-8'))
    coeff=[mpf(r['rational'])+mpf(r['log2_coefficient'])*mp.log(2) for r in raw['coefficients']]
    K=len(coeff)-1;points=[mp.mpf('.5'),mp.mpf('1.5'),mp.mpf('2.5'),mp.mpc('.5','.25'),mp.mpc('1.5','.5')]
    rows=[]
    for x in points:
        total=mp.mpc(0);bj=mp.mpc(1)
        for j,a in enumerate(coeff):
            if j:bj*= (x-j+1)/j
            total+=bj*a
            if j in [16,32,64,96,128,192,256,384,512]:
                rows.append({'x_re':str(x.real),'x_im':str(x.imag),'cutoff':j,
                    'q_re':mp.nstr(total.real,70),'q_im':mp.nstr(total.imag,70),
                    'last_term_abs':mp.nstr(abs(bj*a),30)})
                print(x,j,mp.nstr(total,23),flush=True)
    result={'source':source.name,'dps':dps,'K':K,'seconds':time.monotonic()-start,'samples':rows}
    (HERE/f'data/evaluation_K{K}_{dps}.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--dps',type=int,default=400);a=p.parse_args();run(a.source,a.dps)
