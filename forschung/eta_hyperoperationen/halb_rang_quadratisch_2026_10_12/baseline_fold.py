"""Tangency of a declared convex interpolation, not an actual half rank.

F_theta=(1-theta)*b**x+theta*T_b(x). This gives a comparator only.
Each independent GP worker is closed and each run has a wall-time budget.
"""
import argparse,json,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).absolute().parent
sys.path.insert(0,str(HERE/'runtime'))
from common import engine,mp,stringify

p=argparse.ArgumentParser();p.add_argument('--theta',default='.5');p.add_argument('--dps',type=int,default=25)
p.add_argument('--samples',type=int,default=48);p.add_argument('--radius',default='.04')
p.add_argument('--max-seconds',type=int,default=600);a=p.parse_args()
mp.mp.dps=a.dps+35;start=time.monotonic();records=[]
out={'model':'convex interpolation, NOT H_(3+theta)', 'theta':a.theta,'requested_dps':a.dps,
     'samples':a.samples,'radius':a.radius,'records':records,'status':'running'}
target=HERE/f'data/baseline_theta_{a.theta.replace(".","p")}_{a.dps}.json'
def fmt(x):return mp.nstr(x,a.dps+20)
def save():
    out['seconds']=time.monotonic()-start
    target.write_text(json.dumps(out,indent=2),encoding='utf-8')
def minimum(b,x0):
    if time.monotonic()-start>a.max_seconds-60:raise TimeoutError('run budget reached')
    gp=engine(a.dps,HERE/'runtime/hyper.gp',timeout=90,even_extension=HERE/'runtime/baseline.gp')
    bs=fmt(b);t=time.monotonic()
    try:
        meta=gp.eval_batch(bs,['precis'])
        x=mp.mpf(x0);path=[]
        for k in range(30):
            vals=gp.eval_batch(bs,[f'mixjetset({a.theta},{fmt(x)},{a.radius},{a.samples})',
                                  f'polcoef(Pmix,0)-({fmt(x)})','polcoef(Pmix,1)-1','2*polcoef(Pmix,2)'])[1:]
            if max(abs(z.imag) for z in vals)>mp.mpf('1e-16'):raise ValueError('nonreal jet')
            g,d,d2=[z.real for z in vals]
            if d2<=0:raise ValueError('non-positive curvature')
            dx=max(-mp.mpf('.2'),min(mp.mpf('.2'),d/d2))
            path.append([fmt(x),fmt(g),fmt(d),fmt(d2)])
            if abs(dx)<mp.mpf(10)**(-min(a.dps-6,32)):break
            x-=dx
        else:raise ValueError('minimum not converged')
        g=gp.eval_batch(bs,[f'mix({a.theta},{fmt(x)})-({fmt(x)})'])[0].real
        r={'base':fmt(b),'height':fmt(x),'minimum':fmt(g),'derivative_minus_one':fmt(d),
           'curvature':fmt(d2),'actual_precis':stringify(meta[-1]),'trajectory':path,
           'seconds':time.monotonic()-t}
        records.append(r);save();print(a.theta,fmt(b),fmt(x),fmt(g),flush=True)
        return g,x
    finally:gp.close()

try:
    b0=mp.mpf('1.50');b1=mp.mpf('1.59');g0,x0=minimum(b0,'2.9');g1,x1=minimum(b1,x0)
    for j in range(12):
        b=b1-g1*(b1-b0)/(g1-g0);g,x=minimum(b,x1)
        b0,g0,b1,g1,x1=b1,g1,b,g,x
        if abs(g)<mp.mpf(10)**(-min(a.dps-6,30)):break
    out.update(status='completed',eta=fmt(b1),fold_height=fmt(x1),residual=fmt(g1))
except Exception as exc:
    out.update(status='failed',error=f'{type(exc).__name__}: {exc}');raise
finally:save()
