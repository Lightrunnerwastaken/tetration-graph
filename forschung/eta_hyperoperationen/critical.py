"""Selected-hierarchy positive tangencies; midpoint Cauchy jets, separate sessions."""
import argparse,json,sys,time
from common import engine,HERE,mp,stringify
sys.stdout.reconfigure(encoding="utf-8")
p=argparse.ArgumentParser()
p.add_argument("--rank",type=int,required=True);p.add_argument("--b0",required=True);p.add_argument("--b1",required=True)
p.add_argument("--dps",type=int,default=25);p.add_argument("--samples",type=int,default=48)
p.add_argument("--radius",default=".04");p.add_argument("--x0",default="3.55")
p.add_argument("--steps",type=int,default=8);p.add_argument("--output",required=True)
p.add_argument("--root-digits",type=int,default=0)
p.add_argument("--even",default="even8.gp")
p.add_argument("--resume",default="")
p.add_argument("--eval-seconds",type=int,default=300)
p.add_argument("--max-seconds",type=int,default=1100)
a=p.parse_args();mp.mp.dps=a.dps+25;start=time.time();records=[]
goal=a.root_digits or min(a.dps-7,24)
out={"status":"numerical local positive tangency, not interval or universal hierarchy",
     "rank":a.rank,"requested_dps":a.dps,"samples":a.samples,"radius":a.radius,"root_digits":goal,"records":records}
def save():
    out["seconds"]=time.time()-start
    (HERE/a.output).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
def fmt(v):return mp.nstr(v,a.dps+15)
def minimum(b,x0):
    remaining=a.max_seconds-(time.time()-start)
    if remaining<200:raise TimeoutError("session budget exceeded; preceding values saved")
    gp=engine(a.dps,HERE/"hyper.gp",logged=a.rank>=8,timeout=min(a.eval_seconds,remaining-170),even_extension=HERE/a.even);ts=time.time()
    try:
        bs=fmt(b);meta=gp.eval_batch(bs,[f"hinit({a.rank})","precis","hexreg"])
        state=gp.eval_batch(bs,[f"hrkind[{a.rank}]",f"hrp[{a.rank}]",f"hrlam[{a.rank}]"]) if a.rank>=7 else []
        x=mp.mpf(x0);trajectory=[]
        for k in range(14):
            expr=f"hjetset({a.rank},{fmt(x)},{a.radius},{a.samples})"
            vals=gp.eval_batch(bs,[expr,*[f"hlastjet[{i}]" for i in [1,2,3]]])[1:]
            if max(abs(v.imag) for v in vals)>mp.mpf("1e-18"):raise ValueError("nonreal Cauchy jet")
            g,d1,d2=[v.real for v in vals];assert d2>0
            dx=max(-mp.mpf(".2"),min(mp.mpf(".2"),d1/d2))
            trajectory.append([fmt(x),fmt(g),fmt(d1),fmt(d2)])
            if abs(dx)<mp.mpf(10)**(-goal):break
            x-=dx
        else:raise RuntimeError("minimum Newton did not converge")
        g=gp.eval_batch(bs,[f"hval({a.rank},{fmt(x)})-({fmt(x)})"])[0].real
        rec={"base":fmt(b),"x":fmt(x),"minimum":fmt(g),"derivative_minus_one":fmt(d1),
             "second_derivative":fmt(d2),"engine_digits":stringify(meta[0]),
             "actual_precis":stringify(meta[1]),"hexreg":stringify(meta[2]),
             "hierarchy_state":list(map(stringify,state)),
             "minimum_newton":trajectory,"seconds":time.time()-ts}
        records.append(rec);save()
        print(f"r={a.rank}, b={fmt(b)}, x={mp.nstr(x,18)}, min={mp.nstr(g,6)}, {time.time()-start:.1f}s",flush=True)
        return g,x
    finally:gp.close()
try:
    if a.resume:
        old=json.loads((HERE/a.resume).read_text(encoding='utf-8'))
        seeds=old['records'][-2:];out['seed_source']=a.resume;out['seed_records']=seeds
        b0,b1=[mp.mpf(z['base']) for z in seeds];g0,g1=[mp.mpf(z['minimum']) for z in seeds]
        x1=mp.mpf(seeds[-1]['x'])
    else:
        b0,b1=map(mp.mpf,[a.b0,a.b1]);g0,x0=minimum(b0,a.x0);g1,x1=minimum(b1,x0)
    for j in range(a.steps):
        b=b1-g1*(b1-b0)/(g1-g0);g,x=minimum(b,x1)
        b0,g0,b1,g1,x1=b1,g1,b,g,x
        if abs(g)<mp.mpf(10)**(-goal):break
    out.update(completed=bool(abs(g1)<mp.mpf(10)**(-goal)),eta=fmt(b1),euler=fmt(x1),minimum=fmt(g1));save()
except Exception as exc:
    out.update(completed=False,error=f'{type(exc).__name__}: {exc}');raise
finally:save()
