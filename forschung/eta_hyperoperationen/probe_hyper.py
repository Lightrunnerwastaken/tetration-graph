import argparse,json,sys,time
from common import engine,HERE,mp,stringify
sys.stdout.reconfigure(encoding="utf-8")
p=argparse.ArgumentParser();p.add_argument("--base",default="2");p.add_argument("--rank",type=int,default=7)
p.add_argument("--dps",type=int,default=25);p.add_argument("--timeout",type=int,default=300)
p.add_argument("--tag",default="");p.add_argument("--even",default="even8.gp");a=p.parse_args();mp.mp.dps=a.dps+20
start=time.time();gp=engine(a.dps,HERE/"hyper.gp",logged=True,timeout=a.timeout,even_extension=HERE/a.even);log=[];status="failed";error=None
try:
    expr=["hrdebug=1",f"hinit({a.rank})"]
    values=gp.eval_batch(a.base,expr)
    print("init",list(map(stringify,values)),"seconds",time.time()-start,flush=True)
    expr=["precis",f"hrkind[{a.rank}]",f"hrerror[{a.rank}]",f"hrp[{a.rank}]",f"hrlam[{a.rank}]",f"hrtarget[{a.rank}]",f"hrdepth[{a.rank}]",
          *[f"hval({a.rank},{x})" for x in ["-2","-1","0",".5","1","1.5","2"]],
          *[f"hval({a.rank},{x}+1)-hval({a.rank-1},hval({a.rank},{x}))" for x in [".25",".5"]],
          f"hval({a.rank},.5)-hval({a.rank},.5,2)"]
    values=gp.eval_batch(a.base,expr)
    log.append({"expressions":expr,"values":list(map(stringify,values))})
    print(json.dumps(log[-1],indent=2),flush=True);status="completed"
except Exception as exc:
    error=f"{type(exc).__name__}: {exc}";print(error,flush=True)
finally:
    gp.close()
    (HERE/f"probe_hyper_r{a.rank}_{a.base.replace('.','p')}_{a.dps}{a.tag}.json").write_text(
        json.dumps({"base":a.base,"rank":a.rank,"dps":a.dps,"status":status,
                    "seconds":time.time()-start,"error":error,"log":log},indent=2),encoding="utf-8")
