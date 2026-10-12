import json,sys,time
from common import engine,HERE,mp,stringify
sys.stdout.reconfigure(encoding="utf-8");mp.mp.dps=45
base=sys.argv[1] if len(sys.argv)>1 else "2"
dps=int(sys.argv[2]) if len(sys.argv)>2 else 25
start=time.time();gp=engine(dps);log=[]
try:
    expressions=["hexensure()","hexreg","hexerr","hex(-3.5)","ihex(-3.5)"]
    values=gp.eval_batch(base,expressions)
    log.append({"expressions":expressions,"values":list(map(stringify,values))})
    print(json.dumps(log[-1],indent=2),flush=True)
    values=gp.eval_batch(base,["my(y=1);for(j=1,20,y=ihex(y));y","knc","kngood",
                              "hex(-3.7)","hex(-3.8)"])
    log.append({"values":list(map(stringify,values))})
    print(json.dumps(log[-1],indent=2),flush=True)
finally:
    gp.close()
    (HERE/f"probe_hex_{base.replace('.','p')}_{dps}.json").write_text(json.dumps(
        {"base":base,"dps":dps,"seconds":time.time()-start,"log":log},indent=2),encoding="utf-8")
