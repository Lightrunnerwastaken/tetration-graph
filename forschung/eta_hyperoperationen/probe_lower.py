"""Uniform values at base 2, ranks 3..6, from the frozen engine."""
import json,time
from common import engine,HERE,mp,stringify
mp.mp.dps=50;gp=engine(25,HERE/'hyper.gp');start=time.time();out={}
try:
    expr=['hinit(6)','precis','hexdigits']+[f'hval({r},{x})' for r in range(3,7) for x in ['.5','1.5']]
    vals=gp.eval_batch('2',expr)
    out={'base':'2','expressions':expr,'values':list(map(stringify,vals)),'status':'completed'}
except Exception as exc:out={'status':'failed','error':str(exc)}
finally:
    gp.close();out['seconds']=time.time()-start;(HERE/'probe_lower_2.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(out)
