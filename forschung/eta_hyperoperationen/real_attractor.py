"""Real Koenigs approximants for bounded decation, one bounded session."""
import json,time
from common import HERE,engine,mp,stringify
mp.mp.dps=50
ext=HERE/'hyper_real.gp';ext.write_bytes((HERE/'hyper.gp').read_bytes()+b'\n'+(HERE/'real_attractor.gp').read_bytes())
gp=engine(25,ext,logged=True,timeout=750);start=time.time();out={'rank':10,'base':'1.83','status':'running','samples':[]}
def save():
    out['seconds']=time.time()-start;(HERE/'real_decation.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
try:
    vals=gp.eval_batch('1.83',['hrdebug=1','haready(10)','hla','hlaerr'])
    out['attractor_multiplier_checks']=list(map(stringify,vals));save();print(out,flush=True)
    for k in [5,7,9]:
        if time.time()-start>650:raise TimeoutError('session limit; preceding values saved')
        gp._worker_for('1.83').eval_timeout=min(300,1000-(time.time()-start))
        vals=gp.eval_batch('1.83',[f'hrealset(10,.5,{k})','halflast','hval(9,halflast)',f'hiter(9,1,{k})-ha'])[1:]
        out['samples'].append({'k':k,'values':list(map(stringify,vals))});save();print(out['samples'][-1],flush=True)
    out['status']='completed finite real Koenigs approximants, not a complex construction'
except Exception as exc:out.update(status='failed',error=f'{type(exc).__name__}: {exc}')
finally:gp.close();save()
