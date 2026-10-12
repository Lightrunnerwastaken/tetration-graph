"""Check the required rank recursion for the natural convex family."""
import json,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).absolute().parent
sys.path.insert(0,str(HERE/'runtime'))
from common import engine,mp,stringify
mp.mp.dps=60;start=time.monotonic();records=[]
gp=engine(30,HERE/'runtime/hyper.gp',timeout=120,even_extension=HERE/'runtime/baseline.gp')
try:
    for b in ['2','exp(1)']:
        gp.eval_batch(b,['hinit(5)'])
        for theta in ['0','.5','1']:
            for x in ['.5']:
                # G_theta=(1-theta) H4 + theta H5 is the next-rank comparator.
                y=f'((1-({theta}))*hval(4,{x})+({theta})*hval(5,{x}))'
                left=f'(1-({theta}))*hval(4,({x})+1)+({theta})*hval(5,({x})+1)'
                right=f'mix({theta},{y})'
                vals=gp.eval_batch(b,[y,left,right,f'({left})-({right})'])
                records.append({'base':b,'theta':theta,'x':x,'G_theta_x':stringify(vals[0],45),
                                'left':stringify(vals[1],45),'right':stringify(vals[2],45),
                                'rank_recursion_residual':stringify(vals[3],45)})
finally:
    gp.close()
    (HERE/'data/baseline_recursion.json').write_text(json.dumps({'status':'comparator, not half rank',
            'seconds':time.monotonic()-start,'records':records},indent=2),encoding='utf-8')
print(json.dumps(records,indent=2),flush=True)
