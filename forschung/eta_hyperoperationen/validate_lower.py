"""Check inherited tangencies with complex Cauchy jets, without modifying sources."""
import json,time
from common import HERE,engine,mp,stringify
mp.mp.dps=50
rows=json.loads((HERE/'prior_values.json').read_text())['values'][1:]
records=[];start=time.time()
for row in rows:
    r=row['rank'];b=mp.mpf(row['eta']);x=mp.mpf(row['e']);gp=engine(25,HERE/'hyper.gp')
    try:
        vals=gp.eval_batch(row['eta'],[f'hinit({r})',f'hjetset({r},{row["e"]},.025,48)',
                     *[f'hlastjet[{j}]' for j in [1,2,3]],f'hval({r},{row["e"]})-({row["e"]})'])
        records.append({'rank':r,'base':row['eta'],'x':row['e'],
             'expressions':'hinit, jetset, F-x, Fprime-1, Fsecond, direct F-x',
             'values':list(map(stringify,vals)), 'seconds':time.time()-start})
        print(r,list(map(stringify,vals)),flush=True)
    except Exception as exc:
        records.append({'rank':r,'error':str(exc),'seconds':time.time()-start})
    finally:
        gp.close();(HERE/'lower_validation.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
# Measure the effect of a permissible periodic change of the height parameter.
r=4;x=mp.mpf(rows[0]['e']);b=mp.mpf(rows[0]['eta']);step=mp.mpf('1e-7');vals=[]
for sign in [-1,1]:
    gp=engine(30)
    try:vals.append(gp.eval_batch(mp.nstr(b+sign*step,45),[f'sexp({mp.nstr(x,45)})'])[0].real)
    finally:gp.close()
alpha=(vals[1]-vals[0])/(2*step)
(HERE/'selection_sensitivity.json').write_text(json.dumps({'rank':4,'base_step':str(step),
    'alpha_base_derivative':str(alpha),'eta_prime_periodic_warp':str(-mp.sin(2*mp.pi*x)/alpha),
    'status':'numerical central base difference, no interval guarantee'},indent=2),encoding='utf-8')
