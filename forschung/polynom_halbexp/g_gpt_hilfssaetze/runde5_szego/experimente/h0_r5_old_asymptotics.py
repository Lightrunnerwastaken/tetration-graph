"""Stage 2 preliminary analysis of already completed every-degree runs.

All fits are numerical and windows are reported; no inferred all-n error bound.
"""
import json,math,sys
from pathlib import Path
from h0_r5_paths import locations
import numpy as np
import mpmath as mp
sys.stdout.reconfigure(encoding='utf-8');mp.mp.dps=70
OUT,ROUND4,ROOT=locations(__file__);OUT.mkdir(exist_ok=True)
report={}
for model in ['kubisch','kneser_z']:
    d=json.loads((ROUND4/f'h0_schur_{model}.json').read_text(encoding='utf-8'))
    cases=d['cases'];N=cases[-1]['n'];h={c['n']:mp.mpf(c['h']) for c in cases}
    ratio={n:mp.log(h[n]/h[n-1]) for n in h if n-1 in h}
    fits=[]
    for low in ([200,400,600,800] if N==1000 else [150,250,350,450]):
        for powers in [[0,1],[0,2],[0,1,2],[0,1,2,3],[0,2,3,4],[0,1,2,3,4]]:
            ns=list(range(low,N+1));A=mp.matrix([[mp.mpf(n)**(-p) for p in powers] for n in ns]);b=mp.matrix([ratio[n] for n in ns])
            x,res=mp.qr_solve(A,b);residual=max(abs(sum(x[k]*mp.mpf(n)**(-p) for k,p in enumerate(powers))-ratio[n]) for n in ns)
            fits.append({'low':low,'high':N,'powers':powers,'log_coefficients':[str(v) for v in x],
                         'lambda':str(mp.exp(x[0])),'max_log_ratio_residual':str(residual)})
    # Log Q and log(Q/h), model alpha*n + b*log(n) + c + d/n.
    qfits=[]
    for low in ([300,500,700] if N==1000 else [200,300,400]):
        chosen=[c for c in cases if c['n']>=low];ns=np.array([c['n'] for c in chosen],float)
        X=np.column_stack([ns,np.log(ns),np.ones_like(ns),1/ns])
        for name,vals in [('Q',[math.log(c['Q_norm']) for c in chosen]),('Q_over_h',[math.log(c['Q_norm'])-float(mp.log(h[c['n']])) for c in chosen])]:
            co,*_=np.linalg.lstsq(X,np.array(vals),rcond=None)
            qfits.append({'name':name,'low':low,'high':N,'coefficients_n_log_n_1_inverse_n':co.tolist(),
                          'max_log_residual':float(np.max(abs(X@co-vals)))})
    report[model]={'N':N,'ratio_fits':fits,'norm_fits':qfits,'tail_ratios':{n:str(mp.exp(ratio[n])) for n in [N-100,N-50,N]},
                   'status':'numerical fits, not proof; window dependence retained'}
(OUT/'old_asymptotics.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
for name,v in report.items():
    print(name, v['tail_ratios'])
    for fit in v['ratio_fits']:
        if fit['powers']==[0,1,2,3]:print(fit)
