"""Held-out precision fits and fit-free ratio convergence diagnostics."""
import json,math,sys
from pathlib import Path
from h0_r5_paths import locations
import numpy as np
import mpmath as mp
sys.stdout.reconfigure(encoding='utf-8');mp.mp.dps=65
folder,ROUND4,ROOT=locations(__file__);out={}
for model in ['kubisch','kneser_z']:
    data=json.loads((ROUND4/f'h0_schur_{model}.json').read_text());by={c['n']:c for c in data['cases']};N=max(by)
    hs={n:mp.mpf(c['h']) for n,c in by.items()};rs={n:mp.log(hs[n]/hs[n-1]) for n in by if n-1 in by}
    ratio_tests=[]
    for n in [100,150,200,250]:
        if 4*n<=N:
            p=mp.log((rs[2*n]-rs[n])/(rs[4*n]-rs[2*n]))/mp.log(2)
            c2=(rs[n]-rs[2*n])*4*n*n/3
            ratio_tests.append({'n':n,'effective_order':str(p),'c2_Richardson':str(c2)})
    fitmax=700 if N==1000 else 450;low=400 if N==1000 else 250
    powers=[0,2,3,4];ns=list(range(low,fitmax+1));A=mp.matrix([[mp.mpf(n)**(-p) for p in powers] for n in ns]);b=mp.matrix([rs[n] for n in ns]);co,_=mp.qr_solve(A,b)
    err=max(abs(sum(co[j]*mp.mpf(n)**(-p) for j,p in enumerate(powers))-rs[n]) for n in range(fitmax+1,N+1))
    # Norm fits use scaled inverse powers to keep conditioning reasonable.
    normfits=[]
    for lowq in ([300,500,700] if N==1000 else [200,300,400]):
        chosen=[by[n] for n in range(lowq,N+1)];nv=np.array([c['n'] for c in chosen],float);t=nv/N
        for order in [1,2,3]:
            X=np.column_stack([t,np.log(t),np.ones_like(t)]+[t**(-k) for k in range(1,order+1)])
            for name,vals in [('Q',[math.log(c['Q_norm']) for c in chosen]),('Q_over_h',[math.log(c['Q_norm'])-float(mp.log(hs[c['n']])) for c in chosen])]:
                cf,*_=np.linalg.lstsq(X,np.array(vals),rcond=None)
                normfits.append({'name':name,'low':lowq,'high':N,'inverse_order':order,'exponential_rate':cf[0]/N,
                                 'log_n_power':cf[1],'max_log_residual':float(max(abs(X@cf-vals))),
                                 'scaled_fit_coefficients':cf.tolist()})
    out[model]={'ratio_fit_free_tests':ratio_tests,'heldout':{'train':[low,fitmax],'test':[fitmax+1,N],
                'powers':powers,'coefficients':[str(v) for v in co],'lambda':str(mp.exp(co[0])),
                'max_heldout_log_ratio_error':str(err)},'norm_fits':normfits,'status':'numerical, not all-n estimates'}
    print(model,ratio_tests,'heldout error=',str(err),flush=True)
    for f in normfits:
        if f['name']=='Q_over_h' and f['inverse_order']==3:print(f)
(folder/'second_order.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
