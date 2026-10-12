"""Deterministic model comparisons, not statistical confidence intervals."""
import json,math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import HERE,mp
mp.mp.dps=65
rows=json.loads((HERE/'prior_values.json').read_text())['values']
rows[0].update(eta=str(mp.exp(1/mp.e)),e=str(mp.e))
for r,filename in [(7,'eta7_45.json'),(8,'eta8_pilot.json'),(8,'eta8_refined.json'),(8,'eta8_control.json'),(8,'eta8_control_retry.json')]:
    f=HERE/filename
    if f.exists():
        dat=json.loads(f.read_text())
        if dat.get('completed'):
            rows=[z for z in rows if z['rank']!=r]
            rows.append({'rank':r,'eta':dat['eta'],'e':dat['euler'],'source':filename})
rows.sort(key=lambda z:z['rank'])
for i,z in enumerate(rows):
    gap=2-mp.mpf(z['eta']);z['gap_to_2']=str(gap)
    z['gap_times_r_to_3_over_2']=str(gap*mp.mpf(z['rank'])**mp.mpf('1.5'))
    if i:
        prev=rows[i-1];z['effective_power']=str(mp.log(mp.mpf(prev['gap_to_2'])/gap)/mp.log(mp.mpf(z['rank'])/prev['rank']))
        z['e_previous_minus_eta_squared']=str(mp.mpf(prev['e'])-mp.mpf(z['eta'])**2)
n=np.array([z['rank'] for z in rows],float);eta=np.array([float(z['eta']) for z in rows]);e=np.array([float(z['e']) for z in rows])
train=n<=6;slope,intercept=np.polyfit(np.log(n[train]),np.log(2-eta[train]),1)
A,p=float(np.exp(intercept)),float(-slope)
heldout=[{'rank':int(k),'prediction':float(2-A*k**(-p)),'measured':float(v),
          'error_prediction_minus_measured':float(2-A*k**(-p)-v)} for k,v in zip(n[~train],eta[~train])]
models=[]
funcs=[('L-A*r^(-p)',lambda r,L,A,p:L-A*r**(-p),[2,2.8,1.47],([1.84,.001,.01],[3,100,8])),
       ('L-A*exp(-k*r)',lambda r,L,A,k:L-A*np.exp(-k*r),[1.95,1,.4],([1.84,.001,.001],[3,100,5]))]
for name,fn,guess,bounds in funcs:
    par,_=curve_fit(fn,n,eta,p0=guess,bounds=bounds,maxfev=20000)
    models.append({'model':name,'parameters':par.tolist(),'limit':float(par[0]),
                   'max_residual':float(np.max(abs(fn(n,*par)-eta))),
                   'status':'descriptive fitted family, not a limit proof'})
out={'values':rows,'fit_to_old_ranks_3_to_6':{'A':A,'p':p,'limit_fixed':2,'heldout_predictions':heldout},'models_free_limit':models}
(HERE/'analysis.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
axes[0].plot(n,eta,'o-',label='ηᵣ (numerisch / η₃ exakt)',color='#235587')
axes[0].axhline(2,color='#c97930',linestyle='--',label='Vermuteter Grenzwert 2')
axes[0].set(xlabel='Rang r der berührenden Abbildung Hᵣ',ylabel='Berührbasis ηᵣ',ylim=(1.4,2.025),xticks=n)
axes[0].grid(alpha=.2);axes[0].legend(fontsize=9)
axes[1].loglog(n,2-eta,'o',label='Gemessener Abstand 2−ηᵣ',color='#235587')
nn=np.linspace(3,10,200)
axes[1].loglog(nn,A*nn**(-p),'--',label=f'Fit nur r=3…6: {A:.4f} r⁻{p:.4f}',color='#c97930')
axes[1].set(xlabel='Rang r',ylabel='Abstand zu 2',xticks=np.arange(3,11));axes[1].set_xticklabels(range(3,11))
axes[1].grid(alpha=.2);axes[1].legend(fontsize=9)
fig.suptitle('Berührbasen der gewählten Hyperoperations-Hierarchie\nηᵣ ist Schwellenkandidat für die nächste Stufe Hᵣ₊₁',fontsize=13)
fig.savefig(HERE/'eta_folge.png',dpi=180);fig.savefig(HERE/'eta_folge.svg');plt.close(fig)
print(json.dumps(out,indent=2))
