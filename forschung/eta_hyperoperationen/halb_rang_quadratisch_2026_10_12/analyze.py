"""Summarize saved results without starting computational workers."""
import csv,json,math,time
from pathlib import Path
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).absolute().parent
mp.mp.dps=130
def fmt(x):return mp.nstr(x,70)
raw=[];run_times=[]
for f in sorted((HERE/'data').glob('quadratic_*.json')):
    file_results=json.loads(f.read_text(encoding='utf-8'))
    run_times.append(sum(r['seconds'] for r in file_results))
    for r in file_results:
        r['source']=f.name;raw.append(r)
selected={}
for r in raw:
    if r['K'] not in selected or r['dps']>selected[r['K']]['dps']:selected[r['K']]=r
exact=mp.mpf(341)/12-35*mp.log(2)
rows=[]
for K,r in sorted(selected.items()):
    q=mp.mpf(r['half_quadratic'][2]);p=r['noninteger_values']
    rows.append({'K':K,'dps':r['dps'],'Q_half_x3':fmt(q),'x3_error':fmt(exact-q),
      'x3_tail_bound':fmt(8/(3*mp.sqrt(mp.pi))*(K-3)**(-mp.mpf('4.5'))) if K>=4 else '',
      'Q_half_x0p5':p[0]['quadratic_interpolation'],'Q_half_x1p5':p[1]['quadratic_interpolation'],
      'Q_half_x2p5':p[2]['quadratic_interpolation'],'linear_x0p5_error':p[0]['linear_error'],
      'cohomology_relative_error':r['cohomology_relative_error'],'source':r['source']})
with (HERE/'data/model_size_audit.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
controls=[]
for K in [4,8,16,32,64]:
    vals=sorted([r for r in raw if r['K']==K],key=lambda r:r['dps'])
    if len(vals)<2:continue
    a,b=vals[0],vals[-1]
    differences=[abs(mp.mpf(u)-mp.mpf(v)) for u,v in zip(a['half_quadratic'],b['half_quadratic'])]
    frac=[abs(mp.mpf(u['quadratic_interpolation'])-mp.mpf(v['quadratic_interpolation'])) for u,v in zip(a['noninteger_values'],b['noninteger_values'])]
    controls.append({'K':K,'dps_low':a['dps'],'dps_high':b['dps'],'added_bits':(b['dps']-a['dps'])*math.log2(10),
                     'max_integer_absolute_difference':fmt(max(differences)),
                     'max_noninteger_absolute_difference':fmt(max(frac)),
                     'printed_significant_digits':65})
drift=[r for r in rows if r['K']>=16]
ks=np.array([r['K'] for r in drift],dtype=float)
qs=np.array([float(r['Q_half_x0p5']) for r in drift])
slope,intercept=np.polyfit(np.sqrt(ks),qs,1)
newton=list(csv.DictReader((HERE/'data/exact_scalar.csv').open(encoding='utf-8')))
scalar_errors=[]
for r in rows:
    matched=next(x for x in newton if int(x['K'])==r['K'])
    scalar_errors.append(abs(mp.mpf(r['Q_half_x3'])-mp.mpf(matched['newton_partial'])))
baselines=[]
for f in sorted((HERE/'data').glob('baseline_theta_*.json')):
    r=json.loads(f.read_text(encoding='utf-8'))
    if r['status']=='completed':r['source']=f.name;baselines.append(r)
fit=2-mp.mpf('2.7849974381025855')/mp.mpf('3.5')**mp.mpf('1.467173061332697')
mid=next(r for r in baselines if r['theta']=='.5' and r['requested_dps']==45)
midlow=next(r for r in baselines if r['theta']=='.5' and r['requested_dps']==25)
epsilon=[]
for K in [4,6,8]:
    r=next(r for r in raw if r['K']==K and r['epsilon_checks'])
    checks=r['epsilon_checks'];ratios=[mp.mpf(checks[i]['error'])/mp.mpf(checks[i+1]['error']) for i in range(len(checks)-1)]
    epsilon.append({'K':K,'error_ratios_on_halving':[fmt(t) for t in ratios],
                    'last_error_over_eps_cubed':checks[-1]['error_over_epsilon_cubed']})
out={'meaning':'finite interpolation jet, no established H_3.5',
     'exact_Q_half_x3':fmt(exact),'exact_generator_quadratic_x3':fmt(2*mp.pi**2-13),
     'max_scalar_vs_matrix_absolute_difference':fmt(max(scalar_errors)),
     'precision_controls':controls,'epsilon_checks':epsilon,
     'drift_sqrt_K_fit':{'K_min':16,'K_max':max(ks),'intercept':intercept,'slope':slope,
                        'maximum_residual':float(max(abs(qs-(intercept+slope*np.sqrt(ks))))),
                        'status':'descriptive, not a divergence proof'},
     'eta_fit_at_3p5':fmt(fit),'eta_convex_half':mid['eta'],'fold_convex_half':mid['fold_height'],
     'convex_minus_fit':fmt(mp.mpf(mid['eta'])-fit),
     'relative_difference_to_fit':fmt((mp.mpf(mid['eta'])-fit)/fit),
     'baseline_eta_precision_difference':fmt(abs(mp.mpf(mid['eta'])-mp.mpf(midlow['eta']))),
     'baseline_height_precision_difference':fmt(abs(mp.mpf(mid['fold_height'])-mp.mpf(midlow['fold_height']))),
     'longest_run_compute_seconds':max(run_times+[r['seconds'] for r in baselines]),
     'baseline_roots':[{'theta':r['theta'],'dps':r['requested_dps'],'eta':r['eta'],'height':r['fold_height']} for r in baselines]}
(HERE/'data/analysis.json').write_text(json.dumps(out,indent=2),encoding='utf-8')

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
valid=[r for r in rows if r['K']>=4]
ax[0].loglog([r['K'] for r in valid],[float(r['x3_error']) for r in valid],'o-',label='Fehler bei x=3')
ax[0].loglog([r['K'] for r in valid],[float(r['x3_tail_bound']) for r in valid],'--',label='bewiesene obere Schranke')
ax[0].set(xlabel='Modellgrösse K',ylabel='Absoluter Fehler',title='Ein kontrollierter quadratischer Koeffizient');ax[0].legend();ax[0].grid(alpha=.2)
ax[1].plot(np.sqrt(ks),qs,'o-',label='interpolierter Koeffizient bei x=0.5')
ax[1].plot(np.sqrt(ks),intercept+slope*np.sqrt(ks),'--',label='beschreibender Fit a + c√K')
ax[1].axhline(0,color='grey',lw=.7);ax[1].set(xlabel='√K',ylabel='Quadratischer Koeffizient',title='Drift bei nichtganzzahliger Höhe');ax[1].legend();ax[1].grid(alpha=.2)
fig.suptitle('Halber Rang bis zur zweiten Ordnung — endliches Interpolationsmodell',fontsize=13)
for ext in ['png','svg']:fig.savefig(HERE/f'figures/quadratischer_halbschritt.{ext}',dpi=180)
plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
unique={float(r['theta']):r for r in sorted(baselines,key=lambda r:r['requested_dps'])}
points=[(0,1.4446678610097661),*[(t,float(r['eta'])) for t,r in sorted(unique.items())],(1,1.6353244967152764)]
ax.plot([p[0] for p in points],[p[1] for p in points],'o-',label='Tangenz der konvexen Interpolation')
ax.scatter([.5],[float(fit)],marker='x',s=85,color='#c64924',label='alter η-Fit bei Rang 3.5')
ax.set(xlabel='Mischgewicht θ',ylabel='Tangenzbasis',title='Vergleichskurve (1−θ)bˣ + θT_b(x) — kein nachgewiesener Halb-Rang')
ax.grid(alpha=.2);ax.legend()
for ext in ['png','svg']:fig.savefig(HERE/f'figures/eta_vergleich.{ext}',dpi=180)
plt.close(fig)
print(json.dumps(out,indent=2),flush=True)
