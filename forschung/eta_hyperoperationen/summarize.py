"""Tabulate experiments, precision checks, and runtimes without re-running them."""
import csv,json,math
from pathlib import Path
from common import HERE,mp
mp.mp.dps=70
def read(name):return json.loads((HERE/name).read_text(encoding='utf-8'))
def probe(name):
    dat=read(name)
    return dat,dict(zip(dat['log'][0]['expressions'],dat['log'][0]['values']))
def value(v):return mp.mpc(*map(mp.mpf,v))
low=read('eta7_25.json');high=read('eta7_45.json');controls={}
controls['eta7']={k:{'absolute_difference':str(abs(mp.mpf(low[k])-mp.mpf(high[k]))),
    'decimal_agreement_indicator':float(-mp.log10(abs(mp.mpf(low[k])-mp.mpf(high[k]))))} for k in ['eta','euler']}
p25,v25=probe('probe_hyper_r7_2_25.json');p45,v45=probe('probe_hyper_r7_2_45.json')
controls['H7_base2']={k:str(abs(value(v25[k])-value(v45[k]))) for k in ['hrp[7]','hrlam[7]','hval(7,.5)','hval(7,1.5)']}
if (HERE/'eta8_control_retry.json').exists() and read('eta8_control_retry.json').get('completed'):
    l=read('eta8_refined.json');h=read('eta8_control_retry.json')
    controls['eta8']={k:{'absolute_difference':str(abs(mp.mpf(l[k])-mp.mpf(h[k]))),
         'decimal_agreement_indicator':float(-mp.log10(abs(mp.mpf(l[k])-mp.mpf(h[k]))))} for k in ['eta','euler']}
    controls['eta8'].update(high_root_residual=h['minimum'],low_root_residual=l['minimum'],
                           caution='different root stopping tolerances; agreement is not an absolute error bound')
(HERE/'precision_controls.json').write_text(json.dumps(controls,indent=2),encoding='utf-8')
audit=[]
for f in sorted(HERE.glob('*.json')):
    if f.name in ['runtime_audit.json']:continue
    dat=json.loads(f.read_text(encoding='utf-8'))
    if isinstance(dat,dict) and 'seconds' in dat:
        audit.append({'file':f.name,'seconds':dat['seconds'],'status':dat.get('status'),
                      'completed':dat.get('completed'),'below_20_minutes':dat['seconds']<1200})
(HERE/'runtime_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
analysis=read('analysis.json');names={3:'Exponentiation',4:'Tetration',5:'Pentation',6:'Hexation',7:'Heptation',8:'Octation',9:'Nonation',10:'Decation'}
with (HERE/'eta_values.csv').open('w',newline='',encoding='utf-8-sig') as file:
    writer=csv.writer(file);writer.writerow(['map_rank','map','next_rank','next_operation','eta','e','status'])
    for row in analysis['values']:
        r=row['rank'];writer.writerow([r,names[r],r+1,names[r+1],row['eta'],row['e'],'exact' if r==3 else 'numerical selected local tangency'])
lower=read('probe_lower_2.json');vv=dict(zip(lower['expressions'],lower['values']))
halves=[{'rank':r,'half':vv[f'hval({r},.5)'][0],'one_and_half':vv[f'hval({r},1.5)'][0]} for r in range(3,7)]
for r,file in [(7,'probe_hyper_r7_2_45.json'),(8,'probe_hyper_r8_2_25_pilot.json'),(9,'probe_hyper_r9_2_45_retry.json')]:
    dat,v=probe(file);halves.append({'rank':r,'half':v[f'hval({r},.5)'][0],
                                   'one_and_half':v[f'hval({r},1.5)'][0],
                                   'source':file,'estimated_digits':v.get(f'hrerror[{r}]',[None])[0]})
(HERE/'half_height_values.json').write_text(json.dumps(halves,indent=2),encoding='utf-8')
if (HERE/'real_decation.json').exists():
    real=read('real_decation.json')
    if len(real.get('samples',[]))==3:
        q=mp.mpf(real['attractor_multiplier_checks'][2][0])**2
        estimates=[[str((mp.mpf(real['samples'][i+1]['values'][j][0])-q*mp.mpf(real['samples'][i]['values'][j][0]))/(1-q)) for j in [0,1]] for i in [0,1]]
        diff=[str(abs(mp.mpf(estimates[1][j])-mp.mpf(estimates[0][j]))) for j in [0,1]]
        (HERE/'decation_extrapolation.json').write_text(json.dumps({'lambda_squared':str(q),
            'heights':['0.5','1.5'],'estimates_from_k5_k7':estimates[0],
            'estimates_from_k7_k9':estimates[1],'absolute_differences':diff,
            'status':'observed real convergence; no uniform tail or interval bound'},indent=2),encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
rs=[z['rank'] for z in halves]
for axis,key,line,title in zip(ax,['half','one_and_half'],[1.5,2.5],['Hᵣ(2, 0.5)','Hᵣ(2, 1.5)']):
    axis.plot(rs,[float(z[key]) for z in halves],'o-',color='#235587')
    axis.axhline(line,color='#c97930',ls='--');axis.set(title=title,xlabel='Rang r',xticks=rs);axis.grid(alpha=.2)
fig.suptitle('Zwischenwerte bei Basis 2\nGestrichelte Linien: Bezugswerte, keine bewiesenen Grenzwerte',fontsize=12)
fig.savefig(HERE/'zwischenwerte.png',dpi=180);plt.close(fig)
print(json.dumps({'precision':controls,'max_recorded_session':max(z['seconds'] for z in audit),'half_heights':halves},indent=2))
