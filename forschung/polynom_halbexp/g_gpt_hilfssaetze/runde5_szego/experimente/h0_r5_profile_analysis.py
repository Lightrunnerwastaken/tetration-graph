"""Complex-profile checks, including phase-invariant lower bounds."""
import json,sys,math
from pathlib import Path
from h0_r5_paths import locations
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
folder,ROUND4,ROOT=locations(__file__)
def cv(v):return float(v[0])+1j*float(v[1])
out={}
for model in ['kubisch','kneser_z']:
    by={}
    for file in folder.glob(f'profiles_{model}_*.json'):
        d=json.loads(file.read_text(encoding='utf-8'))
        for c in d['cases']:by[c['n']]=c
    if not by:continue
    ns=sorted(by);comparisons=[]
    for n,m in [(n,m) for n,m in zip(ns,ns[1:])]+[(200,400),(400,800),(300,600),(500,1000)]:
        if n not in by or m not in by:continue
        a=np.array(list(map(cv,by[n]['w_circle_taylor_scaled'])));b=np.array(list(map(cv,by[m]['w_circle_taylor_scaled'])))
        a=np.pad(a,(0,len(b)-len(a)))
        chi_n=cv(by[n]['phase_corner']);chi_m=cv(by[m]['phase_corner'])
        ca=a/chi_n;cb=b/chi_m
        lower=math.sqrt(max(0,float(np.vdot(a,a).real+np.vdot(b,b).real-2*abs(np.vdot(a,b)))))
        phase=np.exp(1j*np.angle(np.vdot(a,b)))
        x=np.array(list(map(cv,by[n]['x_circle_scaled'])))/chi_n
        y=np.array(list(map(cv,by[m]['x_circle_scaled'])))/chi_m
        fixed_n={k:cv(v)/chi_n for k,v in by[n]['fixed_scaled'].items()}
        fixed_m={k:cv(v)/chi_m for k,v in by[m]['fixed_scaled'].items()}
        fixed_d={k:abs(fixed_n[k]-fixed_m[k]) for k in fixed_n}
        comparisons.append({'n':n,'m':m,'W_difference_corner_phase':float(sum(abs(ca-cb))),
                    'W_difference_l2_optimal_phase':float(sum(abs(phase*a-b))),
                    'phase_invariant_l2_lower_bound':lower,'x_circle_sup_difference_corner_phase':float(max(abs(x-y))),
                    'fixed_sup_difference_corner_phase':max(fixed_d.values()),
                    'fixed_argmax_difference':max(fixed_d,key=fixed_d.get),
                    'minimum_possible_fixed_sup_difference_magnitudes':max(abs(abs(fixed_n[k])-abs(fixed_m[k])) for k in fixed_n)})
    # Nearby degrees determine the fast carrier without aliasing by 100-degree samples.
    adjacent=[]
    for n in ns:
        if n+1 in by:
            quotient=cv(by[n+1]['phase_corner'])/cv(by[n]['phase_corner'])
            adjacent.append({'n':n,'phase_increment':float(np.angle(quotient))})
    mag_corner={n:abs(cv(by[n]['corner_scaled'])) for n in ns}
    coeff_table={n:[cv(v) for v in by[n]['w_circle_taylor_scaled'][:7]] for n in ns}
    out[model]={'degrees':ns,'comparisons':comparisons,'adjacent_phase_increments':adjacent,
                'corner_scaled_absolute':mag_corner,
                'phase_corner_definition':'chi_n=D_n(L)/abs(D_n(L)); data-based explicit phase, not a proved analytic phase',
                'status':'finite numerical comparisons; no inferred uniform all-n error bound'}
    print(model,'corner:',mag_corner)
    for c in comparisons:
        if (c['n'],c['m']) in [(200,400),(400,800),(300,600),(500,1000)]:print(c)
(folder/'profile_analysis.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
