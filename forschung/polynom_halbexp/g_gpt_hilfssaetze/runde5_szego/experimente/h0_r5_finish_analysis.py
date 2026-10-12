"""Finalize the numerical stage-1 holdout and precision diagnostics."""
import json,sys,math
from pathlib import Path
from h0_r5_paths import locations
import mpmath as mp
mp.mp.dps=75;sys.stdout.reconfigure(encoding='utf-8')
folder,ROUND4,ROOT=locations(__file__)
analysis=json.loads((folder/'family_analysis.json').read_text());by={c['name']:c for c in analysis['families']};checks=[]
for file in (folder/'holdout_600').glob('family_*.json'):
    d=json.loads(file.read_text());fit=next(c for c in by[d['name']]['fits'] if c['low']==250 and c['powers']==[0,2,3,4])
    co=list(map(mp.mpf,fit['coefficients']));powers=fit['powers']
    errors=[abs(mp.log(mp.mpf(c['ratio']))-sum(co[j]*mp.mpf(c['n'])**(-p) for j,p in enumerate(powers))) for c in d['cases'] if c['n']>400]
    basefile=folder/file.name;old=json.loads(basefile.read_text());overlap=[]
    for a,b in zip(old['cases'],d['cases']):
        assert a['n']==b['n'];overlap.append(abs(mp.mpf(a['ratio'])/mp.mpf(b['ratio'])-1))
    assert max(overlap)<mp.mpf('1e-60')
    checks.append({'model':d['name'],'trained':[250,400],'heldout':[401,600],
                   'max_log_ratio_prediction_error':str(max(errors)),
                   'maximum_overlap_relative_ratio_difference':str(max(overlap)),'seconds':d['seconds']})
phase=json.loads((folder/'phase_analysis.json').read_text());spiral=json.loads((folder/'spiral_analysis.json').read_text())
prediction=2*float(spiral['inverse_kappa'][1])
out={'status':'finite numerical controls only','holdout_families':checks,
     'local_phase_prediction_2Im_inverse_kappa':prediction,
     'corner_phase_fits':phase['corner_higher_order_fits'],
     'phase_candidate_status':'nu=2 Im(1/kappa) is conjectural for the Schur amplitude phase; spatial twist follows algebraically from kappa'}
(folder/'final_controls.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
for name,d in spiral['models'].items():
 print(name,'natural shape norm at last degree',d['profiles'][-1]['natural_majorant_shape_R105'])
 if d['heldout_shape_fits']:print('shape holdout:',d['heldout_shape_fits'])
