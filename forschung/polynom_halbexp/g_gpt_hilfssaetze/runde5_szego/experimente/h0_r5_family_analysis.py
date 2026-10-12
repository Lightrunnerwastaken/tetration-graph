"""Fits for stage 1. All extrapolations stay labelled numerical."""
import json,sys
from pathlib import Path
from h0_r5_paths import locations
import mpmath as mp
sys.stdout.reconfigure(encoding='utf-8');mp.mp.dps=75
folder,ROUND4,ROOT=locations(__file__);result=[]
for file in sorted(folder.glob('family_*.json')):
    if file.name=='family_analysis.json':continue
    d=json.loads(file.read_text(encoding='utf-8'));cases=d['cases'];N=d['N'];fits=[]
    for low in [100,150,200,250]:
        selected=[c for c in cases if c['n']>=low]
        for powers in [[0,2,3,4],[0,1,2,3,4]]:
            A=mp.matrix([[mp.mpf(c['n'])**(-k) for k in powers] for c in selected]);b=mp.matrix([mp.log(mp.mpf(c['ratio'])) for c in selected]);co,_=mp.qr_solve(A,b)
            err=max(abs(sum(co[j]*mp.mpf(c['n'])**(-k) for j,k in enumerate(powers))-mp.log(mp.mpf(c['ratio']))) for c in selected)
            fits.append({'low':low,'high':N,'powers':powers,'lambda':str(mp.exp(co[0])),
                         'coefficients':[str(v) for v in co],'max_log_residual':str(err)})
    item={'name':d['name'],'N':N,'fits':fits,'candidates':d['candidates'],'description':d['description'],'seconds':d['seconds']}
    result.append(item);f=fits[-2]
    print(d['name'],'lambda=',f['lambda'][:22],'c2=',f['coefficients'][1][:22],'res=',f['max_log_residual'][:12])
checks=[]
for file in (folder/'precision_checks').glob('family_*.json'):
    high=json.loads(file.read_text());low=json.loads((folder/file.name).read_text())
    maxerr=max(abs(mp.mpf(a['ratio'])/mp.mpf(b['ratio'])-1) for a,b in zip(high['cases'],low['cases']))
    assert maxerr<mp.mpf('1e-60')
    checks.append({'model':high['name'],'extra_bits':high['extra_bits'],'max_relative_ratio_difference':str(maxerr)})
out={'status':'numerical extrapolation, not certified asymptotic constants','families':result,'precision_checks':checks}
(folder/'family_analysis.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
