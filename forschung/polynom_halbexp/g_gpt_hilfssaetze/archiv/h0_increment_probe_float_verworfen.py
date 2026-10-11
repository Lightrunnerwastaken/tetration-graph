import json
from pathlib import Path
import numpy as np

root = Path(r'D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp')
ref = json.loads((root/'g_gpt_hilfssaetze/ergebnisse/h0_cauchy_kubisch.json').read_text(encoding='utf-8'))
L = complex(.31813150520476413,1.3372357014306895)
rho = abs(L)
for N in [200,300,400,500,600,700,800,900,1000]:
    G = np.load(root/f'w_grenzschicht/w01_N{N}.npz')['G']
    q = G[:,-1]/G[-1,-1]
    points = [('L',L)]
    source = next(x for x in ref['cases'] if x['N']==200)
    points += [(name,complex(*x['x'])) for name,x in source['points'].items() if name.startswith(('bay_','lens_'))]
    vals=[]
    for name,x in points:
        v = G @ ((rho/x)**np.arange(N+1))
        inc = np.linalg.norm(q)*abs(v[-1])
        vals.append(inc)
    print(N,'Qnorm',np.linalg.norm(q),'dL',vals[0],'dmax',max(vals),'sqrtN*dmax',np.sqrt(N)*max(vals),'where',points[np.argmax(vals)][0],flush=True)
