"""Independent scalar Newton-series audit for the rigorously solved x=3 jet."""
import json,csv,time
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).absolute().parent
mp.mp.dps=170
def F(s):
    s=mp.mpf(s)
    return 2*(s+2)*(s+3)*(mp.digamma(s+3)+mp.euler)-(s+2)*(3*s+7)-4
def coefficient(p):
    if p==0:return mp.mpf(0)
    if p==1:return mp.mpf(10)
    if p==2:return mp.mpf(22)/3
    return mp.mpf(24*(-1)**(p+1))/((p-2)*(p-1)*p*(p+1))
def partial(K,s):return mp.fsum(mp.binomial(s,p)*coefficient(p) for p in range(K))
def bound(K):return 8/(3*mp.sqrt(mp.pi))*(K-3)**(-mp.mpf('4.5'))
def fmt(x):return mp.nstr(x,80)
start=time.monotonic(); rows=[]
for K in [3,4,6,8,10,12,16,20,24,32,48,64,80,128,256,512,1024]:
    q=partial(K,'.5');error=F('.5')-q
    rows.append({'K':K,'newton_partial':fmt(q),'exact_limit':fmt(F('.5')),'positive_tail':fmt(error),
                 'tail_bound':fmt(bound(K)) if K>=4 else '',
                 'rescaled_tail_K_4p5':fmt(error*K**mp.mpf('4.5'))})
with (HERE/'data/exact_scalar.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
checks=[]
for p in range(1,15):
    difference=mp.fsum((-1)**(p-j)*mp.binomial(p,j)*F(j) for j in range(p+1))
    checks.append({'p':p,'forward_difference':fmt(difference),'formula':fmt(coefficient(p)),
                   'error':fmt(abs(difference-coefficient(p)))})
out={'exact_Q_half_x3':'341/12 - 35 log(2)','value':fmt(F('.5')),
     'closed_formula_error':fmt(abs(F('.5')-(mp.mpf(341)/12-35*mp.log(2)))),
     'newton_checks':checks,'asymptotic_tail_constant':fmt(8/(3*mp.sqrt(mp.pi))),
     'seconds':time.monotonic()-start}
(HERE/'data/exact_scalar_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(out['value'],out['seconds'],flush=True)
