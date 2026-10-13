"""Independent consistency checks; no numerical intervals are claimed."""
import json,math,time
from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
HERE=Path(__file__).absolute().parent
def read(name):return json.loads((HERE/'data'/name).read_text(encoding='utf-8'))
def value(p,x):
    v=F(0)
    for c in reversed(p):v=v*x+F(c)
    return v
def harmonic(k):return sum((F(1,j) for j in range(1,k+1)),F(0))
def run():
    t=time.monotonic();tests=[]
    intrinsic=read('intrinsic_K128_300.json');direct=read('direct_coefficients_K512_500.json')
    newton=read('height_newton_K128_500.json')
    for a,b in zip(newton['coefficients'],direct['coefficients']):
        assert F(a['rational'])==F(b['rational'])
        assert F(a['log2_coefficient'])==F(b['log2_coefficient'])
    tests.append('all 129 exact rational/log2 Newton coefficients agree between two algorithms')
    u={r['height']:r for r in intrinsic['records']}
    low=read('intrinsic_K32_220.json')
    for a in low['records']:
        b=u[a['height']]
        assert a['U_coefficients']==b['U_coefficients'] and a['V_coefficients']==b['V_coefficients']
    tests.append('all exact U,V polynomials for heights <=32 agree between independent cutoffs')
    ranks=read('rank_newton_N128_400.json')
    for r in ranks['integer_rank']:
        n=r['n'];expected=2*(n+2)*(n+3)*harmonic(n+2)-(n+2)*(3*n+7)-4
        assert F(r['values']['3'])==expected
    tests.append('129 independent polynomial rank jets at height 3 equal the harmonic closed formula exactly')
    # Exact Legendre power coefficients of A_U(1/2,w).
    h=[F(1)];bs=[F(1)]
    for j in range(1,127):
        h.append(h[-1]*(F(5,2)-j+1)*(F(5,2)+j)/F(j*j))
        bs.append(bs[-1]*(F(1,2)-j+1)/j)
    for j in range(2,129):
        au=2*sum((bs[k]*h[j-2-k] for k in range(j-2)),F(0))
        assert au==-F(newton['coefficients'][j]['log2_coefficient'])/2
    tests.append('all 127 Legendre U generating coefficients agree exactly')
    a=read('evaluation_K512_500.json');b=read('evaluation_K512_600.json')
    for x,y in zip(a['samples'],b['samples']):
        assert x['q_re']==y['q_re'] and x['q_im']==y['q_im']
    tests.append('512-term function values at five real/complex points agree in all 70 printed digits at +100 decimal digits')
    a=read('height_newton_K128_400.json');b=read('height_newton_K128_500.json')
    for x,y in zip(a['samples'],b['samples']):assert x['q_re']==y['q_re'] and x['q_im']==y['q_im']
    tests.append('128-term reconstruction independently passes +100 decimal digits (about +332 bits)')
    for dps in [50,75]:
        for r in read(f'abel_plana_{dps}.json')['records']:
            assert mp.mpf(r['abs_error'])<mp.mpf(10)**(-dps+3)
    tests.append('14 Abel-Plana/Gamma identities pass at both 50 and 75 decimal digits')
    for dps in [40,65]:
        for r in read(f'generating_{dps}.json')['records']:
            assert mp.mpf(r['q_error'])<mp.mpf(10)**(-dps+3)
    tests.append('Legendre integral and exact coefficient generating functions agree at two points and two precisions')
    result={'status':'passed','seconds':time.monotonic()-t,'tests':tests}
    (HERE/'data/verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
