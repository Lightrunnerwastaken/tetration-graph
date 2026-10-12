"""High precision phase fits from adjacent late degrees; explicit carrier.

Phase fits are numerical, not a proof of a global conformal g-function.
"""
import json,sys
from pathlib import Path
from h0_r5_paths import locations
import mpmath as mp
mp.mp.dps=65;sys.stdout.reconfigure(encoding='utf-8')
folder,ROUND4,ROOT=locations(__file__);by={}
for file in folder.glob('profiles_kubisch_*.json'):
    for c in json.loads(file.read_text())['cases']:by[c['n']]=c
def cv(v):return mp.mpc(v[0],v[1])
ns=[n for n in range(950,960) if n in by and n+1 in by]
cols=lambda n:[mp.mpf(1),mp.log(mp.mpf(n+1)/n),mp.mpf(1)/(n+1)-mp.mpf(1)/n,mp.mpf(1)/(n+1)**2-mp.mpf(1)/n**2]
A=mp.matrix([cols(n) for n in ns]);freqs={};anchor=1000
for point in ['corner','lens_2','lens_4','bay_2','bay_4']:
    b=mp.matrix([mp.arg(cv(by[n+1]['fixed_scaled'][point])/cv(by[n]['fixed_scaled'][point])) for n in ns])
    cf,_=mp.qr_solve(A,b)
    checks=[];base=cv(by[anchor]['fixed_scaled'][point]);base/=abs(base)
    for n in [200,300,400,500,600,700,800,900,950,960,970,980,990,1000]:
        if n not in by:continue
        actual=cv(by[n]['fixed_scaled'][point]);actual/=abs(actual)
        delta=cf[0]*(n-anchor)+cf[1]*mp.log(mp.mpf(n)/anchor)
        main=base*mp.exp(1j*delta)
        deltafull=delta+cf[2]*(mp.mpf(1)/n-mp.mpf(1)/anchor)+cf[3]*(mp.mpf(1)/n**2-mp.mpf(1)/anchor**2)
        full=base*mp.exp(1j*deltafull)
        checks.append({'n':n,'phase_error_leading':str(abs(main-actual)),'phase_error_corrected':str(abs(full-actual))})
    freqs[point]={'omega':str(cf[0]),'nu_log_n':str(cf[1]),'c_inverse_n':str(cf[2]),'c_inverse_n_squared':str(cf[3]),
                   'phase_anchor_degree':anchor,'phase_anchor':[str(base.real),str(base.imag)],'checks':checks}
    print(point,'omega=',str(cf[0])[:28],'nu=',str(cf[1])[:28])
    if point=='corner':print('checks:',checks)
out={'status':'numerical exact-basis phase fits; no all-n phase theorem','training_degrees':[950,960],
     'phase_formula':'chi_n=chi_1000 exp(i[omega(n-1000)+nu log(n/1000)+c1(1/n-1/1000)+c2(1/n^2-1/1000^2)])',
     'points':freqs}
additional=[]
for order in [1,3,4]:
    rows=[[mp.mpf(1),mp.log(mp.mpf(n+1)/n)]+[mp.mpf(n+1)**(-k)-mp.mpf(n)**(-k) for k in range(1,order+1)] for n in ns]
    b=mp.matrix([mp.arg(cv(by[n+1]['phase_corner'])/cv(by[n]['phase_corner'])) for n in ns])
    co,_=mp.qr_solve(mp.matrix(rows),b)
    additional.append({'inverse_order':order,'omega':str(co[0]),'nu_log_n':str(co[1])})
out['corner_higher_order_fits']=additional
broad=[];omega0=mp.mpf(freqs['corner']['omega']);phase0=cv(by[anchor]['phase_corner'])
for low in [200,400,600]:
    regular=[n for n in [200,300,400,500,600,700,800,900,1000] if n>=low and n!=anchor]
    for order in [1,2,3]:
        if len(regular)<order+2:continue
        rows=[];values=[]
        for n in regular:
            values.append(mp.arg(cv(by[n]['phase_corner'])/phase0*mp.exp(-1j*omega0*(n-anchor))))
            rows.append([mp.mpf(n-anchor),mp.log(mp.mpf(n)/anchor)]+[mp.mpf(n)**(-k)-mp.mpf(anchor)**(-k) for k in range(1,order+1)])
        co,_=mp.qr_solve(mp.matrix(rows),mp.matrix(values))
        broad.append({'low':low,'inverse_order':order,'omega':str(omega0+co[0]),'nu_log_n':str(co[1])})
out['corner_broad_window_fits']=broad
(folder/'phase_analysis.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
