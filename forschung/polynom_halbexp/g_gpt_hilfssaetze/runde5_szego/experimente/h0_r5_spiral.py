"""Compare a fixed-angle corner coordinate with the exact complex-cone scale.

No claim of all-n convergence. Coefficients were shifted before being saved.
No evaluation of original A coefficients in the bay with truncated precision.
"""
import json,sys,math
from pathlib import Path
from h0_r5_paths import locations
import numpy as np
from flint import arb,acb,acb_poly,ctx
sys.stdout.reconfigure(encoding='utf-8');sys.dont_write_bytecode=True
folder,ROUND4,ROOT=locations(__file__)
sys.path.insert(0,str(ROOT/'o_operator'));import o_kern as K
ctx.prec=256
mu=K.Lacb().sqrt();kappa=acb(0,2*arb.pi())/(acb(0,2*arb.pi())-mu.log());inv=1/kappa
rho=abs(K.Lacb());L=acb(K.Lc.real,K.Lc.imag);wL=rho/L
def scalar(v):return acb(arb(v[0]),arb(v[1]))
def vector(values):return np.array([complex(float(v.real),float(v.imag)) for v in values])
def comp(v):return [v.real.mid().str(40,radius=False),v.imag.mid().str(40,radius=False)]
out={'kappa':comp(kappa),'inverse_kappa':comp(inv),'mu':comp(mu),
     'status':'finite numerical profiles in the complex-cone coordinate; no global parametrix theorem','models':{}}
for model in ['kubisch','kneser_z']:
    by={}
    for file in folder.glob(f'profiles_{model}_*.json'):
        for c in json.loads(file.read_text())['cases']:by[c['n']]=c
    if not by:continue
    shifted={};original_shape={};records=[]
    for n,c in sorted(by.items()):
        coefficients=list(map(scalar,c['w_circle_taylor_scaled']));P=acb_poly(coefficients)
        center=scalar(c['w_center']);R=arb(c['w_radius']);dn=acb(n)**(-inv)
        T=P(acb_poly([(wL-center)/R,dn/R]));corner=scalar(c['corner_scaled'])
        assert abs(T[0]-corner)<arb('1e-35')
        shape=T*(1/corner);shifted[n]=vector([shape[k] for k in range(n+1)])
        original_shape[n]=vector(coefficients)/complex(float(corner.real),float(corner.imag))
        records.append({'n':n,'natural_shape_coefficients_first_40':[comp(shape[k]) for k in range(40)],
                        'natural_majorant_shape':float(sum((abs(shape[k]) for k in range(n+1)),arb(0))),
                        'natural_majorant_shape_R105':float(sum((abs(shape[k])*arb('1.05')**k for k in range(n+1)),arb(0))),
                        'shape_tail_after_40_R105':float(sum((abs(shape[k])*arb('1.05')**k for k in range(40,n+1)),arb(0))),
                        'coordinate_derivative':comp(dn),'corner_abs_scaled':float(abs(corner))})
    diffs=[]
    for n,m in [(200,400),(400,800),(500,1000),(300,600),(600,900),(900,1000)]:
        if n not in by or m not in by:continue
        a=shifted[n];b=shifted[m];a=np.pad(a,(0,len(b)-len(a)))
        aa=original_shape[n];bb=original_shape[m];aa=np.pad(aa,(0,len(bb)-len(aa)))
        diffs.append({'n':n,'m':m,'natural_normalized_W_difference':float(sum(abs(a-b))),
                      'fixed_angle_normalized_W_difference':float(sum(abs(aa-bb)))})
    fits=[]
    fit_ns=[n for n in sorted(by) if n in [200,300,400,500,600,700,800,900,1000]]
    for low in [200,400,600]:
        ns=[n for n in fit_ns if n>=low]
        if len(ns)<4:continue
        X=np.array([[1,1/n,1/n**2] for n in ns])
        components=[]
        for k in range(5):
            cf,*_=np.linalg.lstsq(X,np.array([shifted[n][k] for n in ns]),rcond=None)
            components.append([cf[0].real,cf[0].imag])
        amp=np.array([float(abs(scalar(by[n]['corner_scaled']))) for n in ns]);cf,*_=np.linalg.lstsq(X,amp,rcond=None)
        fits.append({'low':low,'high':max(ns),'shape_limit_first_5':components,'amplitude_limit_1_over_n_fit':cf[0],
                     'max_amplitude_residual':float(max(abs(X@cf-amp)))})
    heldout=[]
    train=[n for n in [200,300,400,500,600] if n in by]
    test=[n for n in [700,800,900,1000] if n in by]
    if len(train)==5 and test:
        for powers in [[0,1,2],[0,float(inv.real),1],[0,float(inv.real),1,2*float(inv.real)]]:
            X=np.array([[n**(-p) for p in powers] for n in train]);Y=np.array([shifted[n][:40] for n in train])
            co,*_=np.linalg.lstsq(X,Y,rcond=None)
            errs={n:float(sum(abs(np.array([n**(-p) for p in powers])@co-shifted[n][:40])*1.05**np.arange(40))) for n in test}
            heldout.append({'train':train,'test':test,'powers':powers,'weighted_W_errors_R105':errs,
                            'extrapolated_shape_coefficients_40':[[float(c.real),float(c.imag)] for c in co[0]],
                            'extrapolated_shape_weighted_norm_R105':float(sum(abs(co[0])*1.05**np.arange(40))),
                            'status':'finite numerical extrapolated polynomial; no limiting tail certificate'})
    out['models'][model]={'profiles':records,'differences':diffs,'finite_extrapolations':fits,'heldout_shape_fits':heldout}
    print(model,'differences',diffs,'fits',fits,flush=True)
(folder/'spiral_analysis.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
