"""Independent direct +300-bit checks of new complex profile evaluations."""
import json,sys,time
from pathlib import Path
from h0_r5_paths import locations
sys.dont_write_bytecode=True;sys.stdout.reconfigure(encoding='utf-8')
from flint import arb,acb,arb_mat,acb_poly,ctx
folder,ROUND4,ROOT=locations(__file__)
root=ROOT
sys.path.insert(0,str(root/'o_operator'));import o_kern as K
start=time.time();cases=[]
for model in ['kubisch','kneser_z']:
    n=400;ctx.prec=1500+n;J,_,_=K.matrizen(model,n,ctx.prec)
    e=arb_mat(n+1,1);e[n,0]=1
    c=J.solve(e,algorithm='approx');a=J.transpose().solve(e,algorithm='approx')
    q=(sum((c[k,0].mid()**2 for k in range(n+1)),arb(0)).sqrt()/abs(c[n,0])).mid()
    pol=acb_poly([a[k,0].mid()*q*arb(n).sqrt() for k in range(n+1)])
    refs=[]
    for file in folder.glob(f'profiles_{model}_*.json'):
        d=json.loads(file.read_text());refs.extend([(d,z) for z in d['cases'] if z['n']==n])
    assert len(refs)==1
    d,ref=refs[0];rho=abs(K.Lacb()).mid();err=arb(0)
    for point,x in d['fixed_points'].items():
        old=acb(arb(ref['fixed_scaled'][point][0]),arb(ref['fixed_scaled'][point][1]))
        now=pol(rho/acb(*x));err=max(err,abs(now-old).mid()/(1+abs(now).mid()))
    assert err<arb('1e-40')
    cases.append({'model':model,'n':n,'extra_bits':300,'maximum_relative_complex_fixed_error':err.str(50,radius=False)})
    print(cases[-1],flush=True)
out={'status':'numerical independent direct check; baseline stored 45 digits, no interval inverse certificate',
     'cases':cases,'seconds':time.time()-start}
(folder/'profile_direct_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
