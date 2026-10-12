"""Round 5, stage 1: real polynomial families, native every-degree pivots.

Numerical midpoint diagnostic, 1200+n bits. Each model is saved immediately.
No SVD, no Kneser input for construction, no nonlinear overdetermined systems.
"""
import argparse,json,math,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from flint import arb,acb,arb_poly,arb_mat,ctx
from h0_schur_native import Border

def build(F,tau,N):
    M=F.derivative()(F)*tau
    pw=arb_poly([1]);J=arb_mat(N+1,N+1)
    for j in range(N+1):
        for k in range(N+1):
            J[k,j]=pw[k]+(M[k-j] if k>=j else arb(0))
        pw=(pw*F).truncate(N+1)
    return J,M

def candidate(F,l):
    import numpy as np
    ff=np.array([float(v) for v in F]);z=np.exp(2j*np.pi*(np.arange(131072)+.5)/131072)
    image=np.polynomial.polynomial.polyval(z,ff)
    avg=float(np.maximum(0,np.log(abs(image))).mean())
    a=ff[0];R=abs(l-a)
    lo=2*np.angle(a+R-l);hi=np.angle(l)
    p=np.pi/(hi-lo)
    cap=p*l.imag/np.sin(-p*lo)
    return {'F1':ff[1],'folded_circle_mean_log':avg,
            'exp_mean_log':float(np.exp(avg)),'exp_twice_mean_log':float(np.exp(2*avg)),
            'W_capacity_normalized':float(cap),'W_capacity_squared':float(cap*cap)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--max-n',type=int,default=400);p.add_argument('--models',required=True)
    p.add_argument('--extra-bits',type=int,default=0);p.add_argument('--deadline',type=float,default=900)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);start=time.time()
    sys.path.insert(0,str(a.project_root/'o_operator'));import o_kern as K
    ctx.prec=1300+a.max_n+a.extra_bits
    L=K.Lacb();rho=abs(L);l=L/rho;P=arb_poly([1,-2*l.real,1])
    base=arb_poly([c*rho**(j-1) for j,c in enumerate(K.kubik_koeff())])
    previous_seconds=None
    for name in a.models.split(','):
        if previous_seconds is not None and time.time()-start+1.4*previous_seconds>a.deadline:
            print('Session deadline: start remaining models in a new session.',flush=True);break
        ts=time.time();tau=arb(1);desc={}
        if name=='base':F=base
        elif name.startswith('jet_'):
            epsilon=arb(name[4:]);F=base+epsilon*P*P;desc={'epsilon':str(epsilon),'fixed_L_mu':True,'degree':4}
        elif name.startswith('mu_'):
            change=arb(name[3:]);mu=L.sqrt()*(1+change)
            F=arb_poly([c*rho**(j-1) for j,c in enumerate(K.kubik_koeff(mu))]);desc={'mu_relative_change':str(change),'fixed_L_mu':False,'degree':3}
        elif name.startswith('quad_'):
            kappa=arb(name[5:]);F=arb_poly([0,1])+kappa*P;desc={'kappa':str(kappa),'degree':2}
        elif name.startswith('tau_'):
            tau=arb(name[4:]);F=base;desc={'M_amplitude':str(tau)}
        else:raise ValueError(name)
        J,M=build(F,tau,a.max_n);B=Border(J,ctx.prec);del J
        result={'status':'numerical midpoint diagnostic, not an all-n theorem','name':name,
                'N':a.max_n,'description':desc,'F_coefficients':[v.str(60,radius=False) for v in F],
                'M_coefficients':[v.str(60,radius=False) for v in M],
                'F_at_l':[str(v) for v in [F(l).real,F(l).imag]],
                'derivative_at_l':[str(v) for v in [F.derivative()(l).real,F.derivative()(l).imag]],
                'candidates':candidate(F,complex(float(l.real),float(l.imag))),
                'bits':'1200+n+extra_bits','extra_bits':a.extra_bits,'cases':[]}
        previous=arb(B.api.string(B.api.entry(B.J,0,0),80)).mid()
        for n in range(1,a.max_n+1):
            ctx.prec=1200+n+a.extra_bits;B.step(n,ctx.prec);api=B.api
            h=arb(api.string(B.h,90)).mid();q=arb(api.string(B.qnorm,90)).mid()
            if n>=25:
                row=[arb(api.string(api.entry(B.inv,n,j),int(ctx.prec*.30103)+12)).mid() for j in range(n+1)]
                an=abs(h)*sum((v*v for v in row),arb(0)).sqrt()
                value=q*abs(arb_poly(row)(1/l))
                result['cases'].append({'n':n,'h':h.str(70,radius=False),'ratio':(h/previous).str(65,radius=False),
                    'Q_norm':q.str(65,radius=False),'A_norm':an.str(65,radius=False),
                    'D_corner_abs':float(value),'sqrt_n_D_corner_abs':float(value)*math.sqrt(n)})
            previous=h
        B.close();result['seconds']=time.time()-ts;previous_seconds=result['seconds']
        out=a.output_dir/('family_'+name.replace('-','m').replace('.','p')+'.json')
        out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
        print(f'{name}: ratio({a.max_n})={result["cases"][-1]["ratio"][:23]}; {result["seconds"]:.1f}s',flush=True)
        if time.time()-start>a.deadline:
            print('Session deadline: remaining models must be launched in a new session.',flush=True);break

if __name__=='__main__':main()
