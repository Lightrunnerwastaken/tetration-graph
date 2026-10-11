"""Round 4: every-degree numerical Schur bordering, with a growing precision.

No interval claim: native matrix updates discard all ball radii. An independent
+300-bit run and direct solves diagnose roundoff and recursion errors.
Each invocation has a 900-second soft deadline and writes diagnostic checkpoints.
"""
import argparse,json,math,sys,time,gzip,pickle
from pathlib import Path
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from flint import arb,acb,acb_poly,ctx
from h0_schur_native import Border

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--model',choices=['kubisch','kneser_z'],default='kubisch')
    p.add_argument('--max-n',type=int,default=1000)
    p.add_argument('--min-n',type=int,default=50)
    p.add_argument('--extra-bits',type=int,default=0)
    p.add_argument('--deadline',type=float,default=900)
    p.add_argument('--tag',default='')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();start=time.time();M=a.max_n
    a.output_dir.mkdir(parents=True,exist_ok=True)
    out=a.output_dir/f'h0_schur_{a.model}{a.tag}.json'
    sys.path.insert(0,str(a.project_root/'o_operator'));import o_kern as K
    # The entire input matrix is created once with additional guard precision.
    # A generous initial precision avoids irrevocable early-step rounding;
    # state arithmetic then grows by one bit per degree.
    input_bits=1200+M+a.extra_bits+100
    J,_,_=K.matrizen(a.model,M,input_bits)
    rh=(arb(K.L_RE)**2+arb(K.L_IM)**2).sqrt().mid()
    B=Border(J,input_bits);del J
    print(f'{a.model}: input and native buffers ready after {time.time()-start:.1f}s',flush=True)
    ref=json.loads((a.project_root/'g_gpt_hilfssaetze/ergebnisse/h0_cauchy_kubisch.json').read_text(encoding='utf-8'))
    fixed=[('corner',K.Lc)]+[(k,complex(*v['x'])) for k,v in ref['cases'][0]['points'].items() if k.startswith(('bay_','lens_'))]
    record={'status':'high precision numerical midpoint scan; no interval or all-n certificate',
            'model':a.model,'max_n_requested':M,'input_bits':input_bits,
            'arithmetic_bits':'1200+n+extra_bits; midpoint after every matrix update',
            'extra_bits':a.extra_bits,'beta':ref['beta_used'],'circle_angles_deg':list(range(-180,180,5)),
            'fixed_points':{k:[x.real,x.imag] for k,x in fixed},'cases':[],'completed':False}
    first=1
    checkpoint=out.with_suffix('.checkpoint.pkl.gz')
    if a.resume:
        record=json.loads(out.read_text(encoding='utf-8'))
        with gzip.open(checkpoint,'rb') as f:state=pickle.load(f)
        assert state['model']==a.model and state['max_n']==M and state['extra_bits']==a.extra_bits
        last=state['n'];prevh=arb(state['h']).mid()
        for i,row in enumerate(state['inverse']):
            for j,v in enumerate(row):B.api.set_string(B.api.entry(B.inv,i,j),v,input_bits)
        first=last+1
        record.setdefault('sessions_seconds',[]).append(record['seconds'])
        print(f'resumed without a new inverse solve at n={last}',flush=True)
        del state
    else:
        prevh=arb(B.api.string(B.api.entry(B.J,0,0),int(input_bits*.30103)+12)).mid()
    for n in range(first,M+1):
        bits=1200+n+a.extra_bits;ctx.prec=bits
        B.step(n,bits)
        F=B.api;digits=int(bits*.30103)+12
        h=arb(F.string(B.h,digits)).mid();qn=arb(F.string(B.qnorm,digits)).mid()
        if n>=a.min_n:
            row=[arb(F.string(F.entry(B.inv,n,k),digits)).mid() for k in range(n+1)]
            pol=acb_poly(row)
            pts={}
            points=list(fixed)
            radius=K.RHO*n**(-ref['beta_used'])
            for ps in range(-180,180,5):
                angle=math.atan2(K.Lc.imag,K.Lc.real)+math.radians(ps)
                x=K.Lc+radius*complex(math.cos(angle),math.sin(angle))
                points.append((f'kb_{ps}',x))
            for name,x in points:
                val=qn*abs(pol(rh/acb(x.real,x.imag)))
                if not val.is_finite():raise ArithmeticError(f'nonfinite increment n={n} point={name}')
                pts[name]=float(val.mid())
            fmax=max(((k,v) for k,v in pts.items() if not k.startswith('kb_')),key=lambda kv:kv[1])
            cmax=max(((k,v) for k,v in pts.items() if k.startswith('kb_')),key=lambda kv:kv[1])
            # Whole-circle Taylor majorant, entirely finite but numerical.
            rad=arb(radius);L=acb(K.Lc.real,K.Lc.imag)
            denom=abs(L)**2-rad**2
            center=rh*L.conjugate()/denom;wrad=rh*rad/denom
            shifted=pol(acb_poly([center,1]))
            maj=qn*sum((abs(shifted[j])*wrad**j for j in range(n+1)),arb(0))
            # Relative pivot diagnostic: cancellation against the 2-norm
            # bound for the Schur contraction, not a singular-value bound.
            rr=[arb(F.string(F.entry(B.J,n,k),digits)).mid() for k in range(n)]
            rnorm=sum((v*v for v in rr),arb(0)).sqrt()
            d=arb(F.string(F.entry(B.J,n,n),digits)).mid()
            scale=abs(d)+rnorm*(qn*qn-1).sqrt()
            case={'n':n,'bits':bits,'h':h.str(60,radius=False),'h_ratio':float((h/prevh).mid()),
                  'Q_norm':float(qn.mid()),'relative_pivot_2norm':float((abs(h)/scale).mid()),
                  'h_over_diagonal':float((abs(h)/abs(d)).mid()),
                  'fixed_max':fmax[1],'fixed_argmax':fmax[0],
                  'circle_sample_max':cmax[1],'circle_argmax':cmax[0],
                  'circle_taylor_majorant':float(maj.mid()),'points':pts}
            record['cases'].append(case)
        prevh=h
        if n%25==0 or n==M:
            record['last_n']=n;record['seconds']=time.time()-start
            out.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
            if n>=a.min_n:
                c=record['cases'][-1]
                print(f"{a.model} n={n}: sqrt(n) D fixed={math.sqrt(n)*c['fixed_max']:.9g}, circle={math.sqrt(n)*c['circle_sample_max']:.9g}, majorant={math.sqrt(n)*c['circle_taylor_majorant']:.9g}; h ratio={c['h_ratio']:.9g}; {record['seconds']:.1f}s",flush=True)
        if time.time()-start>a.deadline and n<M:
            # Do not launch a long laptop calculation. Keep diagnostics,
            # but a stopped run is explicitly incomplete.
            record['stopped_at_deadline']=True
            digits=int(bits*.30103)+12
            state={'model':a.model,'max_n':M,'extra_bits':a.extra_bits,'n':n,
                   'h':h.str(digits,radius=False),
                   'inverse':[[B.api.string(B.api.entry(B.inv,i,j),digits) for j in range(n+1)] for i in range(n+1)]}
            with gzip.open(checkpoint,'wb',compresslevel=1) as f:pickle.dump(state,f,protocol=5)
            record['last_n']=n;break
    record['completed']=n==M;record['seconds']=time.time()-start
    if record['completed']:record.pop('stopped_at_deadline',None)
    out.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    B.close()

if __name__=='__main__':main()
