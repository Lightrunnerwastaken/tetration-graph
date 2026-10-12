"""Stage 3: complex degree-increment profiles and both Schur vectors.

900 s soft time limit; saving a full native checkpoint takes additional time.
Use an existing round-4 n=914 checkpoint only for n>914, never for earlier n.
All outputs are numerical midpoint diagnostics.
"""
import argparse,gzip,json,math,pickle,sys,time
from pathlib import Path
sys.dont_write_bytecode=True;sys.stdout.reconfigure(encoding='utf-8')
from flint import arb,acb,acb_poly,ctx
from h0_schur_native import Border

def comp(z):return [z.real.mid().str(45,radius=False),z.imag.mid().str(45,radius=False)]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--model',choices=['kubisch','kneser_z'],required=True);p.add_argument('--max-n',type=int,required=True)
    p.add_argument('--degrees',required=True);p.add_argument('--checkpoint',type=Path)
    p.add_argument('--checkpoint-out',type=Path,
                   help='Optional large intermediate inverse checkpoint for a separate bounded continuation.')
    p.add_argument('--extra-bits',type=int,default=0);p.add_argument('--deadline',type=float,default=850)
    a=p.parse_args();start=time.time();a.output_dir.mkdir(parents=True,exist_ok=True)
    degrees=set(map(int,a.degrees.split(',')));sys.path.insert(0,str(a.project_root/'o_operator'));import o_kern as K
    input_bits=1300+a.max_n+a.extra_bits;J,_,_=K.matrizen(a.model,a.max_n,input_bits)
    B=Border(J,input_bits);del J;first=1
    if a.checkpoint:
        with gzip.open(a.checkpoint,'rb') as f:state=pickle.load(f)
        assert state['model']==a.model and state['max_n']<=a.max_n and state['extra_bits']==a.extra_bits
        first=state['n']+1;assert all(n>=first for n in degrees)
        for i,row in enumerate(state['inverse']):
            for j,v in enumerate(row):B.api.set_string(B.api.entry(B.inv,i,j),v,input_bits)
        del state
    ref=json.loads((a.project_root/'g_gpt_hilfssaetze/ergebnisse/h0_cauchy_kubisch.json').read_text(encoding='utf-8'))
    fixed=[('corner',K.Lc)]+[(k,complex(*v['x'])) for k,v in ref['cases'][0]['points'].items() if k.startswith(('bay_','lens_'))]
    cases=[];output=a.output_dir/f'profiles_{a.model}_{first}_{a.max_n}.json'
    print(f'{a.model}: input prepared, begin n={first}, {time.time()-start:.1f}s',flush=True)
    for n in range(first,a.max_n+1):
        ctx.prec=1200+n+a.extra_bits;B.step(n,ctx.prec)
        if n in degrees:
            api=B.api;digits=int(ctx.prec*.30103)+12
            h=arb(api.string(B.h,digits)).mid();q=arb(api.string(B.qnorm,digits)).mid()
            row=[arb(api.string(api.entry(B.inv,n,j),digits)).mid() for j in range(n+1)]
            col=[arb(api.string(api.entry(B.inv,j,n),digits)).mid() for j in range(n+1)]
            pol=acb_poly(row);rho=(arb(K.L_RE)**2+arb(K.L_IM)**2).sqrt().mid();L=acb(K.Lc.real,K.Lc.imag)
            rad=arb(K.RHO*n**(-ref['beta_used']));den=abs(L)**2-rad**2
            center=rho*L.conjugate()/den;wrad=rho*rad/den
            scale=q*arb(n).sqrt();shift=pol(acb_poly([center,wrad]))*scale
            corner=pol(rho/L)*scale;chi=corner/abs(corner)
            vals={name:comp(pol(rho/acb(x.real,x.imag))*scale) for name,x in fixed}
            circles=[]
            for deg in range(0,360,3):
                s=complex(math.cos(math.radians(deg)),math.sin(math.radians(deg)))
                x=K.Lc+float(rad)*s
                circles.append(comp(pol(rho/acb(x.real,x.imag))*scale))
            case={'n':n,'h':h.str(90,radius=False),'Q_norm':q.str(90,radius=False),
                  'A_norm':(abs(h)*sum((v*v for v in row),arb(0)).sqrt()).str(90,radius=False),
                  'Q_coefficients':[(h*v).str(70,radius=False) for v in col],
                  'A_coefficients':[(h*v).str(70,radius=False) for v in row],
                  'phase_corner':comp(chi),'corner_scaled':comp(corner),'fixed_scaled':vals,
                  'x_circle_scaled_angles_deg':list(range(0,360,3)),'x_circle_scaled':circles,
                  'w_center':comp(center),'w_radius':wrad.str(45,radius=False),
                  'w_circle_taylor_scaled':[comp(shift[k]) for k in range(n+1)],
                  'scaled_circle_majorant':float(sum((abs(shift[k]) for k in range(n+1)),arb(0)).mid())}
            cases.append(case)
            record={'model':a.model,'status':'high precision numerical complex profiles; no all-n or input-error certificate',
                    'bits':'1200+n+extra_bits','extra_bits':a.extra_bits,'start_n':first,'max_n':a.max_n,
                    'beta':ref['beta_used'],'fixed_points':{name:[x.real,x.imag] for name,x in fixed},
                    'cases':cases,'seconds':time.time()-start}
            output.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
            print(f'{a.model} n={n}, sqrt(n) circle majorant={case["scaled_circle_majorant"]:.9f}, {time.time()-start:.1f}s',flush=True)
        if time.time()-start>a.deadline and n<a.max_n:
            raise TimeoutError(f'Soft deadline at n={n}; already saved profiles remain valid. Use separate bounded batches.')
    if a.checkpoint_out:
        # Large scratch artifact, intentionally separate from deliverable JSON.
        digits=int(ctx.prec*.30103)+12
        state={'model':a.model,'max_n':a.max_n,'n':a.max_n,'extra_bits':a.extra_bits,'inverse':[]}
        for i in range(a.max_n+1):
            if time.time()-start>1050:
                B.close()
                raise TimeoutError('Checkpoint export budget exceeded; saved profiles remain valid.')
            state['inverse'].append([B.api.string(B.api.entry(B.inv,i,j),digits)
                                     for j in range(a.max_n+1)])
        a.checkpoint_out.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(a.checkpoint_out,'wb',compresslevel=1) as f:pickle.dump(state,f,protocol=5)
        del state
    B.close();record['completed']=True;record['seconds']=time.time()-start
    output.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'completed {a.model}: {record["seconds"]:.1f}s',flush=True)

if __name__=='__main__':main()
