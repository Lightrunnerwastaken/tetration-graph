"""Independent direct solves at scan precision +300 bits, for selected n.

This bypasses the native bridge and every earlier Schur update. Numerical
midpoint solves only; accuracy counters cover evaluation, not input errors.
"""
import argparse,json,sys,time,math
from pathlib import Path
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from flint import arb,arb_mat,acb,acb_poly,ctx

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--model',choices=['kubisch','kneser_z'],default='kubisch')
    p.add_argument('--degrees',default='50,200,600,1000')
    a=p.parse_args();sys.path.insert(0,str(a.project_root/'o_operator'));import o_kern as K
    scan=json.loads((a.output_dir/f'h0_schur_{a.model}.json').read_text(encoding='utf-8'))
    ns=[int(s) for s in a.degrees.split(',')];by={c['n']:c for c in scan['cases']}
    target=a.output_dir/f'h0_schur_direct_{a.model}.json'
    output=json.loads(target.read_text(encoding='utf-8')) if target.exists() else {'status':'numerical independent direct solve, +300 bits; no interval claim', 'model':a.model,'cases':[]}
    totalstart=time.time()
    executed=[]
    for n in ns:
        if executed:
            previous=next(c for c in output['cases'] if c['n']==executed[-1])
            estimate=previous['seconds']*(n/previous['n'])**3*((1500+n)/(1500+previous['n']))**1.5
            if time.time()-totalstart+estimate>1050:
                print(f'Estimated session limit: run n={n} and later degrees separately.',flush=True)
                break
        start=time.time();bits=1500+n
        J,_,_=K.matrizen(a.model,n,bits)
        e=arb_mat(n+1,1,[arb(k==n) for k in range(n+1)])
        col=J.solve(e,algorithm='approx');row=J.transpose().solve(e,algorithm='approx')
        h=(1/col[n,0]).mid();qn=(sum((col[k,0]**2 for k in range(n+1)),arb(0)).sqrt()/abs(col[n,0])).mid()
        pol=acb_poly([row[k,0].mid() for k in range(n+1)])
        rh=(arb(K.L_RE)**2+arb(K.L_IM)**2).sqrt()
        pts={};accuracy=[]
        for name,coords in scan['fixed_points'].items():
            val=qn*abs(pol(rh/acb(*coords)));pts[name]=float(val.mid());accuracy.append(val.rel_accuracy_bits())
        rad=K.RHO*n**(-scan['beta'])
        for ps in scan['circle_angles_deg']:
            angle=math.atan2(K.Lc.imag,K.Lc.real)+math.radians(ps)
            x=K.Lc+rad*complex(math.cos(angle),math.sin(angle))
            val=qn*abs(pol(rh/acb(x.real,x.imag)));pts[f'kb_{ps}']=float(val.mid());accuracy.append(val.rel_accuracy_bits())
        baseline=by[n]
        L=acb(K.Lc.real,K.Lc.imag);radarb=arb(rad)
        denom=abs(L)**2-radarb**2
        center=rh*L.conjugate()/denom;wrad=rh*radarb/denom
        shifted=pol(acb_poly([center,1]))
        majorant=qn*sum((abs(shifted[j])*wrad**j for j in range(n+1)),arb(0))
        majorant_float=float(majorant.mid())
        differences={name:abs(val-baseline['points'][name])/(1+abs(val)) for name,val in pts.items()}
        hs=arb(baseline['h']);relh=abs((h-hs)/h)
        relq=abs(float(qn)/baseline['Q_norm']-1)
        case={'n':n,'bits':bits,'h':h.str(60,radius=False),'Q_norm':float(qn),
              'h_relative_difference':float(relh.mid()),'Q_relative_difference_float':relq,
              'max_increment_difference_scaled':max(differences.values()),
              'circle_taylor_majorant':majorant_float,
              'majorant_difference_scaled':abs(majorant_float-baseline['circle_taylor_majorant'])/(1+abs(majorant_float)),
              'argmax_difference':max(differences,key=differences.get),
              'minimum_evaluation_accuracy_bits':min(accuracy),'points':pts,'seconds':time.time()-start}
        output['cases']=[c for c in output['cases'] if c['n']!=n]+[case]
        output['cases'].sort(key=lambda c:c['n'])
        executed.append(n)
        target.write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f"direct {a.model} n={n}, difference={case['max_increment_difference_scaled']:.3g}, h rel={case['h_relative_difference']:.3g}, time={case['seconds']:.1f}s",flush=True)
        if time.time()-totalstart>900:
            print('900-second soft limit: remaining degrees can be run separately.',flush=True);break
    output.setdefault('sessions',[]).append({'degrees':executed,'seconds':time.time()-totalstart})
    target.write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':main()
