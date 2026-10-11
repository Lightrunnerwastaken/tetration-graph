"""H0 round 3: exact source-preserving degree increments, numerical probes.

Uses W25.1's existing biorthogonal expansion as a bordered inverse identity.
No SVD, no K2, no nonlinear branch construction. Approximate arb solves;
the finite point sample and decimal model parameters are not certificates.
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.stdout.reconfigure(encoding='utf-8')
from flint import acb, arb, arb_mat, ctx


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--model', choices=['kubisch','kneser_z'], default='kubisch')
    p.add_argument('--degrees', default='200,300,400,600,800,1000')
    p.add_argument('--extra-bits', type=int, default=0)
    p.add_argument('--tag', default='')
    a=p.parse_args()
    ns=[int(x) for x in a.degrees.split(',')]
    if a.model=='kneser_z' and max(ns)>639:
        raise ValueError('J_N requires coefficients through N+1; certified source ends at 640.')
    sys.path.insert(0,str(a.project_root/'o_operator'))
    import o_kern as K
    ref=json.loads((a.project_root/'g_gpt_hilfssaetze/ergebnisse/h0_cauchy_kubisch.json').read_text(encoding='utf-8'))
    refpoints=ref['cases'][0]['points']
    fixed=[('corner',K.Lc)]+[(k,complex(*v['x'])) for k,v in refpoints.items() if k.startswith(('bay_','lens_'))]
    a.output_dir.mkdir(parents=True,exist_ok=True)
    out=a.output_dir/f'h0_border_{a.model}{a.tag}.json'
    record={'status':'numerical diagnostics; approximate arb solves, no interval/continuous-curve certificate',
            'model':a.model,'fixed_contour_source':'round 2 fixed cubic gamma* points; includes L',
            'increment':'||Q_N|| * |(J_N^-1 e_x)_N|; Q_N=(J_N^-1 e_N)/(J_N^-1 e_N)_N',
            'beta_used':ref['beta_used'],'cases':[]}
    for N in ns:
        start=time.time()
        bits=int(1.6*N)+400+a.extra_bits
        J,_,_=K.matrizen(a.model,N,bits)
        rh=(arb(K.L_RE)**2+arb(K.L_IM)**2).sqrt()
        points=list(fixed)
        for ps in range(-180,180,15):
            angle=math.atan2(K.Lc.imag,K.Lc.real)+math.radians(ps)
            x=K.Lc+K.RHO*N**(-ref['beta_used'])*complex(math.cos(angle),math.sin(angle))
            points.append((f'kb_{ps}',x))
        columns=[]
        for _,x in points:
            z=acb(x.real,x.imag)
            ee=[(rh/z)**k for k in range(N+1)]
            columns += [[x.real for x in ee],[x.imag for x in ee]]
        columns.append([arb(int(k==N)) for k in range(N+1)])
        rhs=arb_mat(N+1,len(columns),[columns[j][k] for k in range(N+1) for j in range(len(columns))])
        V=J.solve(rhs,algorithm='approx')
        lastcol=len(columns)-1
        pivot=V[N,lastcol]
        if pivot.contains(0):
            raise ArithmeticError('zero/indeterminate bordered inverse pivot')
        qnorm=(sum((V[k,lastcol]**2 for k in range(N+1)),arb(0)).sqrt()/abs(pivot))
        case={'N':N,'bits':bits,'Q_N_coefficient_norm':float(qnorm.mid()),
              'inverse_pivot':str(pivot),'Schur_complement':str(1/pivot),'points':{}}
        for j,(name,x) in enumerate(points):
            last=acb(V[N,2*j],V[N,2*j+1])
            inc=qnorm*abs(last)
            vnorm=sum((V[k,2*j]**2+V[k,2*j+1]**2 for k in range(N+1)),arb(0)).sqrt()
            assert inc.is_finite() and vnorm.is_finite()
            case['points'][name]={'x':[x.real,x.imag],'G':float(vnorm.mid()),
                                   'degree_increment_norm':float(inc.mid()),
                                   'sqrtN_times_increment':float((inc*arb(N).sqrt()).mid()),
                                   'last_solution_coefficient_abs':str(abs(last))}
        fvals=[(k,v['sqrtN_times_increment']) for k,v in case['points'].items() if not k.startswith('kb_')]
        kvals=[(k,v['sqrtN_times_increment']) for k,v in case['points'].items() if k.startswith('kb_')]
        case['fixed_gamma_max_sqrtN_increment']=max(v for _,v in fvals)
        case['fixed_gamma_argmax']=max(fvals,key=lambda v:v[1])[0]
        case['circle_max_sqrtN_increment']=max(v for _,v in kvals)
        case['seconds']=time.time()-start
        record['cases'].append(case)
        out.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f"{a.model} N={N}: sqrtN*increment L={case['points']['corner']['sqrtN_times_increment']:.9g}, gamma={case['fixed_gamma_max_sqrtN_increment']:.9g} ({case['fixed_gamma_argmax']}), circle={case['circle_max_sqrtN_increment']:.9g}; {case['seconds']:.1f}s",flush=True)


if __name__=='__main__':
    main()
