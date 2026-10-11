"""Rigorous finite-N degree-increment bound on a WHOLE shrinking circle.

Exact cubic parameters are computed from Lambert W with acb enclosures.
The matrix inverse uses arb's checked algorithm. The inverse-circle
image and a finite Taylor polynomial give a continuous-circle majorant.
This certifies only listed N; it is not an all-N H0 certificate.
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8')
sys.dont_write_bytecode=True
from flint import arb, acb, arb_mat, arb_poly, acb_poly, ctx


def exact_model(N,bits):
    ctx.prec=bits
    L=-acb(-1).lambertw(-1)
    rho=abs(L)
    mu=L.sqrt()
    q=(mu-1)/(L-L.conjugate())
    beta=q.imag/L.imag
    alpha=q.real-beta*L.real
    f=arb_poly([rho**2*alpha,1+rho**2*beta-2*L.real*alpha,alpha-2*L.real*beta,beta])
    fp=f.derivative()
    m=arb_poly([fp[0]])+fp[1]*f+fp[2]*f*f
    powers=[arb_poly([1])]
    for j in range(N):
        powers.append((powers[-1]*f).truncate(N+1))
    rp=[rho**j for j in range(N+1)]
    J=arb_mat(N+1,N+1,[((powers[j][k] if k<=powers[j].degree() else arb(0))+
                          (m[k-j] if 0<=k-j<=6 else arb(0)))*rp[k]/rp[j]
                         for k in range(N+1) for j in range(N+1)])
    return J,L,rho,mu


def inverse_circle_norm(Ji,center,wrad,progress=False):
    """Finite Taylor majorant of ||Ji (1,w,...,w^N)^T|| on |w-c|=R."""
    N=Ji.nrows()-1
    sums=[arb(0) for j in range(N+1)]
    for i in range(N+1):
        shifted=acb_poly([acb(Ji[i,k]) for k in range(N+1)])(acb_poly([center,1]))
        for j in range(N+1):
            sums[j]+=abs(shifted[j])**2
        if progress and i>0 and i%200==0:
            print(f'  vector Taylor rows {i}/{N}',flush=True)
    terms=[sums[j].sqrt()*wrad**j for j in range(N+1)]
    return sum(terms,arb(0)),terms


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--degrees',default='200,300,400,600,800,1000')
    p.add_argument('--extra-bits',type=int,default=0)
    p.add_argument('--tag',default='')
    p.add_argument('--certify-G',action='store_true',help='Also certify the full G on the whole circle using this same inverse.')
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    record={'status':'rigorous finite-N whole-circle upper bounds; checked arb solves; no all-N conclusion',
            'parameter_definition':'L=-W_{-1}(-1), mu=sqrt(L), real cubic fixes L and its conjugate',
            'radius':'rho*N^(-beta), beta=1/1.1182 exactly (comparison convention from O38)',
            'majorant':'||Q_N|| * sum_j |[w^j] a_N(c+w)| R^j; a_N(w)=sum_k (J_N^-1)_{Nk} w^k = A_N(w)/h_N',
            'cases':[]}
    for N in [int(x) for x in a.degrees.split(',')]:
        start=time.time()
        bits=int(2.0*N)+500+a.extra_bits
        J,L,rho,mu=exact_model(N,bits)
        # One checked full inverse reuses the same preconditioner for the
        # required last column and last row. Two separate checked solves
        # would repeat that expensive preconditioner.
        Ji=J.inv()
        W=[Ji[k,N] for k in range(N+1)]
        A=[Ji[N,k] for k in range(N+1)]
        del J
        print(f'N={N}: checked inverse, {time.time()-start:.1f}s',flush=True)
        if W[N].contains(0):
            raise ArithmeticError('Inverse pivot enclosure contains zero')
        qnorm=sum((W[k]**2 for k in range(N+1)),arb(0)).sqrt()/abs(W[N])
        apol=acb_poly([acb(A[k]) for k in range(N+1)])
        beta=arb(1)/arb('1.1182')
        rad=rho*(-beta*arb(N).log()).exp()
        denom=rho**2-rad**2
        center=rho*L.conjugate()/denom
        wrad=rho*rad/denom
        shifted=apol(acb_poly([center,1]))
        terms=[qnorm*abs(shifted[j])*wrad**j for j in range(N+1)]
        maj=sum(terms,arb(0))
        scaled=maj*arb(N).sqrt()
        exact_corner=qnorm*abs(apol(rho/L))*arb(N).sqrt()
        assert scaled.is_finite()
        case={'N':N,'bits':bits,'method':'arb_mat.inv checked full inverse',
              'sqrtN_increment_at_L_enclosure':str(exact_corner),
              'sqrtN_whole_circle_increment_upper':str(scaled.upper()),
              'whole_circle_majorant_enclosure':str(scaled),
              'certified_less_than_4':bool(scaled<4),
              'certified_less_than_3p13':bool(scaled<arb('3.13')),
              'Q_N_norm':str(qnorm),'inverse_pivot':str(W[N]),
              'Taylor_terms_scaled_first_12':[str(x*arb(N).sqrt()) for x in terms[:12]],
              'Taylor_tail_after_12_scaled':str(sum(terms[12:],arb(0))*arb(N).sqrt()),
              'seconds':time.time()-start}
        record['cases'].append(case)
        # Save the cheaper scalar certificate before the optional long
        # vector trace calculation; interruption must not lose it.
        (a.output_dir/f'h0_border_circle_cert{a.tag}.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
        if a.certify_G:
            gmajor,gterms=inverse_circle_norm(Ji,center,wrad,progress=True)
            assert gmajor.is_finite()
            case['G_whole_circle_over_sqrtN_upper']=str((gmajor/arb(N).sqrt()).upper())
            case['G_whole_circle_over_sqrtN_enclosure']=str(gmajor/arb(N).sqrt())
            case['G_whole_circle_over_sqrtN_less_than_1']=bool(gmajor/arb(N).sqrt()<1)
            case['G_Taylor_terms_first_12_over_sqrtN']=[(x/arb(N).sqrt()).str(25) for x in gterms[:12]]
            case['seconds']=time.time()-start
            print(f'N={N}: checked whole-circle G/sqrtN {(gmajor/arb(N).sqrt()).str(18)}',flush=True)
        del Ji
        (a.output_dir/f'h0_border_circle_cert{a.tag}.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f'N={N}: sqrtN*increment whole-circle {scaled.str(18)}, at L {exact_corner.str(18)}; {case["seconds"]:.1f}s',flush=True)


if __name__=='__main__':
    main()
