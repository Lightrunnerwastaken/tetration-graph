"""Round 3 exact algebra and fixed cubic transport geometry."""
import argparse
import json
import sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import sympy as S
from flint import arb, acb, ctx


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    x,w=S.symbols('x w')
    # A rational cubic is sufficient to check the exact finite identities.
    f=S.Rational(1,2)+S.Rational(7,8)*x+S.Rational(1,4)*x**2+S.Rational(1,40)*x**3
    m=S.diff(f,x).subs(x,f).expand()
    N=4
    J=S.Matrix(N+1,N+1,lambda k,j:(f**j).expand().coeff(x,k)+(m.coeff(x,k-j) if k>=j else 0))
    Ji=J.inv()
    A=J[:N,:N]
    qi=Ji*S.Matrix([w**k for k in range(N+1)])
    prev=A.inv()*S.Matrix([w**k for k in range(N)])
    pivot=Ji[N,N]
    Q=Ji[:,N]/pivot
    Aa=Ji[N,:].T/pivot
    h=1/pivot
    dif=qi-S.Matrix(list(prev)+[0])
    assert all(S.cancel(v)==0 for v in dif-Q*qi[N])
    assert all(S.cancel(v)==0 for v in dif-Q*((Aa.T*S.Matrix([w**k for k in range(N+1)]))[0]/h))
    assert all(v==0 for v in (J*Q)[:N])
    assert all(v==0 for v in (Aa.T*J)[:N])
    assert S.cancel((Aa.T*J*Q)[0]-h)==0

    # Companion eigenvectors and the separate geometric branch.
    t,aa,b,c,d,lam=S.symbols('t a b c d lam')
    T=S.Matrix([[b/(t-aa),c/(t-aa),d/(t-aa),0],[1,0,0,0],[0,1,0,0],[0,0,0,1/t]])
    eigen=S.Matrix([lam**2,lam,1,0])
    char=(t-aa)*lam**3-b*lam**2-c*lam-d
    defect=T*eigen-lam*eigen
    assert S.cancel(defect[0]+char/(t-aa))==0
    assert all(v==0 for v in defect[1:])
    assert T*S.Matrix([0,0,0,1])==S.Matrix([0,0,0,1/t])
    ls=S.symbols('l1:4')
    V=S.Matrix([[q**2 for q in ls]+[0],list(ls)+[0],[1,1,1,0],[0,0,0,1]])
    vdet=S.factor(V.det())
    assert S.expand(vdet-(ls[0]-ls[1])*(ls[0]-ls[2])*(ls[1]-ls[2]))==0

    ctx.prec=512
    L=-acb(-1).lambertw(-1)
    rho=abs(L)
    mu=L.sqrt()
    q=(mu-1)/(L-L.conjugate())
    beta=q.imag/L.imag
    alpha=q.real-beta*L.real
    fc=[rho**2*alpha,1+rho**2*beta-2*L.real*alpha,alpha-2*L.real*beta,beta]
    disc=4*fc[2]**2-12*fc[3]*fc[1]
    assert disc>0
    critical=[(-2*fc[2]+s*disc.sqrt())/(6*fc[3]) for s in [-1,1]]
    cv=[sum((fc[k]*r**k for k in range(4)),arb(0)) for r in critical]
    assert all(abs(v)<arb('0.476') for v in cv)
    assert rho>1
    # At t=L the other two inverse roots come from the quadratic quotient.
    qb=acb(fc[2])+fc[3]*L
    qc=acb(fc[1])+fc[2]*L+fc[3]*L**2
    qdisc=qb**2-4*fc[3]*qc
    others=[(-qb+s*qdisc.sqrt())/(2*fc[3]) for s in [-1,1]]
    invroots=[1/L]+[1/r for r in others]
    sep=[abs(invroots[i]-invroots[j]) for i in range(3) for j in range(i+1,3)]
    assert all(s>0 for s in sep)
    # Track roots on |t|=rho only as a numerical cross-check of the analytic
    # three-cycle monodromy (both critical values lie inside this circle).
    coeff=np.array([float(v.mid()) for v in fc],dtype=complex)
    roots=None
    first=None
    from itertools import permutations
    for phi in np.linspace(0,2*np.pi,721):
        pol=coeff.copy()
        pol[0]-=float(rho.mid())*np.exp(1j*phi)
        rr=np.roots(pol[::-1])
        if roots is None:
            roots=rr
            first=rr.copy()
        else:
            perm=min(permutations(range(3)),key=lambda p:sum(abs(rr[p[k]]-roots[k])**2 for k in range(3)))
            roots=rr[list(perm)]
    endpoint=[int(np.argmin(abs(first-r))) for r in roots]
    assert sorted(endpoint)==[0,1,2] and all(endpoint[k]!=k for k in range(3))
    result={
        'status':'exact symbolic finite identities; rigorous 512-bit geometry; monodromy tracking separately numerical',
        'bordered_increment_identity_exact':True,
        'biorthogonal_rank_one_identity_exact':True,
        'companion_eigenbasis_exact':True,
        'eigenbasis_determinant':str(vdet),
        'critical_points':[r.str(30) for r in critical],
        'critical_values':[r.str(30) for r in cv],
        'both_critical_values_in_abs_t_lt_0p476':True,
        'inverse_eigenvalue_separations_at_L':[r.str(30) for r in sep],
        'other_preimages_of_L':[[r.real.str(30),r.imag.str(30)] for r in others],
        'root_monodromy_numerical_permutation':endpoint,
        'root_monodromy_tracking_nodes':721,
        'analytic_monodromy_reason':'circle encloses both critical values; deformation to infinity gives the cubic three-cycle',
    }
    a.output_dir.mkdir(parents=True,exist_ok=True)
    (a.output_dir/'h0_border_exact_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
