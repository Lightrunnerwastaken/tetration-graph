"""Exact structural identities and explicit rational norm constants."""
import json,sys,time
from pathlib import Path
from h0_r5_paths import locations
import sympy as s
sys.stdout.reconfigure(encoding='utf-8')
start=time.time();z,t=s.symbols('z t');F=s.Rational(1,3)+s.Rational(7,8)*z+z*z/7+z**3/13
M=1+z/3+z*z/5;checks=[]
for n in range(2,9):
    p=[s.Poly(s.expand(F**j+M*z**j),z) for j in range(n+2)]
    J=s.Matrix(n+1,n+1,lambda k,j:p[j].nth(k))
    S=s.Matrix(n+1,n+1,lambda k,j:int(k==j+1))
    TF=s.Matrix(n+1,n+1,lambda k,j:s.Poly(F,z).nth(k-j) if k>=j else 0)
    L=TF+S;B=S*TF;R=J*S*S-L*J*S+B*J
    U=s.Matrix(n+1,2,lambda k,j:-p[n+1].nth(k) if j==0 else s.Poly(z*F*p[n].as_expr(),z).nth(k))
    V=s.Matrix(n+1,2,lambda k,j:int(k==n-1+j))
    assert R==U*V.T and R.rank()<=2
    companion=L.row_join(-B).col_join(s.eye(n+1).row_join(s.zeros(n+1)))
    H=(J*S).col_join(J)
    assert companion*H-H*S==(-R).col_join(s.zeros(n+1))
    checks.append(n)
assert s.simplify(1/(1-t*F)+M/(1-t*z)-(1+M-t*(z+M*F))/((1-t*F)*(1-t*z)))==0
# Bidiagonal exact inverses.
b=s.Rational(2);m=s.Rational(2);epsilon=s.Rational(1,4)
for n in range(9):
    J=s.diag(*[b**j+m for j in range(n+1)])
    for k in range(1,n+1):J[k,k-1]=epsilon
    inv=J.inv()
    for k in range(n+1):
        assert inv[n,k]==(-epsilon)**(n-k)/s.prod(b**j+m for j in range(k,n+1))
# Exact forward/backward jet equivalence in the nonlinear quadratic family.
quadratic_jet_checks=[]
for eta_check in [s.Rational(-1,20),s.Rational(0),s.Rational(1,20)]:
    degree=6;Fq=2*z+eta_check*z*z;Mq=2+2*eta_check*Fq
    H=z/2 if eta_check==0 else s.series((-1+s.sqrt(1+eta_check*z))/eta_check,z,0,degree+1).removeO()
    CH=s.Matrix(degree+1,degree+1,lambda k,j:s.expand(H**j).coeff(z,k))
    backward=s.Matrix(degree+1,degree+1,lambda k,j:
        s.expand((2+2*eta_check*z)*H**j).coeff(z,k))+s.eye(degree+1)
    forward=s.Matrix(degree+1,degree+1,lambda k,j:s.expand(Fq**j+Mq*z**j).coeff(z,k))
    assert backward.inv()*CH==forward.inv()
    quadratic_jet_checks.append(str(eta_check))
# Global inverse-branch / head-tail model: exact rational inequalities.
R=s.Rational(3,2);r=s.Rational(4,5);eta=s.Rational(1,20);W=s.Rational(43,20)
assert 2*r-eta*r*r>R
tail=W*(r/R)**4
assert tail<s.Rational(1,5)
# Orthonormal Hardy basis z^j/R^j: entries gain R^(row-column).
off_frobenius_squared=(2*eta*R)**2+(s.Rational(3,4)*eta*R)**2+(eta**2*R**2/8)**2+(eta*R/4)**2
assert off_frobenius_squared<s.Rational(1,16)
full_inverse_bound_squared=(1+s.Rational(5,4)*W)**2+s.Rational(25,16)
assert full_inverse_bound_squared<16
out={'status':'exact finite identity checks plus exact rational bounds; full proofs in ETAPPE4',
     'rank_two_displacement_degrees':checks,'two_variable_generating_identity':True,
     'bidiagonal_inverse_degrees':list(range(9)),
     'quadratic_forward_backward_jet_degree':6,'quadratic_eta_checks':quadratic_jet_checks,
     'quadratic_branch_margin':str(2*r-eta*r*r-R),
     'tail_operator_bound':str(tail),'head_off_frobenius_squared':str(off_frobenius_squared),
     'full_inverse_bound_squared':str(full_inverse_bound_squared),'seconds':time.time()-start}
folder,ROUND4,ROOT=locations(__file__)
(folder/'structure_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
