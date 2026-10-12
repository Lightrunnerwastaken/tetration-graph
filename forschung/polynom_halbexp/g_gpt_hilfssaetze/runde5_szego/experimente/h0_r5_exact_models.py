"""Exact finite algebra checks, no asymptotic proof by testing."""
import json,sys,time
from pathlib import Path
from h0_r5_paths import locations
import sympy as s
sys.stdout.reconfigure(encoding='utf-8')
start=time.time();z,w=s.symbols('z w');checks=[]
families=[('pure_cubic',s.Rational(1,3)+s.Rational(3,2)*z+z*z/7+z**3/11,0),
          ('triangular_quadratic',s.Rational(3,2)*z+z*z/7,2+z/3+z*z/5),
          ('affine_M',2*z-1,2),
          ('triangular_boundary',(z+z*z)/2,s.Rational(1,2)+(z+z*z)/2)]
for name,F,M in families:
    for n in range(7):
        J=s.Matrix(n+1,n+1,lambda k,j:s.Poly(s.expand(F**j+M*z**j),z).nth(k))
        det=J.det();b=s.diff(F,z).subs(z,0)
        expected=b**(n*(n+1)//2) if name=='pure_cubic' else s.prod(b**j+s.sympify(M).subs(z,0) for j in range(n+1))
        assert det==expected
        if name=='pure_cubic':
            inv=J.inv();h=1/inv[n,n];Q=inv[:,n]*h
            assert sum(Q[k]*z**k for k in range(n+1))==s.Poly((z-F.subs(z,0))**n,z).as_expr()
        if name=='affine_M':
            inv=J.inv();h=1/inv[n,n];A=inv[n,:].T*h
            assert A==s.eye(n+1)[:,n]
        checks.append([name,n])
# Andreief / rowwise branch expansion, exactly for n=1,2.
F=s.Rational(1,3)+2*z+z*z/5;M=1+z/7
for n in [1,2]:
    zs=s.symbols('z0:'+str(n+1));pol=[s.expand(F**j+M*z**j) for j in range(n+1)]
    J=s.Matrix(n+1,n+1,lambda k,j:s.Poly(pol[j],z).nth(k))
    minus=s.Matrix(n+1,n+1,lambda i,k:zs[i]**(-k)).det()
    pp=s.Matrix(n+1,n+1,lambda i,j:pol[j].subs(z,zs[i])).det()
    integ=s.expand(minus*pp)
    for zz in zs:integ=s.expand(integ).coeff(zz,0)
    assert s.simplify(integ/s.factorial(n+1)-J.det())==0
out={'status':'exact finite checks; written proofs are separate','determinants_and_pairs':checks,
     'Andreief_degrees':[1,2],'seconds':time.time()-start}
folder,ROUND4,ROOT=locations(__file__);folder.mkdir(exist_ok=True)
(folder/'exact_models.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(out)
