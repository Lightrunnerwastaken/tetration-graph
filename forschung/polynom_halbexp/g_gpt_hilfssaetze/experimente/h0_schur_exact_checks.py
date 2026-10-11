"""Exact rational checks of the paired recurrence and the kernel displacement.

Also checks the native matrix bridge against a direct checked small inverse.
These are finite algebra/API checks, not a proof of G19.
"""
import json,sys,time,argparse
from pathlib import Path
sys.dont_write_bytecode=True
import sympy as s
from flint import arb,arb_mat,ctx
from h0_schur_native import Border

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);start=time.time()
    z=s.symbols('z');rho=s.Rational(3,2)
    f=s.Rational(1,2)+s.Rational(3,4)*z+s.Rational(1,4)*z**2+s.Rational(1,32)*z**3
    m=s.expand(s.diff(f,z).subs(z,f));F=s.expand(f.subs(z,rho*z)/rho);M=s.expand(m.subs(z,rho*z))
    cols=[s.Poly(s.expand(F**j+M*z**j),z) for j in range(9)]
    J=s.Matrix(9,9,lambda k,j:cols[j].nth(k))
    As=[];Qs=[];hs=[];checks=[]
    for n in range(9):
        block=J[:n+1,:n+1];inv=block.inv();h=1/inv[n,n]
        Q=inv[:,n]*h;A=inv[n,:].T*h
        if n:
            b=J[:n,n];r=J[n,:n];d=J[n,n]
            predictedQ=s.eye(n+1)[:,n];predictedA=s.eye(n+1)[:,n];predictedh=d
            for j in range(n):
                qpad=Qs[j].col_join(s.zeros(n+1-len(Qs[j]),1))
                apad=As[j].col_join(s.zeros(n+1-len(As[j]),1))
                gamma=(r*qpad[:n,0])[0];eta=(apad[:n,0].T*b)[0]
                predictedQ-=qpad*eta/hs[j]
                predictedA-=apad*gamma/hs[j]
                predictedh-=gamma*eta/hs[j]
            assert Q==predictedQ and A==predictedA and h==predictedh
            # Evaluate the source-preserving functional identity, exactly.
            w=s.Rational(5,3)+s.I*s.Rational(2,7)
            ev=s.Matrix([w**k for k in range(n+1)])
            increment=inv*ev-(J[:n,:n].inv()*ev[:n,0]).col_join(s.zeros(1,1))
            assert all(s.expand(t)==0 for t in increment-Q*(A.T*ev)[0]/h)
        As.append(A);Qs.append(Q);hs.append(h);checks.append(n)
    for j in range(7):
        assert s.expand(cols[j+2].as_expr()-(F+z)*cols[j+1].as_expr()+z*F*cols[j].as_expr())==0
    R,u,v=s.symbols('R u v',positive=True)
    center_gap=R**2*(u**2-v**2)/((R**2-u**2)*(R**2-v**2))
    old_radius=R*u/(R**2-u**2);new_radius=R*v/(R**2-v**2)
    margin=(u-v)*(R-u)*(R-v)/(u*(R**2-v**2))
    assert s.factor(1-(center_gap+new_radius)/old_radius-margin)==0
    ctx.prec=384
    Ja=arb_mat(9,9,[arb(str(v.p))/arb(str(v.q)) for v in J])
    B=Border(Ja,384)
    error=arb(0)
    for n in range(1,9):
        B.step(n,384)
        exact=J[:n+1,:n+1].inv()
        for i in range(n+1):
            for j in range(n+1):
                got=arb(B.api.string(B.api.entry(B.inv,i,j),130))
                ex=arb(str(exact[i,j].p))/arb(str(exact[i,j].q))
                error=max(error,abs(got-ex).upper())
    B.close();assert error<arb('1e-100')
    out={'status':'exact rational finite checks and numerical ABI comparison; no all-n claim',
         'pair_recursion_degrees':checks,'evaluated_source_identity':True,
         'two_state_column_identity_degrees':list(range(7)),
         'pullback_margin_rational_identity':True,
         'native_bridge_max_error_upper':str(error),'native_bits':384,'seconds':time.time()-start}
    (a.output_dir/'h0_schur_exact_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
