"""Exact algebraic checks and high-precision examples for the research note."""
import sys, json
from pathlib import Path
import sympy as sp
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(__file__).resolve().parents[1] / "outputs"

def exact_adjoint():
    z, t, w = sp.symbols("z t w")
    n = 5
    a = sp.Rational(2, 5)
    F = a + z + z*z/7 + sum(z**k/sp.factorial(k)/13 for k in range(3,n+1))
    m = sum(z**k/sp.factorial(k) for k in range(n+1))
    ys = [sp.Rational((-1)**k,k+1) for k in range(n+1)]
    A = sp.Matrix(n+1,n+1,lambda k,j:sp.expand(F**j+m*z**j).coeff(z,k))
    us = list(A.T * sp.Matrix(ys))
    cs = [sum(ys[k]*sp.expand((F-a)**j).coeff(z,k) for k in range(n+1)) for j in range(n+1)]
    ds = [sum(ys[k]*m.coeff(z,k-j) for k in range(j,n+1)) for j in range(n+1)]
    R = sum(cs[j]/(t-a)**(j+1) for j in range(n+1))
    M = sum(ds[j]/t**(j+1) for j in range(n+1))
    U = sum(us[j]/t**(j+1) for j in range(n+1))
    residual = sp.cancel(R+M-U)
    numerator = sp.cancel(residual*t**(n+1)*(t-a)**(n+1))
    assert sp.denom(numerator) == 1
    assert sp.degree(numerator,t) <= n
    series = sp.series(residual.subs(t,1/w),w,0,n+2).removeO()
    assert series == 0
    return {"N":n,"exact_check":True,"numerator_degree":int(sp.degree(numerator,t)),
            "note":"F and m are the N-jets of nonpolynomial analytic germs; the test is exact rational algebra."}

def multiply_series(a,b,n):
    return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(n+1)]

def composition(p,n):
    inner = p[:n+1]
    v = [mp.mpf(0)]*(n+1)
    for coefficient in reversed(p):
        v = multiply_series(v,inner,n)
        v[0] += coefficient
    return v

def polynomial_at_series(p,inner,n):
    v = [mp.mpf(0)]*(n+1)
    for coefficient in reversed(p):
        v = multiply_series(v,inner,n)
        v[0] += coefficient
    return v

def finite_hermite_correction(b,degree,n,rho):
    alpha = b[0]
    inv = [mp.mpf(0)]*(n+1)
    germ = [mp.mpf(0)] + b[1:n+1]
    for j in range(1,n+1):
        value = polynomial_at_series(germ,inv,n)[j]
        inv[j] = ((1 if j==1 else 0)-value)/b[1]
    exponential = [mp.mpf(1)] + [mp.mpf(0)]*n
    for j in range(1,n+1):
        exponential[j] = sum(k*rho*inv[k]*exponential[j-k] for k in range(1,j+1))/j
    outer = [exponential[j]/rho for j in range(n+1)]
    translated = [sum(b[k]*mp.binomial(k,j)*alpha**(k-j)
                      for k in range(j,degree+1)) for j in range(n+1)]
    inverse_power = [(-1)**j*mp.binomial(n+j,j)*alpha**(-n-1-j) for j in range(n+1)]
    return multiply_series([outer[j]-translated[j] for j in range(n+1)],inverse_power,n)

def oversampling():
    mp.mp.dps = 160
    base = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam")
    data = json.loads((base/"kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))
    L = -mp.lambertw(-1,-1)
    rho = abs(L)
    # F(t)=f(rho*t)/rho, hence F(F(t))=exp(rho*t)/rho.
    b = [mp.mpf(row[0])*rho**(k-1) for k,row in enumerate(data["a"])]
    alpha = b[0]
    cases = []
    for n in [5,10,20,40]:
        degree = 5*n
        c = finite_hermite_correction(b,degree,n,rho)
        delta = [mp.mpf(0)]*(degree+1)
        for j in range(n+1):
            for ell in range(j+1):
                delta[n+1+ell] += c[j]*mp.binomial(j,ell)*(-alpha)**(j-ell)
        p = [b[k]+delta[k] for k in range(degree+1)]
        comp = composition(p,n)
        defect = max(abs(comp[k]-rho**(k-1)/mp.factorial(k)) for k in range(n+1))
        corr = rho*sum(abs(x) for x in delta)
        tail = rho*sum(abs(x) for x in b[degree+1:])
        cases.append({"matched_order":n,"degree":degree,
                      "maximum_jet_defect":mp.nstr(defect,12),
                      "correction_l1_rho":mp.nstr(corr,12),
                      "known_tail_l1_rho_through_640":mp.nstr(tail,12)})
    return {"status":"numerical, not an interval certificate; coefficient midpoints used",
            "precision_decimal_digits":160,"alpha":mp.nstr(alpha,30),
            "geometric_base_bound":str(mp.mpf("0.53")**4*mp.mpf("1.37")/mp.mpf("0.16")),
            "cases":cases}

if __name__ == "__main__":
    result={"exact_adjoint":exact_adjoint(),"oversampling":oversampling()}
    (OUT/"lemma_checks.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
