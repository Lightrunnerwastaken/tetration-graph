"""Exact kernel recurrence, finite-rank formula and contour geometry checks."""
import argparse
import json
from pathlib import Path
import sys
import sympy as sp
from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    outdir = (args.output_dir or package / "ergebnisse").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    z, t = sp.symbols("z t")
    a, b, c, d, mt = sp.symbols("a b c d mt")
    F = a + b*z + c*z**2 + d*z**3
    denominator = t*(t-a) - ((t-a)+b*t)*z + (b-c*t)*z**2 + (c-d*t)*z**3 + d*z**4
    assert sp.expand(denominator-(t-F)*(t-z)) == 0
    numerator = t-z+mt*(t-F)
    assert sp.cancel((1/(t-F)+mt/(t-z))*denominator-numerator) == 0
    # The degree-three numerator proves the homogeneous recurrence from order 4.
    assert sp.degree(numerator, z) <= 3
    # Generic nonpolynomial m: the common integral kernel must use m(z), not m(t).
    # Finite-rank projection for a degree-six polynomial m, tested symbolically.
    ms = sp.symbols("m0:7")
    ys = sp.symbols("y0:10")
    m = sum(ms[j]*z**j for j in range(7))
    Y = sum(ys[k]*z**(-k) for k in range(10))
    B = sum(sum(ys[k]*ms[j+k+1] for k in range(6-j))*z**j for j in range(6))
    negative = sum(sum(ys[k]*ms[k-j] for k in range(j,min(9,j+6)+1))*z**(-j-1) for j in range(10))
    assert sp.expand(m*Y/z - negative - B) == 0
    ctx.prec = 512
    L = -acb(-1).lambertw(-1)
    rho = abs(L)
    mu = L.sqrt()
    quotient = (mu-1)/(L-L.conjugate())
    beta = quotient.imag / L.imag
    alpha = quotient.real - beta*L.real
    rr = L.real**2+L.imag**2
    coefficients = [rr*alpha, 1+rr*beta-2*L.real*alpha, alpha-2*L.real*beta, beta]
    radius = 16*rho
    margin = abs(coefficients[3])*radius**3 - sum((abs(coefficients[k])*radius**k for k in range(3)), arb(0)) - arb("1.5")*rho
    assert margin.is_finite() and margin > 0
    kappa = 2*acb(0,1)*arb.pi() / (2*acb(0,1)*arb.pi() - mu.log())
    result = {
        "status": "exact rational-symbolic checks and rigorous 512-bit geometric margin",
        "cubic_forward_kernel_order_four_recurrence": True,
        "degree_six_multiplier_B_depends_on_first_six_y_coefficients": True,
        "large_contour_radius_over_rho": 16,
        "target_radius_over_rho": "1.5",
        "rouche_margin": margin.str(30),
        "kappa_interval": kappa.str(30),
        "beta_interval": (1/kappa.real).str(30),
        "note": "No interval certificate of the inverse matrices or H0 is claimed.",
    }
    (outdir / "h0_exact_checks.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__ == "__main__":
    main()
