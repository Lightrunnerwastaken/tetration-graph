"""Common Cauchy kernel and weighted trace diagnostic for the cubic.

The Laurent polynomial is taken from W's saved approximate corner vector.
Its stored coefficients are treated as exact dyadic input for this diagnostic.
All contour evaluations use 2400-bit acb arithmetic. Not a certificate of H0.
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from flint import acb, acb_poly, arb, arb_poly, ctx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--degrees", default="200,600,1000")
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    project = (args.project_root or package.parent).resolve()
    outdir = (args.output_dir or package / "ergebnisse").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(project / "w_grenzschicht"))
    import w_kern as W
    from w40_lade import lade
    ctx.prec = 2400
    H, Lam, _, _ = lade()
    rho = arb(float(W.RHO))
    first = np.load(project / "w_grenzschicht" / "w01_N200.npz")
    coefficients = [arb(float(x)) for x in first["coef"]]
    f = acb_poly([acb(x) for x in coefficients])
    fp = f.derivative()
    m = acb_poly([0])
    power = acb_poly([1])
    for k in range(3):
        m += fp[k] * power
        power *= f
    L = acb(float(W.L.real), float(W.L.imag))
    a = coefficients[0]
    etaA = acb(a + abs(L - a))
    target = acb("1.14")
    for _ in range(40):
        target -= (f(target) - etaA) / fp(target)
        target = acb(target.real.mid(), target.imag.mid())
    raw = [complex(float(x)) for x in coefficients]
    raw[0] -= complex(float(target.real), float(target.imag))
    roots = []
    for initial in np.roots(list(reversed(raw))):
        root = acb(float(initial.real), float(initial.imag))
        for _ in range(40):
            root -= (f(root) - target) / fp(root)
            root = acb(root.real.mid(), root.imag.mid())
        roots.append(root)
    far = [root for root in roots if abs(root) > rho]
    assert len(far) == 2
    quadrature = {}
    for count in [2048, 4096]:
        ph = 2 * np.pi * np.arange(count) / count
        tau = np.degrees(np.abs(np.angle(np.exp(1j * ph)))) - np.degrees(np.angle(W.L))
        rate = np.zeros(count)
        rate[tau > 0] = np.maximum(Lam(tau[tau > 0]), 0.0)
        rows = []
        for j in range(count):
            z = rho * (acb(0, 1) * (2 * arb.pi() * j / count)).exp()
            kernel = 1 / (target - f(z)) + m(z) / (target - z)
            rows.append((rho / z, kernel, arb(float(rate[j]))))
        quadrature[count] = rows
    record = {
        "status": "2400-bit numerical contour diagnostic, no interval quadrature certificate",
        "input": "w01_N*.npz: approximate corner vector sigma_1*Vt[0]; stored floats treated as dyadic coefficients",
        "contour": "fixed |z|=rho",
        "target_x_over_rho": str((target.real / rho).mid()),
        "kernel": "1/(t-f(z)) + m(z)/(t-z)",
        "weight": "exp(-N*Lambda_reference(angle)) on the fixed source circle",
        "cases": [],
    }
    for N in [int(x) for x in args.degrees.split(",")]:
        start = time.time()
        cache = np.load(project / "w_grenzschicht" / f"w01_N{N}.npz")
        assert np.array_equal(cache["coef"], first["coef"])
        sigma = arb(float(cache["S"][0]))
        y = [sigma * arb(float(v)) for v in cache["Vt"][0]]
        Y = acb_poly([acb(v) for v in y])
        B = sum((y[k] * rho**k * sum((m[j+k+1] * target**j for j in range(6-k)), acb(0)) for k in range(6)), acb(0))
        expected = sum((Y(rho / root) / (root * fp(root)) for root in far), acb(0)) - B
        approximations = []
        diagnostics = None
        for count, rows in quadrature.items():
            integral = acb(0)
            absolute = arb(0)
            weighted_y2 = arb(0)
            weighted_kernel2 = arb(0)
            for w, kernel, rate in rows:
                value = Y(w)
                integral += value * kernel
                absolute += abs(value) * abs(kernel)
                exponential = (N * rate).exp()
                weighted_y2 += (abs(value) / exponential)**2
                weighted_kernel2 += (abs(kernel) * exponential)**2
            integral /= count
            absolute /= count
            ny = (weighted_y2 / count).sqrt()
            nk = (weighted_kernel2 / count).sqrt()
            approximations.append(integral)
            if count == 4096:
                diagnostics = {
                    "absolute_integral": float(absolute.mid()),
                    "weighted_trace_L2": float(ny.mid()),
                    "weighted_kernel_dual_L2": float(nk.mid()),
                    "weighted_CS_bound": float((ny * nk).mid()),
                }
        defect = abs(approximations[-1] - expected)
        difference = abs(approximations[-1] - approximations[0])
        assert defect < arb("1e-65") and difference < arb("1e-65")
        case = {
            "N": N,
            "E_abs": float(abs(expected).mid()),
            "B_abs": float(abs(B).mid()),
            "far_root_contribution_abs": float(abs(expected + B).mid()),
            "quadrature_residue_defect": defect.str(8),
            "quadrature_doubling_difference": difference.str(8),
            **diagnostics,
            "seconds": time.time() - start,
        }
        case["absolute_bound_over_E"] = case["absolute_integral"] / case["E_abs"]
        record["cases"].append(case)
        (outdir / "h0_kernel_probe.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps(case, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
