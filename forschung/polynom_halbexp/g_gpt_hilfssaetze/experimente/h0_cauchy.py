"""H0 round 2: high-precision, SVD-free Cauchy-kernel probes.

Numerical diagnostics, not interval certificates. Reads the other strands;
writes only to --output-dir. No K2 calculation and no nonlinear iteration.
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
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--model", choices=["kubisch", "kneser_z"], default="kubisch")
    parser.add_argument("--degrees", default="200,300,400,600,800,1000")
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    project = (args.project_root or package.parent).resolve()
    outdir = (args.output_dir or package / "ergebnisse").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(project / "o_operator"))
    sys.path.insert(0, str(project / "w_grenzschicht"))
    import o_kern as K
    from w40_lade import lade

    degrees = [int(x) for x in args.degrees.split(",")]
    if args.model == "kneser_z" and max(degrees) > 639:
        raise ValueError("Certified input ends at 640; J_N needs coefficients through N+1.")
    rho, L = K.RHO, K.Lc
    theta = np.angle(L)
    beta = 1 / 1.1182  # The same rounded exponent as O38, for comparison.
    tau_values = [0.25, 0.5, 1, 2, 4, 8, 15, 30, 60]

    # Fixed gamma* sample, constructed from the cubic and W's reference model.
    # This is a numerical contour construction, not a geometry certificate.
    ctx.prec = 300
    cubic = [float(x) for x in K.kubik_koeff()]
    a = cubic[0]
    R0 = abs(L - a)
    phi0 = np.angle(L - a)
    gamma = []
    previous = L
    for tau in tau_values + [float(np.degrees(phi0))]:
        eta = a + R0 * np.exp(1j * (phi0 - np.deg2rad(tau)))
        coefficients = cubic.copy()
        coefficients[0] -= eta
        roots = np.roots(list(reversed(coefficients)))
        root = roots[np.argmin(abs(roots - previous))]
        previous = root
        gamma.append((f"bay_{tau:.8g}", complex(root)))
    H, Lam, _, _ = lade()
    ph = np.linspace(-np.pi, np.pi, 8192, endpoint=False)
    taus = np.degrees(abs(ph)) - np.degrees(theta)
    lv = np.zeros_like(ph)
    lv[taus > 0] = np.maximum(Lam(taus[taus > 0]), 0.0)

    def poisson(s, angle):
        radius = np.exp(s)
        kernel = (radius * radius - 1) / (radius * radius - 2 * radius * np.cos(angle - ph) + 1)
        return float(np.mean(lv * kernel))

    for tau in tau_values + [180 - float(np.degrees(theta))]:
        angle = theta + np.deg2rad(tau)
        lo, hi = 0.0, 0.3
        if poisson(1e-7, angle) <= 0.7e-7:
            hi = 0.0
        else:
            for _ in range(45):
                mid = (lo + hi) / 2
                if poisson(mid, angle) > 0.7 * mid:
                    lo = mid
                else:
                    hi = mid
        gamma.append((f"lens_{tau:.8g}", complex(rho * np.exp(hi + 1j * angle))))

    previous_results = json.loads((project / "o_operator" / f"o38_{args.model}.json").read_text(encoding="utf-8"))
    record = {
        "status": "high-precision numerical solves (algorithm=approx), not interval certificates",
        "model": args.model,
        "model_source": "o_operator/o_kern.py; the cubic uses its stored decimal fixed-point constants",
        "source_geometry": "fixed cubic gamma* samples from w40_modell_ref.npz; the same points are used for both models",
        "kb_convention": "Y/G on the inner half of the circle, T/GT on the outer half, physical x units",
        "beta_used": beta,
        "reference_rate": float(H([-1])[0]),
        "cases": [],
    }
    for N in degrees:
        start = time.time()
        bits = int(1.6 * N) + 400
        J, _, _ = K.matrizen(args.model, N, bits)
        ctx.cap = N + 12
        fk, mk, _ = K.modell(args.model, N)
        fpol, mpol = acb_poly([acb(c) for c in fk]), acb_poly([acb(c) for c in mk])
        rh = (arb(K.L_RE) ** 2 + arb(K.L_IM) ** 2).sqrt()
        rpowers = [rh ** k for k in range(N + 1)]
        points = [("corner", L)]
        for ps in range(-180, 180, 10):
            points.append((f"kb_{ps}", complex(L + rho * N ** (-beta) * np.exp(1j * (theta + np.deg2rad(ps))))))
        points.extend(gamma)
        columns = []
        metadata = []
        for name, x in points:
            xc = acb(float(x.real), float(x.imag))
            vector = [(rh / xc) ** k for k in range(N + 1)]
            columns += [[v.real for v in vector], [v.imag for v in vector]]
            metadata.append((name, x, len(columns) - 2, "G", None))
            if name.startswith("kb_") or name == "corner":
                s1 = list(acb_series((acb_poly([xc]) - fpol).coeffs(), prec=N + 1).inv().coeffs())
                s1 += [acb(0)] * (N + 1 - len(s1))
                delta = [xc ** (-k - 1) * rpowers[k] for k in range(N + 1)]
                kz = acb_poly([xc ** (-k - 1) for k in range(N + 1)])
                s2 = list((mpol * kz).truncate(N + 1).coeffs())
                s2 += [acb(0)] * (N + 1 - len(s2))
                vector = [(s1[k] + s2[k]) * rpowers[k] for k in range(N + 1)]
                columns += [[v.real for v in vector], [v.imag for v in vector]]
                metadata.append((name, x, len(columns) - 2, "GT", delta))
        low_start = len(columns)
        for k in range(6):
            columns.append([arb(int(j == k)) for j in range(N + 1)])
        rhs = arb_mat(N + 1, len(columns), [columns[j][k] for k in range(N + 1) for j in range(len(columns))])
        solution = J.solve(rhs, algorithm="approx")
        print(f"{args.model} N={N}: solved {len(columns)} real right-hand sides at {bits} bits ({time.time()-start:.1f}s)", flush=True)
        point_records = {}
        for name, x, col, field, delta in metadata:
            values = [acb(solution[k, col], solution[k, col + 1]) for k in range(N + 1)]
            if delta is not None:
                values = [values[k] - delta[k] for k in range(N + 1)]
            squared = sum((v.real ** 2 + v.imag ** 2 for v in values), arb(0))
            norm = squared.sqrt()
            assert norm.is_finite()
            item = point_records.setdefault(name, {"x": [x.real, x.imag], "outside_D": bool(abs(x) > rho)})
            item[field] = float(norm.mid())
        low_norms = []
        for j in range(6):
            squared = sum((solution[k, low_start+j] ** 2 for k in range(N+1)), arb(0))
            low_norms.append(float(squared.sqrt().mid()))
        case = {"N": N, "bits": bits, "points": point_records, "inverse_first_six_column_norms": low_norms}
        circle = [v for k,v in point_records.items() if k.startswith("kb_")]
        case["KB_mixed_over_sqrtN"] = max(v["GT"] if v["outside_D"] else v["G"] for v in circle) / math.sqrt(N)
        case["G_whole_circle_over_sqrtN"] = max(v["G"] for v in circle) / math.sqrt(N)
        case["GT_whole_circle_over_sqrtN"] = max(v["GT"] for v in circle) / math.sqrt(N)
        source = previous_results[str(N)]
        source_points = source["kreis"]["1"]
        reference_kb = max(math.exp(v["GT"] if v["aussen"] else v["G"]) for v in source_points) / math.sqrt(N)
        case["reference_O38_mixed_over_sqrtN"] = reference_kb
        case["relative_difference_from_O38"] = case["KB_mixed_over_sqrtN"] / reference_kb - 1
        case["reference_inverse_top_singular_value"] = source["s"][0]
        g_gamma = [v["G"] for k,v in point_records.items() if k.startswith(("bay_", "lens_"))]
        case["gamma_sample_G_over_sqrtN"] = max(g_gamma) / math.sqrt(N)
        if args.model == "kubisch":
            # E = I_F-B = -B on a large contour enclosing all three roots.
            # B has degree 5 and uses only y_0,...,y_5.
            mm = [mk[j] if j < len(mk) else arb(0) for j in range(7)]
            tmax = arb("1.5") * rh
            bound = sum((arb(low_norms[k]) * rpowers[k] *
                         sum((abs(mm[j+k+1]) * tmax**j for j in range(6-k)), arb(0))
                         for k in range(6)), arb(0))
            case["B_majorant_from_column_norms_for_abs_x_le_1p5rho"] = float(bound.mid())
            outer_radius = 16 * rh
            rouche_margin = abs(fk[3]) * outer_radius**3 - sum((abs(fk[j]) * outer_radius**j for j in range(3)), arb(0)) - tmax
            assert rouche_margin > 0
            case["large_contour_radius_over_rho"] = 16
            case["rouche_margin_lower"] = str(rouche_margin.lower())
        case["seconds"] = time.time() - start
        record["cases"].append(case)
        (outdir / f"h0_cauchy_{args.model}.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"{args.model} N={N}: KB/sqrtN={case['KB_mixed_over_sqrtN']:.9g}; Gwhole/sqrtN={case['G_whole_circle_over_sqrtN']:.9g}; GTwhole/sqrtN={case['GT_whole_circle_over_sqrtN']:.9g}; delta_O38={case['relative_difference_from_O38']:.3g}", flush=True)
    if len(record["cases"]) >= 4:
        ns = np.array([x["N"] for x in record["cases"]], dtype=float)
        logs = np.log([x["reference_inverse_top_singular_value"] for x in record["cases"]])
        coefficients = np.linalg.lstsq(np.column_stack([ns, np.log(ns), np.ones(len(ns))]), logs, rcond=None)[0]
        record["rate_fit_of_existing_O38_data"] = {"c": float(coefficients[0]), "logN_coefficient": float(coefficients[1]), "constant": float(coefficients[2]), "note": "Fit of existing O38 singular values, not a newly computed SVD or an asymptotic proof"}
        (outdir / f"h0_cauchy_{args.model}.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
