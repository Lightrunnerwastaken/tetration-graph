"""Reproduce the forward and inverse quadratures in the accompanying note.

Requires mpmath, the gp-tetration fatou_backend package, and PARI/GP.
If fatou_backend is not installed, use --backend-source PATH_TO_GP_TETRATION/src.
The program prints JSON. It disables the backend's on-disk state cache.
"""

import argparse
import json
import sys

import mpmath as mp


def real_batch(function, base, points):
    result = []
    for start in range(0, len(points), 64):
        values = function(base, points[start:start + 64])
        if any(abs(mp.im(value)) > mp.mpf("1e-280") for value in values):
            raise ArithmeticError("A real Kneser evaluation produced a nonreal value")
        result.extend(map(mp.re, values))
    return result


def calculate_forward(backend, count, modes):
    mp.mp.dps = 345
    log_two = mp.log(2)
    p = mp.log(log_two)
    grid = [mp.mpf(j) / count for j in range(count)]
    with backend(dps=335, fatou_gp="fork", state_cache=False, init_timeout=1500) as gp:
        values = real_batch(gp.sexp_batch, "2", [2 + t for t in grid])
        top = [mp.exp(log_two * value + p) + p for value in values]
        reduced = [mp.log(v + mp.log1p(p * mp.exp(-v))) for v in top]
        heights = real_batch(gp.slog_batch, "e", reduced)
        phase = [value - 2 - t for value, t in zip(heights, grid)]

        step = mp.mpf("1e-42")
        probes = [(mp.log(k - p) + sign * step - p) / log_two
                  for k in modes for sign in (-1, 1)]
        measured = real_batch(gp.slog_batch, "2", probes)
        scales = [(right - left) / (2 * step)
                  for left, right in zip(measured[::2], measured[1::2])]

    mean = mp.fsum(phase) / count
    coefficients = []
    for k in modes:
        terms = ((value - mean) * mp.expjpi(-2 * k * t)
                 for value, t in zip(phase, grid))
        coefficients.append(mp.fsum(terms) / count)
    return mean, coefficients, scales


def calculate_inverse(backend, count, modes):
    mp.mp.dps = 335
    log_two = mp.log(2)
    q = mp.log(log_two)
    grid = [(mp.mpf(j) + mp.mpf("0.5")) / count for j in range(count)]
    with backend(dps=325, fatou_gp="fork", state_cache=False, init_timeout=1500) as gp:
        values = real_batch(gp.sexp_batch, "e", [1 + t for t in grid])
        top = [mp.exp(value) for value in values]
        reduced = [(mp.log(v + mp.log1p(-q * mp.exp(-v)) - q) - q) / log_two
                   for v in top]
        heights = real_batch(gp.slog_batch, "2", reduced)
        inverse_lift = [height - 1 for height in heights]
    mean = mp.fsum(t - value for t, value in zip(grid, inverse_lift)) / count
    coefficients = []
    for k in modes:
        integral = mp.fsum(mp.expjpi(-2 * k * value) for value in inverse_lift) / count
        coefficients.append(integral / (2 * mp.pi * mp.j * k))
    return mean, coefficients


def finite_endpoints(backend):
    """Upper endpoint from the logarithmic continuation, at height x = 0."""
    mp.mp.dps = 100
    log_two = mp.log(2)
    p = mp.log(log_two)
    step = mp.mpf("1e-22")
    results = []
    with backend(dps=85, fatou_gp="fork", state_cache=False) as gp:
        minus, plus = gp.sexp_batch("2", [-step, step])
        base_derivative = log_two * mp.re(plus - minus) / (2 * step)
        for level in (6, 7, 8):
            depth = level - 4
            # Selecting this integer at level 8 needs nearly 20,000 digits.
            with mp.workdps(20000 if level == 8 else 100):
                source_log = mp.log(2)
                source_offset = mp.log(source_log)
                w = source_log + source_offset
                multiplier = mp.mpf(1)
                for _ in range(depth):
                    multiplier *= mp.exp(w)
                    w = mp.exp(w) + source_offset
                angle = 2 * mp.pi * mp.nint(mp.exp(w) / (2 * mp.pi))
                endpoint = mp.log(mp.log(-source_offset) - source_offset + 1j * angle)
            source_value = +endpoint
            for _ in range(depth):
                source_value = mp.log(source_value - p)
            point = gp.slog_batch("2", [(source_value - p) / log_two])[0]
            results.append({
                "level": level,
                "real": mp.nstr(point.real, 18),
                "imag": mp.nstr(point.imag, 18),
                "scaled_imag": mp.nstr(point.imag * base_derivative * multiplier, 18),
            })
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend-source", help="directory containing fatou_backend")
    args = parser.parse_args()
    if args.backend_source:
        sys.path.insert(0, args.backend_source)
    from fatou_backend import FatouGP

    modes = [10, 20, 50, 100, 200, 350, 500, 600, 700]
    count = 2048
    mean, forward, scales = calculate_forward(FatouGP, count, modes)
    reverse_mean, reverse = calculate_inverse(FatouGP, count, modes)
    rows = []
    for k, a, other, scale in zip(modes, forward, reverse, scales):
        rows.append({
            "k": k,
            "log10_abs": mp.nstr(mp.log10(abs(a)), 24),
            "relative_complex_difference": mp.nstr(abs(a - other) / abs(a), 6),
            "r_d": mp.nstr(scale, 24),
            "damping_over_k_r": mp.nstr(-mp.log(abs(a)) / (k * scale), 24),
        })
    report = {
        "sample_count_each": count,
        "forward_gp_digits": 335,
        "inverse_gp_digits": 325,
        "forward_mean": mp.nstr(mean, 210),
        "inverse_mean": mp.nstr(reverse_mean, 210),
        "mean_difference": mp.nstr(abs(mean - reverse_mean), 6),
        "modes": rows,
    }
    report["finite_endpoints"] = finite_endpoints(FatouGP)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

