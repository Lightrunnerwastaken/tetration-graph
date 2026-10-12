"""Lightweight audit of supplied rank-asymptotic proposals; no new eta solve.

All critical transformations use old fold data and are consequently not
independent evidence for a limit. Gamma shift checks concern a linear model.
"""
import csv
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).absolute().parent
SOURCE = HERE / "inputs/eta_values.csv"
mp.mp.dps = 90

def bfun(nu, x):
    if x == 0:
        return mp.mpf(0)
    return mp.gamma(x + nu) * mp.rgamma(x) / mp.gamma(nu + 1)

def fmt(x):
    return mp.nstr(x, 45)

def evaluate(dps):
    mp.mp.dps = dps
    rows = []
    with SOURCE.open(encoding="utf-8-sig", newline="") as f:
        raw = list(csv.DictReader(f))
    previous = None
    for row in raw:
        r = int(row["map_rank"])
        m = r - 2  # number of forced successor zeros for H_3=b**x indexing
        eta, e = mp.mpf(row["eta"]), mp.mpf(row["e"])
        ac = -1 / bfun(m, e)
        uc = m**(e - 1) * ac
        slopes = {}
        if previous:
            rp, ep, xp = previous
            for index, index_old, label in [(r, rp, "rank_r"), (m, rp-2, "zero_count_m")]:
                slopes[f"eta_exponent_{label}"] = mp.log((2-ep)/(2-eta)) / mp.log(mp.mpf(index)/index_old)
                slopes[f"e_exponent_{label}"] = mp.log((4-xp)/(4-e)) / mp.log(mp.mpf(index)/index_old)
        rows.append({
            "rank": r, "m": m, "eta": eta, "e": e,
            "A_at_fold": ac, "m_cubed_A_at_fold": m**3*ac,
            "U_at_fold": uc, "U_over_minus_gamma_e": uc/(-mp.gamma(e)),
            "U_log_derivative_at_fold": mp.log(m)-mp.digamma(m+e)+mp.digamma(e),
            "e_gap_times_log_m": (4-e)*mp.log(m), **slopes,
            "unshifted_r_cubed_A_at_fold": -r**3/bfun(r,e),
        })
        previous = (r, eta, e)
    checks = []
    xs = [mp.mpf("-0.4"), mp.mpf("0.25"), mp.mpf("0.7"), mp.mpf("1.5"), mp.mpc("0.3", "0.2")]
    for nu in [mp.mpf(1), mp.mpf("1.5"), mp.mpf(2), mp.mpf("2.5"), mp.mpf("5.25")]:
        for x in xs:
            lhs = bfun(nu+1,x+1)-bfun(nu+1,x)
            rhs = bfun(nu,x+1)
            checks.append(abs(lhs-rhs)/(1+abs(rhs)))
    # v=x(x-1)=2 B_2-2 B_1; a fractional linear rank step sends B_nu to B_(nu+s).
    toy = []
    for j in range(81):
        x=mp.mpf(j)/40
        vals=[2*(bfun(2+s,x)-bfun(1+s,x)) for s in [mp.mpf(0),mp.mpf(".5"),mp.mpf(1)]]
        toy.append({"x":x, "v":vals[0], "half_step":vals[1], "full_step":vals[2]})
    return rows, max(checks), toy

def fractional_sum_checks():
    """Actual semigroup on finite height grids, not just interpolated numbers."""
    def kernel(s,k):
        return mp.rf(s,k)/mp.factorial(k)
    def matrix(s,n):
        return mp.matrix(n,n,lambda i,j:kernel(s,i-j) if i>=j else 0)
    half=matrix(mp.mpf(".5"),16)
    full=matrix(mp.mpf(1),16)
    semigroup=half*half-full
    max_error=max(abs(z) for z in semigroup)
    gamma_errors=[]
    for nu in [mp.mpf(1),mp.mpf(2),mp.mpf("2.5"),mp.mpf("5.25")]:
        for k in range(1,17):
            value=sum(kernel(mp.mpf(".5"),k-j)*bfun(nu,mp.mpf(j)) for j in range(1,k+1))
            gamma_errors.append(abs(value-bfun(nu+mp.mpf(".5"),mp.mpf(k)))/(1+abs(value)))
    return max_error,max(gamma_errors)

def main():
    start=time.monotonic()
    low, err_low, toy = evaluate(90)
    high, err_high, _ = evaluate(170)
    semigroup_error, sum_gamma_error=fractional_sum_checks()
    control=max(abs(x[k]-y[k])/(1+abs(y[k])) for x,y in zip(low, high) for k in x)
    with (HERE / "critical_gamma_audit.csv").open("w", encoding="utf-8", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(high[-1]));writer.writeheader()
        for row in high: writer.writerow({k:fmt(v) if isinstance(v,mp.mpf) else v for k,v in row.items()})
    with (HERE / "linear_half_step.csv").open("w", encoding="utf-8", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(toy[0]));writer.writeheader()
        writer.writerows({k:fmt(v) for k,v in row.items()} for row in toy)
    out={"source":str(SOURCE),"source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         "precision_dps":[90,170],"max_relative_control_difference":fmt(control),
         "max_gamma_difference_identity_error_90dps":fmt(err_low),
         "max_gamma_difference_identity_error_170dps":fmt(err_high),
         "max_16_by_16_half_sum_squared_minus_full_sum_error":fmt(semigroup_error),
         "max_relative_fractional_sum_vs_gamma_shift_error":fmt(sum_gamma_error),
         "critical_table_is_only_algebraic_reexpression_of_existing_fold_data":True,
         "half_step_is_linearized_model_not_an_exact_fractional_hyperoperation":True,
         "seconds":time.monotonic()-start}
    (HERE / "audit_metadata.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
    for row in high:
        print(row["rank"],"m^3 A",mp.nstr(row["m_cubed_A_at_fold"],12),"U",mp.nstr(row["U_at_fold"],12),
              "eta local exponent",mp.nstr(row.get("eta_exponent_rank_r",mp.nan),8),
              "e local exponent",mp.nstr(row.get("e_exponent_rank_r",mp.nan),8))

if __name__ == "__main__":
    main()
