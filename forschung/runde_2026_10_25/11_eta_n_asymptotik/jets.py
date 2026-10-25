"""Taylor-Jets der kritischen Abbildung F_N(x)=H_N(eta_N,x) am Beruehrpunkt e_N.

Cauchy-Abtastung auf einem Kreis um e_N (M Punkte, Radius rho). Liefert a_k = F^(k)(e)/k!
fuer k=0..K und die Invarianten beta=a_2, iteratives Residuum iota = 1 - a_3/a_2^2.
Nutzt die eingefrorene Engine-Kopie unter runtime/ (keine geteilten Zwischenspeicher).
Aufruf: python -B jets.py --rank 6 --dps 25 --M 32 --rho 0.1 --max-seconds 280
"""
import argparse, csv, json, sys, time
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE / "runtime"))
from common import engine, mp, stringify

p = argparse.ArgumentParser()
p.add_argument("--rank", type=int, required=True)
p.add_argument("--dps", type=int, default=25)
p.add_argument("--M", type=int, default=32)
p.add_argument("--rho", type=str, default="0.1")
p.add_argument("--K", type=int, default=8)
p.add_argument("--max-seconds", type=int, default=280)
a = p.parse_args()
mp.mp.dps = a.dps + 20
with (HERE / "eta_values.csv").open(encoding="utf-8-sig", newline="") as f:
    row = next(z for z in csv.DictReader(f) if int(z["map_rank"]) == a.rank)
b, e = row["eta"], row["e"]
E = mp.mpf(e); rho = mp.mpf(a.rho)
start = time.monotonic()
out = {"rank": a.rank, "dps": a.dps, "M": a.M, "rho": a.rho, "base": b, "e": e, "status": "running"}
target = HERE / "cache" / f"jet_r{a.rank}_d{a.dps}_M{a.M}_rho{a.rho}.json"
gp = engine(a.dps, HERE / "runtime" / "hyper.gp", logged=False, timeout=a.max_seconds - 20)
try:
    gp.eval_batch(b, [f"hrmaxdigits[7]={a.dps-6}", f"hinit({a.rank})"])
    print("init", time.monotonic() - start, flush=True)
    zs = [E + rho * mp.expjpi(2 * mp.mpf(j) / a.M) for j in range(a.M)]
    exprs = [f"hval({a.rank},{mp.nstr(z.real, a.dps+8)}+({mp.nstr(z.imag, a.dps+8)})*I)" for z in zs]
    vals = []
    for k in range(0, a.M, 4):
        if time.monotonic() - start > a.max_seconds - 15:
            raise TimeoutError("budget")
        vals += gp.eval_batch(b, exprs[k:k+4])
    coef = []
    for k in range(a.K + 1):
        s = mp.fsum(mp.mpc(vals[j]) * mp.expjpi(-2 * mp.mpf(j) * k / a.M) for j in range(a.M)) / a.M
        coef.append(s / rho**k)
    out["coef"] = [stringify(c, a.dps) for c in coef]
    a2 = coef[2].real; a3 = coef[3].real
    out["beta"] = mp.nstr(a2, a.dps); out["iota"] = mp.nstr(1 - a3 / a2**2, a.dps)
    out["status"] = "completed"
    print("coef", [mp.nstr(c.real, 15) for c in coef], flush=True)
    print("beta", out["beta"], "iota", out["iota"], flush=True)
except Exception as exc:
    out["status"] = "failed"; out["error"] = f"{type(exc).__name__}: {exc}"; print(out["error"], flush=True)
finally:
    out["seconds"] = time.monotonic() - start
    gp.close()
    target.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("seconds", out["seconds"], flush=True)
