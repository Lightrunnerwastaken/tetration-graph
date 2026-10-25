"""Eigener Motor-Zugang (Schnittstelle wie runde_2026_10_18/02, aber state_cache=False wie in 03).
Zwischenspeicher nur in cache/ dieses Ordners.  gp-Fork wird nur gelesen."""
import json, sys, time
from pathlib import Path
import mpmath as mp
HIER = Path(__file__).resolve().parent
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
CDIR = HIER / "cache"
_motoren = {}

def _cfile(basis, dps):
    tag = "e" if basis == "exp(1)" else basis.replace("(", "").replace(")", "")
    return CDIR / f"{tag}_{dps}.json"

def motor(basis, dps):
    key = (basis, dps)
    if key not in _motoren:
        if str(GP_SRC) not in sys.path:
            sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=280.0, init_timeout=280.0, state_cache=False)
        t = time.time()
        s = m.session(basis)
        s.eval_batch(["1"])
        print(f"[motor] {basis} dps={dps} bereit nach {time.time()-t:.1f}s", flush=True)
        _motoren[key] = (m, s)
    return _motoren[key][1]

def schliessen():
    for m, s in _motoren.values():
        try:
            m.close()
        except Exception:
            pass
    _motoren.clear()

def gpzahl(x, digits):
    x = mp.mpc(x)
    f = lambda v: mp.nstr(v, digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    if mp.im(x) == 0:
        return f"({f(mp.re(x))})"
    return f"({f(mp.re(x))}+({f(mp.im(x))})*I)"

def werte(basis, ausdruecke, dps, chunk=200, speichern=True):
    CDIR.mkdir(exist_ok=True)
    cf = _cfile(basis, dps)
    c = json.loads(cf.read_text(encoding="utf-8")) if cf.is_file() else {}
    fehlt = [a for a in dict.fromkeys(ausdruecke) if a not in c]
    for i in range(0, len(fehlt), chunk):
        teil = fehlt[i:i + chunk]
        res = motor(basis, dps).eval_batch(teil)
        for a, w in zip(teil, res):
            c[a] = [mp.nstr(mp.re(w), dps + 10), mp.nstr(mp.im(w), dps + 10)]
        if speichern:
            cf.write_text(json.dumps(c, indent=0), encoding="utf-8")
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[a][0]), mp.mpf(c[a][1])) for a in ausdruecke]

def roh(basis, ausdruecke, dps):
    """ohne Zwischenspeicher"""
    return motor(basis, dps).eval_batch(ausdruecke)
