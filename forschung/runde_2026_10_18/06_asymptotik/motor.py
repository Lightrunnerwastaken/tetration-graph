"""Runde 18.10.2026, 06_asymptotik: Zugang zum Kneser-Motor (fatou.gp-Fork) mit EIGENEM Zwischenspeicher cache/.

Vorlage: forschung/identitaeten_2026_10_14/etappe2/kern.py (nur kopiert, nicht verändert).
"""
import json
import sys
import time
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
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=550.0, init_timeout=550.0)
        t = time.time()
        s = m.session(basis)
        s.eval_batch(["1"])
        print(f"[motor] {basis} dps={dps} bereit nach {time.time()-t:.1f}s", flush=True)
        _motoren[key] = s
    return _motoren[key]


def gpzahl(x, digits):
    x = mp.mpc(x)
    f = lambda v: mp.nstr(v, digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    if mp.im(x) == 0:
        return f"({f(mp.re(x))})"
    return f"({f(mp.re(x))}+({f(mp.im(x))})*I)"


def werte(basis, ausdruecke, dps, chunk=100):
    CDIR.mkdir(exist_ok=True)
    cf = _cfile(basis, dps)
    c = json.loads(cf.read_text(encoding="utf-8")) if cf.is_file() else {}
    fehlt = [a for a in dict.fromkeys(ausdruecke) if a not in c]
    for i in range(0, len(fehlt), chunk):
        teil = fehlt[i:i + chunk]
        res = motor(basis, dps).eval_batch(teil)
        for a, w in zip(teil, res):
            c[a] = [mp.nstr(mp.re(w), dps + 10), mp.nstr(mp.im(w), dps + 10)]
        cf.write_text(json.dumps(c, indent=0), encoding="utf-8")
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[a][0]), mp.mpf(c[a][1])) for a in ausdruecke]


def sexp(zs, dps, basis="exp(1)"):
    return werte(basis, [f"sexp({gpzahl(z, dps + 10)})" for z in zs], dps)


def slog(ws, dps, basis="exp(1)"):
    return werte(basis, [f"slog({gpzahl(w, dps + 10)})" for w in ws], dps)
