"""Eigener Zugang zum fatou.gp-Fork (nur lesend), Zwischenspeicher NUR in diesem Ordner (cache/).

werte(basis, ausdruecke, dps) -> Liste mpc. Ausdrücke sind GP-Ausdrücke, z. B. "sexp(0.5)", "pentx(0.5)".
Der Zustands-Zwischenspeicher des Motors (FATOU_CACHE_DIR) liegt ebenfalls in cache/gpstate.
"""
import json
import os
import sys
import time
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
CDIR = HIER / "cache"
os.environ["FATOU_CACHE_DIR"] = str(CDIR / "gpstate")
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
_motoren = {}


def _cfile(basis, dps):
    tag = basis.replace("(", "").replace(")", "").replace("/", "_").replace("*", "x")
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
        _motoren[key] = (m, s)
    return _motoren[key][1]


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


def schliessen():
    for m, s in _motoren.values():
        try:
            m.close()
        except Exception:
            pass
    _motoren.clear()
