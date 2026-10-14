"""Eigener dünner Zugang zum fatou.gp-Fork (gp-tetration), mit eigenem Zwischenspeicher in diesem Ordner.

Bewusst NICHT experimente/gemeinsam.kneser: dessen 190-MB-Cache wird von anderen Sitzungen geteilt.
Aufruf: werte(basis, ausdruecke, dps) -> Liste mpc; Ausdrücke sind GP-Ausdrücke wie "sexp(0.5)".
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
CACHE = HIER / "daten" / "motor_cache.json"
_motoren = {}


def _cache():
    if CACHE.is_file():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    return {}


def _speichern(c):
    CACHE.parent.mkdir(exist_ok=True)
    CACHE.write_text(json.dumps(c, indent=0), encoding="utf-8")


def motor(basis, dps):
    key = (basis, dps)
    if key not in _motoren:
        if str(GP_SRC) not in sys.path:
            sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=900.0, init_timeout=1100.0)
        t = time.time()
        s = m.session(basis)
        s.eval_batch(["1"])
        print(f"[motor] {basis} dps={dps} bereit nach {time.time()-t:.1f}s", flush=True)
        _motoren[key] = s
    return _motoren[key]


def werte(basis, ausdruecke, dps):
    c = _cache()
    keys = [f"{basis}|{dps}|{a}" for a in ausdruecke]
    fehlt = [i for i, k in enumerate(keys) if k not in c]
    if fehlt:
        s = motor(basis, dps)
        res = s.eval_batch([ausdruecke[i] for i in fehlt])
        for i, w in zip(fehlt, res):
            c[keys[i]] = [mp.nstr(mp.re(w), dps + 10), mp.nstr(mp.im(w), dps + 10)]
        _speichern(c)
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[k][0]), mp.mpf(c[k][1])) for k in keys]


def gpzahl(x, digits):
    """mpmath-Zahl als GP-Literal."""
    x = mp.mpc(x)
    re_, im_ = mp.nstr(mp.re(x), digits, min_fixed=-mp.inf, max_fixed=mp.inf), mp.nstr(mp.im(x), digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    if mp.im(x) == 0:
        return f"({re_})"
    return f"({re_}+({im_})*I)"
