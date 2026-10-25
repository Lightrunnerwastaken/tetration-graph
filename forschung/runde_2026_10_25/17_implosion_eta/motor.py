"""Eigener Zugang zum fatou.gp-Fork (nur lesend), eigener Zwischenspeicher in diesem Ordner, ohne geteilten state_cache.

werte(basis, ausdruecke, dps) -> Liste mpc.  basis: GP-Ausdruck als Text ("2", "exp(1)", "(3/2)+(1/10)*I").
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
CACHE = HIER / "cache" / "motor_cache.json"
_motoren = {}


def _cache():
    if CACHE.is_file():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    return {}


def _speichern(c):
    CACHE.parent.mkdir(exist_ok=True)
    tmp = CACHE.with_suffix(".tmp")
    tmp.write_text(json.dumps(c, indent=0), encoding="utf-8")
    tmp.replace(CACHE)


def motor(basis, dps):
    key = (basis, dps)
    if key not in _motoren:
        if str(GP_SRC) not in sys.path:
            sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=240.0, init_timeout=240.0, state_cache=False)
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


def werte(basis, ausdruecke, dps):
    c = _cache()
    keys = [f"{basis}|{dps}|{a}" for a in ausdruecke]
    fehlt = [i for i, k in enumerate(keys) if k not in c]
    if fehlt:
        s = motor(basis, dps)
        res = s.eval_batch([ausdruecke[i] for i in fehlt])
        for i, w in zip(fehlt, res):
            c[keys[i]] = [mp.nstr(mp.re(w), dps + 5), mp.nstr(mp.im(w), dps + 5)]
        _speichern(c)
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[k][0]), mp.mpf(c[k][1])) for k in keys]


def gpzahl(x, digits):
    x = mp.mpc(x)
    f = lambda v: mp.nstr(v, digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    if mp.im(x) == 0:
        return f"({f(mp.re(x))})"
    return f"({f(mp.re(x))}+({f(mp.im(x))})*I)"
