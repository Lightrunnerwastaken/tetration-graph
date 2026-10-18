"""Katalog der getesteten Störungen φ = x + p(x)."""
import json
from pathlib import Path
import mpmath as mp
E4 = Path(__file__).resolve().parent.parent.parent / "identitaeten_2026_10_14" / "etappe4" / "daten"

def gegenbeispiel(f):
    d = json.loads((E4 / f"gegenbeispiel_{f}_40.json").read_text(encoding="utf-8"))
    return ("trig", [(int(k), mp.mpf(a), mp.mpf(b)) for k, a, b in d["p"]])

def katalog(name):
    if name.startswith("sin"):        # sin:eps
        return ("trig", [(1, mp.mpf(name.split(":")[1]), mp.mpf(0))])
    if name.startswith("cos"):        # cos:eps   eps (cos 2πx − 1)
        return ("trig", [(1, mp.mpf(0), mp.mpf(name.split(":")[1]))])
    if name.startswith("zz_unbenutzt"):
        return ("trig", [(2, mp.mpf(name.split(":")[1]), mp.mpf(0))])
    if name.startswith("moeb"):       # moeb:re:im
        _, re, im = name.split(":")
        return ("moeb", mp.mpc(mp.mpf(re), mp.mpf(im)))
    if name.startswith("gb"):         # gb:c
        return gegenbeispiel(name.split(":")[1])
    if name == "kneser":
        return ("kneser", None)
    raise ValueError(name)

def p_komplex(stoer, z):
    """p(z) und p'(z) für komplexes z (für die Singularitätensuche)."""
    art, par = stoer
    if art == "kneser":
        return mp.mpc(0), mp.mpc(0)
    if art == "trig":
        p = mp.mpc(0); dp = mp.mpc(0)
        for k, a, b in par:
            w = 2 * mp.pi * k
            p += a * mp.sin(w * z) + b * (mp.cos(w * z) - 1)
            dp += w * (a * mp.cos(w * z) - b * mp.sin(w * z))
        return p, dp
    if art == "moeb":
        a = mp.mpc(par); ab = mp.conj(a)
        q = mp.exp(2j * mp.pi * z); qb = mp.exp(-2j * mp.pi * z)  # conj(q) auf R, analytisch fortgesetzt
        # arg(1 − ā q) = [log(1 − ā q) − log(1 − a q̄)]/(2i)
        A = (mp.log(1 - ab * q) - mp.log(1 - a * qb)) / 2j
        A0 = mp.im(mp.log(1 - ab))
        p = -(A - A0) / mp.pi
        dA = (-ab * 2j * mp.pi * q / (1 - ab * q) - a * 2j * mp.pi * qb / (1 - a * qb)) / 2j
        return p, -dA / mp.pi
    raise ValueError(art)
