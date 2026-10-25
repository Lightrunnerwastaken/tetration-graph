"""Rohdaten für die Zertifizierung: Motorwerte von sexp auf den Cauchy-Kreisen |z - c| = r.

Gelesen werden nur die Zwischenspeicher von Etappe 3 (etappe3/cache/<basis>_<dps>.json, Schlüssel "sexp(<z>)").
Zu jedem Zentrum c werden die M Knoten z_j = c + r e^{2πij/M} gesucht; der gespeicherte Knoten z~_j ist die
Dezimaldarstellung, die der Motor tatsächlich bekommen hat. Rückgabe als Dezimalzeichenketten (keine Rundung hier).
"""
import json
import re
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E3 = HIER
CACHE = E3 / "cache"
_RX = re.compile(r"^sexp\(\(([^()]+?)(?:\+\(([^()]+)\)\*I)?\)\)$")
_geladen = {}


def _lade(tag, dps):
    key = (tag, dps)
    if key not in _geladen:
        d = json.loads((CACHE / f"{tag}_{dps}.json").read_text(encoding="utf-8"))
        eintraege = []
        for k, v in d.items():
            m = _RX.match(k)
            if not m:
                continue
            eintraege.append((m.group(1), m.group(2) or "0", v[0], v[1]))
        _geladen[key] = eintraege
    return _geladen[key]


def kreis(tag, dps, c, r, M):
    """Liste von M Tupeln (re z~, im z~, re s, im s) als Strings, Index j = Knotennummer."""
    with mp.workdps(40):
        c = mp.mpf(c); r = mp.mpf(r)
        out = [None] * M
        for zr, zi, sr, si in _lade(tag, dps):
            z = mp.mpc(mp.mpf(zr), mp.mpf(zi))
            if abs(abs(z - c) - r) > mp.mpf("1e-30"):
                continue
            t = mp.arg(z - c) / (2 * mp.pi) * M
            j = int(mp.nint(t)) % M
            if abs(t - mp.nint(t)) > mp.mpf("1e-20"):
                continue
            out[j] = (zr, zi, sr, si)
    fehlend = [j for j in range(M) if out[j] is None]
    if fehlend:
        raise ValueError(f"{tag}_{dps} c={c} r={r} M={M}: {len(fehlend)} Knoten fehlen, z.B. {fehlend[:5]}")
    return out


if __name__ == "__main__":
    for tag, dps, cs, r, M in [("e", 130, ["-0.4", "-0.2", "0", "0.2", "0.4"], "0.2", 256),
                               ("e", 100, ["-0.4", "-0.2", "0", "0.2", "0.4"], "0.2", 128),
                               ("e", 70, ["-0.4", "-0.2", "0", "0.2", "0.4"], "0.2", 256),
                               ("e", 80, ["-0.6", "-0.8", "-1.0", "-1.2", "-1.4", "-1.6"], "0.2", 256),
                               ("2", 80, ["-0.4", "-0.2", "0", "0.2", "0.4"], "0.2", 256)]:
        for c in cs:
            try:
                k = kreis(tag, dps, c, r, M)
                print(tag, dps, c, "ok", k[0][2][:20])
            except ValueError as e:
                print(e)
