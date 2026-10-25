"""Schneller Kneser-sexp-Auswerter (Basis e) fuer Bilder: Taylor-Reihen um mehrere Zentren 0.5+i*yc,
Koeffizienten aus Cauchy-Abtastung mit dem gp-Motor (motor.py, eigener Zwischenspeicher in cache/).

feld(z)  -> F(z) fuer 0<=Re z<=1, |Im z|<=7.5 (numpy complex128), unten per Spiegelung F(conj z)=conj F(z).
Aufruf als Skript: baut die Koeffizienten und prueft gegen direkte gp-Werte (Gegenprobe).
"""
import json
from pathlib import Path
import numpy as np
import mpmath as mp
import motor

HIER = Path(__file__).resolve().parent
KOEF = HIER / "cache" / "taylor_koef.json"
ZENTREN = [0.0, 1.5, 3.0, 4.5, 6.0]
R = 2.0
M = 96


def bauen():
    if KOEF.is_file():
        return json.loads(KOEF.read_text())
    out = {}
    for yc in ZENTREN:
        c = mp.mpc(0.5, yc)
        Rc = 1.4 if yc == 0 else R
        pts = [c + Rc * mp.expj(2 * mp.pi * k / M) for k in range(M)]
        ausd = [f"sexp({motor.gpzahl(p, 25)})" for p in pts]
        v = motor.werte("exp(1)", ausd, 20)
        a = np.array([complex(x) for x in v])
        co = np.fft.fft(a) / M  # a_k R^k
        co = co / (Rc ** np.arange(M))
        out[str(yc)] = [[float(x.real), float(x.imag)] for x in co]
    KOEF.write_text(json.dumps(out))
    return out


_K = None


def feld(z):
    global _K
    if _K is None:
        raw = bauen()
        _K = {float(k): np.array([complex(a, b) for a, b in v]) for k, v in raw.items()}
    z = np.asarray(z, dtype=complex)
    unten = z.imag < 0
    zz = np.where(unten, np.conj(z), z)
    res = np.empty_like(zz)
    yc = np.array(ZENTREN)
    idx = np.argmin(np.abs(zz.imag[..., None] - yc), axis=-1)
    for i, y in enumerate(ZENTREN):
        m = idx == i
        if not m.any():
            continue
        co = _K[y][:60]
        d = zz[m] - complex(0.5, y)
        acc = np.zeros_like(d)
        for a in co[::-1]:
            acc = acc * d + a
        res[m] = acc
    return np.where(unten, np.conj(res), res)


if __name__ == "__main__":
    bauen()
    test = [0.0, 1.0, 0.3 + 0.4j, 0.9 + 1.1j, 0.1 + 2.2j, 0.7 + 3.8j, 0.2 + 5.3j, 0.95 + 7.4j, 0.5 - 2.0j]
    ausd = [f"sexp({motor.gpzahl(mp.mpc(t), 20)})" for t in test]
    ref = motor.werte("exp(1)", ausd, 20)
    for t, r in zip(test, ref):
        f = feld(np.array([t]))[0]
        print(t, f, "abw", abs(f - complex(r)))
    motor.schliessen()
