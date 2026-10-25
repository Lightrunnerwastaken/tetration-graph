"""S1: Kennzahlen der Rang-4-Welle P4 = Phi4 - mu4, die über Glattheit / Singularität der Rang-5-Phase entscheiden.

Aus den Fourier-Moden von Phi4 (Runde 18, e1_*.json, nur gelesen) für mehrere Basispaare (b, d):
  pp      = max P4 - min P4                         (Spannweite)
  dmax/dmin = max / min P4'                         (Schranken für den Faktor je Ebene)
  var     = int_0^1 P4'^2                           (Varianz eines Schritts der Irrfahrt log Pi_n)
  chi     = -int_0^1 log(1 + P4'(t)) dt             (Drift von log Pi_n unter Lebesgue; > 0 nach Jensen)
  kl      =  int_0^1 (1+P4') log(1 + P4') dt        (Drift unter dem Mass dH; > 0)
  hell    = 1 - int_0^1 sqrt(1 + P4') dt            (Hellinger-Abstand je Ebene; Kakutani)
  c0      = min_u max_{s in [1,2]} |P4(u+s) - P4(u)|/s   (Untergrenze der Wackel-Amplitude, Satz A)
  n_stern = 1/chi                                   (Ebenen, bis Pi_n um e^-1 gefallen ist)
Gegenprobe: zwei Gitter (Quadratur mit 4000 und 8000 Punkten).
Aufruf: python s1_konstanten.py
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40
ALT = Path(__file__).resolve().parents[2] / "runde_2026_10_18" / "01_basiswechsel_raenge" / "daten"
PAARE = [("exp1", "2", 128), ("2", "exp1", 128), ("exp1", "3", 64), ("2", "3", 64), ("exp1", "10", 64), ("exp1", "4", 64),
         ("exp1", "9_5", 64), ("exp1", "14_5", 32), ("exp1", "27_10", 32)]


def nst(x, n):
    return mp.nstr(mp.mpf(float(x)), n)


def lade(b, d, M):
    f = ALT / f"e1_{b}_{d}_M{M}.json"
    dat = json.loads(f.read_text(encoding="utf-8"))
    a = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re, im in dat["f4"]]
    return a, dat


def P(a, t):        # P4(t) = 2 Re sum_{k>=1} a_k e^{2 pi i k t}
    return 2 * mp.re(mp.fsum(a[k] * mp.expjpi(2 * k * t) for k in range(1, len(a) - 1)))


def dP(a, t):
    return 2 * mp.re(mp.fsum(a[k] * 2j * mp.pi * k * mp.expjpi(2 * k * t) for k in range(1, len(a) - 1)))


def kennzahlen(a, G):
    # numpy (float64 reicht: P4' ~ 1e-3, gebraucht werden ~1e-10 relativ); Moden k = 1 .. len(a)-2
    ak = np.array([complex(x) for x in a[1:len(a) - 1]])
    k = np.arange(1, len(ak) + 1)
    t = np.arange(G) / G
    E = np.exp(2j * np.pi * np.outer(t, k))
    p = 2 * np.real(E @ ak)
    q = 2 * np.real(E @ (ak * 2j * np.pi * k))
    var = np.mean(q * q)
    chi = -np.mean(np.log1p(q))
    kl = np.mean((1 + q) * np.log1p(q))
    hell = 1 - np.mean(np.sqrt(1 + q))
    s = 1 + np.arange(0, 201) / 200
    c0 = np.inf
    for i in range(0, G, max(1, G // 1000)):
        w = np.abs(2 * np.real(np.exp(2j * np.pi * np.outer(t[i] + s, k)) @ ak) - p[i]) / s
        c0 = min(c0, w.max())
    return dict(pp=p.max() - p.min(), dmax=q.max(), dmin=q.min(), var=var, chi=chi, kl=kl, hell=hell, c0=c0)


def main():
    zeilen = []
    for b, d, M in PAARE:
        try:
            a, dat = lade(b, d, M)
        except FileNotFoundError:
            print("fehlt", b, d, M)
            continue
        k1 = kennzahlen(a, 4000)
        k2 = kennzahlen(a, 8000)
        rel = abs(k2["chi"] - k1["chi"]) / k1["chi"]
        print(f"({b},{d}): pp={nst(k1['pp'], 6)} P4' in [{nst(k1['dmin'], 6)}, {nst(k1['dmax'], 6)}] "
              f"var={nst(k1['var'], 6)} chi={nst(k1['chi'], 8)} (Gegenprobe rel {nst(rel, 2)}) var/2={nst(k1['var']/2, 8)} "
              f"kl={nst(k1['kl'], 6)} hell={nst(k1['hell'], 6)} var/8={nst(k1['var']/8, 6)} c0={nst(k1['c0'], 6)} "
              f"c0/pp={nst(k1['c0']/k1['pp'], 4)} n*=1/chi={nst(1/k1['chi'], 4)}", flush=True)
        zeilen.append({"b": b, "d": d, "M": M, **{k: nst(v, 12) for k, v in k1.items()}, "chi_8000": nst(k2["chi"], 12)})
    out = Path(__file__).resolve().parent / "daten" / "s1_konstanten.json"
    out.write_text(json.dumps(zeilen, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
