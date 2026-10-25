"""S2: Dichtheit der Verzweigungspunkte, vorgefuehrt.
Zu einem regulaeren Punkt z* des Hauptblatts werden Verzweigungspunkte (auf Nebenblaettern) gesucht, die immer naeher liegen.
Satz: z ist Projektion eines Verzweigungspunkts <=> sexp(z+M) = sexp(M-n) fuer ein M und ein n >= 2 (Hauptblatt).
Konstruktion: Bahn o_i = exp^[i](c*), c* = sexp(z*); Ziel T in exp^{-1}(E_j) = E_{j-1} + 2 pi i Z (j = 1..s-1) nahe o_{s-1};
dann rueckwaerts u_i = log(u_{i+1}) + 2 pi i m (Zweig naechst o_i). Dann exp^[s](u_0) = E_j exakt, also z ~ j-1-s.
Tiefe d = s - j: die Singularitaet ist vom Typ log^[d] (d = 1: reiner Logarithmus).
Aufruf: python s2_dicht.py DPS
"""
import sys, json
import mpmath as mp
from motor import werte, gpzahl, schliessen

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = 3 * DPS

def sexp(zs):
    return werte("exp(1)", [f"sexp({gpzahl(z, DPS + 10)})" for z in zs], DPS)

E = [mp.mpf(0)]
for _ in range(4):
    E.append(mp.exp(E[-1]))  # E_k = exp^[k](0) = sexp(k-1)

def naechster_log(u, ziel):
    v = mp.log(u)
    m = mp.nint((mp.im(ziel) - mp.im(v)) / (2 * mp.pi))
    return v + 2j * mp.pi * m

def suche(cstern, smax=7):
    o = [mp.mpc(cstern)]
    for _ in range(smax):
        if abs(o[-1]) > 1e8:
            break
        o.append(mp.exp(o[-1]))
    erg = []
    for s in range(2, len(o) + 1):
        best = None
        for j in range(1, s):
            if j - 1 >= len(E):
                break
            os_ = o[s - 1]
            for n in set([int(mp.nint(mp.im(os_) / (2 * mp.pi))) + dd for dd in (-1, 0, 1)] +
                         [sg * int(mp.nint(abs(os_) / (2 * mp.pi))) + dd for sg in (-1, 1) for dd in (-1, 0, 1)]):
                T = E[j - 1] + 2j * mp.pi * n
                u = T
                for i in range(s - 2, -1, -1):
                    u = naechster_log(u, o[i])
                t = u
                for _ in range(s - 1):
                    t = mp.exp(t)
                fehler = abs(t - T) / max(1, abs(T))  # exp^[s-1](u) = T, exp(T) = E_j exakt
                dist = abs(u - cstern)
                if best is None or dist < best[0]:
                    best = (dist, j, u, fehler)
        dist, j, u, fehler = best
        erg.append(dict(s=s, j=j, tiefe=s - j, u=u, dist=dist, fehler=fehler))
    return o, erg

aus = []
for zs in [mp.mpc("0.3", "0"), mp.mpc("1.3", "0.1"), mp.mpc("-0.5", "0.25")]:
    h = mp.mpf(10) ** (-(DPS // 3))
    a, b, c0 = sexp([zs + h, zs - h, zs])
    d0 = (a - b) / (2 * h)
    o, erg = suche(c0)
    print(f"z* = {mp.nstr(zs, 6)}: c* = sexp(z*) = {mp.nstr(c0, 12)}, |Bahn| = {[mp.nstr(abs(x), 4) for x in o]}", flush=True)
    for e in erg:
        # Newton im z-Raum: sexp(z) = u
        z = zs + (e["u"] - c0) / d0
        if mp.mpf(10) ** (-DPS + 12) < e["dist"] < 0.05:
            for _ in range(12):
                try:
                    a, b, f = sexp([z + h, z - h, z])
                except Exception:
                    z = mp.nan; break
                st = (f - e["u"]) / ((a - b) / (2 * h))
                z -= st
                if abs(st) < mp.mpf(10) ** (-DPS + 8):
                    break
        print(f"  s={e['s']} Tiefe={e['tiefe']} (Klasse von {e['j']-1-e['s']}): |c-c*|={mp.nstr(e['dist'], 4)}  "
              f"z={mp.nstr(z, 18)}  |z-z*|={mp.nstr(abs(z - zs), 4)}  Probe={mp.nstr(e['fehler'], 3)}", flush=True)
        aus.append([mp.nstr(zs, 10), e["s"], e["tiefe"], mp.nstr(e["dist"], 6), mp.nstr(mp.re(z), DPS - 8), mp.nstr(mp.im(z), DPS - 8),
                    mp.nstr(abs(z - zs), 6)])
json.dump(aus, open(f"daten/dicht_{DPS}.json", "w"), indent=1)
schliessen()
