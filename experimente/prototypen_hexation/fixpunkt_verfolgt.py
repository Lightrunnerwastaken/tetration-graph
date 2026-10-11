"""Fixpunkte der von der reellen Achse aus fortgesetzten Pentation zur Basis 2 (Blattverfolgung entlang fester Wegarten)."""
import sys, time
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from blatt import BlattPentation, verfolge
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

BP = BlattPentation("2")
q = mp.mpc("-1.992948435932322332933652", "1.445128496639205463533289")


def weg_gerade(z):
    return lambda t: 1 + t * (z - 1)


def weg_oben(z):
    """1 → 1 + i·Im z (senkrecht) → z (waagrecht); t ∈ [0, 1], Knick bei t = 1/2."""
    ecke = 1 + 1j * mp.im(z)
    return lambda t: 1 + 2 * t * (ecke - 1) if t <= mp.mpf(1) / 2 else ecke + (2 * t - 1) * (z - ecke)


def P_bis(weg, z, n=80):
    w = weg(z)
    ts = [mp.mpf(k) / n for k in range(n + 1)]
    param, werte, smax = verfolge(ts, lambda tt: BP.P_pfad([w(t) for t in tt]))
    return werte[-1], smax


def fixpunkt(weg, z0, name):
    z = mp.mpc(z0)
    h = mp.mpf("1e-10")
    for it in range(25):
        f0, s = P_bis(weg, z)
        f1, _ = P_bis(weg, z + h)
        d = (f1 - f0) / h
        dz = (f0 - z) / (d - 1)
        log(f"  {name} Newton {it}: z = {mp.nstr(z, 15)}, |P(z) − z| = {mp.nstr(abs(f0 - z), 3)}, P′ = {mp.nstr(d, 8)}, Winkelsprung {mp.nstr(s, 3)}")
        z -= dz
        if abs(dz) < mp.mpf(10) ** -20:
            break
    return z, d


for name, weg in (("gerade", weg_gerade), ("oben herum", weg_oben)):
    v, s = P_bis(weg, q)
    log(f"Weg {name}: P(q) fortgesetzt = {mp.nstr(v, 15)} (Hauptzweig: q), Abstand {mp.nstr(abs(v - q), 3)}")
    z, mu = fixpunkt(weg, q, name)
    log(f"Weg {name}: Fixpunkt z* = {mp.nstr(z, 22)}, Multiplikator {mp.nstr(mu, 12)}, |μ| = {mp.nstr(abs(mu), 8)}")
