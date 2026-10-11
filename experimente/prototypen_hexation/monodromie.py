"""Monodromie von F_q (reguläre Superfunktion der Pentation P_2 am Fixpunkt q) um die vermuteten Verzweigungspunkte.

 F_q ist exakt periodisch mit Π = 2πi/ln μ_q. Vermutung: F_q ist verzweigt bei ζ_c − 2 + m + kΠ (m ≥ 0), ζ_c = Ecke der reellen Kurve
 (F_q(ζ_c) = P(L_T)); bei m = 0 logarithmisch (F_q → −∞, Sprung 2πi/ln 2 um den Punkt), sonst Potenz-artig (Exponent ln λ_T/ln 2).
 Test: Weg von oben (Hauptzweig dort richtig) herunter, ein Kreis um den Punkt, wieder hinauf; F_q mit Blattverfolgung entlang des Weges.
 Kehrt der Wert oben nicht zurück, ist der Punkt ein Verzweigungspunkt. Kontrolle: Kreis um einen Punkt dazwischen.
"""
import sys, time
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from blatt import BlattPentation, verfolge
from hexq import HexQ
from stufen import psi_wert
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
BP = BlattPentation("2")
PI_ = 2j * mp.pi / HQ.lnmu
zeta_c = mp.mpc("0.40312930", "-1.70978835")
LT = BP.Pe.L
log("Periode Π =", mp.nstr(PI_, 10), " L_T =", mp.nstr(LT, 12), " P(L_T) =", mp.nstr(BP.Pe.P([LT])[0], 12))


def F_weg(zetas, n):
    s = [mp.exp((z - n) * HQ.lnmu) for z in zetas]
    if max(abs(x) for x in s) > HQ.rlok:
        raise ValueError("n zu klein")
    w = [psi_wert(HQ.q, HQ.psi, x) for x in s]
    sp = [mp.mpf(0)] * len(zetas)
    for _ in range(n):
        w, s_ = BP.P_pfad(w)
        sp = [max(a, b) for a, b in zip(sp, s_)]
    return w, sp


def schleife(c, r=mp.mpf("0.2"), hoehe=mp.mpf("2.5"), N=160):
    oben = c + 1j * hoehe
    unten = c + 1j * r
    runter = [oben + (unten - oben) * mp.mpf(k) / 60 for k in range(61)]
    kreis = [c + r * mp.expjpi(mp.mpf(1) / 2 + 2 * mp.mpf(k) / N) for k in range(1, N + 1)]
    hoch = runter[::-1][1:]
    return runter + kreis + hoch


def test(name, c):
    pts = schleife(c)
    n = max(0, max(int(mp.ceil(mp.re((mp.log(abs(mp.exp(z * HQ.lnmu))) - mp.log(HQ.rlok)) / mp.re(HQ.lnmu)))) for z in pts)) + 1
    def weg(tt):
        out = []
        for t in tt:
            j = min(int(mp.floor(t)), len(pts) - 2)
            out.append(pts[j] + (t - j) * (pts[j + 1] - pts[j]))
        return F_weg(out, n)
    p, w, sm = verfolge([mp.mpf(j) for j in range(len(pts))], weg, runden=14)
    haupt = HQ.F([pts[0]])[0]
    log(f"{name}: Mitte {mp.nstr(c, 8)}, Kette n = {n}, {len(p)} Punkte, Sprung {mp.nstr(sm, 3)}")
    log(f"     Start {mp.nstr(w[0], 12)} (Hauptzweig {mp.nstr(haupt, 12)}), nach dem Umlauf {mp.nstr(w[-1], 12)}, Differenz {mp.nstr(w[-1] - w[0], 6)}")
    return w[-1] - w[0]


for name, c in (("Log-Punkt  m = 0 (ζ_c − 2 + Π)", zeta_c - 2 + PI_),
                ("Kontrolle zwischen m = 0 und 1", zeta_c - mp.mpf("1.5") + PI_),
                ("Potenz-Punkt m = 1 (ζ_c − 1 + Π)", zeta_c - 1 + PI_),
                ("Potenz-Punkt m = 2 (ζ_c + Π)", zeta_c + PI_),
                ("Kontrolle zwischen m = 1 und 2", zeta_c - mp.mpf("0.5") + PI_)):
    test(name, c)
log(f"zum Vergleich 2πi/ln 2 = {mp.nstr(2j * mp.pi / mp.log(2), 10)}")
