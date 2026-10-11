"""Exponent der Potenz-Verzweigung von F_q bei ζ_c + Π (Wert P(L_T)) und ζ_c − 1 + Π (Wert L_T): |F_q(ζ₀ + i r) − w₀| ~ C r^α.
 Vorhersage α = ln λ_T / ln 2 (λ_T = Multiplikator der Tetration am Fixpunkt L_T). Weg senkrecht von oben, mit Blattverfolgung."""
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
alpha = mp.log(BP.Pe.lam) / mp.log(2)
log("Vorhersage α = ln λ_T / ln 2 =", mp.nstr(alpha, 10))


def F_weg(zetas, n):
    w = [psi_wert(HQ.q, HQ.psi, mp.exp((z - n) * HQ.lnmu)) for z in zetas]
    sp = [mp.mpf(0)] * len(zetas)
    for _ in range(n):
        w, s_ = BP.P_pfad(w)
        sp = [max(a, b) for a, b in zip(sp, s_)]
    return w, sp


for name, z0, w0, n in (("ζ_c + Π (Wert P(L_T))", zeta_c + PI_, BP.Pe.P([LT])[0], 7), ("ζ_c − 1 + Π (Wert L_T)", zeta_c - 1 + PI_, LT, 6)):
    rs = [mp.mpf(2) ** -k / 5 for k in range(0, 7)]          # 0.2 … 0.003125
    oben = z0 + 2.5j
    ts = [mp.mpf(k) / 400 for k in range(401)]
    ziel = [z0 + 1j * r for r in rs]
    pts = sorted(set([oben + (z0 + 1j * rs[-1] - oben) * t for t in ts] + ziel), key=lambda z: -mp.im(z))
    p, w, sm = verfolge(pts, lambda zz: F_weg(zz, n), runden=10)
    werte = {mp.nstr(mp.im(z - z0), 5): abs(v - w0) for z, v in zip(p, w) if any(abs(z - zz) < mp.mpf(10) ** -30 for zz in ziel)}
    log(f"{name}: Sprung {mp.nstr(sm, 3)}; |F − w₀| je r:", {k: mp.nstr(v, 6) for k, v in werte.items()})
    ks = list(werte.keys())
    steig = [mp.log(werte[ks[i]] / werte[ks[i + 1]]) / mp.log(mp.mpf(ks[i]) / mp.mpf(ks[i + 1])) for i in range(len(ks) - 1)]
    log("   lokale Exponenten:", [mp.nstr(s, 6) for s in steig])


# Umlauf-Faktor: (F_nach − w₀)/(F_vor − w₀) = e^{2πiα} an derselben Stelle (unabhängig vom genauen Mittelpunkt)
log("Vorhersage e^{2πiα} =", mp.nstr(mp.expjpi(2 * alpha), 10))
for name, z0, w0, n in (("ζ_c + Π", zeta_c + PI_, BP.Pe.P([LT])[0], 7), ("ζ_c − 1 + Π", zeta_c - 1 + PI_, LT, 6)):
    for r in (mp.mpf("0.05"), mp.mpf("0.02")):
        oben = z0 + 2.5j
        runter = [oben + (z0 + 1j * r - oben) * mp.mpf(k) / 200 for k in range(201)]
        kreis = [z0 + r * mp.expjpi(mp.mpf(1) / 2 + 2 * mp.mpf(k) / 200) for k in range(1, 201)]
        p, w, sm = verfolge(runter + kreis, lambda zz: F_weg(zz, n), runden=10)
        i0 = next(i for i, z in enumerate(p) if abs(z - (z0 + 1j * r)) < mp.mpf(10) ** -30)
        q_ = (w[-1] - w0) / (w[i0] - w0)
        a_ = mp.re(mp.log(q_) / (2j * mp.pi))
        log(f"{name}, r = {r}: Faktor {mp.nstr(q_, 8)}, also α ≡ {mp.nstr(a_ % 1, 8)} (mod 1); Sprung {mp.nstr(sm, 3)}")
