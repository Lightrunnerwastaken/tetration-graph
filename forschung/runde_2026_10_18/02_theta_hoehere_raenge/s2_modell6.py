"""Rang 6: Singularitaetsmodell fuer theta6 und Vergleich mit den gemessenen a_n (s1_theta6.py).
Vorhersage: theta6 ist genau bei x* + Z singulaer, x* = ihex(p), p = reeller Fixpunkt von sexp (= Wert, den die Pentation
fuer Re z -> -inf annimmt und ausl.). Nahe x*:  theta6(z) = 1 - z + h(l(z)),  l = [log sigma_p(hex z - p)]/ln lam_p - pcent.
Koeffizienten: a_n = e^{-2 pi i n x*} b_n + ...,  b_n = [Q^n] h((c6 + log(1-Q))/ln lam_p),
c6 = ln(A6/2pi) + i pi/2 - pcent ln lam_p,  A6 = hex'(x*);  zweite Ordnung mit delta1 = 1/2 + i kappa6/2pi.
Aufruf: python s2_modell6.py daten/theta6_40_0.04_512.json daten/theta6_30_0.05_512.json"""
import json, sys, time
from pathlib import Path
import mpmath as mp
import motor
from lib6 import Rang6
DPS = 40
t0 = time.time()
R = Rang6("2", DPS, start=mp.mpc(-1.99294843593232, 1.44512849663921), rho=mp.mpf("0.2"), rad=mp.mpf("0.04"))
mp.mp.dps = DPS + 10
motor.roh("2", ["pentx(0.5)", "hexx(0.5)"], DPS)
k = motor.werte("2", ["tfixed", "log(lambdat)", "pcent", "ihexx(tfixed)", "polcoeff(ipentz,1)", "polcoeff(ipentz,2)", "pentdigits", "hexdigits"], DPS)
p, llp, pcent, xs, s1, s2 = [mp.re(v) for v in k[:6]]
print("p =", mp.nstr(p, 30), " ln lam_p =", mp.nstr(llp, 25), " pcent =", mp.nstr(pcent, 25), " x* =", mp.nstr(xs, 30))
print("sigma_p'(0) =", mp.nstr(s1, 10), " s2 =", mp.nstr(s2, 20), " Stellen pent/hex:", mp.nstr(k[6].real, 5), mp.nstr(k[7].real, 5))
# ihex'(p), ihex''(p) per Cauchy
M, r = 32, mp.mpf("0.02")
pts = [p + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
v = motor.werte("2", [f"ihexx({motor.gpzahl(z, DPS + 10)})" for z in pts], DPS)
c = [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * kk) / M) for j in range(M)) / M / r ** kk for kk in range(4)]
d1, d2 = mp.re(c[1]), mp.re(2 * c[2])
A6 = 1 / d1; hpp = -d2 / d1 ** 3
kappa6 = hpp / (2 * A6) + s2 / s1 * A6
print("A6 = hex'(x*) =", mp.nstr(A6, 25), " hex''(x*) =", mp.nstr(hpp, 20), " kappa6 =", mp.nstr(kappa6, 20))
c6 = mp.log(A6 / (2 * mp.pi)) + 1j * mp.pi / 2 - pcent * llp
d1c = mp.mpf(1) / 2 + 1j * kappa6 / (2 * mp.pi)
# b_n per Cauchy-Integral auf |Q| = rq
MQ = 512; yq = mp.mpf("0.05"); rq = mp.exp(-2 * mp.pi * yq)
Qs = [rq * mp.expjpi(2 * mp.mpf(j) / MQ) for j in range(MQ)]
vs = [(c6 + mp.log(1 - Q)) / llp for Q in Qs]
hv, ns = R.h6(vs)
per = 2j * mp.pi / R.logmu
for j in range(1, MQ):
    m = min(range(-5, 6), key=lambda m: abs(hv[j] + m * per - hv[j - 1])); hv[j] += m * per
print(f"[{time.time()-t0:.0f}s] h an {MQ} Punkten (ipent-Schritte {min(ns)}..{max(ns)}), Randschritt {mp.nstr(min(abs(hv[0]+m*per-hv[-1]) for m in range(-5,6)),3)}")
def fft(a):
    n = len(a)
    if n == 1: return a[:]
    ev, od = fft(a[0::2]), fft(a[1::2]); out = [None] * n
    for kk in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(kk) / n) * od[kk]; out[kk] = ev[kk] + tw; out[kk + n // 2] = ev[kk] - tw
    return out
Fb = fft(hv)
b = [Fb[n] / MQ / rq ** n for n in range(MQ // 2)]
lade = lambda f: [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in json.loads(Path(f).read_text())["a"]]
a1, a2 = lade(sys.argv[1]), lade(sys.argv[2])
rot = lambda n: mp.expjpi(2 * n * xs)        # e^{+2 pi i n x*}
zeilen = []
print(" n   |a_n|          |a-a'|/|a|   |a e^{2pi i n x*} - b|/|a|  n*rel1   |..-b2|/|a|  n^2*rel2   arg(a e^{2pi i n x*})")
for n in [1, 2, 4, 8, 16, 32, 48, 64, 96, 128, 160]:
    an = a1[n]; prec = abs(an - a2[n]) / abs(an)
    bb = an * rot(n)
    r1 = abs(bb - b[n]) / abs(an)
    b2 = b[n] - d1c * ((n + 1) * b[n + 1] - 2 * n * b[n] + (n - 1) * b[n - 1])
    r2 = abs(bb - b2) / abs(an)
    print(f"{n:4d} {mp.nstr(abs(an),8):>13} {mp.nstr(prec,2):>9}  {mp.nstr(r1,3):>9}  {mp.nstr(n*r1,3):>8}  {mp.nstr(r2,3):>9} {mp.nstr(n*n*r2,3):>8}  {mp.nstr(mp.arg(bb),8)}")
    zeilen.append([n, float(abs(an)), float(prec), float(r1), float(r2), float(mp.arg(bb))])
# fuehrender Term
def hprime(vv):
    return mp.diff(lambda t: R.h6([t])[0][0], vv)
fuehr = []
for n in [16, 64, 128]:
    vn = (c6 - mp.log(n)) / llp
    Tn = -hprime(vn) / (n * llp)
    print(f"fuehrend n={n}: T_n/(a_n e^(2pi i n x*)) = {mp.nstr(Tn/(a1[n]*rot(n)),8)}")
    fuehr.append([n, str(Tn / (a1[n] * rot(n)))])
out = {"p": str(p), "ln_lam_p": str(llp), "pcent": str(pcent), "x_stern": str(xs), "A6": str(A6), "kappa6": str(kappa6),
       "c6": str(c6), "zeilen": zeilen, "fuehrend": fuehr}
(Path(__file__).resolve().parent / "daten" / "modell6.json").write_text(json.dumps(out, indent=1))
print(f"[{time.time()-t0:.0f}s] fertig")
