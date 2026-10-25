"""F1: Welle P(u) = u + 1 - sexp(u) auf [-1,0] (1-periodisch fortgesetzt), Basis e, Kneser.
- Spruenge J_k = sexp^(k)(0) - sexp^(k)(-1) der k-ten Ableitung von P an den ganzen Zahlen (J_0 = J_1 = 0 fuer P).
- Fourier-Koeffizienten a_n = int_{-1}^0 P(u) e^{-2 pi i n u} du, direkt per Quadratur, gegen Reihe sum_k J_k/(2 pi i n)^{k+1}.
- Bernoulli-Entwicklung P(u) = a_0 - sum_{k>=3} B_k(u+1) J_{k-1}/k! + Rest.
- Ableitungsschranken: Q(u) = 1/sexp'(u) = slog'/slog_lin'.
Zwei Taylor-Saetze um 0 (taylor_0_1_512_60, taylor_0_1.5_400_50) als Gegenprobe; Taylor um -1 aus taylor_-1_0.6_400_60.
"""
import json, sys
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
mp.mp.dps = 32
satz = sys.argv[1] if len(sys.argv) > 1 else "taylor_0_1_512_60.json"
C0 = [mp.mpf(s) for s in json.load(open(HIER / "daten" / satz))][:200]
Cm = [mp.mpf(s) for s in json.load(open(HIER / "daten" / "taylor_-1_0.6_400_60.json"))][:120]
S = lambda u: mp.polyval(C0[::-1], u)
d1 = [k * C0[k] for k in range(1, len(C0))]
d2 = [k * (k - 1) * C0[k] for k in range(2, len(C0))]
d3 = [k * (k - 1) * (k - 2) * C0[k] for k in range(3, len(C0))]
S1 = lambda u: mp.polyval(d1[::-1], u)
S2 = lambda u: mp.polyval(d2[::-1], u)
S3 = lambda u: mp.polyval(d3[::-1], u)
P = lambda u: u + 1 - S(u)

print("Satz:", satz)
c = C0[1]
print("c = sexp'(0) = sexp'(-1) =", mp.nstr(c, 30), " Taylor(-1) gibt", mp.nstr(Cm[1], 30))
# Spruenge J_k
K = 40
J = [mp.factorial(k) * (C0[k] - Cm[k]) for k in range(K)]
print("J_0..J_3:", [mp.nstr(x, 20) for x in J[:4]], " c^2 =", mp.nstr(c**2, 20))
print("J_3 Formel c^3+3c*sexp''(-1) =", mp.nstr(c**3 + 3 * c * 2 * Cm[2], 20))
for k in range(2, 26):
    print(f"  J_{k:2d} = {mp.nstr(J[k], 18):>26}   J_k/k! = {mp.nstr(J[k]/mp.factorial(k), 12)}")

# Fourier-Koeffizienten direkt
# Gauss-Legendre: 8 Teilintervalle zu 96 Knoten, P einmal pro Knoten
_gl = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)
KN, GW = [], []
TEIL = 8
for j in range(TEIL):
    lo, hi = mp.mpf(-1) + mp.mpf(j) / TEIL, mp.mpf(-1) + mp.mpf(j + 1) / TEIL
    for x, w in _gl:
        KN.append((lo + hi) / 2 + (hi - lo) / 2 * x); GW.append(w * (hi - lo) / 2)
PK = [P(u) for u in KN]
def a_direkt(n):
    return mp.fsum(w * p * mp.expj(-2 * mp.pi * n * u) for u, w, p in zip(KN, GW, PK))

a0 = a_direkt(0)
print("a_0 =", mp.nstr(a0, 30), " (1/2 - int sexp)")
def a_reihe(n, kmax):
    s = mp.mpc(0)
    for k in range(2, kmax + 1):
        s += J[k] / (2j * mp.pi * n) ** (k + 1)
    return s

tab = []
for n in [1, 2, 3, 4, 5, 6, 8, 10, 15, 20]:
    an = a_direkt(n)
    # beste Kuerzung der Reihe
    best = None
    for kmax in range(2, K):
        e = abs(an - a_reihe(n, kmax))
        if best is None or e < best[1]:
            best = (kmax, e)
    lead = J[2] / (2j * mp.pi * n) ** 3
    tab.append((n, an, best))
    print(f"n={n:2d} a_n={mp.nstr(an.real,12):>18} {mp.nstr(an.imag,15):>20}i  |a_n/fuehrend|={mp.nstr(abs(an/lead),10)}  beste Kuerzung k<={best[0]} Fehler {mp.nstr(best[1],3)}  e^(-2pi n)={mp.nstr(mp.exp(-2*mp.pi*n),3)}")

# Rekonstruktion von P aus Fourier-Reihe (Probe)
N = 60
an_list = {n: a_direkt(n) for n in range(1, N + 1)}
def P_fourier(u):
    return a0 + sum(2 * mp.re(an_list[n] * mp.expj(2 * mp.pi * n * u)) for n in an_list)
for u in [mp.mpf(-0.8), mp.mpf(-0.5), mp.mpf(-0.2)]:
    print(f"  P({u}) = {mp.nstr(P(u), 15)}   Fourier N=60: {mp.nstr(P_fourier(u), 15)}")

# Bernoulli-Entwicklung
gitter = mp.linspace(-1, 0, 101)
print("Bernoulli: max |P - (a0 - sum_{k=3}^{K} B_k(u+1) J_{k-1}/k!)| auf [-1,0]")
for Kb in [3, 4, 5, 6, 7, 8, 10, 12, 15, 20]:
    def PB(u):
        return a0 - sum(mp.bernpoly(k, u + 1) * J[k - 1] / mp.factorial(k) for k in range(3, Kb + 1))
    err = max(abs(P(u) - PB(u)) for u in gitter)
    print(f"  K={Kb:2d}: {mp.nstr(err, 5)}")
print("max|P| =", mp.nstr(max(abs(P(u)) for u in gitter), 8), " c^2/6 * max|B3| =", mp.nstr(c**2 / 6 * mp.sqrt(3) / 36, 8))

# Ableitungsschranken / Konkavitaet
xw = mp.findroot(S2, -0.52)
m = S1(xw)
print("Wendepunkt sexp in [-1,0]: x_w =", mp.nstr(xw, 30), " sexp(x_w) =", mp.nstr(S(xw), 30))
print("m = sexp'(x_w) = min sexp' auf [-1,0] =", mp.nstr(m, 30))
print("Schranken Q = slog'/slog_lin' in [1/c, 1/m] = [", mp.nstr(1 / c, 25), ",", mp.nstr(1 / m, 25), "]")
print("Probe: Vorzeichenwechsel von sexp'' auf [-1,0]:", sum(1 for i in range(len(gitter) - 1) if S2(gitter[i]) * S2(gitter[i + 1]) < 0))
ex = [mp.findroot(lambda u: S1(u) - 1, u0) for u0 in (-0.8, -0.2)]
print("sexp'=1 bei u =", [mp.nstr(u, 25) for u in ex], " P dort", [mp.nstr(P(u), 25) for u in ex])
print("  zugehoerige x = sexp(u) (slog'(x)=1 in (0,1)):", [mp.nstr(S(u), 25) for u in ex])
print("Wendestelle der Welle P (P''=0 <-> sexp''=0): u = x_w, P(x_w) =", mp.nstr(P(xw), 25))
print("P'(x_w) = 1 - m =", mp.nstr(1 - m, 25), "(max Steigung der Welle);  min Steigung 1-c =", mp.nstr(1 - c, 25))
print("sexp'''(x_w) =", mp.nstr(S3(xw), 20))
