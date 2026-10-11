"""Strang A: Amplitude c_N des Störmodus für alle N = 30..194 und Anpassung c_N = C·N^{1−2γ}(1 + b₁/N + b₂/N²) + Schwingung.

c_N := c_{−1}(N) (Koeffizient von e^{−2πiw} in θ_N, Mode singulär bei L), gemessen wie in k3b (Linie y = 0, T(w) vom Motor).
Normierte Amplitude ρ_N = c_N · L^{3/2} e^{−2πiκ₀} / (d₁ N^{1−2γ}): Verhältnis der t^{1−γ}-Amplitude von c_N δ_{−1} an L
(δ_{−1} ≈ L^{3/2} e^{−2πiκ₀} t^{1−γ}) zur Kneser-Leitamplitude d₁ (t^{1+γ}).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 60
ctx.prec = 260
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
sing = json.loads((GEM / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
d1 = cx(sing["d_n"]["1"])
k0 = cx(sing["kappa"]["0"])
q = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
P = {int(k): [mp.mpf(x) for x in v] for k, v in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8")).items()}
kd = json.loads((DATEN / "kreis_r1.34_M4096.json").read_text(encoding="utf-8"))
M = kd["M"]
Ah = [mp.mpc(*x) for x in kd["A"]]
Av = Ah + [mp.conj(x) for x in Ah[1:-1][::-1]]
dd = acb.dft([acb(arb(mp.nstr(mp.re(x), 58)), arb(mp.nstr(mp.im(x), 58))) for x in Av])
Ak = [mp.mpf(dd[k].real.mid().str(58, radius=False)) / M / mp.mpf("1.34") ** k for k in range(1200)]


def reihe(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def abl(c, z):
    return reihe([k * c[k] for k in range(1, len(c))], z)


m = 48
X0 = mp.mpf("-1.6")
ws = [X0 + mp.mpf(j) / m for j in range(m)]
gp = FatouGP(dps=57, fatou_gp="fork")
zs = [mp.re(z) for z in gp.session("exp(1)").sexp_batch(ws)]
gp.close()
Tp = [abl(q, z) / abl(Ak, z) for z in zs]
fz = [reihe(q, z) for z in zs]
norm = L * SL * mp.exp(-2j * mp.pi * k0) / d1
cN, rho = {}, {}
for N in range(30, 195):
    th = [(reihe(P[N], z) - f) / t for z, f, t in zip(zs, fz, Tp)]
    c = mp.fsum(t * mp.expjpi(2 * w) for t, w in zip(th, ws)) / m          # Koeffizient von e^{−2πiw}
    cN[N] = c
    rho[N] = c * norm / mp.power(N, 1 - 2 * G)
print(" N    |ρ_N|        arg ρ_N     (ρ_N = c_N L^{3/2} e^{−2πiκ₀} / (d₁ N^{1−2γ}))")
for N in range(30, 195, 8):
    print(f"{N:4d}  {mp.nstr(abs(rho[N]), 6):>10}  {mp.nstr(mp.arg(rho[N]), 6):>9}")
# Anpassung log ρ_N = log C' + b1/N + b2/N^2 + Schwingung (cos/sin(N arg L)/N) über N = 60..194 (kleinste Quadrate, komplex)
Ns = list(range(60, 195))
th = mp.im(L) and mp.arg(L)
basis = lambda N: [mp.mpf(1), mp.mpf(1) / N, mp.mpf(1) / N ** 2, mp.cos(N * th) / N, mp.sin(N * th) / N,  # noqa: E731
                   mp.cos(N * th) / N ** 2, mp.sin(N * th) / N ** 2]
A = mp.matrix([basis(N) for N in Ns])
yre = mp.matrix([mp.log(abs(rho[N])) for N in Ns])
# Phase stetig machen
ph = []
for N in Ns:
    p = mp.arg(rho[N])
    if ph:
        while p - ph[-1] > mp.pi:
            p -= 2 * mp.pi
        while p - ph[-1] < -mp.pi:
            p += 2 * mp.pi
    ph.append(p)
yim = mp.matrix(ph)
xr = mp.qr_solve(A, yre)[0]
xi = mp.qr_solve(A, yim)[0]
C = mp.exp(xr[0] + 1j * xi[0])
print(f"\nAnpassung N = 60..194: ρ_∞ = {mp.nstr(C, 8)} (|ρ_∞| = {mp.nstr(abs(C), 6)}), b₁ = {mp.nstr(xr[1] + 1j * xi[1], 5)}, "
      f"b₂ = {mp.nstr(xr[2] + 1j * xi[2], 5)}")
res = max(abs(mp.fsum(a * b for a, b in zip(basis(N), xr)) - mp.log(abs(rho[N]))) for N in Ns)
print(f"   max Rest in log|ρ|: {mp.nstr(res, 3)}")
# freier Exponent: log c_N = a + β log N + b/N + Schwingung
B2 = lambda N: [mp.mpf(1), mp.log(N), mp.mpf(1) / N, mp.cos(N * th) / N, mp.sin(N * th) / N]  # noqa: E731
A2 = mp.matrix([B2(N) for N in Ns])
phc = []
for N in Ns:
    p = mp.arg(cN[N])
    if phc:
        while p - phc[-1] > mp.pi:
            p -= 2 * mp.pi
        while p - phc[-1] < -mp.pi:
            p += 2 * mp.pi
    phc.append(p)
br = mp.qr_solve(A2, mp.matrix([mp.log(abs(cN[N])) for N in Ns]))[0]
bi = mp.qr_solve(A2, mp.matrix(phc))[0]
print(f"freier Exponent (mit b/N): c_N ~ N^({mp.nstr(br[1], 5)} + {mp.nstr(bi[1], 5)} i);  Soll 1 − 2γ = {mp.nstr(1 - 2 * G, 6)}")
A3 = mp.matrix([[mp.mpf(1), mp.log(N), mp.cos(N * th) / N, mp.sin(N * th) / N] for N in Ns])
br3 = mp.qr_solve(A3, mp.matrix([mp.log(abs(cN[N])) for N in Ns]))[0]
bi3 = mp.qr_solve(A3, mp.matrix(phc))[0]
print(f"freier Exponent (ohne b/N): c_N ~ N^({mp.nstr(br3[1], 5)} + {mp.nstr(bi3[1], 5)} i)")
json.dump({str(N): {"c": [mp.nstr(mp.re(cN[N]), 20), mp.nstr(mp.im(cN[N]), 20)], "rho": [mp.nstr(mp.re(rho[N]), 12), mp.nstr(mp.im(rho[N]), 12)]}
           for N in cN}, open(DATEN / "amplitude_cN.json", "w", encoding="utf-8"), indent=0)
