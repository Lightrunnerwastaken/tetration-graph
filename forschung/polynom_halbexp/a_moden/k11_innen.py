"""Strang A: Innenschicht an L. E_N(z) = p_N(z) − f_K(z) bei z = L(1 − v/N) (v = 0, 1, 2, 4, 8): skaliert E_N wie N^{−γ}?

f_K nahe L aus der exakten lokalen Formel f_K = φ(√L Ψ e^{Lε(Ψ^γ)}) (k5/k6; Ψ als Reihe, κ_m aus dem Motor).
Kontrolle von f_K an einem Punkt mit der Taylor-Reihe um 0 (1000 Glieder).
Ausgabe: Φ_N(v) = E_N(L(1 − v/N)) · N^{γ}; konvergiert Φ_N(v) für N → ∞, ist die Innenamplitude N·(1/N)^{1+γ}·const.
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 60
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
d = json.loads((DATEN / "singularitaet_theta.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
phi = [cx(p) for p in d["phi"]]
psi = [cx(p) for p in d["psi"]]
eps = [cx(p) for p in d["eps_q"]]
q = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
P = {int(k): [mp.mpf(x) for x in v] for k, v in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8")).items()}
pab = DATEN / "polynome_ab_195.json"
if pab.is_file():
    P.update({int(k): [mp.mpf(x) for x in v] for k, v in json.loads(pab.read_text(encoding="utf-8")).items()})


def pv(c, x):
    s = mp.mpc(0)
    for a in reversed(c):
        s = s * x + a
    return s


def fK(z):
    s = pv(psi, z - L)
    qq = mp.exp(G * mp.log(s))            # Hauptzweig: auf der Seite von 0 richtig (k5)
    e = mp.fsum(eps[k] * qq ** k for k in range(1, len(eps)))
    return pv(phi, SL * s * mp.exp(L * e))


z0 = L * mp.mpf("0.93")
print(f"Kontrolle f_K bei z = 0.93 L: lokal − Taylor(0) = {mp.nstr(fK(z0) - pv(q, z0), 3)}")
print(f"Kontrolle f_K(L) − L = {mp.nstr(fK(L * (1 - mp.mpf(10) ** -30)) - L, 3)}")
VS = [0, 1, 2, 4, 8]
Ns = sorted(n for n in P if n >= 40)
print("\n N    " + "   ".join(f"Φ_N({v})" + " " * 22 for v in VS))
werte = {}
for N in Ns:
    row = []
    for v in VS:
        z = L * (1 - mp.mpf(v) / N)
        E = pv(P[N], z) - (L if v == 0 else fK(z))
        row.append(E * mp.power(N, G))
    werte[N] = row
    if N % 10 == 0 or N == max(Ns):
        print(f"{N:4d}  " + "  ".join(f"{mp.nstr(x, 6):>28}" for x in row))
print("\nfreier Exponent von E_N(L) = p_N(L) − L (Anpassung log E = a + b log N + c/N, N ≥ 80):")
sel = [N for N in Ns if N >= 80]
A = mp.matrix([[1, mp.log(N), mp.mpf(1) / N] for N in sel])
ph, last = [], None
for N in sel:
    p = mp.arg(werte[N][0] / mp.power(N, G))
    if last is not None:
        while p - last > mp.pi:
            p -= 2 * mp.pi
        while p - last < -mp.pi:
            p += 2 * mp.pi
    ph.append(p)
    last = p
br = mp.qr_solve(A, mp.matrix([mp.log(abs(werte[N][0] / mp.power(N, G))) for N in sel]))[0]
bi = mp.qr_solve(A, mp.matrix(ph))[0]
print(f"   p_N(L) − L ~ N^({mp.nstr(br[1], 5)} + {mp.nstr(bi[1], 5)} i);  −γ = {mp.nstr(-G, 5)},  −1−γ = {mp.nstr(-1 - G, 5)}")
