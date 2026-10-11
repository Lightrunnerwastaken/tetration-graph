"""Strang A: universelle Konstante der Amplitudenregel und Vorhersage p_N(1) − T(½) für N > 194.

Regel (Strang B, generisch): c_N ≈ K · N · u_N / (δ_{−1})_N, u_N = [z^N] D₁(1 − z/L)^{1+γ}, (δ_{−1})_N = [z^N] D_δ(1 − z/L)^{1−γ},
D_δ = L^{3/2} e^{−2πiκ₀} (−L)^{1−γ} (lokale Form von δ_{−1} an L). Mit [z^N](1−z/L)^α ≈ L^{−N} N^{−α−1}/Γ(−α):
K = c_N Γ(−1−γ) D_δ / (N^{1−2γ} D₁ Γ(−1+γ)) = ρ_N (−L)^{−2γ} Γ(−1−γ)/Γ(−1+γ)  (ρ_N aus k7, exakte Γ-Quotienten zusätzlich).
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
DATEN = HIER / "daten"
mp.mp.dps = 40
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
LML = mp.log(-L)
sing = json.loads((GEM / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
d1, k0 = cx(sing["d_n"]["1"]), cx(sing["kappa"]["0"])
D1 = d1 * mp.exp((1 + G) * LML)
Dd = L * SL * mp.exp(-2j * mp.pi * k0) * mp.exp((1 - G) * LML)
am = json.loads((DATEN / "amplitude_cN.json").read_text(encoding="utf-8"))


def koef(al, N):
    return mp.exp(mp.loggamma(N - al) - mp.loggamma(-al) - mp.loggamma(N + 1)) * L ** (-N)


print(" N    K_N = c_N (δ_{−1})_N / (N u_N)   (exakte Γ)        |K_N|")
Ks = {}
for N in range(40, 195):
    c = cx(am[str(N)]["c"])
    uN = D1 * koef(1 + G, N)
    dN = Dd * koef(1 - G, N)
    Ks[N] = c * dN / (N * uN)
    if N % 15 == 0 or N == 194:
        print(f"{N:4d}  {mp.nstr(Ks[N], 10):>34}  {mp.nstr(abs(Ks[N]), 6)}")
# Mittel über eine Schwingungsperiode (≈ 4.7) am Ende, und 1/N-Extrapolation aus Mitteln
def mittel(a, b):
    return mp.fsum(Ks[N] for N in range(a, b + 1)) / (b - a + 1)


m1, m2 = mittel(100, 128), mittel(166, 194)
x1, x2 = mp.mpf(1) / 114, mp.mpf(1) / 180
Kinf = m2 + (m2 - m1) * (0 - x2) / (x2 - x1)
print(f"Mittel N = 100..128: {mp.nstr(m1, 8)}; N = 166..194: {mp.nstr(m2, 8)}; linear in 1/N extrapoliert: K_∞ ≈ {mp.nstr(Kinf, 6)}")
# Vorhersage p_N(1) − T(½) für N = 195..240 aus c_N = K_N-Extrapolation: c_N = K·N·u_N/(δ_{−1})_N, K ≈ m2 (+ Schwingung ignoriert)
T1 = mp.mpf("1.5751579377843008124")
print("\nVorhersage (K = Mittel 166..194, Schwingung ±3 % ignoriert): p_N(1) − T(½) ≈ 2 T'(½) Re c_N")
vor = {}
for N in list(range(196, 241, 2)):
    c = m2 * N * D1 * koef(1 + G, N) / (Dd * koef(1 - G, N))
    vor[N] = 2 * T1 * mp.re(c)
    print(f"  N = {N}: {mp.nstr(vor[N], 4)}   (arg c = {mp.nstr(mp.arg(c), 4)})")
json.dump({"K_N": {str(N): [mp.nstr(mp.re(v), 12), mp.nstr(mp.im(v), 12)] for N, v in Ks.items()},
           "vorhersage_pN1": {str(N): mp.nstr(v, 6) for N, v in vor.items()}}, open(DATEN / "konstante_K.json", "w"), indent=0)
