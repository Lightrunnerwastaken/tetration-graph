"""b10: Amplitude c_N des Modus e^{2πiw} von θ_N (Kneser-Koordinaten) für viele N; Test des Gesetzes
c_N ≈ C · N^{1 − 2 conj(γ)},  γ = 2πi/L = 4.447 + 1.058i.  (Strecke y = 0, w0 = −1.55, M = 64 Punkte.)"""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
Ac = [mp.mpf(c) for c in json.load(open(G + "kneser_slog_koeff_300.json"))["A"]]
dAc = [k * Ac[k] for k in range(1, len(Ac))]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
gam = 2j * mp.pi / L

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

def T(w, z0):
    z = mp.mpc(z0)
    for _ in range(60):
        dz = (horner(Ac, z) - w) / horner(dAc, z)
        z -= dz
        if abs(dz) < mp.mpf(10) ** -36:
            return z
    raise RuntimeError

M, w0 = 64, mp.mpf("-1.55")
pts, z = [], mp.mpf(0)
for j in range(M):
    w = w0 + mp.mpf(j) / M
    z = T(w, z if j else T(w, mp.mpf(-0.5)))
    pts.append((w, z))
ex = [mp.exp(-2j * mp.pi * w) for w, _ in pts]
out = {}
print("γ =", mp.nstr(gam, 8), "  1 − 2·conj(γ) =", mp.nstr(1 - 2 * mp.conj(gam), 8))
for N in list(range(10, 195, 4)) + [194]:
    th = [horner(Ac, horner(P[N], z)) - w - mp.mpf(1) / 2 for w, z in pts]
    c1 = mp.fsum(t * e for t, e in zip(th, ex)) / M
    C = c1 / mp.power(N, 1 - 2 * mp.conj(gam))
    out[N] = (c1, C)
    print(f"N={N:3d}: c_N = {mp.nstr(c1, 6):>32}  |c_N| = {mp.nstr(abs(c1), 4):>10}  C_N = c_N/N^(1-2γ̄) = {mp.nstr(C, 6)}")
json.dump({N: [str(mp.re(c)), str(mp.im(c))] for N, (c, C) in out.items()}, open("b10_amplitude.json", "w"))
