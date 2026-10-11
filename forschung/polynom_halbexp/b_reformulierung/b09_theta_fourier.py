"""b09: Idee 3 numerisch — p_N in Kneser-Koordinaten. θ_N(w) = A(p_N(T(w))) − w − 1/2 auf waagrechten Strecken der Länge 1
im Korridor B, Fourier-Koeffizienten θ̂_n(y) = ∫_0^1 θ_N(w0+t+iy) e^{−2πin(w0+t+iy)} dt.
A = Knesers slog aus 300 Taylor-Koeffizienten (|z| ≤ 1.25), T = A^{-1} per Newton. Kein Kneser-Motor.
Erwartung: θ_N ≈ Σ_{n ungerade} c_n e^{2πinw} (antiperiodisch), c_{±1} dominiert, c_{-n} = conj(c_n)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]]
Ac = [mp.mpf(c) for c in json.load(open(G + "kneser_slog_koeff_300.json"))["A"]]
dAc = [k * Ac[k] for k in range(1, len(Ac))]

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

def A(z):
    if abs(z) > 1.27:
        raise ValueError(f"|z| = {abs(z)} zu gross")
    return horner(Ac, z)

def T(w, z0):
    z = mp.mpc(z0)
    for _ in range(60):
        dz = (horner(Ac, z) - w) / horner(dAc, z)
        z -= dz
        if abs(dz) < mp.mpf(10) ** -35:
            return z
    raise RuntimeError("Newton T")

M = 64
def strecke(y, w0):
    pts = []
    # Startwert: Newton von w = -1 (z = 0) aus entlang einer Strecke fortsetzen
    z = T(mp.mpc(-1, 0), mp.mpc(0))
    wa = mp.mpc(-1, 0)
    ziel = mp.mpc(w0, y)
    for s in range(1, 41):
        z = T(wa + (ziel - wa) * s / 40, z)
    for j in range(M):
        w = mp.mpc(w0 + mp.mpf(j) / M, y)
        z = T(w, z)
        pts.append((w, z))
    return pts

def fourier(vals, pts, n):
    return mp.fsum(v * mp.exp(-2j * mp.pi * n * w) for v, (w, z) in zip(vals, pts)) / M

for y, w0 in ((0.0, -1.55), (0.4, -1.65), (0.8, -1.75)):
    pts = strecke(y, w0)
    zmax = max(abs(z) for w, z in pts)
    th_K = [A(horner(K, z)) - w - mp.mpf(1) / 2 for w, z in pts]
    print(f"\nHöhe y={y}, w0={w0}: max|T(w)| = {mp.nstr(zmax, 4)}, Kneser-Kontrolle max|θ_K| = {mp.nstr(max(abs(t) for t in th_K), 3)}")
    for N in (20, 40, 80, 120, 194):
        try:
            th = [A(horner(P[N], z)) - w - mp.mpf(1) / 2 for w, z in pts]
        except ValueError as e:
            print(f"  N={N}: {e}"); continue
        co = {n: fourier(th, pts, n) for n in (-3, -2, -1, 0, 1, 2, 3)}
        per = max(abs(th[j] + th[(j + M // 2) % M]) for j in range(M // 2))  # Antiperiodizität (nur Kontrolle, Strecke Länge 1)
        print(f"  N={N}: max|θ_N|={mp.nstr(max(abs(t) for t in th), 3)}  ĉ_1={mp.nstr(co[1], 4)}  |ĉ_-1|={mp.nstr(abs(co[-1]), 3)}"
              f"  |ĉ_0|={mp.nstr(abs(co[0]), 3)} |ĉ_±2|={mp.nstr(abs(co[2]), 2)},{mp.nstr(abs(co[-2]), 2)} |ĉ_±3|={mp.nstr(abs(co[3]), 2)},{mp.nstr(abs(co[-3]), 2)}")
