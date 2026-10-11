"""Strang A, Aufgabe 3a/4: Wie sieht d_N = p_N − q_N aus? Zerlegung in abgeschnittene Störrichtungen δ_n^{(N)}.

Für N = 20..194: (1) Profil |d_{N,k}|·|L|^k über k; (2) Kleinste-Quadrate-Anpassung d_N ≈ Σ_n (α_n Re δ_n + β_n Im δ_n)_{≤N}
(n = 1; 1,3; 1,3,5) mit Gewichten |L|^k bzw. 1.3^k; relativer Rest; (3) Werte bei x = 1: d_N(1) gegen Modell, und
p_N(1) − f_K(1) = d_N(1) − r_N(1) (r_N = Kneser-Schwanz).
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
mp.mp.dps = 70
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)

q = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
P = {int(k): [mp.mpf(x) for x in v] for k, v in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8")).items()}
D = json.loads((GEM / "stoerrichtungen_delta.json").read_text(encoding="utf-8"))["delta"]
delta = {int(n): [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in zip(v["re"], v["im"])] for n, v in D.items()}


def lstsq(spalten, ziel, w):
    """min |Σ c_i spalten_i − ziel|_2 mit Gewichten w (Normalgleichungen in mpmath, kleine Dimension)."""
    m = len(spalten)
    S = [[x * wk for x, wk in zip(s, w)] for s in spalten]
    z = [x * wk for x, wk in zip(ziel, w)]
    G = mp.matrix(m, m)
    b = mp.matrix(m, 1)
    for i in range(m):
        b[i] = mp.fsum(a * c for a, c in zip(S[i], z))
        for j in range(m):
            G[i, j] = mp.fsum(a * c for a, c in zip(S[i], S[j]))
    c = mp.lu_solve(G, b)
    rest = [zz - mp.fsum(c[i] * S[i][k] for i in range(m)) for k, zz in enumerate(z)]
    return [c[i] for i in range(m)], mp.sqrt(mp.fsum(x * x for x in rest)) / mp.sqrt(mp.fsum(x * x for x in z))


def horner(c, x):
    s = mp.mpf(0)
    for a in reversed(c):
        s = s * x + a
    return s


print("Profil |d_{N,k}|·|L|^k (k = 0, N/4, N/2, 3N/4, N−1, N) und Vergleich mit |a_N|·|L|^N:")
for N in (20, 40, 60, 80, 100, 140, 194):
    d = [P[N][k] - q[k] for k in range(N + 1)]
    ks = [0, N // 4, N // 2, 3 * N // 4, N - 1, N]
    print(f"  N = {N:3d}: " + ", ".join(f"k={k}: {mp.nstr(abs(d[k]) * RL ** k, 3)}" for k in ks)
          + f";  |a_N||L|^N = {mp.nstr(abs(q[N]) * RL ** N, 3)};  a_N^(N)/a_N = {mp.nstr(P[N][N] / q[N], 4)}")

print("\nAnpassung d_N ≈ Σ (α_n Re δ_n + β_n Im δ_n), Gewicht |L|^k (ℓ²); relativer Rest und Koeffizienten:")
zeilen = []
for N in list(range(20, 200, 10)) + [194]:
    if N > 194:
        continue
    d = [P[N][k] - q[k] for k in range(N + 1)]
    w = [RL ** k for k in range(N + 1)]
    out = {}
    for ns in ((1,), (1, 3), (1, 3, 5)):
        sp = []
        for n in ns:
            sp += [[mp.re(x) for x in delta[n][:N + 1]], [mp.im(x) for x in delta[n][:N + 1]]]
        c, rest = lstsq(sp, d, w)
        out[ns] = (c, rest)
    c1, r1 = out[(1,)]
    c3, r3 = out[(1, 3)]
    c5, r5 = out[(1, 3, 5)]
    # Werte bei x = 1
    d1 = horner(d, mp.mpf(1))
    mod1 = c1[0] * mp.re(horner(delta[1][:N + 1], mp.mpf(1))) + c1[1] * mp.im(horner(delta[1][:N + 1], mp.mpf(1)))
    voll1 = c1[0] * mp.re(horner(delta[1], mp.mpf(1))) + c1[1] * mp.im(horner(delta[1], mp.mpf(1)))
    rN1 = mp.fsum(q[N + 1:])
    zeilen.append((N, c1, r1, r3, r5, d1, mod1, voll1, rN1))
    print(f"  N = {N:3d}: Rest(δ1) = {mp.nstr(r1, 3)}, Rest(δ1,δ3) = {mp.nstr(r3, 3)}, Rest(δ1,3,5) = {mp.nstr(r5, 3)};  "
          f"|c_1| = {mp.nstr(mp.sqrt(c1[0] ** 2 + c1[1] ** 2), 4)}, arg = {mp.nstr(mp.atan2(c1[1], c1[0]), 4)};  "
          f"d_N(1) = {mp.nstr(d1, 4)}, Modell(δ1 abgeschn.) = {mp.nstr(mod1, 4)}, Modell(δ1 voll) = {mp.nstr(voll1, 4)}, r_N(1) = {mp.nstr(rN1, 4)}")

print("\nSkalierung |c_1(N)|: lokale Exponenten log(c(N2)/c(N1))/log(N2/N1):")
for (N1, c1a, *_), (N2, c1b, *_) in zip(zeilen[:-1], zeilen[1:]):
    a = mp.sqrt(c1a[0] ** 2 + c1a[1] ** 2)
    b = mp.sqrt(c1b[0] ** 2 + c1b[1] ** 2)
    print(f"  N = {N1}→{N2}: {mp.nstr(mp.log(b / a) / mp.log(mp.mpf(N2) / N1), 4)}")
