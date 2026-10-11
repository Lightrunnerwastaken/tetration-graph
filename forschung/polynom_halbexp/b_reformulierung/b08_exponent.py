"""b08: Abklingexponent von p_N - f_K (Kreis |x|=1, Punkt x=1, Randpunkt L) für N bis 194.
Vorhersage aus dem Miller-Bild: Kneser-Schwanz ~ k^-6.447, Störrichtung δ_{±1} ~ k^{+2.447} ⇒ c_N ~ N^-8.894,
im Inneren |p_N - f_K| ~ N^-8.894, am Rand L ~ N^-8.894 · N^3.447 = N^-5.447."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s
pts = [mp.expjpi(2 * mp.mpf(j) / 48) for j in range(48)]
ref = [horner(K, z) for z in pts]
Ns = list(range(20, 195, 10)) + [194]
res = {}
for N in Ns:
    m = max(abs(horner(P[N], z) - r) for z, r in zip(pts, ref))
    e1 = horner(P[N], 1) - horner(K, 1)
    eL = horner(P[N], L) - L
    res[N] = (m, e1, eL)
prev = None
for N in Ns:
    m, e1, eL = res[N]
    s = ""
    if prev:
        Np, (mp_, e1p, eLp) = prev
        g = lambda a, b: mp.log(a / b) / mp.log(mp.mpf(N) / Np)
        s = f" lokale Exponenten: Kreis {mp.nstr(g(m, mp_), 3)}, L {mp.nstr(g(abs(eL), abs(eLp)), 3)}"
    print(f"N={N}: max_|x|=1 |p_N-f_K| = {mp.nstr(m,3)}, p_N(1)-f_K(1) = {mp.nstr(mp.re(e1),3)}, |p_N(L)-L| = {mp.nstr(abs(eL),3)}" + s)
    prev = (N, res[N])
# globaler Fit über N=100..194
import math
xs = [math.log(N) for N in Ns if N >= 80]
for idx, name in ((0, "Kreis"), (2, "L")):
    ys = [math.log(float(abs(res[N][idx]))) for N in Ns if N >= 80]
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    print(f"Fit N>=80 {name}: Exponent {sl:.3f}")
