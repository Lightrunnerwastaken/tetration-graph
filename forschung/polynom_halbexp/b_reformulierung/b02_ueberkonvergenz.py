"""b02: Approximiert p_N die Fortsetzung von f_K auch ausserhalb von D (z. B. reell bei x > |L|)?
f_K(x) für x > |L| über f_K(x) = exp(f_K(log x)) (f(e^y) = e^{f(y)}), f_K(log x) aus der Taylorreihe (|log x| < 1)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_80.json"))]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

def fK(z):
    # Taylorreihe innerhalb |z| <= 1.0 gut; sonst über exp∘f∘log bzw. log∘f∘exp
    if abs(z) <= 1.0:
        return horner(K, z)
    return mp.exp(fK(mp.log(z)))

for x in (1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, -1.2, -1.5, -2.0):
    x = mp.mpf(x)
    try:
        ref = fK(x) if x > 0 else None
    except Exception:
        ref = None
    row = []
    for N in (20, 40, 80, 120, 160, 194):
        v = horner(P[N], x)
        row.append(mp.nstr(mp.re(v - ref), 3) if ref is not None else mp.nstr(mp.re(v), 6))
    print(f"x={mp.nstr(x,3)}: f_K={mp.nstr(ref,10) if ref is not None else '-'}  p_N - f_K (N=20,40,80,120,160,194):", row)

# max |p_N| auf Kreisen r > |L|, verglichen mit (r/|L|)^N
for r in (1.5, 2.0, 3.0):
    row = []
    for N in (20, 40, 80, 160):
        m = max(abs(horner(P[N], r * mp.expjpi(2 * mp.mpf(j) / 256))) for j in range(256))
        row.append(mp.nstr(m, 4))
    print(f"r={r}: max|p_N| (N=20,40,80,160):", row, " (r/|L|)^N:", [mp.nstr((r / abs(L)) ** N, 3) for N in (20, 40, 80, 160)])
