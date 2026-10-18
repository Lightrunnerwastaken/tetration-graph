"""E4: Taylor-Koeffizienten a_n = (1/n!) d^n/db^n T_b(x) bei b = η, aus der regulären Familie:
Cauchy-Koeffizienten um b0 = η − s (Kreis q·s, ganz im Shell–Thron-Inneren), dann Umentwicklung nach (b − η).
Vergleich mit Kneser-Differenzenquotienten von oben (Motor).  Aufruf: python e4_ableitungen_eta.py DPS M s q x"""
import sys, json, mpmath as mp
from regulaer import T_reg
from motor import HIER
dps = int(sys.argv[1]); M = int(sys.argv[2]); mp.mp.dps = dps
s = mp.mpf(sys.argv[3]); q = mp.mpf(sys.argv[4]); x = mp.mpf(sys.argv[5])
NMAX = int(sys.argv[6]) if len(sys.argv) > 6 else 40
eta = mp.exp(mp.exp(-1)); b0 = eta - s; rho = q * s
v = [T_reg(b0 + rho * mp.expjpi(2 * mp.mpf(j) / M), x) for j in range(M)]
c = [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * k) / M) for j in range(M)) / M / rho**k for k in range(M // 2)]
# Umentwicklung: T = Σ c_k (b−b0)^k = Σ c_k (u + s)^k, u = b−η
a = []
for n in range(int(sys.argv[7]) if len(sys.argv)>7 else 9):
    a.append(mp.fsum(c[k] * mp.binomial(k, n) * s**(k - n) for k in range(n, NMAX)))
print(f"x={x} s={s} q={q} M={M} NMAX={NMAX}")
for n, an in enumerate(a):
    print(n, mp.nstr(mp.re(an), 22), mp.nstr(mp.im(an), 3))
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / f"e4_x{sys.argv[5]}_s{sys.argv[3]}_{dps}_{M}_{NMAX}.json").write_text(
    json.dumps([mp.nstr(mp.re(z), 30) for z in a]))
