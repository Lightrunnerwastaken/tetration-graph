"""Modell 4b: g = b^x (b = 1.5, 2, e, 10) und f∘f∘f = exp gegen Kneser (Werte von Strang A, gemeinsam/kneser_andere_basen.json).
Fehler p_N(0) − sexp_b(−1/(ord)·…), p_N(1) − sexp_b(1/ord); Fit der Abklingrate (geometrisch vs. polynomial)."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_exp_basis, horner

mp.mp.dps = 60
HIER = Path(__file__).parent
K = json.loads((HIER.parent / "gemeinsam" / "kneser_andere_basen.json").read_text(encoding="utf-8"))["werte"]
NMAX = int(sys.argv[1])
FAELLE = sys.argv[2].split(",")          # z.B. 1.5,2,10,e,e3
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 900
aus = {}
for fall in FAELLE:
    t0 = time.time()
    if fall == "e3":
        b, ordn, k0, k1 = mp.e, 3, "sexp_exp(1)(-0.666667)", "sexp_exp(1)(0.333333)"
    else:
        b = mp.e if fall == "e" else mp.mpf(fall)
        ordn = 2
        nm = "exp(1)" if fall == "e" else fall
        k0, k1 = f"sexp_{nm}(-0.5)", f"sexp_{nm}(0.5)"
    r0, r1 = mp.mpf(K[k0]), mp.mpf(K[k1])
    lb = mp.log(b)
    L = -mp.lambertw(-lb) / lb
    L = mp.conj(L) if mp.im(L) < 0 else L
    mu = L * lb
    s = mp.re(2j * mp.pi / mp.log(mu))
    loes = zweig(g_exp_basis(b), NMAX, bits=BITS, ordnung=ordn)
    d0 = {N: loes[N][0] - r0 for N in loes}
    d1 = {N: horner(loes[N], mp.mpf(1)) - r1 for N in loes}
    print(f"\n== {fall} (Ordnung {ordn}): |L_b| = {mp.nstr(abs(L), 6)}, arg L_b = {mp.nstr(mp.arg(L) * 180 / mp.pi, 5)}°, "
          f"μ = L ln b = {mp.nstr(mu, 6)}, Re(2πi/log μ) = {mp.nstr(s, 5)}, f(0) = {mp.nstr(r0, 8)}, 1/|L_b| = {mp.nstr(1 / abs(L), 5)}, "
          f"f(0)/|L_b| = {mp.nstr(r0 / abs(L), 5)}  ({time.time() - t0:.0f}s)")
    for N in (5, 10, 20, 30, 40, 50, 60, 80, 100, 120, 140, 160):
        if N in loes:
            print(f"   N={N:3d}  p(0)-f(0)={mp.nstr(d0[N], 4):>11}  p(1)-f(1)={mp.nstr(d1[N], 4):>11}")
    hull = lambda d, N: max(abs(d[n]) for n in range(N - 4, N + 1))
    Nm = max(loes)
    for N1, N2 in ((20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 140), (140, 160)):
        if N2 <= Nm:
            q = hull(d0, N2) / hull(d0, N1)
            print(f"   Einhüllende |p(0)-f(0)| {N1}->{N2}: Rate/Grad {mp.nstr(q ** (mp.mpf(1) / (N2 - N1)), 4)}, "
                  f"Exponent (falls N^-β) {mp.nstr(-mp.log(q) / mp.log(mp.mpf(N2) / N1), 4)}")
    aus[fall] = {"d0": {N: mp.nstr(d0[N], 5) for N in loes}, "d1": {N: mp.nstr(d1[N], 5) for N in loes},
                 "L": mp.nstr(L, 15), "staerke": mp.nstr(s, 8)}
json.dump(aus, open(HIER / "daten" / f"m4b_{'_'.join(FAELLE)}_N{NMAX}.json", "w"), indent=0)
