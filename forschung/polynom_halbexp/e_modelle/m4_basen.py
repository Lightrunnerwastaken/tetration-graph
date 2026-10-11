"""Modell 4: g = b^x für verschiedene Basen b (auch b < e^{1/e} mit reellen Fixpunkten).
Ausgabe: p_N(0), p_N(1), Selbstkonvergenz |p_N − p_Nmax| bei x = 0 und 1, lokale Raten; Fixpunkt L_b, Multiplikator μ = L_b ln b,
Singularitätsstärke Re(2πi/log μ). Für b < e^{1/e}: Vergleich mit der regulären halben Iterierten am anziehenden Fixpunkt a
(f = g^{-n} ∘ F_a ∘ g^n, F_a = Taylor-Reihe der regulären Wurzel bei a)."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_exp_basis, horner

mp.mp.dps = 80
HIER = Path(__file__).parent
NMAX = int(sys.argv[1])
BS = sys.argv[2].split(",")
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 900
ORD = int(sys.argv[4]) if len(sys.argv) > 4 else 2


def basis(s):
    return {"e": mp.e, "eta": mp.e ** (1 / mp.e), "sqrt2": mp.sqrt(2)}.get(s, None) or mp.mpf(s)


def regulaer_halb(b, a, M=60):
    """Taylor-Reihe der regulären halben Iterierten von b^x am Fixpunkt a (Koeffizienten in (y−a)^k)."""
    lb = mp.log(b)
    G = [a * lb ** k / mp.factorial(k) for k in range(M + 1)]   # b^(a+t) = a·Σ (lb t)^k/k!
    G[0] = 0                                                       # Abweichung g(a+t) − a
    F = [mp.mpf(0), mp.sqrt(G[1])]
    for n in range(2, M + 1):
        # [t^n] F∘F = G_n; F∘F = Σ F_m F(t)^m ; unbekannt F_n tritt als F_n(F_1^n + F_1) auf
        F.append(mp.mpf(0))
        Ft = F[:]
        # Koeffizient n von Σ_m F_m (F(t))^m mit F_n = 0
        pot = [mp.mpf(1)] + [mp.mpf(0)] * n
        s = mp.mpf(0)
        Fs = Ft[:n + 1]
        pw = [mp.mpf(1)] + [mp.mpf(0)] * n
        for m in range(1, n + 1):
            neu = [mp.fsum(pw[i] * Fs[j - i] for i in range(j + 1)) for j in range(n + 1)]
            pw = neu
            s += Fs[m] * pw[n]
        F[n] = (G[n] - s) / (F[1] ** n + F[1])
    return F


tabelle = {}
for bs in BS:
    b = basis(bs)
    lb = mp.log(b)
    t0 = time.time()
    try:
        loes = zweig(g_exp_basis(b), NMAX, bits=BITS, ordnung=ORD)
    except RuntimeError as e:
        print("b =", bs, e)
        continue
    Nm = max(loes)
    w0 = {N: loes[N][0] for N in loes}
    w1 = {N: horner(loes[N], mp.mpf(1)) for N in loes}
    # Fixpunkt(e)
    L = -mp.lambertw(-lb) / lb
    L = mp.conj(L) if mp.im(L) < 0 else L
    mu = L * lb
    staerke = mp.re(2j * mp.pi / mp.log(mu)) if abs(mp.im(L)) > 1e-30 else None
    zeile = {"b": mp.nstr(b, 12), "L": mp.nstr(L, 12), "|L|": mp.nstr(abs(L), 8), "mu": mp.nstr(mu, 8),
             "p(0)": mp.nstr(w0[Nm], 30), "p(1)": mp.nstr(w1[Nm], 30)}
    print(f"\n== b = {bs} ({mp.nstr(b, 10)}), Ordnung {ORD}: L_b = {mp.nstr(L, 10)}, |L_b| = {mp.nstr(abs(L), 8)}, μ = {mp.nstr(mu, 8)}"
          + (f", Re(2πi/log μ) = {mp.nstr(staerke, 6)}" if staerke is not None else "") + f"  ({time.time() - t0:.0f}s)")
    print(f"   p_Nmax(0) = {mp.nstr(w0[Nm], 25)}, p_Nmax(1) = {mp.nstr(w1[Nm], 25)}")
    if abs(mp.im(L)) < 1e-30 and ORD == 2:   # reelle Fixpunkte: anziehender = untere Lösung von b^x = x
        a = mp.re(L)
        F = regulaer_halb(b, a, 60)
        def freg(x, n=60):
            y = x
            for _ in range(n):
                y = b ** y
            z = a + mp.fsum(F[k] * (y - a) ** k for k in range(1, len(F)))
            for _ in range(n):
                z = mp.log(z) / lb
            return z
        r0, r1 = freg(mp.mpf(0)), freg(mp.mpf(1))
        print(f"   reguläre Wurzel am anziehenden Fixpunkt a = {mp.nstr(a, 12)}: f(0) = {mp.nstr(r0, 25)}, f(1) = {mp.nstr(r1, 25)}")
        zeile["freg(0)"], zeile["freg(1)"] = mp.nstr(r0, 30), mp.nstr(r1, 30)
        ref0, ref1 = r0, r1
    else:
        ref0, ref1 = w0[Nm], w1[Nm]
    for N in (1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160):
        if N in loes and N <= Nm:
            print(f"   N={N:3d}  p(0)-ref={mp.nstr(w0[N] - ref0, 4):>11}  p(1)-ref={mp.nstr(w1[N] - ref1, 4):>11}")
    zeile["d0"] = {N: mp.nstr(w0[N] - ref0, 5) for N in loes}
    zeile["d1"] = {N: mp.nstr(w1[N] - ref1, 5) for N in loes}
    tabelle[bs] = zeile
json.dump(tabelle, open(HIER / "daten" / f"m4_basen_ord{ORD}_{'_'.join(BS)}.json", "w"), indent=0)
