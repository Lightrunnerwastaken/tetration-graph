"""Modell 1d: Konvergenzrate p_N → f_M (elliptische Möbius-Abbildung) als Funktion von c.
Fehlerfolge e_N = p_N(x0) − f_M(x0) und d_N = a_0^(N) − τ; Prony-Anpassung e_N ≈ Σ A_j z_j^N (zwei bis vier Moden)."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_moebius_ell, horner

mp.mp.dps = 80
HIER = Path(__file__).parent
NMAX = int(sys.argv[1])
CS = [mp.mpf(x) for x in sys.argv[2].split(",")]
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 900
aus = {}


def prony(e, n0, m):
    """Löse e_{n+m} = Σ_k c_k e_{n+k} (n = n0..n0+m−1) und gib die Wurzeln des Polynoms zurück."""
    A = mp.matrix(m, m)
    b = mp.matrix(m, 1)
    for i in range(m):
        for k in range(m):
            A[i, k] = e[n0 + i + k]
        b[i] = e[n0 + i + m]
    co = mp.lu_solve(A, b)
    poly = [1] + [-co[m - 1 - k] for k in range(m)]
    return sorted(mp.polyroots(poly, maxsteps=200, extraprec=200), key=lambda z: -abs(z))


for c in CS:
    t0 = time.time()
    tau = mp.tan(c / 2)
    try:
        loes = zweig(g_moebius_ell(c), NMAX, bits=BITS)
    except RuntimeError as err:
        print("c =", c, err)
        continue
    d0 = {N: loes[N][0] - tau for N in loes}
    fM = lambda z: (z + tau) / (1 - tau * z)
    e5 = {N: horner(loes[N], mp.mpf("0.5")) - fM(mp.mpf("0.5")) for N in loes}
    Nmax = max(loes)
    zeile = {"tau": mp.nstr(tau, 10), "Nmax": Nmax, "zeit": round(time.time() - t0)}
    for m in (2, 4):
        n0 = Nmax - 2 * m - 2
        try:
            z = prony(d0, n0, m)
            zeile[f"prony{m}_a0"] = [(mp.nstr(abs(w), 5), mp.nstr(mp.arg(w) * 180 / mp.pi, 4)) for w in z]
        except Exception as err:
            zeile[f"prony{m}_a0"] = str(err)
    zeile["e_a0"] = {N: mp.nstr(d0[N], 4) for N in loes}
    zeile["e_x05"] = {N: mp.nstr(e5[N], 4) for N in loes}
    # gemittelte Rate der Einhüllenden über die letzten 20 Grade
    hull = lambda d, N: max(abs(d[n]) for n in range(N - 5, N + 1))
    zeile["rate20"] = mp.nstr((hull(d0, Nmax) / hull(d0, Nmax - 20)) ** (mp.mpf(1) / 20), 5)
    aus[mp.nstr(c, 6)] = zeile
    print(f"c={mp.nstr(c, 4)} tau={mp.nstr(tau, 5)} Nmax={Nmax} rate20={zeile['rate20']} prony2={zeile['prony2_a0']} prony4={zeile['prony4_a0']} ({zeile['zeit']}s)", flush=True)
json.dump(aus, open(HIER / "daten" / f"m1d_raten_{sys.argv[2].replace(',', '_')}.json", "w"), indent=0)
