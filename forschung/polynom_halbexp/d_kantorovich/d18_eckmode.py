"""D18 — Eckmode bis N = 1000: J_N am Punkt x = Π_N f_K (Strang A, gemeinsam/lin_PiN_fK_bis_1000.json; ≈ p_N bis auf
2. Ordnung). Rescaled J̃_kj = |L|^{k−j} J_kj, Inverse wahlweise in arb ("arb") oder float64 ("float", nach Test).
Ausgabe je N: Singulärwerte σ1, σ2, σ3 der umskalierten Inversen (ℓ²_{|L|}), Lage/Profil der Singulärvektoren,
Überlapp mit Eingaben (e^x, Kneser-Typ, J∞f_K-Typ), gespeichert in d18_eckmode.json und d18_vek_<N>.npz.
"""
import json
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
MODUS = sys.argv[2] if len(sys.argv) > 2 else "float"
PI = json.load(open(K.GEMEINSAM + "/lin_PiN_fK_bis_1000.json"))
L = 1.3745570107437075
datei = K.HIER + f"/d18_eckmode_{MODUS}.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
mp.mp.dps = 40
Lc = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")


def eingabe_kneser_typ(n, sig=mp.mpc("5.4468", "1.0578")):
    c, out = mp.mpc(1), []
    for k in range(n + 1):
        out.append(float(2 * mp.re(c)) * L ** k)
        c = c * (k - sig) / (k + 1) / Lc
    return np.array(out)


t0 = time.time()
for N in NS:
    a_str = PI[str(N)]["b"][:N + 1]
    ctx.prec = int(0.72 * N) + 250
    a = [arb(s) for s in a_str]
    _, J = K.FJ(a, N)
    n = N + 1
    k = np.arange(n)
    if MODUS == "arb":
        ctx.prec = int(0.95 * N) + 300
        Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
        G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    else:
        # umskaliert in arb, dann float64
        Lp = [arb(L) ** i for i in range(n)]
        Jt = np.array([[float((J[i, j] * Lp[i] / Lp[j]).mid()) for j in range(n)] for i in range(n)])
        G = np.linalg.inv(Jt)
    U, S, Vt = np.linalg.svd(G)
    ex = np.array([float(1 / mp.factorial(i)) for i in range(n)]) * L ** k
    kt = eingabe_kneser_typ(N)
    z = {"s1": float(S[0]), "s2": float(S[1]), "s3": float(S[2]), "s4": float(S[3])}
    for name, v in (("ex", ex), ("kneser", kt)):
        vn = v / np.linalg.norm(v)
        z[f"ov1_{name}"] = float(abs(Vt[0] @ vn))
        z[f"ov2_{name}"] = float(abs(Vt[1] @ vn))
        z[f"kopf_{name}"] = float(np.linalg.norm(vn[N - 30:]))
    # Profil der Eingabe-Singulärvektoren (Abstand zur Spitze)
    for i in (0, 1):
        prof = np.abs(Vt[i])
        z[f"v{i+1}_masse_letzte30"] = float(np.sum(prof[N - 30:] ** 2))
        z[f"v{i+1}_profil"] = [float(prof[N - d]) for d in (0, 1, 2, 5, 10, 20, 40, 80) if N - d >= 0]
    # Ausgabe-Vektor auf dem Kreis
    M = 8192
    Fu = np.abs(np.fft.ifft(np.concatenate([U[:, 0], np.zeros(M - n)])) * M)
    z["u1_phi"] = 360.0 * int(np.argmax(Fu)) / M
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    np.savez_compressed(K.HIER + f"/d18_vek_{MODUS}_{N}.npz", u1=U[:, 0], u2=U[:, 1], v1=Vt[0], v2=Vt[1], S=S[:6])
    print(f"N={N} ({time.time()-t0:.0f}s) σ = {S[0]:.4f}, {S[1]:.4f}, {S[2]:.4f}, {S[3]:.4f} | Überlapp v1·e^x {z['ov1_ex']:.2e}, "
          f"v1·Kneser {z['ov1_kneser']:.2e} (Kopf {z['kopf_kneser']:.2e}) | v1 letzte30 {z['v1_masse_letzte30']:.3f} | u1 bei {z['u1_phi']:.2f}°",
          flush=True)
