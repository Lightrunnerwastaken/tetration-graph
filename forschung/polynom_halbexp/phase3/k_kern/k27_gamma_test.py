"""k27 — Ist das Eckprofil eine Mellin–Barnes-Funktion vom Malyuzhinets-Typ?

Kandidat (BERICHT §18): y(ζ) = Σ_k Γ(β − iαk)·(−λζ)^k / k!,  α = log μ/(2π) = L/(4π)  (μ = √L, log μ = L/2).
Begründung: Mellin y = (1/2πi)∫Γ(s)Γ(β + iαs)λ^{−s}ζ^{−s}ds; die Pole von Γ(β + iαs) rechts sind genau die Nullstellen von 1 + μ^{s−2}
(für β = ½ − 2iα), also hat H = μy + μ^{−1}y(·/μ) keine Potenzasymptotik — der adjungierte Kopf ist im homogenen Keil ~ ζ^{−s₀}.
Skaleninvarianter Test: r_k = k c_k/c_{k−1} = −λ Γ(β − iαk)/Γ(β − iα(k−1)),  q_k := r_k/r_{k−1} hängt nur von β ab.
Daten: Gs normierte Profile S_n (spiral_analysis.json, nur gelesen), Gs extrapoliertes S_∞, Os Ψ_Y/Ψ_T-Taylor (o39_*.npz, nur gelesen).
Aufruf: python k27_gamma_test.py
"""
import json
import os
import sys

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
HIER = os.path.dirname(os.path.abspath(__file__))
WURZEL = os.path.join(HIER, "..", "..")
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)


def lade_G():
    d = json.load(open(os.path.join(WURZEL, "g_gpt_hilfssaetze", "runde5_szego", "spiral_analysis.json"), encoding="utf-8"))
    out = {}
    for p in d["models"]["kubisch"]["profiles"]:
        c = [complex(float(a), float(b)) for a, b in p["natural_shape_coefficients_first_40"]]
        out[p["n"]] = np.array(c)
    return out


def q_aus_koeff(c):
    r = np.array([k * c[k] / c[k - 1] for k in range(1, len(c))])  # r[k-1] = r_k
    q = r[1:] / r[:-1]  # q[k-2] = q_k, k >= 2
    return r, q


def q_modell(beta, kmax, vorz=-1):
    """q_k für c_k ∝ Γ(β + vorz·iαk)/k!."""
    out = []
    for k in range(2, kmax + 1):
        g = lambda j: mp.loggamma(beta + vorz * 1j * ALPHA * j)
        out.append(complex(mp.exp(g(k) + g(k - 2) - 2 * g(k - 1))))
    return np.array(out)


def fit_beta(q, kmin, kmax, vorz):
    """Kleinste Quadrate über β (komplex) für q_k, k in [kmin, kmax]."""
    ks = range(kmin, kmax + 1)

    def res(b):
        qm = q_modell(b, kmax, vorz)
        return [qm[k - 2] - q[k - 2] for k in ks]

    best = None
    for b0 in [0.5 - 2j * complex(ALPHA), 0.5, 1.0, 0.5 + 2j * complex(ALPHA), 0.2 + 0.3j, 1.5 - 0.5j]:
        b = complex(b0)
        for _ in range(40):
            r0 = np.array(res(b))
            h = 1e-6
            J1 = (np.array(res(b + h)) - r0) / h
            # komplexe Ableitung (holomorph in b)
            step = -np.vdot(J1, r0) / np.vdot(J1, J1)
            b = b + step
            if abs(step) < 1e-12:
                break
        err = np.linalg.norm(res(b)) / np.linalg.norm(q[kmin - 2:kmax - 1])
        if best is None or err < best[1]:
            best = (b, err)
    return best


def main():
    G = lade_G()
    ns = sorted(G)
    print("α = L/(4π) = %s;  Kandidat β₀ = ½ − 2iα = %s" % (mp.nstr(ALPHA, 10), mp.nstr(0.5 - 2j * ALPHA, 10)))
    for n in [ns[0], 200, 400, 600, 800, ns[-1]]:
        if n not in G:
            continue
        c = G[n]
        r, q = q_aus_koeff(c)
        print("\nn = %d: r_k = k c_k/c_(k-1) (k = 1..8): %s" % (n, " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in r[:8])))
        print("        q_k = r_k/r_(k-1) (k = 2..10): %s" % " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in q[:9]))
        for vorz in (-1, +1):
            for (k0, k1) in [(2, 12), (10, 30), (20, 38)]:
                b, err = fit_beta(q, k0, k1, vorz)
                print("        Familie Γ(β %s iαk), k ∈ [%d,%d]: β = %.5f%+.5fi, rel. Fehler %.2e" % ("−" if vorz < 0 else "+", k0, k1, b.real, b.imag, err))
    # Modell-q für β₀ zum Vergleich
    qm = q_modell(0.5 - 2j * ALPHA, 12, -1)
    print("\nModell β₀ = ½ − 2iα, Familie −: q_k (k = 2..10): %s" % " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in qm[:9]))
    # O: Ψ_Y, Ψ_T Taylor
    for fn in ["o39_kubisch_800.npz", "o39_kubisch_400.npz"]:
        p = os.path.join(WURZEL, "o_operator", fn)
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        print("\n%s: Felder %s" % (fn, list(z.keys())))
        if "tay" in z:
            t = z["tay"]
            print("   tay shape", t.shape, t.dtype)


if __name__ == "__main__":
    main()
