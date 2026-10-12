"""k13 — Erster Test zu Idee R12 (Quadranten-Wiener–Hopf): sind die künftigen Köpfe, auf Grade ≤ N projiziert, skalenfrei in (ξ, λ)?

Alle Köpfe bis N' aus EINER Inversen durch Rang-1-Abbau (Satz K5): X^{(n)} = J̃_n⁻¹, c_n = X^{(n)}e_n, ℓ_nᵀ = e_nᵀX^{(n)}/X^{(n)}_nn,
X^{(n−1)} = Block_n(X^{(n)} − c_n ℓ_nᵀ). In w: r̂_nᵀ = ℓ_nᵀC⁻¹. D_N = −Σ_{n>N}(P_Nc_n)(P_Nr̂_n)ᵀ.
Gemessen: e_N(λ) = ‖P_Nc_{λN}‖·‖P_Nr̂_{λN}‖ für λ ∈ (1, 2] und mehrere N — Skalenfreiheit ⇔ N·e_N(λ) fällt auf eine Kurve; dazu das Indexprofil.
Aufruf: python k13_koepfe_2d.py kubisch 300 80,100,120,150
"""
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde  # noqa: E402
import k04_konform as K4  # noqa: E402
from k05_widom import inv  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    name = sys.argv[1]
    Nr = int(sys.argv[2])
    Ns = [int(s) for s in sys.argv[3].split(",")]
    t0 = time.time()
    bits = int(1.2 * Nr) + 400
    J = J_tilde(name, Nr, bits=bits)
    ctx.prec = bits
    X = inv(J).mid()
    F, d = K4.theodorsen()
    Ci = K4.float_mat(inv(K4.C_matrix(F, Nr)))
    print("Inverse und C⁻¹: %.0f s" % (time.time() - t0), flush=True)
    heads_c = [None] * (Nr + 1)
    heads_l = [None] * (Nr + 1)
    cur = X
    for n in range(Nr, -1, -1):
        piv = cur[n, n]
        c = [cur[i, n] for i in range(n + 1)]
        r = [cur[n, j] for j in range(n + 1)]
        heads_c[n] = np.array([float(v.mid()) for v in c])
        heads_l[n] = np.array([float((v / piv).mid()) for v in r])
        if n == 0:
            break
        nxt = arb_mat(n, n)
        for i in range(n):
            ci = c[i] / piv
            for j in range(n):
                nxt[i, j] = cur[i, j] - ci * r[j]
        cur = nxt.mid()   # nur Mittelpunkte weiterführen (Kugelradien würden über 300 Abbauschritte explodieren)
    print("Köpfe abgebaut: %.0f s" % (time.time() - t0), flush=True)
    rhat = [heads_l[n] @ Ci[:n + 1, :n + 1] for n in range(Nr + 1)]
    # Kopf-Daten für spätere Analysen speichern (gepolstert auf Länge Nr+1)
    Cm = np.zeros((Nr + 1, Nr + 1)); Rm = np.zeros((Nr + 1, Nr + 1))
    for n in range(Nr + 1):
        Cm[n, :n + 1] = heads_c[n]; Rm[n, :n + 1] = rhat[n]
    np.savez(__file__.replace(".py", "_%s_%d.npz" % (name, Nr)), c=Cm, r=Rm)
    # Kontrolle: Kopfterm-Normen in x und w bei einigen n
    for n in [100, 200, Nr]:
        if n <= Nr:
            print("n = %d: ‖c_n‖ = %.3e, ‖ℓ_n‖ = %.3e, Kopf x = %.4g, Kopf w = %.4g"
                  % (n, np.linalg.norm(heads_c[n]), np.linalg.norm(heads_l[n]),
                     np.linalg.norm(heads_c[n]) * np.linalg.norm(heads_l[n]), np.linalg.norm(heads_c[n]) * np.linalg.norm(rhat[n])), flush=True)
    lams = [1.05, 1.1, 1.2, 1.3, 1.5, 1.7, 2.0]
    for N in Ns:
        row = []
        for lam in lams:
            n = int(round(lam * N))
            if n > Nr:
                row.append(float("nan"))
                continue
            pc = np.linalg.norm(heads_c[n][:N + 1])
            pr = np.linalg.norm(rhat[n][:N + 1])
            row.append(N * pc * pr)
        print("N = %3d  N·‖P_Nc_{λN}‖·‖P_Nr̂_{λN}‖:  " % N + "  ".join("λ=%.2f: %.4g" % (l, v) for l, v in zip(lams, row)), flush=True)
    # Indexprofil (Schwerpunkt und Maximum in ξ = i/N) der projizierten Köpfe bei λ = 1.2 und 1.5
    for N in Ns:
        out = []
        for lam in [1.2, 1.5]:
            n = int(round(lam * N))
            if n > Nr:
                continue
            v = np.abs(heads_c[n][:N + 1]) ** 2
            u = np.abs(rhat[n][:N + 1]) ** 2
            xi = np.arange(N + 1) / N
            out.append("λ=%.1f: c-Schwerpunkt ξ=%.3f (Max %.3f), r̂-Schwerpunkt ξ=%.3f (Max %.3f)"
                       % (lam, (v * xi).sum() / v.sum(), xi[np.argmax(v)], (u * xi).sum() / u.sum(), xi[np.argmax(u)]))
        print("N = %3d  " % N + " | ".join(out), flush=True)
    # D_N aus den Köpfen (Summe n = N+1 … N') — Singulärwerte zum Vergleich mit k05/k09
    for N in Ns:
        Dm = np.zeros((N + 1, N + 1))
        for n in range(N + 1, Nr + 1):
            Dm -= np.outer(heads_c[n][:N + 1], rhat[n][:N + 1])
        s = np.linalg.svd(Dm, compute_uv=False)
        print("N = %3d  σ(D_N aus Köpfen bis %d) = %s" % (N, Nr, " ".join("%.4f" % x for x in s[:6])), flush=True)
    print("fertig %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
