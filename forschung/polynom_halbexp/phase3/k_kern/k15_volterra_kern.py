"""k15 — R2 analytisch: Gs Volterra-Rekursion G22b als Renormierungsgleichung. Ist ihr Kern skalenfrei?

Exakt (bewiesen, BERICHT §14): G22b ⇔ 𝒜 = L⁻¹(w^m)_m mit dem unteren LU-Faktor L von J̃ (γ_nj/h_j = L_nj, Lⱼ e_j = J̃c_j).
In den Köpfen 𝒟_n = q_n𝒜_n/h_n: 𝒟_n(w) = (q_n/h_n)w^n − Σ_{j<n} K_nj 𝒟_j(w), K_nj := (q_n h_j)/(h_n q_j)·L_nj — DERSELBE Kern für alle w.
Test: (i) Grösse und Lage von K_nj; (ii) Skalenfreiheit n·|K_{n,tn}| für n = 150, 200, 300; (iii) Anteil der Quelle; (iv) schnelle Phase von K.
Aufruf: python k15_volterra_kern.py kubisch 300 [|μ|]
"""
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, Lc, RHO  # noqa: E402
from k05_widom import inv  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    name = sys.argv[1]
    Nr = int(sys.argv[2])
    t0 = time.time()
    bits = int(1.2 * Nr) + 400
    mu = None
    if len(sys.argv) > 3:   # optional: |μ| ändern bei festem arg μ (kubisches Modell), z. B. 1.30
        from flint import acb, arb
        from k_kern import Lacb
        ctx.prec = bits
        m0 = Lacb().sqrt()
        mu = m0 * (arb(sys.argv[3]) / abs(m0))
        print("Modell %s mit |μ| = %s (1/|μ| = %.4f)" % (name, sys.argv[3], 1 / float(sys.argv[3])), flush=True)
    J = J_tilde(name, Nr, bits=bits, mu=mu)
    ctx.prec = bits
    X = inv(J).mid()
    Jf = None
    heads_c = [None] * (Nr + 1)
    heads_l = [None] * (Nr + 1)
    piv = [None] * (Nr + 1)
    cur = X
    for n in range(Nr, -1, -1):
        p = cur[n, n]
        c = [cur[i, n] for i in range(n + 1)]
        r = [cur[n, j] for j in range(n + 1)]
        heads_c[n] = c
        heads_l[n] = [v / p for v in r]
        piv[n] = p
        if n == 0:
            break
        nxt = arb_mat(n, n)
        for i in range(n):
            ci = c[i] / p
            for j in range(n):
                nxt[i, j] = cur[i, j] - ci * r[j]
        cur = nxt.mid()
    print("Köpfe: %.0f s" % (time.time() - t0), flush=True)
    # h_n = 1/p_n, Q_n = h_n c_n, q_n = ‖Q_n‖, 𝒜_n(w_L) mit w_L = ρ/L
    wL = RHO / Lc
    h = np.array([float((1 / piv[n]).mid()) for n in range(Nr + 1)])
    q = np.array([np.linalg.norm([float((v / piv[n]).mid()) for v in heads_c[n]]) for n in range(Nr + 1)])
    A_L = np.array([sum(complex(float(heads_l[n][k].mid())) * wL ** k for k in range(n + 1)) for n in range(Nr + 1)])
    D_L = q * A_L / h
    # L_nj = (J̃ c_j)_n für n > j — in arb, dann float
    Lmat = np.zeros((Nr + 1, Nr + 1))
    for j in range(Nr + 1):
        cj = heads_c[j]
        for n in range(j, Nr + 1):
            s = 0
            for k in range(j + 1):
                s += J[n, k] * cj[k]
            Lmat[n, j] = float(s.mid())
    print("L: %.0f s; Diagonale (sollte 1 sein): min %.6f, max %.6f" % (time.time() - t0, Lmat.diagonal().min(), Lmat.diagonal().max()), flush=True)
    K = np.zeros((Nr + 1, Nr + 1))
    for n in range(Nr + 1):
        for j in range(n):
            K[n, j] = q[n] * h[j] / (h[n] * q[j]) * Lmat[n, j]
    # Kontrolle der Rekursion in w_L
    for n in [100, 200, Nr]:
        src = q[n] / h[n] * wL ** n
        rec = src - np.dot(K[n, :n], D_L[:n])
        print("n = %d: 𝒟_n(L) = %.6e%+.6ei, Rekursion %.6e%+.6ei, |Quelle| = %.2e, √n|𝒟_n(L)| = %.5f"
              % (n, D_L[n].real, D_L[n].imag, rec.real, rec.imag, abs(src), np.sqrt(n) * abs(D_L[n])), flush=True)
    # Lage und Skalenfreiheit des Kerns
    for n in [150, 200, 250, Nr]:
        row = np.abs(K[n, :n])
        tot = row.sum()
        jmax = np.argmax(row)
        ts = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]
        vals = ["t=%.2f: %.3e" % (t, n * row[min(int(t * n), n - 1)]) for t in ts]
        near = "  ".join("m=%d: %.3e" % (m, row[n - m]) for m in [1, 2, 3, 5, 10, 20])
        print("n = %d: Σ_j|K_nj| = %.4g, Maximum bei j = %d (t = %.3f); n·|K_{n,tn}|: %s" % (n, tot, jmax, jmax / n, "  ".join(vals)), flush=True)
        print("          nahe Diagonale |K_{n,n−m}|: " + near, flush=True)
    # Phase: arg K_{n,n-m} − arg K_{n-1,n-1-m} (schneller Träger)
    for m in [1, 2, 5]:
        ph = [np.angle(K[n, n - m]) for n in range(Nr - 5, Nr + 1)]
        print("Vorzeichen/Phase K_{n,n−%d}, n = %d…%d: %s" % (m, Nr - 5, Nr, " ".join("%+.3f" % p for p in ph)), flush=True)
    np.savez(__file__.replace(".py", "_%s_%d%s.npz" % (name, Nr, ("_mu" + sys.argv[3]) if len(sys.argv) > 3 else "")), h=h, q=q, D_L=D_L, L=Lmat, K=K)
    print("fertig %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
