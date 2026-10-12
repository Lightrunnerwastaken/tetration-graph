"""k59c — Liegt der dominante Eigenvektor der linearisierten Dyson-Abbildung (Variante c-rü, k59) in der Tangente der m-Familie?
Die Riccati-Gleichung [S, X] = X E X hängt nur von f ab; jedes X_m = (C_f + M_m)⁻¹ löst sie. Tangente: δX = −X·T_k·X (T_k = S̃^k, Toeplitz),
in die Eichung (letzte Zeile 0) gebracht durch Abzug der unteren Toeplitz-Matrix mit derselben letzten Zeile.
Ausgabe: Anteil der dominanten Eigenvektoren (|λ| ≥ 0.9·ρ) im Spann dieser n+1 Richtungen (kleinste Quadrate).
Aufruf: python k59c_m_familie.py 30 40
"""
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, modell  # noqa: E402
from k59_dyson_spektrum import to_np, S_inv_rueck  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def toep_letzte_zeile(row):
    """Untere Dreiecks-Toeplitz-Matrix T mit vorgegebener letzter Zeile (T_{n,j} = t_{n−j})."""
    n1 = len(row); n = n1 - 1
    t = [row[n - k] for k in range(n1)]      # t_k = T_{n, n−k}
    T = np.zeros((n1, n1))
    for i in range(n1):
        for j in range(i + 1):
            T[i, j] = t[i - j]
    return T


def lauf(n):
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    Xr = Jr.inv()
    rho_a = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    S = arb_mat(n + 1, n + 1)
    for k in range(1, n + 1):
        S[k, k - 1] = rho_a
    Er = Jr * S - S * Jr
    fk, _ = modell("kubisch", n)
    f = arb_poly(fk[:4])
    pw = arb_poly([1])
    for _ in range(n + 1):
        pw = arb_poly([pw[k] if k <= pw.degree() else arb(0) for k in range(n + 1)]) * f
    for k in range(n + 1):
        ckk = pw[k] if k <= pw.degree() else arb(0)
        Er[k, n] = Er[k, n] + ckk * rho_a ** k / rho_a ** n
    rho = float(rho_a.mid())
    X = to_np(Xr); Ec = to_np(Er); Sn = to_np(S)
    A = Ec @ X; B = X @ Ec
    idx = [(i, j) for j in range(n + 1) for i in range(n)]
    M = np.zeros((len(idx), len(idx)))
    for k, (i, j) in enumerate(idx):
        R = np.zeros((n + 1, n + 1)); R[i, :] += A[j, :]; R[:, j] += B[:, i]
        Y = S_inv_rueck(R, rho)
        M[:, k] = [Y[p] for p in idx]
    ev, W = np.linalg.eig(M)
    o = np.argsort(-np.abs(ev)); ev = ev[o]; W = W[:, o]
    # Tangente der m-Familie in Eichung
    basis = []
    Tk = np.eye(n + 1)
    for k in range(n + 1):
        d = -X @ Tk @ X
        d = d - toep_letzte_zeile(d[n, :])
        basis.append([d[p] for p in idx])
        Tk = Tk @ (Sn / rho)
    Bm = np.array(basis).T
    Q, _ = np.linalg.qr(Bm)
    rho_L = abs(ev[0])
    print("n = %d: ρ(L) = %.4f" % (n, rho_L))
    for q in range(len(ev)):
        if abs(ev[q]) < 0.9 * rho_L or q >= 12:
            break
        w = W[:, q]
        anteil = np.linalg.norm(Q.conj().T @ w) / np.linalg.norm(w)
        print("   λ = %.4f%+.4fi: Anteil in der m-Tangente = %.3f" % (ev[q].real, ev[q].imag, anteil))
    # Kontrolle: Wirkung von L auf die m-Tangente (sollte ≈ Identität modulo Eichung sein)
    w = Bm[:, 1] / np.linalg.norm(Bm[:, 1])
    Lw = M @ w
    print("   Kontrolle: ‖L·t₁ − t₁‖/‖t₁‖ = %.3e (t₁ = m-Tangente zu x¹)" % (np.linalg.norm(Lw - w)))


def main():
    for a in sys.argv[1:]:
        lauf(int(a))


if __name__ == "__main__":
    main()
