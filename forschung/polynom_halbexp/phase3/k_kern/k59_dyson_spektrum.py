"""k59 — Dyson-Abbildung aus (CD1) linearisiert: Spektralradius (BERICHT §41). Kubisch, skaliert (J̃ = DJD⁻¹, S̃ = ρS).
Riccati [S̃, X] = X Ẽ X (X = J̃_n⁻¹, Ẽ = J̃S̃ − S̃J̃); Eichung: erste Spalte von X gegeben.
𝒮₀⁻¹: [S̃, Y] = R ⇔ Y_{i,j+1} = Y_{i−1,j} − R_{ij}/ρ (j = 0 … n−1), Y_{·,0} = 0.
L[δX] = 𝒮₀⁻¹(δX·E·X + X·E·δX) auf {δX : δX e_0 = 0};  Variante „c“: E → Ẽ_c (nur Kompositionskanal, Rang-1-Teil −(T_nf^{n+1})e_nᵀ entfernt).
Kontrollen: Residuum der Riccati-Gleichung; Rekonstruktion X = 𝒮_ν⁻¹(XẼX) aus der ersten Spalte.
Aufruf: python k59_dyson_spektrum.py 10 20 30 40
"""
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, modell  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def to_np(A):
    return np.array([[float(A[i, j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])


def S_inv(R, rho, Y0=None):
    n1 = R.shape[0]
    Y = np.zeros_like(R)
    if Y0 is not None:
        Y[:, 0] = Y0
    for j in range(n1 - 1):
        Y[1:, j + 1] = Y[:-1, j] - R[1:, j] / rho
        Y[0, j + 1] = -R[0, j] / rho
    return Y


def S_inv_rueck(R, rho, Yn=None):
    """[S̃, Y] = R rückwärts: Y_{i−1,j} = Y_{i,j+1} + R_{ij}/ρ, letzte Zeile gegeben (0), Y_{·,n+1} = 0."""
    n1 = R.shape[0]
    Y = np.zeros_like(R)
    if Yn is not None:
        Y[n1 - 1, :] = Yn
    for i in range(n1 - 1, 0, -1):
        Y[i - 1, :-1] = Y[i, 1:] + R[i, :-1] / rho
        Y[i - 1, -1] = R[i, -1] / rho
    return Y


def lauf(n):
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    Xr = Jr.inv()
    rho_a = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    S = arb_mat(n + 1, n + 1)
    for k in range(1, n + 1):
        S[k, k - 1] = rho_a
    Er = Jr * S - S * Jr
    # Rang-1-Teil: −D T_n f^{n+1} ρ^{−n} e_nᵀ  (skaliert)
    fk, _ = modell("kubisch", n)
    f = arb_poly(fk[:4])
    pw = arb_poly([1])
    for _ in range(n + 1):
        pw = arb_poly([pw[k] if k <= pw.degree() else arb(0) for k in range(n + 1)]) * f
    Ecr = arb_mat(Er)
    for k in range(n + 1):
        ckk = pw[k] if k <= pw.degree() else arb(0)
        Ecr[k, n] = Ecr[k, n] + ckk * rho_a ** k / rho_a ** n
    rho = float(rho_a.mid())
    X = to_np(Xr); E = to_np(Er); Ec = to_np(Ecr)
    # Kontrollen in hoher Genauigkeit
    res = Xr * Er * Xr - (S * Xr - Xr * S)
    rmax = max(abs(float(res[i, j].mid())) for i in range(n + 1) for j in range(n + 1))
    xmax = max(abs(float(Xr[i, j].mid())) for i in range(n + 1) for j in range(n + 1))
    Yrec = S_inv(X @ E @ X, rho, X[:, 0])
    rec = np.max(np.abs(Yrec - X)) / np.max(np.abs(X))
    Yrec2 = S_inv_rueck(X @ E @ X, rho, X[n, :])
    rec2 = np.max(np.abs(Yrec2 - X)) / np.max(np.abs(X))
    out = []
    varianten = (("voll", E, False), ("c", Ec, False), ("c-rü", Ec, True), ("v-rü", E, True))
    wahl = os.environ.get("K59_VARIANTEN")
    if wahl:
        varianten = tuple(v for v in varianten if v[0] in wahl.split(","))
    for name, EE, rueck in varianten:
        A = EE @ X          # δX·A
        B = X @ EE          # B·δX
        if rueck:
            idx = [(i, j) for j in range(n + 1) for i in range(n)]
        else:
            idx = [(i, j) for j in range(1, n + 1) for i in range(n + 1)]
        pos = {p: k for k, p in enumerate(idx)}
        M = np.zeros((len(idx), len(idx)))
        for k, (i, j) in enumerate(idx):
            R = np.zeros((n + 1, n + 1))
            R[i, :] += A[j, :]
            R[:, j] += B[:, i]
            Y = S_inv_rueck(R, rho) if rueck else S_inv(R, rho)
            M[:, k] = [Y[p] for p in idx]
        ev = np.linalg.eigvals(M)
        ev = ev[np.argsort(-np.abs(ev))]
        nrm = np.linalg.norm(M, 2)
        out.append((name, ev[:6], nrm, np.sum(np.abs(ev) > 1)))
    print("n = %d (%.0f s): Riccati-Residuum %.1e (max|X| %.1e), Rekonstruktion aus erster Spalte %.1e, aus letzter Zeile %.1e" % (n, time.time() - t0, rmax, xmax, rec, rec2))
    for name, ev, nrm, cnt in out:
        print("   %-4s ρ(L) = %.4g, ‖L‖₂ = %.3g, #|λ|>1 = %d; grösste: %s" % (name, abs(ev[0]), nrm, cnt, ", ".join("%.4g%+.4gi" % (z.real, z.imag) for z in ev[:6])))


def main():
    for a in sys.argv[1:]:
        lauf(int(a))


if __name__ == "__main__":
    main()
