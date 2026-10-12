"""k11 — Familie: ist die ECHTE Newton-Matrix T_N𝒥_P|_{P_N} (P = T_N f_K bzw. P = p_N) in der konformen Koordinate ebenso stabil wie J_N?

𝒥_P δ = δ∘P + (P'∘P)·δ.  Wichtig: P'∘P mod x^{N+1} hängt von ALLEN Koeffizienten von P' ab (P(0) = a₀ ≠ 0), nicht nur von P mod x^N.
Für P = T_N f_K ist P'∘P − f'∘f = −(Σ_{k>N} k a_k x^{k−1})∘f + … : auf Γout (|f| > ρ) exponentiell gross — „Überlauf im Multiplikationskanal“.
Gemessen je N: ‖Ĵ_N⁻¹‖ (konform, Linse S = 0.0166 wie k04), σ₁(J̃_N⁻¹) (x-Koordinate), max_d ρ^d|[x^d](m_P − m)|.
Varianten: lin = J_N (Linearisierung bei f_K, m = e^x/f'), taylor = Newton bei P = T_N f_K, pN = Newton bei P = p_N (gemeinsam/pN_N100_bis_1000.json).
Aufruf: python k11_newton.py 200,300,400 [lin,taylor,pN,PiN,seg1e-09,seg1e-06,seg0.001,seg0.5]
3. Argument (optional): Linsenhöhe S (Standard 0.0166).
seg λ: P = p_N + λ(T_N f_K − p_N) (Strecke Lösung → Taylor-Start); PiN: P = Π_N f_K (gemeinsam/lin_PiN_fK_bis_1000.json, Strang-A-Newton-Start).
"""
import json
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
from flint import arb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import modell, L_RE, L_IM, GEMEINSAM  # noqa: E402
import k04_konform as K4  # noqa: E402
from k05_widom import inv  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def J_aus(Pk, mk, N, bits):
    """J̃ (ρ-skaliert) aus Kompositionspolynom P (Koeff. Pk) und Multiplikator mk; liefert auch die Potenzen."""
    ctx.prec = bits
    f = arb_poly(Pk)
    pw = [arb_poly([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(N + 1))
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    rp = [rho ** k for k in range(N + 1)]
    if mk is None:  # Newton: m_P = Σ_k k P_k P^{k−1} mod x^{N+1}
        m = arb_poly([0])
        for k in range(1, N + 1):
            m += pw[k - 1] * (Pk[k] * k)
        mk = [m[d] if d <= m.degree() else arb(0) for d in range(N + 1)]
    J = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = pw[j]
        for k in range(N + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            J[k, j] = (ck + mm) * rp[k] / rp[j]
    return J, mk, rp


LAMS = [1e-9, 1e-6, 1e-3, 0.5]
WAHL = None


def main():
    global WAHL
    Ns = [int(s) for s in sys.argv[1].split(",")]
    if len(sys.argv) > 2:
        WAHL = sys.argv[2].split(",")
    S = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0166
    F, d = K4.theodorsen(S=S)
    print("Linsenhöhe S = %g, F_0 = %.5f (mittlere Linsenhöhe unter harmonischem Mass)" % (S, F[0]), flush=True)
    phis = np.linspace(-np.pi, np.pi, 4001)
    xs = abs(complex(float(L_RE), float(L_IM))) * np.exp(K4.s_profil(phis, S=S) + 1j * phis)
    pN = json.load(open(os.path.join(GEMEINSAM, "pN_N100_bis_1000.json"), encoding="utf-8"))["p"]
    PiN = json.load(open(os.path.join(GEMEINSAM, "lin_PiN_fK_bis_1000.json"), encoding="utf-8"))
    for N in Ns:
        bits = int(1.2 * N) + 300
        ctx.prec = bits
        C = K4.C_matrix(F, N)
        fk, mk_lin = modell("kneser", N)
        varianten = [("lin", fk, mk_lin), ("taylor", fk, None)]
        if str(N) in pN:
            ctx.prec = bits
            pk = [arb(s) for s in pN[str(N)]]
            varianten.append(("pN", pk, None))
            dk = np.array([float((pk[k] - fk[k]).mid()) for k in range(N + 1)])
            vals = np.polyval(dk[::-1], xs)
            rr = np.abs(xs) <= abs(complex(float(L_RE), float(L_IM))) * (1 + 1e-12)
            print("N = %4d  |p_N − T_N f_K|: max auf ∂D-Teil von γ = %.3e, max auf der Linse = %.3e (bei arg %.1f°)"
                  % (N, np.max(np.abs(vals[rr])), np.max(np.abs(vals[~rr])), np.degrees(phis[~rr][np.argmax(np.abs(vals[~rr]))])), flush=True)
            for lam in LAMS:
                varianten.append(("seg%g" % lam, [pk[k] + (fk[k] - pk[k]) * arb(lam) for k in range(N + 1)], None))
        if str(N) in PiN:
            ctx.prec = bits
            varianten.append(("PiN", [arb(s) for s in PiN[str(N)]["b"]], None))
        varianten = [v for v in varianten if WAHL is None or v[0] in WAHL]
        for name, Pk, mk in varianten:
            t0 = time.time()
            J, mk_used, rp = J_aus(Pk, mk, N, bits)
            dm = [abs(float((mk_used[k] - mk_lin[k]).mid()) * float(rp[k].mid())) for k in range(N + 1)]
            kmax = int(np.argmax(dm))
            Jh = C * J
            X = K4.float_mat(inv(Jh))
            Y = K4.float_mat(inv(J))
            sX = np.linalg.svd(X, compute_uv=False)
            cNN = float(C[N, N].mid())
            sens = cNN * np.linalg.norm(X[:, N]) * np.linalg.norm(X[0, :])
            kopp = cNN * abs(X[0, N])
            sY = np.linalg.svd(Y, compute_uv=False)
            print("N = %4d  %-6s ‖Ĵ⁻¹‖ = %.6f (2. %.4f)   σ₁(J̃⁻¹) = %.5g (2. %.4g)   max_d ρ^d|Δm_d| = %.3e bei d = %d, Δm_N·ρ^N = %.3e  (%.0f s)"
                  % (N, name, sX[0], sX[1], sY[0], sY[1], dm[kmax], kmax, dm[N], time.time() - t0), flush=True)
            print("          Empfindlichkeit der (N,0)-Stelle: s_N = C_NN·‖Ĵ⁻¹e_N‖·‖e_0ᵀĴ⁻¹‖ = %.3e, C_NN·|(Ĵ⁻¹)_0N| = %.3e;  Indikator Δm_N ρ^N · s_N = %.3e"
                  % (sens, kopp, dm[N] * sens), flush=True)


if __name__ == "__main__":
    main()
