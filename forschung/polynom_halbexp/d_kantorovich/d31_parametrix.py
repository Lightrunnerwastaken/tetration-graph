"""D31 — Wo sitzt Y − T_N(χ) (Y = (I + K_N)⁻¹ = J_N⁻¹ M_N, χ = 1_Γin)? Eckschicht (Flip-Teil W_N X W_N) oder Inneres?
Kubisches Modell. Frobenius-Anteile in Blöcken; Vergleich der Eckblöcke in Flip-Koordinaten für zwei N."""
import sys
import numpy as np
from flint import arb, arb_mat, ctx
import d_kern as K
sys.stdout.reconfigure(encoding="utf-8")
L = complex(0.31813150520476413, 1.3372357014306895); R = abs(L); mu = np.sqrt(L)
rL = (mu - 1) / (L - L.conjugate()); beta = rL.imag / L.imag; alpha = rL.real - beta * L.real
coef = [R**2*alpha, 1 + R**2*beta - 2*L.real*alpha, alpha - 2*L.real*beta, beta]
f = np.polynomial.Polynomial(coef); m = f.deriv()(f)
Mq = 1 << 15; phi = 2*np.pi*np.arange(Mq)/Mq; x = R*np.exp(1j*phi)
chi = np.fft.fft((np.abs(f(x)) < R).astype(float)) / Mq
ecken = {}
for N in [int(s) for s in sys.argv[1].split(",")]:
    n = N + 1; ctx.prec = int(0.9*N) + 300
    a = [arb(float(c)) for c in coef] + [arb(0)]*(n - 4)
    _, J = K.FJ(a, N)
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n); sc = R**k
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * sc[:, None] / sc[None, :]
    Mt = np.zeros((n, n))
    for d in range(min(n, len(m.coef))):
        Mt += np.diag(np.full(n - d, m.coef[d]*R**d), -d)
    Y = G @ Mt
    U, S, Vt = np.linalg.svd(Y)
    Yr = Y - S[0]*np.outer(U[:, 0], Vt[0]) - S[1]*np.outer(U[:, 1], Vt[1])     # Eckmoden abgezogen
    T = np.array([[chi[(i - j) % Mq].real for j in range(n)] for i in range(n)])
    D = Yr - T
    fro = (D**2).sum()
    for B in (20, 50, 100):
        ecke = (D[N - B:, :]**2).sum() + (D[:, N - B:]**2).sum() - (D[N - B:, N - B:]**2).sum()
        print(f"N={N}: Anteil der Frobenius-Masse von (Y − Eckmoden − T_N(χ)) in den letzten {B} Zeilen/Spalten: {ecke/fro:.3f}")
    sD = np.linalg.svd(D, compute_uv=False)
    sI = np.linalg.svd(D[:N - 100, :N - 100], compute_uv=False)
    print(f"   σ(D) = {np.round(sD[:6], 4)}, ‖D_innen (ohne letzte 100)‖ = {sI[0]:.4f}, Frobenius² gesamt {fro:.2f}", flush=True)
    ecken[N] = D[N - 40:, N - 40:][::-1, ::-1]
Ns = sorted(ecken)
if len(Ns) > 1:
    a_, b_ = ecken[Ns[-2]], ecken[Ns[-1]]
    print(f"Eckblock 40×40 in Flip-Koordinaten: ‖D(N={Ns[-2]}) − D(N={Ns[-1]})‖_F / ‖D(N={Ns[-1]})‖_F = {np.linalg.norm(a_ - b_)/np.linalg.norm(b_):.3f}")
