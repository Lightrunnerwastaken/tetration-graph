"""Runde 2, Schritt 1: Möbius-Drehung in Cayley-Koordinaten w = (1+ix)/(1−ix). Dort ist f_M: w ↦ λw (λ = e^{ic}),
J_∞ ≅ 𝓛Δ(w) = Δ(λw) + λΔ(w), P_N ≅ 𝒫_N = {(i/2)(1+λw)² b(w)/(1+w)^N : deg b ≤ N}, Abschneiden = N-Jet bei w = 1,
und f_M selbst ≅ U(w) = (λ²w² − 1)/2.  Π_N f_M ≅ Δ ∈ 𝒫_N mit 𝓛(Δ − U) = O((w−1)^{N+1}).
Teil A: numerische Gegenprobe gegen die x-Formulierung.  Teil B: symbolisch (sympy) für kleine N — Struktur von b."""
import sys
import mpmath as mp
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 50


def pi_N_w(N, c):
    """Löse für b (Koeffizienten) numerisch; gib Δ als Funktion zurück."""
    lam = mp.expj(c)
    # Unbekannte b_0..b_N; Bedingung: Taylor-Koeffizienten (in t = w − 1) von 𝓛(Δ − U) bis Ordnung N verschwinden
    M = N + 1
    def taylor_of(fun):
        return mp.taylor(fun, 1, N)
    A = mp.matrix(M, M)
    for j in range(M):
        def Dj(w, j=j):
            return 1j / 2 * (1 + lam * w) ** 2 * w ** j / (1 + w) ** N
        col = taylor_of(lambda w: Dj(lam * w) + lam * Dj(w))
        for k in range(M):
            A[k, j] = col[k]
    U = lambda w: (lam ** 2 * w ** 2 - 1) / 2
    rhs = taylor_of(lambda w: U(lam * w) + lam * U(w))
    b = mp.lu_solve(A, mp.matrix(rhs))
    return lam, [b[j] for j in range(M)]


def x_wert(lam, b, N, x):
    w = (1 + 1j * x) / (1 - 1j * x)
    Delta = 1j / 2 * (1 + lam * w) ** 2 * sum(b[j] * w ** j for j in range(N + 1)) / (1 + w) ** N
    phi_strich = lambda v: 2 / (1j * (1 + v) ** 2)
    return phi_strich(lam * w) * Delta          # δ(x) = φ'(λw)Δ(w)


if __name__ == "__main__":
    # Teil A: Vergleich mit x-Formulierung (Π_N f_M direkt über J_N an f_M)
    from flint import arb, ctx
    from kern import F_J
    c = mp.mpf("0.8")
    tau = mp.tan(c / 2)
    for N in (4, 8, 12):
        lam, b = pi_N_w(N, c)
        # x-Formulierung: Π_N f_M = J_N^{-1} T_N J_∞ f_M, J_∞ f_M = f∘f + (f'∘f) f = g + (f'∘f) f
        ctx.prec = 300
        f = lambda y: (y + tau) / (1 - tau * y)
        fk = mp.taylor(f, 0, N)
        gk = mp.taylor(lambda y: f(f(y)), 0, N)
        a = [arb(mp.nstr(v, 60)) for v in fk]
        _, J = F_J(a, [arb(mp.nstr(v, 60)) for v in gk], N)      # J an T_N f (≈ J an f bis Ordnung N: J hängt nur von T_N f ab)
        rhs = mp.taylor(lambda y: f(f(y)) + (1 + tau ** 2) / (1 - tau * f(y)) ** 2 * f(y), 0, N)
        import flint
        R = flint.arb_mat([[arb(mp.nstr(v, 60))] for v in rhs])
        d = J.solve(R)
        px = [mp.mpf(d[k, 0].mid().str(40, radius=False)) for k in range(N + 1)]
        for x in (mp.mpf("0.3"), mp.mpf("-0.5")):
            v1 = mp.polyval(px[::-1], x)
            v2 = x_wert(lam, b, N, x)
            print(f"N={N} x={x}: x-Form {mp.nstr(v1, 15)}  w-Form {mp.nstr(v2, 15)}  f_M(x) {mp.nstr(f(x), 15)}")
