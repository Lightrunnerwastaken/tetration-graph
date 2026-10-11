"""Runde 2, Schritt 2: Struktur der Lösung S (Cayley-Koordinate) des linearen Problems Π_N f_M für die Möbius-Drehung.
E = Δ − U = (1+λw) S(w) / (2(1+w)^N), S Polynom vom Grad ≤ N+1 mit S(−1/λ) = 2(1 − 1/λ)^N und
Φ(w) = (1+λ²w)(1+w)^N S(λw) + λ(1+λw)^{N+1} S(w) ≡ 0 mod (w−1)^{N+1}.
Ausgabe: Nullstellen von S (Betrag, Argument), |S| auf dem Einheitskreis, Vergleich mit Π_N f_M − f_M."""
import sys
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 120


def loese_S(N, c):
    lam = mp.expj(c)
    M = N + 2
    A = mp.matrix(M, M)
    rhs = mp.matrix(M, 1)
    # Taylor bei w = 1 von Φ für Basis S = w^j
    for j in range(M):
        col = mp.taylor(lambda w: (1 + lam ** 2 * w) * (1 + w) ** N * (lam * w) ** j + lam * (1 + lam * w) ** (N + 1) * w ** j, 1, N)
        for k in range(N + 1):
            A[k, j] = col[k]
        A[N + 1, j] = (-1 / lam) ** j
    rhs[N + 1] = 2 * (1 - 1 / lam) ** N
    s = mp.lu_solve(A, rhs)
    return lam, [s[j] for j in range(M)]


if __name__ == "__main__":
    c = mp.mpf(sys.argv[1])
    N = int(sys.argv[2])
    lam, s = loese_S(N, c)
    roots = mp.polyroots(s[::-1], maxsteps=400, extraprec=600)
    roots = sorted(roots, key=lambda r: mp.arg(r))
    print(f"c = {c}, N = {N}: λ = e^(ic), w_p = −1/λ (arg {mp.nstr(mp.arg(-1 / lam) * 180 / mp.pi, 5)}°)")
    print("Nullstellen von S (|r|, arg in Grad):")
    print("  " + ", ".join(f"({mp.nstr(abs(r), 4)}, {mp.nstr(mp.arg(r) * 180 / mp.pi, 4)})" for r in roots))
    # |E| auf dem Einheitskreis und in w = ∞-Richtung
    for th in (0, 30, 60, 90, 120, 150, 179):
        w = mp.expjpi(mp.mpf(th) / 180)
        E = (1 + lam * w) * mp.polyval(s[::-1], w) / (2 * (1 + w) ** N)
        print(f"  w = e^(i{th}°): |E| = {mp.nstr(abs(E), 4)}")
