"""Runde 3, Schritt 1: J_N als Resolvente eines „komprimierten Drehoperators“.
Cayley-Bild: E(λw) + λE(w) = O((w−1)^{N+1}) ⇔ j_λE (in Potenzen von w/λ − 1) = −λ·j_1E (in Potenzen von w − 1).
Für Δ ∈ 𝒫_N bestimmt j_1Δ das Element; 𝒮_N := (j_1Δ ↦ j_λΔ) (beide Jets in den gedrehten Koordinaten).
Dann gilt exakt: J_N ≅ 𝒮_N + λ·I (bis auf invertierbare Dreiecks-Faktoren). Frage: wo liegt das Spektrum von 𝒮_N relativ zu −λ?"""
import sys
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 80


def S_matrix(c, N):
    lam = mp.expj(c)
    n = N + 1
    A = mp.matrix(n, n)
    B = mp.matrix(n, n)
    for j in range(n):
        D = lambda w, j=j: (1 + lam * w) ** 2 * w ** j / (1 + w) ** N
        a = mp.taylor(lambda t: D(1 + t), 0, N)
        b = mp.taylor(lambda t: D(lam * (1 + t)), 0, N)
        for k in range(n):
            A[k, j] = a[k]
            B[k, j] = b[k]
    return lam, B * mp.inverse(A)


if __name__ == "__main__":
    cs = [mp.mpf(x) for x in sys.argv[1].split(",")]
    Ns = [int(x) for x in sys.argv[2].split(",")]
    for c in cs:
        for N in Ns:
            lam, S = S_matrix(c, N)
            ev = mp.eig(S, left=False, right=False)
            ev = sorted(ev, key=lambda z: mp.arg(z))
            d = min(abs(e + lam) for e in ev)
            mods = [abs(e) for e in ev]
            print(f"c={mp.nstr(c, 3)} N={N}: |Eigenwerte| ∈ [{mp.nstr(min(mods), 4)}, {mp.nstr(max(mods), 4)}], "
                  f"arg-Bereich [{mp.nstr(mp.arg(ev[0]) * 180 / mp.pi, 4)}°, {mp.nstr(mp.arg(ev[-1]) * 180 / mp.pi, 4)}°], "
                  f"−λ bei {mp.nstr(mp.arg(-lam) * 180 / mp.pi, 4)}°, dist(−λ, Spektrum) = {mp.nstr(d, 3)}", flush=True)
            print("   Eigenwerte (|μ|, arg°): " + ", ".join(f"({mp.nstr(abs(e), 3)},{mp.nstr(mp.arg(e) * 180 / mp.pi, 4)})" for e in ev), flush=True)
