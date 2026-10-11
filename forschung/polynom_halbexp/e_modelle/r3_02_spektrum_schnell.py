"""Runde 3, Schritt 2: Spektrum des komprimierten Drehoperators λ⁻¹𝒮_N (J_N ≅ λ(I + λ⁻¹𝒮_N)) — schnell und streng (acb).
Basis von 𝒫_N: Δ_j(w) = (1+λw)² w^j/(1+w)^N; A = Jets bei w = 1 (Potenzen von t = w − 1), B = Jets bei λ in gedrehten Koordinaten
(Δ_j(λ(1+t))). 𝒮 = B A⁻¹.  Ausgabe: nicht-reelle bzw. nicht-positive Eigenwerte von λ⁻¹𝒮 (ausser dem Paar nahe λ^{∓1})."""
import sys
import mpmath as mp
from flint import acb, acb_mat, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")


def binom_reihe(a, z, n):
    """Koeffizienten von (1 + z t)^a bis t^n (a ganz, auch negativ)."""
    out, c = [acb(1)], acb(1)
    for k in range(1, n + 1):
        c = c * (a - k + 1) / k * z
        out.append(c)
    return out


def mult(p, q, n):
    return [sum((p[i] * q[k - i] for i in range(k + 1)), acb(0)) for k in range(n + 1)]


def S_matrix(c, N, bits):
    ctx.prec = bits
    lam = (acb(0, 1) * arb(mp.nstr(c, 50))).exp()
    n = N
    # Jets bei 1: Δ_j(1+t) = (1+λ+λt)² (1+t)^j (2+t)^{−N}
    q1 = [(1 + lam) ** 2, 2 * lam * (1 + lam), lam ** 2] + [acb(0)] * n
    inv2 = [x * arb(2) ** (-N) for x in binom_reihe(-N, acb("0.5"), n)]
    base1 = mult(q1[:n + 1], inv2, n)
    # Jets bei λ (gedreht): Δ_j(λ(1+t)) = (1+λ²+λ²t)² λ^j (1+t)^j (1+λ+λt)^{−N}
    q2 = [(1 + lam ** 2) ** 2, 2 * lam ** 2 * (1 + lam ** 2), lam ** 4] + [acb(0)] * n
    inv3 = [x * (1 + lam) ** (-N) for x in binom_reihe(-N, lam / (1 + lam), n)]
    base2 = mult(q2[:n + 1], inv3, n)
    A = acb_mat(n + 1, n + 1)
    B = acb_mat(n + 1, n + 1)
    for j in range(n + 1):
        pj = binom_reihe(j, acb(1), n)
        a = mult(base1, pj, n)
        b = mult(base2, pj, n)
        for k in range(n + 1):
            A[k, j] = a[k]
            B[k, j] = b[k] * lam ** j
    return lam, B * A.inv()


if __name__ == "__main__":
    cs = [mp.mpf(x) for x in sys.argv[1].split(",")]
    Ns = [int(x) for x in sys.argv[2].split(",")]
    for c in cs:
        for N in Ns:
            lam, S = S_matrix(c, N, 40 * N + 400)
            T = S * (1 / lam)
            ev = T.eig()
            ev = [mp.mpc(e.real.mid().str(30, radius=False), e.imag.mid().str(30, radius=False)) for e in ev]
            unimod = sorted(ev, key=lambda e: abs(abs(e) - 1))[:2]
            rest = [e for e in ev if e not in unimod]
            schlecht = [e for e in rest if abs(mp.im(e)) > 1e-12 * abs(e) or mp.re(e) <= 0]
            print(f"c={mp.nstr(c, 4)} N={N}: Paar nahe Einheitskreis {[mp.nstr(e, 5) for e in unimod]}; "
                  f"übrige {len(rest)}: reell-positiv bis auf {len(schlecht)}"
                  + (f" — Ausreisser (|μ|, arg°): {[(mp.nstr(abs(e), 3), mp.nstr(mp.arg(e) * 180 / mp.pi, 4)) for e in schlecht[:6]]}" if schlecht else "")
                  + f"; Spanne {mp.nstr(min(abs(e) for e in rest), 3)} … {mp.nstr(max(abs(e) for e in rest), 3)}", flush=True)
