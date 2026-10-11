"""Strang P, Prüfung 25 (Runde 3, Auftrag 1, Baustein K1 von Satz P-V): ZERTIFIZIERTE Poincaré-Funktion an L.

P ganz, P(Lw) = e^{P(w)}, P(0) = L, P'(0) = 1  (L = e^L der obere Fixpunkt; e^{P} = L·exp(P − L)).
(1) Koeffizienten p_n (n ≤ n₀) exakt per Rekursion p_n(L^n − L) = L·[wⁿ]exp(P − L)|_{p_n = 0}  (arb/acb).
(2) Rest R = P − T_{n₀} auf |w| ≤ r₀ als Fixpunkt von 𝒯[R](w) = L·exp(T(w/L) − L + R(w/L)) − T(w) in der Kugel
    {sup_{|w|≤r₀} |R(w)|/|w|^{n₀+1} ≤ C}:  𝒯[R] = δ + L·e^{T(w/L)−L}(e^{R(w/L)} − 1),  δ := L·e^{T(w/L)−L} − T (Defekt, O(w^{n₀+1})).
    Selbstabbildung: D + M₁·C·|L|^{−n₀−1}·e^{C(r₀/|L|)^{n₀+1}} ≤ C; Kontraktion: M₁|L|^{−n₀−1}e^{…} < 1. Der Fixpunkt löst die Funktional-
    gleichung mit P(0) = L, P'(0) = 1 und ist daher (formale Eindeutigkeit, L^n ≠ L) die Poincaré-Funktion.
(3) Für |w| > r₀: P(w) = exp^k(P(w/L^k)) in Intervallarithmetik.
Ausgabe: n₀, r₀, C (Restkonstante), Prüfwerte. Aufruf: python p25_poincare.py
"""
import json
import sys
from pathlib import Path

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
ctx.prec = 400
N0 = 80
R0 = arb("0.5")
M_EXTRA = 240
R1 = arb(2)


def exp_reihe(a, n):
    """Koeffizienten von exp(A) bis Ordnung n, A(0) = 0 (Rekursion n e_n = Σ k a_k e_{n−k})."""
    e = [acb(1)] + [acb(0)] * n
    for m in range(1, n + 1):
        s = acb(0)
        for k in range(1, min(m, len(a) - 1) + 1):
            s += k * a[k] * e[m - k]
        e[m] = s / m
    return e


def main():
    L = -acb(-1).lambertw(-1)                      # L = e^L, Im L > 0
    assert abs(L - L.exp()) < arb("1e-100"), L
    absL = abs(L)
    # (1) Koeffizienten
    p = [L, acb(1)] + [acb(0)] * (N0 - 1)
    for n in range(2, N0 + 1):
        a = [acb(0)] + p[1:n] + [acb(0)]            # A = P − L bis Ordnung n, p_n vorläufig 0
        Q = exp_reihe(a, n)[n]
        p[n] = L * Q / (L ** n - L)
    # Prüfung der Rekursion
    a = [acb(0)] + p[1:]
    E = exp_reihe(a, N0)
    fehler = max(float(abs(p[n] * L ** n - L * E[n]).mid()) for n in range(1, N0 + 1))
    # (2) Defekt δ = L·exp(T(w/L) − L) − T(w): Koeffizienten bis N0 + M_EXTRA exakt
    nmax = N0 + M_EXTRA
    b = [acb(0)] + [p[n] / L ** n for n in range(1, N0 + 1)] + [acb(0)] * M_EXTRA
    Eb = exp_reihe(b, nmax)
    delta = [L * Eb[k] - (p[k] if k <= N0 else 0) for k in range(nmax + 1)]
    unten = max(float(abs(delta[k]).mid() + abs(delta[k]).rad()) for k in range(1, N0 + 1))
    # grobe Schranke von |δ| auf |w| = R1 über Majoranten
    S1 = sum((abs(p[n]) * (R1 / absL) ** n for n in range(1, N0 + 1)), arb(0))
    S2 = absL + sum((abs(p[n]) * R1 ** n for n in range(1, N0 + 1)), arb(0))
    B1 = absL * S1.exp() + S2
    D = sum((abs(delta[k]) * R0 ** (k - N0 - 1) for k in range(N0 + 1, nmax + 1)), arb(0))
    q = R0 / R1
    D += B1 * R1 ** (-(N0 + 1)) * q ** (M_EXTRA) / (1 - q)      # Σ_{k>nmax} B1 R1^{-k} r0^{k-N0-1}
    # M1 = sup_{|w|≤r0} |L e^{T(w/L) − L}|
    S1r = sum((abs(p[n]) * (R0 / absL) ** n for n in range(1, N0 + 1)), arb(0))
    M1 = absL * S1r.exp()
    # C aus D + M1 C |L|^{-N0-1} e^{C (r0/|L|)^{N0+1}} ≤ C
    fakt = M1 * absL ** (-(N0 + 1))
    C = D / (1 - fakt) * 2
    for _ in range(20):
        lhs = D + fakt * C * (C * (R0 / absL) ** (N0 + 1)).exp()
        if lhs < C:
            break
        C *= 2
    lip = fakt * (C * (R0 / absL) ** (N0 + 1)).exp()
    ok = bool(lhs < C) and bool(lip < 1)
    print(f"L = {L.str(20)}; Rekursionsfehler {fehler:.2e}; Defekt-Koeff. ≤ N0 (sollen 0 sein): {unten:.2e}")
    print(f"n₀ = {N0}, r₀ = 0.5: D = {float(D.mid()):.3e}, M₁ = {float(M1.mid()):.3f}, C = {float(C.mid()):.3e}, Lipschitz = {float(lip.mid()):.3e}, ok = {ok}")
    print(f"   ⇒ |P(w) − T_{N0}(w)| ≤ C|w|^{N0 + 1}; bei |w| = 0.5: {float((C * R0 ** (N0 + 1)).mid()):.3e}")

    def P_ball(w):
        k = 0
        while abs(w) > R0 * arb("0.999"):
            w = w / L
            k += 1
        s = acb(0)
        for c in reversed(p):
            s = s * w + c
        s += acb(0, 0) + (C * abs(w) ** (N0 + 1)) * acb(arb(0, 1), arb(0, 1))
        for _ in range(k):
            s = s.exp()
        return s

    # Prüfwerte: Funktionalgleichung an Punkten mit |w| > r0, und P(Ψ)-Probe über log-Iteration
    pruef = []
    for w in [acb("0.3", "0.2"), acb("1.5", "-0.7"), acb("-2", "1"), acb(3, 0)]:
        lhs_ = P_ball(L * w)
        rhs_ = P_ball(w).exp()
        pruef.append((str(w), float(abs(lhs_ - rhs_).mid() + abs(lhs_ - rhs_).rad()), float(abs(lhs_).mid()),
                      float(lhs_.real.rad() + lhs_.imag.rad())))
    for w, d, betrag, rad in pruef:
        print(f"   w = {w}: |P(Lw) − e^{{P(w)}}| ≤ {d:.2e}  (|P(Lw)| = {betrag:.4f}, Radius {rad:.1e})")
    (HIER / "p25_poincare.json").write_text(json.dumps({
        "n0": N0, "r0": 0.5, "C": float(C.mid() + C.rad()), "ok": ok, "lipschitz": float(lip.mid()),
        "koeff": [[str(c.real.mid()), str(c.imag.mid())] for c in p[:N0 + 1]],
        "koeff_radius": max(float(c.real.rad() + c.imag.rad()) for c in p)}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
