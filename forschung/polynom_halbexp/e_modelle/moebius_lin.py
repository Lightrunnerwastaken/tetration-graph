"""Runde 2 — Werkzeuge für das lineare Problem im Möbius-Drehungsmodell (Fixpunkte ±i, symmetrisch).
f(x) = (x + τ)/(1 − τx), τ = tan(c/2); m = f'∘f = (1+τ²)(1−τx)²/(1−τ²−2τx)²; J_∞u = u∘f + m·u.
J_N (an der EXAKTEN Lösung f) als arb-Matrix; Π_N u = J_N⁻¹ T_N J_∞ u; Pol-Antwort E_N(ζ, z) = Π_N k_ζ(z) − k_ζ(z)."""
import mpmath as mp
from flint import arb, acb, arb_mat, acb_mat, arb_poly, ctx


def setup(c, N, bits=800):
    ctx.prec = bits
    c = arb(mp.nstr(c, 60)) if not isinstance(c, arb) else c
    tau = (c / 2).tan()
    # Taylor-Koeffizienten von f (bis N)
    fk = [tau] + [(1 + tau * tau) * tau ** (k - 1) for k in range(1, N + 1)]
    f = arb_poly(fk)
    # Potenzen von f (abgeschnitten): T_N(f^j)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    # m = (1+τ²)(1−τx)²·(1 − τ² − 2τx)^{−2}
    a = 1 - tau * tau
    inv = [(2 * tau / a) ** k / a for k in range(N + 1)]            # 1/(a − 2τx) = Σ (2τ/a)^k x^k / a
    inv = arb_poly(inv)
    m = ((1 + tau * tau) * arb_poly([1, -tau]) * arb_poly([1, -tau]) * inv * inv).truncate(N + 1)
    J = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = P[j] + (m * arb_poly([0] * j + [1])).truncate(N + 1)
        for k in range(N + 1):
            J[k, j] = col[k]
    return {"c": c, "tau": tau, "N": N, "P": P, "m": m, "J": J, "f": f}


def rhs_pol(S, zeta):
    """T_N J_∞ k_ζ für komplexes ζ (acb): Koeffizienten von 1/(ζ − f(x)) + m(x)/(ζ − x) bis N."""
    N, P, m = S["N"], S["P"], S["m"]
    zeta = acb(zeta)
    a0 = S["tau"]
    # 1/(ζ − f(x)): f(x) = a0 + h(x), h = O(x): 1/(ζ − a0 − h) = Σ_n h^n/(ζ − a0)^{n+1}
    h = S["f"] - arb_poly([a0])
    hp = [arb_poly([1])]
    for _ in range(N):
        hp.append((hp[-1] * h).truncate(N + 1))
    v = [acb(0)] * (N + 1)
    for n in range(N + 1):
        fac = 1 / (zeta - a0) ** (n + 1)
        for k in range(N + 1):
            v[k] += hp[n][k] * fac
    # m(x)/(ζ − x) = m(x) Σ x^k/ζ^{k+1}
    for k in range(N + 1):
        s = acb(0)
        for i in range(k + 1):
            s += m[i] / zeta ** (k - i + 1)
        v[k] += s
    return v


def pol_antwort(S, zeta, Jinv=None):
    """Koeffizienten von Π_N k_ζ (acb-Liste)."""
    N = S["N"]
    v = rhs_pol(S, zeta)
    if Jinv is None:
        Jinv = S["J"].inv()
        S["Jinv"] = Jinv
    out = []
    for k in range(N + 1):
        s = acb(0)
        for j in range(N + 1):
            s += Jinv[k, j] * v[j]
        out.append(s)
    return out


def auswerten(coef, z):
    z = acb(z)
    s = acb(0)
    for cf in reversed(coef):
        s = s * z + cf
    return s


def setup_allg(A, B, c, N, bits=800):
    """Allgemeine elliptische Drehung um A ± iB (Winkel c): f(x) = A + B·tan(arctan((x−A)/B) + c/2) als Möbius-Abbildung
    f = (αx + β)/(γx + δ); J_N an der exakten Lösung."""
    ctx.prec = bits
    A, B = arb(mp.nstr(A, 60)), arb(mp.nstr(B, 60))
    h = (arb(mp.nstr(c, 60)) / 2).tan()
    # u = (x−A)/B, tan(arctan u + c/2) = (u + h)/(1 − h u)  ⇒  f(x) = A + B(x − A + Bh)/(B − h(x − A))
    al, be = B - A * h + h * A, (B * (B * h - A) + A * B) - 0      # Zähler: A(B − h(x−A)) + B(x − A + Bh)
    # direkt: Zähler = A·B − A h x + A² h + B x − A B + B² h = (B − A h) x + h (A² + B²);  Nenner = (B + h A) − h x
    al, be = B - A * h, h * (A * A + B * B)
    ga, de = -h, B + h * A
    a0 = be / de
    # Taylor von f: (αx + β)/(δ + γx) = (αx + β)/δ · Σ (−γx/δ)^k
    geo = [(-ga / de) ** k / de for k in range(N + 1)]
    fk = [be * geo[0]] + [be * geo[k] + al * geo[k - 1] for k in range(1, N + 1)]
    f = arb_poly(fk)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    # m = f'(f(x)) = (αδ − βγ)/(γ f(x) + δ)²;  γf + δ = (γ a0 + δ) + γ (f − a0)
    det = al * de - be * ga
    q0 = ga * a0 + de
    hser = f - arb_poly([a0])
    inv, hp = arb_poly([0]), arb_poly([1])
    for n in range(N + 1):
        inv += hp * ((-ga) ** n / q0 ** (n + 1))
        hp = (hp * hser).truncate(N + 1)
    inv = inv.truncate(N + 1)
    m = (det * inv * inv).truncate(N + 1)
    J = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = P[j] + (m * arb_poly([0] * j + [1])).truncate(N + 1)
        for k in range(N + 1):
            J[k, j] = col[k]
    S = {"c": c, "tau": a0, "N": N, "P": P, "m": m, "J": J, "f": f, "a0": a0, "pol": de / h}
    return S
