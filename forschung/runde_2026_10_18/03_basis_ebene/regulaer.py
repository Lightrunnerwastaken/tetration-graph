"""Reguläre Tetration am anziehenden Fixpunkt (Koenigs), unabhängig vom Motor.
Koordinate: Multiplikator λ (|λ|<1), Basis b = exp(λ e^{−λ}), anziehender Fixpunkt a = e^λ (fp.satz.multiplikator_rezept).
g(λ, x) = σ^{-1}(λ^x σ(1)),  σ Koenigs-Funktion mit σ'(a)=1; ausgewertet als log_b^n(a + λ^{x+n}σ(1) + ...)."""
import mpmath as mp

def fix_anz(b):
    lb = mp.log(b)
    best=None
    for k in (0,-1,1):
        L = -mp.lambertw(-lb,k)/lb
        lam = L*lb
        if abs(lam)<1 and (best is None or abs(lam)<abs(best[1])): best=(L,lam)
    return best

def g_lam(lam, x, ordnung=3):
    """T^reg_b(x) mit b = exp(λ e^{−λ})."""
    lam = mp.mpc(lam)
    lb = lam*mp.exp(-lam); a = mp.exp(lam)
    return _T(lb, a, lam, x, ordnung)

def T_reg(b, x, ordnung=3):
    b = mp.mpc(b); lb = mp.log(b)
    a, lam = fix_anz(b)
    return _T(lb, a, lam, x, ordnung)

def _T(lb, a, lam, x, ordnung):
    # Schröder-Reihe σ(a+u)=u+s2 u^2+s3 u^3: f(a+u)=a+λ(u+lb u^2/2+lb^2 u^3/6)
    f1, f2, f3 = lam, lam*lb/2, lam*lb**2/6
    s2 = f2/(lam - lam**2)
    s3 = (f3 + 2*s2*f1*f2)/(lam - lam**3)
    tol = mp.mpf(10)**(-mp.mp.dps)
    m = 0; w = mp.mpc(1)
    # bis |u|^4 < tol
    while abs(w-a)**4 > tol:
        w = mp.exp(lb*w); m += 1
        if m > 10**6: raise RuntimeError("keine Konvergenz")
    u = w-a; sig1 = (u + s2*u**2 + s3*u**3)/lam**m
    n = m
    s = mp.exp((x+n)*mp.log(lam))*sig1
    # Umkehrreihe: u = s - s2 s^2 + (2 s2^2 - s3) s^3
    w = a + s - s2*s**2 + (2*s2**2 - s3)*s**3
    for _ in range(n):
        w = mp.log(w)/lb
    return w
