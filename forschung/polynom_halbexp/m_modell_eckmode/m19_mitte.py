"""M19 — Mittelstück streng: c_m ≤ λ_min(I + Re K_m) auf L²([−T2, T2])² (beide Kanten), radiales Modell, Streifen.
Verfahren (rechnergestützt, arb + verifizierte Gleitkomma-Cholesky):
 1. Paneele der Länge h, je n Gauss-Legendre-Knoten (arb.legendre_p_root). Endlicher Rang: k_N(t,t') = Σ k(t_i,t'_j) ℓ_i(t) ℓ_j(t')
    (Tensor-Lagrange je Paneelpaar). In der Orthonormalbasis ℓ_i/√ω_i ist K_N die Matrix 𝐀_ij = √ω_i k(t_i,t'_j) √ω_j (exakt).
 2. 𝐀 in arb (Kugeln), B = I + Re 𝐀 als Gleitkommamatrix + Radien; λ_min(B) ≥ s − γ‖R̂‖²_F, falls Cholesky von B − sI gelingt
    (Higham, Satz 10.3: R̂*R̂ = B − sI + ΔB, |ΔB| ≤ γ_{N+1}|R̂*||R̂|; komplex mit Faktor 4).
 3. ‖K_m − K_N‖ ≤ ‖k − k_N‖_HS; je Paneelpaar sup|k − k_N| ≤ C_int (M_t' + Λ_n M_t) (Hermite-Restformel auf der Bernstein-Ellipse
    E_ρ, Lebesgue-Konstante Λ_n), M_t' = sup |k| für t ∈ P reell, t' auf ∂E(Q); M_t = sup |k^c| für t auf ∂E(P), t' ∈ Q reell
    (k^c: analytische Fortsetzung in t über P0 = (1/2πi)[φ(Φ−t') − φ(Φ*−t')], φ(z) = 1/(1−e^z), Φ*(t) = conj Φ(conj t)).
 Ergebnis: c_m ≥ min(1, λ_min(B) − ‖B − B_wahr‖) − ‖k − k_N‖_HS.
Aufruf: python m19_mitte.py T2 h n rho [nell nsub fern]
"""
import sys
import time

import numpy as np
from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 120
T2 = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(13)
h = arb(sys.argv[2]) if len(sys.argv) > 2 else arb("0.5")
n = int(sys.argv[3]) if len(sys.argv) > 3 else 16
rho = arb(sys.argv[4]) if len(sys.argv) > 4 else arb("2.8")
nell = int(sys.argv[5]) if len(sys.argv) > 5 else 24
nsub = int(sys.argv[6]) if len(sys.argv) > 6 else 4
fern = int(sys.argv[7]) if len(sys.argv) > 7 else 12    # Paneelabstand, ab dem grob geschätzt wird

pi = arb.pi()
eps = arb("0.31382")
th = arb("76.618") * pi / 180
c0 = 1 - 2 * eps * th.cos()
L = acb(th.cos(), th.sin())
Lb = L.conjugate()
sinth, costh = th.sin(), th.cos()
t0 = time.time()


# Formen mit EINMALIGEM Auftreten der Variablen (gegen Abhängigkeitsverlust der Kugelarithmetik)
def f(z):
    return eps * (z + c0 / (2 * eps)) ** 2 + eps - c0 * c0 / (4 * eps)


def fp(z):
    return c0 + 2 * eps * z


def fpf(x):      # f'(f(x)) = c0 + 2ε f(x)
    return c0 + 2 * eps * f(x)


def xi(M):
    return Lb + (Lb - L) / (M - 1)


def hinv(x):
    return (-c0 + (c0 * c0 - 4 * eps * eps + 4 * eps * x).sqrt()) / (2 * eps)


def Lam(x):
    return (1 + eps * (L - Lb) / (c0 + eps * (x + Lb))).log()


def Phi(t, e):
    """Φ_e(t) für reelles oder komplexes t (acb)."""
    t = acb(t)
    M = L * t.exp() if e == 0 else -L * t.exp()
    x = xi(M)
    if e == 0:
        return t + Lam(x)
    return t + acb(0, pi) - Lam(hinv(x))


def Wf(t, e):
    t = acb(t)
    M = L * t.exp() if e == 0 else -L * t.exp()
    x = xi(M)
    return 1 / fpf(x) if e == 0 else fp(x)


def sqw(t, e):
    """√w_e(t), w_0 = sinθ/(2π(cosh t − cos θ)), w_1 = sinθ/(2π(cosh t + cos θ)); t reell oder komplex."""
    t = acb(t)
    den = t.cosh() - costh if e == 0 else t.cosh() + costh
    return (sinth / (2 * pi) / den).sqrt()


# --- Gauss-Legendre-Knoten, Lebesgue-Konstante, Hermite-Konstante
xs, ws = [], []
for k in range(n):
    x, w = arb.legendre_p_root(n, k, weight=True)
    xs.append(x)
    ws.append(w)


def lag_sum(xb):
    s = arb(0)
    for i in range(n):
        p = arb(1)
        for j in range(n):
            if j != i:
                p *= (xb - xs[j]) / (xs[i] - xs[j])
        s += abs(p)
    return s


nl = 4000
Lam_n = arb(0)
om_max = arb(0)
for k in range(nl):
    xb = arb(-1 + (2 * k + 1) / nl, 1 / nl * (1 + 1e-12))
    Lam_n = arb(max(Lam_n.upper(), lag_sum(xb).upper()))
    om = arb(1)
    for j in range(n):
        om *= (xb - xs[j])
    om_max = arb(max(om_max.upper(), abs(om).upper()))
a_e = (rho + 1 / rho) / 2
b_e = (rho - 1 / rho) / 2
om_min = None
ne = 720
for k in range(ne):
    ph = 2 * pi * (arb(k) + arb("0.5")) / ne
    z = acb(a_e * ph.cos(), b_e * ph.sin())
    z = acb(arb(z.real.mid(), (a_e * pi / ne).upper()), arb(z.imag.mid(), (a_e * pi / ne).upper()))
    om = arb(1)
    for j in range(n):
        om *= abs(z - xs[j])          # Produkt der Beträge (kein Rechteck-Umwickeln)
    v = om.lower()
    om_min = v if om_min is None else arb(min(om_min, v))
dist = a_e - 1
Cint = a_e * om_max / (arb(om_min) * dist)
print(f"n = {n}, ρ = {rho}: Λ_n ≤ {float(Lam_n.upper()):.4f}, max|ω| ≤ {float(om_max.upper()):.3e}, min_E|ω| ≥ {float(om_min):.3e}, "
      f"C_int ≤ {float(Cint.upper()):.3e} ({time.time()-t0:.0f}s)", flush=True)

# --- Paneele und Knoten
npan = int(float((2 * T2 / h).mid()) + 0.5)
pan_a = [-T2 + h * k for k in range(npan)]
tn = []      # Knoten (arb)
om_ = []     # Gewichte
for a in pan_a:
    for i in range(n):
        tn.append(a + h / 2 * (1 + xs[i]))
        om_.append(h / 2 * ws[i])
N = len(tn)
print(f"Paneele: {npan} je Kante, N = {N} je Kante, Matrix {2*N}×{2*N}", flush=True)

# --- Matrix 𝐀 in arb
rowdat = {}
for e in (0, 1):
    rd = []
    for t, w in zip(tn, om_):
        p = Phi(t, e)
        u, y = p.real, p.imag
        fac = w.sqrt() * sqw(t, e) * Wf(t, e) * y.sin() / (2 * pi)
        rd.append((u.cosh(), u.sinh(), y.cos(), fac))
    rowdat[e] = rd
coldat = {}
for e in (0, 1):
    coldat[e] = [(t.cosh(), t.sinh(), w.sqrt() / sqw(t, e).real) for t, w in zip(tn, om_)]
Amid = np.zeros((2 * N, 2 * N), complex)
Arad = np.zeros((2 * N, 2 * N))
for e in (0, 1):
    for i, (cu, su, cy, fac) in enumerate(rowdat[e]):
        for ep in (0, 1):
            sgn = -1 if ep == 0 else 1
            for j, (ct, st, cf) in enumerate(coldat[ep]):
                den = ct * cu - st * su + sgn * cy
                v = fac * (cf / den)
                Amid[e * N + i, ep * N + j] = complex(float(v.real.mid()), float(v.imag.mid()))
                Arad[e * N + i, ep * N + j] = float(v.real.rad()) + float(v.imag.rad())
    print(f"  Matrixzeilen Kante {e} fertig ({time.time()-t0:.0f}s)", flush=True)
Arad += 2.3e-16 * np.abs(Amid)                       # Rundung mid → float
Bm = np.eye(2 * N) + 0.5 * (Amid + Amid.conj().T)
Brad = 0.5 * (Arad + Arad.T) + 2.3e-16 * (np.abs(Bm) + 1)
fehlerB = float(np.sqrt(np.sum(Brad ** 2))) * (1 + 1e-6)
lam_est = np.linalg.eigvalsh(Bm)[0]
s = lam_est - 1e-7
R = np.linalg.cholesky(Bm - s * np.eye(2 * N))   # wirft, falls nicht positiv definit
gam = 4 * (2 * N + 1) * 2.0 ** -53 / (1 - 4 * (2 * N + 1) * 2.0 ** -53)
fr2 = float(np.sum(np.abs(R) ** 2)) * (1 + 1e-6)
lam_lo = s - gam * fr2
print(f"λ_min(B) Schätzung {lam_est:.8f}; verifiziert ≥ {lam_lo:.8f} (Cholesky-Rest {gam*fr2:.2e}); ‖B − B_wahr‖ ≤ {fehlerB:.2e} "
      f"({time.time()-t0:.0f}s)", flush=True)


# --- Analytizität im INNERN der Ellipsen (Voraussetzung der Hermite-Restformel)
# t'-Richtung: Pole von P_{e'}(u+iy, t') liegen bei |Im t'| = y bzw. π − y; 1/√w analytisch für |Im t'| ≤ 0.31 (Re cosh t' ≥ cos 0.31 > cos θ).
#   Nötig also: Ellipsen-Halbachse (h/2)b_e < min(y, π − y) über alle Zeilen — y aus den Einschlüssen unten.
# t-Richtung: Kästchenüberdeckung von {|Im t| ≤ (h/2)b_e, Re t ∈ [−T2 − (h/2)a_e, T2 + (h/2)a_e]}: Log-/Wurzelargumente meiden (−∞, 0],
#   f'∘f ≠ 0, 0 < Im Φ_0 < 2π (Pole von φ), 0 < Im Φ_1 < π (Pole von φ und φ1), entsprechend für Φ* = conj Φ(conj t).


def ok_achse(z):
    return bool(z.real > 0) or bool(z.imag > 0) or bool(z.imag < 0)


def analytisch_t():
    hb = (h / 2 * b_e).upper() * (1 + 1e-9)
    ha = (h / 2 * a_e).upper()
    dx = arb("0.02")
    nx = int(float(((2 * (T2 + ha)) / dx).upper())) + 1
    ny = int(float((2 * arb(hb) / dx).upper())) + 1
    schlecht = 0
    ymin = {0: None, 1: None}
    ymax = {0: None, 1: None}
    for i in range(nx):
        xr = -T2 - ha + dx * (arb(i) + arb("0.5"))
        for j in range(ny):
            yi = -arb(hb) + dx * (arb(j) + arb("0.5"))
            z = acb(arb(xr.mid(), (dx / 2).upper() * 1.000001), arb(yi.mid(), (dx / 2).upper() * 1.000001))
            for zz, konj in ((z, False), (z.conjugate(), True)):
                Mz = L * zz.exp()
                x0 = xi(Mz)
                r0 = 1 + eps * (L - Lb) / (c0 + eps * (x0 + Lb))
                g = acb(zz) + r0.log()
                c1 = ok_achse(r0) and bool(abs(fpf(x0)) > 0) and bool(g.imag > 0) and bool(g.imag < 2 * pi)
                M1 = -L * zz.exp()
                x1 = xi(M1)
                rad = c0 * c0 - 4 * eps * eps + 4 * eps * x1
                c2 = ok_achse(rad)
                if c2:
                    v = (-c0 + rad.sqrt()) / (2 * eps)
                    r1 = 1 + eps * (L - Lb) / (c0 + eps * (v + Lb))
                    Hh = zz + acb(0, pi) - r1.log()
                    c2 = ok_achse(r1) and bool(Hh.imag > 0) and bool(Hh.imag < pi)
                if not (c1 and c2):
                    schlecht += 1
    return schlecht


nschlecht = analytisch_t()
print(f"Analytizität im Innern der t-Ellipsen: {'geprüft' if nschlecht == 0 else 'NICHT geprüft'} ({nschlecht} kritische Kästchen) "
      f"({time.time()-t0:.0f}s)", flush=True)


# --- Interpolationsfehler (HS)
def ellipse_boxes(mid):
    """Kästchen auf mid + (h/2)·∂E_ρ."""
    out = []
    for k in range(nell):
        ph = 2 * pi * (arb(k) + arb("0.5")) / nell
        zr = mid + h / 2 * a_e * ph.cos()
        zi = h / 2 * b_e * ph.sin()
        r = (h / 2 * a_e * pi / nell).upper() * (1 + 1e-12)
        out.append(acb(arb(zr.mid(), r), arb(zi.mid(), r)))
    return out


def sub_balls(a):
    out = []
    for k in range(nsub):
        m = a + h * (arb(k) + arb("0.5")) / nsub
        out.append(arb(m.mid(), (h / (2 * nsub)).upper() * (1 + 1e-12)))
    return out


def phi0(z):
    return 1 / (1 - z.exp())


def phi1(z):
    return 1 / (1 + z.exp())


# vorberechnet je Paneel: reelle Teilkugeln (Zeilendaten), Ellipsenkästchen (für t komplex) und (für t' komplex)
pan_mid = [a + h / 2 for a in pan_a]
row_sub = {e: [] for e in (0, 1)}      # je Paneel: Liste (u, y, sqw·W) für reelle Teilkugeln
row_ell = {e: [] for e in (0, 1)}      # je Paneel: Liste (Φ, Φ*, sqw·W) für komplexe t auf der Ellipse
col_sub = {e: [] for e in (0, 1)}      # je Paneel: Liste (t' reell, 1/sqw)
col_ell = {e: [] for e in (0, 1)}      # je Paneel: Liste (t' komplex, 1/sqw)
for k, a in enumerate(pan_a):
    sb = sub_balls(a)
    eb = ellipse_boxes(pan_mid[k])
    for e in (0, 1):
        rs = []
        for t in sb:
            p = Phi(t, e)
            rs.append((p.real, p.imag, sqw(t, e) * Wf(t, e)))
        row_sub[e].append(rs)
        re_ = []
        for z in eb:
            p = Phi(z, e)
            ps = Phi(z.conjugate(), e).conjugate()
            re_.append((p, ps, sqw(z, e) * Wf(z, e)))
        row_ell[e].append(re_)
        col_sub[e].append([(t, 1 / sqw(t, e)) for t in sb])
        col_ell[e].append([(z, 1 / sqw(z, e)) for z in eb])


def hull(lst):
    out = lst[0]
    for v in lst[1:]:
        out = tuple(a.union(b) for a, b in zip(out, v))
    return out


row_sub_h = {e: [hull(v) for v in row_sub[e]] for e in (0, 1)}
row_ell_h = {e: [hull(v) for v in row_ell[e]] for e in (0, 1)}
col_sub_h = {e: [hull(v) for v in col_sub[e]] for e in (0, 1)}
col_ell_h = {e: [hull(v) for v in col_ell[e]] for e in (0, 1)}
ymin_all = min(float(rs[1].lower()) for e in (0, 1) for P_ in row_sub[e] for rs in P_)
pimy_all = min(float((pi - rs[1]).lower()) for e in (0, 1) for P_ in row_sub[e] for rs in P_)
halb = float((h / 2 * b_e).upper())
print(f"t'-Richtung: Halbachse {halb:.4f} < min y = {ymin_all:.4f} und < min(π − y) = {pimy_all:.4f}: "
      f"{'geprüft' if halb < min(ymin_all, pimy_all) else 'NICHT geprüft'}", flush=True)
print(f"Vorberechnung Paneeldaten fertig ({time.time()-t0:.0f}s)", flush=True)


def M_tp(e, ep, P, Q, rows=None, cols=None):
    """sup |k(t,t')|, t ∈ P reell (Teilkugeln), t' ∈ ∂E(Q)."""
    m = arb(0)
    sg = -1 if ep == 0 else 1
    rows = row_sub[e][P] if rows is None else rows
    cols = col_ell[ep][Q] if cols is None else cols
    for (u, y, fw) in rows:
        sy, cy = y.sin(), y.cos()
        for (z, iw) in cols:
            v = fw * iw * sy / (2 * pi) / ((z - u).cosh() + sg * cy)
            m = arb(max(m.upper(), abs(v).upper()))
    return m


def M_t(e, ep, P, Q, rows=None, cols=None):
    """sup |k^c(t,t')|, t ∈ ∂E(P) komplex, t' ∈ Q reell."""
    m = arb(0)
    rows = row_ell[e][P] if rows is None else rows
    cols = col_sub[ep][Q] if cols is None else cols
    for (p, ps, fw) in rows:
        for (tq, iw) in cols:
            if ep == 0:
                pc = (phi0(p - tq) - phi0(ps - tq)) / acb(0, 2 * pi)
            else:
                pc = -(phi1(p - tq) - phi1(ps - tq)) / acb(0, 2 * pi)
            v = fw * iw * pc
            m = arb(max(m.upper(), abs(v).upper()))
    return m


def M_grob(e, ep, P, Q):
    """Fernpaare: ein Kästchen um ganz E(P) bzw. E(Q) — grobe obere Schranke für beide Richtungen."""
    # t' komplex über das ganze Kästchen von E(Q), t reell über P
    zq = acb(arb(pan_mid[Q].mid(), (h / 2 * a_e).upper()), arb(0, (h / 2 * b_e).upper()))
    tp_ = arb(pan_mid[P].mid(), (h / 2).upper())
    p = Phi(tp_, e)
    u, y = p.real, p.imag
    sg = -1 if ep == 0 else 1
    v1 = sqw(tp_, e) * Wf(tp_, e) / sqw(zq, ep) * y.sin() / (2 * pi) / ((zq - u).cosh() + sg * y.cos())
    zp = acb(arb(pan_mid[P].mid(), (h / 2 * a_e).upper()), arb(0, (h / 2 * b_e).upper()))
    tq = arb(pan_mid[Q].mid(), (h / 2).upper())
    pz = Phi(zp, e)
    pzs = Phi(zp.conjugate(), e).conjugate()
    if ep == 0:
        pc = (phi0(pz - tq) - phi0(pzs - tq)) / acb(0, 2 * pi)
    else:
        pc = -(phi1(pz - tq) - phi1(pzs - tq)) / acb(0, 2 * pi)
    v2 = sqw(zp, e) * Wf(zp, e) / sqw(tq, ep) * pc
    return abs(v1), abs(v2)


def phi_schranke(re_lo, re_hi):
    """Für z mit Re z ∈ [re_lo, re_hi]: Schranke für |φ(z)| (Re > 0) bzw. |φ(z) − 1| (Re < 0), gleich für φ1 = 1/(1+e^z)."""
    if bool(re_lo > 0):
        return 1 / (re_lo.exp() - 1)
    if bool(re_hi < 0):
        return re_hi.exp() / (1 - re_hi.exp())
    return None


def M_fern(e, ep, P, Q):
    """Fernpaare mit reeller Intervallrechnung (keine weiten komplexen Kästchen):
    M_t' über |cosh(α+iβ) ∓ cos y| ≥ sinh|α| − 1 und |1/√w(z)| ≤ √(2π(cosh Re z + 1)/sin θ);
    M_t über |P^c| ≤ (1/2π)(φ-Schranke(Φ−t') + φ-Schranke(Φ*−t')) (die Einsen heben sich bei Re < 0 weg)."""
    inf = arb("1e300")
    # M_t': t reell in P (Hülle), t' auf ∂E(Q)
    u, y, fw = row_sub_h[e][P]
    zr = arb(pan_mid[Q].mid(), (h / 2 * a_e).upper())
    al = abs(zr - u)
    almin = arb(al.lower())
    if bool(almin.sinh() > 1):
        iw = (2 * pi * (zr.cosh() + 1) / sinth).sqrt()
        m1 = arb((abs(fw) * iw / (2 * pi) / (almin.sinh() - 1)).upper())
    else:
        m1 = M_tp(e, ep, P, Q)
    # M_t: t komplex auf ∂E(P) (Hülle der Feinwerte), t' reell in Q (Hülle)
    p_, ps, fw2 = row_ell_h[e][P]
    tq, iw2 = col_sub_h[ep][Q]
    b1 = phi_schranke(arb((p_.real - tq).lower()), arb((p_.real - tq).upper()))
    b2 = phi_schranke(arb((ps.real - tq).lower()), arb((ps.real - tq).upper()))
    same = (bool(p_.real - tq > 0) and bool(ps.real - tq > 0)) or (bool(p_.real - tq < 0) and bool(ps.real - tq < 0))
    if b1 is not None and b2 is not None and same:
        m2 = arb((abs(fw2) * abs(iw2) * (b1 + b2) / (2 * pi)).upper())
    else:
        m2 = M_t(e, ep, P, Q)
    return m1, m2


hs2 = arb(0)
maxerr = arb(0)
for e in (0, 1):
    for ep in (0, 1):
        for P in range(npan):
            for Q in range(npan):
                if abs(P - Q) >= fern:
                    m1, m2 = M_fern(e, ep, P, Q)
                else:
                    m1, m2 = M_tp(e, ep, P, Q), M_t(e, ep, P, Q)
                err = Cint * (m1 + Lam_n * m2)
                if not err.is_finite():
                    print("  NICHT ENDLICH:", e, ep, P, Q, m1, m2)
                hs2 += h * h * err * err
                maxerr = arb(max(maxerr.upper(), err.upper()))
        print(f"  Block ({e},{ep}) fertig, HS² bisher ≤ {float(hs2.upper()):.3e} ({time.time()-t0:.0f}s)", flush=True)
hs = hs2.sqrt()
print(f"‖k − k_N‖_HS ≤ {float(hs.upper()):.3e} (max Paneelfehler {float(maxerr.upper()):.2e})")
cm = min(1.0, lam_lo - fehlerB) - float(hs.upper())
print(f"ERGEBNIS: c_m ≥ {cm:.6f}  (T2 = {T2}, h = {h}, n = {n}, ρ = {rho}) ({time.time()-t0:.0f}s)")
