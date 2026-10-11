"""b43 (Runde 4): Zweig-Schritt N → N+1 als Kopfantwort.
Exakt: F_{N+1}(pad p_N) = c_{N+1} x^{N+1},  c_{N+1} := [x^{N+1}](p_N∘p_N) − 1/(N+1)!  (alle tieferen Gleichungen gelten schon).
Homotopie T_{N+1}[p∘p − eˣ] = σ x^{N+1}, σ: c → 0, ergibt p_{N+1} − pad p_N = −c ∫₀¹ J_{N+1}(p(s))⁻¹ x^{N+1} ds ≈ −c·t_{N+1}.
Ebenso exakt: x₁(N+1) − pad x₁(N) = −ρ^f_N · t_{N+1},  ρ^f_N := [x^{N+1}] J_∞(x₁(N) − f_K),  t_M := J_M(f_K)⁻¹ x^M.
Gemessen: |c||L|^{N+1}, T_{N+1} = sup_D̄|t_{N+1}|/|L|^{N+1}, Abstand der Vorhersage −c t_{N+1} (mit J(f) bzw. J(pad p_N))
zum echten Schritt p_{N+1} − pad p_N; Prüfung der x₁-Identität bei N = 296…299 (Daten von Strang A)."""
import json, sys, time
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 1000
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
P = json.load(open(G + "polynome_195_bis_300.json"))
lin = json.load(open(G + "lin_PiN_fK_bis_1000.json"))
mp.mp.dps = 60
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
aL = abs(Lm)
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 240) for j in range(240)] + [Lm]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def supL(c):
    cm = [mp.mpf(x.mid().str(50, radius=False)) if hasattr(x, "mid") else mp.mpf(x) for x in c]
    return max(abs(horner(cm, z)) for z in kreis)


def coef(p, k):
    c = p.coeffs()
    return c[k] if k < len(c) else arb(0)


def jac(a, M):
    """J_M(a)[k,j] = [x^k] a^j + [x^{k−j}] a'∘a, k,j ≤ M; a: Liste arb (Länge ≤ M+1)."""
    f = arb_poly(a)
    pw = [arb_poly([1])]
    for _ in range(M):
        pw.append((pw[-1] * f).truncate(M + 1))
    fp = f.derivative()
    fpf = arb_poly([0])
    for m in range(len(fp.coeffs())):
        fpf += coef(fp, m) * pw[m] if m <= M else 0
    Jm = arb_mat(M + 1, M + 1)
    for k in range(M + 1):
        for j in range(M + 1):
            Jm[k, j] = coef(pw[j], k) + (coef(fpf, k - j) if j <= k else arb(0))
    return Jm, pw


def fakul(n):
    r = arb(1)
    for k in range(2, n + 1):
        r *= k
    return r


fK = [arb(c) for c in K[:400]]
ergebnis = {}
for N in [int(s) for s in sys.argv[1:]] or [200, 250, 299]:
    t0 = time.time()
    M = N + 1
    JfM, pwf = jac(fK[:M + 1], M)
    Jfi = JfM.inv()
    t = [Jfi[k, M] for k in range(M + 1)]
    TM = supL(t) / aL ** M
    # c_{N+1} aus p_N
    pN = [arb(s) for s in P[str(N)]]
    pN1 = [arb(s) for s in P[str(N + 1)]]
    fp_ = arb_poly(pN)
    s, pw = arb(0), arb_poly([1])
    for j in range(N + 1):
        s += pN[j] * coef(pw, M)
        pw = (pw * fp_).truncate(M + 1)
    c = s - 1 / fakul(M)
    c_norm = abs(float((c * arb(mp.nstr(aL, 50)) ** M).mid()))
    vorh_f = [-c * x for x in t]
    pad = pN + [arb(0)]
    Jp, _ = jac(pad, M)
    d1 = Jp.solve(arb_mat([[arb(0)]] * M + [[c]]))
    vorh_p = [-d1[k, 0] for k in range(M + 1)]
    echt = [pN1[k] - pad[k] for k in range(M + 1)]
    s_echt = supL(echt); s_f = supL(vorh_f); s_p = supL(vorh_p)
    abw_f = supL([echt[k] - vorh_f[k] for k in range(M + 1)]); abw_p = supL([echt[k] - vorh_p[k] for k in range(M + 1)])
    print(f"N={N}→{M} ({time.time()-t0:.0f}s): |c_(N+1)||L|^(N+1) = {c_norm:.3e}, T_(N+1) = {mp.nstr(TM,4)};"
          f" sup|p_(N+1) − pad p_N| = {mp.nstr(s_echt,3)}; Vorhersage −c·t (J(f)): {mp.nstr(s_f,3)} (Abw. {mp.nstr(abw_f,2)}),"
          f" −c·J(pad p_N)⁻¹x^(N+1): {mp.nstr(s_p,3)} (Abw. {mp.nstr(abw_p,2)})", flush=True)
    # x₁-Identität, falls Daten für N und N+1 vorliegen
    if str(N) in lin and str(M) in lin:
        x1N = [arb(s) for s in lin[str(N)]["b"]]; x1M = [arb(s) for s in lin[str(M)]["b"]]
        # ρ^f_N = [x^{N+1}] (x₁∘f + g x₁ − eˣ − g f)
        Mm = M + 1
        f = arb_poly(fK[:Mm + 1])
        pwm = [arb_poly([1])]
        for _ in range(N):
            pwm.append((pwm[-1] * f).truncate(Mm))
        x1f = sum((x1N[j] * coef(pwm[j], M) for j in range(N + 1)), arb(0))
        ex = [arb(1)]
        for k in range(1, Mm + 1):
            ex.append(ex[-1] / k)
        fpc = [coef(f.derivative(), k) for k in range(Mm)]
        gk = []
        for n in range(Mm):
            gk.append((ex[n] - sum((fpc[k] * gk[n - k] for k in range(1, n + 1)), arb(0))) / fpc[0])
        u = [x1N[k] - fK[k] if k <= N else -fK[k] for k in range(Mm)]
        gu = sum((gk[M - j] * u[j] for j in range(M + 1)), arb(0))
        rho = x1f + gu - 1 / fakul(M)            # ρ^f = [x^{N+1}](x₁∘f) + [x^{N+1}](g·(x₁ − f)) − [x^{N+1}](f∘f)
        pred = [-rho * x for x in t]
        echt1 = [x1M[k] - (x1N[k] if k <= N else arb(0)) for k in range(M + 1)]
        print(f"      x₁-Identität: sup|x₁(N+1) − pad x₁(N)| = {mp.nstr(supL(echt1),3)}, −ρ^f t: {mp.nstr(supL(pred),3)},"
              f" Abweichung {mp.nstr(supL([echt1[k]-pred[k] for k in range(M+1)]),2)}; |ρ^f||L|^(N+1) = {abs(float((rho*arb(mp.nstr(aL,50))**M).mid())):.3e}", flush=True)
