# Entfaltung unter eta: b = exp(lam e^{-lam}) < eta, Fixpunkte a = e^lam (anziehend, Multiplikator lam = e^{-s})
# und r = e^mu (abstossend, mu > 1). Uebergang theta(w) = alpha_r(Psi_a(w)) - w zwischen den Koenigs-Abelfunktionen
# auf dem reellen Stueck (a, r). Fourier-Koeffizienten t_n vergleichen mit Hornabbildung a_n:
#   Vorhersage |t_n| ~ |a_n| e^{-2 pi^2 n/s} * (Potenzfaktor)^n, also t_2/t_1^2 -> a_2/a_1^2 (normierungsfrei).
# Aufruf: e6_entfaltung.py dps M s1 s2 ...
import mpmath as mp, sys, json

def setup(s):
    lam = mp.exp(-s)
    lnb = lam * mp.exp(-lam)
    mu = -mp.lambertw(-lnb, -1).real
    return lam, mu, lnb, mp.exp(lam), mp.exp(mu)

def alpha_a(z, P):
    lam, mu, lnb, a, r = P
    n = 0; tol = mp.mpf(10) ** (-(mp.mp.dps // 3))
    v = z - a
    while abs(v) > tol:
        z = mp.exp(lnb * z); n += 1; v = z - a
    kap = lnb / (2 * (1 - lam))
    sig = (v + kap * v**2) * mp.power(lam, -n)
    return mp.log(sig) / mp.log(lam)

def alpha_r(z, P):
    lam, mu, lnb, a, r = P
    n = 0; tol = mp.mpf(10) ** (-(mp.mp.dps // 3))
    v = z - r
    while abs(v) > tol:
        z = mp.log(z) / lnb; n += 1; v = z - r
    kap = lnb / (2 * (1 - mu))
    sig = (v + kap * v**2) * mp.power(mu, n)
    return mp.log(-sig) / mp.log(mu)

def psi_a(w, P):
    lam, mu, lnb, a, r = P
    zeta = mp.power(lam, w)
    m = 0; tol = mp.mpf(10) ** (-(mp.mp.dps // 3))
    while abs(zeta) > tol:
        zeta *= lam; m += 1
    kap = lnb / (2 * (1 - lam))
    v = zeta - kap * zeta**2 + 2 * kap**2 * zeta**3
    z = a + v
    for _ in range(m):
        z = mp.log(z) / lnb
    return z

if __name__ == '__main__':
    mp.mp.dps = int(sys.argv[1]); M = int(sys.argv[2])
    a12 = mp.mpc('-0.0253204930997575512295974960086', '1.57421906523195162665651202605')
    a13 = mp.mpc('-2.25479912452380311104015274723', '-0.0244087142968831711717887164663')
    absa1 = mp.mpf('64.15318451854875769074032856744156932425')
    out = {}
    for sarg in sys.argv[3:]:
        s = mp.mpf(sarg); P = setup(s)
        lam, mu, lnb, a, r = P
        z0 = (a + r) / 2
        w0 = alpha_a(z0, P)
        th = []
        for j in range(M):
            w = w0 + mp.mpf(j) / M
            z = psi_a(w, P)
            th.append(alpha_r(z, P) - w)
        t = {}
        for n in range(-3, 4):
            t[n] = mp.fsum(th[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M * mp.expjpi(-2 * n * w0)
        b = mp.exp(lnb)
        sr = mp.log(mu)
        print(f"s={sarg} b={mp.nstr(b,12)} eta-b={mp.nstr(mp.e**(1/mp.e)-b,5)} ln mu={mp.nstr(sr,10)}")
        print("   t0 =", mp.nstr(t[0], 20))
        for n in [1, 2, 3]:
            print(f"   t{n} =", mp.nstr(t[n], 15), " |t-n|/|tn|=", mp.nstr(abs(t[-n]) / abs(t[n]), 8))
        x1 = abs(t[1]) * mp.exp(2 * mp.pi**2 / s) / absa1
        print("   |t1| e^{2pi^2/s}/|a1| =", mp.nstr(x1, 12), "  mit s_r:", mp.nstr(abs(t[1]) * mp.exp(2 * mp.pi**2 / sr) / absa1, 12),
              "  mit (s+s_r)/2:", mp.nstr(abs(t[1]) * mp.exp(4 * mp.pi**2 / (s + sr)) / absa1, 12))
        print("   t2/t1^2 =", mp.nstr(t[2] / t[1]**2, 15), "  (a2/a1^2 =", mp.nstr(a12, 8), ")")
        print("   t3/t1^3 =", mp.nstr(t[3] / t[1]**3, 15), "  (a3/a1^3 =", mp.nstr(a13, 8), ")")
        out[sarg] = {str(n): [str(mp.re(v)), str(mp.im(v))] for n, v in t.items()}
        sys.stdout.flush()
    json.dump(out, open(f"cache/e6_{mp.mp.dps}_{M}_{'_'.join(sys.argv[3:])}.json", "w"), indent=1)
