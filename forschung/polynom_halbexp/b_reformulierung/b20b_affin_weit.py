"""b20: Affines Modell der Grenzschicht: f(z) = L + μ(z − L), μ = e^{L/2} (f∘f = L + L(z−L) = Linearisierung von exp bei L).
J_∞u = u∘f + μu; für u_s = (1 − z/L)^s gilt J_∞u_s = (μ^s + μ)u_s exakt, also Π_N u_s = (μ^s + μ) J_N⁻¹ T_N u_s.
Kernfunktionen: (z − L)^{1+(2n+1)γ}. Modus-Inhalt von e = u − Π_N u über θ(w) = e(z)/(z − L), z = L + r0·e^{L(w − w0)}.
Frage: Amplitude von (z−L)^{1−γ} ~ N^{1−γ−s} (naiver Koeffizientenquotient) oder N^{2−γ−s} (wie im echten Problem)?"""
import sys, math
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 110
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
mu = mp.exp(L / 2)
gam = 2j * mp.pi / L
a0 = L * (1 - mu)
print("μ =", mp.nstr(mu, 8), " f(0) =", mp.nstr(a0, 8), " |f(0)| =", mp.nstr(abs(a0), 6))

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

def Jmat(N):
    # Spalte j: Koeffizienten von (a0 + μx)^j + μ x^j
    J = mp.matrix(N + 1, N + 1)
    for j in range(N + 1):
        for k in range(j + 1):
            J[k, j] = mp.binomial(j, k) * a0 ** (j - k) * mu ** k
        J[j, j] += mu
    return J

M, r0 = 64, mp.mpf("0.55")
def moden(e_fun):
    out = {}
    pts = []
    for j in range(M):
        w = mp.mpf(j) / M
        z = L + r0 * mp.exp(L * w) * mp.exp(1j * mp.pi)     # Spirale um L, nach innen in die Scheibe
        pts.append((w, z))
    th = [e_fun(z) / (z - L) for w, z in pts]
    for n in (-1, 1):
        out[n] = mp.fsum(t * mp.exp(-2j * mp.pi * n * w) for t, (w, z) in zip(th, pts)) / M
    return out, max(abs(z) for w, z in pts)

for s in (1 + gam, mp.mpf("5.447"), mp.mpf("2.5")):
    res = {}
    for N in (40, 60, 80, 100, 120):
        J = Jmat(N)
        uk = [mp.binomial(s, k) * (-1 / L) ** k for k in range(N + 1)]
        b = mp.matrix([(mu ** s + mu) * c for c in uk])
        P = [mp.mpc(0)] * (N + 1)          # J ist obere Dreiecksmatrix: Rückwärtseinsetzen
        for k in range(N, -1, -1):
            s_ = b[k] - mp.fsum(J[k, j] * P[j] for j in range(k + 1, N + 1))
            P[k] = s_ / J[k, k]
        efun = lambda z: (1 - z / L) ** s - horner(P, z)
        m, zmax = moden(efun)
        res[N] = m[-1]
    Ns = sorted(res)
    ex = [complex(mp.log(res[b2] / res[a2]) / mp.log(mp.mpf(b2) / a2)) for a2, b2 in zip(Ns, Ns[1:])]
    print(f"s = {mp.nstr(s,4)} (max|z| auf Spirale {mp.nstr(zmax,3)}): c_-1 = {[mp.nstr(res[N],3) for N in Ns]}")
    print(f"   lokale Exponenten {[f'{e.real:+.2f}{e.imag:+.2f}i' for e in ex]};  1−γ−s = {complex(1 - gam - s):.3f},  2−γ−s = {complex(2 - gam - s):.3f}")
