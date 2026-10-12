"""k36 — Asymptotik des adjungierten Kopfes längs der Linse γ nahe L für grosse n, aus dem Eckprofil Y (Kegelskala; gültig für 1 ≪ A ≪ n^{Re a}
über die Anpassung an die äussere Lösung — Annahme!).  Frage: bleibt max_γ |Λ|/|Λ(L)| (und die lokale L²-Masse) beschränkt, wenn n → ∞?
Hintergrund (BERICHT §20): das Eckmuster dreht sich gegenüber festen Richtungen um Im(a)·ln n (Spirale, komplexes κ₁); die Γin-Tangente wandert in
den Wachstumssektor von Y. Linsen: tangential (S, τ₁, p) und Keil (Höhe c·τ).
Aufruf: python k36_linse_asymptotik.py
"""
import math
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 140
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
RHO = abs(L)
ALPHA = L / (4 * mp.pi)
A_EXP = 1 + 1j * ALPHA
LAM_A = mp.mpc(-0.203297, -1.406008)
BA = 0.5 - 2j * ALPHA
K = 900
CY = [mp.gamma(BA - 1j * ALPHA * k) / (mp.factorial(k) * mp.gamma(BA)) for k in range(K)]
ARGL = mp.arg(L)


def Y(u):
    s = mp.mpc(0); p = mp.mpc(1)
    for c in CY:
        s += c * p; p *= -u
    return s


def s_tang(tau, S=0.0166, tau1=math.radians(12), p=1.3):
    return S * (1 - mp.exp(-(tau / tau1) ** p)) if tau > 0 else 0


def s_keil(tau, c):
    return c * tau if tau > 0 else 0


def profil(n, sfun, Amax, nA=60):
    dn = mp.power(n, -A_EXP)
    scale = RHO * abs(dn)
    wL = RHO / L
    out = []
    for i in range(1, nA + 1):
        A = Amax * (i / nA) ** 1.5
        tau = A * scale / RHO
        x = RHO * mp.exp(sfun(tau)) * mp.expjpi((ARGL + tau) / mp.pi)
        u = LAM_A * (RHO / x - wL) / dn
        out.append((float(A), float(abs(Y(u)))))
    return out


def main():
    lenses = [("tangential p=1.3, τ₁=12°", lambda t: s_tang(t)), ("tangential p=1.15, τ₁=5°", lambda t: s_tang(t, tau1=math.radians(5), p=1.15)),
              ("Keil c=0.25", lambda t: s_keil(t, 0.25)), ("Keil c=0.5", lambda t: s_keil(t, 0.5)), ("ohne Linse (Kreis)", lambda t: 0)]
    for n in [300, 600, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8]:
        print("n = %g: Kegelskala ρ|n^{−a}| = %.3e" % (n, float(RHO * abs(mp.power(n, -A_EXP)))), flush=True)
        for name, sf in lenses:
            pr = profil(n, sf, 90.0)
            m = max(pr, key=lambda t: t[1])
            l2 = math.sqrt(sum(((pr[i][0] - pr[i - 1][0]) * pr[i][1] ** 2) for i in range(1, len(pr))))
            tail = pr[-1][1]
            print("   %-26s max|Y| = %.3g bei A = %.1f; |Y(A=90)| = %.3g; (∫_0^90|Y|²dA)^{1/2} = %.3g" % (name, m[1], m[0], tail, l2), flush=True)


if __name__ == "__main__":
    main()
