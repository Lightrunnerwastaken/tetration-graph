"""m1: Monodromie des regulaeren Erzeugers J_L = L*sigma/sigma' von exp (Basis e) ueber C \\ P, P = {0,1,e,e^e,...}.

Methode (unabhaengig von der Iterierte-Log-Methode aus Runde 18/10):
  Psi = sigma^{-1} ist ganz, Psi(L u) = exp(Psi(u)), Psi(0) = L, Psi'(0) = 1.
  Ein Blatt von sigma an der Stelle z ist ein Punkt u der Faser Psi^{-1}(z).
  Fortsetzung von sigma laengs eines Weges z(t) = Hochheben: Newton-Verfolgung von Psi(u(t)) = z(t).
  J_L auf dem Blatt u: J = L * u * Psi'(u).

Kodierung der Faser ueber L (Satz M2 im Bericht):
  code(u) = (m_1, code(u')), Psi(u/L) = L + 2 pi i m_1, u' = Hebung von ell_{m_1}^{-1} ab u/L.
  ell_m: L -> -1 + i Im L -> -1 + i(Im L + 2 pi m) -> L + 2 pi i m   (laeuft links an 0 vorbei).
Vorhersage: a = Schleife um 0:  (m_1, m_2, ...) -> (m_1 + 1, m_2, ...)
            b_k = Schleife um p_k: m_{k+1} += 1, falls m_1 = ... = m_k = 0, sonst nichts.

Aufruf: python m1_monodromie.py [A|B]   A: Gleitkomma (complex), B: mpmath dps 30, feinere Schritte.
"""
import sys
import cmath
import math

EINST = sys.argv[1] if len(sys.argv) > 1 else "A"

TEIL2 = EINST == "B2"            # B2 = Einstellung B, zweite Haelfte der Teilmenge
if TEIL2:
    EINST = "B"

if EINST == "A":
    EXP, LOG, PI = cmath.exp, cmath.log, math.pi
    C = complex
    def absv(x): return abs(x)
    NK, R0, H = 40, 0.08, 0.02      # Reihenglieder, Radius fuer Psi_loc, Grundschritt
    TOL = 1e-13
    def L_berechnen():
        z = 0.3 + 1.3j
        for _ in range(60):
            z = z - (cmath.exp(z) - z) / (cmath.exp(z) - 1)
        return z
else:
    import mpmath as mp
    mp.mp.dps = 30
    EXP, LOG, PI = mp.exp, mp.log, mp.pi
    def C(x, y=0):
        return mp.mpc(x, y)
    def absv(x): return abs(x)
    NK, R0, H = 50, 0.05, 0.012
    TOL = mp.mpf(10) ** -25
    def L_berechnen():
        return mp.findroot(lambda z: mp.exp(z) - z, mp.mpc("0.318", "1.337"))

L = L_berechnen()
I = C(0, 1)
P = [C(0)]
for _ in range(4):
    P.append(EXP(P[-1]))           # 0, 1, e, e^e, e^e^e
P = P[:4]                          # nur 0,1,e,e^e als Hindernisse benutzen

# Koeffizienten von Psi(u) = L + sum c_n u^n
c = [C(0)] * (NK + 1)
E = [C(0)] * (NK + 1)
c[1] = C(1); E[0] = C(1); E[1] = C(1)
for n in range(2, NK + 1):
    T = sum(k * c[k] * E[n - k] for k in range(1, n)) / n
    c[n] = L * T / (L ** n - L)
    E[n] = c[n] + T


def psi(u):
    """Psi(u), Psi'(u)."""
    n = 0
    v = u
    while absv(v) > R0:
        v = v / L
        n += 1
    w = L; dw = C(0)
    vp = C(1)
    for k in range(1, NK + 1):
        dw += k * c[k] * vp
        vp = vp * v
        w += c[k] * vp
    dw = dw / L ** n
    for _ in range(n):
        w = EXP(w)
        dw = w * dw
    return w, dw


def mindist(z):
    return min(absv(z - p) for p in P)


def verfolge(u, punkte):
    """Hebt den Polygonzug punkte (z-Werte) ab dem Fasernpunkt u (Psi(u)=punkte[0]) hoch."""
    w, dw = psi(u)
    assert absv(w - punkte[0]) < 1e-8, (w, punkte[0])
    for z0, z1 in zip(punkte[:-1], punkte[1:]):
        laenge = absv(z1 - z0)
        t = 0.0
        while t < 1.0:
            z = z0 + (z1 - z0) * t
            h = min(H, 0.03 * mindist(z)) / max(laenge, 1e-30)
            t1 = min(1.0, t + max(h, 1e-6))
            zn = z0 + (z1 - z0) * t1
            u = u + (zn - z) / dw          # Praediktor
            for _ in range(10):
                w, dw = psi(u)
                d = (w - zn) / dw
                u = u - d
                if absv(d) < TOL * (1 + absv(u)):
                    break
            else:
                raise RuntimeError("Newton konvergiert nicht")
            t = t1
        w, dw = psi(u)
    return u


def kreis(m, r, start_winkel, n=64):
    """Polygon fuer Kreis um m, Radius r, gegen den Uhrzeigersinn, ab Winkel start_winkel."""
    return [m + r * EXP(I * (start_winkel + 2 * PI * j / n)) for j in range(n + 1)]


def ell(m):
    """Verbindungsweg ell_m von L nach L + 2 pi i m (links an 0 vorbei)."""
    if m == 0:
        return [L, L]
    return [L, C(-1) + I * L.imag, C(-1) + I * (L.imag + 2 * PI * m), L + 2 * PI * I * m]


def lasso(k, r):
    """Schleife b_k: von L geradlinig Richtung p_k, Kreis Radius r gegen Uhrzeigersinn, zurueck."""
    p = P[k]
    richt = (L - p) / absv(L - p)
    q = p + r * richt
    ang = cmath.phase(complex(richt)) if EINST == "A" else mp.arg(richt)
    return [L, q] + kreis(p, r, ang)[1:] + [L]


# Schleife a um 0: Rechteck [-1,0.5]x[-1, Im L] durch L, gegen Uhrzeigersinn (homotop zum Lasso laengs [L,0])
A_WEG = [L, C(-1) + I * L.imag, C(-1, -1), C(0.5, -1), C(0.5) + I * L.imag, L]
WEGE = {"a": A_WEG, "b1": lasso(1, 0.3), "b2": lasso(2, 0.5), "b3": lasso(3, 1.0)}


def umkehr(weg):
    return list(reversed(weg))


def kodiere(u, maxlen=12):
    code = []
    for _ in range(maxlen):
        if absv(u) < 1e-9:
            return tuple(code)
        v = u / L
        w, _ = psi(v)
        m = int(round(float(((w - L) / (2 * PI)).imag)))
        assert absv(w - L - 2 * PI * I * m) < 1e-6, w
        u = verfolge(v, umkehr(ell(m)))
        code.append(m)
    return tuple(code) + ("...",)


def dekodiere(code):
    u = C(0)
    for m in reversed(code):
        v = verfolge(u, ell(m))       # Fasernpunkt ueber L + 2 pi i m
        u = L * v
    return u


def vorhersage(name, code):
    code = list(code) + [0] * 6
    if name == "a":
        code[0] += 1
    else:
        k = int(name[1:])
        if all(x == 0 for x in code[:k]):
            code[k] += 1
    while code and code[-1] == 0:
        code.pop()
    return tuple(code)


def norm(code):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    return tuple(code)


if __name__ == "__main__":
    print(f"Einstellung {EINST}: L = {L}")
    print("Koeffizienten |c_n| (n=10,20,30,40):", [float(absv(c[n])) for n in (10, 20, 30, 40)])
    # Gegenprobe mit Runde 18: J_L(L+2 pi i) laengs der Senkrechten = -14.3223459699 + 10.7051679582 i
    s1 = verfolge(C(0), [L, L + 2 * PI * I])
    w, dw = psi(s1)
    print("s1 = sigma(L+2pi i) laengs Senkrechter =", s1, "  J_L = L s1 Psi'(s1) =", L * s1 * dw)
    s1b = verfolge(C(0), ell(1))
    print("   ueber ell_1 (links herum):", s1b, " Differenz", float(absv(s1 - s1b)))
    tests = [(), (1,), (-1,), (2,), (0, 1), (0, -1), (3, -2), (0, 0, 1), (1, 0, 1), (0, 1, -1), (-1, 0, 0, 1), (0, 0, 0, 2)]
    if EINST == "B":   # Gegenprobe: Teilmenge (mpmath ist ~100-mal langsamer)
        tests = [(0, -1), (0, 0, 1)] if TEIL2 else [(), (1,), (0, -1), (0, 0, 1)]
        WEGE.pop("b3")
    fehler = 0
    gesamt = 0
    for code in tests:
        u = dekodiere(code)
        kc = norm(kodiere(u))
        zeile = [f"Blatt {code}: u = {complex(u):.10f}, code zurueck {kc}"]
        if kc != norm(code):
            fehler += 1
        for name, weg in WEGE.items():
            u2 = verfolge(u, weg)
            k2 = norm(kodiere(u2))
            erw = vorhersage(name, code)
            ok = (k2 == erw)
            gesamt += 1
            if not ok:
                fehler += 1
            zeile.append(f"  {name}: -> {k2} (erwartet {erw}) {'ok' if ok else 'FALSCH'}  |du|={float(absv(u2-u)):.3e}")
        print("\n".join(zeile))
        sys.stdout.flush()
    print(f"Fehler: {fehler} von {gesamt} Schleifenproben (+ Kodierproben)")
