"""b31: Explizite Green-Funktion g von C ∖ Ḡ', Ḡ' = D̄(0,|L|) ∪ D̄(a0,|L−a0|) (beide Kreise durch L und L̄).
M(ζ) = (ζ − L)/(ζ − L̄) bildet beide Kreise auf Geraden durch 0 ab, das Äussere von Ḡ' auf einen Keil der Öffnung ω = 158.93°;
W = (e^{−iβ}M)^{π/ω} bildet den Keil auf die obere Halbebene ab; g = log|(W − conj W∞)/(W − W∞)|, W∞ = W(∞).
Prüfungen: g = 0 auf ∂G', g ~ log|ζ| für |ζ| → ∞, g(L + u e^{iψ}) ~ c·u^{π/ω}, Vergleich e^{−g} mit Cs 1/r_W."""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
aL = abs(L); a0 = mp.mpf("0.498563287941114434679619092493133305042"); R0 = abs(L - a0)
M = lambda z: (z - L) / (z - mp.conj(L))
# Bildgeraden: Richtung des Bildes eines weiteren Punktes jedes Kreises; Äusseres enthält M(∞) = 1
d1 = mp.arg(M(-aL)); d2 = mp.arg(M(a0 + R0))
# Randstrahlen = Bilder der äusseren Bögen: Richtung arg M(−|L|) (∂D ausserhalb D(a0)) und arg M(a0+R0) (∂D(a0) ausserhalb D)
lo, hi = sorted([d1, d2])
if lo < 0 < hi:          # Keil (lo, hi) enthält arg M(∞) = 0
    beta, omega = lo, hi - lo
else:                    # sonst der Komplementkeil (hi, lo + 2π)
    beta, omega = hi, lo + 2 * mp.pi - hi
def W(z):
    m = M(z) * mp.expj(-beta)
    a = mp.arg(m)
    if a < -mp.mpf(10) ** -20: a += 2 * mp.pi
    elif a < 0: a = mp.mpf(0)
    return abs(m) ** (mp.pi / omega) * mp.expj(a * mp.pi / omega)
Winf = W(mp.mpf(10) ** 12)
def g(z):
    w = W(z)
    return mp.log(abs((w - mp.conj(Winf)) / (w - Winf)))
# Prüfungen
print("g auf ∂D (Punkte ausserhalb D(a0)):", [mp.nstr(g(aL * mp.expjpi(t)), 3) for t in (0.6, 0.8, 1.0, 1.2, 1.4)])
print("g auf ∂D(a0) (ausserhalb D):", [mp.nstr(g(a0 + R0 * mp.expjpi(t)), 3) for t in (-0.3, -0.1, 0.0, 0.1, 0.3)])
for R in (10, 100, 1000):
    print(f"g({R}) − log {R} =", mp.nstr(g(mp.mpf(R)) - mp.log(R), 8), " (→ −log cap)")
cap = mp.exp(-(g(mp.mpf(10) ** 6) - mp.log(mp.mpf(10) ** 6)))
print("Kapazität von Ḡ' ≈", mp.nstr(cap, 8), "  (|L| =", mp.nstr(aL, 6), ")")
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
for u in (1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2):
    z = L + mp.mpf(u) * mp.expj(psi)
    print(f"u={u}: g = {mp.nstr(g(z), 6)}, g/u^(π/ω) = {mp.nstr(g(z) / mp.mpf(u) ** (mp.pi / omega), 6)},  e^(−g) = {mp.nstr(mp.exp(-g(z)),6)}, 1/r_W = {mp.nstr(max(aL/abs(z), R0/abs(z-a0)),6)}")
for z in (mp.mpf(2.2), mp.mpf(-2.2), mp.mpc(2.4, 2.4), mp.mpf(1.9), mp.mpf(-1.6)):
    print(f"ζ={mp.nstr(z,3)}: e^(−g) = {mp.nstr(mp.exp(-g(z)),5)}  vs  1/r_W = {mp.nstr(max(aL/abs(z), R0/abs(z-a0)),5)}  (C misst Zerfall nach 1/r_W)")
