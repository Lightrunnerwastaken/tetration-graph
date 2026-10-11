"""b35 (Runde 3): Strenge Prüfungen (arb-Intervallrechnung) für Lemma 7.4.
(a) Auf dem äusseren Bogen C_a = {a0 + R0 e^{iφ}} ausserhalb D und ausserhalb D(L, r): |Log z| < |L|  (⇒ z ∈ exp(D)).
    Zusammen mit Subharmonizität von |Log z|² und |Log z| ≤ |L| auf ∂D (Gleichheit nur in L, L̄; da Log L = L) folgt
    D(a0,R0) ∖ (D̄ ∪ D(L,r) ∪ D(L̄,r)) ⊂ exp(D) kompakt.  Prüfung auch auf den Bögen ∂D(L,r) im Halbmond.
(b) Linker Bogen ∂D ausserhalb D(a0) und ausserhalb D(L,r): Re z ≤ log|L| − c(r) (⇒ e^z ∈ D kompakt) — streng;
    f_K(e^z) ∉ (−∞, 0] — numerisch mit Knesers Reihe (4000 Koeffizienten)."""
import sys, json
from flint import arb, acb, ctx
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 200
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
L = acb(arb(mp.nstr(mp.re(Lm), 50)), arb(mp.nstr(mp.im(Lm), 50)))
# Fehlerkugel für L (Fixpunkt e^L = L auf 1e-45 genau)
L = acb(arb(L.real.mid(), 1e-40), arb(L.imag.mid(), 1e-40))
aL = abs(L)
a0 = arb("0.498563287941114434679619092493133305042", 1e-38)
R0 = abs(L - a0)
phiL = (L - a0).arg()
r = arb("0.02")
print("|L| =", aL.str(15), " R0 =", R0.str(15), " φ_L =", (phiL * 180 / arb.pi()).str(10), "°")

def log_ok_on_arc(phi_lo, phi_hi, n):
    """prüft |Log z|² < |L|² für z = a0 + R0 e^{iφ}, φ ∈ [phi_lo, phi_hi], zerlegt in n Intervalle; ausgenommen |z − L| < r."""
    schlecht, maxq = 0, arb(0)
    for i in range(n):
        lo = phi_lo + (phi_hi - phi_lo) * i / n
        hi = phi_lo + (phi_hi - phi_lo) * (i + 1) / n
        ph = arb((lo + hi) / 2, abs(hi - lo) / 2)
        z = a0 + R0 * acb(ph.cos(), ph.sin())
        if abs(z - L).upper() < r.lower() or abs(z - L.conjugate()).upper() < r.lower():
            continue
        if abs(z).upper() < aL.lower():          # Punkt in D: nichts zu prüfen
            continue
        lz = z.log()
        q = (lz.real ** 2 + lz.imag ** 2) / aL ** 2
        maxq = max(maxq, q.upper()) if isinstance(maxq, arb) else q.upper()
        if not (q < 1):
            schlecht += 1
    return schlecht, maxq

pi = arb.pi()
s, mq = log_ok_on_arc(-phiL.upper(), phiL.upper(), 20000)
print(f"(a) äusserer Bogen C_a (|z−L| ≥ {r}): Intervalle mit möglicher Verletzung: {s}; max |Log z|²/|L|² ≤ {arb(mq).str(8)}")
# Bögen ∂D(L, r) im Halbmond D(a0,R0) \ D̄:
schlecht = 0; maxq = 0
n = 4000
for i in range(n):
    ph = arb((2 * i + 1) * pi / n - pi, pi / n)      # φ ∈ [−π, π]
    z = L + r * acb(ph.cos(), ph.sin())
    if abs(z).upper() < aL.lower() or abs(z - a0).lower() > R0.upper():
        continue                                     # sicher in D oder sicher ausserhalb D(a0)
    lz = z.log(); q = (lz.real ** 2 + lz.imag ** 2) / aL ** 2
    maxq = max(maxq, float(q.upper().str(10, radius=False)))
    if not (q < 1):
        schlecht += 1
print(f"(a') Bogen |z − L| = {r} im Halbmond: mögliche Verletzungen {schlecht}, max |Log z|²/|L|² ≤ {maxq:.6f}")
# (b) linker Bogen: z = |L| e^{iθ}, θ ∈ [θ_L, 2π − θ_L] ausserhalb D(L,r): Re z ≤ log|L| − c
thL = L.arg()
n = 20000; maxre = None
for i in range(n):
    lo = thL + (2 * pi - 2 * thL) * i / n; hi = thL + (2 * pi - 2 * thL) * (i + 1) / n
    th = arb((lo + hi) / 2, abs(hi - lo) / 2)
    z = aL * acb(th.cos(), th.sin())
    if abs(z - L).upper() < r.lower() or abs(z - L.conjugate()).upper() < r.lower():
        continue
    d = z.real - aL.log()
    maxre = d.upper() if maxre is None else max(maxre, d.upper())
print(f"(b) linker Bogen ausserhalb D(L,{r}): max (Re z − log|L|) ≤ {arb(maxre).str(8)}  (< 0 ⇒ e^z ∈ D kompakt)")
# (b') numerisch: f_K(e^z) weg von (−∞, 0]
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
Kc = [mp.mpf(c) for c in (d4["a"] if isinstance(d4, dict) else d4)][:2000]
mp.mp.dps = 30
def fK(w):
    s = mp.mpc(0)
    for c in reversed(Kc):
        s = s * w + c
    return s
minabst = 10; minarg = 10
for i in range(2001):
    th = mp.arg(Lm) + (2 * mp.pi - 2 * mp.arg(Lm)) * i / 2000
    z = abs(Lm) * mp.expj(th)
    if abs(z - Lm) < 0.02 or abs(z - mp.conj(Lm)) < 0.02:
        continue
    v = fK(mp.exp(z))
    minabst = min(minabst, abs(v)); minarg = min(minarg, mp.pi - abs(mp.arg(v)))
print(f"(b') numerisch: min |f_K(e^z)| = {mp.nstr(minabst, 5)}, min Abstand von arg zu ±π = {mp.nstr(minarg, 5)} rad (⇒ Log holomorph)")
