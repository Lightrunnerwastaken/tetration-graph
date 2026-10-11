"""O15 — Aufgabe 2: Anregung der Eckmoden-Paare durch Pole (kubisches Modell), aus gespeicherten Singulärvektoren (o07b_kubisch_N.npz).

Y_N = J_N^{-1} M_N = Σ σ_i u_i v_i^* (ℓ²_ρ).  Pol-Antwort Π_N k_η = J_N^{-1} w_η = Y_N M_N^{-1} w_η,
M_N^{-1} w_η = P_N[ 1/(m(x)(η − f(x))) + 1/(η − x) ]  (exakt: M^{-1}P_N[m g] = P_N g).
Eckanteil von Modus i:  e_i(η) = σ_i |⟨v_i, c(η)⟩|,  c(η) = ρ-skalierte Koeffizienten von M^{-1}w_η.
Pole η: (a) Γin-Seite von ∂K' (= ∂D) bei arg L + τ; (b) f(Γout)-Seite von ∂K': η = f(ρ e^{i(arg L − τ)}); (c) Schlitz η = L(1+u);
(d) Lemma-D6-Kontrolle: Profil der adjungierten Moden σ_i|Ṽ_i(arg L + τ)| auf |x| = ρ (beidseitig).
Aufruf: python o15_anregung.py N1,N2
"""
import sys

import mpmath as mp
import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
Ns = [int(s) for s in sys.argv[1].split(",")]
Lc = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
rho = abs(Lc)
thL = mp.arg(Lc)
from flint import ctx
ctx.prec = 200
fc = [mp.mpf(x.mid().str(50, radius=False)) for x in K.kubik_koeff()]          # f-Koeffizienten (reell)
f = lambda x: sum(c * x ** i for i, c in enumerate(fc))
fp_c = [i * c for i, c in enumerate(fc)][1:]
# m(x) = f'(f(x)) als Polynom (Grad 6)
def polmul(a, b):
    r = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r
mpoly = [mp.mpf(0)]
pw = [mp.mpf(1)]
for i, c in enumerate(fp_c):
    term = [c * t for t in pw]
    mpoly = [ (mpoly[k] if k < len(mpoly) else 0) + (term[k] if k < len(term) else 0) for k in range(max(len(mpoly), len(term)))]
    pw = polmul(pw, fc)


def koeff(eta, N):
    """ρ-skalierte Taylor-Koeffizienten (k ≤ N) von 1/(m(x)(η − f(x))) + 1/(η − x)."""
    q = polmul(mpoly, [eta - fc[0]] + [-c for c in fc[1:]])
    n = N + 1
    c = [mp.mpc(0)] * n
    c[0] = 1 / q[0]
    for k in range(1, n):
        s = 0
        for j in range(1, min(k, len(q) - 1) + 1):
            s += q[j] * c[k - j]
        c[k] = -s / q[0]
    out = np.array([complex(c[k] + 1 / eta ** (k + 1)) * float(rho) ** k for k in range(n)])
    return out


pole = []
for tau in (0.5, 2, 5, 10, 20, 40, 90):
    pole.append((f"Γin τ=+{tau}°", rho * mp.expj(thL + mp.radians(tau))))
for tau in (0.5, 2, 5, 10, 20, 40):
    pole.append((f"f(Γout) τ=−{tau}°", f(rho * mp.expj(thL - mp.radians(tau)))))
for u in (0.001, 0.003, 0.01, 0.03, 0.1):
    pole.append((f"Schlitz u={u}", Lc * (1 + u)))
pole.append(("Ecke L", Lc * (1 + mp.mpf("1e-9"))))

daten = {}
for N in Ns:
    d = np.load(f"o07b_kubisch_{N}.npz")
    s, V = d["s"], d["V"]
    nm = min(6, V.shape[1])
    daten[N] = {}
    for name, eta in pole:
        c = koeff(eta, N)
        e = [s[i] * abs(np.vdot(V[:, i], c)) for i in range(nm)]
        daten[N][name] = e
    # Profil der adjungierten Moden auf |x| = ρ: Ṽ_i(φ) = Σ_k v_{i,k} e^{−ikφ}
    daten[N]["_s"] = s[:nm]
    daten[N]["_V"] = V[:, :nm]

N1, N2 = Ns[0], Ns[-1]
print(f"σ (N={N1}): {np.array2string(daten[N1]['_s'], precision=4)}\nσ (N={N2}): {np.array2string(daten[N2]['_s'], precision=4)}")
print(f"Raten von σ_i ({N1}→{N2}): {np.array2string(np.log(daten[N2]['_s'] / daten[N1]['_s']) / (N2 - N1), precision=4)}")
print(f"\nEckanteil e_i(η) = σ_i|⟨v_i, M⁻¹w_η⟩| bei N = {N2} und Rate ({N1}→{N2}) je Modus i = 1..6:")
for name, eta in pole:
    e1, e2 = np.array(daten[N1][name]), np.array(daten[N2][name])
    rate = np.log(e2 / e1) / (N2 - N1)
    print(f"  {name:18s}: " + "  ".join(f"{a:9.3g} ({r:+.4f})" for a, r in zip(e2, rate)))
print("\nProfil σ_i|Ṽ_i(arg L + τ)| (|x| = ρ), N = {0}, und Rate ({1}→{0}) — τ > 0 Γin-Seite, τ < 0 Γout-Seite:".format(N2, N1))
for tau in (-30, -8, -2, -0.5, 0, 0.5, 2, 8, 30, 103.4):
    zeile = []
    for i in range(6):
        vals = []
        for N in (N1, N2):
            V = daten[N]["_V"][:, i]
            k = np.arange(len(V))
            phi = float(thL) + np.radians(tau)
            vals.append(daten[N]["_s"][i] * abs(np.sum(V * np.exp(-1j * k * phi))))
        zeile.append(f"{vals[1]:9.3g} ({np.log(vals[1]/vals[0])/(N2-N1):+.4f})")
    print(f"  τ={tau:+7.1f}°: " + "  ".join(zeile))
