"""k18 — Herkunft der Dekorrelation der adjungierten w-Köpfe r̂_n: Phase und Amplitude ihrer Erzeugenden auf |w| = 1.
G_n(θ) := Σ_i r̂_{n,i} e^{−iiθ} (w-Seite; bis auf den glatten Faktor k'w/k die adjungierte Kopffunktion 𝒜_n auf γ, K1″(b)).
Gemessen: Phaseninkrement Φ_n(θ) = arg(G_n/G_{n−1}) über θ (stationäre Punkte?), Amplitudenprofil |G_n(θ)|, Abklingen von |cos(r̂_n, r̂_{n−k})|.
Daten: k13_koepfe_2d_kubisch_300.npz."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = np.load("k13_koepfe_2d_kubisch_300.npz"); C = d["c"]; R = d["r"]; Nr = C.shape[0] - 1
th = np.linspace(-np.pi, np.pi, 721)
E = np.exp(-1j * np.outer(np.arange(Nr + 1), th))       # e^{−iiθ}
G = R @ E                                                 # Zeilen n: G_n(θ)
argL = np.degrees(np.angle(0.31813150520476413 + 1.3372357014306895j))
for n in [200, 250, 300]:
    a = np.abs(G[n]); ph = np.angle(G[n] / G[n - 1])
    w = a ** 2 / (a ** 2).sum()
    order = np.argsort(-w)
    top = order[:5]
    print("n = %d: Amplitude |G_n| — Masse in |θ| ≤ 30°: %.3f, um ±arg L (±15°): %.3f; Maxima bei θ = %s°"
          % (n, w[np.abs(th) <= np.radians(30)].sum(), w[np.abs(np.abs(np.degrees(th)) - argL) <= 15].sum(), ", ".join("%.0f" % np.degrees(th[i]) for i in top)))
    # Phaseninkrement entlang θ an Stellen mit nennenswerter Amplitude
    sel = w > 1e-4
    dph = np.unwrap(ph[sel])
    print("          Phaseninkrement Φ_n(θ) auf dem Träger: von %.3f bis %.3f rad; dΦ/dθ (Median |.|) = %.3f; Anteil |dΦ/dθ| < 0.05: %.3f"
          % (dph.min(), dph.max(), np.median(np.abs(np.gradient(dph, th[sel]))), np.mean(np.abs(np.gradient(dph, th[sel])) < 0.05)))
nr = np.linalg.norm(R, axis=1); cosr = (R @ R.T) / np.outer(nr, nr)
nc = np.linalg.norm(C, axis=1); tau = nc * nr
ks = [1, 2, 3, 4, 5, 7, 10, 14, 20, 28, 40, 56, 80, 113, 160]
mx = [max(abs(cosr[n, n - k]) for n in range(200, Nr + 1)) for k in ks]
print("max_{n∈[200,300]} |cos(r̂_n, r̂_{n−k})|: " + "  ".join("k=%d: %.3f" % (k, v) for k, v in zip(ks, mx)))
B0 = max(sum(np.sqrt(np.abs(cosr[n, :Nr + 1]))) for n in range(Nr + 1))
print("B₀ = sup_n Σ_m |cos^r_nm|^{1/2} = %.3f;  τ̄ = sup_n τ_n = %.4f (bei n = %d), τ_300 = %.4f" % (B0, tau.max(), np.argmax(tau), tau[Nr]))
for N in [100, 200, 300]:
    A2 = max(sum(np.sqrt(tau[n] * tau[m] * abs(cosr[n, m])) for m in range(N + 1)) for n in range(N + 1))
    print("N = %d: Schranke K7  τ̄·√((N+1)·B₀) = %.3f;  feiner √(A₁A₂) mit A₁ ≤ sup √τ_n Σ√τ_m: A₁ ≤ %.3f, A₂ = %.3f ⇒ %.3f"
          % (N, tau[:N + 1].max() * np.sqrt((N + 1) * B0), max(np.sqrt(tau[n]) * np.sqrt(tau[:N + 1]).sum() for n in range(N + 1)), A2,
             np.sqrt(max(np.sqrt(tau[n]) * np.sqrt(tau[:N + 1]).sum() for n in range(N + 1)) * A2)))
