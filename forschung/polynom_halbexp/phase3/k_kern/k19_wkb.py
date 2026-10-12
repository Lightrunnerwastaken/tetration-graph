"""k19 — WKB-Form der adjungierten w-Köpfe: G_n(θ) = a_n(θ)·e^{inΦ(θ)} mit FESTER Phase Φ? Glattheit von a_n, Lage der Eckanteile.
G_n(θ) = Σ_i r̂_{n,i}e^{−iiθ}. Φ aus dem Phaseninkrement arg(G_n/G_{n−1}), gemittelt über zwei getrennte n-Fenster (n-Unabhängigkeit prüfen).
a_n = G_n e^{−inΦ}: Fourier-Abfall von a_n (glatt ⇔ schneller Abfall), Vergleich a_250 vs a_300 (n-Unabhängigkeit der Form).
Vorhersage Lemma K8: |⟨r̂_n, r̂_m⟩| = |(1/2π)∫ a_n ā_m e^{i(n−m)Φ}dθ| — mit den gemessenen a_n nachgerechnet und mit cos^r verglichen.
Daten: k13_koepfe_2d_kubisch_300.npz."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "k13_koepfe_2d_kubisch_300.npz"); R = d["r"]; Nr = R.shape[0] - 1
M = 2048
th = 2 * np.pi * np.arange(M) / M - np.pi
E = np.exp(-1j * np.outer(np.arange(Nr + 1), th))
G = R @ E
def phase_incr(lo, hi):
    z = sum(G[n] * np.conj(G[n - 1]) for n in range(lo, hi))
    return np.angle(z), np.abs(z)
P1, w1 = phase_incr(200, 231); P2, w2 = phase_incr(270, 301)
amp = np.abs(G[300]) ** 2; amp /= amp.sum()
sel = amp > 1e-5
dP = np.angle(np.exp(1j * (P1 - P2)))
print("Φ aus n∈[200,230] vs [270,300]: max|Δ| auf dem Träger = %.4f rad, gewichtet mittel %.4f rad" % (np.abs(dP[sel]).max(), (np.abs(dP) * amp).sum()))
Phi = np.unwrap(P2)
dPhi = np.gradient(Phi, th)
print("dΦ/dθ auf dem Träger: min %.3f, max %.3f, gewichteter Mittelwert %.3f; Anteil mit |dΦ/dθ| < 0.3: %.4f (gewichtet)"
      % (dPhi[sel].min(), dPhi[sel].max(), (dPhi * amp).sum(), amp[np.abs(dPhi) < 0.3].sum()))
# Stellen kleiner |dΦ/dθ| (fast stationär) — wo?
idx = np.where(sel & (np.abs(dPhi) < 0.5))[0]
if len(idx):
    print("  fast stationäre Phase bei θ ≈ %s° (Massenanteil %.4f)" % (", ".join("%.0f" % np.degrees(th[i]) for i in idx[:: max(1, len(idx) // 8)]), amp[idx].sum()))
a = {n: G[n] * np.exp(-1j * n * P2) for n in [200, 250, 300]}
for n in [200, 250, 300]:
    f = np.fft.fft(a[n]) / M
    p = np.abs(f) ** 2; p /= p.sum()
    k = np.fft.fftfreq(M, 1 / M)
    print("a_%d: Fourier-Masse bei |k| ≤ 5: %.4f, ≤ 20: %.4f, ≤ 60: %.4f, > 60: %.2e" % (n, p[np.abs(k) <= 5].sum(), p[np.abs(k) <= 20].sum(), p[np.abs(k) <= 60].sum(), p[np.abs(k) > 60].sum()))
def rel(u, v): return np.linalg.norm(u / np.linalg.norm(u) - v / np.linalg.norm(v) * (np.vdot(v, u) / abs(np.vdot(v, u))))
print("Formvergleich ‖â_250 − â_300‖ (normiert, Phase angeglichen) = %.4f, ‖â_200 − â_300‖ = %.4f" % (rel(a[250], a[300]), rel(a[200], a[300])))
# Eckanteil: Masse von |G_n|² in Fenstern um die Maxima, Breite c·n^{-Re a}
ra = 0.8935861640191741
pk = th[np.argmax(np.abs(G[300]) * (th > 0))]
for n in [150, 200, 250, 300]:
    g = np.abs(G[n]) ** 2; g /= g.sum()
    for cc in [1.0, 3.0]:
        win = cc * n ** (-ra)
        m = (np.abs(th - pk) <= win) | (np.abs(th + pk) <= win)
        print("n = %d: Eckfenster ±%.1f·n^{-Re a} = ±%.2f° um ±%.1f°: Massenanteil %.4f" % (n, cc, np.degrees(win), np.degrees(pk), g[m].sum()))
