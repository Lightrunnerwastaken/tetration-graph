"""Strang V — Geometrie und Green-Funktion von K_γ = (D̄ ∪ ℬ̄) ∪ f(D̄ ∪ ℬ̄) für Satz K1 (Auftrag Leitung 23:20, Punkt 1), Gleitkomma (mk_kneser).
Normierte Variable z = x/ρ. Schritte:
  1. Rand: Linse z = e^{s(φ)}e^{iφ} über Γin, f(Γout) = F(e^{iφ}), |φ| < θ; Prüfungen |F| > 1 auf Γout, |F| < 1 auf Γin und auf der Linse.
  2. Sternförmigkeit (Randwinkel monoton) ⇒ Radiusfunktion R(ϑ).
  3. Äussere Theodorsen-Abbildung ψ: {|w| > 1} → Ĉ∖K_γ, ψ(w) = w·exp(H(w)), H(w) = Σ h_n w^{−n}; g_γ = log|ψ⁻¹|.
  4. Lip(g_γ) = max_∂ 1/|ψ'| (Maximumprinzip für |Φ'/Φ|), in x-Einheiten /ρ.
  5. Schlitz η = L(1+u): g_γ(η) gegen c_g u^{1/α'} (α' = 0.787172), dist(η, K_γ)/(ρu) =: c_d(u).
  6. M_γ = max |m| auf ∂(D ∪ ℬ) (Maximumprinzip), m = e^{ρz}/F'(z).
  7. Niveaulinie Γ_d = ψ(e^d·∂𝔻): Länge (x-Einheiten), max|x|, Abstand zu D̄ (für den Konturteil von Satz P-I).
Aufruf: python v30_geometrie_Kgamma.py S [tau1] [p] [M]
"""
import math
import sys
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
import mk_kneser as MK  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
S = float(sys.argv[1])
TAU1 = math.radians(float(sys.argv[2])) if len(sys.argv) > 2 else math.radians(12)
PEXP = float(sys.argv[3]) if len(sys.argv) > 3 else 1.3
M = int(sys.argv[4]) if len(sys.argv) > 4 else 2 ** 15
rho, th = MK.RHO, MK.TH
ALPHA = 0.787172


def s_lin(phi):
    t = np.minimum(np.mod(phi - th, 2 * np.pi), np.mod(-phi - th, 2 * np.pi))      # Abstand zu ±θ in Γin-Richtung
    inn = (np.mod(phi, 2 * np.pi) > th) & (np.mod(phi, 2 * np.pi) < 2 * np.pi - th)
    return np.where(inn, S * (1 - np.exp(-(t / TAU1) ** PEXP)), 0.0)


def F(z):
    return np.asarray(MK.F(np.atleast_1d(np.asarray(z, dtype=complex))))


def Fp(z):
    return np.asarray(MK.Fp(np.atleast_1d(np.asarray(z, dtype=complex))))


# 1. Rand
nb = 6000
phi_out = np.linspace(-th, th, nb)[1:-1]
f_out = F(np.exp(1j * phi_out))
phi_in = np.linspace(th, 2 * np.pi - th, nb)[1:-1]
z_lin = np.exp(s_lin(phi_in)) * np.exp(1j * phi_in)
f_in = F(np.exp(1j * phi_in))
f_lin = F(z_lin)
print(f"[Linse S = {S}, τ₁ = {math.degrees(TAU1):.0f}°, p = {PEXP}] min|F| auf Γout = {np.min(np.abs(f_out)):.4f} (> 1?), "
      f"max|F| auf Γin = {np.max(np.abs(f_in)):.4f}, max|F| auf der Linse = {np.max(np.abs(f_lin)):.4f} (< 1 ⇒ f(ℬ) ⊂ D)", flush=True)
# Abstand von f(ℬ) zum Rand ∂D, aufgelöst nach τ
tau = np.degrees(phi_in - th)
for tg in (1, 3, 10, 30, 103):
    i = np.argmin(np.abs(tau - tg))
    print(f"   τ = {tg}°: 1 − |F(Linse)| = {1 - abs(f_lin[i]):.4e}, Linsenhöhe s = {s_lin(phi_in[i]):.4e}", flush=True)
# 2. Rand von K_γ: Linse (Γin-Seite) + f(Γout); Sternform
rand = np.concatenate([f_out, z_lin])
ang = np.unwrap(np.angle(rand))
order = np.argsort(np.mod(np.angle(rand), 2 * np.pi))
a_sorted = np.mod(np.angle(rand), 2 * np.pi)[order]
r_sorted = np.abs(rand)[order]
mono_out = np.all(np.diff(np.unwrap(np.angle(f_out))) > 0)
mono_in = np.all(np.diff(np.unwrap(np.angle(z_lin))) > 0)
print(f"Sternform: Randwinkel monoton auf f(Γout): {mono_out}, auf der Linse: {mono_in}; Winkel f(Γout) von {math.degrees(np.angle(f_out[0])):.2f}° "
      f"bis {math.degrees(np.angle(f_out[-1])):.2f}°", flush=True)


def R_of(theta):
    return np.interp(np.mod(theta, 2 * np.pi), a_sorted, r_sorted, period=2 * np.pi)


# 3. äussere Theodorsen
t = 2 * np.pi * np.arange(M) / M
n = np.fft.fftfreq(M, d=1.0 / M)
theta = t.copy()
for it in range(400):
    lr = np.log(R_of(theta))
    gh = np.fft.fft(lr)
    conj = np.real(np.fft.ifft(1j * np.sign(n) * gh))          # Aussenproblem: Vorzeichen umgekehrt
    neu = t + conj
    dd = np.max(np.abs(neu - theta))
    theta = 0.5 * theta + 0.5 * neu                               # gedämpft (Ecken)
    if dd < 1e-13:
        break
lr = np.log(R_of(theta))
hh = np.fft.fft(lr) / M                                           # Re H auf dem Rand: Σ (h_n e^{−int} + konj.)
# H(w) = h0 + Σ_{n≥1} 2·ĥ_{−n}... direkt: H(e^{it}) = lr(t) + i·conj(t)
Hb = lr + 1j * np.real(np.fft.ifft(1j * np.sign(n) * np.fft.fft(lr)))
psi_b = np.exp(1j * t) * np.exp(Hb)
res = np.max(np.abs(np.abs(psi_b) - R_of(np.angle(psi_b))))
cap = math.exp(hh[0].real)
print(f"Theodorsen aussen: {it} Iterationen, letzte Änderung {dd:.1e}, Randfehler max||ψ| − R(arg ψ)| = {res:.1e}; Kapazität (z-Einheiten) = {cap:.6f}", flush=True)
# Koeffizienten von H in w^{−n}: H(w) = c_0 + Σ_{n≥1} c_n w^{−n}
ch = np.fft.fft(Hb) / M
c = np.zeros(M // 2, dtype=complex)
c[0] = ch[0]
c[1:] = ch[-1:-(M // 2):-1]                                       # Fourier-Index −n ↔ w^{−n}


def psi(w):
    w = np.atleast_1d(w)
    H = np.polynomial.polynomial.polyval(1 / w, c)
    return w * np.exp(H)


def psip(w):
    w = np.atleast_1d(w)
    iw = 1 / w
    H = np.polynomial.polynomial.polyval(iw, c)
    dH = -np.polynomial.polynomial.polyval(iw, c[1:] * np.arange(1, len(c))) * iw ** 2
    return np.exp(H) * (1 + w * dH)


# 4. Lipschitz (Rand: Werte auf |w| = 1 + 1e-6 statt genau 1, Ecken haben |ψ'| → ∞ bzw. 0)
wb = (1 + 1e-9) * np.exp(1j * t[::4])
pp = np.abs(psip(wb))
lip_z = np.max(1 / pp)
print(f"Lip(g_γ) = max 1/|ψ'| = {lip_z:.4f} (z-Einheiten) = {lip_z / rho:.4f} pro x-Einheit; Stelle arg ψ = "
      f"{math.degrees(np.angle(psi(wb[np.argmax(1 / pp)])[0])):.1f}°", flush=True)


# 5. Schlitz
def u_von_g(g):
    """Schnitt der Niveaulinie {g_γ = g} = ψ(e^g ∂𝔻) mit dem Strahl arg z = θ: liefert u mit z = e^{iθ}(1 + u)."""
    tt = 2 * np.pi * np.arange(8192) / 8192
    zl = psi(math.exp(g) * np.exp(1j * tt))
    a_ = np.angle(zl * np.exp(-1j * th))
    i = np.where((a_[:-1] < 0) & (a_[1:] >= 0))[0]
    i = i[np.argmin(np.abs(np.abs(zl[i]) - 1))]
    lam = -a_[i] / (a_[i + 1] - a_[i])
    zs = zl[i] + lam * (zl[i + 1] - zl[i])
    return abs(zs) - 1


GG = np.exp(np.linspace(math.log(1e-5), math.log(0.45), 60))
UU = np.array([u_von_g(g) for g in GG])


def g_von(z):
    u = abs(z) - 1
    return float(np.exp(np.interp(math.log(u), np.log(UU), np.log(GG)))), None


def dist_K(z):
    return np.min(np.abs(rand - z))


zeile = []
cg_min, cd_min = 1e9, 1e9
for u in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3):
    z = np.exp(1j * th) * (1 + u)
    g, w = g_von(z)
    cg = g / u ** (1 / ALPHA)
    cd = dist_K(z) / u
    cg_min, cd_min = min(cg_min, cg), min(cd_min, cd)
    zeile.append(f"u = {u}: g = {g:.4e}, g/u^(1/α') = {cg:.4f}, dist/(ρu) = {cd:.4f}")
print("Schlitz L(1+u):\n   " + "\n   ".join(zeile) + f"\n   ⇒ c_g ≥ {cg_min:.4f}, c_d ≥ {cd_min:.4f} (u ∈ [0.001, 0.3]); u-Bereich der Tabelle [{UU.min():.1e}, {UU.max():.2f}]", flush=True)
# 6. M_γ
zz = np.concatenate([np.exp(1j * phi_out), z_lin])
mm = np.abs(np.exp(rho * zz) / Fp(zz))
print(f"M_γ = max|m| auf ∂(D ∪ ℬ) = {np.max(mm):.4f} (auf Γout {np.max(mm[:len(phi_out)]):.4f}, Linse {np.max(mm[len(phi_out):]):.4f})", flush=True)
# 7. Niveaulinien
for d in (0.1, 0.15, 0.2, 0.3):
    wl = math.exp(d) * np.exp(1j * t[::2])
    zl = psi(wl)
    laenge = np.sum(np.abs(np.diff(np.concatenate([zl, zl[:1]])))) * rho
    dD = np.min(np.abs(zl)) - 1
    print(f"Γ_d, d = {d}: Länge {laenge:.3f} (x), max|x| = {np.max(np.abs(zl)) * rho:.3f}, min|x| − ρ = {dD * rho:.4f}, "
          f"Schnitt mit Schlitz bei u₀ = {np.min(np.abs(zl[np.argmin(np.abs(np.angle(zl) - th))]) - 1):.4f}", flush=True)
# 8. Speichern für die N₀-Bilanz (v31)
import json  # noqa: E402
lv = {}
for d in (0.05, 0.07, 0.09, 0.1, 0.12, 0.15, 0.2, 0.3):
    wl = math.exp(d) * np.exp(1j * t[::2])
    zl = psi(wl)
    i0 = np.argmin(np.abs(np.angle(zl) - th))
    lv[str(d)] = dict(laenge=float(np.sum(np.abs(np.diff(np.concatenate([zl, zl[:1]])))) * rho), max_abs=float(np.max(np.abs(zl)) * rho),
                      min_abs_minus_rho=float((np.min(np.abs(zl)) - 1) * rho), u0=float(abs(zl[i0]) - 1),
                      punkte=[[float(v.real * rho), float(v.imag * rho)] for v in zl[::64]])
out = dict(S=S, tau1=math.degrees(TAU1), p=PEXP, lip_x=float(lip_z / rho), M_gamma=float(np.max(mm)), c_d=float(cd_min),
           g_tabelle=dict(u=[float(x) for x in UU], g=[float(x) for x in GG]), niveaulinien=lv, kapazitaet_z=cap)
(HIER / f"v30_Kgamma_S{S}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(f"gespeichert: v30_Kgamma_S{S}.json", flush=True)
