"""exp.fixpunkte — Fixpunkte von b^z über Lambert W, Multiplikatoren, Minimum von b^x − x."""
import mpmath as mp

from gemeinsam import Lauf, plt

mp.mp.dps = 60
lauf = Lauf("exp.fixpunkte", __file__)
eta = mp.e ** (1 / mp.e)
tabelle = []
for b in [mp.mpf("1.2"), mp.mpf("1.3"), mp.sqrt(2), mp.mpf("1.44"), mp.mpf("1.5"), mp.mpf(2), mp.e, mp.mpf(3), mp.mpf(10)]:
    lb = mp.log(b)
    for k in (0, -1):
        L = -mp.lambertw(-lb, k) / lb
        lam = L * lb
        res = abs(mp.power(b, L) - L)
        tabelle.append({"b": mp.nstr(b, 8), "k": k, "L": mp.nstr(L, 12), "lambda": mp.nstr(lam, 12),
                        "abs_lambda": float(abs(lam)), "residuum": mp.nstr(res, 3)})
        lauf.pruefe(f"b={mp.nstr(b, 5)}, W_{k}: b^L = L", res < mp.mpf(10) ** -50, mp.nstr(res, 3))
        if b > eta:
            lauf.pruefe(f"b={mp.nstr(b, 5)}, W_{k}: |λ| > 1 (abstossend)", abs(lam) > 1, mp.nstr(abs(lam), 8))
        elif k == 0:
            lauf.pruefe(f"b={mp.nstr(b, 5)}: unterer Fixpunkt reell, 0 < λ < 1", abs(mp.im(L)) < 1e-40 and 0 < mp.re(lam) < 1,
                        mp.nstr(mp.re(lam), 8))
        else:
            lauf.pruefe(f"b={mp.nstr(b, 5)}: oberer Fixpunkt reell, λ > 1", abs(mp.im(L)) < 1e-40 and mp.re(lam) > 1,
                        mp.nstr(mp.re(lam), 8))

# Minimum von b^x − x gegen Formel (1 + ln ln b)/ln b
for b in [mp.mpf("1.2"), eta, mp.mpf(2), mp.mpf(10)]:
    lb = mp.log(b)
    xs = -mp.log(lb) / lb
    num = mp.findroot(lambda x: mp.power(b, x) * lb - 1, xs + mp.mpf("0.1"))
    wert = mp.power(b, num) - num
    formel = (1 + mp.log(lb)) / lb
    lauf.pruefe(f"b={mp.nstr(b, 6)}: min(b^x − x) = (1+ln ln b)/ln b", abs(wert - formel) < mp.mpf(10) ** -50,
                f"{mp.nstr(wert, 10)}")
L_oben = -mp.lambertw(-1, -1)
L_unten = -mp.lambertw(-1, 0)
lauf.pruefe("L_e = −W_{−1}(−1) = 0.3181315… + 1.3372357… i (Im > 0)",
            abs(L_oben - mp.mpc("0.31813150520476413531", "1.33723570143068940890")) < 1e-18, mp.nstr(L_oben, 15))
lauf.pruefe("−W_0(−1) ist der konjugierte Fixpunkt (Im < 0)", abs(L_unten - mp.conj(L_oben)) < 1e-40, mp.nstr(L_unten, 15))

# Bild: |λ| über b und Fixpunktkurve
bs = [1.01 + i * 0.005 for i in range(1, 90)] + [1.46 + 0.02 * i for i in range(1, 200)]
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
unten, oben, komplex = [], [], []
for b in bs:
    lb = mp.log(b)
    l0 = -mp.lambertw(-lb, 0)
    l1 = -mp.lambertw(-lb, -1)
    if b < float(eta):
        unten.append((b, float(mp.re(l0))))
        oben.append((b, float(mp.re(l1))))
    else:
        komplex.append((b, float(abs(l0)), complex(-mp.lambertw(-lb, -1) / lb)))
ax[0].plot([u[0] for u in unten], [u[1] for u in unten], label="unterer Fixpunkt: λ (anziehend)")
ax[0].plot([u[0] for u in oben], [u[1] for u in oben], label="oberer Fixpunkt: λ (abstossend)")
ax[0].plot([c[0] for c in komplex], [c[1] for c in komplex], label="primärer Fixpunkt L: |λ|")
ax[0].axvline(float(eta), ls="--", color="#888")
ax[0].axhline(1, ls=":", color="#888")
ax[0].set_xscale("log")
ax[0].set_xlabel("Basis b")
ax[0].set_ylabel("Multiplikator")
ax[0].set_title("Multiplikatoren; η = e^{1/e} gestrichelt")
ax[0].legend(fontsize=8)
ax[1].plot([c[2].real for c in komplex], [c[2].imag for c in komplex], color="#8c4fa3")
ax[1].plot([c[2].real for c in komplex], [-c[2].imag for c in komplex], color="#8c4fa3", alpha=0.4)
for b, n in [(2, "2"), (float(mp.e), "e"), (10, "10")]:
    L = complex(-mp.lambertw(-mp.log(b), -1) / mp.log(b))
    ax[1].plot([L.real], [L.imag], "o", color="#b5523b")
    ax[1].annotate(f"L_{n}", (L.real, L.imag), textcoords="offset points", xytext=(5, 4))
ax[1].plot([float(mp.e)], [0], "s", color="#2f8f83")
ax[1].annotate("e (b = η)", (float(mp.e), 0), textcoords="offset points", xytext=(5, -12))
ax[1].set_xlabel("Re L")
ax[1].set_ylabel("Im L")
ax[1].set_title("Primärer Fixpunkt L_b für b > η")
lauf.bild(fig)
lauf.daten = {"tabelle": tabelle}
lauf.zusammenfassung = "Fixpunkte über W_0/W_{−1} mit Residuum < 1e-50; für b > η alle geprüften |λ| > 1; Minimumsformel exakt bestätigt."
lauf.schreibe()
