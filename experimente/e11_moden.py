"""exp.moden — Modenprofile der Phasen Φ_(e,b) aus gp-tetration; Universalität und m(2)-Fassungen.

(1) Normierte Amplituden |a_k|/|a_1| für alle Basen: nahezu gleiche Form (Universalität, bw.beob.modenabfall).
(2) Drei gespeicherte Fassungen von m(2) = µ_(e,2) vergleichen (Versionskonflikt-Suche).
"""
import mpmath as mp

from gemeinsam import Lauf, gp_json, plt

mp.mp.dps = 60
lauf = Lauf("exp.moden", __file__)
pm = gp_json("phi_modes.json")
hub = gp_json("mu_hub.json")["values"]
basen = pm["bases"]
profile = {}
for b, dat in basen.items():
    amp = [abs(mp.mpc(m["re"], m["im"])) for m in dat["modes"]]
    profile[b] = [a / amp[0] for a in amp]
k_max = min(len(p) for p in profile.values())
spreiz = []
for k in range(1, 6):
    w = [mp.log(p[k]) for p in profile.values()]
    spreiz.append(max(w) - min(w))
lauf.pruefe("normierte Profile k = 2..6: Spannweite von ln(|a_k|/|a_1|) < 0.25 über alle Basen", max(spreiz) < 0.25,
            [mp.nstr(s, 3) for s in spreiz])

m2_log = mp.mpf("-1.12840377662892400987921790203642626711169249485385417303588")   # Runde 1, 180 Stellen (gekürzt)
m2_modes = mp.mpf(basen["2"]["mu"])
m2_hub = mp.mpf(hub["2"])
d1 = abs(m2_modes - m2_log)
d2 = abs(m2_hub - m2_log)
lauf.pruefe("phi_modes.json: m(2) stimmt mit dem 180-Stellen-Wert auf 1e-30", d1 < 1e-30, mp.nstr(d1, 3))
lauf.pruefe("mu_hub.json: m(2) stimmt mit dem 180-Stellen-Wert auf 1e-30", d2 < 1e-30, mp.nstr(d2, 3))
lauf.daten = {"spannweite_ln": [mp.nstr(s, 4) for s in spreiz], "m2_abweichung_phi_modes": d1, "m2_abweichung_mu_hub": d2,
              "basen": list(basen)}

fig, ax = plt.subplots(figsize=(6.5, 4))
for b, p in profile.items():
    ax.semilogy(range(1, len(p) + 1), [float(x) for x in p], "o-", ms=3, label=f"b = {b}")
ax.set_xlabel("Mode k")
ax.set_ylabel("|a_k| / |a_1|")
ax.set_title("Normierte Modenprofile von Φ_(e,b) (gp-tetration)")
ax.legend(fontsize=7, ncol=2)
lauf.bild(fig)
lauf.zusammenfassung = (f"Profilform über {len(profile)} Basen nahezu gleich (Spannweite ln ≤ {mp.nstr(max(spreiz), 2)} für k ≤ 6). "
                        f"m(2): phi_modes weicht {mp.nstr(d1, 2)}, mu_hub {mp.nstr(d2, 2)} vom 180-Stellen-Wert ab.")
lauf.schreibe()
