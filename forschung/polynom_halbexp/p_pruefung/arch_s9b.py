"""Hilfsskript (Strang P): fügt in ARCHITEKTUR.md §9.3 den gemischten Weg (Rest in K', Eckteil in W̄) als bevorzugte Zielaussage ein."""
from pathlib import Path

p = Path(__file__).resolve().parent.parent / "ARCHITEKTUR.md"
s = p.read_text(encoding="utf-8")
anker = "**Zielaussagen (das ist zu beweisen; die Kette dahinter ist in 9.4 bewiesen):**"
neu = """**BEVORZUGT (P 00:58, nach W §10/§13): gemischter Weg — Rest in K', Eckteil in W̄ — braucht KEIN (E_A).** Die Zerlegung
y_ζ = Q_N y_ζ + (I−Q_N)y_ζ besteht aus zwei rationalen Funktionen (Pole 0, a₀, O(ζ^{−N−2}), C 13.0); jede darf ihr EIGENES Gebiet benutzen:
- Rest (I−Q_N)y: Bernstein–Walsh in Ĉ∖K' mit (a′)_r ⇒ (E_R') (Satz P-II, ganz ∂K' liegt ausserhalb f_K(D)) ⇒ |(I−Q_N)y_ζ| ≤ CN^{a+2}e^{−(N+2)g'(ζ)}.
  Keine Eckbedingung nötig (volle Green-Rate).
- Eckteil Q_N y: Poisson–Jensen in Ĉ∖W̄ mit Randdaten auf ∂W̄ aus O6/D6 + **(V1)–(V3)** (W §10): (V1) Profil auf dem D-Bogen, λ ≤ Aτ^p, p > 1/α = 1.1326,
  λ < 0.065 fern; (V2) Aussenwerte ausserhalb D̄ nahe Γin (Linse) polynomial; (V3) **Werte der adjungierten Moden in der Bucht h(A-Bogen) ⊂ D
  polynomial** (σ|V| ≤ CN^κ; numerisch W: Sockel). Dazu (E_sing) und Kantenschwanz für A_f. ⇒ |Q_N y_ζ| ≤ CN^κ e^{−ε(N+2)g_W(ζ)}.
- Zusammen auf Schlitz ∪ Γ'_d: **(H₀^mix)** |E_N(ζ,z)| ≤ CN^κ(e^{−Ng'(ζ)} + e^{−εNg_W(ζ)})/dist(ζ,K'). Satz P-I(i) mit dem langsameren Teil ⇒
  sup_{D̄}|Π_N f_K − f_K| = O(N^{κ−4.2877}); Schwelle κ < 4.29 (Rest: κ ≈ a + 2 ≈ 2 bei σ_{r+1} ≈ 1).
(V3) ist eine Aussage nur über die r(N) Eckmoden (endlicher Rang, Grenzschicht-Theorie W), nicht über die volle Pol-Antwort — deutlich
schwächer als (E_A). Der K*-Weg unten bleibt als Alternative stehen.

"""
assert anker in s
s = s.replace(anker, neu + "**Alternative (K*-Weg, P 00:42) — Zielaussagen:**")
s = s.replace("(E_R*) + (E_K*) ⇒ (H₀*) [P-I] · (H₀*) ⇒ (H*) [Familien-Stabilität, offen] ·",
              "(E_R*) + (E_K*) ⇒ (H₀*) [P-I] · BEVORZUGT: (a′)_r ⇒ Rest in K' [P-II] + (V1)–(V3) + (E_sing) + Kantenschwanz ⇒ Eckteil in W̄ [P-IV, C 13.2] ⇒ (H₀^mix) ⇒ (LP) mit N^{κ−4.29} [P-I] ·\n(H₀*) bzw. (H₀^mix) ⇒ (H) [Familien-Stabilität, offen] ·")
p.write_text(s, encoding="utf-8")
print("ok")
