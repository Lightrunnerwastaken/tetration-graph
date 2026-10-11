"""Hilfsskript (Strang P): ARCHITEKTUR §9.8 durch das Ergebnis von Satz P-V ersetzen (Runde 4). Aufruf: python arch_s98.py HH:MM"""
import sys
from pathlib import Path

ZEIT = sys.argv[1]
p = Path(__file__).resolve().parent.parent / "ARCHITEKTUR.md"
s = p.read_text(encoding="utf-8")
i = s.index("### 9.8 Zertifizierte Kneser-Werte")
neu = f"""### 9.8 Zertifizierte Kneser-Werte — Satz P-V BEWIESEN (rigoros, arb; P, Stand 11.10. {ZEIT})

**Ergebnis.** Knesers sexp zur Basis e ist auf dem Streifen −½ ≤ Re z ≤ ½ (alle Höhen) mit Fehler **ε ≤ 1.3e-168** durch explizite Daten
(p26_iter6.json) eingeschlossen: sexp_K = T̃ ± ε für |Im z| ≤ 0.45 (T̃ reelles Polynom vom Grad 519 um 0) und sexp_K = Ũ ± ε für Im z ≥ 0.35
(Ũ(z) = P(L^{{z+θ̃(z)}}), P die Poincaré-Funktion an L, θ̃ trigonometrisches Polynom mit 241 Moden), unten gespiegelt, sonst über exp/Log.
Daraus **zertifiziert: Knesers f_K-Koeffizienten a_0 … a_640 mit |Fehler|·ρⁿ ≤ 2.9e-127** (≥ 126 gewichtete Stellen; a_0 auf 8.8e-143,
a_600 auf 3.8e-211 absolut) und |a_n| ≤ 2.128·1.30^{{−n}} für alle n — Datei **gemeinsam/kneser_halbexp_koeff_zert_640.json** (Os Anfrage).
Gegenprobe: Strang As 50-stellige Datei stimmt damit auf 1.1e-43 (gewichtet) überein.

**Beweis (Fassung 2 — quasikonforme Glättung; das Cousin-Netz der Fassung 1 hatte am Tripelpunkt einen log-Verlust und ist verworfen):**
1. *Daten (Numerik, p26):* Cousin-Newton auf dem Zylinder (Schleife Im z = 0.25, Naht Re z = −½, Knotenkonstante per Rampe), Start aus
   dem Motor (dps 50; komplexe Werte des Forks versagen ab dps 60), sechs Schritte, Residuen 2e-31 → 2e-59 → 4e-95 → 5e-130 → 5e-169.
2. *P zertifiziert* (kc_rig, p25): p_n per Rekursion, |P − T_{{n0}}| ≤ C|w|^{{n0+1}} auf |w| ≤ ½ (Fixpunktargument); Ũ ist ganz mit Ũ(z+1) = e^{{Ũ(z)}} EXAKT.
3. *Defekte (p27):* e mit Ũ(ζ+e) = T̃(ζ) auf G_B = [−0.55, 0.55]×[0.30, 0.50], d mit T̃(ζ+1+d) = e^{{T̃(ζ)}} auf G_S = [−0.65, −0.45]×[−0.55, 0.55]:
   Maximumprinzip + Taylor-Modelle (Ordnung 240, Cauchy-Schwanz, 900 Bit) ⇒ sup|Ũ − T̃| ≤ 2.7e-173, sup|T̃(·+1) − e^{{T̃}}| ≤ 2.1e-173; |Ũ'| ≥ 0.83,
   |T̃'| ≥ 1.14 (Taylor-Modelle auf Kästchen) ⇒ **|e| ≤ 6.5e-173, |d| ≤ 3.7e-173** (skalares Kantorovich, eindeutig im Radius 0.02).
4. *Exakte Lipschitz-Lösung F₀ (p28):* W₀ = χ(Re z)·d(z−1) (χ: 0 → 1 auf [0.4, 0.467]), Y = W₀ + e(·+W₀), V₀ = ψ(Im z)·Y (ψ: 1 → 0 auf
   [0.40, 0.45]); F₀ = T̃∘(id+W₀) für |Im z| ≤ 0.40, F₀ = Ũ∘(id+V₀) für Im z ≥ 0.35 (dort gleich). **F₀(z+1) = e^{{F₀(z)}} gilt exakt**: an der Naht
   nach Definition von d, oben weil Y periodisch ist — wegen der exakten Identität e(ζ) − e(ζ+1+d(ζ)) = d(ζ) (folgt aus der Exaktheit von Ũ und
   der Eindeutigkeit von e). Damit gibt es KEINE Tripelpunkt-Verträglichkeitsbedingung mehr (die Fassung-1-Rampe wird nicht gebraucht).
5. *Beltrami (p28):* |μ₀| = |∂̄F₀/∂F₀| ≤ 1.3e-171, Träger in |Im z| ≤ 0.45. Messbarer Riemannscher Abbildungssatz auf dem Zylinder (q = e^{{2πiz}},
   Hauptlösung ψ₀ = q + C[h], h = (I − μS)⁻¹μ in L⁴, ‖S‖₄ ≤ 12) ⇒ qk φ mit φ(z+1) = φ(z)+1, **|φ − id| ≤ 1.3e-169**; Stoilow: F = F₀∘φ⁻¹ ist
   analytisch, F(z+1) = e^{{F(z)}}, reell auf ℝ, F → L (Im → +∞), L̄ (Im → −∞).
6. *Identifikation (p28):* Prüfung A (F ∉ (−∞, 0] auf dem Streifen) und B (F ∉ (0, 1] für 0 < Re ζ < ½; nahe 0 analytisch über F(0)=1,
   F'(0) ≥ 1.0918, |F''| ≤ 1.71) auf 2958 Scheiben für 0 ≤ Im ≤ 3, Schwanz Im ≥ 3: |F − L| ≤ 0.065 < Im L ⇒ F holomorph auf Re z > −2 mit
   beiden Grenzwerten für jedes x > −2 ⇒ **Paulsen–Cowgill Prop. 2** (Fassung Paper VI Thm 11.11): F(· + τ) = sexp_K, |τ| ≤ 5.2e-169.
7. *Fehler (p29):* ε_T, ε_U ≤ 1.27e-168 (mit sup|T̃'| ≤ 1.95, sup|Ũ'| ≤ 1.69).
*Zitierte Sätze:* Ahlfors–Bers (messbarer RMT, Hauptlösung), Stoilow-Faktorisierung, Beurling-Transformation ‖S‖_p ≤ 4(p*−1)
(Bañuelos–Wang 1995; bei |μ| ~ 1e-171 genügt jede endliche Schranke), Paulsen–Cowgill, Adv. Comput. Math. 43 (2017), Prop. 2.
**f_K-Koeffizienten (p29):** f_K(x) = sexp(slog x + ½) auf |x| = 1.2 (4096 Punkte; Krawczyk 1D für slog auf D(w̃, r₀) mit sexp-Bällen, dann
sexp-Ball auf D + ½; Fehler je Wert ≲ 3e-167), DFT, Aliasfehler ≤ M r'^{{−n}}(1.2/1.3)^{{4096}} mit M = sup_{{|x|=1.30}}|f_K| ≤ 2.128 (adaptive
Bogen-Krawczyk, 1950 Scheiben).
**gp-tetration (nur gelesen):** dort gab es für b = e nur Vorzertifikate (Paper VI setzt den Anker T_e voraus); Satz P-V liefert genau
dieses Anker-Zertifikat (Def. 11.3) — mit 168 Stellen.
**Als Nächstes (P):** (i) f_K und f_K' auf ∂D inkl. Umgebung von L (fK_zert vorhanden) ⇒ Os Pfadbedingungen P1–P4 streng; (ii) Startbedingung (A):
Π_{{N₀}}f_K aus den zertifizierten Koeffizienten (T_N J_∞f_K mit Schwanzschranke) und sup_{{D̄}}|p_{{N₀}} − Π_{{N₀}}f_K| ≤ ‖·‖_{{ℓ¹_ρ}} gegen die Box aus p24;
(iii) o22 mit diesen Koeffizienten (Sache von O).
"""
s = s[:i] + neu
p.write_text(s, encoding="utf-8")
print("§9.8 ersetzt")
