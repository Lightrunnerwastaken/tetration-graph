### 9.3 Korrigierte Zielaussage (Stand 11.10. 00:45) — massgeblich für W, O, M, X

**Rahmen.** Für ein Kompaktum K ⊃ D̄ (0, a₀ ∈ int K, Ĉ∖K zusammenhängend und regulär) sei g_K die Green-Funktion von Ĉ∖K mit Pol ∞,
ω^K_ζ das harmonische Mass, α_K·π der Aussenwinkel bei L, S_L der radiale Schlitz L(1+u) (Richtung 76.62°) bis zur Niveaulinie Γ^K_d.
Drei Rahmen (Zahlen p01, p13, p13b):

| Rahmen | K | Aussenwinkel | α_K | Schwelle α_K(1+Re γ) | Eckbedingung an das Γin-Profil λ ≤ Aτ^p | cap | κ₀ |
|---|---|---|---|---|---|---|---|
| W̄ (alt) | D̄ ∪ D̄(a₀,R₀) | 158.934° | 0.882964 | 4.8094 | p > 1/α = 1.1326 | 1.5169 | 6.2732 |
| K' | D̄ ∪ f_K(D̄) | 141.691° = π − arg μ | 0.787172 | 4.2877 | p > 1/α' = 1.2704 | 1.6556 | 5.8170 |
| K* (Hybrid, empfohlen) | = W̄ in D(L,r₁) ∪ D(L̄,r₁), = K' ausserhalb D(L,2r₁) ∪ D(L̄,2r₁) | 158.934° | 0.882964 | 4.8094 | p > 1.1326 | — | — |

K': g'(L(1+u))/u^{1/α'} = 0.839 (u → 0) … 0.754 (u = 0.204), u₀'(0.05) = 0.114, u₀'(0.1) = 0.204, g' ≤ g_W geprüft (p13b).

**Warum der Hybrid.** In K' liegen alle Urbilder der Randpunkte auf oder ausserhalb ∂D ⇒ (a′) trägt per Norm (Satz P-II); aber die
Eckbedingung verlangt p > 1.2704, und das gemessene Γin-Profil hat p ≈ 1.26–1.29 (W 00:24; Vermutung W1: p = 1/α' genau) ⇒ K' ist
**grenzwertig**: sup_Schlitz ω'/g' wächst logarithmisch (0.33/0.53/0.73/0.93 bei u = 1e-4/-6/-8/-10, p17), also formal nicht „für alle N“
(Kreuzung bei N ~ 1e14). In W̄ konvergiert dasselbe Profil (0.22), aber der A-Bogen liegt in f_K(D) (Eingaben ~(ρ/|h|)^N, p02). K* nimmt
die Ecke von W̄ und den Rest von K'; ∂K* ∩ f_K(D) liegt dann in A_r := D(L,2r₁) ∪ D(L̄,2r₁) (A-Bogen-Stück nahe L plus Verbindungsbogen zu
f_K(Γout)).

**Zielaussagen (das ist zu beweisen; die Kette dahinter ist in 9.4 bewiesen):**
- **(a′)_r** (O 00:32): Es gibt r(N) und a < ∞ mit σ_{r(N)+1}((I + K_N)⁻¹) ≤ C·N^a (gemessen 1 + 2.2e-5 bei Kneser, O; r(N) ≈ N/250).
- **(E_L⁺)** (gemeinsames Eckprofil aller r(N) Moden, ZWEIDIMENSIONAL nahe L): Φ_N(p) := ‖Q_N Y_N e(p)‖ (e(p) = ρ-skalierte Koeffizienten
  von 1/(p − x), |p| ≥ ρ), λ_N(p) := (1/N)·log(Φ_N(p)/Φ_N(L)). Verlangt: λ_N ≤ A|τ|^p mit p > 1.1326 auf der Γin-Seite von ∂D nahe L;
  λ_N ≤ 0 auf der Γout-Seite nahe L (W: = 0⁻) **und für |p| > ρ nahe L** (Aussenwerte — sie treten als Urbild- und Polterme auf;
  kubisch p10: −0.004 … −0.014); λ_N ≤ λ_max < 0.065 sonst; Φ_N(L) ≤ CN^κ.
- **(E_sing)** (nur Kneser): der Eckanteil des singulären Teils von f_K bei L ist polynomial: ‖Q_N Y_N P_N[φ_sing]‖ ≤ CN^κ für die
  Singulärterme φ_sing ~ (x−L)^{1+γ}(…) des O6-Rests von 1/(m(η − f)). Gestützt: Eckanteil von Π_N f_K = 8e-9 (D d19, N = 800).
- **(E_A)** (die einzige echte Auslöschungsaussage): Für η ∈ ∂K* ∩ f_K(D) ist R_N(η) = |y_η|_{ℓ¹_ρ} ≤ C N^κ e^{Nλ_A(η)} mit λ_A ≤ 0 nahe L
  (obwohl dort |w_η| ~ e^{cN|η−L|}, c ≈ 0.18). Gestützt: R_N auf dem A-Bogen fällt (C c27: A0.02: 14.2 → 2.1 für N = 200 → 1000).
- **(Fam)** Familien-Stabilität (H₀*) → (H*) (wie bisher offen).

**Satz (Hauptsatz im K*-Rahmen; Kette bewiesen, 9.4):** (K) + Univalenz (num.) + b35b′ + (a′)_r + (E_L⁺) + (E_sing) + (E_A) + Bulk ⇒ (H₀*) ⇒
sup_{D̄}|Π_N f_K − f_K| = O(N^{κ−4.81}); mit (Fam) und (A) ⇒ p_N → f_K. Ohne Zerlegung genügt die messbare einteilige Form:
R_N(L) ≤ CN^κ und Θ*_N := sup_{S_L ∪ Γ*_d} ω*_ζ(λ_N)/g*(ζ) ≤ 1 − ε. Numerisch: Θ_N(W̄) ≤ 0.105 (C, N ≤ 1000), Θ'_N(K') ≤ 0.110 (p15/p15b,
konservativ, N ≤ 1000); K* liegt dazwischen.

### 9.4 Bewiesene Reduktionen (P, 11.10. 00:45; Prüfskripte p13–p17)

**Satz P-I (rahmenunabhängig, bewiesen).** Für jedes K wie oben mit (F_K) [f_K holomorph um K bis auf S_L, S_L̄, stetig bis an die Ufer,
Sprung ≤ C_J|ζ−L|^{5.447}] und g_K(L(1+u)) ≥ c_K u^{1/α_K} auf dem Schlitz gilt:
(i) (H₀^K) ⇒ sup_{D̄}|Π_N f_K − f_K| ≤ C(ε^{−5.447α_K}N^{κ−5.447α_K} + N^κ e^{−εNd}); ebenso für Π^P_N ⇒ H1/H2 mit dieser Schwelle (Beweis = §4.1).
(ii) Einteilig (= C 13.4 in Ĉ∖K): R_N(L) ≤ CN^κ und Θ^K_N ≤ 1 − ε ⇒ (H₀^K). (iii) Zerlegt (= C 13.1/13.2): (E_R^K) + (E_K^{∂K}) ⇒ (H₀^K).
Beweis von (ii)/(iii) ⇒ (H₀^K): E_N(ζ,z) = y_ζ(z) − [k_ζ]_{>N}(z), |y_ζ(z)| ≤ R_N(ζ) auf D̄, |[k_ζ]_{>N}(z)| ≤ e^{−(N+1)g_K(ζ)}/(|ζ|−ρ) da
g_K ≤ log(|ζ|/ρ) für K ⊃ D̄; y_ζ, Q_N y_ζ, (I−Q_N)y_ζ sind rational mit Polen 0, a₀ ∈ int K und O(ζ^{−N−2}) ⇒ Bernstein–Walsh bzw.
Poisson–Jensen in Ĉ∖K. ∎
(F_{K'}) und (F_{K*}) folgen aus L8.2–8.4 und f_K = exp∘(f_K|_V)⁻¹ auf f_K(V), V Umgebung von D̄ ohne kleine Scheiben um L, L̄, sofern f_K
dort univalent ist (num.: Windungszahl 1, min|f(x)−f(y)|/|x−y| = 0.3456 auf ∂D, p16; f_K = exp∘h auf f_K(Γout) vom Motor bestätigt, p04).
Die Ecken-Asymptotik g_K ≍ u^{1/α_K} ist Standard (Ecke mit Dini-glatten Schenkeln); für K' numerisch: Exponent 1.27018 gegen 1.27037.

**Satz P-II ((a′)_r ⇒ Rest-Schranke; bewiesen bis auf num. Univalenz).** Sei Γ ⊂ ∂K ein Randstück mit Γ ∩ f_K(D) = ∅ und |η| ≥ ρ auf Γ.
Dann |w_η|_{ℓ²_ρ} ≤ C N^{3/2} gleichmässig auf Γ, also unter (a′)_r: |(I−Q_N)y_η|_{ℓ¹_ρ} ≤ C N^{a+2} auf Γ.
Beweis: w_η = P_N[1/(η−f)] + P_N[m·k_η]; zweiter Teil ≤ ‖m‖_∞√(N+1)/ρ. Erster: Koeffizienten über |x| = r_N = ρ(1−1/N):
|c_k|ρ^k ≤ e/dist(η, f_K(D̄_{r_N})) ≤ eN/(cρ) — denn η ∉ f_K(D), f_K bi-Lipschitz auf D̄ (für η = f(x₀), x₀ ∈ ∂D: |η − f(x)| ≥ c_f|x₀ − x|
≥ c_f ρ/N; nahe L trennt der Winkel 38.3° zwischen Γin und f_K(D) zusätzlich). Dann (I−Q_N)y_η = 𝓡_N w_η − (I−Q_N)[k_η]_{≤N} (C 13.0) und
‖𝓡_N‖ ≤ σ_{r+1}((I+K_N)⁻¹)·max_{∂D}|1/m| (D8; max = 1.3663, p12). ∎ Numerisch |[1/(η−f)]_{≤N}| ≤ 12 auf ∂K' bis N = 400 (A-Scheitel 1.5e33), p16.
Für K' ist Γ = ∂K' (ganz); für K* ist Γ = ∂K* ∖ A_r; auf ∂K* ∩ A_r liefert (E_A) beide Teile (orthogonale Projektion: |Qy|, |(I−Q)y| ≤ |y|).

**Satz P-III (Eck-Kriterium, bewiesen).** Randdaten λ ≤ Aτ^p auf einer Seite der Ecke, ≤ 0 auf der anderen:
lim_{u→0} ω^K_{L(1+u)}(λ)/g_K = (1/κ₀^K)∫λ s^{−2}ds ist endlich genau dann, wenn p·α_K > 1 (Randkoordinate s ≍ τ^{1/α_K}); bei p·α_K = 1 wächst
das Verhältnis wie (A'/κ₀)·log(1/u) (p17).

**Satz P-IV (Lokalisierung für A_f = O6 + Rest; bewiesen).** Mit O6: Q_NΠ_N k_η = F_Q(η) + Σ_j r_j F_Q(p_j) + Q_N Y_N c_hol. Für Kneser ist c_hol
nicht geometrisch abfallend (f_K ist in L, L̄ auf |x| = ρ singulär); zerlegt man c_hol = c_reg + c_sing (c_reg holomorph auf D_R, R > ρ; c_sing =
Singulärterme bei L, L̄, Koeffizienten ~k^{−6.45}), so gilt ‖Q_NY_N c_reg‖ ≤ ‖Y_N P_{N'}‖C + σ₁C(ρ/R)^{N'} (O6), und ‖Q_NY_N c_sing‖ ist (E_sing).
Mit (E_L⁺) an η und an den Urbildern p_j (Γin-Punkte: p_j ausserhalb D̄ nahe L; f_K(Γout)-Punkte: η aussen, p_j ∈ Γout) ist die Randrate des
Eckteils auf ∂K* ∖ A_r höchstens max(Γin-Profil, 0) + O(log N/N), auf ∂K* ∩ A_r nach (E_A) ≤ 0; Satz P-III mit α = 0.883 und p ≈ 1.27 ⇒ (E_K^{∂*})
mit ε > 0, sobald der Fernteil λ_max·sup ω/g < 1 (W̄-Konstante 14.65 ⇒ λ_max < 0.068). ∎ (Zusatz aus O6: Bulk ‖Y_N P_{0.7N}‖ ≤ C, gemessen 1 + δ₀.)

**Offen und präzise benannt:** (a′)_r [O] · (E_L⁺) inkl. Aussenwerte nahe L [W] · (E_sing) [W/O] · (E_A) [neu: Auslöschung auf dem A-Bogen nahe L;
Kandidat Os Dualbild O2b — Kontrolle von I − Y_N auf H²(D(a₀,R₀))] · Bulk ‖Y_N P_{0.7N}‖ ≤ C [O] · (Fam) · (A) [P, Zertifikat in Arbeit].

===96===
### 9.6 Abhängigkeiten (Stand 00:45)

(K) ⇒ (F) [L8.2–8.4, b35b′ num.] ⇒ (F_{K*}) [+ Univalenz num.] · (F_{K*}) + (H₀*) ⇒ (LP), Rate N^{κ−4.81} [P-I] · (F_{K*}) + (H*) ⇒ H1* ⇒ (mit (A)) H2* ⇒
Vermutung · (a′)_r ⇒ (E_R*) auf ∂K*∖A_r [P-II] · (E_A) ⇒ beide Teile auf ∂K*∩A_r · (E_L⁺) + (E_sing) + Bulk ⇒ (E_K^{∂*}) [P-IV + P-III] ⇒ (E_K*) [C 13.2] ·
(E_R*) + (E_K*) ⇒ (H₀*) [P-I] · (H₀*) ⇒ (H*) [Familien-Stabilität, offen] · alternativ einteilig: R_N(L) polynomial + Θ*_N ≤ 1−ε ⇒ (H₀*) [P-I(ii)].
