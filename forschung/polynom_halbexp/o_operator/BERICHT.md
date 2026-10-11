# Strang O — Operatortheorie: (a′) Rest-Kontraktion von (I + K_N)⁻¹

Stand: 11.10.2026, 04:46 (Beginn 10.10. 23:29; Runde 2 ab §9, Runde 3 ab §13, Runde 4 ab §19, Runde 5 ab §23, Runde 6 ab §26, Runde 7 ab §32, Runde 8 ab §33). Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.
Bezeichnungen wie ARCHITEKTUR.md und D §11: H = H²(|x| < ρ) ≅ ℓ²_ρ, P_N = Polynome vom Grad ≤ N (orthogonale Projektion in H),
Kδ = (δ∘f)/m, m = f'∘f, K_N = P_N K P_N, **Y_N := (I + K_N)⁻¹ = J_N⁻¹M_N**, a₀ = f(0), R₀ = |L − a₀|, h = f⁻¹,
T⁽⁰⁾_N, T^{(a₀)}_N = Taylor-Polynom vom Grad N bei 0 bzw. a₀, S = Taylor-Verschiebung 0 → a₀ (Koeffizientenwechsel), U_N = ([x^k](f−a₀)^j).
Modelle (Code `o_kern.py`): `kubisch` (D §9.3, μ = √L), `moeb_ell` (Strang Es „Möbius-exp“: Drehung um L, L̄ mit g(0) = 1, c = 0.705104),
`moeb_sym` (Drehung um ±i, ρ = 1, Winkel c), `moeb_lox` (f = M⁻¹(μM), D d24), `affin` (f = L + μ(x − L)).

## 0. Kurzfassung

1. **Satz O1 (bewiesen, exakt; geprüft 1e-136):** J_N = M_N(I + K_N) = U_N(I + K'_N)S mit K'_N = T^{(a₀)}_N K⁻¹|_{P_N}, K⁻¹δ = f'·(δ∘h);
   (I + K'_N)⁻¹ = S(I − Y_N)S⁻¹; und die **Mischformel**  Π_N g = Y_N T⁽⁰⁾_N g + (I − Y_N) T^{(a₀)}_N g.
   Die J-Projektion mischt die beiden Taylor-Polynome; Y_N ≈ Projektor auf den Γin-Anteil. Dualbild: (a′) bei a₀ heisst
   „I − Y_N beschränkt in H²(D(a₀,R₀))“.
2. **Satz O2 (bewiesen): Herglotz-Kriterium.** Für Möbius-f mit σ_D∘f Spiegelung (z. B. elliptisch, Fixpunkte auf ∂D) ist K auf allen
   Polynomen akkretiv ⇔ Re b ≥ 0 auf einer expliziten Scheibe D_H (b explizit rational). Dann **‖Y_N‖ ≤ 1 für alle N** ((a′) mit r = 0,
   Konstante exakt 1). Symmetrisches Modell: ⇔ sin²(c/2) ≤ (4√2 − 5)/4, d. h. c ≤ c* = 0.834467. **Möbius-exp:** Re b ≥ 0.607 (arb),
   Dualbild Re b ≥ 1.204 (arb) ⇒ ‖Y_N‖ ≤ 1 und ‖I − Y_N‖_{H²(D(a₀,R₀))} ≤ 1 für alle N, ‖J_N⁻¹‖ ≤ max_{∂D}|1/m| = 1.3225 (scharf).
3. **Satz O5 (bewiesen; „Riemann–Hilbert-Test“ ohne Riemann–Hilbert):** Im Möbius-exp-Modell gilt UNBEDINGT
   sup_{D̄}|Π_N f_M − f_M| ≤ C N^κ·0.53495^N (Geometrie K' = D̄ ∪ f(D̄) nach Strang P; wahre Rate 0.40918). Symmetrisch: für c < π/4.
   Damit ist Es offene Voraussetzung (H*)_M für die lineare Aussage umgangen (nicht scharf, aber vollständig).
4. **Satz O3 (bewiesen, lokales Mellin-Modell am Randfixpunkt):** lokal ist K Mellin-Multiplikator μ^{s−1} auf Re s = −½.
   ‖(I + K)⁻¹‖_lok = sup|1 + μ^{−3/2+iσ}|⁻¹: elliptisch (μ = e^{ic}) ≤ 1 ⇔ c ≤ π/3 (= Es Schwelle!); μ = √L: 1 + 1.368e-6.
   **K ist an einem loxodromischen Fixpunkt (|μ| ≠ 1) in keiner mellin-diagonalen Norm akkretiv** — Ds Akkretivitäts-Befund ist
   strukturell, kein Normwahl-Problem. Für Kneser braucht (a′) Normargumente, keine Wertebereichsargumente.
5. **Numerisch, wichtig für die Formulierung von (a′) und (E):** (i) Der Rest ist N-unabhängig **1 + 5.304e-5** (N = 400…1000), nicht
   1 + o(1) — harmlos. (ii) **Die Zahl der Eckmoden wächst:** Paare setzen bei N ≈ 60, 450, 950 ein; σ(Y_N) kubisch bei N = 1000:
   4212, 949, 5.84, 2.89, 1.0030, 1.0004, dann 1.00005 ⇒ „höchstens drei Eckmoden“ ist für N ≳ 1000 falsch; r(N) wächst etwa wie N/250.
   (iii) Phasenraum: Y_N ≈ χ_{Γin} im Inneren, Mellin-Übergang bei L, dunkles Band k/N ∈ [1/|μ|, 1] bei L und Deckschicht k ≈ N über Γin.

6. **Runde 2/3 (§9–§18):** (a′)_r mit Restkonstante δ₀ = Niederfrequenz-Norm des unendlichen Operators (Kneser 2.2e-5); die Kette verträgt
   σ_{r+1} ≤ CN^a, a < 3.15. Lemma O6 (gemeinsame Eckformel, Ersatz für D6 bei A_f; braucht nur den Kantenschwanz). O5 nicht auf |μ| > 1 schiebbar.
   **Satz O7:** T = I + K: E²(Ω) → H²(D) Fredholm (Lemma O7a lokal bewiesen, Kleben zitiert), Index 0 per Homotopie von Möbius-exp, Kneser
   invertierbar ⇐ D4|(G2) + (S) + Pfad. (a″)_∞ ⇔ ‖G‖ < ∞ erklärt Xs 0.7328 = min_{Γin}|m|/(1+δ₀).

## 0a. Ausrichtung auf Strang Ps Korrektur (K' = D̄ ∪ f_K(D̄), BRETT 23:46–00:11; Auftrag der Leitung, 00:13)

- Die Zielaussage (a′) bleibt **im Bild bei 0** (H²(D_ρ)) genau die richtige Eingabe für Ps korrigierte Kette: für η ausserhalb K' liegen ζ und
  alle Urbilder ausserhalb D̄, also ist J_∞k_η ∈ H²(D_ρ) mit Norm ≤ C·dist(η, K')^{−1/2}, und P_N ist orthogonal ⇒ (a′) ⇒ (E_R') wie D7, ohne
  den A-Bogen. Satz O5 (§3) ist diese korrigierte Kette **vollständig durchgeführt** im Möbius-exp-Modell (dort r = 0): Kontraktion (O2) + K' +
  Bernstein–Walsh ⇒ unbedingte Konvergenz mit Rate e^{−g_{K'}(ζ_p)} = 0.535 (statt 0.409 in der W-Geometrie).
- Was die K'-Geometrie kostet, sieht man dort explizit: Rate g_{K'} ≤ g_W. Für Kneser heisst das Schwelle α'(1+Re γ) − 1 statt α(1+Re γ) − 1
  (P: 4.29 statt 4.81 an den Schlitzen) — weiterhin weit über der gemessenen polynomialen Pol-Antwort.
- Die W-Geometrie zurückzugewinnen wäre Sache des **Dualbilds** (Satz O1(b), O2b): in Möbius-Modellen ist W = K' ∩ K'' mit
  K'' = D̄(a₀,R₀) ∪ h(D̄(a₀,R₀)); das Dualbild kontrolliert Pole ausserhalb K'' — aber nur auf D(a₀,R₀) in z. Ein Zwei-Variablen-Argument,
  das beides verbindet, habe ich nicht (§6).
- Konsequenz für die Parametrix: kein neues (a′)', aber die **wachsende Eckmodenzahl r(N)** (§5.1) muss in (H₀')/(E_K') eingehen.

## 1. Exakte Struktur (Satz O1, bewiesen)

**Satz O1.** Sei f holomorph nahe D̄ ∪ {a₀}, m = f'∘f ohne Nullstelle, u = f − a₀, U_N = ([x^k]u^j)_{k,j≤N} (untere Dreiecksmatrix),
S = Koeffizientenwechsel eines Polynoms vom Grad ≤ N von der Basis x^k zur Basis (x − a₀)^k. Dann
(a) J_N = U_N S + M_N (D3) = M_N(I + K_N) = U_N(I + K'_N)S, wobei K'_N = U_N⁻¹M_N S⁻¹ = S K_N⁻¹ S⁻¹ der Abschnitt bei a₀ von
    K⁻¹δ = f'·(δ∘h) ist: K'_N η = Taylor-Koeffizienten bei a₀ von f'·(δ∘h), δ = S⁻¹η.
(b) (I + K'_N)⁻¹ = S(I − Y_N)S⁻¹ =: Y'_N.
(c) (Mischformel) Für g holomorph nahe 0 und a₀: **Π_N g = Y_N T⁽⁰⁾_N g + (I − Y_N) T^{(a₀)}_N g.**
(d) Für Pole: Π_N k_ζ − T^{(a₀)}_N k_ζ = Y_N D_ζ, Cs y_ζ = Π_N k_ζ − T⁽⁰⁾_N k_ζ = −(I − Y_N)D_ζ mit dem Polynom
    D_ζ(z) = T⁽⁰⁾_N k_ζ − T^{(a₀)}_N k_ζ = [((z−a₀)/(ζ−a₀))^{N+1} − (z/ζ)^{N+1}]/(ζ − z).
*Beweis.* (a) P_N[δ∘f] = U_N S δ (δ∘f = (δ∘τ_{a₀})∘u, u(0) = 0) und P_N[mδ] = M_N δ. Für K'_N: U_N⁻¹ bildet die 0-Koeffizienten von g auf die
a₀-Koeffizienten von g∘h ab (Jets: (g∘h)∘f = g, h(a₀) = 0); also U_N⁻¹M_N S⁻¹η = a₀-Jet von (mδ)∘h = (m∘h)(δ∘h) = f'·(δ∘h).
(b) Y'_N = S J_N⁻¹U_N und I − Y_N = J_N⁻¹(J_N − M_N) = J_N⁻¹U_N S. (c) Π_N g = J_N⁻¹(P_N[g∘f] + P_N[mg]) = J_N⁻¹U_N·(a₀-Jet von g) + Y_N T⁽⁰⁾_N g
und J_N⁻¹U_N = (I − Y_N)S⁻¹, S⁻¹(a₀-Jet) = T^{(a₀)}_N g. (d) Einsetzen. ∎
Probe (o10_mischformel.py, kubisch, N = 60, 1500 Bit): g = eˣ und g = 1/(2.3 − x): Abweichung 2e-144 bzw. 3e-136. Für Kneser gilt alles
wörtlich mit Cs A_f (J = U S + M mit der Potenzreihe f).
*Bemerkungen.* (i) Für g ∈ P_N gibt (c) Π_N g = g. (ii) Für D ⊄ D(a₀,R₀) ist ∂D ∖ D̄(a₀,R₀) = Γin (beide Kreise gehen durch L, L̄): T^{(a₀)} divergiert
im Halbmond neben Γin, wo Y_N ≈ I den T⁽⁰⁾-Anteil nimmt; auf Γout ⊂ D(a₀,R₀) ist Y_N ≈ 0 und T^{(a₀)} trägt. Das ist der ganze Mechanismus der
Zwei-Kreis-Menge W in einer Zeile. (iii) Strang X hat die a₀-Identität für g = f_K (X §2.2); (c) ist die allgemeine lineare Form mit Y_N.

## 2. Satz O2 — Herglotz-Kriterium für Akkretivität (bewiesen)

**Lemma (Akkretivität ⇒ Kontraktion, trivial).** Ist Re⟨Kδ, δ⟩_H ≥ 0 für alle Polynome δ, so ‖(I + K_N)⁻¹‖ ≤ 1 für jedes N.
(⟨K_Nδ,δ⟩ = ⟨Kδ,δ⟩ für δ ∈ P_N, da P_N orthogonal; also Re⟨(I+K_N)δ,δ⟩ ≥ ‖δ‖².) Allgemeiner: ist Re⟨Kδ,δ⟩ ≥ −ε‖δ‖² auf einem Unterraum
der Kodimension r, so σ_{r+1}((I+K_N)⁻¹) ≤ 1/(1−ε).

**Satz O2.** Sei Ω = D(0,r), φ(y) = (αy+β)/(γy+δ) Möbius, a holomorph auf einer Umgebung von 0, und
Δ₁(y,z̄) := r²(γy+δ) − (αy+β)z̄ = e^{iκ}Δ_H(y,z̄) mit hermitescher Form Δ_H (Signatur (1,1)), Δ_H(0,0) > 0. D_H := Zusammenhangskomponente von
{y : Δ_H(y,ȳ) > 0} mit 0 (eine Kreisscheibe/Halbebene), b(y) := r²e^{−iκ}(γy+δ)a(y). Ist b holomorph auf D_H mit Re b ≥ 0, so gilt
Re⟨M_aC_φ p, p⟩_{H²(Ω)} ≥ 0 für alle Polynome p, also ‖(I + P_N M_aC_φ P_N)⁻¹‖ ≤ 1 für alle N. Umgekehrt folgt aus der Akkretivität Re b ≥ 0
auf D_H (Carathéodory–Toeplitz).
*Beweis.* Mit dem Szegő-Kern s_z(x) = 1/(1 − x z̄/r²) ist ⟨M_aC_φ x^j, x^k⟩ + conj⟨…x^k, x^j⟩ = r^{2(j+k)}·[y^k z̄^j]𝔎(y,z) mit
𝔎(y,z) := a(y)s_z(φ(y)) + conj(a(z)s_y(φ(z))). Also: Form ≥ 0 auf Polynomen ⇔ Taylor-Matrix von 𝔎 in (0,0) positiv semidefinit. Wegen
1 − φ(y)z̄/r² = Δ₁(y,z̄)/(r²(γy+δ)) und der Hermitezität ist 𝔎 = [b(y) + conj b(z)]/Δ_H(y,z̄). Schreibe Δ_H = c₀ℓ(y)conj(ℓ(z))(1 − Y(y)conj Y(z))
mit c₀ > 0, affinem ℓ ≠ 0 auf D_H und einer Möbius-Abbildung Y: D_H → 𝔻. Dann 𝔎 = (1/c₀)(1/ℓ(y))conj(1/ℓ(z))·[B(Y) + conj B(Z)]/(1 − YZ̄),
B = b∘Y⁻¹. Re B ≥ 0 ⇒ Herglotz: [B(Y) + conj B(Z)]/(1 − YZ̄) = ∫2/((1 − Ye^{−iθ})(1 − Z̄e^{iθ}))dν(θ) ≥ 0; Schur-Produkt mit dem Rang-1-Faktor ⇒ 𝔎 ≥ 0
auf D_H × D_H ∋ (0,0) ⇒ Taylor-Matrix psd. Umkehrung: Carathéodory–Toeplitz (Taylor-Matrix psd in einem Punkt ⇔ Re B ≥ 0 auf 𝔻). ∎
*Wann ist Δ₁ hermitesch?* Δ₁(y,ȳ) = 0 ⇔ φ(y) = r²/ȳ, d. h. Fixpunkte von σ∘φ (σ = Spiegelung an ∂Ω). Ist φ elliptisch mit beiden Fixpunkten auf ∂Ω,
so ist σφσ = φ⁻¹, σ∘φ eine antikonforme Involution = Spiegelung an einem Kreis durch die Fixpunkte ⇒ Δ₁ hermitesch bis auf eine Phase.

**Korollar O2a (symmetrisch, bewiesen).** f = (x+τ)/(1−τx), τ = tan(c/2), r = 1: Δ_H ∝ 1 − YZ̄ mit Y = (y+τ)/√(1+τ²), und
B(Y) ∝ (1 − 2tY)²/(1 − tY), t = sin(c/2). min_{|Y|=1} Re B = 1 − (5 + 4t²)²/32 (falls (5+4t²)/(16t) ≤ 1) ⇒ **akkretiv ⇔ t² ≤ (4√2−5)/4 ⇔
c ≤ c* = 0.834467**. Probe (o08, mpmath 120 Stellen, kleinster Eigenwert von Re K_N, N = 10/20/30/40): c = 0.80: +0.013/8e-5/4e-7/2e-9;
0.83: +0.011/6e-5/3e-7/1e-9; 0.84: +0.009/−0.49/−204/−4.5e4; 0.86: negativ ab N = 10. ✓
**Korollar O2b (Möbius-exp, bewiesen mit arb-Einschluss, o09c).** Bild 0: Δ₁ exakt hermitesch (ρ²γ = −β), D_H = D(−0.569713, 1.605138) (Rand durch
L, L̄), b = ρ²(γ(αy+β) + δ(γy+δ))²/(det(γy+δ)), einziger Pol 3.9527 ∉ D̄_H; **Re b ≥ 0.6074 auf ∂D_H** (4096 Bögen, Ball-Arithmetik) ⇒ Minimumprinzip
⇒ K akkretiv ⇒ ‖Y_N‖ ≤ 1 ∀N. Dualbild (K⁻¹ = M_{f'}C_h auf H²(D(a₀,R₀)), verschoben): D_H' = D(0.521994, 1.501048), Pol von b' bei 3.454 ∉ D̄_H',
**Re b' ≥ 1.2038** ⇒ ‖I − Y_N‖_{H²(D(a₀,R₀))} ≤ 1 ∀N. Mit D8: ‖J_N⁻¹‖ ≤ max_{∂D}|1/m| = 1.32254 — genau Ds gemessener Wert (d21), also scharf.
Numerisch zur Kontrolle (o01b, o06): σ₁(Y_N) − 1 = −1e-13 (N = 50), ≤ 4e-16 (100, 200); Dualbild ≤ 7e-16 (200).

**Numerik jenseits des Kriteriums (o03, symmetrisch, ‖Y_N‖ − 1 für N = 2…400):** c ≤ 0.85: ≤ 1e-15 bis N = 400 (also Kontraktion auch etwas
jenseits von c*; Akkretivität ist hinreichend, nicht notwendig). c = 0.90/0.95/1.0: kleiner, langsam wachsender Überschuss (2e-7/3e-5/4e-4 bei
N = 128). c = 1.04 (< π/3): O(1–8). c > π/3: exponentiell (Versagen des Verfahrens).

## 3. Satz O5 — unbedingte lineare Konvergenz im Möbius-exp-Modell (bewiesen)

**Satz O5.** Möbius-exp (Strang E §5, R3.5): f_M = Drehung um L, L̄ mit f_M∘f_M(0) = 1, Pol ζ_p = 3.95268. Dann gibt es C, κ mit
sup_{|z|≤ρ}|Π_N f_M(z) − f_M(z)| ≤ C N^κ·e^{−N g_{K'}(ζ_p)},  e^{−g_{K'}(ζ_p)} = 0.534953,
wobei K' = D̄ ∪ f_M(D̄) und g_{K'} die Green-Funktion von Ĉ ∖ K' mit Pol ∞ ist. (Wahre Rate nach E: e^{−g_W(ζ_p)} = 0.409185.)
*Beweis.* (1) ‖J_N⁻¹‖_{H→H} ≤ ‖Y_N‖‖M_N⁻¹‖ ≤ 1·sup_D|1/m| = 1.32254 (Kor. O2b, D8; m = f'∘f hat Pole nur bei h(ζ_p) = 1.8894, |·| > ρ).
(2) Für ζ ∉ K' liegen ζ und h(ζ) ausserhalb D̄, also J_∞k_ζ = 1/(ζ − f(x)) + m(x)/(ζ − x) ∈ H mit ‖J_∞k_ζ‖_H ≤ C·dist(ζ,K')^{−1/2}; P_N ist orthogonal,
daher ‖Π_N k_ζ‖_H ≤ 1.3226·C·dist^{−1/2} und sup_{D̄}|Π_N k_ζ| ≤ √(N+1)·‖Π_N k_ζ‖_H. Also |E_N(ζ,z)| ≤ C'√N·dist(ζ,K')^{−1/2} + 1/dist(ζ,D̄).
(3) E_N(·,z) ist rational mit Polen nur in 0, a₀ ∈ K' und O(ζ^{−N−2}) (B L7.1); Maximumprinzip auf {g_{K'} > 1/N} (dort dist ≥ cN^{−1/α'}) wie
B L7.2 ⇒ |E_N(ζ_p,z)| ≤ C N^κ e^{−N g_{K'}(ζ_p)}. (4) Es M1: Π_N f_M − f_M = r·E_N(ζ_p,·). ∎
Geometrie (o11): in M = (x−L)/(x−L̄) ist D eine Halbebene, f_M die Drehung um c = 0.70510, K' der Sektor der Öffnung π + c, W der Sektor π + c/2
(Probe: g_W gibt Es 0.409185 exakt). Symmetrisches Modell: Bedingung ζ_p ∉ K' ⇔ |cot c| > 1 ⇔ c < π/4; Rate sin(πc/2ω')/sin(π(π−3c)/2ω'),
ω' = π − c (c = 0.6: 0.537 statt 0.3705).
*Einordnung.* Es Satz M-lin brauchte (H*)_M auf der W-Niveaulinie (nicht bewiesen). O5 gibt dieselbe Aussage mit schwächerer Rate ohne jede
Annahme: Kontraktion (O2) + Ps Kompaktum K' (alle Urbilder ausserhalb D̄) + Bernstein–Walsh. Für die scharfe W-Rate fehlt die Kontrolle für
ζ ∈ K' ∖ W (A-Bogen-Seite); das Dualbild (O2b) liefert dort nur Kontrolle auf D(a₀,R₀), nicht auf D̄. In der Sektorsprache: W = K' ∩ K'',
K'' = D̄(a₀,R₀) ∪ h(D̄(a₀,R₀)).
*Nicht erfasst:* der nichtlineare Teil (Es M-voll braucht (H*_sek) für die Sekanten-Familie — dort ist φ = p kein Möbius mehr).

## 4. Satz O3 — lokales Mellin-Modell am Randfixpunkt (bewiesen) und warum Akkretivität bei Kneser scheitert

**Satz O3.** H_φ = {ε : φ < arg ε < φ + π}, H²(H_φ) (Randwerte in L²). K_μδ(ε) := δ(με)/μ (Linearisierung von Kδ = δ∘f/m mit f(L+ε) ≈ L + με, m ≈ μ).
Mellin: δ = ∫δ̂(σ)ε^{−1/2+iσ}dσ ⇒ ‖δ‖² = ∫|δ̂|²W(σ)dσ, W(σ) = 2π(e^{−2σφ} + e^{−2σ(φ+π)}), und K_μ ist Multiplikation mit μ^{s−1}, s = −½ + iσ
(gleiche Gewichte links und rechts). Folglich
(i) ‖(I + K_μ)⁻¹‖ = sup_σ |1 + μ^{−3/2+iσ}|⁻¹, (I+K_μ)⁻¹ bildet den Übergang von 0 (σ → −∞, Γout-Seite) nach 1 (σ → +∞, Γin-Seite);
(ii) K_μ akkretiv ⇔ Re μ^{−3/2+iσ} ≥ 0 ∀σ ⇔ |μ| = 1 und cos(3/2·arg μ) ≥ 0.
Für |μ| ≠ 1 durchläuft arg μ^{−3/2+iσ} = σ ln|μ| − 3/2·arg μ ganz ℝ, während |μ^{s−1}| = |μ|^{−3/2}e^{−σ arg μ} nach ∞ geht ⇒ Re K_μ ist nach unten
unbeschränkt. Gewichte |ε|^{2a} verschieben nur die Gerade (dasselbe Argument). **Folgerung:** An Knesers Fixpunkt (μ = √L, |μ| = 1.1724,
arg μ = 0.6686) gibt es keine mellin-diagonale Norm, in der K akkretiv ist — Ds negative Eigenwerte (−2.8e9) sind strukturell. Werte:
μ = √L: sup = **1 + 1.368e-6** (bei σ = 17.65); elliptisch μ = e^{ic}: Strahl in Richtung −1.5c, Norm ≤ 1 ⇔ **c ≤ π/3** — genau Es Schwelle
„f∘f∘f(0) = ∞“ (beide Bedingungen sind „drei halbe Drehungen = Viertelkreis“; ob das mehr als eine Koinzidenz ist: offen).
**Affines Modell** (f = L + μ(x−L), P_N invariant, kein Abschneide-Effekt, o01): ‖Y_N‖ = 1.2e15 (N = 50), 2e30 (100) — die globale Bedingung
„f(Γin) ⊂ D und h(Γout) ⊂ D“ ist dort verletzt; (I+K)⁻¹ ist unbeschränkt. Lokale Symbole allein reichen also nicht; man braucht die globale
Struktur (S): f(Γin) ⊂ D, h(Γout) ⊂ D (gilt kubisch/Kneser).

## 5. Numerik zum kubischen Modell (Kneser auf 1–3 % gleich, D)

**5.1 Spektrum von Y_N (o01, o07, o07b; volle Matrix in arb, 0.85N + 200 Bit):**

| N | σ₁ | σ₂ | σ₃ | σ₄ | σ₅ | σ₆ | nächster | #σ > 0.5 |
|---|---|---|---|---|---|---|---|---|
| 100 | 1.0484 | 1.0017 | 1.000053 | 1.000028 | | | | 0.554(N+1) |
| 200 | 1.492 | 1.0999 | 1.000053 | 1.000028 | | | | |
| 400 | 6.913 | 2.848 | 1.00017 | 1.00005 | | | | |
| 600 | 52.61 | 15.90 | 1.0432 | 1.0065 | 1.00005 | | | |
| 700 | 153.7 | 42.28 | 1.2099 | 1.0500 | 1.00005 | | | |
| 800 | 457.6 | 116.6 | 1.694 | 1.223 | 1.00005 | | | |
| 1000 | 4212 | 948.6 | 5.843 | 2.887 | 1.00301 | 1.00038 | 1.00005 | |
| 1200 | 40044 | 8177 | 30.65 | 12.35 | 1.0955 | 1.0226 | 1.00005 | |

Befunde: (a) Der Wert nach den Eckmoden ist **N-unabhängig 1 + 5.304e-5** (bis N = 1200) (dann 1 + 2.81e-5, …) — eine Eigenschaft des unendlichen
Operators (I + K)⁻¹, grösser als die lokale Mellin-Norm 1 + 1.37e-6; (a′) gilt also als „σ_{r+1} ≤ 1 + 5.3e-5 + o(1)“, nicht als 1 + o(1)
(für den Beweis unerheblich: polynomiale Schranken genügen). (b) **Eckmoden kommen paarweise und die Paare setzen nacheinander ein**:
Paar 1 bei N ≈ 60, Paar 2 bei ≈ 450, Paar 3 bei ≈ 950 (σ − 1 ≈ 1e-3). r(N) = 2 (N ≤ 400), 4 (600–800), 6 (1000–1200).
   Raten 1000→1200: Paar 1 0.0113, Paar 2 0.0083 (steigend), Paar 3 im Anlauf; Eingangstiefe (Schwerpunkt N − k) bei N = 1200:
   56–64 / 77–90 / 174–194, also ∝ N (0.05N, 0.07N, 0.15N). Ds „r ≤ 3“ beruhte auf
den letzten 100 Spalten. Folge: (a′) und Cs (E_K) müssen mit wachsendem r(N) formuliert werden (Q_N = Projektor auf r(N) Eck-Ausgänge),
oder man zeigt, dass spätere Paare von kohärenten Eingaben noch schwächer angeregt werden. (c) #σ > 0.5 ≈ |Γin|/2π = 0.574: Y_N ≈ Projektor
auf Γin-Inhalt.

**5.2 Phasenraum (o04, N = 300; Wellenpakete g bei (θ, k), Breite √N in k):** ⟨Yg,g⟩ | ‖Yg‖:
k/N ≤ 0.7: Γout (θ < arg L) 0 | 0, Γin 1.00 | 1.00, Übergang bei L in einer Schicht ~0.05 rad (L−0.05: 0.1–0.3, L+0.05: 0.7–0.9) = Mellin-Bereich;
k/N = 0.85: L+0.2: 0.55 | 1.00; k/N = 0.95: L+0.2: 0.01 | 0.26, L+0.6: 0.31 | 0.93, θ = π: 0.65 | 1.10; k/N = 1: θ = π: 0.30 | 0.94.
⇒ ein **dunkles Band** (Y ≈ 0) bei L für k/N ∈ [≈ 1/|μ| = 0.853, 1] (Bs Streckband) und eine **Deckschicht** k ≈ N über ganz Γin, in der Y
nicht ≈ I ist (dort leben die Eingangsvektoren v_i der Eckmoden, D §9.2). Y − T_N(χ) − Ecke hat Singulärwerte ≈ 1 nur mit Vektoren, die zu
99–100 % auf k > N − 30 sitzen (Ds 90 %-Befund präzisiert): die Toeplitz-Parametrix stimmt bis auf die Deckschicht.

**5.3 Dualbild (o06):** Singulärwerte von Y'_N = S(I − Y_N)S⁻¹ in H²(D(a₀,R₀)), kubisch:
N = 300: 12.9, 4.05 | 1.0066, 1.0063; N = 400: 85.1, 21.1, 1.163, 1.028 | 1.00663, 1.00660; N = 500: 620.7, 133.5, 2.41, 1.42 | 1.00688, 1.00688.
Möbius-exp ≤ 1 + 7e-16. Gleiches Bild wie bei 0 (Eckmoden + beschränkter Rest), aber die Eckmoden wachsen im Dualbild schneller (≈ 0.020/Grad)
und das zweite Paar setzt früher ein (N ≈ 350); Rest-Überschuss ≈ 6.6–6.9e-3 (ein Paar, langsam steigend — beobachten).

**5.4 Struktur der Paare (o12, N = 700):** Paar 1 (σ = 154, 42): Ausgang u bei arg L ± 0.1°, radial (x−L)^{≈−3.3…−3.96} (Konkurrenz-Exponent
Re(1−γ) = −3.447), Eingang v zu 92–99 % auf k > N − 60, Schwerpunkt N − k = 34–43, V auf Γin (≤ 0.08 auf Γout). Paar 2 (σ = 1.21, 1.05):
Ausgang bei ±79.8° (3° innerhalb Γin neben L), radial steiler (−3.6…−4.4), Eingang tiefer (Schwerpunkt N − k = 92–111, nur 22–35 % in den
letzten 60), V auch auf Γout merklich (0.7–0.9). Lesart (Vermutung): eine Familie von WKB-Eckmoden mit wachsender Knotenzahl in der Tiefe
(N − k)/N; ihre Zahl wächst linear in N.

## 6. Sackgassen (mit Grund)

- Akkretivität für Kneser/kubisch in jeder mellin-diagonalen oder Diagonal-Gewichtsnorm: unmöglich (Satz O3(ii)); Multiplikatoren B mit
  Re⟨B(I+K)δ,δ⟩ ≥ c‖δ‖² ebenfalls (lokal dreht μ^{iσ} jede Konstante durch alle Richtungen).
- Störung um ein elliptisches Modell: |μ| = 1.17 ist keine kleine Störung von |μ| = 1 (Spirale statt Strahl).
- Naive Parametrix P_N G P_N (G = (I+K)⁻¹): Fehler −P_N K Q_N G P_N ist exponentiell gross (Schwanz von Gg wird über Γout verstärkt).
- Transmissionsnorm H²(D ∪ f(D)) im Urbild: inf-sup scheitert (Γout-Pakete hoher Frequenz haben riesige Norm, werden aber fast identisch abgebildet).
- Kontrolle für ζ auf dem A-Bogen mit z ∈ D̄ aus O1 allein: jedes der beiden Bilder hat dort einen exponentiell grossen Summanden (Mischformel),
  die Auslöschung ist mikrolokal (Y ≈ I auf dem Halbmond-Anteil) — normweise nicht zu haben (deckt sich mit Ps Lücke 23:46/00:04).

## 7. Nächste Schritte

1. (a′) für das kubische Modell/Kneser: Normargument mit lokalem Prinzip; konkret Ms radiales Modell (K = Op(s(φ)^k)/m, Kohn–Nirenberg) als
   Testbett: Bulk + Mellin-Teil beweisen (dort ist alles explizit), Deckschicht/Ecke als endlich-dimensionales bzw. WKB-Problem abtrennen (mit W, M).
2. Wachstum r(N) klären: N = 1400–1600 (Paar 4 erwartet bei ≈ 1500), asymptotische Raten der Paare, Anregung der späteren Paare durch
   Pole nahe L — entscheidet, ob (E_K') mit wachsendem r(N) trägt (Abstimmung mit W/M).
3. O5 schärfen (W-Rate): Zwei-Variablen-Maximumprinzip für E_N(ζ,z) mit beiden Kontraktionen (O2b); nichtlinear: Herglotz für Sekanten-Operatoren.
4. Herglotz-Kriterium für nicht-Möbius-f: Positivität von 𝔎 = a(y)s_z(f(y)) + conj als Kernkriterium (kubisch: scheitert lokal nach O3; aber
   „Akkretivität modulo Ecke“ in einer nicht-diagonalen Norm prüfen).

## 8. Dateien (alle in o_operator/)

o_kern.py (Modelle, J̃, M̃, Y_N, Dualbild Y'_N) · o01_spektrum.py, o01b.py, o01c.py (Spektren) · o02_moebius_exakt.py (mpmath, A*A − I) ·
o03_scan_c.py (‖Y_N‖ über c, N) · o04_phasenraum.py · o05_akkretiv.py · o06_dual.py · o07_eckzahl.py (+o07_kubisch.log/json), o07b_eckzahl_vek.py
(+Vektoren npz) · o08_herglotz_test.py · o09_herglotz_moebexp.py, o09b_herglotz_dual.py, o09c_arb_dual.py (arb-Einschluss) · o10_mischformel.py ·
o11_moebexp_LP.py (Geometrie/Raten zu O5).


---

# Runde 2 (ab 00:17, Auftrag der Leitung): (a′)_r im radialen Modell, Anregung späterer Paare, Störung |μ| = 1 → 1.17

## 9. Korrigierte Formulierung (a′)_r und was die Kette verträgt

**9.1 Die Restkonstante ist eine Niederfrequenz-Eigenschaft des unendlichen Operators (numerisch, o16, o17).**
‖P_n Y_N P_n‖ (Block der ersten n+1 Koeffizienten) erreicht σ_{r+1}(Y_N) schon für n = 50 auf 7 Stellen; der zugehörige Singulärvektor
hat k-Schwerpunkt 9–13, 88–98 % Masse auf k < 20 und sitzt auf dem Kreis bei arg ≈ ±125…145° (tief in Γin, NICHT an L).
δ₀ := lim σ_{r+1} − 1 = ‖(I+K)⁻¹‖ − 1:

| Modell | δ₀ | ‖P_nY_NP_n‖ − 1 für n = 20/30/40/50 (N = 140) | lokale Mellin-Norm − 1 |
|---|---|---|---|
| Kneser (A_f aus gemeinsam/kneser_halbexp_koeff_1000.json, m = eˣ/f', kein Motor) | **2.2e-5** (σ₃(Y_140) − 1 = 2.17e-5) | 3.58e-6 / 1.32e-5 / 1.80e-5 / 2.00e-5 | 1.37e-6 |
| kubisch (D) | 5.30e-5 | 3.44e-5 / 4.79e-5 / 5.14e-5 / 5.24e-5 | 1.37e-6 |
| radial (M, α = 0.31382/ρ) | 1.046e-3 | 9.92e-4 / 1.039e-3 / 1.045e-3 / 1.046e-3 | 2.10e-5 |
| Möbius-exp | ≤ 0 (bewiesen, O2) | — | 0 |

**9.2 Formulierung (für P, ARCHITEKTUR §9.3).** (a′)_r: Es gibt δ₀ < ∞ und eine Folge r(N) mit σ_{r(N)+1}((I+K_N)⁻¹) ≤ 1 + δ₀ + o(1).
Vorschlag r(N) := #{i : σ_i(Y_N) > 1 + 2δ₀} = 2 × (Zahl der Eckmoden-Paare); gemessen r(N) ≈ N/250 (§5.1), ebenso im radialen Modell
(o07b_radial.log: N = 400: 5.43, 2.36 | 1.00105; N = 800: 285.8, 75.0, 1.449, 1.129 | 1.00105; N = 1000: 2374, 550, 4.17, 2.21, 1.00105, 1.00100 —
Paar 2 bei N ≈ 550, Paar 3 bei ≈ 1000).
**Was die Kette verträgt (Buchhaltung vorläufig — massgeblich ist Ps §9.4):** (E_R') braucht nur eine Randschranke auf {g' = 1/N} (danach
Bernstein–Walsh; die Rationalität bleibt bei der Zerlegung erhalten, C 13.0). Dort ist
‖(I−Q_N)Π_N k_ζ‖_{∞,D̄} ≤ √(N+1)·σ_{r+1}(Y_N)·‖1/m‖_∞·‖J_∞k_ζ‖_{H²} mit ‖J_∞k_ζ‖ ≲ dist(ζ,K')^{−1/2} ≲ N^{1/(2α')} (Ecke).
Mit σ_{r+1} ≤ C N^a also κ ≤ a + ½ + 1/(2α') = a + 1.135; Schwelle (4.1') κ < α'(1+Re γ) = 4.2877.
**⇒ σ_{r+1}(Y_N) ≤ C N^a mit a < 3.15 genügt sicher, C beliebig** (rechnet man den Faktor 1/dist aus (H₀') mit, sogar etwa a < 4.4). Die gemessene
Konstante 1 + δ₀ (a = 0) ist also kein Engpass; „1 + o(1)“ wird nicht gebraucht. Der Engpass ist allein der Eckteil mit wachsendem r(N).

## 10. Lemma O6 (gemeinsame Eckformel, bewiesen) und die Anregung der späteren Paare

**Lemma O6.** f mit J_N = U S + M (Polynom oder Cs A_f), Y_N = J_N⁻¹M_N, Q eine orthogonale Projektion in ℓ²_ρ. Für |p| > ρ sei e(p) der Vektor
der ρ-skalierten Taylor-Koeffizienten (k ≤ N) von 1/(p − x) und F_Q(p) := Q Y_N e(p). Dann:
(a) p ↦ p·F_Q(p) ist ein ℂ^{N+1}-wertiges Polynom in w = ρ/p vom Grad ≤ N; log‖F_Q‖ ist subharmonisch in w, also für |p| ≥ ρ
    log(|p|·‖F_Q(p)‖) ≤ ∫ log(ρ‖F_Q(ρe^{it})‖)·P_w(t)dt (Poisson-Kern in |w| < 1).
(b) Für η ∉ K' sei 1/(m(x)(η − f(x))) meromorph auf D_R (R > ρ) mit einfachen Polen p_j (Residuen −r_j) und holomorphem Rest c_hol
    (|c_hol,k| ≤ C_R(ρ/R)^k). Dann QΠ_N k_η = F_Q(η) + Σ_j r_j F_Q(p_j) + QY_N c_hol und, für jedes N' ≤ N,
    ‖QY_N c_hol‖ ≤ ‖Y_N P_{N'}‖·C_R/(1 − ρ/R) + σ₁(Y_N)·C_R(ρ/R)^{N'}/(1 − ρ/R).
*Beweis.* M_N⁻¹w_η = P_N[1/(m(η−f))] + P_N[1/(η−x)] (weil M_N⁻¹P_N[m g] = P_N g); Partialbrüche auf D_R; Linearität; (a) da
Y_N e(p) = Σ_k ρ^k p^{−k−1}Y_N e_k und log‖·‖ eines vektorwertigen holomorphen Polynoms subharmonisch ist. ∎
*Bedeutung.* Für polynomiales f ist das Ds D6 in koordinatenfreier Form, für alle r(N) Eckmoden gemeinsam (kein √r-Verlust); für Kneser
(A_f, unendlich viele Urbilder) ersetzt es D6 (Ps offener Punkt „Lokalisierung für A_f“, §9.2(2)/9.6): in D_R liegen nur endlich viele Pole
(η nahe L: das nahe Urbild h(η)), m = eˣ/f' hat keine Nullstellen, und der holomorphe Rest kostet nur die **Bulk-Bedingung**
(Bulk) ‖Y_N P_{0.7N}‖ ≤ C, zusammen mit ln(R/ρ) > c₁/0.7 = 0.016 (c₁ = 0.0113 Eckrate). Gemessen (o07b-Vektoren): der Eckanteil von Y_N P_{cN} ist
für c ≤ 0.7 nicht sichtbar (‖Y_N P_{0.7N}‖ = 1 + δ₀ bis N = 1200 kubisch, 1000 radial); bei c = 0.8: 1.31/2.57/4.40 (N = 700/1000/1200).
Die σ₁-gewichteten Eingangsvektoren fallen in die Tiefe: σ₁|v_{1,k}| = 0.036–0.038 bei k = N/2 (N = 700 und 1200 gleich), 0.18 → 0.35 bei
k = 0.75N — die adjungierte Mode hat einen WKB-Ausläufer ins Innere, ab Tiefe ≈ 0.3N ist er O(1).
⇒ Der Eckteil auf ∂K' ist durch das **gemeinsame Randprofil** ‖Q_N Y_N e(ρe^{iφ})‖ = (Σ_{i≤r}σ_i²|Ṽ_i(φ)|²)^{1/2}/ρ beschränkt — die richtige
Form von (E_L)' für wachsendes r(N) (Brücke zu W).

**10.1 Messung (o15, kubisch, Singulärvektoren aus o07b, N = 700/1000/1200; Eckanteil e_i(η) = σ_i|⟨v_i, M⁻¹w_η⟩|, Pol-Koeffizienten mit mpmath).**
Werte bei N = 1200, Raten 1000 → 1200 (Paar 3 war bei N = 1000 noch nicht gebildet; seine Raten sind Anlauf):

| Pol η | Paar 1 (σ-Rate 0.0113) | Paar 2 (σ-Rate 0.0083/0.0073) | Paar 3 |
|---|---|---|---|
| Ecke L | 4.49 (+0.0005) | 2.21 (+0.0007) | 0.84 |
| Γin τ = 2° / 5° / 10° | 8.7 (+0.0014) / 30 (+0.0032) / 171 (+0.0052) | 2.8 (+0.0008) / 3.4 (+0.0016) / 6.6 (+0.0033) | 1.9 / 1.8 / 1.8 |
| Γin τ = 20° / 40° / 90° | 1.5e3 (+0.0076) / 1.3e4 (+0.0098) / 6.3e4 (+0.0115) | 20 (+0.0059) / 48 (+0.0084) / 24 (+0.0085) | 1.7 / 0.4 / 0.5 |
| f(Γout) τ = −0.5° / −5° / −40° | 0.38 (−0.0019) / 7.7e-5 (−0.0023) / 1.8e-8 (−0.0023) | 0.23 (−0.0013) / 5e-5 (−0.0018) / 1.2e-8 (−0.0017) | 0.06 / 4e-6 / 3e-10 |
| Schlitz u = 0.001 / 0.01 / 0.1 | 3.05 (+0.0002) / 0.10 (−0.0020) / 3.8e-5 (−0.0023) | 1.54 (+0.0004) / 0.062 (−0.0016) / 2.4e-5 (−0.0019) | 0.59 / 0.017 / 2.4e-6 |

Zusatzrate λ_p(τ) = Rate − Eckrate (Γin, τ = 2/5/10/20/40/90°): Paar 1: 0.0009/0.0027/0.0047/0.0071/0.0093/0.0110, Paar 2: 0.0001/0.0009/0.0026/
0.0052/0.0077/0.0078 ⇒ **λ₂ ≤ λ₁ an jeder Stelle**; Eckanregung beider Paare polynomial (≈ N^{0.5…0.9}); auf dem neuen Randstück f(Γout) von ∂K'
und auf dem Schlitz (u ≥ 0.003–0.01) sind beide Paare klein UND fallend. Das Profil auf |x| = ρ (σ_i|Ṽ_i|, gleiche Datei) stimmt auf der Γin-Seite
mit den Eckanteilen überein und fällt auf der Γout-Seite (Paar 1 −0.0023, Paar 2 −0.0018…−0.0010) — zweiseitig, vorzeichenbehaftet, wie P verlangt.
**Vermutung (numerisch gestützt):** jedes spätere Paar p hat ein Zusatzprofil λ_p ≤ λ₁ und eine polynomiale Eckanregung; dann trägt (E_K') für das
ganze Q_N, sobald es für Paar 1 trägt.
*Achtung (W, BRETT „00:40“):* Werte von V_i im Inneren von D (A-Bogen-Urbilder) nicht mit float-Vektoren auswerten — hier nicht verwendet (alle Pole
auf ∂K' haben Urbilder auf/ausserhalb ∂D, die Koeffizienten von M⁻¹w_η sind O(1)).

## 11. Lässt sich O5 von |μ| = 1 auf |μ| = 1.1717 schieben? — Nein, nicht gleichmässig in N (numerisch + O3)

o14 (Möbius M⁻¹(μM), arg μ = arg√L; Einträge σ₁ / r(N) / Rest − 1):

| \|μ\| | N = 100 | N = 200 | N = 300 |
|---|---|---|---|
| 1.00 | 1 / 0 / 2e-16 | 1 / 0 / 7e-16 | 1 / 0 / 9e-16 |
| 1.01 | 1 / 0 / 2e-16 | 1 / 0 / 4e-5 | 1.004 / 0 / 4e-3 |
| 1.03 | 1.003 / 0 / 3e-3 | 1.18 / 1 / 2e-9 | 2.07 / 1 / 8e-9 |
| 1.06 | 1.16 / 1 / 7e-6 | 4.92 / 1 / 2e-5 | 54.4 / 2 / 2e-5 |
| 1.10 | 2.44 / 1 / 4e-4 | 111 / 2 / 6e-4 | 9990 / 2 / 5e-3 |
| 1.1717 | 27.2 / 1 / 7e-3 | 4.8e4 / 3 / 6e-3 | 1.3e8 / 4 / 6e-3 |

Befund: Für jedes |μ| > 1 entsteht eine Eckmode; ihr Einsatzpunkt N₁(|μ|) → ∞ für |μ| → 1 (≈ 300 bei 1.01, ≈ 120 bei 1.03), die Rate wächst mit |μ|
(1.1717: 0.075 wie D d24). Der Rest nach Abspalten bleibt 1 + (2e-9 … 7e-3). **Begründung (bewiesen, soweit O3):** für |μ| ≠ 1 ist das lokale Symbol
μ^{s−1} eine Spirale — K ist nicht akkretiv (O3(ii)), das Herglotz-Argument O2 bricht für jedes |μ| > 1 sofort; zugleich ist die Trefethen–Chapman-
Twistbedingung erfüllt (X §1.3), die exponentiell kleine Quasimoden an der Kante erzwingt (sup_N‖Y_N‖ = ∞ für jedes |μ| > 1 ist hier nur numerisch).
**Was sich schieben lässt:** nur der Rest σ_{r+1} und die lineare Kette mit abgespaltenen Eckmoden — also genau (a′)_r + (E_K'). Eine gleichmässige
Störungsrechnung um |μ| = 1 ist unmöglich (die Eckmoden sind eine singuläre Störung mit Rate c(μ) → 0 für |μ| → 1).

## 12. Stand Runde 2 und nächste Schritte

- (a′)_r im radialen Modell: **nicht bewiesen.** Gezeigt (numerisch): dieselbe Struktur wie kubisch/Kneser (Paare, Rest = Niederfrequenz-δ₀).
  Bewiesen: Lemma O6; die Restkonstante ist kein Engpass (a < 3.15 genügt). Offen: (i) ‖(I+K)⁻¹‖ < ∞ für den unendlichen Operator — Skizze:
  Transmissionsproblem auf ∂K' = Γin ∪ f(Γout) (auf Γin wertet K im Inneren aus, auf Γout wertet K⁻¹ = f'C_h im Inneren aus), also eine
  Randintegralgleichung zweiter Art auf ∂K' mit Mellin-Ecken bei L, L̄, deren lokales Symbol 1 + μ^{s−1} auf Re s = −½ nicht verschwindet
  (Fredholm); Injektivität = D4; (ii) die Deckschicht: nach Abspalten von r(N) Moden nichts über 1 + δ₀ — das ist der eigentliche Kern.
- Nächster Schritt (Vorschlag): (i) als Satz ausarbeiten (Fredholm-Theorie auf Kurven mit Ecken; Literatur über X); (ii) mit W: WKB der Paarfamilie
  (Tiefe c_p N, Einsatz N_p), Ziel „λ_p ≤ λ₁“ und polynomiale Eckanregung für alle p.

Dateien Runde 2: o07b_eckzahl_vek.py mit Modell `radial` (o07b_radial.log, npz), o14_lox_familie.py, o15_anregung.py (+o15_700_1200.log,
o15_1000_1200.log), o16_rest.py, o17_kneser_rest.py; o_kern.py um die Modelle `radial` und `kneser` (A_f aus gemeinsam/) erweitert.


---

# Runde 3 (ab 00:35, Auftrag der Leitung): unendlicher Operator (Fredholm über ∂K'), Bulk-Bedingung, Aufteilung mit W

**Aufteilung (BRETT 00:42):** O = unendlicher Operator T = I + K und Bulk; W = Kante (Eckmoden-Familie, Deckschicht, (E_L)').
Schnittstelle: (Kantenschwanz) ‖Q_N Y_N P_{cN}‖ ≤ C N^b für ein c ≥ 0.7 (Q_N = Projektor auf die r(N) Eck-Ausgänge).

## 13. Der unendliche Operator als Abbildung E²(Ω) → H²(D)

**13.1 Richtiger Definitionsbereich.** Ω := D ∪ f(D) = int K' (Ps Kompaktum). Für δ ∈ H²(D) ist Kδ = (δ∘f)/m nur sinnvoll, wenn δ auf f(D)
fortsetzbar ist; aus (I + K)δ = g folgt umgekehrt δ = (m∘h)(g − δ)∘h auf f(D), d. h. jede Lösung lebt auf Ω. Daher
  **T: E²(Ω) → H²(D),  Tδ := δ|_D + (δ∘f)/m**   (E² = Smirnov-Klasse, Norm = L²(∂Ω, |dz|) der Randwerte),
beschränkt (‖δ|_D‖_{H²(D)} ≲ ‖δ‖_{E²(Ω)}: Carleson für den Bogen Γout, der ∂Ω in L, L̄ unter dem Winkel arg μ trifft; ‖δ∘f‖_{H²(D)} ≍ ‖δ‖_{E²(f(D))}
wegen 0 < c ≤ |f'| ≤ C auf D̄). G := (I+K)⁻¹ ist dann T⁻¹ gefolgt von der Einschränkung auf D; ‖G‖_{H²→H²} ≤ C‖T⁻¹‖.

**13.2 Strukturbedingungen (S)** (für f holomorph um D̄; für Kneser mit der schwachen Singularität in L, L̄, siehe 13.6):
(S1) f univalent auf D̄, m = f'∘f holomorph und nullstellenfrei um D̄; (S2) f(L) = L, f(L̄) = L̄, f'(L) = μ, 0 < β := arg μ < π;
(S3) ∂D ∖ {L, L̄} = Γin ⊔ Γout mit f(Γin) ⊂ D, f(Γout) ⊂ ℂ ∖ D̄; dann ist Ω ein Jordangebiet mit ∂Ω = Γ̄in ∪ f(Γ̄out) und Ecken in L, L̄ mit
Innenwinkel π + β (Aussenwinkel π − β = 141.69° bei Kneser, = Ps α'π); (S4) 1 + μ^{s−1} ≠ 0 auf Re s = −½.

## 14. Lemma O7a — das lokale Modell in der Ecke (bewiesen)

H := {θ₀ < arg ε < θ₀ + π} (D nahe L, ε = x − L), Γin = Strahl θ₀, Γout = Strahl θ₀ + π; f(L+ε) ≈ L + με dreht um β, also
S := H ∪ μH = {θ₀ < arg ε < θ₀ + π + β} (Ω nahe L), ∂S = Γin ∪ μΓout. Lokaler Operator **T_L: E²(S) → H²(H), F ↦ F|_H + F(μ·)/μ.**
**Lemma O7a.** Mit ω(σ) := |1 + μ^{−3/2+iσ}|²·(1 + e^{−2πσ})/(1 + e^{−2(π+β)σ}) gilt: T_L ist genau dann ein Isomorphismus, wenn (S4) gilt, und dann
‖T_L⁻¹‖ = (inf_σ ω)^{−1/2}, ‖T_L‖ = (sup_σ ω)^{1/2}; ω(σ) → 1 (σ → +∞) und → |μ|^{−3} (σ → −∞).
*Beweis.* Mellin–Paley–Wiener im Sektor: F ∈ E²(S) ⇔ F(ε) = ∫F̂(σ)ε^{−1/2+iσ}dσ mit ‖F‖² = 2π∫|F̂|²(e^{−2σθ₀} + e^{−2σ(θ₀+π+β)})dσ
(auf dem Strahl arg ε = θ ist F(re^{iθ}) = ∫F̂ e^{−iθ/2}e^{−σθ}r^{−1/2+iσ}dσ, Mellin–Plancherel in r); analog H²(H) mit Gewicht e^{−2σθ₀} + e^{−2σ(θ₀+π)}.
Für ε ∈ H ist με ∈ S und (με)^s = μ^s ε^s (μ^s := e^{s(ln|μ| + iβ)}, verträglich mit dem Zweig auf S), also (T_LF)^ = (1 + μ^{s−1})F̂. T_L ist damit
Multiplikation mit 1 + μ^{s−1} zwischen zwei gewichteten L²-Räumen; Normen = sup/inf des Quotienten ω (Faktor e^{−2σθ₀} kürzt sich). Grenzwerte:
σ → +∞: Gewichtsquotient → 1, μ^{s−1} → 0; σ → −∞: Quotient ~ e^{2βσ}, |μ^{s−1}|² = |μ|^{−3}e^{−2βσ}. Surjektiv: F̂ := Ĝ/(1 + μ^{s−1}). ∎
*Werte (o20_lokal.py, o20.log):* μ = √L (Kneser, kubisch): inf ω = 0.6090, sup ω = 2.478 ⇒ ‖T_L⁻¹‖ = 1.2814, ‖T_L‖ = 1.5741; radial (μ = 1 + 0.6106i):
1.2804 / 1.6440; Möbius-exp: ‖T_L⁻¹‖ = 1.0000. (S4) ist für alle drei weit erfüllt (min|1 + μ^{s−1}| = 1 − 1e-6.)
*Bemerkung.* Für die Einschränkung auf H (= lokales G) ist die Norm sup|1 + μ^{s−1}|⁻¹ = 1 + 1.37e-6 (Satz O3); der Unterschied zu 1.28 liegt in der
E²(S)-Norm, die auch die Randwerte auf μΓout = f(Γout) misst. Die Oszillation |μ|^{iσ} bei σ → −∞ (Komposition = Dilatation) stört nicht, weil
T_L exakt diagonal ist — sie verhindert nur, das Symbol als stetige Funktion auf [−∞, ∞] zu lesen (wichtig für den Index, 15.2).

## 15. Satz O7 — Fredholm, Index, Invertierbarkeit

**Satz O7 (Fredholm).** Unter (S1)–(S4) ist T: E²(Ω) → H²(D) ein Fredholm-Operator.
*Beweis (Gerüst; die drei lokalen Bausteine sind bewiesen, das Zusammensetzen ist das Standard-Lokalprinzip).*
(i) *Punkte von Γin (fern von L, L̄):* dort ist ∂Ω = ∂D mit derselben Gebietsseite; F ↦ F|_D ist lokal die Identität, und F(f(x))/m(x) ist eine
Auswertung im Innern (f(x₀) ∈ D) ⇒ lokal kompakt (Montel: gleichmässige Konvergenz auf Kompakta). Lokal: T ≡ I.
(ii) *Punkte von Γout:* x₀ liegt im Innern von Ω (Γout ⊂ f(D)), F(x) ist dort Innenauswertung ⇒ lokal kompakt; F ↦ (F∘f)/m bildet die
E²-Randwerte auf f(Γout) ⊂ ∂Ω lokal konform auf H²-Randwerte auf Γout ab (f bildet die D-Seite bei x₀ in Ω ab), mal m⁻¹ ≠ 0 ⇒ lokal invertierbar.
(iii) *Ecken L, L̄:* nach Abbildung auf Tangentialkoordinaten ist T lokal T_L + (Störungen, die in der Ecke verschwinden: f(L+ε) − L − με = O(ε²),
Krümmung von ∂D und f(∂D), m − μ = O(ε)); solche Störungen haben beliebig kleine lokale Norm. Lemma O7a ⇒ lokal invertierbar.
(iv) Lokalprinzip (Simonenko / Allan–Douglas) in der Fassung für gewichtete Shift-Operatoren: Lokalisierung im Bildraum bei x ∈ ∂D gepaart mit
Lokalisierung im Urbild bei x (Γin) bzw. f(x) (Γout) — mit dem Rückzug Φ: E²(Ω) → L²(∂D), (ΦF)|_{Γin} = F, (ΦF)|_{Γout} = F∘f wird das ein
Lokalprinzip für Multiplikationen und Cauchy-Projektionen auf ∂D plus die Mellin-Ecken; lokal invertierbar überall ⇒ Fredholm.
**Status:** (i)–(iii) bewiesen (iii über Lemma O7a), (iv) zitiert (Roch–Santos–Silbermann 2011 Kap. 4–5; Antonevich, Linear Functional Equations,
1996, für Shift-Operatoren) — genaue Fundstellen bei X angefragt (BRETT 00:42).

**15.2 Index 0 per Homotopie (Satz O7b).** Der Index ist längs stetiger Familien von Fredholm-Operatoren konstant (nach Rückzug auf feste
Räume mit Φ_t; die Eckwinkel π + arg μ_t ändern sich stetig). Basis: Möbius-exp f_M, dort T invertierbar (Satz O2: K akkretiv ⇒ ‖(I+K)δ‖ ≥ ‖δ‖
⇒ injektiv mit abgeschlossenem Bild; Surjektivität aus den gleichmässig beschränkten Abschnitten Y_N über schwache Grenzwerte, da die
Zeilenfunktionale (I+K)*e_k ∈ H²(D) sind, weil die k-te Taylorzahl von Kδ nur vom Jet von δ in a₀ ∈ D abhängt) ⇒ Index 0.
Pfad f_t = (1−t)f_M + t f_Z (alle f_t fixieren L, L̄). **Pfadbedingungen numerisch geprüft (o19, t = 0, 0.1, …, 1, Gitter 20000 Randpunkte,
120×480 Innenpunkte):** f_Z = kubisch: min|f_t'| ≥ 0.312, Selbstabstand f_t(∂D) ≥ 0.070 (einfache Kurve ⇒ univalent), |f_t| − ρ wechselt das
Vorzeichen genau in ±76.63° (Γin/Γout unverändert), min|m_t| ≥ 0.733, Pol 3.95 von f_M ≥ 1.45 von D̄ ∪ f_t(D̄), μ_t von e^{0.705i} nach √L mit
min_σ|1 + μ_t^{s−1}| = 1.0000; f_Z = radial: ebenso (min|f_t'| ≥ 0.227, Selbstabstand ≥ 0.051, min|m_t| ≥ 0.691). ⇒ **ind T = 0 für das kubische
und das radiale Modell** (bewiesen bis auf 15.1(iv) und die numerische Pfadprüfung, die auf Gittern und nicht mit Einschluss erfolgte).
*Warum nicht einfacher:* die Homotopie I + λK (λ ∈ (0,1]) erreicht keinen Basisfall: bei λ → 0 verliert T die Fredholm-Eigenschaft (ω → 0 für
σ → −∞) und es kreuzen unendlich viele Nullstellen von 1 + λμ^{s−1} die Gerade Re s = −½ (bei ln λ = (π(2k+1)β + 1.5|ln μ|²)/ln|μ|, k = −1: λ = 1.6e-4);
für λ ∈ (1.6e-4, 4.9e7) ist T_λ Fredholm mit konstantem Index. Ebenso ist das kubische Modell mit reellem μ (Γout leer) kein Basisfall (singuläre Deformation).

**15.3 Invertierbarkeit (Korollar O7c).** Fredholm mit Index 0 + injektiv ⇒ T invertierbar ⇒ ‖G‖ < ∞.
- *Kneser:* Injektivität ist Ds Lemma D4 (Kern von J_∞ trivial unter |δ| ≤ C|x−L|^{−β}|x−L̄|^{−β}, β < 3.447; bewiesen|(G2)): ein Kernelement
  δ ∈ E²(Ω) ist analytisch über Γin (δ = −(δ∘f)/m, f(Γin) ⊂ D) und über Γout (Inneres von Ω) und hat in L nach Lemma O7a eine Entwicklung nach den
  Nullstellen s_k = 1 + (2k+1)γ mit Re s_k > −½, also Re s_k ≥ 5.447 — D4 greift. Offen für Kneser: (S1)–(S3) sind numerisch (P p01: f_K univalent auf D̄,
  A-Bogen ⊂ f_K(D), D-Bogen ∩ f_K(D̄) = ∅), der Homotopie-Pfad braucht f_K ausserhalb D (nur mit Motor: → P), und die Eckenstörung ist nicht
  O(ε²), sondern O(ε^{min(2, 1+Re γ)}) = O(ε²) (Re γ = 4.45) mit Schlitz — (iii) bleibt gültig.
- *Kubisch/radial:* Injektivität **nicht bewiesen** (ein globales Abel-Koordinatensystem wie in D4 fehlt; die lokale Analyse allein kann es nicht
  geben, weil μ^{1−s_k} = −1 exakt ist — lokale Eigenfunktionen existieren). Numerisch: Abschnitte konvergieren eintragsweise, Niederfrequenzblock
  ‖P_nY_NP_n‖ für N = 60…140 bis auf 1e-9 stabil (o17), Rest σ_{r+1} N-unabhängig ⇒ kein Hinweis auf einen Kern.
- *Wert:* ‖G‖ = 1 + δ₀ (numerisch; Kneser 1 + 2.2e-5) — keine explizite Schranke aus der Fredholm-Theorie.

## 16. Bulk-Bedingung: Rückführung auf den Kantenschwanz

**16.1 Korrektur von Lemma O6.** In O6(b) zerfällt QY_N c_hol mit Q = Q_N; nur der Eckteil braucht den Bulk:
‖Q_N Y_N c_hol‖ ≤ ‖Q_N Y_N P_{cN}‖·‖c_hol‖ + σ₁(Y_N)·‖(I − P_{cN})c_hol‖, und ‖c_hol‖ ≤ C_R, ‖(I−P_{cN})c_hol‖ ≤ C_R(ρ/R)^{cN}.
Der Restteil (I − Q_N)Y_N c_hol ist ≤ σ_{r+1}(Y_N)C_R und läuft über (a′)_r. **O6 braucht also nur (Kantenschwanz)** ‖Q_N Y_N P_{cN}‖ ≤ CN^b plus
ln(R/ρ) > c₁/c; die volle Bulk-Bedingung ‖Y_N P_{cN}‖ ≤ ‖Q_NY_NP_{cN}‖ + σ_{r+1}(Y_N) folgt aus (Kantenschwanz) + (a′)_r.
**16.2 Messung (o07b-Vektoren, nur die r Eckmoden):** ‖Q_N Y_N P_{cN}‖ für c = 0.5/0.6/0.7/0.8/0.9:
kubisch N = 700 (r = 4): 0.237/0.422/0.683/1.31/10.2; N = 1000 (r = 4): 0.266/0.463/0.797/2.55/88; N = 1200 (r = 6): 0.318/0.552/0.954/4.40/390;
radial N = 600/800/1000: 0.42/0.52/0.56 bei c = 0.7, 0.78/0.99/1.25 bei c = 0.8. ⇒ für c ≤ 0.7 wächst der Kantenschwanz nur wie ≈ N^{0.5}
(Raten 0.0006–0.0007/Grad bei c = 0.5–0.7), ab c ≈ 0.8 exponentiell (Deckschicht). Das ist eine Aussage über das Abklingen der adjungierten
Eckmoden in die Tiefe (WKB-Ausläufer) — Kante, also W.
**16.3 Was die Fredholm-Theorie zum Bulk beiträgt (Skizze, nicht bewiesen):** ist T invertierbar und konvergieren die Abschnitte stark, so gilt
P_nY_NP_n → P_nGP_n für festes n; eine gleichmässige Aussage für n = cN braucht genau die Trennung „Kante weit weg“, d. h. dass die von der
Abschneidekante erzeugten Moden (Eingänge auf k ≳ (1−c₀)N) in die Tiefe schneller abklingen, als σ_i wächst — der Kantenschwanz von 16.2.
Ein rein operatortheoretischer Beweis ohne die Kantenanalyse ist mir nicht gelungen (Versuch: Y_N g − P_N G g = Y_N P_N K Q_N G g — der Schwanz
Q_N G g wird über Γout exponentiell verstärkt; Kohärenz geht dabei verloren, wie bei der naiven Parametrix §6).

## 17. Stand Runde 3

- Bewiesen: Lemma O7a (lokales Eckmodell, mit Konstanten ‖T_L⁻¹‖ = 1.281 für μ = √L); Satz O7 bis auf das zitierte Lokalprinzip;
  Index 0 für kubisch/radial bis auf das Lokalprinzip und die (numerische) Pfadprüfung; Invertierbarkeit für Kneser zusätzlich bis auf D4|(G2),
  (S) (numerisch, P) und den Kneser-Pfad (→ P).
- Korrigiert: O6 braucht nur (Kantenschwanz) (W) + (a′)_r.
- Offen: Injektivität kubisch/radial (globales Argument fehlt); explizite Schranke für ‖G‖; Bulk = Kantenschwanz (W).

Dateien Runde 3: o18_kubisch_ellip.py (Basispunkt-Tests: kubisch mit |μ| = 1 und arg μ = 0.3 akkretiv bis N = 40, ‖Y_200‖ − 1 = 3e-14; arg μ = 0.5/
0.6686/0.705: Re K ≥ −2.9e-5/−2.9e-3/−5.1e-3, sättigend ⇒ ‖Y_N‖ ≤ 1/(1−ε)), o19_pfad.py (Pfadbedingungen), o20_lokal.py (+o20.log, lokale Konstanten).

## 18. Abgleich mit Strang X (Rechteck-Abschnitte, (a″); BRETT 00:38/00:45, Leitung 00:47)

X: mit d Zusatzzeilen R_d = T_{N+d}J_∞|_{P_N} verschwinden die Eckmoden, σ_min(R_d) = 0.7328 = min_{Γin}|m| ab d ≈ N/60…N/100, und Lemma 3.1
(σ_{d+1}(J_N⁻¹) ≤ 1/σ_min(R_d)). Einordnung:
- **(a″)_∞ = unendlicher Operator.** Für d ≥ 2N (kubisch; dann ist T_{N+d}J_∞δ = J_∞δ) ist (a″) die Aussage ‖J_∞δ‖_{H²(D)} ≥ c‖δ‖_{H²(D)} auf Polynomen.
  Aus Satz O7 + Invertierbarkeit: ‖J_∞δ‖ = ‖m·(I+K)δ‖ ≥ min_{∂D}|m|·‖(I+K)δ‖ ≥ (min_{∂D}|m|/‖G‖)‖δ‖, und min_{∂D}|m| = min_{Γin}|m| = 0.7328 (auf
  Γout ist |m| ≥ 1.17); mit ‖G‖ = 1 + δ₀ erklärt das Xs Wert 0.7328 auf 4 Stellen. Umgekehrt ist (a″)_∞ mit dichtem Bild gleichwertig mit ‖G‖ < ∞.
- **d = ∞ → d(N):** wie viele Überlaufzeilen die Eckmoden einfangen, ist Kantenanalyse (S-Matrix, Streckband) — vorgeschlagen für X/W.
  Eine naive Reduktion über ‖(I − T_{N+d})J_∞δ‖ scheitert (der ferne Überlauf C_fδ jenseits N+d ist für Γout-Inhalt riesig, obwohl dort auch
  ‖J_∞δ‖ riesig ist — die Ungleichung ist dort trivial, aber nicht über eine Normdifferenz).
- **Kette:** (a″)_∞ [O, Satz O7] + (d(N)-Reduktion) [Kante] ⇒ (a″)_{d(N)} ⇒ (Lemma 3.1) (a′)_{r=d(N)} mit Konstante 1/σ_min(R_d) — das ersetzt
  „r(N) Eckmoden abspalten“ durch „d(N) Zeilen hinzunehmen“, ohne Singulärvektoren.


---

# Runde 4 (ab 00:49, Auftrag der Leitung): Lokalprinzip ausgeführt, Index streng, Schnittstellen

## 19. Beweis von Satz O7 (Fassung 2: das Kleben elementar, die klassischen Zutaten mit Fundstelle)

**Bezeichnungen.** Γ = ∂D; Ω = D ∪ f(D) mit ∂Ω = Γ̄in ∪ f(Γ̄out); ψ: Γ → ∂Ω Homöomorphismus (ψ = id auf Γ̄in, ψ = f auf Γ̄out). B: E²(Ω) → L²(∂Ω)
Randspur, P_Ω: L²(∂Ω) → E²(Ω) Cauchy-Projektion, P_D: L²(Γ) → H²(D) Szegő/Cauchy-Projektion. Zerlegung von T (exakt):
  TF = P_D[ ℛ(BF) + 𝒱F ],  ℛu := a·(u∘ψ) mit a = 1 auf Γin, a = 1/m auf Γout;  𝒱F := (F∘f)/m auf Γin, 𝒱F := F auf Γout (Innenwerte),
und ℛ(BF) + 𝒱F ist bereits die Randfunktion der H²(D)-Funktion F|_D + (F∘f)/m (P_D wirkt als Identität).
**Klassische Zutaten (Fundstellen):**
(K1) ∂Ω ist eine Lipschitz-Jordankurve (stückweise analytisch, zwei Ecken mit Innenwinkel π + β ∈ (π, 2π)); der Cauchy-Operator ist auf
     L²(∂Ω) beschränkt [Coifman–McIntosh–Meyer, Ann. of Math. 116 (1982) 361–387] und E²(Ω) = P_ΩL²(∂Ω) (Smirnov-Gebiet) [Duren, Theory of H^p
     Spaces, Kap. 10].
(K2) Für c ∈ C(Γ') (Γ' Carleson-Kurve, insbesondere Lipschitz) ist der Kommutator S_{Γ'}cI − cS_{Γ'} kompakt [A. Yu. Karlovich, J. Integral
     Equations Appl. 15 (2003) 263–320, Lemma 6.5; zitiert nach Karlovich, arXiv:0810.3110, §5.2 — Fundstelle über ar5iv geprüft, Originalbeweis
     nicht gelesen = S]. Für den Kreis klassisch (Hartman).
(K3) Montel: Innenauswertung entlang einer kompakten Teilmenge von Ω ist kompakt E²(Ω) → L²; Cauchy-Integrale mit glatten Kernen sind kompakt.

**Lemma O7d (Kleben, bewiesen).** X, Y Banachräume, T ∈ L(X,Y). Es gebe Φ_1…Φ_n ∈ L(X), Ψ_1…Ψ_n ∈ L(Y) mit ΣΦ_j = I, ΣΨ_j = I,
TΦ_j − Ψ_jT kompakt, und R_j, R'_j ∈ L(Y,X) mit R_jTΦ_j = Φ_j + C_j + S_j, Ψ_jTR'_j = Ψ_j + C'_j + S'_j (C kompakt, Σ‖S_j‖ < 1, Σ‖S'_j‖ < 1).
Dann ist T Fredholm.
*Beweis.* R := ΣR_jΨ_j: RT = ΣR_j(TΦ_j + kompakt) = I + C + S mit ‖S‖ < 1 ⇒ (I+S)⁻¹RT = I + kompakt. R' := ΣΦ_jR'_j: TR' = Σ(Ψ_jT + kompakt)R'_j
= I + C' + S' ⇒ TR'(I+S')⁻¹ = I + kompakt. Links- und Rechtsregularisator ⇒ Fredholm. ∎

**Lokalisierer.** Stetige Zerlegung der Eins auf Γ: χ_L (≡ 1 nahe L), χ_L̄, χ_in (Träger in Γin), χ_out (Träger in Γout), Σχ = 1.
Ψ_j := P_D M_{χ_j} auf H²(D) (Toeplitz), Φ_j := P_Ω M_{χ_j∘ψ⁻¹} auf E²(Ω). Dann ΣΨ_j = I, ΣΦ_j = I exakt.

**Lemma O7e (Vertauschung, bewiesen mit K2, K3).** TΦ_j − Ψ_jT ist kompakt.
*Beweis.* Nach (K2) ist BΦ_jF = (χ_j∘ψ⁻¹)BF + C₁F (C₁ kompakt), also ℛBΦ_jF = χ_jℛBF + ℛC₁F (exakt: ℛ((χ∘ψ⁻¹)u) = χ·ℛu). Mit der Zerlegung
TΦ_jF = P_D[χ_jℛBF] + P_D𝒱Φ_jF + kompakt = Ψ_jTF − P_D[χ_j𝒱F] + P_D𝒱Φ_jF + kompakt. Bleibt 𝒱Φ_j − χ_j𝒱 kompakt: fern von L, L̄ ist 𝒱
selbst kompakt (K3: die Auswertungskurven f(Γin), Γout haben dort positiven Abstand zu ∂Ω); nahe L ist χ_j ≡ c (const) und
Φ_jF = cF + P_Ω((χ_j∘ψ⁻¹ − c)F), wobei χ_j∘ψ⁻¹ − c nahe L verschwindet ⇒ die Innenauswertung dieses Cauchy-Integrals nahe L hat glatten
Kern ⇒ kompakt (K3). ∎

**Lemma O7f (lokale Inverse, bewiesen bis auf K1–K3).**
(in) Für u mit Träger in einem Teilbogen von Γin sind P_D u und P_Ω u dasselbe Cauchy-Integral (gleicher Bogen, gleiche Formel), also ist
     T lokal die Identität: R_in g := P_Ω(χ'·g) (χ' etwas grössere Abschneidung) erfüllt R_inTΦ_in − Φ_in kompakt, Ψ_inTR_in − Ψ_in kompakt.
(out) R_out g := P_Ω(χ'·(m g)∘h): Konforme Invarianz modulo glatter Kerne — f'(x)/(f(x)−f(y)) − 1/(x−y) ist glatt (f univalent, analytisch) ⇒
     R_outTΦ_out − Φ_out und Ψ_outTR_out − Ψ_out kompakt (𝒱 ist auf Γout fern von L kompakt).
(L) **Koenigs-Koordinaten:** für |μ| ≠ 1 und f holomorph bei L gibt es κ holomorph nahe L, κ(L) = 0, κ'(L) = 1, κ∘f = μκ (Linearisierung
     des anziehenden Zweigs h = f⁻¹; für Möbius-f mit |μ| = 1 leistet das die Möbius-Abbildung M(x) = (x−L)/(x−L̄); für Kneser leistet es
     Knesers eigene Koordinate e^{L·slog} — dort ist f_K exakt die Multiplikation mit √L, holomorph auf Ω nahe L, da der Schlitz ausserhalb Ω liegt).
     In κ-Koordinaten ist die Komposition EXAKT die Dilatation z ↦ μz; die Kurven κ(Γin), μκ(Γout) sind analytische Bögen durch 0 mit den
     Tangentenstrahlen des Modells (Lemma O7a). Abweichung Kurve ↔ Strahl: die Cauchy-Kerne unterscheiden sich lokal um k(s,t) = O((|s|+|t|)^{α−1})
     (α = 1: analytische Bögen; bei Kneser ist κ(x) = C(x−L)(1 + a₁(x−L) + … + b(x−L)^γ + …) wegen Re γ = 4.45 > 1 ebenfalls C^{1,1}) ⇒ Hilbert–Schmidt-Norm auf {|s|,|t| < r} ist O(r^α) → 0; ebenso die Innenauswertungen
     (Auswertungskurven im Abstand ≳ |z|·sin β vom Rand, Abweichung O(|z|^{1+α})) und der Faktor 1/m − 1/μ = O(|z|). Also T lokal bei L =
     (Konjugation mit dem unitären Rückzug F ↦ (F∘κ⁻¹)((κ⁻¹)')^{1/2}) ∘ T_L ∘ (…) + (Norm ≤ C r^α) + kompakt. Lemma O7a ⇒ lokale Inverse mit
     beliebig kleinem S-Rest (r klein wählen).
⇒ **Satz O7 (Fredholm) ist bewiesen**, gestützt nur auf K1–K3 (klassisch, zitiert) und die eigenen Lemmata O7a, O7d–O7f.
*Ehrlicher Rest:* (i) K2 habe ich nur über die Sekundärquelle geprüft (ar5iv-Volltext von arXiv:0810.3110, dort „[17, Lemma 6.5]“); (ii) die
Carleson-Einbettung für die Spur von E²(Ω) auf Γout (Beschränktheit von T, Schritt 0) ist Standard, aber hier nicht ausgeschrieben.

## 20. Index 0 streng (Satz O7b): Pfadbedingungen mit arb zertifiziert

**Reduktion auf Randbedingungen (bewiesen).** Für f_t = (1−t)f_M + t f_Z (alle reell, f_t(L) = L exakt):
 (P1) Re f_t' > 0 auf ∂D ⇒ (harmonisches Minimum) auf D̄ ⇒ f_t univalent auf einer konvexen Umgebung von D̄ (Noshiro–Warschawski);
 (P2) Re m_t > 0 auf ∂D (m_t = f_t'∘f_t holomorph auf D̄ wegen P3) ⇒ m_t ≠ 0 auf D̄;
 (P3) Re f_t < p (Pol von f_M) auf ∂D ⇒ auf D̄ ⇒ f_t(D̄) ∌ p ⇒ f_t, m_t holomorph um D̄;
 (P4) g_t = |f_t|² − ρ² > 0 auf [0, θ_L−η], < 0 auf [θ_L+η, π], g_t' < 0 auf [θ_L−η, θ_L+η] (g_t(θ_L) = 0 exakt) ⇒ mit der Symmetrie f_t reell:
      Γin/Γout genau die Bögen jenseits bzw. diesseits ±arg L (arg L = Im L = 1.33724 exakt, da L = e^L);
 (P5) X_t = π·Im Log μ_t/|Log μ_t|² > 3/2 und 0 < arg μ_t < π ⇒ die Nullstellen s_k = 1 + (2k+1)X_t + i(…) von 1 + μ_t^{s−1} liegen nie auf Re s = −½.
**Zertifikat (o21_zertifikat.py, o21.log; Ball-Arithmetik 120 Bit, adaptive Bisektion in (t, φ) ∈ [0,1] × [0,π], Bälle mit Rundungszuschlag):**
kubisch: P1 ✓ (64 Rechtecke), P2 ✓ (584), P3 ✓ (64), P4a ✓ (1272), P4b ✓ (690), P4c ✓ (64, η = 0.02), P5 ✓ (min X_t ≥ 4.369);
radial: P1 ✓, P2 ✓ (514), P3 ✓, P4a ✓ (1302), P4b ✓ (784), P4c ✓, P5 ✓ (min X_t ≥ 4.381). Gegenprobe: die falsche Bedingung Re f_t' > 0.6 wird
nicht zertifiziert (Abbruch an einer konkreten Stelle) — das Verfahren ist nicht trivial wahr.
**Index-Konstanz (Skizze mit Fundstelle):** Rückzug Φ_t: E²(Ω_t) → 𝓗_t := Ran P_t ⊂ L²(Γ), P_t = Φ_tP_{Ω_t}Φ_t⁻¹; P_t und T̃_tP_t sind in t normstetig
(Kerne der zurückgezogenen Cauchy-Operatoren hängen glatt von t ab, Differenzen sind Calderón-Kommutatoren mit Lipschitz-Symbol O(|t−t'|);
an den Ecken Mellin-Faltungen mit stetig vom Winkel π + arg μ_t abhängendem Symbol). Mit U_t = P_tP_{t₀} + (I−P_t)(I−P_{t₀}) (für t nahe t₀
invertierbar, 𝓗_{t₀} → 𝓗_t) ist T̃_tU_t eine normstetige Familie zwischen festen Räumen ⇒ Index lokal konstant [Kato, Perturbation Theory for Linear Operators, Kap. IV §5 (Stabilität des Index), U].
Basis t = 0: Möbius-exp ist invertierbar (O2: ‖(I+K)δ‖ ≥ ‖δ‖ auf Polynomen ⇒ injektiv mit abgeschlossenem Bild, weil
‖δ‖_{E²(Ω)} ≲ ‖δ‖_{H²(D)} + ‖m(Tδ − δ)‖ ≤ C‖Tδ‖; surjektiv per schwachem Grenzwert der gleichmässig beschränkten Abschnitte Y_N g, da die
Zeilenfunktionale δ ↦ [x^k]Tδ stetig auf H²(D) sind — sie hängen nur vom Jet in a₀ ∈ D ab). ⇒ **ind T = 0 für kubisch und radial.**
*Status:* Pfad (S) streng (arb); Normstetigkeit der Familie: Skizze (Calderón-Kommutatoren) — der einzige nicht ausgeschriebene Schritt.

## 21. Injektivität und Schranke für ‖G‖ — Stand

- **Kneser:** injektiv nach D4 (bewiesen|(G2)); ein Kernelement δ ∈ E²(Ω) ist analytisch über Γin und Γout und hat nach Lemma O7a in L nur
  Exponenten mit Re s ≥ 5.447 (keine Nullstelle von 1 + μ^{s−1} in −½ < Re s < 5.447), erfüllt also Ds Wachstumsbedingung. Für ind = 0 bei
  Kneser fehlt der Pfad Möbius-exp → f_K (braucht f_K ausserhalb D) **→ P (Motor):** bitte P1–P5 für f_t = (1−t)f_M + t f_K prüfen (P1–P4 auf ∂D,
  P5 analytisch: μ_t = (1−t)e^{0.705104i} + t√L, X_t ≥ 4.36 wie kubisch). Dann: T_Kneser invertierbar ⇐ (G2) + diese Prüfung.
- **Kubisch/radial:** Injektivität **offen** — die lokale Analyse kann sie nicht liefern (μ^{1−s_k} = −1 exakt: lokale Eigenfunktionen existieren),
  und ein globales Abel-Koordinatensystem wie bei Kneser (D4) fehlt für die Modelle. Numerisch kein Kern (Niederfrequenzblock stabil, o17).
  Folge: für die Modelle ist „T Fredholm mit Index 0“ bewiesen, „T invertierbar“ hängt an der Injektivität.
- **Explizite Schranke für ‖G‖:** nicht erreicht. Grund: Akkretivitäts-/Wertebereichsschranken sind bei |μ| ≠ 1 unmöglich (O3), die
  Fredholm-Theorie liefert keine Konstante, und eine rechnergestützte Unterschranke inf ‖J_∞δ‖/‖δ‖ über ALLE δ braucht eine quantitative
  Fassung von Lemma O7f für hohe Frequenzen (endlich viele Frequenzen allein geben nur obere Schranken der Konstante). Plan für später:
  δ = δ_{≤n} + δ_{>n}; für δ_{>n} lokale Schranken mit Fehlern O(r) + O(1/(nr)) (Mellin bei L, Dämpfung auf Γin, Dominanz auf Γout), für δ_{≤n}
  eine arb-Rechnung; Kopplungsterme über (K2)-Kommutatoren quantitativ.

## 22. Schnittstellen

**An X: (a″)_∞ ⇔ ‖G‖ < ∞ (Satz O8, bewiesen).** Sei c_∞ := inf_{δ Polynom} ‖J_∞δ‖_{H²(D)}/‖δ‖_{H²(D)}. Dann
(i) ist T invertierbar, so c_∞ ≥ min_{∂D}|m|/‖G‖, G = (I+K)⁻¹ = Einschränkung∘T⁻¹;
(ii) ist c_∞ > 0 und ind T = 0 (Satz O7b), so ist T invertierbar und ‖G‖ ≤ max_{∂D}|m|/c_∞.
*Beweis.* (i) Für Polynome δ ∈ E²(Ω): ‖J_∞δ‖ = ‖m·Tδ‖_{L²(Γ)} ≥ min|m|·‖Tδ‖ und ‖δ‖_{H²(D)} = ‖G Tδ‖ ≤ ‖G‖‖Tδ‖. (ii) ‖Tδ‖ = ‖J_∞δ/m‖ ≥ c_∞‖δ‖/max|m|
auf Polynomen; Polynome sind dicht in E²(Ω) (Smirnov-Gebiet, K1) und T stetig ⇒ ‖Tδ‖ ≥ (c_∞/max|m|)‖δ‖_{H²(D)} auf E²(Ω); mit
‖δ‖_{E²(Ω)} ≲ ‖δ‖_{H²(D)} + ‖m(Tδ − δ)‖_{H²(D)} folgt Injektivität mit abgeschlossenem Bild; Index 0 ⇒ surjektiv; ‖Gg‖ ≤ (max|m|/c_∞)‖g‖. ∎
Kubisch: min_{∂D}|m| = min_{Γin}|m| = 0.7328, max_{∂D}|m| = ‖m‖_∞ (auf Γout); Xs Messung c_∞ = 0.7328 heisst ‖G‖ = 1 + O(δ₀) (verträglich).
Die Reduktion d = ∞ → d(N) (Rechteck-Abschnitte) ist Xs Teil (Kante).

**An W: präzises Ziel (Kantenschwanz).** Für jedes N sei r(N) := #{i : σ_i(Y_N) > 1 + 2δ₀} und Q_N der Orthogonalprojektor auf span{u_1,…,u_{r(N)}}
(linke Singulärvektoren von Y_N = (I+K_N)⁻¹ in ℓ²_ρ). **(KS)** Es gibt c₀ ∈ [0.7, 1), b < ∞, C mit ‖Q_N Y_N P_{c₀N}‖_{ℓ²_ρ→ℓ²_ρ} ≤ C N^b für alle N.
Äquivalent: Σ_{i≤r(N)} σ_i² Σ_{k≤c₀N}|v_{i,k}|² ≲ N^{2b} (die adjungierten Eckmoden fallen bis zur Tiefe (1−c₀)N mindestens wie 1/σ_i ab).
Gemessen b ≈ 0.5 bei c₀ = 0.7 (kubisch bis N = 1200, radial bis 1000), exponentiell bei c₀ = 0.9. Verwendet in Lemma O6 (Lokalisierung für A_f)
zusammen mit (a′)_r. Kanonische (SVD-freie) Variante passend zu Xs Rechteck-Abschnitten: ‖R⁺_bot S⁻¹ B R⁺_top P_{c₀N}‖ ≤ CN^b.


---

# Runde 5 (ab 01:01, Auftrag der Leitung): Lücken in O7 geschlossen; (a″)_∞ für Kneser

## 23. Die drei Lücken in Satz O7

**23.1 Lemma O7g (Kommutator; eigener Beweis, ersetzt das Karlovich-Zitat).** Γ' rektifizierbare Jordankurve mit Sehnen-Bogen-Bedingung
|ζ − z| ≥ κ·ℓ(ζ,z) (ℓ = kürzere Bogenlänge), S_{Γ'} beschränkt auf L²(Γ'). Dann ist cS_{Γ'} − S_{Γ'}c kompakt für jedes c ∈ C(Γ').
*Beweis.* Für c Lipschitz bzgl. Bogenlänge (Konstante λ) hat cS − Sc den Kern (1/πi)(c(z) − c(ζ))/(ζ − z) mit |·| ≤ λℓ(ζ,z)/(π|ζ − z|) ≤ λ/(πκ):
beschränkter Kern auf einem Raum endlichen Masses ⇒ Hilbert–Schmidt ⇒ kompakt. Lipschitz-Funktionen liegen gleichmässig dicht in C(Γ'),
‖[S, c] − [S, c_n]‖ ≤ 2‖S‖‖c − c_n‖_∞ → 0, und die kompakten Operatoren sind abgeschlossen. ∎
∂Ω ist stückweise analytisch mit zwei Ecken der Innenwinkel π + β ∈ (π, 2π) ⇒ Sehnen-Bogen; S_{∂Ω} beschränkt nach Coifman–McIntosh–Meyer
(Ann. of Math. 116 (1982)) bzw. David (Ann. Sci. ÉNS 17 (1984)) für Sehnen-Bogen-Kurven [U: aus dem Gedächtnis, Standard].

**23.2 Lemma O7h (Spuren im Innern; T beschränkt; bewiesen bis auf die Eckenasymptotik konformer Abbildungen).** Für jeden analytischen Bogen γ ⊂ Ω̄,
der ∂Ω nur in L, L̄ und dort nichttangential trifft (γ = Γout, γ = f(Γin)), gilt ∫_γ|F|²|dz| ≤ C_γ‖F‖²_{E²(Ω)} für F ∈ E²(Ω). Daraus
‖Tδ‖_{H²(D)} ≤ (1 + C_{Γout}^{1/2} + ‖1/m‖_∞‖1/f'‖_∞^{1/2}(1 + C_{f(Γin)}^{1/2}))‖δ‖_{E²(Ω)}.
*Beweis.* (i) Mittelwert: |F(z)|² ≤ (πr_z²)⁻¹∫_{B(z,r_z)}|F|²dA mit r_z = c_γ|z − L|/2 (nichttangential: dist(z, ∂Ω) ≥ c_γ|z − L| nahe L; fern der Ecken
ist dist ≥ d₀ > 0 und alles trivial). Fubini: ∫_{γ∩U}|F|²|dz| ≤ C∫_{Ω∩U}|F(w)|²dA(w)/|w − L|. (ii) Mit der Riemann-Abbildung Φ: Ω → 𝔻,
ζ_L = Φ(L): F = H(Φ)·(Φ')^{1/2}, H ∈ H²(𝔻), ‖H‖ = ‖F‖_E. An der Ecke (Innenwinkel π + β) ist w − L ≍ (ζ − ζ_L)^a, a = (π+β)/π, |dw/dζ| ≍ |ζ − ζ_L|^{a−1}
[Pommerenke, Boundary Behaviour of Conformal Maps, Thm 3.9, U], also ∫|F|²dA(w)/|w − L| ≍ ∫|H|²·|ζ − ζ_L|^{1−a+2(a−1)−a}dA(ζ) = ∫|H|²|ζ − ζ_L|⁻¹dA(ζ).
(iii) μ = |ζ − ζ_L|⁻¹dA ist ein Carleson-Mass (Box der Kantenlänge h an ζ_L: μ ≍ ∫_0^h r⁻¹·r dr·π = πh; fern von ζ_L: μ ≤ h²/dist ≤ h) ⇒
∫|H|²dμ ≤ C‖H‖²_{H²} [Carleson]. ∎
*Bemerkung.* Im exakten Sektor S geht (ii)–(iii) ohne Riemann-Abbildung: ∫_S|F|²dA/|w| = ∫_θ∫_r|F(re^{iθ})|²dr dθ ≤ (π+β)·2π∫|F̂(σ)|²max_θ e^{−2σθ}dσ
≤ (π+β)‖F‖²_{E²(S)} (Mellin wie in O7a; max über den Sektor am Rand).

**23.3 Lemma O7i (Normstetigkeit der Familie; bewiesen bis auf den Calderón-Kommutatorsatz).** ψ_t: Γ → ∂Ω_t (id auf Γin, f_t auf Γout),
Φ_t v := (v∘ψ_t)|ψ_t'|^{1/2} (unitär L²(∂Ω_t) → L²(Γ)), P_t := Φ_tP_{Ω_t}Φ_t⁻¹. Dann ‖P_t − P_s‖ ≤ C|t − s| und ebenso für T̃_tP_t.
*Beweis.* P_Ω = ½(I + S_∂Ω). Der Kern von Φ_tS_tΦ_t⁻¹ ist k_t(x,y) = (πi)⁻¹|ψ_t'(x)|^{1/2}ψ_t'(y)|ψ_t'(y)|^{−1/2}/(ψ_t(y) − ψ_t(x)).
Mit A := ψ_t − ψ_s (= 0 auf Γin, = f_t − f_s = (t−s)(f_Z − f_M) auf Γout; A(L) = A(L̄) = 0, also global Lipschitz mit ‖A'‖ ≤ C|t − s|):
1/(ψ_t(y) − ψ_t(x)) − 1/(ψ_s(y) − ψ_s(x)) = (A(x) − A(y))/((ψ_t(y)−ψ_t(x))(ψ_s(y)−ψ_s(x))) — ein erster Calderón-Kommutator auf einer
Sehnen-Bogen-Kurve, Norm ≤ C‖A'‖_∞ [Calderón, PNAS 53 (1965); auf Sehnen-Bogen-Kurven Coifman–McIntosh–Meyer 1982, U]; die Faktoren |ψ'|^{±1/2}, ψ'
sind Lipschitz in t und multiplizieren einen gleichmässig beschränkten Cauchy-Operator (die Sehnen-Bogen-Konstanten sind auf [0,1] gleichmässig,
weil arg μ_t ∈ [0.6686, 0.7051] die Eckwinkel von 0 und 2π fernhält). Für T̃_t: ℛ_t ist Multiplikation mit a_t (Lipschitz in t);
𝒱_t sind Cauchy-Integrale an Innenpunkten w_t(x) mit Kern ψ_t'(y)/(ψ_t(y) − w_t(x)); deren t-Differenz ist O(|t−s|)/|ψ(y) − w(x)| und
an den Ecken ein homogener Kern vom Grad −1 auf zwei Strahlen, dessen L²-Norm (Mellin/Hardy) endlich ist. ∎
**Index-Konstanz.** U_t := P_tP_{t₀} + (I − P_t)(I − P_{t₀}) ist für |t − t₀| klein invertierbar und bildet Ran P_{t₀} auf Ran P_t ab; T̃_tU_t|_{Ran P_{t₀}} ist
normstetig zwischen festen Räumen ⇒ Index lokal konstant [Kato, Kap. IV §5, U]. Mit O7d–O7h und der arb-zertifizierten Pfadbedingung (§20):
**ind T = 0 entlang des ganzen Pfades.**

**Stand Satz O7 nach Runde 5:** bewiesen; an Fremdergebnissen nur noch (U, Standard): Beschränktheit des Cauchy-Operators und des ersten
Calderón-Kommutators auf Sehnen-Bogen-Kurven (CMM 1982/David 1984), Eckenasymptotik konformer Abbildungen (Pommerenke Thm 3.9), Carleson-
Einbettung, Index-Stabilität (Kato). Das Karlovich-Zitat ist durch O7g ersetzt.

## 24. (a″)_∞ für Kneser — qualitativ bewiesen, rechnergestützt bis N = 100

**Korollar O7j (Kneser).** Unter (G2) (Ds/Bs Geometrie von slog(D)), den Strukturbedingungen (S) für f_K (P p01, numerisch) und P1–P4 für den Pfad
f_t = (1−t)f_M + t f_K (→ P, Motor; P5 ist zertifiziert, §20) ist T: E²(Ω) → H²(D) invertierbar (O7: Fredholm; O7b: Index 0; D4: injektiv), also
**c_∞ = inf_{δ Polynom} ‖J_∞δ‖_{H²(D)}/‖δ‖_{H²(D)} ≥ min_{∂D}|m|/‖G‖ > 0** (Satz O8). Das ist (a″) mit d = ∞, Xs Bulk-Baustein in Satz X1. Für Kneser
lebt (O7f) in Knesers eigener Koordinate e^{L·slog}: dort ist f_K exakt die Multiplikation mit √L, κ(x) = C(x − L)(1 + O(x−L)) ist C^{1,1} (Re γ > 1),
und der Schlitz liegt ausserhalb Ω. Konstante dabei nicht explizit.
**Rechnergestützt (o22, o22.log):** Rechteck-Abschnitt R_N = T_{3N}J_∞|_{P_N} mit Knesers Koeffizienten als arb-Bällen (Radius 10× Strang As
Fehlerschätzung — die Koeffizientenfehler selbst sind nicht zertifiziert), Zertifikat: genäherte Eigenbasis Q (mpmath, 90 Stellen), M = Qᵀ(RᵀR − c²I)Q
in arb strikt diagonaldominant ⇒ σ_min(R_N) > c. Da Abschneiden die Norm verkleinert, ist σ_min(R_N) ≤ c_N := inf_{P_N}‖J_∞δ‖/‖δ‖:
| N | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|
| σ_min(R_N) (Mittelwert) | 0.732274 | 0.731958 | 0.731925 | 0.731915 | 0.731911 |
| zertifiziert σ_min(R_N) > | 0.7315 | 0.7315 | 0.7315 | 0.7315 | 0.7315 |
(Max. Radius in RᵀR: 1e-36 / 8e-28 / 3e-19 / 1e-10 / 4e-2; N = 120: Mittelwert 0.731909, Radius 1e7 — nicht mehr zertifizierbar. Die 46-stelligen
Kneser-Daten tragen bis N = 100; darüber braucht es Koeffizienten mit mehr Stellen und zertifizierten Fehlern (P).)
**Struktur von c_∞ (numerisch):** c_N fällt gegen ≈ 0.73189 = min_{∂D}|m|/(1 + δ₀) (min_{∂D}|m| = 1/1.36630 = 0.731903, angenommen bei φ = π auf Γin;
δ₀ = 2.2e-5 Niederfrequenz-Effekt). Hohe Frequenzen haben den Wert min_{Γin}|m| (Wellenpakete bei φ = π, δ∘f dort gedämpft), die Ecke L den grösseren
Wert |μ|·inf|1 + μ^{s−1}| = 1.1724·(1 − 1.4e-6), Γout-Inhalt beliebig grosse Werte.

**Plan für die explizite Schranke c_∞ ≥ c* > 0 (Zertifikat, Bausteine):**
B1 Niederfrequenzteil P_n: arb wie o22 (erledigt bis n = 100 mit Strang-A-Daten; mit zertifizierten Koeffizienten von P bis n ≈ 200).
B2 Eckmodell: Lemma O7a mit exakten Konstanten (erledigt: inf ω = 0.609, Werte in o20.log).
B3 Eckstörung in Knesers Koordinate κ: explizite Schranke für die Kernabweichung Kurve ↔ Strahl auf |z| < r (O(r) mit Konstante aus κ'' — braucht
   rigorose Werte von slog nahe L → P/A).
B4 Bulk hoher Frequenzen (δ = x^{n+1}δ̃): auf Γin Dämpfung |δ(f(x))| ≤ ‖δ‖(|f(x)|/ρ)^{n+1}(1 − |f(x)|²/ρ²)^{−1/2}; auf Γout Zwei-Konstanten für δ̃ in Ω
   (Randwerte auf f(Γout) gedämpft um (ρ/|y|)^{n+1}) — braucht explizite harmonische Masse von Ω (Riemann-Abbildung mit Einschluss).
B5 Kleben: Kosinus zwischen J_∞(P_n) und J_∞(P_n^⊥) explizit (Glattheit von J_∞e_j, j ≤ n, gegen Oszillation hoher Frequenzen; Stationärphasen-Schranken).
Ehrlich: B3–B5 sind zusammen ein eigenes Projekt (mehrere Tage); qualitativ ist (a″)_∞ für Kneser mit Korollar O7j erledigt.

## 25. Injektivität in den Modellen (nachrangig)

Offen wie in §21. Für Ms radiales Modell (Gårding-Weg) frage ich M über das BRETT, ob die Injektivität dort gebraucht wird; die Fredholm-Eigenschaft
mit Index 0 ist für das radiale Modell bewiesen (§20 zertifiziert, §23).


---

# Runde 6 (ab 01:40, Auftrag der Leitung): explizite Schranke c_∞ ≥ c* — Umformulierung als Formungleichung

## 26. Satz O9 — Reduktion von c_∞ auf eine hermitesche Formungleichung (FA)

**Bezeichnungen.** Für Polynome δ auf ∂D: Kδ := (δ∘f)/m, Tδ := δ + Kδ (= J_∞δ/m), und **Bδ := δ auf Γin, Kδ auf Γout** — das sind die Randwerte
von δ auf ∂Ω = Γin ∪ f(Γout), auf ∂D zurückgezogen (‖Bδ‖ ≍ ‖δ‖_{E²(Ω)}). Gewicht p = |m|² auf ∂D, ‖g‖²_p := ∫_{∂D} p|g|² dφ/2π; dann ist
‖Tδ‖_p = ‖J_∞δ‖_{H²(D)}.
**(FA_λ):** Re⟨Tδ, Bδ⟩_p ≥ λ‖Bδ‖²_p für alle Polynome δ.
Wegen Re⟨Tδ, Bδ⟩_p = ‖Bδ‖²_p + X(δ) mit **X(δ) := Re⟨Kδ, δ⟩_p** (Integral über den ganzen Kreis; auf Γin ist Kδ die Innenspur, auf Γout δ)
heisst (FA_λ): **X ≥ −(1−λ)‖B·‖²_p** — der negative Teil des Wertebereichs von K ist klein, gemessen in der E²(Ω)-Norm (nicht in H²(D)).
**Satz O9 (bewiesen).** Gilt ‖Tδ‖_p ≥ λ‖Bδ‖_p für alle Polynome δ (z. B. als Folge von (FA_λ) per Cauchy–Schwarz), λ > 0, so ist c_∞ ≥ c(λ) := [max_{0≤b≤1/λ} ((λ⁻² − b²)/m_in² + (1 + b)²/m_out²)]^{−1/2},
m_in := min_{Γin}|m|, m_out := min_{Γout}|m|; grob c(λ) ≥ min_{∂D}|m|·(λ⁻² + (1 + λ⁻¹)²)^{−1/2}.
*Beweis.* t := ‖J_∞δ‖ = ‖Tδ‖_p ≥ λ‖Bδ‖_p (unter (FA_λ): t‖Bδ‖_p ≥ Re⟨Tδ, Bδ⟩_p ≥ λ‖Bδ‖²_p). Mit a := ‖mδ‖_{L²(Γin)} = ‖Bδ‖_{p,Γin},
b := ‖δ∘f‖_{L²(Γout)} = ‖Bδ‖_{p,Γout} ist a² + b² ≤ t²/λ²; auf Γin |δ| ≤ |mδ|/m_in; auf Γout mδ = J_∞δ − δ∘f ⇒ ‖mδ‖_{Γout} ≤ t + b.
Also ‖δ‖² ≤ a²/m_in² + (t + b)²/m_out² ≤ t²·(Maximum oben). ∎ (Kein Carleson-Mass nötig: die Innenspur auf Γout wird über mδ = J_∞δ − δ∘f kontrolliert.)
**Kneser:** m_in = 0.731905 (φ = π), m_out = |μ| = 1.172416 (in den Ecken), max_{Γout}|m| = 2.338. c(0.98988) = 0.568, c(0.9) = 0.535,
c(0.75) = 0.471, **c(0.5) = 0.340**, c(0.25) = 0.179 (grobe Fassung: 0.325/0.307/0.272/0.203/0.114).

## 27. Lemma O10 — die Ecke erfüllt (FA) mit λ_L = 0.98988 (bewiesen, arb)

Im Eckmodell von Lemma O7a (Sektor S, Mellin F̂(σ), s = −½ + iσ, z := μ^{s−1}, Strahlgewichte a = e^{−2σθ₀}, b = e^{−2σ(θ₀+π)}; p ist in der
Ecke konstant und kürzt sich) ist Bδ ↔ (F̂ auf Γin, zF̂ auf Γout), X ↔ Re z·(a + b)|F̂|², also
  Re⟨T_LF, BF⟩ = 2π∫|F̂|²[a(1 + Re z) + b(|z|² + Re z)]dσ,  ‖BF‖² = 2π∫|F̂|²(a + b|z|²)dσ.
**Lemma O10.** Das Modell erfüllt (FA_λ) genau für λ ≤ λ_L := inf_σ q(σ), q = 1 + Re z·(1 + r)/(1 + r|z|²), r = e^{−2πσ}, und
**λ_L = 0.9898822 für μ = √L** (Minimum bei σ = −5.04; Zertifikat o25: Intervall-Bisektion auf [−60, 60], Schwänze |q − 1| ≤ 1e-17). sup q = 1.546.
*Bedeutung.* Das Hindernis aus Satz O3 (Re z < −1 für σ ∈ (−23, −3.6): K ist in keiner Mellin-diagonalen H²-Norm akkretiv) verschwindet,
weil ‖B‖² den Anteil b|z|² (Randwerte auf f(Γout)) enthält: dort dominiert |z|² + Re z ≥ |z|² − |z|. Im elliptischen Fall (|μ| = 1) liegt z auf
einem festen Strahl mit Re z ≥ 0 ⇒ λ_L = 1. Folgerung (Skizze, Weyl-Folgen in der Ecke wie in O7): λ_∞ := sup{λ : (FA_λ)} ≤ λ_L.

## 28. Numerik (o24_testform.py, o24_kneser.log): (FA) global

λ_N := min über P_N von Re⟨Tδ, Bδ⟩_p/‖Bδ‖²_p (verallgemeinertes Eigenwertproblem; Gram-Matrizen aus Gauss-Quadratur auf Γout und Γin
mit Randverdichtung, in arb assembliert, Eigenwerte mpmath), p = |m|²:
| N | 10 | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|---|
| λ_N Kneser | 1.00170 | 0.999843 | 0.997973 | 0.996755 | 0.995974 | 0.995425 |
| λ_N kubisch | 1.00132 | 0.999829 | 0.997896 | 0.996658 | 0.995870 | 0.995319 |
X/‖B‖² ∈ [−0.0046, 1.5015] (N = 100). Kontrolle: min ‖J_∞δ‖/‖δ‖ über P_N = 0.731958 (N = 40), 0.7319113 (N = 100) — gleich Os σ_min(R_N) aus o22 ✓.
Mit p = 1: λ_40 = 0.99934 und min ‖(I + K)δ‖/‖δ‖ = 0.999986 (= 1/‖G‖, vgl. ‖G‖ = 1 + 2.2e-5).
**Befund:** λ_N fällt langsam gegen den Eckwert λ_L = 0.98988 (feinere Eckskalen werden aufgelöst); kein Hinweis auf einen diskreten Wert unterhalb.
**Vermutung:** λ_∞ = λ_L, also c_∞ ≥ c(λ_L) = 0.568 (wahr: 0.7319). Rigoros wird c(λ) für jedes bewiesene λ; Ziel λ = ½ ⇒ c* = 0.340.

## 29. Was für einen Beweis von (FA_λ) mit explizitem λ fehlt (Plan B3'–B5', ersetzt B3–B5 aus §24)

Gegenüber §24 ist die Lage besser: statt einer Unterschranke für einen nicht-normalen Operator ist eine **hermitesche Form** nach unten
abzuschätzen, mit Spielraum 1 − λ (bis 0.99 statt 0), und ihr wesentliches Spektrum ist exakt bekannt (O10).
- **Gute Koordinate:** M(x) := (x − L)/(x − L̄) bildet ∂D auf eine Gerade durch 0 ab; in ζ = log M (bis auf Drehung) ist **D ein EXAKTER Streifen
  {θ₀ < Im ζ < θ₀ + π}** — für jedes Modell (Kneser, kubisch, radial), Γin = untere, Γout = obere Kante. f wird ζ ↦ g(ζ) = ζ + ℓ + O(e^{Re ζ})
  am L-Ende (ℓ = log μ), = ζ − ℓ̄ + O(e^{−Re ζ}) am L̄-Ende; nur die obere Kante von Ω (= f(Γout)) ist gekrümmt (Linse der Höhe → β an beiden Enden).
  Die B-Daten liegen auf den zwei GERADEN Kanten des D-Streifens (δ unten, Kδ oben); Fourier längs des Streifens ist exakt.
- **B3' (Ecke, quantitativ):** am L-Ende Begradigung der oberen Ω-Kante durch ψ = id + O(e^{Re ζ}) (Streifen-Poissonkern, explizit), danach ist
  (FA) das Modell O10 bis auf Störungen O(e^{Re ζ}) in Innenauswertungen (Cauchy-Abschätzungen im Streifen: Abstand ≥ min(β, π)) und Jacobi-
  Faktoren 1 + O(e^{Re ζ}). Liefert Q(φ) ≥ (λ_L − λ − C₃r)‖Bφ‖² für φ, die in |x − L| < r leben (Lokalisierung siehe B4').
- **B4' (Kleben):** IMS-Zerlegung mit analytischen Abschneidern in log-Skala (χ_j analytisch ⇒ χ_jδ ∈ E²(Ω)); weil Σ_j χ_j(ζ+ℓ)χ̄_j(ζ) − 1 =
  O(|ℓ|²/W²) + O(|ℓ|H/W²) (erste Ordnung hebt sich weg), ist der Fehler zweiter Ordnung in 1/W; er geht als η(‖Tδ‖_p + 2‖Bδ‖_p)² ein und wird
  absorbiert (bei η = 0.03, λ = 0.9 bleibt ‖Tδ‖ ≥ 0.68‖Bδ‖). Kosten: Übergangslänge W ≈ 5–8 in log-Skala.
- **B5' (Mitte):** endliche Rechnung auf dem Mittelstück. Die Polynom-Basis reicht dafür NICHT (sie müsste Skalen r·e^{−W} ≈ 1e-4 an den Ecken
  auflösen, Grad ~10⁴); nötig ist die Fourier-/Gabor-Basis im (begradigten) Streifen — dort ist das Komplement exakt „hohe Frequenz längs des
  Streifens“, randkonzentriert, mit Innenspuren ≤ e^{−|σ|·Abstand} (explizit). Dafür braucht es die konforme Begradigung der oberen Ω-Kante in
  ζ = log M (Theodorsen-Iteration mit a-posteriori-Schranke) — rechnerisch gut machbar, für Kneser aus f_K allein (kein slog nötig).
Ehrlich: B3' ist ein sauberes Lemma (nächste Etappe); B5' ist die eigentliche Arbeit (mehrtägig). **[Überholt durch §30 (02:50): die
Konturverschiebung auf Mittellinien macht B3'–B5' überflüssig.]**


## 30. Satz O11 — (FA) durch Konturverschiebung auf Mittellinien: c_∞ ≥ 0.54 (Konstanten numerisch, Strenge-Rechnung folgt)

**Die Idee in einem Satz.** X(δ) = Re⟨Kδ, δ⟩_p paart δ in einem Randpunkt mit δ in einem Innenpunkt (Γin: x ↔ f(x); Γout: x ↔ f(x) auf
dem Linsendeckel). Schreibt man conj(δ) als analytische Funktion der gespiegelten Variablen und verschiebt die Kontur bis zur „Mitte“
zwischen den beiden Punkten, so wird aus der Paarung fast ein Betragsquadrat |u(z)|² mit positivem Gewicht — der Rest ist eine kleine
reelle Verschiebung τ (die loxodromische Drehung, 0.159 an den Enden) und wird durch eine POSITIVE Form abgeschätzt. Das ist die
Verallgemeinerung der Drei-Geraden-Identität ⟨u(·+iβ), u⟩ = ‖u(·+iβ/2)‖² auf variable Verschiebungen.

**Koordinaten (o26).** ζ = log M(x) − iθ₀, M = (x − L)/(x − L̄): D = {0 < Im ζ < π} exakt, Γin = {Im ζ = 0}, Γout = {Im ζ = π};
f ↔ g(ζ) = log M(f(x)) − iθ₀; Ω̃ := D ∪ g(D) = {0 < Im ζ < h(Re ζ)}, h ∈ [π + 0.585, π + 0.669] (Linse in der Mitte dünner).
δ = u·J^{−1/2} (J = dx/dζ; H²-isometrisch), Kδ ↔ R·u∘g mit R = (g'/f')^{1/2}/m. Dann (Normierung 1/(2πρ) weggelassen)
  ‖Bδ‖²_p = ∫ p_b|u|² dξ (Boden) + ∫ p_t|R_t|²|u∘g|² dξ (über Γout = Linsendeckel),
  X = Re∫ p_bR_b·u(g(ξ))·conj(u(ξ)) dξ + Re∫ p_tR̄_t·u(ξ+iπ)·conj(u(g(ξ+iπ))) dξ.
**Mittellinien (o29).** Unten: conj(u(ξ)) = ũ(ξ), ũ(ζ) := conj(u(ζ̄)) analytisch; p(ζ) := m(x(ζ))·conj(m(x(ζ̄))); Kontur ζ = c(ξ) = ξ − iε(ξ)
mit Im g(c) = ε („symmetrischer Punkt“). Dann ist g(c) = z + τ, z := ξ + iε, τ := Re g(c) − ξ **reell**, und
  X_b = Re∫ W_b(ξ)·u(z + τ)·conj(u(z)) dξ,  W_b = p(c)R(c)c'(ξ).
Oben analog (Kontur nach oben in die Linse, ε mit Im g(ξ + i(π − ε)) = π + ε): X_t = Re∫ W_t(ξ)·u(z)·conj(u(z + τ)) dξ, z = ξ + i(π + ε).
f_K ausserhalb D: jenseits Γin über f(x) = log f(eˣ) (eˣ ∈ D), in der Linse über f(x) = e^{f⁻¹(x)}.
**Probe (o31):** X direkt auf ∂D = X über die Mittellinien für zufällige Polynome bis auf 1e-5 (0.85098/0.85098, 5.16369/5.16370, …) ✓.

**Satz O11 (bewiesen unter den unten genannten, rechnerisch prüfbaren Bedingungen).** Ist Re W_b, Re W_t > 0 längs der Kurven, so gilt
für alle Polynome δ:  X(δ) ≥ −D(u) ≥ −Λ‖Bδ‖²_p  mit  D(u) := Σ_Kurven ∫ (|W|²/(4 Re W))·|u(z + τ) − u(z)|² dξ,
also (FA_{1−Λ}) und c_∞ ≥ c(1 − Λ) (Satz O9).
*Beweis.* Punktweise Re(W·u(z+τ)ū(z)) = Re W|u|² + Re(W·Δu·ū) ≥ Re W|u|² − |W||Δu||u| ≥ −|W|²|Δu|²/(4 Re W). Dann
|Δu|² ≤ τ²∫₀¹|u'(z + rτ)|²dr und u'(z) = (1/2πi)∮_{∂Ω̃} u(ζ)dζ/(ζ − z)² (Cauchy; ∂Ω̃ = Boden ∪ Linsendeckel, genau die Träger der B-Norm).
Schur-Test für den Kern k(ξ, ζ) = √ρ(ξ)/(2π|ζ − z_r(ξ)|²√w(ζ)), ρ := |W|²τ²/(4 Re W), w = B-Gewicht (Boden p_b, Deckel p_t|R_t|²/|g'|):
Λ ≤ sup_r [sup_ξ ∫k|dζ|]·[sup_ζ ∫k dξ]. ∎
*Bedingungen (für Polynome δ genügt Analytizität von f_K, f_K' ≠ 0 auf den überstrichenen Gebieten — jenseits Γin bis Tiefe ε ≤ 0.585 und in der
Linse — sowie Abklingen an den Enden, das |J|^{1/2} ~ e^{−|ξ|/2} liefert).*
**Zahlen (Kneser, o29/o32, Gleitkomma):**
| Grösse | unten (Γin) | oben (Linse) |
|---|---|---|
| ε (Höhe der Mittellinie) | 0.334 (Enden) … 0.585 (Mitte) | 0.334 … 0.307 |
| τ (Restdrehung) | ±0.159 (Enden), max 0.27 bei ξ ≈ ±1, 0 in der Mitte | ±0.159 … 0 |
| min Re W | 0.563 | 0.582 |
| max ρ = |W|²τ²/(4 Re W) | 0.0395 | 0.0127 |
Schur: **Λ ≤ 0.086** (sup über r ∈ {0, ½, 1}: 0.078/0.082/0.086) ⇒ **λ ≥ 0.914 ⇒ c_∞ ≥ c(0.914) = 0.540.**
Galerkin (o30, Wellenpakete, drei Basen, Integration ξ ∈ [−40, 40]): Λ_wahr ≈ 0.0296 ⇒ λ ≈ 0.970, c ≈ 0.561. (Mit Integration nur über
[−12, 12] entstanden Scheinwerte 0.14/0.28 — abgeschnittene Deckelnorm; behoben.) Exakter Endwert λ_L = 0.98988 (O10), wahres c_∞ = 0.7319.

**Was das ersetzt.** Plan B3'–B5' (§29: Ecke, IMS-Kleben, Mittelrechnung) ist überflüssig: die Konturverschiebung behandelt Ecken (ξ → ±∞,
dort konstante Modellwerte) und Mitte in EINEM Argument, ohne Lokalisierung, ohne Fredholm-Theorie, ohne Injektivität. Die einzige
Rechnung ist die Auswertung expliziter Funktionen längs zweier Kurven und zweier Kurvenintegrale (Schur).
**Für die Strenge (nächste Etappe):** (1) einfache explizite Kurven (stückweise linear an die numerischen angelehnt; τ darf dann komplex sein,
|Δu| ≤ |τ|·…); (2) arb-Einschlüsse von W, τ, ρ, w mit Ableitungsschranken zwischen Gitterpunkten; (3) Enden |ξ| ≥ ξ₁ über die lokale
Entwicklung an L (B3: Abweichungen O(e^{ξ})); (4) Schur-Integrale mit Fehlerschranken; (5) Analytizität von f_K auf den überstrichenen
Gebieten (eˣ ∈ D, f(eˣ) ≠ 0; Linse ⊂ f(D)); (6) Kneser-Koeffizienten mit zertifizierten Fehlern + Schwanzschranke (→ P).


## 31. Abgleich mit Strang M (Satz M6) und neue Aufteilung (03:05, Auftrag der Leitung)

- **Waagrechte Konturen (o33/o34).** Satz O11 gilt mit jeder Kontur; mit waagrechten Geraden (unten Im = −0.40, oben π + 0.334) wird τ komplex,
  Schranke unverändert (|Δu|² ≤ |τ|²∫₀¹|u'(z + rτ)|²dr). Kneser: min Re W 0.574/0.582, Schur Λ ≤ 0.092 ⇒ λ ≥ 0.908, **c_∞ ≥ 0.538**.
- **Prüfstein Ms radiales Modell (o34, Gleitkomma):** min Re W 0.66/0.74, Λ ≤ 0.345 ⇒ λ ≥ 0.655, **c_∞ ≥ 0.408** — Ms strenger Satz M6 liefert 0.4249
  (wahr 0.6905). Kubisch: Λ ≤ 0.138, c_∞ ≥ 0.520. Die beiden Wege (Ms 𝒯 mit IMS + Nyström; Os Konturverschiebung + Schur) sind unabhängig und
  stimmen in der Grössenordnung überein. Ms 𝒯-Weg vermeidet die Linse ganz (Γout-Daten mit h = f⁻¹ zurückgezogen, nur Werte in D̄).
- **Aufteilung (Vorschlag, BRETT 03:05):** (a″)_∞ für Kneser: Hauptweg M (𝒯), O11 als Gegenprobe (nur noch 1D-Einschlüsse + Schur-Integrale nötig,
  sobald P zertifizierte f_K-Werte liefert). **O neu: Kante** — Rechteck-Abschnitte mit kleinem d = αN (für Xs Θ^c-Seite nötig): Eckmodell im
  Doppel-Skalenlimes x = L·e^{−w/N} (D → Re w > 0, P_N → Laplace-Spektrum [0, 1], f → w ↦ μw, T_{(1+α)N} → Spektrum ≤ 1 + α), kritisches α* und Kleben
  mit dem Bulk. Antwort an W zu e^{Nφ}-Gewichten: Restverschiebung τ kostet exp(N|∇φ||τ|); gleichmässig nur bis r ~ N^{−0.9}.


---

# Runde 7 (ab 03:07, Auftrag der Leitung): Kante — Rechteck-Abschnitte mit kleinem d, Innenschicht an L

## 32. Abstimmung und Befunde zur Kante (o35–o37; kubisch exakt, Kneser mit Strang-A-Daten)

**Aufteilung (BRETT 03:05/03:32/03:41):** X behandelt (K_d) von der Rechteck-Seite (Gram-Bedingung der Überlaufprofile), W die äussere Lösung
(Kegel, Raten, Gonchar–Rakhmanov–Stahl), **O die Innenschicht an L** (Skala 1/N, Amplitudenebene) — Ws Programm §22.6(C) und die j/N-Struktur für X.

**32.1 Die Eckmode ist asymptotisch ganz eine Innenschicht (numerisch, kubisch N = 200…600).** Für die normierte Eckmode v (grösster
Singulärwert von J_N⁻¹, Eingang; v aus der arb-Lösung J_N v = s·u, sonst Gleitkomma-Müll) und y := N(θ_L − φ):
| N | 200 | 300 | 400 | 500 | 600 |
|---|---|---|---|---|---|
| s₁ = σ₁(J_N⁻¹) | 1.852 | 3.698 | 9.100 | 24.77 | 70.71 |
| \|v(L)\|/√N | 0.296 | 0.368 | 0.398 | 0.410 | 0.416 |
| Masse in \|y\| ≤ 20 (beide Ecken) | 0.51 | 0.73 | 0.83 | 0.88 | 0.90 |
| Masse in \|y\| ≤ 50 | 0.66 | 0.86 | 0.94 | 0.97 | 0.98 |
⇒ v(L e^{−iy/N}) ≈ √N·V(y) mit einem Grenzprofil V (Γout-Seite schon konvergent: |V(10)|/|V(0)| = 0.22 → 0.39, |V(40)|/|V(0)| ≈ 6e-4;
Γin-Seite fällt wie |y|^{−1.8…−1.9}). **Folgerung: der Eckvorfaktor ist N^{1/2}** — eine Funktion der Breite 1/N mit N-unabhängigem Massenanteil;
das passt zu Ws gemessenem κ(0) ≈ 0.51 (Normierungen noch abzugleichen: W misst y = σv).
**32.2 Universalität (numerisch):** Kneser und kubisches Modell (gleiches L, μ = √L, sonst verschieden) stimmen in der Innenschicht auf 3 Stellen
überein: N = 200/300: |v(L)|/√N = 0.2950/0.3668 (kubisch 0.2964/0.3678), Profil y = −50/+10/+20: 0.321/0.221/0.0158 (0.318/0.222/0.0159),
s₁ = 1.841/3.649 (1.852/3.698), Kantenwerte identisch. ⇒ Die Eckmode hängt in führender Ordnung nur von μ ab (lokales Objekt).
Aber: im affinen Modell f = L + μ(x − L) gibt es keinen Abschneide-Effekt (Komposition erhält den Grad, P_N invariant, Rechteck = quadratisch). Die Innenschicht ist also
nicht das affine Innenproblem allein; die Abschneidung wirkt über die äussere Lösung (Ws Kegel, ebenfalls nur von μ abhängig) — die Eckmode
ist ein Anpassungsproblem innen/aussen. Der Kern des affinen Innenoperators μV(w) + V(μw) = 0 sind genau meine Mellin-Moden w^{s_k}
(s_k = 1 + (2k+1)πi/log μ, Re s_k = 5.45 + 8.90k); das gemessene Γin-Abklingen |y|^{−1.9} ist KEINE reine Mellin-Mode.
**32.3 Kante in Koeffizienten (Ratenebene, universell):** log10|[J_∞v]_k|/N an der letzten behaltenen Zeile k = N: −0.0349/−0.0349/−0.0354/
−0.0357/−0.0361; an der ersten Überlaufzeile k = N+1: +0.0309/+0.0308/+0.0304/+0.0301/+0.0299; oberster Eingangskoeffizient |v_N|: −0.0383/
−0.0387/−0.0397/−0.0394/−0.0397 (N = 200…600). Die Eckmode löscht also N+1 Zeilen „natürlicher Grösse“ 10^{+0.03N} bis auf 10^{−0.035N};
ab Zeile N+1 keine Auslöschung mehr. Überlaufrate Φ(j/N) = log10|[J_∞v]_{N+j}|/N konvergiert: 0.042 (0.01), 0.051 (0.02), 0.071 (0.05).
Eine einzelne Eckmode wird also schon von Zeile N+1 erfasst; (K_d) ist die Frage, ob r Eckmoden + Volumenvektoren die ersten d Überlaufzeilen
gemeinsam auslöschen können (Xs lineare Unabhängigkeit) — die Antwort liegt im Vorfaktor (Amplitudenebene) der Überlaufprofile.
**32.4 Warum kein naives Innenproblem:** Im Grenzwert N → ∞ ist V bandbegrenzt (Spektrum ⊂ [0, 1]), V(μ·) wächst auf der Γout-Seite wie
e^{Im μ·y}; die Projektion P_{≤1} ist dafür nicht erklärt (adjungiert ebenso). Das Innenproblem braucht also die Anpassungsbedingung aus dem
Kegel (W) — es ist die lokale Parametrix eines Riemann–Hilbert-/Orthogonalitätsproblems an L, kein selbständiges Wiener–Hopf-Problem.
**Stand (ehrlich):** Ein Beweis der Kanten-Aussage (K_d) ist nicht in Reichweite dieser Runde. Belastbar sind die Befunde 32.1–32.3 (Grenzschicht,
Vorfaktor N^{1/2}, Universalität, Kantenraten) — Eingänge für Ws lokale Parametrix und Xs Gram-Bedingung.
Dateien: o35_ueberlauf.py, o36_innenschicht.py (+ o36_kubisch.log, o36_kneser.log, o36_*.npz), o37_profilvergleich.py.

**32.5 Schnittstelle zu W (BRETT 03:49), übernommen:** Ws Aussenteil braucht von O genau die **Kreisbedingung (KB)**:
max über |t − L| = A·N^{−1/Re κ₁} (= A·N^{−0.894}) von |Y_i| und |T_i| ≤ C_E·N^{κ_E} (gemessen κ_E ≈ 0.9; Y_i, T_i analytisch um L ⇒ dieselbe Schranke
auf der ganzen Scheibe). Anschluss bei |w| → ∞: Wachstum ≤ exp(|c₁|N^{1−κ₁}|Lw|^{Re κ₁}) mit Naht und zweitem Knoten 161°.
*Was ein Beweis von (KB) bräuchte (Einschätzung O):* die Innenschicht 32.1 als Lösung einer Anpassungsaufgabe — innen das affine Modell
(Mellin-Kern w^{s_k}), aussen Ws Kegel — mit Fehlerschranken; die triviale Bernstein-Schranke |v(z)| ≤ √(N+1)(|z|/ρ)^N gibt auf dem (KB)-Kreis nur
exp(A·N^{0.106}/ρ), also nicht polynomial, und sieht den Faktor σ_i nicht. Das ist ein Forschungsprogramm (lokale Parametrix), keine Etappe dieser Runde.

---

# Runde 8 (ab 04:03): lokale Parametrix an L — Kreisbedingung (KB), festes Grenzobjekt

## 33. (KB) mit κ_E = 1/2 für alle Eckmoden zugleich; das feste Grenzobjekt Ψ und warum es noch nicht zertifiziert ist

**33.1 Objekte und eine SVD-freie Majorante (bewiesen, elementar).** J̃^{−1} = U S Vᵀ, u_i = U[:, i] (Vorwärts-Eckmode, Einheitsvektor),
y_i = σ_i V[:, i] = J̃^{−T}u_i (Ws Normierung J̃ᵀy_i = u_i). Mit ẽ_x := ((ρ/x)^k)_k (abgeschnittener Cauchy-Kern bei x),
c = T_N J_∞(1/(η − ·)), δ_η = T_N 1/(η − ·):
  Y_u(x) = uᵀ J̃^{−1} ẽ_x,   T_u(η) = uᵀ(J̃^{−1}c̃_η − δ̃_η)   (für u = u_i genau Ws Y_i, T_i; Abgleich mit w46 bei N = 200: 1.230/2.354/4.796 gegen 1.22/2.33/4.73,
  ln|T₁(L)| = 0.571 gegen 0.57).
Cauchy–Schwarz gibt für JEDEN Einheitsvektor u: |Y_u(x)| ≤ G(x) := ‖J̃^{−1}ẽ_x‖, |T_u(η)| ≤ G^T(η) := ‖J̃^{−1}c̃_η − δ̃_η‖.
G ist die Norm der Abschnittslösung J_N δ = T_N[x/(x − ·)] — eine einzige lineare Gleichung, ohne SVD, ohne Modenzahl r(N), und sie deckt alle
Eckmoden und ihre Linearkombinationen ab. Ausserdem exakt: Y_i(x) = σ_i²·(J̃u_i)(ρ²/x), d. h. (KB) sagt, dass das Residuum der Vorwärts-Eckmode
bei L um den Faktor σ_i² kleiner ist, als seine ℓ²-Norm σ_i^{−1} erwarten liesse.
**Lemma O12 (Eckformel bei L; W §17(2), hier mit voller Reihenlänge geprüft: Rest 1e-14 absolut, N = 300).** T(ζ) = Y(ζ)m(ζ)/ζ
+ Σ_{f(x)=ζ} Y(x)/(x f'(x)) + B_ζ − 𝒰(ζ), 𝒰(ζ) = Σ u_kρ^kζ^{−k−1}. In ζ = L (f(L) = L, m(L) = f'(L) = μ):
**Y(L)·(μ + 1/μ)/L = T(L) + 𝒰(L) − B_L − (ferne Urbilder)**, und |𝒰(ζ)| ≤ √(N+1)/ρ für |ζ| ≥ ρ (bewiesen). Der Eckwert von Y ist also bis auf
√N-Terme der von T (Lemma W9 in Fixpunktform).
*Werkzeug-Fehler (behoben):* python-flint kappt Reihen bei ctx.cap = 10 Termen; ohne ctx.cap = n + 10 war c_η falsch (T-Werte der ersten
o38-Fassung verworfen und neu gerechnet).

**33.2 Messung (o38, o40; arb-Inverse, Gleitkomma-Auswertung; alle Werte MINUS ½·ln N).** Kreise |t − L| = A·ρ·N^{−1/Re κ₁}, Ws Konvention
(innen Y, aussen T); KB₁ = Eckmode 1, KB_r = Wurzel der Quadratsumme über die r Moden mit σ > 1.4 (Ws Grösse), KB_G = Majorante (alle u).

| N | ln σ₁ | Ecke: Y₁ / T₁ / G / G^T | A = 0.5: KB₁ / KB_G | A = 1: KB₁ / KB_r / KB_G | A = 2: KB₁ / KB_G | A = 4: KB₁ / KB_G |
|---|---|---|---|---|---|---|
| 200 | 0.62 | −2.47 / −2.08 / −0.85 / −0.82 | −1.96 / −0.54 | −1.42 / −1.42 / **−0.17** | −0.30 / 0.73 | 2.15 / 2.96 |
| 300 | 1.31 | −2.30 / −1.86 / −0.88 / −0.78 | −1.80 / −0.57 | −1.27 / −0.98 / −0.20 | −0.16 / 0.70 | 2.24 / 2.90 |
| 400 | 2.21 | −2.24 / −1.79 / −0.89 / −0.76 | −1.74 / −0.58 | −1.21 / −0.89 / −0.22 | −0.11 / 0.67 | 2.26 / 2.86 |
| 600 | 4.26 | −2.21 / −1.74 / −0.92 / −0.73 | −1.71 / −0.61 | −1.19 / −0.84 / −0.24 | −0.10 / 0.64 | 2.26 / 2.81 |
| 800 | 6.43 | −2.20 / −1.73 / −0.93 / −0.71 | −1.70 / −0.62 | −1.18 / −0.79 / −0.26 | −0.10 / 0.62 | 2.25 / 2.78 |
| 1000 | 8.65 | −2.20 / −1.73 / −0.94 / −0.70 | −1.70 / −0.64 | −1.18 / −0.75 / −0.27 | −0.10 / 0.61 | 2.24 / 2.76 |

Kneser (Ps zertifizierte Koeffizienten, Modell `kneser_z`, N = 200/300/400/600): **jede Zelle gleich der kubischen auf ≤ 0.01.**
**Ergebnis (numerisch, N = 200…1000): (KB) gilt mit κ_E = 1/2 — für jede Eckmode einzeln und, über die Majorante, für alle zugleich:**
  max_{|t−L| = Aρ N^{−0.894}} (|Y_u| innen, |T_u| aussen) ≤ C_E(A)·√N,  C_E(0.5) = 0.58, **C_E(1) = 0.84**, C_E(2) = 2.1, C_E(4) = 19
(Höchstwerte über die gemessenen N; die Majoranten-Spalten FALLEN mit N, der Höchstwert liegt bei N = 200). Für die Einzelmoden allein:
C₁(1) = 0.31, C₁(2) = 0.91. Ws κ_E ≈ 0.9 kam von der Quadratsumme über r(N) ~ N/250 Moden (KB_r steigt, weil r steigt); einzeln ist es N^{1/2}.
Eckwerte: |Y₁(L)| → 0.11√N, |T₁(L)| → 0.18√N, G(L) → 0.39√N; Konvergenz an der Ecke schnell (Abstände zum Grenzwert 0.03/0.01/0.002 bei
N = 600/800/1000), auf dem A = 4-Kreis langsam (Drift ≈ 0.2 von N = 200 bis 1000, Ordnung N^{−a} mit a ≈ 0.1–0.5).

**33.3 Das feste Grenzobjekt (numerisch).** Auf der Skala z = (t − L)/(ρN^{−1/Re κ₁}) kollabieren die Profile Ψ_Y(z) = Y_u(t)/Y_u(L),
Ψ_T(z) = T_u(t)/T_u(L) (o38b): ln|Ψ_Y| auf |z| = 1 je ψ (Grad, von der äusseren Normalen) für N = 200 → 800: ψ = 180°: 1.03 → 0.99,
ψ = 0°: −0.96 → −0.93, Vorzeichenwechsel bei ψ ≈ −100° und +80°; auf |z| = 2 und 4 ebenso (4.50 → 4.30 bei ψ = 180°, A = 4).
Taylor-Koeffizienten um z = 0 (o39, z in Einheiten ρN^{−0.894}e^{i arg L}):
  Ψ_Y = 1 − (0.959 + 0.223i)z + (0.466 + 0.235i)z² − (0.147 + 0.126i)z³ + (0.033 + 0.046i)z⁴ − …  (N = 800; N = 400 auf ≤ 0.015 gleich)
  Ψ_T = 1 − (1.300 − 0.042i)z + (0.839 − 0.045i)z² − (0.359 − 0.022i)z³ + (0.115 − 0.006i)z⁴ − …
Kneser = kubisch auf 3 Stellen (N = 400/600). Ψ_Y ≈ e^{−az}, a ≈ 0.96 + 0.22i: Y wird nahe L von Graden k ~ a·N^{0.894} getragen
(y_k(ρ/L)^k addiert sich für k ≲ 0.7N gleichphasig, die grossen Koeffizienten bei k ≈ 0.95N löschen sich in L aus; N = 400).
**Struktur (numerisch):** (i) alle Eckmoden haben dasselbe Ψ_Y (Moden 1/2 bis 0.05, Mode 3 bei N = 800/1000 bis 0.3 auf A = 4) — der lokale
adjungierte Eckraum ist eindimensional, nur die Amplitude (Eckwert) ist global. (ii) Die reellen Eckmoden kommen paarweise: u₁, u₂ haben in L
gleichen Betrag (0.416/0.411·√N bei N = 600) und Phasenabstand 87.7°/89.6°/90.0° (N = 400/500/600); c = ū₂(L)u₁ − ū₁(L)u₂ ist eine
**komplexe Eckmode, die nur in L sitzt** (bei L̄ 1.9 %/1.0 %/0.6 % des Wertes bei L). (iii) Die Quelle 𝒰 ist auf dieser Skala NICHT klein:
|𝒰|/|Y| = 2.72/2.75 auf |z| = 1 (N = 400/800) — die Vorwärtsmode trägt dieselbe √N-Skalierung bis zur Anpassungsskala.

**33.4 Formulierung des Grenzproblems (Skizze).** Mit y(z) := Y(L + ρN^{−1/κ₁}z)/√N, τ, 𝔲 entsprechend (komplexe Potenz N^{1/κ₁}, damit
auch die langsame Phase e^{0.025i·ln N} verschwindet) gilt nach Lemma O12 exakt bis auf O(N^{−0.894}) und die (gemessen 1e-7 kleinen) fernen Terme:
  **(P1)** τ(z) = (μ/L)·[y(z) + μ^{−2} y(z/μ)] − 𝔲(z),  y, τ, 𝔲 ganz;
  **(P2)** 𝔲(z) = conj(ũ(−z̄))/L mit der Vorwärtsmode ũ (Spiegelung am Kreis), und ũ erfüllt μũ(z) + ũ(μz) = õ(z) mit dem Überlauf õ;
  **(P3)** Wachstum für |z| → ∞ aus Ws Kegel: ln|y| ~ Re(c₁z^{κ₁}) im Scheibensektor, τ rezessiv (Nulllinie = Naht Z_T) im Aussensektor,
       Deckverschiebung 2πi − log μ, zweiter Knoten 161°.
(P1)–(P3) hängt nur von μ und der Kegelkonstanten c₁ (Skalierung) ab; die Messungen 33.2–33.3 zeigen genau diese Universalität. Mellin in z
macht (P1) diagonal: τ̂(s) = (μ/L)(1 + μ^{s−2})ŷ(s) − 𝔲̂(s), und 1 + μ^{s−2} = 2πμ^{(s−2)/2}/[Γ(½ + iα(s−2))·Γ(½ − iα(s−2))], α = log μ/(2π):
(P3) wird damit ein skalares Faktorisierungsproblem vom Malyuzhinets-Typ (Keil mit Klebung) — der natürliche Weg zu einer EXPLIZITEN lokalen
Parametrix. Nicht ausgeführt.

**33.5 Stand (ehrlich).** *Bewiesen:* Majoranten G, G^T (33.1), Lemma O12 samt |𝒰| ≤ √(N+1)/ρ, die Identität Y_i = σ_i²(J̃u_i)(ρ²/·).
*Numerisch, belastbar:* (KB) mit κ_E = 1/2 und den Konstanten aus 33.2 für N = 200…1000, für alle Eckmoden; Universalität (Kneser = kubisch);
das feste Objekt Ψ (Kollaps, Taylor-Koeffizienten), Eindimensionalität des lokalen Eckraums, komplexe Eckmode.
*Nicht erreicht — Zertifikat des Grenzproblems:* (a) die Skalen trennen sich nur wie N^{0.106} (1/N gegen N^{−0.894}), und die Quelle 𝔲 ist auf
der Anpassungsskala so gross wie y (33.3(iii)) — ein quellenfreies Grenzproblem wäre für erreichbare N falsch; (b) (P3) braucht Ws
äussere Daten (c₁, Sektoren), die heute numerisch sind; (c) Existenz und Eindeutigkeit von (P1)–(P3) sind nicht bewiesen. Die endlichen N selbst
sind billig zu zertifizieren (Sekunden bis Minuten je N), aber (KB) für ALLE N ≥ N₀ braucht den asymptotischen Schritt, und der erbt — wie Ws
Eckzone (§24.2) — einen globalen Eingang: die Amplitude der Eckmode (Eckwert) ist nicht lokal bestimmt. Einschätzung: (KB) gehört als Teil von
(H₀) in den Hauptsatz; diese Runde liefert die scharfe Form (κ_E = 1/2, C_E(1) = 0.84, für alle u) und die Gestalt der lokalen Parametrix.

**33.6 (K_d) — wie es aus der Parametrix folgen würde (Notiz, nicht verfolgt).** Xs Gram-Bedingung verlangt lineare Unabhängigkeit der
Überlaufprofile j ↦ [J_∞u_i]_{N+j}, j ≤ εN. Nach 33.3(i)/(ii) unterscheidet die Ecke die Moden NICHT in der Gestalt: ein Paar reeller Moden
ist Real-/Imaginärteil EINER komplexen, nur in L sitzenden Mode c (und c̄ in L̄). Für das Paar folgt (K_d), sobald die Parametrix das Überlaufprofil
von c in der Form [J_∞c]_{N+j} = A(j/N)·e^{ijθ_c}(1 + o(1)) mit θ_c ≠ 0 mod π liefert (die Phase dreht mit j; Real- und Imaginärteil trennen
sich auf j ~ 1/θ_c — passt zu Xs „Trennung auf j ≈ 0.07N“, wenn θ_c ~ 14/N). Zwischen verschiedenen komplexen Moden (r = 4 ab N ≈ 800)
entscheidet nur der Eckwert bzw. die globale Amplitude — dort ist (K_d) keine Eckaussage mehr.
Dateien: o38_adjungiert.py (+ o38_kubisch/kneser_z .json/.log), o38b_kollaps.py, o39_zskala.py (+ o39_*.npz/.log), o40_kb_tabelle.py; o_kern.py
um `kneser_z` (Ps Zertifikatsdatei) erweitert.

## 34. Abschluss (Team wird angehalten, 11.10. 04:48) und Wiederaufnahme-Notiz

**Gesamtstand Strang O.**
- *Bewiesen:* Satz O1 (Faktorisierung), O7/O7b (Index 0, Pfadbedingungen in arb), Satz O9 (Reduktion von c_∞ auf (FA)), Lemma O10
  (Eckwert λ_L = 0.98988, arb), Lemma O12 (Eckformel in L, 33.1), SVD-freie Majoranten G, G^T (33.1), Y_i = σ_i²(J̃u_i)(ρ²/·).
- *Numerisch:* Satz O11 (Konturverschiebung, c_∞ ≥ 0.54 Kneser, 0.408 radial gegen Ms strenge 0.4249), Innenschicht der Vorwärts-Eckmode
  (Breite 1/N, |u(L)| → 0.42√N, §32), **(KB) mit κ_E = 1/2, C_E(1) = 0.84 für alle Eckmoden, N = 200…1000, Kneser = kubisch** (§33.2),
  festes Objekt Ψ und komplexe Eckmode (§33.3).
- *Offen:* (KB) für alle N ≥ N₀ (braucht globalen Eingang, gehört zu (H₀)); Zertifikat des Grenzproblems (P1)–(P3); (K_d); O11 streng.

**Schnittstelle zu W (Kreisbedingung §32.5 → Schweiss-RH-Aufgabe).** W übernimmt (KB) als innere Randbedingung seines Kegels (Maximumprinzip
auf |t − L| ≤ r_N, W §23.2(iii)) und klebt bei r ~ A·N^{−0.894} an die äussere Lösung. Als Riemann–Hilbert-Aufgabe gelesen: äussere Parametrix =
Ws Kegel (Ratenebene, c₁, Naht Z_T, Deckverschiebung 2πi − log μ), lokale Parametrix in L = Lösung von (P1)–(P3) (§33.4), Schweissnaht = Kreis
|z| = A; Fehlerordnung der Anpassung numerisch N^{−a}, a ≈ 0.1–0.5 auf A = 4, an der Ecke schnell. Was W von O übernehmen kann, sofort:
die Zahlen aus §33.2 (Majorante, gilt für jedes Einheits-u) statt der Quadratsumme über r(N).

**O11 als spätere Gegenprobe.** Mit Ps zertifizierten f_K-Koeffizienten (gemeinsam/kneser_halbexp_koeff_zert_640.json, Satz P-V) braucht O11 nur
1D-Einschlüsse längs zweier Geraden (o29/o33: Mittellinien bzw. waagrechte Konturen) und zwei Schur-Integrale (o32); Werkzeug: Ms
mk_kneser_arb.py (Klasse KneserArb). Ergebnis wäre eine strenge zweite Schranke c_∞ ≥ ~0.53 unabhängig von Ms 𝒯-Weg.

**Wiederaufnahme — nächster Schritt.** (1) (P1)–(P3) per Mellin-Faktorisierung (Γ(½ ± iα(s−2)), §33.4) lösen und mit den Taylor-Koeffizienten
von Ψ_Y, Ψ_T (o39_*.npz, Feld `tay`) vergleichen; dazu von W c₁ und die Sektorwinkel holen. (2) Falls das trägt: Fehlerordnung N^{−a} der
Anpassung abschätzen und (KB) als „(H₀)-Teil mit expliziter Parametrix“ formulieren. (3) Unabhängig: O11 streng (siehe oben).
Einstieg: dieser Bericht §32–§34; Daten o38_kubisch.json / o38_kneser_z.json (Kreise, Ecke, alle Moden), o39_<modell>_<N>.npz (Kreise bis |z| = 8,
Taylor, D6-Probe), o36_*.npz (Vorwärts-Eckmoden, Innenschicht); Tabellen neu erzeugen mit `python o40_kb_tabelle.py kubisch kneser_z` und
`python o38b_kollaps.py kubisch Y`. Interpreter C:/Users/janis/miniconda3/python.exe; Läufe bis N = 1000 dauern ≈ 6 min (arb-Inverse).
Es laufen keine Prozesse von Strang O.
