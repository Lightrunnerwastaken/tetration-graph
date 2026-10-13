# Hauptsatz in der K_{r₀}-Fassung — Polynom-Wurzel von exp → Kneser

**Fassung 3**, Stand 12.10.2026 12:49 (Strang K); §2.6 korrigiert 13:07 ((KK) → (KK′), k_kern §39). Sie ersetzt Fassung 2 (11:05, unverändert in Anhang B) und Fassung 1 (10:23, Anhang A).

Grundlagen: ARCHITEKTUR.md §4 und §9 (nichtlineare Kette), k_kern/BERICHT.md §24–§38 (lineare Aussage), KERNPROBLEM.md §4 (der gemeinsame Kern (KK)),
v_pruefstand/BERICHT §1.11, §1.12, §3.19–§3.23 (Vs Gegenlesen und Messungen).

Statuswörter:
- **V**: von einem zweiten Strang Schritt für Schritt geprüft.
- **S**: rigoros gerechnet (arb), vom Autor selbst geprüft.
- **B**: vom Autor bewiesen, noch nicht gegengelesen.
- **num**: numerisch, für endlich viele N.
- **offen**: weder bewiesen noch gerechnet.

Was sich gegenüber Fassung 2 ändert:
- **(P) entfällt** (Vs Anmerkung C). (G-P) wird absolut formuliert, als Schranke für |Δu_n(z)| auf K_{r₀}. Damit braucht man keine Unterschranke für |P_n|.
  Mit Vs Messung E4′ (V §3.22) ist die Summe Σ_n Q_n beschränkt: Das Budget ist γ = 0 statt 0.89.
- **E ⊃ K_{r₀}** statt E ⊃ D̄ (Vs Anmerkung A). Die Aussenöffnung von E in L ist **131.7°**, nicht 99° (Anmerkung B). Für die Konstanten bei N₀ = 400 ist r₁ ≈ 0.2 besser als 0.05.
- **(S_E) hat jetzt den Status „gelesener Satz + numerisch geprüfte Geometrie“.** Satz K22 ruht auf dem gelesenen Satz 2 von Andrievskii–Pritsker–Varga. Die Geometrie von E ist numerisch geprüft (kubisch k54, f_K v55).
  s_E = 5.417·0.732 = 3.96, gemessen 3.9–4.1.
- **Familie:** F3 in der Form „sup_𝒞|P − f_K| ≤ δ₀“ ist **falsch** (V §3.20: auf der Linse e^{+0.015N}). Sie wird ersetzt durch **(Fam-koh)^{r₀}** für Linse UND Σ₁.
  Die Hilfsannahme (EP⁺) ist gestrichen, weil Fenster 1 primal ein Wachstumsfenster ist (k_kern §37.1).
- **Der gemeinsame Kern** ist benannt (§2.6; seit 13:07 als (KK′), weil die Steg-Lücke von (KK) ein Messartefakt war). (EP) und (Fam-koh)^{r₀} bleiben benannte Hypothesen.
- Fassung 1, Anhang A.4: Die Schlitzkonstante für K′ ist **8.2, nicht 16.3** (V §1.11: die „2·“ in ARCH §9.7 (1) ist ein Schreibfehler).

## 1. Aussage

**Bezeichnungen.**
- L = 0.31813 + 1.33724i (e^L = L), ρ = |L|, D = {|x| < ρ}, μ = √L, a = 1 + iL/(4π), r_n = ρn^{−Re a} (Kegelskala).
- f_K ist Knesers halbe e-Funktion. p_N sind die Polynome des Autors (T_N[p_N∘p_N − eˣ] = 0, Zweig der σ-Homotopie ab p₁).
- **K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄))** mit festem r₀ ∈ (r_*, 0.1]. r_* ist der Abstand der Zone, in der der primale Kopf auf der Γout-Seite zu wachsen beginnt
  (Kneser ~ 1e-6 nach der Profilrechnung, Schätzung; Modell M4 ≈ 0.13, Vs T1).
- Köpfe: c_n = J̃_n⁻¹e_n, ℓ_n = J̃_n^{−T}e_n, δ_n(z) = Σ c_{n,k}(z/ρ)^k, Λ_n(t) = Σ ℓ_{n,k}(ρ/t)^k.
  Kern K_n(z, t) = ε_zᵀJ̃_n⁻¹ẽ_t = t·(J_n⁻¹T_nk_t)(z); Zuwachs ΔK_n = K_n − K_{n−1} = δ_n(z)Λ_n(t)/ℓ_nn (Rang 1).
- Die zwei **Eck-Pol-Antworten**: u_n := J_n⁻¹T_nk_L = Π_n s (primal; s = J_∞⁻¹k_L hat einen Pol in L) und w_n(t) := K_n(L, t) (adjungiert).
  Es gilt ΔK_n(z, L) = L·Δu_n(z) und, wegen Rang 1, **ΔK_n(z, t) = L·Δu_n(z)·Λ_n(t)/Λ_n(L)**.
- **Kontur 𝒞** (fest, für alle n dieselbe, geschlossen um 0):
  - die Spiralen Σ₁, Σ₂ in L: θ(r) = u_A + 71.46° + 1.623°·ln(r/ρ), u_A ≈ 80° (ausserhalb D) bzw. ≈ −90° (in D), 0 < r ≤ r₁ (r₁ = 0.05; für die Konstanten besser 0.2);
  - die Linse |t| = ρ(1 + 0.03) über Γin, der Γout-Bogen |t| = ρ, Verbindungsbögen;
  - konjugiert bei L̄.
- **Kompaktum E** := Abschluss von int 𝒞 ∪ f_K(int 𝒞). Es gilt **E ⊃ K_{r₀} ∪ 𝒞 ∪ f_K(𝒞)** (nicht E ⊃ D̄). f_K ist auf E stetig und im Innern holomorph.
  E ist eine Quasischeibe mit Spiralecken der Öffnung Θ ≈ 132° in L, L̄ (num, k54/v55). Der Spiralschnitt von L aus läuft in der Lücke, die E frei lässt.

**Satz A (lineare Konvergenz bei f_K; Beweis §3, Status B, Logik V).** Es gelten (EP), (FA) und (S_E) aus §2. Dann gilt für alle N
  **sup_{K_{r₀}} |Π_N f_K − f_K| ≤ C·N^{γ − s_E},   γ := max(0, 1 + κ_G + κ_A)**; im Grenzfall 1 + κ_G + κ_A = 0 kommt ein Faktor ln N dazu.
Gemessen: κ_G ≈ −3.2 (r₀ = 0.05; R_P ~ n^{−3.1} mal |P_n| ~ n^{−0.106}), κ_A = 0, also γ = 0 und Rate N^{−3.96}.

**Hauptsatz (Status: bedingt).** Es gelten die Annahmen von Satz A für die ganze Familie (F1, F2, F3′, F4 in §2.4). Dann gilt für alle N ≥ N₀ = 400:
- p_N ist die einzige Lösung von T_N[P∘P − eˣ] = 0 in B_N^{r₀} := {P : sup_{K_{r₀}}|P − f_K| ≤ R_N}, R_N = K₀C′N^{γ−s′}, K₀ = 6, mit einem festen s′ ∈ (γ, s_E].
- **sup_{K_{r₀}} |p_N − f_K| ≤ C′·N^{γ − s′} → 0.**
- Für N < 400 sind die p_N streng eingeschlossen (Zertifikat, §3).

**Folgerungen** (B, je eine Zeile):
- (a) |a_k^{(N)} − a_k| ≤ (ρ − r₀)^{−k}·C′N^{γ−s′}, also a_k^{(N)} → a_k für jedes k (Cauchy auf |x| = ρ − r₀ ⊂ K_{r₀}).
- (b) p_N(1) → f_K(1) = sexp_K(½) = T_e(½), denn 1 ∈ K_{r₀} (|1 − L| = 1.50).
- (c) Die Identifikation mit Kneser ist eingebaut, weil die Kugel um f_K liegt.

**Nicht behauptet** ist die gleichmässige Konvergenz bis an L heran (Vs T1: im Modell M4 Plateau von Π_N U − U in einer festen Γout-Zone). Deshalb ist **r₀ > r_* nötig**.

## 2. Alle Annahmen

### 2.1 (EP) Eck-Pol-Hypothese: polynomiale Zuwächse der beiden Eck-Pol-Antworten
| Teil | Aussage (für alle n) | gemessen | Status |
|---|---|---|---|
| (G-P) | **G_n := sup_{z∈K_{r₀}} \|Δu_n(z)\| ≤ C n^{κ_G}** (absolut) | G_n = \|P_n\|·R_P/\|L\|: R_P(0.05) ~ n^{−3.1} bis n = 1800, \|P_n\| = 0.535·n^{−0.106}; Kneser = kubisch auf < 1 % (V §3.22) | num; offen. Rand 120°–180° bis 1800 noch linear fallend (vorübergehend, K20) |
| (A-lok) | sup_{t ∈ Σ₁∪Σ₂, \|t−L\| ≤ r₁}\|Λ_n(t)/Λ_n(L)\| ≤ C n^{κ_A} (ebenso bei L̄) | 1.0082 bis n = 2000 (V §3.22), κ_A = 0 | num; lokal, die Profile F_A sind explizit (K10) |

(P) (|P_n| ≤ Cn^{β₁}) wird nicht mehr gebraucht (V Anmerkung C); das Gesetz P_n = aCL/(nd_n)(1 + O(1/n)) bleibt als Befund (k_kern §25, num 1e-5).

### 2.2 (FA) adjungierte Fernrate mit Reserve (n-te-Wurzel-Ebene)
  **lim sup_n max_{t ∈ 𝒞 ∖ U_{r₁}(L, L̄)} (1/n)·ln|Λ_n(t)/Λ_n(L)| ≤ −η < 0** (gleichmässige Form, V Anmerkung D).
- Gemessen: Linse s = 0.03 −0.011/−0.013, Γout-Bogen −0.053/−0.037 (k52, n = 300/600); R_A fern 0.178 (n = 200) → 4.2e-8 (n = 1800) (V §3.22).
- Status **offen**. (FA) folgt aus W12′ zusammen mit Ws Modellrechnung. Der Vergleichsweg mit Superlösung scheitert (k_kern §33, O1/O2).

### 2.3 (S_E) Approximation von f_K auf E — Satz K22
  **dist_E(f_K, P_N) ≤ C_E·N^{−s_E}, s_E = 5.417·α_E, α_E = (u₁ − u₂ − 38.3°)/180°**; bei (80°, −90°) s_E = 3.96.
- Beweis: Andrievskii–Pritsker–Varga, Satz 2 (gelesen; E ∈ H*, Quasischeiben eingeschlossen), Glattheit von f_K in L aus (K3) (k_kern §35–§36).
- Eingaben: E ist eine Quasischeibe mit Spiralecke Θ ≈ 132° (num: k54 kubisch, v55 f_K; Raster und Stichproben, nicht streng).
- Gemessen s_E = 3.9–4.1 (80°/−90°, r₁ = 0.2) und 3.3–3.4 (60°/−90°) (V §3.19). Status **B + num** (Geometrie), Logik V.

### 2.4 Familie (offen, ausdrücklich)
Für alle N ≥ 400 und alle Fast-Lösungen P ∈ 𝒮_N ∩ B_N^{r₀}, für die Jacobi-Operatoren J(P)δ = δ∘P + (P′∘P)δ und die Sekanten 𝒥_Pδ = δ∘P + M_Pδ:
- **(F1)** Alle führenden Abschnitte sind invertierbar (die Köpfe existieren).
- **(F2)** (EP) und (FA) gelten gleichmässig, mit denselben Konstanten.
- **(F3′) = (Fam-koh)^{r₀}:** Es gibt einen Keim-Repräsentanten P̃ ≡ P mod x^{N+2} mit sup_{t ∈ 𝒞∖K_{r₀}}|P̃(t) − f_K(t)| ≤ δ₀ und M_{P̃} beschränkt auf 𝒞
  (Linse UND Σ₁). Das genügt, weil T_N𝒥_P nur vom Keim abhängt (k_kern §34.1).
- **(F4) Kopf-Stabilität (H-iii)^{r₀}:** sup_{K_{r₀}}|J_N(P)⁻¹xᴺ| ≤ 2·sup_{K_{r₀}}|(T_N𝒥_{P′}|_{P_N})⁻¹xᴺ|.

Stand:
- F3 alt (mit P selbst) ist falsch: sup_Linse|p_N − f_K| = 1.7e-9 (N = 200) → 1.75 (N = 2000); Σ₁ wächst langsam ab N ≈ 1000; Γout-Bogen und Σ₂ konvergieren (V §3.20).
- (F3′): Der beste Repräsentant nach Carathéodory–Fejér trägt nur bis N ≈ 6300–7800 (k_kern §34.2(b)). Für alle N braucht es Kohärenz, also (KK).
- num: (S_γ) am echten p_N bis 1000 mit derselben Konstante wie bei f_K (V §3.4). F1, F4 sind Matrix-Aussagen, numerisch wie bei f_K.

### 2.5 Endliche numerische Eingaben (noch nicht streng)
- Univalenz von f_K auf int 𝒞 (min|f_K′| = 0.336, v55), Geometrie von E (Quasischeibe, Öffnung 132°, Dreipunkt 1.07; k54, v55).
- b35b′ (min|f_K(eᶻ)| = 0.7345) und Radius r > 0.02 der lokalen Formel, aus (K) ⇒ (F).

### 2.6 Der gemeinsame Kern (KK′) (KERNPROBLEM.md §4, Nachtrag 2; korrigiert 13:07)
Die Fassung (KK) von 12:44 („am Steg wächst das Produkt primal × adjungiert mit +0.003 je Grad“) ist **falsch und zurückgezogen**: Gleichartig gemessen ist die Summe
dort ≈ 0, und das Produkt ist explizit, ΔK_n(t,t)·m(t) → κ(t) (k_kern §39, numerisch). Mit Vs Carathéodory-Fejér-Test hat die Datenpaarung der Familie am Steg
sogar eine strikte Reserve von ≈ 0.0024 je Grad (V §3.24). Es bleiben als exponentielle Hindernisse:
1. die **Komposition** der Daten: (F3′) braucht P̃ punktweise nahe f_K auf der Kontur. Je Kopf genügt der Keim mod x^{n+1}; laut V gibt es ein Linsenfenster S ≈ 0.013–0.015;
2. die **Eck-Wachstumszone** ausserhalb D nahe L (dort wächst auch das Produkt; Σ₁ bis r = 0.1 polynomial, bei r = 0.2 nicht) — also r₁ ≤ 0.1;
3. nur für den Parametrix-Weg zu (EP): die Primale an f(Γout) (k_kern §38).
Status: **numerisch**. (EP) und (F3′) bleiben benannte Hypothesen.

## 3. Beweis von Satz A und Kette des Hauptsatzes

**Beweis von Satz A** (B; Logik von V geprüft, V §1.11; k_kern §28.1, §30.2):
1. Für jedes q_N ∈ P_N gilt exakt Π_Nf − q_N = J_N⁻¹T_NJ(f − q_N).
2. Lemma K13 mit der Kontur 𝒞 und die Rang-1-Darstellung geben J_N⁻¹T_NW(z) = Σ_{n≤N}(1/2πi)∮_𝒞 ΔK_n(z,t)W(t)dt/t.
3. Mit W = J(f − q_N) gilt |W| ≤ (1 + M)·sup_E|f − q_N| auf 𝒞, denn 𝒞 ∪ f(𝒞) ⊂ E.
4. Also sup_{K_{r₀}}|Π_Nf − f| ≤ (1 + c_𝒞(1+M)Σ_{n≤N}Q_n)·dist_E(f, P_N), Q_n := max_{z∈K_{r₀}, t∈𝒞}|ΔK_n(z,t)| ≤ |L|·G_n·R_A(n), R_A(n) := max_{t∈𝒞}|Λ_n(t)/Λ_n(L)|.
5. R_A(n) ≤ Cn^{κ_A} + C_η e^{−ηn/2} nach (A-lok) und (FA); G_n ≤ Cn^{κ_G} nach (G-P).
6. Daraus Σ_{n≤N}Q_n ≤ C N^γ, und mit (S_E) folgt die Behauptung. ∎

**Kette des Hauptsatzes:**

| # | Schritt | Fundstelle | Status |
|---|---|---|---|
| 1 | Kneser zertifiziert: sexp_K ± 1.3e-168, a₀ … a₆₄₀, obere Struktur (K3), Schwanz (K2) \|b_k\| ≤ 461.36 k⁻⁶ | ARCH §9.8–9.9; V §2.2 | **V** |
| 2 | f_K holomorph auf E bis auf die Schnitte (aus (K), Univalenz) | ARCH §9.1; v55 | **V\|(K)** + num (§2.5) |
| 3 | Satz A für jedes P der Familie (Sekanten 𝒥_P statt J; Daten über P̃) | §3 oben, F2, F3′ | **B**, bedingt |
| 4 | H1 in der K_{r₀}-Norm: a-priori-Schranke über P − f = λ(Π^Pf − f) und Schritt 3; implizite Funktionen (F1); Kompaktheit; genau ein Blatt | ARCH §4.2; k_kern §23.4 | Logik **V\|(H)**; Kompaktheit **V** (V §1.11) |
| 5 | H2 (σ-Homotopie N → N+1) mit F4, K₀ = 6 | ARCH §4.3 | Logik **V\|(H),(A)**; Fassung ohne sup_{D̄} **B** (V: noch nicht im Wortlaut geprüft) |
| 6a | Zweig p₁ → p₁₉₃ | ARCH §9.7 (p18, p24, p30) | **S** (P) |
| 6b | Zweig p₁₉₃ → p₄₀₀ | V §2.1 | **S** (V, mit eigenem Code nachgeprüft) |
| 6c | Kopplung (A): sup_{D̄}\|p₄₀₀ − f_K\| ≤ 1.75e-12, also auch auf K_{r₀} | V §2.1 (v15) | **S** |
| 7 | Folgerungen (a)–(c) | §1 | **B** |

**Startbedingung bei N₀ = 400:** R₄₀₀ = 6·C′·400^{γ−s′} ≥ 1.75e-12 mit C′ = C_E·(1 + c_𝒞(1+M)C_Q), Σ_{n≤N}Q_n ≤ C_Q N^γ.
- Den Exponenten s′ ∈ (γ, s_E] im Kugelradius darf man kleiner wählen als s_E (die Schranke von Satz A gilt dann erst recht). Mit γ = 0, s′ = 2 genügt C′ ≥ 4.7e-8.
- Gemessen ist Q_n(r₀ = 0.05) = 0.236 / 0.076 / 0.0042 bei n = 200 / 600 / 1800 (V §3.22); C_E ist noch nicht gemessen.

## 4. Bausteine und ihr Status

| Baustein | Inhalt | Fundstelle | Status |
|---|---|---|---|
| P-V, p32a/(K3), (K2) | Kneser streng; f_K(L) = L, f_K′(L) = √L; Schwanz C₆ = 461.36 | ARCH §9.8–9.9; V §2.2 | **V** |
| Zertifikat 1 → 400, Kopplung (A) | Zweig der σ-Homotopie streng | ARCH §9.7; V §2.1 | **S** |
| H1, H2 | nichtlineare Kette (Existenz, Eindeutigkeit, Zweig) | ARCH §4.2–4.3 | **V\|(H)**; Umstellung auf K_{r₀}: Kompaktheit **V**, H2 **B** |
| Satz K5 / Rang-1-Ränderung | J̃_n⁻¹ = Σ c_mℓ_mᵀ/ℓ_mm; ΔK_n = δ_nΛ_n/ℓ_nn | k_kern §11.3, §24.1 | **V** (V §1.12) |
| Lemma K13 | Kopfentwicklung mit beliebigen (kopfweisen oder festen) Konturen | §22.1 | **V** (V §1.11, §1.12) |
| Satz K19 | Datenseite ohne abgeschnittene Pole: (1 + cΣQ_n)·dist_E | §30.2 | **V** (V §1.11) |
| Satz K22 | (S_E) mit s_E = 5.417·α_E | §35–§36 | **B** auf gelesenem Satz (APV 2), Logik **V**; Geometrie num (k54, v55) |
| Lemma K10.1, Satz K10, Satz K11/K11′ | Eckprofile F_β explizit; Eindeutigkeit (Mellin-Liouville) | §19 | **V** (mit Ks Vorbehalten) |
| Integraldarstellung | F_β(u) = Γ(β)⁻¹∫₀^∞t^{β−1}e^{−t−ut^{−iα}}dt; Potenzsektor per Watson-Lemma | §27.3 | **V** (V §1.12) |
| Spiralfenster | Ränder der Nicht-Wachstumsfenster sind n-unabhängige Log-Spiralen | §28.2 | **V** (V §1.12) |
| Lemma K20 | n-te-Wurzel-Rate der Primalen ≡ 0 auf D̄ | §31.1 | **V** (bis auf polare Mengen) |
| Satz K16 | Bootstrap für das Eckprodukt (bedingt) | §27.1 | **V** (wird für Satz A nicht mehr gebraucht) |
| O1 | keine superharmonische Vergleichsfunktion | §33.2 | **V** |
| Identität R_n | R_n = (I − T_n)Ju_n = u_n∘f + m·u_n − T_nk_L | §38.1 | **B** (eine Zeile) |
| T1 | primale Γout-Zone im Modell real ⇒ r₀ > r_* nötig | V §3.14 | num (V) |

## 5. Prüfstellen für V (Fassung 3)
1. **Absolute (G-P) und Budget:** Q_n ≤ |L|·G_n·R_A(n) aus Rang 1; Σ Q_n ≤ CN^γ mit γ = max(0, 1 + κ_G + κ_A). Stimmt die Buchhaltung, und reicht V §3.22 für κ_G ≈ −3.2 bei r₀ = 0.05?
2. **(F3′):** Darf man in K19 für 𝒥_P die Werte von P auf 𝒞 durch die eines Keim-Repräsentanten P̃ ersetzen (nur T_N𝒥_P geht ein)? Gilt das auch für M_P?
3. **H2 ohne sup_{D̄}** im Wortlaut von ARCH §4.3 (aus Fassung 2 offen).
4. **(KK′)** (KERNPROBLEM §4, Nachtrag 2): Die Steg-Zahl +0.003 ist zurückgezogen (13:07). Zu prüfen: h_P + h_A ≈ 0 am Steg und ΔK_n(t,t)·m(t) → κ(t) (k57, k57b) mit eigenem Code.
5. Startbedingung: s′ < s_E im Kugelradius — verträglich mit H1/H2 (R_N ≤ δ₀ der Familienannahmen)?

## 6. Für den Autor
Wenn drei Aussagen über das lineare Problem stimmen, nähern sich Deine Polynome p_N Knesers Funktion überall auf der Scheibe bis |L|, ausser in zwei winzigen
festen Kreisen um L und L̄. Daraus folgt, dass alle Koeffizienten und p_N(1) gegen Knesers Werte gehen.
- Die dritte Aussage (Annäherung durch Polynome) ist jetzt erledigt: ein bekannter Satz aus der Literatur, dazu eine am Computer geprüfte Form des Gebiets.
- Die erste (eine Antwortfolge an der Ecke bleibt klein) und die zweite (fern der Ecke wird eine andere Folge exponentiell klein) gelten in allen Rechnungen bis Grad 2000.
- Dieselben Aussagen müssen auch für die Nachbarpolynome im Newton-Schritt gelten (die „Familie“). Das ist offen.
- Neu verstanden (13:07): Knapp ausserhalb der Scheibe wird eine Hälfte der Rechnung gross und die andere genau gleich schnell klein. Das Produkt bleibt also harmlos,
  und es lässt sich sogar mit einer einfachen Formel beschreiben. Offen bleibt vor allem, ob die Nachbarpolynome dort nahe genug an Knesers Funktion bleiben.

## 7. Neu bewiesen in Phase 3 (für den Wissensgraphen)

| Satz / Lemma | Inhalt (kurz) | Fundstelle | Strang | von V gegengelesen? |
|---|---|---|---|---|
| Satz K1, Lemma K1.0 | Kurvenkriterium: polynomiale G_N-Schranke auf einer festen Kurve ⇒ (H₀); Cauchy-Darstellung der Abschnittslösung | k_kern §1 | K | **ja** (V §1.3) |
| Lemma K1′, K1″ | Stabilität in der konformen Gleichungskoordinate ⇒ (H₀^{K_γ}) | k_kern §3, §5.1 | K | **ja** (V §1.6) |
| Lemma K2 | Kompositions-Vorkonditionierung vertauscht mit dem Abschneiden | k_kern §0 Nr. 3 | K | **ja** (V §1.12; „Phasenraumgewicht“ nur Lesart) |
| Satz K3 | Homogenität exakter linearer Identitäten (verallgemeinert W13) | k_kern §9 | K | **ja** (V §1.12) |
| Satz K4 | Strukturidentität P′∘P = (eˣ + R_P′)/P′ | k_kern §9 | K | **ja** (V §1.12) |
| Satz K5 | exakte Kopf-Zerlegung (LU) der Abschnittsinversen | k_kern §11.3 | K | **ja** (V §1.12; inhaltsgleich G9/G22) |
| Lemma K6 | Linsenpunkte über das Grad-Antwortprofil | k_kern §11.2 | K | **ja** (V §1.12) |
| Lemma K10.1, Satz K10, K11/K11′ | Gamma-Zerlegung; explizite Eckprofile; Eindeutigkeit | k_kern §19 | K | **ja** (V §1.10) |
| Grenzprofil des adjungierten Kopfes | Mellin–Barnes-Form, Test gegen G | k_kern §18 | K | **ja** (V §1.9) |
| Lemma K13 | Kopfentwicklung mit beliebigen Konturen | k_kern §22.1 | K | **ja** (V §1.11, §1.12) |
| Rang-1-Zuwachs / P_n = K_n − K_{n−1} | Eckprodukt = Zuwachs des Rayleigh-Quotienten ẽ_L*J̃_n⁻¹ẽ_L | k_kern §24.1, §28.1 | K | **ja** (V §1.12, numerisch v48) |
| sup_{D̄}\|δ_n\| ≤ √(n+1)‖c_n‖ | primale Randschranke | k_kern §23.1 | K | **ja** (V §1.12) |
| Reduktionssatz, lokale Fassung | lokale Beschränktheit genügt für die Konvergenz in D | k_kern §23.2 | K | **ja** (V §1.12) |
| Satz K16 | Bootstrap für das Eckprodukt (bedingt) | k_kern §27.1 | K | **ja** (V §1.12) |
| Integraldarstellung F_β, Watson-Asymptotik | Eckprofile als Gamma-Integral; Potenzsektor elementar | k_kern §27.3 | K | **ja** (V §1.12; unabhängig numerisch über V(ζ), 1e-12) |
| Spiralfenster | n-unabhängige Log-Spiralen | k_kern §28.2 | K | **ja** (V §1.12) |
| Satz K19 | Datenseite über Approximation auf E | k_kern §30.2 | K | **ja** (V §1.11) |
| Satz K22 | (S_E): dist_E(f_K, P_N) ≤ C_E N^{−5.417α_E} (APV Satz 2, Glattheit aus (K3), δ_T = 0.030) | k_kern §35–§36 | K | **ja** (V §1.12; APV von K gelesen; Geometrie und Rate numerisch V §3.19, §3.23) |
| Identität des Abschneide-Rests | R_n = (I − T_n)Ju_n = u_n∘f + m·u_n − T_nk_L (Galerkin, eine Zeile) | k_kern §38.1 | K | nein |
| Lemma K20 | primale n-te-Wurzel-Rate ≡ 0 auf D̄ | k_kern §31.1 | K | **ja** (V §1.12; für den punktweisen limsup nur bis auf polare Mengen) |
| O1 | keine superharmonische Vergleichsfunktion (Laplace-Masse) | k_kern §33.2 | K | **ja** (V §1.12) |
| (K2) | \|b_k\| ≤ 461.36·k⁻⁶, Σ_{k>640}\|b_k\| ≤ 8.6e-13 | V §2.2 | V (Rechnung) | zweiter Leser K (§7) |
| Zweig 193 → 400 | σ-Homotopie streng | V §2.1 | V | unabhängig nachgerechnet (V) |
| Sprungschranke | \|Δf_K(L(1+u))\| ≤ 1.0378·u^{5.4469}, 0 < u ≤ 0.48 | V §3.15 | V | nein |
| Satz V2 | (S_γ) mit B ≤ 16.4 ab N₀ = 400 ⇒ Hauptsatz (Ks Linse) | V §3.9, §3.15 | V | nein (streng bis auf Xs Lemmata und Fam) |
| (a″)_∞ Kneser | σ_min ≥ 0.641 für alle N | m_modell_eckmode §15, §17 | M | **ja** (V §1.8) |
| G8–G11, G22, G23, G39, G40, G43, G44 | Diagonalisierung, Rang-1-Zuwachs, Summenkriterium, Modelle | g_gpt_hilfssaetze | G | **ja** (V §1.2) |

Zurückgezogen in Phase 3: Satz K17 (k_kern §28.4 → §29.1). Gestrichen: die Hilfsannahme (EP⁺) (k_kern §37.1). Widerlegt: F3 in der Form von Fassung 2 (V §3.20). Skizzen ohne vollständigen Beweis: K14, K15 (Datenseite), Lemma K12, die Herleitung von C (§25), das Vergleichsgerüst (§33.1).

---

## Anhang A: Fassung 1 (12.10. 10:23), unverändert

Stand 12.10.2026 10:23 (Strang K). Zum Gegenlesen als Ganzes durch V. Massgeblich für die Kette sind ARCHITEKTUR.md §4 und §9.

Statuswörter (wie in v_pruefstand/ABHAENGIGKEIT.md):
- **V**: von einem zweiten Strang Schritt für Schritt geprüft.
- **S**: rigoros gerechnet (arb), vom Autor selbst geprüft.
- **B**: vom Autor bewiesen, noch nicht gegengelesen.
- **num**: numerisch, für endlich viele N.
- **offen**: weder bewiesen noch gerechnet.

### A.1. Aussage

**Bezeichnungen.**
- L = 0.31813 + 1.33724i ist der Fixpunkt e^L = L, ρ = |L|, D = {|x| < ρ}.
- f_K ist Knesers halbe e-Funktion.
- p_N sind die Polynome des Autors vom Grad N mit T_N[p_N∘p_N − eˣ] = 0, als Zweig der σ-Homotopie ab p₁.
- Für ein festes r₀ ∈ (0, 0.1] ist **K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄))**.
- B_N^{r₀} := {P : sup_{K_{r₀}}|P − f_K| ≤ R_N}, R_N = K₀C′N^{κ−4.29}, K₀ = 6.

**Hauptsatz.** Es gelte die Annahme (H)^{r₀} aus §2. Dann gilt für alle N ≥ N₀ = 400:
- p_N ist die einzige Lösung von T_N[P∘P − eˣ] = 0 in B_N^{r₀}.
- **sup_{K_{r₀}} |p_N − f_K| ≤ C′·N^{κ−4.29}.** Wegen κ < 4.29 geht das gegen 0.
- Für N < 400 sind die p_N streng eingeschlossen (Zertifikat, §3 Schritt 7).

**Folgerungen** (B, je eine Zeile):
- (a) **Koeffizienten.** Der Kreis |x| = ρ − r₀ liegt in K_{r₀}. Mit Cauchy folgt |a_k^{(N)} − a_k| ≤ (ρ − r₀)^{−k}·C′N^{κ−4.29}, also a_k^{(N)} → a_k für jedes k.
- (b) **p_N(1) → f_K(1) = sexp_K(½) = T_e(½).** Es ist 1 ∈ K_{r₀}, denn |1| < ρ und |1 − L| = 1.50 > r₀.
- (c) **Identifikation mit Kneser.** Sie ist eingebaut, weil die Kugel um f_K liegt. Ein eigener Eindeutigkeitssatz ist nicht nötig.

**Nicht behauptet** ist die gleichmässige Konvergenz bis an L heran. Dafür gibt es die Warnkette in k_kern/BERICHT §22.6, entschieden wird sie durch Test T1 (V). Diese Kette betrifft nur Punkte z nahe L, also nicht K_{r₀}.

### A.2. Alle Annahmen

**Annahme (H)^{r₀} — die eine offene Aussage für alle N.** Es gibt C, κ < 4.29, ε ∈ (0, 1] und d = 0.3, so dass für alle N ≥ 400 und alle Fast-Lösungen P ∈ 𝒮_N ∩ B_N^{r₀} gilt:
- **(H-i)** J_N(P) und T_N𝒥_P|_{P_N} sind invertierbar. Dabei ist J(P)δ = δ∘P + (P′∘P)δ und 𝒥_Pδ = δ∘P + M_Pδ.
- **(H-ii)^{r₀}** Für z ∈ K_{r₀} und ζ ∈ S_L(d) ∪ S_L̄(d) ∪ Γ_d (Rahmen K′ = D̄ ∪ f_K(D̄)) gilt |E^P_N(ζ, z)| ≤ C·N^κ·e^{−εN g_{K′}(ζ)}/dist(ζ, K′). Dabei ist E^P_N(ζ, z) := Π^P_N k_ζ(z) − k_ζ(z) die lokale Pol-Antwort.
- **(H-iii)^{r₀}** sup_{K_{r₀}}|J_N(P)⁻¹xᴺ| ≤ 2·sup_{K_{r₀}}|(T_N𝒥_{P′}|_{P_N})⁻¹xᴺ| für alle P, P′ in derselben Menge (Kopf-Stabilität).

Das ist wörtlich die Annahme (H) aus ARCHITEKTUR §5.1, mit zwei Änderungen: D̄ ist überall durch K_{r₀} ersetzt, und die Kugel ist in sup_{K_{r₀}} gemessen.

Status:
- Für alle N: **offen**.
- num: (H-ii) in der stärkeren D̄-Form bis N = 1000 (Θ ≤ 0.11, ARCH §9.3).
- num: die Linsenform (S_γ) am echten p_N bis N = 1000 mit derselben Konstante wie bei f_K (V §3.4).
- num: (H-iii) mit Quotient 1.0000 bei N = 60 (ARCH §5.3).

**Endliche numerische Eingaben**, die noch nicht streng sind (aus (K) ⇒ (F), ARCH §9.1 und v_pruefstand/ABHAENGIGKEIT.md Nr. 4):
- b35b′: min|f_K(eᶻ)| = 0.7345.
- Radius r > 0.02 der lokalen Formel.
- Univalenz: Gitterminimum c_out = 0.3456.
- Untere Schranke der Green-Funktion von K′ auf dem Schlitz: g′(L(1+u)) ≥ 0.754·u^{1/α′} für u ≤ 0.204 (p13b).
- Abstand dist(L(1+u), K′) ≥ 0.933ρu.

Weitere Annahmen gibt es nicht. Insbesondere braucht diese Fassung **keine Eckbedingung (EB_N)**. Begründung in §5.

### A.3. Beweiskette und Bausteine

| # | Schritt | Fundstelle | Status |
|---|---|---|---|
| 1 | Kneser zertifiziert: sexp_K ± 1.3e-168, a_0 … a_640 (Fehler·ρⁿ ≤ 2.9e-127), obere Struktur (K3), Schwanz (K2) mit \|b_k\| ≤ 461.36·k⁻⁶ | ARCH §9.8, §9.9; V §2.2 | **V** |
| 2 | (K) ⇒ (F): f_K holomorph um K′ bis auf die Schlitze, Sprung ≤ C_J\|ζ − L\|^{5.447} | B L8.2–8.4, ARCH §9.1 | **V\|(K)** + num-Eingaben aus §2 |
| 3 | Sprungschranke C_J = 1.0378 für 0 < u ≤ 0.48 | V §3.15 | **S** (V), kein zweiter Leser |
| 4 | P-I(i), punktweise in z: (H-ii)^{r₀} ⇒ sup_{K_{r₀}}\|Π^P_N f_K − f_K\| ≤ C′N^{κ−4.29}. Beweis: Cauchy-Darstellung über die Schlitzufer und Γ_d. Er ist punktweise in z, die Einschränkung auf K_{r₀} gilt also wörtlich. | ARCH §4.1, §9.4 | **V** (P-I); Einschränkung auf K_{r₀}: **B** (eine Zeile) |
| 5 | H1 in der K_{r₀}-Norm. (1) a-priori-Schranke: Sekanten-Identität und Schritt 4 für Π^P. (2) implizite Funktionen mit (H-i). (3) Kompaktheit: P_N ist endlichdimensional und K_{r₀} hat innere Punkte, also ist sup_{K_{r₀}} eine Norm auf P_N. (4) Es gibt genau ein Blatt. | ARCH §4.2; k_kern §23.4 | Logik **V\|(H)** (P); Umstellung auf K_{r₀}: **B** |
| 6 | H2 (σ-Homotopie N → N+1) mit (H-iii)^{r₀} und K₀ = 6 | ARCH §4.3 | Logik **V\|(H),(A)**; Umstellung: **B** |
| 7a | Zweig p₁ → p₁₉₃ als σ-Homotopie | ARCH §9.7 (p18, p24, reproduziert in p30) | **S** (P) |
| 7b | Zweig p₁₉₃ → p₄₀₀ | V §2.1 | **S** (V, mit eigenem Code unabhängig nachgeprüft) |
| 7c | Kopplung (A): sup_{D̄}\|p₄₀₀ − f_K\| ≤ 1.75e-12. Das gilt erst recht auf K_{r₀}. | V §2.1 (v15) | **S** |
| 8 | Folgerungen (a)–(c) | §1 | **B** |

Ablauf: Aus 1–4 folgt die lineare Konvergenz auf K_{r₀} für jedes P der Familie. Mit 5 und 6 bleibt der Zweig ab N₀ in der Kugel. 7 setzt bei N₀ = 400 an.

### A.4. Explizite Konstanten für N₀ = 400

| Grösse | Wert | Status |
|---|---|---|
| Sprungkonstante C_J (u ≤ 0.48) | 1.0378 | S |
| Green-Faktor von K′ auf dem Schlitz, c_g′ (u ≤ 0.204) | ≥ 0.754 | num (p13b) |
| Abstandsfaktor dist(L(1+u), K)/(ρu) | ≥ 0.933 | num (W̄-Wert; für K′ nachprüfen, Prüfstelle 1) |
| Exponent der Schlitzrate 5.447·α′ (α′ = 0.787172, Aussenwinkel 141.69°) | 4.2877 | B (P-I) |
| Schlitzkonstante 2·C_J/(0.933π)·α′Γ(4.2877)·c_g′^{−4.2877} | 16.3 | K; siehe Prüfstelle 1 |
| Konturkonstante bei d = 0.3: M_d\|Γ_d\|/(2πδ_d), mit \|Γ\| = 12.871, δ = 0.5095, M_d = 4.251 | 17.1 | num (X) |
| Kugelfaktor K₀ | 6 | V (H2) |
| Kopplung bei N = 400 | 1.75e-12 | S |
| Kneser-Schwanz: C₆ = 461.36, Σ_{k>640}\|b_k\| ≤ 8.6e-13 | — | V |

Daraus: C′ = C·(16.3·ε^{−4.29} + 17.1·N^{4.29}e^{−0.3εN}). Der Konturteil ist bei N = 400 für ε ≥ 0.4 kleiner als 1e-8·C (gemessen ist ε ≥ 0.60, W).

Bedingung für den Start bei N₀ = 400: R₄₀₀ = 6·C′·400^{κ−4.29} ≥ 1.75e-12, also (bei ε ≤ 1) **C·400^κ ≥ 2.6e-3**. Das ist schwach; V gibt für (A) eine Reserve von etwa 1e8 gegen R_{N₀} an (ABHAENGIGKEIT Nr. 20).

Die eigentliche Last liegt darin, dass (H)^{r₀} ab N = 400 gilt.

Leitet man den Familienteil statt dessen aus P = f_K her (Weg X4a/X8), so kommt eine Bedingung hinzu: A_θ·N·R_N ≤ ½ mit A_θ = 13.4 (X8′, B). Bei C_E ≤ 1 und κ_E ≤ 2 gibt das N₀ ≈ 318 (x_alternativen/KONSTANTEN.md §000).

### A.5. Verhältnis zu (EB_N) und zu (AH)

**(EB_N) — Antwort auf die Frage der Leitung.** Für den Hauptsatz in der Form von §2 ist (EB_N) **keine Annahme**: P wird nur auf K_{r₀} gebraucht.

(EB_N) wird erst nötig, wenn man (H)^{r₀} für die Familie aus der einen Aussage bei P = f_K herleiten will. Der Störungsschritt (X8/X4a) schätzt 𝒥_P − J_∞ über sup|P − f| auf ganz ∂D ab, also auch nahe L.
- Auf der Kegelscheibe folgt (EB_N) aus der a-priori-Schranke, wenn (H-ii) auch für z in der Kegelscheibe und den Ruhefenstern gilt. Dort greift die Warnkette nicht (A ≤ 16).
- Das Band zwischen Kegelscheibe und U_{r₀}(L) ist damit **nicht** gedeckt. Für diesen Weg wäre das eine zweite Aussage („Bandkontrolle“).
- Gemessen gilt (EB_N) am echten p_N: η_N/r_N ≤ 1.5e-3 ab N = 800 (V §3.18).

**(AH) „Eckdominanz“ (k_kern/BERICHT §28).** (AH) ist die gemessene Strukturaussage über die Ecke:
- Die Köpfe sind auf einer festen Spiralkontur durch L in der Ecke am grössten.
- Ab n ≈ 600 gilt Q_n = 1.008·|P_n| (num, n ≤ 600).

**Korrektur zu §28.4: Satz K17 ist zurückgezogen.** Behauptet war: (AH) ⇒ (H-ii)^{r₀} für alle N. Das ist nicht gezeigt.
- Für Quellen ζ auf dem Schlitz nahe L muss f(Kontur) nahe L in D̄ liegen.
- Das feste Fenster der Spirale verlässt diesen Sektor aber unterhalb von r ≈ 1e-6 bis 1e-11 (je nach Fensterwinkel).
- Das Argument trägt deshalb nur bis N ≈ 1e8 bis 1e11.

(AH) bleibt eine gemessene Eckaussage. Ob sie für alle N hinreicht, ist **offen**.

**Exponentielle Ebene (k50, n = 300, num)** — Ausgangspunkt für den nächsten Block:
- Auf exponentieller Ebene lässt sich die Ecke im festen Abstand δ umgehen. Mit dem Datengewicht (|f(t)|/ρ)ⁿ ist G_n < 0 auf den Bögen |t − L| = δ ∈ [0.003, 0.1] in den Richtungen 130° bis 265°.
- Dort dreht sich nichts: Die Daten fallen linear in δ, der Kopf wächst nur wie δ^{1.118}.
- Auf der Γout-Seite in D ist dieses Gewicht dagegen zu schwach (G ≥ +0.03). Dort braucht es den Rahmen K′ wie in (H-ii).

### A.6. Prüfstellen für V

1. **Schlitzkonstante.** Nach der Formel in ARCH §9.7 (1) ergibt sich mit den W̄-Werten (C_J = 0.95, c_g = 0.829, α = 0.883) 25.5, angegeben ist dort 12.7. Ist der Faktor 2 (beide Ufer) doppelt gezählt? Davon hängt ab, ob es für K′ 16.3 oder 8.2 ist.
   Zudem passen Rahmen und Schlitzlänge noch nicht zusammen: Bei d = 0.3 reicht der Schlitz in K′ über u = 0.204 hinaus, und bis dorthin ist c_g′ nur numerisch belegt. Zwei Auswege:
   - c_g′ bis u₀′(0.3) nachrechnen;
   - oder Rahmen W̄ nehmen (Rate 4.81, Konstanten aus ARCH §9.7 (1), Schlitz bis u₀ = 0.451, durch C_J bis 0.48 gedeckt).
2. **Schritt 5 (3).** Trägt die Kompaktheit von 𝒵 wörtlich, wenn die Kugel nur in sup_{K_{r₀}} gemessen wird? Ist 𝒮_N mit dieser Norm abgeschlossen?
3. **Schritt 4 für die Sekanten.** Die Cauchy-Darstellung von Π^P f_K über Kerne k_ζ braucht nur den Keim von 𝒥_P k_ζ in 0, also P(0) ≠ ζ. Stimmt es, dass P nahe L dabei nirgends eingeht?
4. **Schritt 6.** (H-iii) ist hier in sup_{K_{r₀}} formuliert. Prüfen, ob H2 an keiner Stelle sup_{D̄} benutzt (Homotopie-Integral, Ball-Erhalt).
5. **(H-ii)^{r₀}** ist schwächer als die alte (H-ii) und von der Warnkette nicht betroffen. Ob sie für alle N gilt, ist die offene Frage. Weil Quellen auf dem Schlitz nahe L bleiben, ist sie aber nicht frei von der Ecke.

### A.7. Für den Autor

Wir zeigen: Wenn eine einzige Aussage über das lineare Problem stimmt — die Annahme (H)^{r₀} —, dann nähern sich Deine Polynome p_N Knesers Funktion überall auf der Kreisscheibe bis |x| = |L|, nur nicht in zwei winzigen festen Kreisen um die Fixpunkte L und L̄. Daraus folgt, dass jeder Koeffizient und der Wert p_N(1) gegen Knesers Werte gehen.

Alles andere ist bewiesen oder vom Computer streng nachgerechnet, bis Grad 400. Einige Randkonstanten sind nur numerisch geprüft (§2, §4).

Die offene Aussage gilt in allen Rechnungen bis Grad 1000. Bewiesen ist sie für alle Grade noch nicht.

---

## Anhang B: Fassung 2 (12.10. 11:05), unverändert (ohne §7, die Liste lebt oben weiter)

**Fassung 2**, Stand 12.10.2026 11:05 (Strang K). Sie ersetzt Fassung 1 (10:23, unverändert im Anhang A). Zum Gegenlesen als Ganzes durch V.

Grundlagen: ARCHITEKTUR.md §4 und §9 (nichtlineare Kette), k_kern/BERICHT.md §24–§33 (lineare Aussage).

Statuswörter:
- **V**: von einem zweiten Strang Schritt für Schritt geprüft.
- **S**: rigoros gerechnet (arb), vom Autor selbst geprüft.
- **B**: vom Autor bewiesen, noch nicht gegengelesen.
- **num**: numerisch, für endlich viele N.
- **offen**: weder bewiesen noch gerechnet.

Was sich gegenüber Fassung 1 ändert:
- Die eine Annahme (H)^{r₀} (Pol-Antwort für ALLE Quellen ζ und die ganze Familie) wird zerlegt. Bei P = f_K bleiben (EP) (zwei Funktionenfolgen, polynomial), (FA) (n-te-Wurzel-Ebene) und (S_E) (klassische Approximation).
- Die Datenseite läuft über K19 statt über P-I. Damit gibt es keine abgeschnittenen Pole mehr und keine Eckbedingung (EB_N).
- Die Familie ist als eigener offener Punkt ausgewiesen (F1–F4).
- r₀ muss grösser sein als der Abstand r_* der primalen Γout-Zone (Vs T1).

### B.1. Aussage

**Bezeichnungen.**
- L = 0.31813 + 1.33724i (e^L = L), ρ = |L|, D = {|x| < ρ}, μ = √L, a = 1 + iL/(4π), r_n = ρn^{−Re a} (Kegelskala).
- f_K ist Knesers halbe e-Funktion. p_N sind die Polynome des Autors (T_N[p_N∘p_N − eˣ] = 0, Zweig der σ-Homotopie ab p₁).
- **K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄))** mit festem r₀ ∈ (r_*, 0.1]. Dabei ist r_* der Abstand der Zone, in der der primale Kopf auf der Γout-Seite zu wachsen beginnt:
  - für Kneser nach der Profilrechnung r_* ~ 1e-6 (k_kern/BERICHT §23.3; Schätzung);
  - im Modell M4 von V etwa 0.13 (Vs T1).
- Köpfe: c_n = J̃_n⁻¹e_n, ℓ_n = J̃_n^{−T}e_n, δ_n(z) = Σ c_{n,k}(z/ρ)^k, Λ_n(t) = Σ ℓ_{n,k}(ρ/t)^k.
  Kern K_n(z, t) = ε_zᵀJ̃_n⁻¹ẽ_t = t·(J_n⁻¹T_nk_t)(z); Zuwachs ΔK_n = K_n − K_{n−1} = δ_n(z)Λ_n(t)/ℓ_nn (Rang 1).
- Die zwei **Eck-Pol-Antworten**:
  - u_n := J_n⁻¹T_nk_L = Π_n s (primal; s = J_∞⁻¹k_L hat einen Pol in L);
  - w_n(t) := K_n(L, t) (adjungiert, ausgewertet in L).
  - Ferner gilt **P_n := ΔK_n(L, L)** (Eckprodukt).
- **Kontur 𝒞.** Sie ist geschlossen um 0 und fest, für alle n dieselbe:
  - die Spiralen Σ₁, Σ₂ in L: θ(r) = u_A + 71.46° + 1.623°·ln(r/ρ), mit u_A ≈ 80° (ausserhalb D) bzw. u_A ≈ −90° (in D), 0 < r ≤ r₁ = 0.05;
  - die Linse |t| = ρ(1 + 0.03) über Γin;
  - der Γout-Bogen |t| = ρ;
  - konjugiert bei L̄.
- **Kompaktum E** ⊃ D̄ ∪ 𝒞 ∪ f_K(𝒞) ∪ f_K(Inneres von 𝒞): f_K ist darauf stetig und im Innern holomorph. Ein Spiralschnitt von L aus läuft in der Lücke ≈ (25°, 124°), die E frei lässt; ebenso bei L̄.

**Satz A (lineare Konvergenz bei f_K; Beweis §3, Status B).** Es gelten (EP), (FA) und (S_E) aus §2. Dann gilt für alle N
  **sup_{K_{r₀}} |Π_N f_K − f_K| ≤ C·N^{γ − s_E},   γ := 1 + β₁ + κ_P + κ_A < s_E.**

**Hauptsatz (Status: bedingt).** Es gelten die Annahmen von Satz A für die ganze Familie, also F1–F4 in §2.4. Dann gilt für alle N ≥ N₀ = 400:
- p_N ist die einzige Lösung von T_N[P∘P − eˣ] = 0 in B_N^{r₀} := {P : sup_{K_{r₀}}|P − f_K| ≤ R_N}, R_N = K₀C′N^{γ−s_E}, K₀ = 6.
- **sup_{K_{r₀}} |p_N − f_K| ≤ C′·N^{γ − s_E} → 0.**
- Für N < 400 sind die p_N streng eingeschlossen (Zertifikat, §3).

**Folgerungen** (B, je eine Zeile):
- (a) |a_k^{(N)} − a_k| ≤ (ρ − r₀)^{−k}·C′N^{γ−s_E}, also a_k^{(N)} → a_k für jedes k (Cauchy auf |x| = ρ − r₀ ⊂ K_{r₀}).
- (b) p_N(1) → f_K(1) = sexp_K(½) = T_e(½), denn 1 ∈ K_{r₀} (|1 − L| = 1.50).
- (c) Die Identifikation mit Kneser ist eingebaut, weil die Kugel um f_K liegt.

**Nicht behauptet** ist die gleichmässige Konvergenz bis an L heran. Vs T1 hat das in einem Modell mit schneller Spirale (M4) gezeigt:
- In einer Γout-Zone in festem Abstand von L konvergiert Π_N U − U seit N ≈ 1000 nicht mehr (Plateau bis N = 2800).
- In L selbst, auf Γin und im Innern konvergiert es weiter.
- Deshalb ist **r₀ > r_* nötig** und nicht nur bequem.

### B.2. Alle Annahmen

#### B.2.1 (EP) Eck-Pol-Hypothese: polynomiale Zuwächse der beiden Eck-Pol-Antworten
| Teil | Aussage (für alle n) | gemessen | Status |
|---|---|---|---|
| (P) | \|P_n\| ≤ C n^{β₁} | 0.535·n^{−0.106} (k46); \|K_n(L,L)\|/n^{0.894} → 0.599 bis n = 2000, kubisch = Kneser (V §3.17) | num; Gesetz P_n = aCL/(n d_n)(1 + O(1/n)) mit lokal hergeleitetem C (§25, Skizze + num auf 1e-5); K16 reduziert es auf Lemma ρ |
| (G-P) | R_P(n) := sup_{z∈K_{r₀}}\|Δu_n(z)\|/\|Δu_n(L)\| ≤ C n^{κ_P} | Inneres ~ g(z)·n^{−3.05} = n^{−Re(a s₀^P)}; Rand fern von L vorübergehend exponentiell (k51, k51b, n ≤ 1000) | num; nach K20 hat sie keine n-te-Wurzel-Reserve, ist also echt polynomial |
| (A-lok) | sup_{t ∈ Σ₁∪Σ₂, \|t−L\| ≤ r₁}\|Δw_n(t)\|/\|Δw_n(L)\| ≤ C n^{κ_A} (ebenso bei L̄) | ≤ 1.008 (k49, n ≤ 600) | num; lokal, die Profile F_A sind explizit (K10) |

**Budget:** γ = 1 + β₁ + κ_P + κ_A < s_E. Gemessen ist γ ≈ 0.89.

Bemerkung: (G-P) ist der Eckteil der alten Pol-Antwort-Hypothese an der Quelle ζ = L (k_kern/BERICHT §32). Ein Beweis über die äussere Entwicklung bräuchte eine polynomiale Parametrix; diese prüft V numerisch (Vortest (b), BRETT 11:00).

#### B.2.2 (FA) adjungierte Fernrate mit Reserve (n-te-Wurzel-Ebene)
  **lim sup_n max_{t ∈ 𝒞 ∖ U_{r₁}(L, L̄)} (1/n)·ln|Λ_n(t)/Λ_n(L)| ≤ −η < 0.**

Gemessen (k52, n = 300/600): Linse s = 0.03: −0.011/−0.013; Γout-Bogen (≥ 10° von L, L̄): −0.053/−0.037.

Status **offen**:
- (FA) folgt aus W12′ (Gestalt des Grenzproblems bis auf Konstanten eindeutig) zusammen mit Ws Modellrechnung.
- Der Vergleichsweg mit Superlösung scheitert (§33, Hindernisse O1 und O2).
- W13 stört nicht, weil nur die relative Rate zur Ecke zählt.

#### B.2.3 (S_E) Approximation von f_K auf E
  **dist_E(f_K, P_N) := inf_{q∈P_N} sup_E|f_K − q| ≤ C_E·N^{−s_E}.**
- Erwartet ist s_E ≈ 5.447·α_E, wobei α_Eπ ≈ 99°…109° der Aussenwinkel von E in L ist; also s_E ≈ 3.0…3.6.
- Theorie: Dzjadyk-Typ, für die Spiralecke nicht nachgeprüft. Für Konvergenz ohne Rate genügt Mergelyan; die Rate braucht s_E > γ.
- Status: **offen**; V misst s_E.

#### B.2.4 Familie (offen, ausdrücklich)
Für alle N ≥ 400 und alle Fast-Lösungen P ∈ 𝒮_N ∩ B_N^{r₀}, sowohl für die Jacobi-Operatoren J(P)δ = δ∘P + (P′∘P)δ als auch für die Sekanten 𝒥_Pδ = δ∘P + M_Pδ:
- **(F1)** Alle führenden Abschnitte sind invertierbar (die Köpfe existieren).
- **(F2)** (EP) und (FA) gelten gleichmässig, mit denselben Konstanten.
- **(F3) Datenbedingung auf 𝒞:** sup_{t∈𝒞}|P(t) − f_K(t)| ≤ δ₀ und M_P beschränkt auf 𝒞. Nötig, weil K19 für 𝒥_P den Wert von P auf 𝒞 ∖ K_{r₀} braucht (Linse, Spiralen an L); die K_{r₀}-Kugel kontrolliert ihn nicht.
- **(F4) Kopf-Stabilität (H-iii)^{r₀}:** sup_{K_{r₀}}|J_N(P)⁻¹xᴺ| ≤ 2·sup_{K_{r₀}}|(T_N𝒥_{P′}|_{P_N})⁻¹xᴺ|.

Stand:
- num: (S_γ) am echten p_N bis 1000 mit derselben Konstante wie bei f_K (V §3.4).
- num: (EB_N) am echten p_N: η_N/r_N = 5e-5 bei N = 2000 (V §3.18); das stützt F3 nahe L.
- Zu prüfen (PLAN): ob Keimabhängigkeit (T_N𝒥_P hängt nur vom Keim von P ab) zusammen mit K4 (P′∘P = (eˣ + R_P′)/P′) F2/F3 auf P = f_K zurückführt, wie bei der Linse (§10).

#### B.2.5 Endliche numerische Eingaben (noch nicht streng)
- Univalenz von f_K auf D̄ ∪ Linse (Gitterminimum c_out = 0.3456). Sie wird gebraucht für: f_K holomorph auf E, und m = eˣ/f_K′ beschränkt auf 𝒞.
- b35b′ (min|f_K(eᶻ)| = 0.7345) und Radius r > 0.02 der lokalen Formel, aus (K) ⇒ (F).

Die Green- und Abstandsfaktoren von Fassung 1 (c_g′, 0.933) werden in Fassung 2 **nicht mehr** gebraucht; an ihre Stelle tritt (S_E).

### B.3. Beweis von Satz A und Kette des Hauptsatzes

**Beweis von Satz A** (B; k_kern/BERICHT §28.1, §30.2):
1. Für jedes q_N ∈ P_N gilt exakt Π_Nf − q_N = J_N⁻¹T_NJ(f − q_N).
2. Lemma K13 mit der Kontur 𝒞 für alle Köpfe und die Rang-1-Darstellung geben
   J_N⁻¹T_NW(z) = Σ_{n≤N}(1/2πi)∮_𝒞 ΔK_n(z,t)W(t)dt/t.
3. Mit W = J(f − q_N) gilt |W| ≤ (1 + M)·sup_E|f − q_N| auf 𝒞, denn 𝒞 ∪ f(𝒞) ⊂ E.
4. Also sup_{K_{r₀}}|Π_Nf − f| ≤ (1 + c_𝒞(1+M)Σ_{n≤N}Q_n)·dist_E(f, P_N) mit Q_n := max_{z∈K_{r₀}, t∈𝒞}|ΔK_n(z,t)| = |P_n|·R_P(n)·R_A(n).
5. Aufteilung von R_A:
   - fern (𝒞 ∖ U_{r₁}): ≤ C_η e^{−ηn/2} nach (FA);
   - nah: ≤ Cn^{κ_A} nach (A-lok);
   - R_P ≤ Cn^{κ_P} nach (G-P), |P_n| ≤ Cn^{β₁} nach (P).
6. Daraus Σ_{n≤N}Q_n ≤ C N^γ, und mit (S_E) folgt die Behauptung. ∎

**Kette des Hauptsatzes:**

| # | Schritt | Fundstelle | Status |
|---|---|---|---|
| 1 | Kneser zertifiziert: sexp_K ± 1.3e-168, a₀ … a₆₄₀, obere Struktur (K3), Schwanz (K2) \|b_k\| ≤ 461.36 k⁻⁶ | ARCH §9.8–9.9; V §2.2 | **V** |
| 2 | f_K holomorph auf E bis auf die Schnitte (aus (K), Univalenz) | ARCH §9.1 | **V\|(K)** + num (§2.5) |
| 3 | Satz A für jedes P der Familie (Sekanten 𝒥_P statt J) | §3 oben, F2, F3 | **B**, bedingt |
| 4 | H1 in der K_{r₀}-Norm: a-priori-Schranke über die Sekanten-Identität P − f = λ(Π^Pf − f) und Schritt 3; implizite Funktionen (F1); Kompaktheit (P_N endlichdimensional, sup_{K_{r₀}} ist eine Norm); genau ein Blatt | ARCH §4.2; k_kern §23.4 | Logik **V\|(H)**; Umstellung **B** |
| 5 | H2 (σ-Homotopie N → N+1) mit F4, K₀ = 6 | ARCH §4.3 | Logik **V\|(H),(A)**; Umstellung **B** |
| 6a | Zweig p₁ → p₁₉₃ | ARCH §9.7 (p18, p24, p30) | **S** (P) |
| 6b | Zweig p₁₉₃ → p₄₀₀ | V §2.1 | **S** (V, mit eigenem Code nachgeprüft) |
| 6c | Kopplung (A): sup_{D̄}\|p₄₀₀ − f_K\| ≤ 1.75e-12, also auch auf K_{r₀} | V §2.1 (v15) | **S** |
| 7 | Folgerungen (a)–(c) | §1 | **B** |

**Startbedingung bei N₀ = 400:** R₄₀₀ = 6·C′·400^{γ−s_E} ≥ 1.75e-12.
- C′ = C_E·(1 + c_𝒞(1+M)C_Q), wobei C_Q die Konstante in Σ_{n≤N}Q_n ≤ C_Q N^γ ist.
- Gemessen ist Q_n mit z ∈ K_{0.05}, ohne z = L und mit grober ferner Linse: 0.54 / 0.17 / 0.087 bei n = 200 / 400 / 600 (k49).
- C_E und s_E misst V. Bei s_E − γ ≈ 2 ist 400^{−2} = 6e-6; die Bedingung hat also grosse Reserve, solange C′ ≳ 1e-6.

### B.4. Bausteine und ihr Status

| Baustein | Inhalt | Fundstelle | Status |
|---|---|---|---|
| P-V, p32a/(K3), (K2) | Kneser streng; f_K(L) = L, f_K′(L) = √L; Schwanz C₆ = 461.36 | ARCH §9.8–9.9; V §2.2 | **V** |
| Zertifikat 1 → 400, Kopplung (A) | Zweig der σ-Homotopie streng | ARCH §9.7; V §2.1 | **S** |
| H1, H2 | nichtlineare Kette (Existenz, Eindeutigkeit, Zweig) | ARCH §4.2–4.3 | **V\|(H)**; Umstellung auf K_{r₀} **B** |
| Satz K5 / Rang-1-Ränderung | J̃_n⁻¹ = Σ c_mℓ_mᵀ/ℓ_mm; ΔK_n = δ_nΛ_n/ℓ_nn | k_kern §11.3, §24.1 | **B**; inhaltsgleich mit G9/G22 (**V** 20:25) |
| Lemma K13 | Kopfentwicklung mit beliebigen (kopfweisen oder festen) Konturen | §22.1 | **B** (elementar) |
| Satz K19 | Datenseite ohne abgeschnittene Pole: (1 + cΣQ_n)·dist_E | §30.2 | **B** (elementar) |
| Lemma K10.1, Satz K10, Satz K11/K11′ | Eckprofile F_β explizit; Eindeutigkeit (Mellin-Liouville) | §19 | **V** (08:27, mit Ks Vorbehalten) |
| Integraldarstellung | F_β(u) = Γ(β)⁻¹∫₀^∞t^{β−1}e^{−t−ut^{−iα}}dt; Potenzsektor per Watson-Lemma | §27.3; k48 (1e-45) | **B** + num |
| Spiralfenster | Ränder der Nicht-Wachstumsfenster sind n-unabhängige Log-Spiralen (1.623°/e-Faktor) | §28.2 | **B** (aus der Kegelform; die Kegelform selbst V) |
| Lemma K20 | n-te-Wurzel-Rate der Primalen ≡ 0 auf D̄ (Oberhalbstetigkeit + Innengesetz) | §31.1 | **B** (Folgerung) + num |
| Satz K16 | Bootstrap: (ρ1), (ρ2) + glatte Zuwächse ⇒ Gesetz von P_n | §27.1 | **B** (bedingt) |
| O1 | keine superharmonische Vergleichsfunktion mit der Singularität der Köpfe (Laplace-Masse) | §33.2 | **B** (eine Zeile) |
| Eckkonstante C | C = −λ̃_PΓ(½+iα)/((μ+1/μ)Γ(½+2iα)) | §25 | Skizze + num (1e-5) |
| T1 | primale Γout-Zone im Modell real ⇒ r₀ > r_* nötig | V §3.14 | num (V) |

### B.5. Prüfstellen für V (Fassung 2)
1. **K19 (§30.2):**
   - Gilt die Identität Π_Nf − q_N = J_N⁻¹T_NJ(f − q_N)?
   - Ist Lemma K13 mit der festen Kontur 𝒞 durch L zulässig? Dazu gehören: Analytizität von W = J(f − q_N) im Innern von 𝒞, die Lage des Spiralschnitts und die Zweigfortsetzung über die Blätter bei sehr kleinem r.
2. **Geometrie von E:** Lücke um den Schnitt; f_K(𝒞) ⊂ Holomorphiegebiet; m = eˣ/f_K′ beschränkt auf 𝒞.
3. **Messung s_E** (BRETT 10:48, Punkt (2)): Hängt das Budget γ < s_E?
4. **Umstellung H1/H2 auf K_{r₀}:** Kompaktheit von 𝒵 in der K_{r₀}-Norm; H2 ohne sup_{D̄}.
5. **F3:** Ist die Datenbedingung auf 𝒞 für die Fast-Lösungen messbar? Vorschlag: sup_𝒞|p_N − f_K| bis N = 2000 — Linse, Spiralen, Γout-Bogen.

**Prüfvermerk V (12.10. 11:42; Einzelheiten v_pruefstand/BERICHT §1.11, §3.19, §3.20):**
- 1 (K19): richtig. Identität exakt; K13 mit fester Kontur zulässig (Cauchy für Jordangebiete, Spiralen rektifizierbar, f_K in L stetig); die Fortsetzung über
  den Spiralschnitt bei sehr kleinem r ist mit Ps oberer Struktur nachgerechnet (Im z steigt pro e-Faktor um +0.708 − 0.0048, singulärer Teil bleibt O(|t − L|^{5.447})).
- 2 (E): E ⊃ D̄ ∖ U_{r₀}(L, L̄) statt E ⊃ D̄ (wie Ks Korrektur 11:20). Aussenwinkel bei L ist (u₁ − u₂ − 38.3°) = 131.7°, nicht 99° (die 25°-Kante war f(D̄)).
  Jenseits r₁ = 0.05 verengt der Verbindungsbogen zur Linse die Lücke auf ≈ 97°; für die Konstanten bei N₀ = 400 besser r₁ ≈ 0.2. f_K-Geometrie: v55 läuft.
- 3 (s_E): gemessen 3.9–4.1 (80°/−90°, r₁ = 0.2) und 3.3–3.4 (60°/−90°) — Satz K22 bestätigt; mit r₁ = 0.05 bei N ≤ 300 nur 2.9–3.3. Budget γ ≈ 0.9 trägt.
- 4 (H1/H2): Kompaktheit trägt (P_N endlichdimensional, sup_{K_{r₀}} Norm, 𝒮_N abgeschlossen). H2 ohne sup_{D̄}: im Wortlaut von ARCH §4.3 noch nicht geprüft.
- 5 (F3): gemessen — auf der Linse s = 0.03 wächst |p_N − f_K| wie e^{+0.015N} (N = 2000: 1.75), am Steg 1e-4; Γout-Bogen und Σ₂ konvergieren. **F3 ist mit
  dieser Linse für N ≳ 1500 falsch**; die Familie braucht eine andere Datenbedingung oder eine Kontur ohne Linsenteil, auf dem p_N divergiert.
- Fassung 1, Prüfstelle 1 (Faktor 2): Schreibfehler in ARCH §9.7 (1); 12.7 enthält beide Schlitze; für K′ 8.2, nicht 16.3.

### B.6. Für den Autor
Wenn drei Aussagen über das lineare Problem stimmen, nähern sich Deine Polynome p_N Knesers Funktion überall auf der Scheibe bis |L|, ausser in zwei winzigen festen Kreisen um L und L̄. Daraus folgt, dass alle Koeffizienten und p_N(1) gegen Knesers Werte gehen.
- Die erste betrifft zwei Antwortfolgen an der Ecke, deren Form wir explizit kennen.
- Die zweite ist eine Aussage über exponentielles Wachstum, wie sie die Potentialtheorie behandelt.
- Die dritte ist klassische Approximationstheorie.

Dazu kommt, dass dieselben Aussagen auch für die Nachbarpolynome im Newton-Schritt gelten müssen (die „Familie“). Das ist noch offen.

Alles andere ist bewiesen oder bis Grad 400 vom Computer streng nachgerechnet. Die ersten beiden Aussagen gelten in den Rechnungen (bis Grad 600 bzw. 2000); die dritte misst V gerade.
