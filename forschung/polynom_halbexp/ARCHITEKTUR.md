# Beweisarchitektur — Polynom-Wurzel von exp → Kneser (konsolidiert, alle Stränge)

Stand: 10.10.2026, 22:33, ergänzt 23:03 um C §13, D §11, E Runde 3, A 23:02 und b46 (Zeiten aus `date`). Verfasst von Strang B im Auftrag der Leitung (Runde 5), aus den BERICHTen a–e,
dem BRETT und README. Kennzeichnung: **bew.** = bewiesen · **bew.|X** = bewiesen unter Annahme X · **Sk.** = Skizze (Lücke
genannt) · **num.** = numerisch · **Verm.** = Vermutung. Fundstellen: `A §x` = a_moden/BERICHT.md, `B §x` = b_reformulierung/…,
`C §x` = c_literatur/…, `D §x` = d_kantorovich/…, `E §x` = e_modelle/…, `R` = README.md.

---

## 0. Kurzfassung

**Hauptsatz (bew.|(K),(H),(A); Beweis §4).** Gilt Knesers Struktursatz (K) und die *eine* lineare Familien-Annahme (H) (§5.1),
so hat für alle N ≥ N₀ das abgeschnittene System T_N[p∘p − eˣ] = 0 in der Kugel
B_N = {p : sup_{D̄}|p − f_K| ≤ R_N}, R_N := K₀·C'·N^{κ−4.69}, **genau eine** Lösung p*_N, und sup_{D̄}|p*_N − f_K| ≤ C'N^{κ−4.69}.
Liegt der fortgesetzte Zweig bei einem N₀ in B_{N₀} (Startbedingung (A)), so ist **p_N = p*_N für alle N ≥ N₀**. Insbesondere
gilt **p_N → f_K gleichmässig auf D̄** (also auf D), mit Rate N^{κ−4.69} (gemessen κ ≈ 0: Rate ≈ N^{−4.8} am Rand, N^{−7.9} im Inneren).

**Die minimalen offenen Annahmen** (alles andere in der Kette ist bewiesen):

| | Annahme | Art | numerisch gestützt |
|---|---|---|---|
| (K) | Knesers Struktursatz: T(w) = φ(e^{L(w+θ(w))}), θ 1-periodisch, beschränkt für Im w > c | Literatur (Kneser 1950, Paulsen–Cowgill 2017; nicht im Original geprüft [S]) | A: κ_m, Formel am Motor auf 5e−37 |
| (H) | Pol-Fehler der Linearisierungen an Fast-Lösungen sind auf dem Schlitz bei L, L̄ **polynomial** (mit einem Bruchteil ε des Green-Abfalls) und auf einer festen Niveaulinie Γ_d exponentiell klein (= Cs (E)); dazu Invertierbarkeit und Kopf-Stabilität (§5.1) | **eine lineare Aussage** über eine Familie von Matrizen | P = f_K: Eckprofil N ≤ 1200, (a′) N ≤ 400 (D), Pol-Antworten N ≤ 400 (B, C); Familie: N ≤ 1000 (A), ≤ 500 (B) |
| (A) | p_{N₀} ∈ B_{N₀} für ein N₀ | endliche Rechnung | Krawczyk bis 120 (R); Zweig = Lösung bei Π_N f_K bis 300 (A); p_N nahe Π_N f_K bis 1500 (A) |
| (b35b′) | f_K(eᶻ) ∉ (−∞,0] auf dem linken Bogen von ∂D | endliche Prüfung | min|f_K(eᶻ)| = 0.73 (B) |

**Kern der Sache:** (H) eingeschränkt auf P = f_K heisst (H₀). Das ist Strang Cs **Eck-Bedingung (E) = (E_R) ∧ (E_K)**
(C §13.2), mit der bewiesen ist: (F) + (E) ⇒ Π_N f_K → f_K für alle N (C Satz 13.1 = B-Max in lokaler Form). Spektral zerfällt (E) in
zwei lineare Aussagen über J_N(f_K):
- **(a′)** Rest-Kontraktion σ_{r+1}((I + K_N)⁻¹) ≤ 1 + o(1); daraus folgt (E_R) (D Satz D8). Numerisch auf 1e−4 genau (N ≤ 400).
- **(E_L)** Randprofil der adjungierten Eckmoden: σ_i|Ṽ_i(arg L + τ)| ≤ C N^κ e^{λ(τ)N} mit λ ≤ A|τ|^p, p > 1/α = 1.13, und
  λ < 0.068; daraus folgt (E_K) (D §11.3 + C Lemma 13.2). Numerisch p ≈ 1.44–1.52, ∫λ dν = 0.086 gleichmässig (N ≤ 1200).

**(E_L) ist die gesuchte lineare Aussage über die Eckmoden nahe L.** Die exponentiell wachsenden Eckmoden (σ₁ ~ e^{0.012N}) werden von
Polen nahe L nicht exponentiell angeregt: die Anregungsrate verschwindet zur Ecke hin schneller als |τ|^{1.13}. Eine Rate, die zur Ecke
hin nur polynomial abnimmt, reicht beweisbar nicht (C Satz 13.3).

Die Familien-Erweiterung (H₀) → (H) (früher „(H*) ⇒ (H*_n)“) ist offen. Nach B-Id/B-NF ist jede Ordnung wieder Π_N auf explizite
Eingaben (B §9).

**Nicht mehr gebraucht** (und teils widerlegt): gleichmässige Stabilität (a) [falsch, D §4.5, §9.1]; (S)/(H*) „polynomial im
ganzen Band“ [für grosse N falsch, C §12.9, b45]; (G2) [nur für den Nebenweg B-lin]; Montel/Eindeutigkeitssatz aus R [nur für
den Nebenweg]; separate quadratische Kontrolle (Q) [folgt aus (H)]; Hierarchie (H*_n) [in (H) aufgegangen]; Eindeutigkeitskugeln aus
Operatornormen; Kontraktion auf **fester** Kugel (E) [nicht nötig; numerisch bis N = 400 unproblematisch, Restrisiko extrapoliert ab
N ~ 3000, §6.5]. Kugelargumente für Newton selbst scheitern sicher (A 23:02: Zufallsrichtungen werden ≳ 1.17^N verstärkt) — die
Architektur benutzt nur kohärente Richtungen.

---

## 1. Bezeichnungen

f = f_K (Knesers halbe e-Funktion), a₀ = f(0) = 0.49856, L = 0.31813 + 1.33724i, ρ = |L| = 1.37456, R₀ = |L − a₀| = 1.34935,
γ = 2πi/L = 4.44695 + 1.05794i, μ = f'(L) = √L (|μ| = 1.1717). D = D(0,ρ), Ḡ' = W̄ = D̄ ∪ D̄(a₀,R₀) („Zwei-Kreis-Menge“, Ecken
L, L̄ mit Aussenwinkel απ, **α = 0.882964**). g = Green-Funktion von ℂ∖Ḡ' mit Pol ∞ (explizit: Möbius + Potenz, B Lemma 7.3 =
C §12.2; cap = 1.5169); Γ_η = {g = η}. Exponent der Randrate: α(1 + Re γ) = 4.81, α(2 + Re γ) = 5.69, **Schwelle 4.69** = 5.69 − 1.
P_N = Polynome vom Grad ≤ N, T_N = Abschnitt bei 0. J_∞u = u∘f + (f'∘f)u, J_N = T_N J_∞|_{P_N} (= Cs A_f, Ds A_N).
**Π_N u := J_N⁻¹T_N J_∞u** (J-Projektion, B §1.2); k_ζ(z) = 1/(ζ − z); Pol-Antwort Π_N k_ζ; Pol-Fehler E_N(ζ,z) = Π_N k_ζ(z) − k_ζ(z);
Cs y_ζ = Π_N k_ζ − T_N k_ζ (Pol-Antwort ohne „Kopf“), R_N(ζ) = Σ_j ρ^j|(y_ζ)_j|.
Für einen Keim P nahe f: Sekantenoperator 𝒥_P v = v∘P + M_P v, M_P = (f∘P − f∘f)/(P − f); Tangente J_N(P)δ = T_N[δ∘P + (P'∘P)δ];
Π^P_N := (T_N𝒥_P|_{P_N})⁻¹T_N𝒥_P, E^P_N(ζ,z) := Π^P_N k_ζ(z) − k_ζ(z). Bei P = f: 𝒥_f = J_∞, Π^f_N = Π_N.
x₁ := Π_N f_K; ε_N := sup_{D̄}|x₁ − f_K|; t_N := J_N⁻¹x^N (Kopfantwort).

---

## 2. Katalog aller Aussagen

### 2.1 README (Ausgangslage)

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| R1 | Rückführung: (i) sup_N sup_D|p_N| < ∞, (ii) |p_N| ≤ |L| auf Kompakta des Korridors V₀, (iii) p_N(0) ≤ ρ₀ < |L| ⇒ p_N → f_K | bew. (Montel + Max.-Prinzip + R2) | R §2.1 |
| R2 | Eindeutigkeit: jede auf D holomorphe reelle Lösung von f∘f = exp, die V₀ in D abbildet, ist f_K | bew. bis Form von Ω = slog(D) (num., „(G2)“) | R §2.2 |
| R3 | N ≤ 120: eindeutige Lösung in Krawczyk-Box, Σ|a_k||L|^k ≤ 2.2405 | bew. (arb) | R §2.3, e51 |

### 2.2 Strang A — Moden & Kneser-Daten

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| A1 | Singularitätsformel f_K = φ(√L·Ψ·e^{Lε(Ψ^γ)}), ε = 2κ₁q + 4πiκ₁²q² + …; Leitterm d₁t^{1+γ}, d₁ = 2L^{3/2}κ₁, D₁ = −0.03200 − 0.01124i | Sk. bei A; **bew.|(K)** als B Lemma 8.2; num. (Motor 5e−37, Darboux 1e−18…1e−24) | A §4, B §8.3 |
| A2 | Taylor-Koeffizienten a_k = 2Re[D₁L^{−k}k^{−2−γ}/Γ(−1−γ)](1+O(1/k)); Hüllkurve 71.57·k^{−6.447}ρ^{−k} | num. (folgt aus A1 + Flajolet–Odlyzko, nicht ausgeführt) | A §4 |
| A3 | Ein-Moden-Struktur: p_N − f_K = c_Nδ_{−1} + c̄_Nδ_{+1} im Inneren bis rel. 1e−14 | num. (= B §2.2) | A §3 |
| A4 | Amplitudengesetz c_N = K·N·u_N/(δ_{−1})_N, c_N ~ N^{1−2γ}, K∞ ≈ 0.0017 − 0.0533i; Blindtest Vorzeichenwechsel N = 215–218 | num.; K∞ Verm. | A §4b, §4c, R2.3 |
| A5 | max_D|Π_N f_K| = 2.232651 für alle N ≤ 1000; |Π_N f_K(L) − L| ~ N^{−4.85} | num. | A R2.2, R2.4 |
| A6 | Echte p_N bis N = 1000 (Newton von Π_N f_K, quadratisch); p_N − Π_N f_K ≈ 0.75e−3·N·|E_N(L)|²; kein Exponentialanteil (c < 1e−4/Grad); Zweig bis N = 300 = diese Lösung (6e−61) | num. | A R2.6 |
| A7 | (F)-Konstanten: Sprung |Δf|/u^{5.447} ≤ 0.95 (u ≤ 0.2), max|f_K| auf ∂(W+0.1) = 3.32 | num. | A R2.5 |
| A8 | Vermutung S (strukturierte lineare Abschätzung S1/S2) | Verm. (überholt durch (H)) | A §4d |

### 2.3 Strang B — Reformulierung & Beweisstrategie

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| B1.1 | Jet-/Dreiecksform J_N = H(S + Γ); Fixpunktform = Richardson-Iteration mit Spektrum bei −1 | bew. | B §1.1 |
| B1.2 | J-Projektion Π_N: reproduziert P_N, annulliert ker J_∞, hängt nur von Jets bei 0, a₀ ab | bew. | B §1.2 |
| B1.3 | Lineare Eindeutigkeit: ker J_∞ ∩ {|δ| ≤ K|z−L|^{−a}|z−L̄|^{−a}, a < 3.447} = {0} | bew.|(G2) (= D4) | B §1.3 |
| B1.4–1.6 | θ-Gleichung für p_N; Zwei-Konstanten-Schranke; quantitativer Bôcher | bew. | B §1.4–1.6 |
| B-lin | |Π_N f_K| ≤ B(|z−L||z−L̄|)^{−a}, a < 3.447 ⇒ Π_N f_K → f_K lokal glm. auf D | bew.|(G2) | B §1.7 |
| B-Stab | gleichm. Stabilität in A_s ⇒ (LP) | Implikation bew.; Voraussetzung **num. falsch** | B §7.1 |
| L7.1 | ζ ↦ Π_N k_ζ(z) rational, Pole nur in 0, a₀, E_N = O(ζ^{−N−2}) | bew. (= C 12.1) | B §7.2 |
| L7.2 | Bernstein–Walsh in ζ: |E_N(ζ)| ≤ max_{g=t}|E_N|·e^{−(N+2)(g(ζ)−t)} | bew. (= C 12.2) | B §7.2 |
| L7.3 | Green-Funktion explizit, α = 0.883, cap 1.5169, g(L+ue^{iψ}) ≈ 0.655u^{1.1326} | bew. (+num. Kontrolle) | B §7.2 |
| L7.4/7.5 = L8.2–8.4 | (F): f_K holomorph um Ḡ' bis auf zwei Schlitze, Sprung ≤ C_J|ζ−L|^{5.447} | **bew.|(K)**, eine Prüfung num. (b35b′), Rest arb-streng (b35) | B §8.3 |
| B-Max | (F) + (H*) ⇒ sup_{D̄}|Π_N f_K − f_K| ≤ C N^{κ−4.69} | bew.|(H*); lokale Form (H₀) genügt (§4.1) | B §7.3, §9.3c |
| B-Sek | p Lösung ⇒ p = Π^p_N f (Sekanten-Projektion); (H*) für Π^p ⇒ sup|p − f| ≤ CN^{κ−4.69} | bew. (a posteriori) | B §8.1 |
| L8.1 | multilineare Pol-Kerne rational in jeder Variablen | bew. | B §8.2 |
| B-Q | (H*_n^{loc}) ∀n ⇒ Lösungen p*_N → f_K | bew.|(H*_n^{loc}) (in §4 durch (H) ersetzt) | B §8.2, §9.3c |
| B-Id | D²Φ(f)[a,b] = −J_∞[(ab)'/(2f')] + Terme mit J_∞a, J_∞b (exakt, beliebiges f) | bew. (symbolisch geprüft b40) | B §9.1 |
| Kor. B-Ξ₂ | Ξ₂(ζ,ζ') = ¼Π_N[(E_ζE_ζ')'/f'] − Kopfterm·t_N | bew.; num. (Kopfterm ≤ 1.4e−9) | B §9.1–9.2 |
| L9.1 | J_∞a = 0 ⇒ W = a/(2f') exp-invariant; τ_t⁻¹∘f∘τ_t exakte Lösungen | bew. | B §9.4 |
| B-NF | formale Normalform S(a) = a + aa'/(2f') + …, Φ(f+S(a)) − eˣ = R(a; J_∞a) | Sk. | B §9.4 |
| B-Zweig | λσ-Fortsetzung ⇒ p_N = p*_N | bew.|(A),(B),(C) von B §9.5 (in §4.3 vereinfacht) | B §9.5 |
| B9.5(i–iii) | F_{N+1}(pad p_N) = c_{N+1}x^{N+1}; x₁(N+1) − pad x₁(N) = −ρ^f_N t_{N+1}; Zweigschritt = −c∫J⁻¹x^{N+1} | bew.; num. auf 1e−22 (b43) | B §9.5 |
| Bnum | Ordnungen: |p_N − x₁| = c₂ε², c₂ ~ N^{0.94}; Rest ~ c₃ε³, c₃ ~ N^{1.82} (N ≤ 500); Eckpaar-Kerne ~N^{0.92} (N ≤ 400) | num. | B §8.4, §9.2–9.3 |

### 2.4 Strang C — Literatur & Analogien

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| C-Lit | Kein Konvergenzbeweis für irgendein Abschneideverfahren gegen Kneser bekannt (Walker 1991, Paulsen 2016, Trappmann 2011, TPID 18) | Literatur [O] | C §0, §2 |
| C-Abel | Linearisierung in Abel-Koordinaten ≅ u ↦ u(w+½) + u(w) | bew. (formal) | C §1 |
| C-Trans | Transferoperator-Rahmen bei 0 nicht anwendbar (kein f⁻¹-invariantes U ∋ 0) | bew. | C §2.4 |
| C9.1–9.4 | Schwanzkopplung; F_N(q_N) = −B_N t_N; Pol-Kopplung geschlossen; Cauchy-Darstellung | bew. (9.1 = D1) | C §9.2 |
| C9.5 | (F) + (H_R1) + (H_R2) ⇒ α_N = |J⁻¹F(q_N)| ≤ C N^{−4.81} | bew.|(F),(H_R) | C §9.3 |
| C10.1 | Exakte Zerlegung des nichtlinearen Rests N_q(δ) = [Q_f(δ)] − [Q_t(δ)] | bew. | C §10.1 |
| C10.4 | Sekantensatz p − q = Ā⁻¹B_N t_N | bew. (a posteriori) | C §10.4 |
| C-Q | Vermutung Q (starke Konvergenz auf regulären Eingaben) | Verm.; für W̄-holomorphe v bew.|(S^f) (C12.5) | C §11, §12.7 |
| C12.1 | y_ζ rational in ζ, Pole nur 0, a₀ | bew. (= L7.1) | C §12.1 |
| C12.2 | R_N(ζ) ≤ S_N e^{−(N+2)g(ζ)} | bew. (= L7.2) | C §12.2 |
| C12.3 | (F) + (S) [S_N = max_{∂W}R_N = O(N^b), b < 5.69] ⇒ α_N ≤ C S_N N^{−5.69} | bew.|(F),(S) | C §12.3 |
| C12.5 | Q für auf W̄ holomorphe v, exponentielle Rate | bew.|(S^f) | C §12.7 |
| C12.6 | Poisson–Jensen: log R_N(ζ) ≤ ∫_{∂W}log R_N dω_ζ − (N+2)g(ζ); explizite Schwellen 0.065/0.068, 14.65 | bew. | C §12.9 |
| Cnum | S_N = R_N(L) ≈ 0.33N^{1.12} (N ≤ 400); Spitze auf fernem D-Bogen ~1.0115^N, überholt die Ecke bei N ≈ 420 | num. | C §12.4, §12.9 |
| C13.0 | Eckzerlegung erhält Rationalität und Abfall; richtiges Eckfunktional ist ⟨u_i, y_ζ⟩ (nicht ⟨u_i, Π_N k_ζ⟩) | bew. | C §13.1 |
| C13.1 | (F) + (E) ⇒ α_N = O(N^{κ−5.693}) für alle N; (E) = (E_R) ∧ (E_K), (E_K) nur auf dem Schlitz S_L(d) und der Niveaulinie Γ_d, mit Bruchteil ε der Green-Rate | bew.|(F),(E) | C §13.2–13.3 |
| C13.2 | Randform: ω_ζ(λ_N) ≤ (1−ε)g(ζ) auf Schlitz ∪ Γ_d ⇒ (E_K) (vorzeichenbehaftete Randrate) | bew. | C §13.4 |
| C13.3 | Schärfe: für einen Eckanteil vom Rang 1 ohne Nullstellen ist die Randform auch notwendig; feste Rate mit nur polynomialer Eckunterdrückung reicht nicht | bew. | C §13.4 |
| C13.4 | einteilige Form: R_N(L) ≤ CN^κ und Θ_N := sup ω_ζ(λ_N)/g(ζ) ≤ 1−ε ⇒ α_N → 0 (ohne Zerlegung, aus Randdaten messbar) | bew.|(F) | C §13.5 |
| C13.5 | kohärente Minimalform: (F) + (E_R) + |Q_N δ1| → 0 ⇒ α_N → 0 | bew. (trivial) | C §13.5 |

### 2.5 Strang D — Kantorovich / Induktion / Eckmoden

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| D1 | Schwanz wirkt linear: [p∘p]_{≤N} = F_N(p_lo) + B_N(p_lo)a_hi | bew. | D §1 |
| D2 | Bordering: System N+1 ⇔ {F_N + t·b = 0, Quadratik in t} | bew. | D §1 |
| D3 | J_N = U_N S_N + M_N (Taylor-Verschiebung) | bew. | D §4.1 |
| D4 | Injektivität der Grenz-Linearisierung (ker J_∞ trivial unter Wachstumsschranke) | Sk.|(G2), f' ≠ 0 num. (= B1.3) | D §4.9 |
| D5 | J⁻¹P_N[h∘x₀ + (x₀'∘x₀)h] = P_N h + J⁻¹B_N h_hi ⇒ Eckanteil kohärenter Eingaben = σ_i⟨v_i, B_N h_hi⟩ | bew. | D §9.4 |
| D6 | Eckformel ⟨u_i, Π_N k_ζ⟩ = σ_i[m(ζ)V_i(ζ)/ζ + Σ_j V_i(x_j)/(x_jx₀'(x_j)) + ℓ_i(P_N Q_ζ)] | bew. (polynomiales x₀) | D §10.1 |
| D7 | Zerlegung J⁻¹ = Σ_{i≤r}σ_iu_iv_iᵀ + R_N ⇒ (H*) = Rest-Teil (‖R_N‖·Eingabe) + Eck-Teil | bew. (für Pole mit Urbildern ausserhalb D̄) | D §10.2 |
| Dnum1 | Grenzoperator: im Inneren Toeplitz mit Sprung-Symbol (f'e^{−x})·1_{Γin}; (a) in ℓ¹ falsch (log) | num. + Argument | D §4.2–4.3 |
| Dnum2 | Eckmoden σ₁,₂ ≈ C N^{−a} e^{0.0120N} bis N = 1000; Ursache: abstossender Randfixpunkt mit Drehung; kubisches Modell reproduziert sie | num. | D §9.1–9.3 |
| Dnum3 | (Rest) σ_{r+1} = 1.3663 = ‖σ‖_∞ bis N = 1000 (r ≤ 3) | num. | D §10.3 |
| Dnum4 | (E_D) Eckanregung λ(t) ≈ 0.04t nahe L, ∫λdν ≈ 0.11 (N = 1000) | num. | D §10.4 |
| Dnum5 | Stabilität bei p_N, nicht bei q_N; Kantorovich um q_N hoffnungslos | num. | D §4.7 |
| D8 | J_N = M_N(I + K_N), Kδ = (δ∘f)/(f'∘f); S_N = σ_{r+1}(J_N⁻¹) ≤ σ_{r+1}((I+K_N)⁻¹)·‖σ‖_∞ | bew. | D §11.1 |
| (a′) | σ_{r+1}((I+K_N)⁻¹) ≤ 1 + o(1) | num. (1.0001–1.0002, N ≤ 400, kubisch); Akkretivität und reine Toeplitz-Parametrix scheitern | D §11.2 |
| D-(b) | (E_D) ⇐ Randprofil der adjungierten Eckmode σ_i|Ṽ_i(arg L+τ)| ≤ CN^κe^{λ(τ)N}, λ ≤ A|τ|^p (p > 1.13), λ < 0.068 | Reduktion bew. (Pole mit arg(ζ−L) ∈ (24.9°, 166.6°)); Profil num.: p ≈ 1.44–1.52, ∫λdν = 0.086 (N ≤ 1200) | D §11.3 |

### 2.6 Strang E — Modellprobleme

| ID | Aussage | Status | Fundstelle |
|---|---|---|---|
| E-Ers | Ersatzformen (A)/(B); Abel-Darstellung der Linearisierung | bew. | E §1 |
| E-P | Modell P (Fixpunkt im Entwicklungspunkt): p_N = T_N h; parabolisch Divergenz | bew. | E §2 |
| T | Translationsmodell: u_N = φ_s(D)T_N h; Konvergenz ⇔ Typ < π/s; Fehler bei gleichem Typ ~ N^{1−Lücke} | bew. | E §3 |
| K | kontraktiver Fall: |J_N⁻¹| glm. beschränkt, p_N → f* geometrisch | bew. (N ≥ N₀) | E §4 |
| K′ | g = (3−x)/(8−3x): ganzer Zweig, ‖p_N − T_N f‖₁ ≤ 7.71·6^{−N−1} | **vollständig bew. (arb)** | E §4 |
| E-Möb | Möbius-Drehung: Rate e^{−g(ζ_p)}, Schwelle c = π/3 ⇔ Pol in Ḡ' | num. (+ M2, M3 bew.) | E §5, R2 |
| M1 | Π_N f_M − f_M = r·E_N(ζ_p,·) | bew. | E R2.1 |
| M2 | Cayley-Bild; ζ_p ∈ Ḡ' ⇔ c > π/3 | bew. | E R2.1 |
| M3 | Green-Funktion explizit; e^{−g(ζ_p)} = sin(πc/2ω)/sin(π(π−2c)/2ω) | bew. | E R2.1 |
| M4 | p löst System ⇔ p = Π^p f; p − f_M = r·E^p_N(ζ_p,·) | bew. | E R2.1 |
| M-lin | (H*)_M ⇒ (LP) im Möbius-Modell, Rate e^{−Ng(ζ_p)} | bew.|(H*)_M | E R2.2 |
| M-voll | (H*_sek) auf fester Kugel ⇒ Lösung (Brouwer) | bew.|(H*_sek) | E R2.2 |
| Enum | (H*)_M mit C ≈ 1, κ = 0 bis N = 160; Möbius ohne Eckmode; exp: T(p) = Π^p f_K Kontraktion auf fester Kugel (N ≤ 160, Zufallsrichtungen; B b46 bis N = 400 auch in Eckrichtung) | num. | E R2.3 |
| E-R3 | Möbius: det J_N = (1+τ²)^N Q_N(τ²)/(1−τ²)^{2N+2}, Q_N irreduzibel (N ≤ 5) ⇒ keine geschlossene Form; J_N ≅ λ(I + λ⁻¹𝒮_N) (komprimierte Drehung), Schwelle pseudospektral; Residuen-Lemma für die Negativrichtung | bew. (exakt) bzw. bew. für N ≤ 5; Schwelle num. | E „Runde 3“ |
| E-Prüf | Drehung um L, L̄ (Knesers W-Geometrie, keine Eckmode): (H*) mit C = 1.01, κ = 0 (N ≤ 160; symmetrisch bis 320) — Prüfstein für den regulären Teil (E_R) | num. | E „Runde 3“ |

---

## 3. Abhängigkeitsgraph

Pfeile „X, Y ⟹ Z“ sind **bewiesen**, wo nicht anders vermerkt.

```
(K) [Lit.] ──B L8.2──▶ Singularitätsformel (A1) ──B L8.3──▶ Sprung ≤ C|ζ−L|^{5.447} ┐
(K) + b35 (arb) + b35b′ (num.) ──B L8.4──▶ Holomorphie um Ḡ' bis auf Schlitze ────────┴─▶ (F)

L7.1 = C12.1 (Rationalität)  +  L7.2 = C12.2 (Bernstein–Walsh)  +  L7.3 (Green explizit)
        │                                   │
        │          C12.6 (Poisson–Jensen) ──┤
        ▼                                   ▼
(a′) [Rest-Kontraktion] ──D8──▶ (E_R)          (E_L) [Profil der adj. Eckmoden] ──D6 + Poisson aussen + C13.2──▶ (E_K)
(E) := (E_R) ∧ (E_K)  ≡  (H₀) [lokale Pol-Antwort an P = f, nur auf Schlitz ∪ Γ_d]
C13.4: (einteilige Randform) R_N(L) ≤ CN^κ, Θ_N ≤ 1 − ε ──▶ (H₀) ohne Zerlegung
(S) ⇔ (H*) [polynomial im ganzen Band; num. FALSCH ab N ≈ 420–650] ──▶ (H₀)

(F) + (H₀) ──B-Max (lokal, §4.1) = C12.3/12.6──▶ ‖Π_N f_K − f_K‖_{D̄} ≤ C N^{κ−4.69}   [lineare Vermutung (LP)]
(F) + (H) ──B-Sek + Bootstrap + Überlagerung (§4.2)──▶ ∃! p*_N ∈ B_N, ‖p*_N − f_K‖ ≤ C N^{κ−4.69}
(H) + (A) ──§4.3 (σ-Homotopie, Eindeutigkeit in B_N)──▶ p_N = p*_N ∀ N ≥ N₀  ──▶ VERMUTUNG

Nebenwege:
(H₀) + (H₂^{loc}) + … ──B-Q──▶ p*_N (analytische Reihe)            [ersetzt durch §4.2]
B-Id, B-NF (Sk.) ──▶ „jede Ordnung von (H) = Π_N auf explizite Eingaben“  [Struktur für (H₀) ⇒ (H)]
(F) + (H_R1,2) ──C9.5──▶ α_N → 0 (Konsistenz am Zentrum q_N)
p_N Lösung + Ā invertierbar ──C10.4 / B-Sek / M4──▶ „Vermutung ⇔ lineare Schranke für den Sekantenoperator“
(LP-B) + (G2) ──B-lin──▶ (LP);  R1 + R2 + (G2) ──▶ Vermutung aus Schranken (i)–(iii)
(H*_sek auf fester Kugel) ──M-voll/E (Brouwer)──▶ Existenz      [bis N = 400 auch in Eckrichtung erfüllt, b46; §6.5]
```

---

## 4. Hauptsatz und die zwei neuen Beweisschritte

### 4.1 B-Max in lokaler Form (bew.|(F),(H₀))

Sei d > 0 so, dass f_K ausserhalb Ḡ' bis Γ_d holomorph und beschränkt ist (bis auf die Schlitze; A: max|f_K| auf ∂(W+0.1) = 3.32;
C: für d = 0.1 liegt Γ_d in diesem Bereich), S_L(d) = {L(1+u) : 0 < u ≤ u₀(d)} der Schlitz bis Γ_d.
**(H₀)** (= Cs (E) in Pol-Fehler-Form): es gibt ε ∈ (0,1] mit
  |E_N(ζ,z)| ≤ C N^κ e^{−εN g(ζ)} / dist(ζ,Ḡ')   für z ∈ D̄ und ζ ∈ S_L(d) ∪ S_L̄(d) ∪ Γ_d
(auf Γ_d also ≤ C N^κ e^{−εNd}/dist). Mehr verlangt der Beweis nicht. Eine Schranke im ganzen Band mit e^{cN}-Zuschlag fern der Ecke
(meine frühere Form) impliziert das mit ε = 1 − c/d.
**Aussage:** sup_{D̄}|Π_N f_K − f_K| ≤ C'(ε^{−4.81}N^{κ−4.81} + N^κ e^{−εNd}).
*Beweis.* Cauchy-Darstellung über Γ_d und die Schlitzufer (B §7.3 (1), C §13.3 (2)). Auf den Ufern ist dist ≍ u und
∫₀u^{5.447}N^κ e^{−εNc_g u^{1/α}}u^{−1}du = O(N^κ(εN)^{−α·5.447}) = O(ε^{−4.81}N^{κ−4.81}); auf Γ_d ist |E_N| ≤ CN^κe^{−εNd}/dist. ∎
(Die Schwelle ist hier 4.81, B-Max hatte 4.69 wegen des Umwegs über die Niveaulinie g = 1/N; im Folgenden bleibt 4.69 als sichere
Schranke stehen.) (Dasselbe gilt wörtlich für jeden Sekantenoperator Π^P_N mit
P(0) ∈ int Ḡ', weil L7.1/L7.2 nur Rationalität mit Polen in 0 und P(0) brauchen — B-Sek, E M4.)

### 4.2 Satz H1 — Existenz, Eindeutigkeit, Konvergenz (bew.|(F),(H-i),(H-ii); neu)

Sei 𝒮_N die Menge der „Fast-Lösungen“ P = (1−λ)f + Q (λ ∈ [0,1], Q ∈ P_N) mit T_N[P∘P − eˣ] ∈ ℂ·x^N und
B_N := {P : sup_{D̄}|P − f| ≤ R_N}, R_N = K₀C'N^{κ−4.69}. Unter (H) (§5.1) gilt für N ≥ N₀: die Lösungsmenge
𝒵 := {(λ,Q) : T_N[Φ((1−λ)f + Q) − eˣ] = 0, P ∈ B_N} ist der Graph **einer** stetigen (analytischen) Funktion λ ↦ Q(λ) auf [0,1];
p*_N := Q(1) ist die **einzige** Lösung des abgeschnittenen Systems in B_N, und sup_{D̄}|p*_N − f| ≤ C'N^{κ−4.69}.
*Beweis.* (1) *A-priori-Schranke:* Für (λ,Q) ∈ 𝒵 ist Φ(P) − Φ(f) = 𝒥_P(P − f) (Sekanten-Identität), also T_N𝒥_P(Q − λf) = 0,
d. h. Q = λΠ^P_N f und P − f = λ(Π^P_N f − f). Nach 4.1 (für Π^P) ist sup_{D̄}|P − f| ≤ λC'N^{κ−4.69} ≤ R_N/K₀ — strikt im Inneren
von B_N. (2) *Lokal:* D_Q T_N[Φ(P) − eˣ] = J_N(P) ist invertierbar (H-i) ⇒ Satz über implizite Funktionen: 𝒵 ist nahe jedem Punkt
ein Graph über λ; wegen (1) bleiben die Nachbarlösungen in B_N. (3) *Global:* 𝒵 ist kompakt (abgeschlossen, und Q ist durch
sup_{D̄}|Q| beschränkt), die Projektion π: 𝒵 → [0,1] ein lokaler Homöomorphismus, also eine Überlagerung ihres (offenen und
abgeschlossenen) Bildes [0,1]. (4) *Ein Blatt:* Über λ = 0 ist P = f + Q, und T_N𝒥_P Q = 0 erzwingt Q = 0 (H-i). Also hat die
Überlagerung genau ein Blatt; über λ = 1 liegt genau ein Punkt p*_N. ∎
*Bemerkungen.* (i) Es wird **keine** Operatornorm von J_N⁻¹ benutzt (die Eckspitze stört nicht). (ii) (H) wird nur an Lösungen
in B_N gebraucht — nicht auf einer festen Kugel aller Polynome (anders als E M-voll). (iii) p*_N ist reell (Konjugation erhält 𝒵).
(iv) Die Kugel ist rein durch f_K definiert; die Eindeutigkeit gilt unabhängig davon, wie man p_N berechnet.

### 4.3 Satz H2 — Zweig-Identifikation (bew.|(F),(H),(A); vereinfacht gegenüber B-Zweig)

„Stetig von p₁ fortgesetzt“ wird präzisiert als σ-Homotopie: von pad(p_N) aus löst man T_{N+1}[P∘P − eˣ] = (1−s)c_{N+1}x^{N+1},
s: 0 → 1 (Newton ab pad(p_N) ist das numerische Verfahren dafür; b43: erster Newton-Schritt = linearisierter Homotopieschritt
auf 1e−22). **Aussage:** Ist p_{N₀} ∈ B_{N₀} (A), so gilt p_N = p*_N für alle N ≥ N₀.
*Beweis (Induktion, M = N+1).* P₀ := pad(p_N) ∈ 𝒮_M (Kopf-Rest c_M x^M), sup|P₀ − f| ≤ C'N^{κ−4.69}. Exakt (Sekante an P₀):
c_M·(T_M𝒥_{P₀})⁻¹x^M = P₀ − Π^{P₀}_M f, also |c_M|·sup_{D̄}|(T_M𝒥_{P₀})⁻¹x^M| ≤ 2C'N^{κ−4.69} (4.1 auf Π^{P₀}_M). Längs der
Homotopie ist P(s) − P₀ = −c_M∫₀^s J_M(P(s'))⁻¹x^M ds', und nach (H-iii) |…| ≤ 2·2C'N^{κ−4.69}. Also bleibt P(s) in B_M
(für K₀ ≥ 6), J_M(P(s)) ist invertierbar (H-i), der Pfad erreicht s = 1, und P(1) ∈ B_M löst das System M ⇒ P(1) = p*_M (H1). ∎
(Die λσ-Variante B §9.5 braucht stattdessen Analytizität in (λ,σ); beide ergeben dieselbe Lücke (H).)

### 4.4 Hauptsatz

(K) ⟹ (F) [B L8.2–8.4, mit b35b′]; (F) + (H) ⟹ H1; H1 + (A) ⟹ H2 ⟹ p_N = p*_N → f_K gleichmässig auf D̄ mit Rate N^{κ−4.69}. ∎

---

## 5. Die offenen Annahmen — vereinheitlicht und minimiert

### 5.1 Die Annahme (H)

Es gibt C, κ < 4.69, d > 0, ε ∈ (0,1], K₀ ≥ 6, N₀ mit: für alle N ≥ N₀ und alle P ∈ 𝒮_N ∩ B_N
- **(H-i)** J_N(P) und T_N𝒥_P|_{P_N} sind invertierbar;
- **(H-ii)** (lokale Pol-Antwort = Cs (E) für Π^P_N) für z ∈ D̄ und ζ ∈ S_L(d) ∪ S_L̄(d) ∪ Γ_d:
  |E^P_N(ζ,z)| ≤ C N^κ e^{−εN g(ζ)}/dist(ζ,Ḡ');
- **(H-iii)** (Kopf-Stabilität) sup_{D̄}|J_N(P)⁻¹x^N| ≤ 2·sup_{D̄}|(T_N𝒥_{P'}|_{P_N})⁻¹x^N| für P, P' ∈ 𝒮_N ∩ B_N.

Bei P = f_K ist (H-ii) genau **(H₀)** = die lokale Form von (H*). (H) verlangt dieselbe Schranke für die Linearisierungen an
allen Fast-Lösungen im Ball — das sind **kohärente** Störungen von J_∞ (P − f = λ(Π^P f − f) + Kopf), nicht beliebige.

### 5.2 Wie die Annahmen der Stränge zusammenhängen

| Annahme | Herkunft | Verhältnis zu (H₀)/(H) | Status |
|---|---|---|---|
| (a) sup_N‖J_N⁻¹‖ < ∞ | C R1, D | würde (H₀) geben (B-Stab), ist aber **falsch** (Eckmode) | num. widerlegt |
| (S) = (H*) global polynomial | C §12.3, B §7.3 | ⇒ (H₀); **zu stark** | num. falsch ab N ≈ 420 (ℓ¹, C) bzw. ≈ 650 (sup mit Kopf, aus b45 extrapoliert) |
| (H_R1), (H_R2) | C §9.3 | ⇐ (S) (C12.2); lokale Form ⇐ (H₀) | wie (S) |
| (E) = (E_R) ∧ (E_K) | C §13.2 | **≡ (H₀)** bis auf polynomiale Faktoren (E_N = y_ζ − (k_ζ − T_Nk_ζ), der Schwanz des Kerns ist ≤ e^{−N log(|ζ|/ρ)}/dist); (F) + (E) ⇒ (LP) (C13.1) | num. (über (a′), (E_L)) |
| (E_R) Rest-Teil | C §13.2, D §10.3 | max_{∂W}|(I−Q_N)y_η| ≤ CN^κ; ⇐ (a′) (D8) | σ_{r+1}(J_N⁻¹) = 1.3663 bis N = 1000 |
| (a′) | D §11.1 | σ_{r+1}((I+K_N)⁻¹) ≤ 1 + o(1) ⇒ (E_R) | 1.0001–1.0002 (N ≤ 400) |
| (E_K) Eck-Teil | C §13.2 | |Q_N y_ζ| ≤ CN^κe^{−εNg} auf Schlitz ∪ Γ_d; ⇐ Randform C13.2 ⇐ (E_L) | num. |
| (E_L) = (E_D) Eckmoden-Profil | D §10.4, §11.3 | λ(τ) ≤ A|τ|^p, p > 1.13, sup λ < 0.068 ⇒ (E_K) | p ≈ 1.44–1.52, ∫λdν = 0.086 (N ≤ 1200) |
| C13.4 einteilig | C §13.5 | R_N(L) ≤ CN^κ ∧ Θ_N ≤ 1−ε ⇒ (H₀) ohne Zerlegung | aus Randdaten messbar (C §14 angekündigt) |
| (E_koh) | C §13.5 | schwächste Form, nur Konsistenz, keine Pol-Aussage | |Q_Nδ1| ≤ 3.2e−13 (N = 1000) |
| (F) | C §9.3 | **bewiesen** aus (K) (+ b35b′) | — |
| (K) | Kneser | Literatur | A |
| (H*_n) global | B §8.2 | zu stark (Paare mit linkem Pol wachsen 1.012^N) | num. falsch |
| (H*_n^{loc}), (H*_Π) | B §9.3 | ∀n zusammen ⇔ (H) längs des λ-Pfads (Taylor-Entwicklung von E^P in P − f; B-Id/B-NF: jede Ordnung = Π_N auf explizite Eingabe) | num. n = 2 bis N = 400, n = 3 integriert bis 500 |
| (Q) quadratische Kontrolle | B §0 | folgt aus (H) (H1) | — |
| (H*_sek), (Lip) auf fester Kugel | E R2.2 | stärker als (H) (alle Polynome der Kugel, nicht nur Fast-Lösungen) | num. N ≤ 160 (E), Lip auch in Eckrichtung bis N = 400 (b46); Restrisiko extrapoliert ab N ~ 3000 |
| (G2) | R, B §1.3 | nur für Nebenwege (B-lin, R1/R2, D4) | num. |
| Vermutung Q | C §11 | Spezialfall von (H₀) (holomorphe Eingaben) | bew.|(S^f) |
| Vermutung S | A §4d | durch (H₀) + B-Max ersetzt | — |

**Ergebnis der Minimierung.** Es bleiben **eine lineare Annahme (H)** plus Literatur (K) und endliche Rechnungen ((A), b35b′).
(H) zerfällt sachlich in
1. **(H₀) = (E) an P = f_K** — die lineare Vermutung. Als *eine* messbare Aussage: Cs einteilige Randform C13.4 (Pol-Antwort in der
   Ecke polynomial, und ihr Wachstumsprofil auf ∂W bleibt im harmonischen Mittel um den Faktor 1−ε unter der Green-Funktion). Spektral:
   (a′) [Bulk, D8] + **(E_L)** [Eckmoden]: *die Eckmoden von J_N⁻¹ (abstossender Randfixpunkt L, σ₁ ~ e^{0.012N}) werden von Polen nahe L
   nur polynomial angeregt, die Anregungsrate verschwindet zur Ecke schneller als |τ|^{1/α}* (gemessen |τ|^{1.5}); fern der Ecke genügt
   eine Rate < 0.068 (gemessen 0.012).
2. **Stabilität von (H₀) unter kohärenten Störungen** (P ∈ 𝒮_N ∩ B_N). Ob das aus (H₀) folgt, ist die offene Frage
   „(H*) ⇒ (H*_n)“; B-Id zeigt, dass die erste Ordnung Π_N auf Produkte von Pol-Fehlern ist (keine neue Matrix), und
   die Messung gibt genau eine Eckskalen-Ableitung (N^{≈0.9}) pro Ordnung (B §9.3d).

Eine **einzige** Aussage über Eckmoden nahe L, die alles impliziert, ist mir nicht gelungen: (E_L) liefert (H₀), aber nicht die
Familien-Stabilität. Am nächsten kommt (E_L) **gleichmässig für die Linearisierungen J_N(P), P ∈ 𝒮_N ∩ B_N** (plus (Rest)
für dieselben P); das ist (H) in der Sprache der Eckmoden.

### 5.3 Numerische Stützung im Einzelnen

| Teil | Grösse | Werte (N-Bereich) | Konstante |
|---|---|---|---|
| (H₀) Ecke | dist·sup_{D̄}|Π_N k_ζ| (mit Kopf), t = 0.5/2/8 | 0.17/0.53/0.95 bei N = 400, fallend seit N = 60 (B b25, b41) | C ≈ 1, κ = 0 |
| (H₀) Ecke | R_N(L) = S_N (ohne Kopf, ℓ¹) | 0.33·N^{1.12}, N ≤ 400 (C); sup_{D̄}|y_ζ| = 20.5/33.2/46.6 bei N = 200/300/400 (b45) | κ ≈ 1.1 (ohne Kopf) |
| (H₀) Schlitz | Eckanteil |⟨u₁, Π_N k_η⟩|, η = L + u e^{iψ} | u = 0.001: 1.6 → 2.2; u ≥ 0.003 konstant oder fallend; N = 400 → 1000 (D) | polynomial |
| (H₀) fern | Spitzenanteil (linker Pol) | sup_{D̄}|y_ζ| = 1.5/4.6/15.6 (N = 200/300/400, b45) ⇒ e^{0.012N} | c ≈ 0.012 < η₀ |
| (E_L) | Profil λ(τ) nahe L; ∫λ dν | λ ≈ A τ^{1.44…1.52}, ∫λdν = 0.086 gleichmässig in t → 0 (N ≤ 1200, D §11.3) | Reserve Faktor ≈ 11 |
| (E_R)/(Rest) | σ_{r+1}(J_N⁻¹) | 1.3663 bis N = 1000 (D) | = ‖σ‖_∞ |
| (a′) | σ_{r+1}((I+K_N)⁻¹) − 1 | 1e−4 … 2e−4 (N ≤ 400, kubisch, D §11.2) | — |
| (H) Familie | Kerne 2. Ordnung (Eckpaare) | 1.6e−3 → 9.2e−3 (N = 60 → 400), ~N^{0.92} (b41) | — |
| (H) Familie | Ordnungen am echten p_N | c₂ = 0.09 → 0.41, c₃ = 0.05 → 0.98 (N = 100 → 500, b44); A: |p_N − Π_N f_K|/|E_N(L)|² = 0.75 (N = 1000) | eine Potenz N^{0.9} je Ordnung |
| (H) Familie | Sekanten-Familie fester Kugel | Pol-Antworten unverändert, Lipschitz 1e−9 (N ≤ 160, E); Lip/ε_N = 4.4 (Eckrichtung), 10.6 (zufällig) bei N = 400 (b46) | — |
| (H-iii) | Kopfantwort | T_N = sup|J_N⁻¹x^N|/ρ^N ≈ 0.937^N (4.2e−9 bei N = 300); J(f) ↔ J(pad p_N) gleich auf 1e−10 relativ (b41, b43) | — |
| (A) | Start | Krawczyk bis N = 120 (R3); Zweig = Lösung bei Π_N f_K bis N = 300 (A, 6e−61) | — |
| (H) Familie | echte p_N bis N = 1500 | Newton von Π_N f_K quadratisch (5.7e−33 → 7.1e−62), |p − Π_N f|(L) ≈ 1.1·|E_N(L)|², ω_eff ≈ 1e3 konstant (A 23:02) | — |
| ε_N | sup_{D̄}|Π_N f_K − f_K| | 8.1e−11 (150), 2.7e−12 (300), 6.5e−13 (400), 2.2e−13 (500); |·(L)| 7.7e−15 (1000, A) | C' ≈ 1.3, κ ≈ 0 |

---

## 6. Widersprüche zwischen den Strängen — aufgelöst oder offen

**6.1 Max-Formel (C, „Gesetz R“) gegen Green-Gesetz (B L7.2/7.3, E M3).** C maß R_N^{1/N} → max(ρ/|ζ|, R₀/|ζ−a₀|) mit
„systematischem Vorlauf +0.016“. E zeigt im Möbius-Modell: |E_N(ζ)|·e^{Ng(ζ)} ist konstant (6 Punkte, N = 40–160), die Max-Formel
weicht exponentiell ab (Faktor 411 bei N = 160); die Max-Formel ist die Green-Funktion ohne Eckkorrektur (ω = π statt π − c/2).
Bei exp ist Cs Vorlauf genau log(e^{−g}/Max) (ζ = 2.2: 0.805 gegen 0.793). **Aufgelöst:** richtig ist e^{−Ng}; Cs bewiesene Sätze
(12.2, 12.3, 12.6) benutzen ohnehin g. Nur die numerische „Vermutung R“ ist entsprechend zu korrigieren.

**6.2 Ds Extrapolation „Probleme für Kneser ab N ~ 5000“ (D §4.8).** Beruhte auf „Überlapp ≈ Kopf der Eingabe ~N^{−6.45}“.
Für kohärente Eingaben fällt der Überlapp ∝ 1/σ₁ (D5 + Messung bis N = 1000). **Aufgelöst, von D zurückgezogen** (D §9.5).

**6.3 „(H*) hält mit Konstante 1 bis N = 400“ (B) gegen „Spitze überholt (S) ab N ≈ 420“ (C).** Kein Widerspruch, sondern zwei
Normen. B misst sup_{D̄}|Π_N k_ζ|·dist **mit** dem Kopf 1/(ζ − z), der beim linken Pol in z ≈ ζ den Wert N/2 hat; C misst y_ζ = Π_N k_ζ − T_N k_ζ
**ohne** Kopf (ℓ¹). b45: linker Pol ζ = −(ρ + 2/N): sup_{D̄}|y_ζ| = 1.5, 4.6, 15.6 (N = 200, 300, 400), Maximum bei arg z = 76.5°
(Ds Spitze), Rate 1.0122/Grad; mit Kopf bleibt dist·sup = 0.7517 (Kopf dominiert bis N ≈ 650). **Aufgelöst:** (H*) in globaler
polynomialer Form ist für grosse N falsch (C ab ≈ 420 in ℓ¹, B-Form ab ≈ 650); die lokale Form (H₀) (Ecke polynomial, fern
Rate < η₀) ist das Richtige und reicht (4.1). B §7.4 („κ = 0“) gilt nur bis N ≈ 400 und nur für die gemessene Norm.

**6.4 Cs „p_N konvergiert auf W“ (C §3.5, §9.4) gegen Bs „Π_N f_K konvergiert nicht auf dem Halbmond“ (B §9.2, BRETT 22:22).**
Nachgerechnet (b45, korrigiert): |Π_N f_K − f_K| bei x = 1.6: 1.8e−16, 2.5e−19, 1.5e−20, 2.8e−21, 5.4e−22; bei x = 1.8: 3.0e−12,
6.5e−14, 6.3e−16, 7.7e−17, 7.8e−18 (N = 100…500) — **konvergent im Inneren des Halbmonds**. Am Randpunkt a₀+R₀ = 1.848 wächst der
Fehler dagegen ab N ≈ 300 (b42: 1.5e−11 → 1.4e−10 bei N = 500), ebenso die Pol-Antworten auf dem äusseren Bogen (b41).
**Aufgelöst:** lokal gleichmässige Konvergenz auf dem offenen W, **nicht gleichmässig bis zum äusseren Rand** von D(a₀,R₀). Bs
Formulierung vom 22:22 war zu stark (Korrektur §8). Für den Beweis heisst das: alle Abschätzungen auf D̄ führen, nie Cauchy-Konturen
durch den Halbmondrand (B §9.3a).

**6.5 E: Kontraktion bzw. (H*_sek) auf einer FESTEN Kugel gegen Ds Eckmode.** Befürchtung: Störungen des Operators in Richtung
der Eckmode werden um σ₁ ~ e^{0.012N} verstärkt. Nachgerechnet (b46; Es Sekanten-Abbildung T(P) = Π^P_N f_K exakt, M_P über
f^{(j+1)}∘f ohne Division, Kontrolle T(p_N) = p_N auf 1e−60; P = p_N + ηq, Σ|q_k|²ρ^{2k} = 1, η = 1e−4 … 1e−2):

| N | σ₁ | ε_N | Lip in Eckrichtung u₁ | Lip zufällig | Lip/ε_N (Ecke / zufällig) |
|---|---|---|---|---|---|
| 200 | 1.84 | 1.89e−11 | 2.2e−11 | 9.0e−11 | 1.2 / 4.8 |
| 300 | 3.65 | 2.64e−12 | 7.1e−12 | 1.9e−11 | 2.7 / 7.0 |
| 400 | 8.92 | 6.49e−13 | 2.8e−12 | 6.9e−12 | 4.4 / 10.6 |

Die Ableitung von T ist ∝ ε_N (sie wirkt auf T(P) − f, das kohärent ist, und q tritt nur multiplikativ auf, nie in einer Komposition);
die Eckrichtung ist **nicht** besonders verstärkt, Lip/ε_N wächst langsamer als σ₁ (×1.5 statt ×2.4 je 100 Grade).
**Bis N = 400 aufgelöst:** Es feste Kugel trägt mit Lip ≈ 1e−11. Restrisiko nur extrapoliert: wüchse Lip schlimmstenfalls wie σ₁ε_N,
wäre die Grenze bei N ≈ 3000. Die Formulierung über Fast-Lösungen (§4.2) braucht die feste Kugel nicht. Zum Vergleich A (23:02):
**Newton** verstärkt Zufallsfehler 1e−60 bei N = 1500 auf 1e−24 (≳ 1.17^N, Cs Worst-Case-ω) — dort wirkt F'' mit Komposition auf die
Störung; das betrifft Kugelargumente für Newton, nicht T und nicht §4.

**6.6 Cs frühe „gleichmässige Stabilität in ℓ¹_{|L|,4}“ (C Runde 1) und As „(M4) J_N⁻¹ beschränkt“ (A §4d)** gegen Ds Eckmode.
**Aufgelöst:** von C selbst korrigiert (C §0b); (a) ist falsch, wird aber nicht gebraucht.

**6.7 Bs „(H*_n) gleichmässig im Band“ (B §8.2).** Für n = 2 numerisch falsch (Paare mit linkem Pol ~1.012^N, b41). **Aufgelöst**
durch die lokale Form (B §9.3c) und jetzt durch (H) (§5.1).

**6.8 Faktor N in c_N.** Vier Beschreibungen: Streckband k/N ∈ [1/|μ|,1] (B), Kopfnorm (N+1)/ρ des Pols bei L (C), Schwanzsumme
(E Satz T), Innenschicht N^{−γ}Φ(v) (A). Kein Widerspruch; keine ist bewiesen. Für den Hauptsatz unerheblich (nur κ < 4.69 nötig).

**6.9 Eckenexponent.** Bs früheres β/θ = 4.876 ist korrigiert zu α(1 + Re γ) = 4.81 (B 20:20); gemessen 4.83–4.86 (A, B, C). Offen
ist nur ein Detail: A findet den Imaginärteil des Randexponenten näher an Im γ (1.058) als an α·Im γ (0.934) (A R2.4).

**6.10 Möbius-Modell als Prüfstein.** Es hat |f'(±i)| = 1 und daher **keine** Eckmode (D, E). Es testet Auswahl, Schwelle und
Green-Gesetz, aber nicht die Gefahr aus 6.3/6.5. Das reelle kubische Modell (D §9.3) ist das richtige Modell für (E_L).

**6.11 Definition des Zweigs.** R: „Newton, Start p_{N−1} mit a_N = 0“. Für einen Beweis ist die σ-Homotopie (4.3) die präzise
Form; Newton ist dort das numerische Verfahren (erster Schritt = Homotopie-Linearisierung bis 1e−22, b43; D: linear bis 1e−12).
Offen nur als Definitionsfrage; inhaltlich kein Unterschied in allen Daten (A: Zweig = Lösung bei Π_N f_K bis N = 300).

**6.12 Ds lineares Profil λ ≈ 0.04·t gegen Cs Divergenz-Befund.** C (BRETT 22:39): ein einseitig in t lineares Profil gibt am
Schlitz ~N^{0.12} (formal nicht „für alle N“). D (22:40) misst feiner p ≈ 1.44–1.52 > 1/α (in s: |s|^{1.27…1.34}); damit ist
∫λdν = 0.086 gleichmässig. **Aufgelöst** (numerisch); Cs Satz 13.3 zeigt zugleich, dass die Rate selbst zur Ecke hin verschwinden muss.

**6.13 Eckfunktional.** Ds ⟨u_i, Π_N k_ζ⟩ fällt bei ∞ nur wie ζ^{−1}; für Maximumprinzip/Poisson–Jensen ist ⟨u_i, y_ζ⟩ (Kopf
abgezogen) das richtige Funktional (C13.0). Am Schlitz sind beide fast gleich. **Aufgelöst.**

**6.14 As Rundungswarnung (23:02) gegen „kein Kippen“.** Kein Widerspruch: p_1500 existiert und konvergiert quadratisch, aber nur aus
ungerundetem Start; Zufallsfehler werden exponentiell verstärkt, kohärente nicht (ω_eff ≈ 1e3). Bestätigt die Grundentscheidung der
Architektur (nur kohärente Richtungen, keine Kugel-Normen) und verlangt für (A) als Zertifikat bei grossem N₀ durchgereichte Präzision.

---

## 7. Was als Nächstes (Priorität)

1. **(E_L) beweisen** — die harte Stelle von (H₀). Testbett: das reelle kubische Modell (D §9.3; J_N exakt, Eckmode wie bei Kneser
   auf 3–5 %). Zu zeigen: das Randprofil der adjungierten Eckmode fällt zur Ecke wie e^{−N(c−A|τ|^p)}, p > 1.13. Die transponierte
   Transfergleichung allein ist zirkulär (D §11.3); nötig ist eine Grenzschicht-Beschreibung der Eckmode an der Abschneidekante.
2. **(a′) beweisen:** Parametrix „Toeplitz im Inneren + lokale Mellin-Operatoren an L, L̄ + Eckteil an der Abschneidekante“ für
   I + K_N (D §11.2; die klassische Theorie stückweise stetiger Symbole greift nicht, weil K auf Γout eine Komposition ist).
3. **Familien-Stabilität (H₀) → (H):** mit B-Id/B-NF (Ordnung n = Π_N auf explizite Eingaben) und D5/D6 (Eckanteil kohärenter
   Eingaben). Konkret: Eckanteil von Π_N[(E_ζE_ζ')'/f'] mit der Eckformel D6 abschätzen.
4. **(A) und b35b′ streng:** zertifizierte Kneser-Werte (Intervalle) für f_K auf D̄ bzw. auf dem linken Bogen; dann Kantorovich mit
   Zentrum x₁ bei einem festen N₀ (η = c₂ε² ≈ 1e−24 bei N = 300, also leicht).
5. **Risiko 6.5 weiter beobachten** (Lipschitz der Sekanten-Abbildung in Eckrichtung bei N ~ 1000–2000, mit durchgereichter
   Präzision wegen A 23:02); bis N = 400 unkritisch.

---

## 8. Korrekturen, die diese Konsolidierung an Strang-B-Aussagen vornimmt

- BRETT 22:22 / B §9.0 Punkt 3: „Π_N f_K konvergiert NICHT auf dem Halbmond“ → richtig: konvergiert im Inneren (x = 1.6, 1.8 bis
  N = 500), nicht gleichmässig bis zum äusseren Rand (a₀+R₀).
- B §7.4/BRETT 21:18: „(H*) mit κ = 0“ gilt in der gemessenen Norm (mit Kopf) bis N ≈ 400; der kopffreie Teil wächst beim linken Pol
  exponentiell (b45). Massgeblich ist die lokale Form (H₀).
- Satz B-Q und B-Zweig bleiben richtig, sind aber durch H1/H2 (§4.2–4.3) mit schwächerer Annahme ersetzt.

Dateien zu dieser Konsolidierung: b_reformulierung/b45_abgleich.py (+ b45.log; Teil (2) dort im ersten Lauf mit 15-stelligen
Kneser-Daten — Skript korrigiert, die richtigen Werte in §6.4 stammen aus der Direktrechnung), b46_sekante_eckrichtung.py (+ b46.log),
b41–b44 (Runde 4).

---

## 9. Phase 2: Prüfstand und Korrekturen (Hüter: Strang P; gilt vor §0–§8)

Stand: 11.10.2026, 00:42; Runde 3 eingearbeitet 02:19 (Zeiten aus `date`). Einzelnachweise: `p_pruefung/PRUEFLISTE.md`, Skripte p01–p24. **Massgeblich für die
Phase-2-Stränge ist die korrigierte Zielaussage 9.3.** Statuswörter: **verifiziert** (von P nachvollzogen) · **verifiziert|X** (unter
Annahme X) · **Fehler** · **Lücke** · **offen** · **num.** (nur numerisch gestützt).

### 9.1 Hauptkette (§4) — verifiziert

| Schritt | Status | Bemerkung |
|---|---|---|
| (K) ⇒ (F): L8.2, L8.3, L8.4 | verifiziert|(K) | numerisch bleiben: b35b′ (vom Motor bestätigt, min|f_K(e^z)| = 0.7345) und Radius r > 0.02 der lokalen Formel |
| 4.1 (B-Max lokal) | verifiziert|(F),(H₀) | Bemerkung „ε = 1 − c/d“ falsch auf dem Schlitz (nur via Poisson–Jensen, c < 0.068) |
| H1 (a-priori, Kompaktheit, Überlagerung, ein Blatt) | verifiziert|(H-i),(H-ii) | Identitäten am echten Problem auf 1e-50 (p03) |
| H2 (σ-Homotopie) | verifiziert|(H),(A) | Newton-Definition: jeder Grenzwert in B_M ist p*_M (Eindeutigkeit) |
| (A) | **Lücke (Kopplung), ausgearbeitet in §9.7** | R_N ≈ 240·C·N^{κ−4.81}; (A) = (H) ab N₀ + Zweig-Zertifikat bis N₀ (**rigoros bis N = 193**: p18 bis 45, p24 kohärent ab 45) + zertifizierte f_K auf D̄ auf ≈ R_{N₀} (**jetzt vorhanden: Satz P-V, §9.8**; offen nur noch die Rechnung sup|p_{N₀} − Π_{N₀}f_K|) |
| R3 (Krawczyk N ≤ 120) | verifiziert wie formuliert | zertifiziert isolierte Lösungen, **nicht** den Zweig (Übergänge N → N+1 fehlen) |

### 9.2 Fehler und Lücken in Phase 1

1. **Lücke: (a′) ⇒ (E_R) über D8** (§3, §5.2). Der ganze A-Bogen ∂D(a₀,R₀)∖D liegt in f_K(D); dort hat w_η = [J_∞k_η]_{≤N} den
   Pol h(η) ∈ D und wächst wie (ρ/|h(η)|)^N (1.21^N am Scheitel, p02). D7 setzt selbst „Urbilder ausserhalb D̄“ voraus.
2. **Lücke: (E_L) ⇒ (E_K)** (§3, D-(b)). (E_L) kontrolliert V_i nur auf |x| = ρ; gebraucht werden Werte an Urbildern h(η): auf dem
   A-Bogen liegen sie IN D (dort ist V_i durch (E_L) gar nicht kontrolliert), auf dem D-Bogen knapp ausserhalb — das äussere
   Poisson-Mittel des Profils gibt dort nur eine in |τ| lineare Randrate (nach C 13.4(a) nicht „für alle N“). [Korrigiert 00:37: meine
   p10-Zahl „+0.031“ war float-Rauschen (W 00:40); genau gerechnet sind die Urbildterme klein — die Lücke ist logisch, nicht numerisch.]
   Für Kneser fehlt zudem die Eckformel D6 selbst (bewiesen nur für polynomiales x₀); Os Ersatz O6 braucht bei A_f zusätzlich 9.4 (E_sing).
3. **Fehler (klein):** „(E) ≡ (H₀)“ — bewiesen ist nur (E) ⇒ (H₀).
4. **Fehler (klein):** „höchstens drei Eckmoden“ (D §9.5, §10.3) — r(N) wächst (O, BRETT 00:14: 4 Moden bei N = 800–1000, ~N/250).
5. **Korrektur (X, von P nachgerechnet p08):** Cs „Π_N f_K divergiert im Inneren des Halbmonds ab N ≈ 600–900“ (C 23:05 (5), §6.4) ist
   Rundungsrauschen der 60-Stellen-Datei; echt ist nur das Wachstum am Scheitel a₀+R₀ (b42).
7. **Daten-Fehler (W 00:50, Liste W „01:02“):** Ds d22/d28/d29 benutzten python-flint `acb_series(...).inv()` ohne `ctx.cap` — nur 10 Terme.
   **Falsch:** D §9.4-Tabelle komplett („ganze Antwort“, „Eck-Anteil des Konsistenzfehlers“, Folgerung „~N^{0.7}“; BRETT D 21:21); D §10.4 Zeilen L und
   Schlitz u = 0.001…0.03 (BRETT D 22:20); D §11.3 „Eckanteil in L 2.67/3.00/3.29“ und Schlitzwerte N = 600…1200; §5.3 Zeile „(H₀) Schlitz“ oben.
   Richtig (W w23, kubisch, mit Kopf / kopffrei): L 2.47/3.14/3.65/4.10 (N = 400…1000; kopffrei 3.35…5.62); Schlitz u = 0.003: 1.81/2.00/2.03/2.00;
   u = 0.03: 0.089 → 0.0095 (kopffrei 7.8e-3 → 3.7e-6) — kopffrei fällt der Eckanteil auf dem Schlitz wie e^{−0.73·N·g}. **Nicht betroffen:** D-Bogen-
   Werte (λ(τ)-Tabelle, ∫λdν), d19 (σ, Eckanteil von Π_N f_K 8e-9), d25/d26, d30/d31, Cs c27/c29, Ps p14 (eigene Newton-Inversion, gegen mpmath geprüft).
6. **Numerik-Warnung:** M_P = (f∘P − f∘f)/(P − f) nie per Reihendivision (P − f hat eine Nullstelle sehr nahe 0; p09b) — nur
   M_P = Σ_j (f^{(j)}∘f/j!)(P−f)^{j−1}. Rundung: 60-Stellen-Dateien nur bis N ≈ 1000, darüber durchgereichte Präzision (A 23:02, X 23:53).

### 9.3 Korrigierte Zielaussage (Stand 11.10. 00:42) — massgeblich für W, O, M, X

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
**KORREKTUR (P, Runde 3, nach W 01:06/§14):** Der Exponent 1.26–1.29 ist zurückgezogen (W1 zurückgezogen); das massgebliche Γin-Profil (Eckteil
der kopffreien Pol-Antwort, kubisch UND Kneser) hat **p = 1.5–1.6**. Mit Ws neuem Profil (λ = 0.0001/0.0003/0.0009/0.0016/0.0024 bei
0.5°/1°/2°/3°/4°) und konservativ λ = 0 auf der anderen Seite gibt p17b den Ecken-Grenzwert **K': 0.191 (p = 1.5) / 0.155 (p = 1.6)**, W̄: 0.089/0.080
— beide konvergent (p·α' = 1.18 > 1). **K' ist also NICHT grenzwertig**; der Hybrid K* und (E_A) werden für die Ecke nicht mehr gebraucht.
Der reine K'-Weg (P-II für den Rest, P-III mit p ≈ 1.5 für den Eckteil, bzw. Variante X1) ist damit gleichwertig zum gemischten Weg; X 01:08
misst zusätzlich Abklingen des kanonischen Eckteils auf der f_K(Γout)-Seite (Ecken-Grenzwert −30).

**BEVORZUGT (P 00:54, nach W §10/§13): gemischter Weg — Rest in K', Eckteil in W̄ — braucht KEIN (E_A).** Die Zerlegung
y_ζ = Q_N y_ζ + (I−Q_N)y_ζ besteht aus zwei rationalen Funktionen (Pole 0, a₀, O(ζ^{−N−2}), C 13.0); jede darf ihr EIGENES Gebiet benutzen:
- Rest (I−Q_N)y: Bernstein–Walsh in Ĉ∖K' mit (a′)_r ⇒ (E_R') (Satz P-II, ganz ∂K' liegt ausserhalb f_K(D)) ⇒ |(I−Q_N)y_ζ| ≤ CN^{a+2}e^{−(N+2)g'(ζ)}.
  Keine Eckbedingung nötig (volle Green-Rate).
- Eckteil Q_N y: Poisson–Jensen in Ĉ∖W̄ mit Randdaten auf ∂W̄ aus O6/D6 + **(V1)–(V3)** (W §10): (V1) Profil auf dem D-Bogen, λ ≤ Aτ^p, p > 1/α = 1.1326,
  λ < 0.065 fern; (V2) Aussenwerte ausserhalb D̄ nahe Γin (Linse) UND nahe L auch auf der A-Seite (der A-Bogen verlässt L unter 7.7°, also aussen) polynomial; (V3) **Werte der adjungierten Moden in der Bucht h(A-Bogen) ⊂ D
  polynomial** (σ|V| ≤ CN^κ; numerisch W: Sockel). Dazu (E_sing) und Kantenschwanz für A_f. ⇒ |Q_N y_ζ| ≤ CN^κ e^{−ε(N+2)g_W(ζ)}.
- Zusammen auf Schlitz ∪ Γ'_d: **(H₀^mix)** |E_N(ζ,z)| ≤ CN^κ(e^{−Ng'(ζ)} + e^{−εNg_W(ζ)})/dist(ζ,K'). Satz P-I(i) mit dem langsameren Teil ⇒
  sup_{D̄}|Π_N f_K − f_K| = O(N^{κ−4.2877}); Schwelle κ < 4.29 (Rest: κ ≈ a + 2 ≈ 2 bei σ_{r+1} ≈ 1).
(V3) ist eine Aussage nur über die r(N) Eckmoden (endlicher Rang, Grenzschicht-Theorie W), nicht über die volle Pol-Antwort — deutlich
schwächer als (E_A). Der K*-Weg unten bleibt als Alternative stehen.

**Variante X1 (X 00:58, von P geprüft 01:00 — Reduktion korrekt):** statt SVD/Q_N die Kleinste-Quadrate-Projektion Π̂ = R⁺T_{N+d}J_∞ (jedes
d ≤ 2N) und der kanonische Eckteil c_N(k_ζ) = Π_N k_ζ − Π̂k_ζ (Rang ≤ d; Formel c_N = J_N⁻¹Zᵀ(I+ZZᵀ)⁻¹·o(g), Z = BJ_N⁻¹, geprüft). Annahmen:
**(a″)** σ_min(T_{N+d}J_∞|_{P_N}) ≥ c₀ > 0 [O; numerisch 0.7319 = 1/1.3663] und **(Ẽ)** |c_N(k_ζ)|₁ ≤ CN^κe^{−εNg'(ζ)} auf Schlitz ∪ Γ'_δ [W].
Eckfreier Teil wie P-II in K' (Bernstein–Walsh, |ŷ_η|₁ ≤ CN² auf ∂K'), Rate N^{2−5.075} + N^κ(εN)^{−5.075}. Kein Q_N, kein r(N), keine Bulk-Annahme.
~~Achtung: (Ẽ) in K' grenzwertig~~ — überholt (Runde 3): mit p ≈ 1.5 ist die K'-Ecke nicht grenzwertig (p17b), X 01:08 misst Θ ≈ 0.13 für den
kanonischen Eckteil (N = 200…600, kleine d) bzw. Θ^c < 0 für d = 2N (N ≤ 300, X 01:22 = Satz X3). Für d = 2N zeigt P (PRUEFLISTE, X 01:22) den Schwanz
‖T_{>3N}J_∞(f_K)δ‖ ≤ CN^{−3.95}‖δ‖ (Darboux an L, L̄; unter der Bedingung (T3): max_{|x|=R}|f_K| < (R/ρ)³ für ein R < |L ± 2πi|).

**Alternative (K*-Weg, P 00:42) — Zielaussagen:**
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

**Satz (Hauptsatz im K*-Rahmen; Kette bewiesen, 9.4):** (K) + Univalenz (num.) + b35b′ + (a′)_r + (E_L⁺) + (E_sing) + (E_A) + Kantenschwanz ⇒ (H₀*) ⇒
sup_{D̄}|Π_N f_K − f_K| = O(N^{κ−4.81}); mit (Fam) und (A) ⇒ p_N → f_K. Ohne Zerlegung genügt die messbare einteilige Form:
R_N(L) ≤ CN^κ und Θ*_N := sup_{S_L ∪ Γ*_d} ω*_ζ(λ_N)/g*(ζ) ≤ 1 − ε. Numerisch: Θ_N(W̄) ≤ 0.105 (C, N ≤ 1000); Θ'_N(K') mit echten
R_N auf ganz ∂K' (p14/p15): −0.21/−0.13/−0.09/−0.06/−0.03/+0.03 (N = 100/200/300/400/600/800), konservativ (λ := 0 auf f_K(Γout), Cs
D-Bogen-Daten) ≤ 0.110 bis N = 1000 (p15b); K* liegt dazwischen.

**Gleichwertiger K'-Schluss.** Statt (E_A) genügt im K'-Rahmen, dass der Eckteil auf f_K(Γout) nahe L ABKLINGT:
‖Q_N y_η‖ ≤ CN^κ e^{−bN|η−L|^{q'}} mit q' < 1/α' = 1.2704 (dann ist der Ecken-Grenzwert in P-III −∞ trotz grenzwertigem Γin-Profil).
Die echte R_N auf f_K(Γout) klingt so ab (λ ≈ −g_W(η), weil f_K(Γout) ausserhalb W̄ liegt; p15: Schlitzverhältnisse ≈ −18); Os Eckanteile
auf f_K(Γout) fallen (o15). (a′)_r allein liefert dort nur „polynomial“, nicht „abklingend“. Beide Schlüsse sind dieselbe Auslöschung auf
der A-Seite der Ecke L.

### 9.4 Bewiesene Reduktionen (P, 11.10. 00:42; Prüfskripte p13–p17)

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
Eckteils auf ∂K* ∖ A_r höchstens max(Γin-Profil, 0) + O(log N/N), auf ∂K* ∩ A_r nach (E_A) ≤ 0; Satz P-III mit α = 0.883 und p ≈ 1.5 (W §14; früher 1.27) ⇒ (E_K^{∂*})
mit ε > 0, sobald der Fernteil λ_max·sup ω/g < 1 (W̄-Konstante 14.65 ⇒ λ_max < 0.068). Aussenwerte fern von L sind durch Os Poisson–Jensen O6(a)
automatisch ≤ max des Kreisprofils; nur die Aussenwerte NAHE L brauchen die Zusatzannahme in (E_L⁺). ∎ (Zusatz aus O6, korrigiert O §16.1:
Kantenschwanz ‖Q_N Y_N P_{cN}‖ ≤ CN^b für ein c ≥ 0.7 mit ln(R/ρ) > c₁/c; gemessen ≈ N^{0.5} bei c ≤ 0.7.)

**Offen und präzise benannt:** (a′)_r [O] · (E_L⁺) inkl. Aussenwerte nahe L [W] · (E_sing) [W/O] · (E_A) [neu: Auslöschung auf dem A-Bogen nahe L;
Kandidat Os Dualbild O2b — Kontrolle von I − Y_N auf H²(D(a₀,R₀))] · Kantenschwanz ‖Q_NY_NP_{cN}‖ ≤ CN^b [W, nach O §16] · **(Fam) = OFFENE AUFGABE F (9.6)** · (A) [P: Zweig rigoros bis N = 193; Kneser-Werte zertifiziert (§9.8); offen nur die Kopplungsrechnung und (H) ab N₀].

### 9.5 Phase-2-Ergebnisse der Stränge (Status nach P-Prüfung)

| Eintrag | Aussage | Status |
|---|---|---|
| X 23:53 | Halbmond-„Divergenz“ = Rundungsrauschen; Ersatzform-Identität Π_N f = T^{a₀}_N f − T^{a₀}_N[f'(e_N∘h)] | verifiziert |
| X §1.5, §2.4 | Möbius = Toeplitz × GL(2)-Darstellung; Selektion nur am Rand (exakte f_c ≠ f_K um Ū) | verifiziert (§2.4 |(K)) |
| M 23:55 | radiales Quadrat, einzige Familie | verifiziert |
| O 00:14 | O1 (exakte Faktorisierung, Mischformel Π_N g = Y_N T⁽⁰⁾g + (I−Y_N)T^{(a₀)}g) | verifiziert |
| O 00:14 | O2 (Herglotz-Kriterium ⇒ ‖(I+K_N)⁻¹‖ ≤ 1); O2b-Zahlen (arb) nicht nachgerechnet | verifiziert (Struktur) |
| O 00:14 | O3 (keine mellin-diagonale akkretive Norm bei |μ| ≠ 1) | verifiziert im lokalen Modell |
| O 00:14 | O5 (Möbius-exp unbedingt: Rate e^{−g_{K'}(ζ_p)} = 0.535 über K') — erste vollständige K'-Kette | verifiziert|O2b |
| O 00:32 | (a′)_r-Formulierung (σ_{r(N)+1} ≤ 1 + δ₀, δ₀ = 2.2e-5 Kneser); „a < 3.15 genügt“ (Buchhaltung grob, Richtung stimmt) | übernommen in 9.3 |
| O 00:32 | Lemma O6 (gemeinsame Eckformel) | verifiziert für polynomiales f; **für A_f ist (b) verletzt** (f_K singulär bei L auf |x| = ρ) ⇒ Ergänzung (E_sing), 9.4 P-IV |
| X 00:58 | Satz X1 ((LP) ⇐ (a″) + (Ẽ), kanonischer Eckteil) und Formel §3.6 | **verifiziert** (Cauchy-Darstellung, Zerlegung y = ŷ + c_N, Rationalität, P-II mit N+d Zeilen, Bernstein–Walsh; Kleinste-Quadrate-Formel per Normalgleichung + Push-through) |
| X 00:38/00:45 | Lemma 3.1 (Rechteck-Abschnitt: J_N⁻¹ = R⁺_top + R⁺_bot S⁻¹BR⁺_top, σ_{d+1}(J_N⁻¹) ≤ 1/σ_min(R)) und Korollar 3.2 (Π_N g = Π̂g + R⁺_bot S⁻¹ρ̂(g)) | **verifiziert** (Sylvester-Determinante, Linksinverse, Einsetzen nachgerechnet); Befunde zu Π̂ (keine Spitze) numerisch |
| M 00:43 (2) | adjungierte Gleichung in Zwei-Punkt-Form (J^Tv)_{N−l} = [t^l](mW − PP_{a₀}[W(h)h'(t/h)^{N+1}]) | **verifiziert** — unabhängig an f_K nachgerechnet (p19, N = 14: 1.3e-68) |
| M 00:43 (1) | Ort der Singularität entscheidet (Ecke: konvergent; D-Bogen: divergent) | numerisch (nicht nachgerechnet); passt zu D6 |
| O 00:48 | Lemma O7a (lokales Eckmodell T_L: E²(S) → H²(H), Multiplikator 1 + μ^{s−1}, ‖T_L⁻¹‖ = 1.2814) | **verifiziert** (Mellin–Plancherel im Sektor, Zweig von (με)^s, Gewichtsquotient) |
| O 00:48 | Satz O7 (T: E²(Ω) → H²(D) Fredholm, Ω = int K'; Index 0 per Homotopie) | verifiziert bis auf das zitierte Lokalprinzip (iv) und die numerischen Pfadbedingungen (Gitter, ohne Einschluss). **Kneser-Pfad f_t = (1−t)f_M + t f_K von P geprüft (p20, Gitter 2400, t = 0…1):** min|f_t'| 0.53 → 0.3456, f_t' ohne Nullstelle (Windung 0), Selbstabstand ≥ 0.030, (S3) überall, min|m_t| ≥ 0.7319, |μ_t| 1 → 1.1724, min|1+μ_t^{s−1}| ≥ 0.999999 (f_K ausserhalb D̄ über exp∘f_K∘Log, kein Motor nötig); Os P1–P5 (o21) ebenso erfüllt: min Re f_t' = 0.3456, min Re m_t = 0.7159, min X_t = 4.447 (p20b; Gitter, nicht zertifiziert) |
| O 01:00 | O7 Fassung 2 (Kleben O7d, Vertauschung O7e, lokale Inverse O7f über Koenigs-/Kneser-Koordinate e^{L·slog}) | O7d verifiziert (Links-/Rechts-Regularisator); O7e/O7f verifiziert bis auf die klassischen Zutaten K1–K3 (Coifman–McIntosh–Meyer, Kommutator-Kompaktheit, Montel) und die Normstetigkeit der Familie (Skizze); Kneser-Koordinate κ = Ψ(1 + O(Ψ^γ)) ist auf Ω nahe L eindeutig (Sektor 218° < 360°, Schlitz ausserhalb Ω) und C^{1,1} ✓ |
| O 01:00 | Satz O8: c_∞ = inf‖J_∞δ‖/‖δ‖ ≥ min|m|/‖G‖; c_∞ > 0 + ind 0 ⇒ ‖G‖ ≤ max|m|/c_∞ | verifiziert (J_∞ = M_m(I+K), ‖m g‖ ≥ min|m|‖g‖ auf H²) |
| W 00:24 | Γin-Profil-Exponent 1.26–1.29 (Ds 1.44–1.52 enthielt N^κ); Γout λ = 0⁻ | übernommen; macht K' grenzwertig (9.3) |
| W 00:40 | P's p10-Wert „+0.031“ auf dem A-Bogen war float-Rauschen; genau: Urbildterme klein | **richtig** (bei N = 400 stimmt meine float-Zahl −0.0152 mit Ws genauer −0.0153 überein; bei N = 1000 versagt float) — P korrigiert |
| W 01:06/01:29 (§14) | Exponent 1.27 zurückgezogen, Profil p = 1.5–1.6; (E_L)-Weg trägt in W̄ (Θ_∞ 0.22/0.10; Kneser 0.215/0.103), (E_A) gestützt | **übernommen; Folgerung nachgerechnet** (p17b: K'-Ecke 0.19/0.155, nicht grenzwertig — §9.3 korrigiert) |
| X 01:08 | kanonischer Eckteil in K': Θ ≈ 0.13, Ecken-Grenzwert −30 | plausibel, konsistent mit p17b (nicht nachgerechnet) |
| X 01:22 | Satz X3 (d = 2N: (a″)_∞ + einteilige Randbedingung Θ^c), Kneser Θ^c < 0 (N ≤ 300) | Reduktion wie X1 (verifiziert); Schwanz ‖T_{>3N}J_∞δ‖ ≤ CN^{−3.95}‖δ‖ von P skizziert (Darboux, Bedingung (T3)) |
| O 01:11 | O7g (Kommutator), O7h (Spuren, Carleson), O7i (Normstetigkeit) ⇒ O7 bis auf Standardsätze; O7j (T invertierbar für Kneser, qualitativ); o22: σ_min(T_{3N}J_∞|P_N) > 0.7315, N ≤ 100 | O7g/O7h/O7i **verifiziert** bis auf CMM/David, Calderón, Pommerenke, Carleson, Kato (Standard); O7j **verifiziert|(G2),(S),P1–P4,D4** (D4 nur Skizze); o22 **bedingt** (Kneser-Koeffizienten nicht zertifiziert, §9.8) |
| M 01:04/§8 | keine Gårding-Ungleichung; (a″)_d mit d ≈ r(N) ≡ (a′)_r; Lemma M3 (σ_min ≤ min_Γin|m| + C/N für alle d — 0.690525 scharf); Lemma M2 | M3 **verifiziert bis auf F = 0.8751 < 1 (Gitter)**; M2 **verifiziert** |
| O 02:09 | Satz O9 (c_∞ ≥ c(λ) aus ‖Tδ‖_p ≥ λ‖Bδ‖_p; c(0.98988) = 0.568, c(0.5) = 0.340) und Lemma O10 (Ecke: (FA) bis λ_L) | O9 **verifiziert**; O10 **verifiziert** (Formel), λ_L = 0.98988321 (P, mpmath; O-Text 0.9898822 vermutlich Tippfehler) |
| P Runde 3 | **Kohärentes Zweig-Zertifikat** (p24): jeder Übergang N → N+1, N = 45 … 192, mit EINEM s-Segment | **rigoros** (arb), §9.7 (3b) |
| P Runde 3/4 | Kneser-Werte mit Schranken: in gp-tetration NICHT vorhanden; **Satz P-V bewiesen** (sexp auf 1.3e-168, f_K-Koeffizienten a_0…a_640 mit Fehler·ρⁿ ≤ 2.9e-127) | **rigoros** (arb; zitiert: Ahlfors–Bers, Stoilow, Beurling-Norm, Paulsen–Cowgill), §9.8 |

### 9.6 Abhängigkeiten (Stand 00:42; Runde 3 ergänzt)

(K) ⇒ (F) [L8.2–8.4, b35b′ num.] ⇒ (F_{K*}) [+ Univalenz num.] · (F_{K*}) + (H₀*) ⇒ (LP), Rate N^{κ−4.81} [P-I] · (F_{K*}) + (H*) ⇒ H1* ⇒ (mit (A)) H2* ⇒
Vermutung · (a′)_r ⇒ (E_R*) auf ∂K*∖A_r [P-II] · (E_A) ⇒ beide Teile auf ∂K*∩A_r · (E_L⁺) + (E_sing) + Kantenschwanz ⇒ (E_K^{∂*}) [P-IV + P-III] ⇒ (E_K*) [C 13.2] ·
(E_R*) + (E_K*) ⇒ (H₀*) [P-I] · BEVORZUGT: (a′)_r ⇒ Rest in K' [P-II] + (V1)–(V3) + (E_sing) + Kantenschwanz ⇒ Eckteil in W̄ [P-IV, C 13.2] ⇒ (H₀^mix) ⇒ (LP) mit N^{κ−4.29} [P-I] ·
(H₀*) bzw. (H₀^mix) ⇒ (H) [Familien-Stabilität, offen] · alternativ einteilig: R_N(L) polynomial + Θ*_N ≤ 1−ε ⇒ (H₀*) [P-I(ii)].

**OFFENE AUFGABE F — Familien-Stabilität (H₀) ⇒ (H)** (sichtbar geführt; die Leitung vergibt sie).
*Gegeben:* (H₀) in einer der Formen 9.3 (K', gemischt, X1/X3) für die feste Linearisierung J_∞ = J_∞(f_K).
*Zu zeigen:* (H-i)–(H-iii) aus §5.1 für ALLE Fast-Lösungen P ∈ 𝒮_N ∩ B_N, also für die Sekantenoperatoren 𝒥_P v = v∘P + M_P v mit
P − f = λ(Π^P_N f − f) + Kopfterm, sup_{D̄}|P − f| ≤ R_N = O(N^{κ−4.69}): dieselben Pol-Antwort-Schranken (bis auf einen Faktor) und die Kopf-Stabilität.
*Die Schwierigkeit:* v∘P − v∘f braucht v auf P(D̄) ≈ f_K(D̄) ⊄ D̄ — für Polynome vom Grad N nur mit e^{N g}-Verlust kontrollierbar; tragfähig sind nur
**kohärente** Störungen (Π^P f − f klein auf ganz K', nicht nur auf D̄). Eine Kugel aller Polynome um f in sup_D̄ ist zu gross (p22: die Kopfantworten
J_N⁻¹x^k sind auf P(D̄) exponentiell gross, ~1.3^N).
*Stützen:* p09 (λ-Pfad N = 60: Pol-Antworten unverändert, (H-iii)-Quotient 1.0000); p24 (längs der echten Übergänge ändert sich J_N relativ nur um
θ ≤ 6e-8, N ≤ 193). *Bausteine:* Störungsformel E^P_N − E_N über 𝒥_P − J_∞ = (Komposition) + (M_P − M_f); Os c_∞ > 0 ist offen in P (Kato).
*Entlastung durch p24:* (Fam) wird nur noch für N ≥ N₀ gebraucht, N₀ ≤ 193 frei wählbar.

### 9.7 Startbedingung (A) als endlicher Test — Konstanten und Zweig-Zertifikat (P, Stand 00:54)

**(1) Explizite Konstante in 4.1 / H1** (W̄- bzw. K*-Ecke; Zahlen p01, A: C_J ≤ 0.95 für u ≤ 0.2, dist(L(1+u), K) ≥ 0.933ρu, c_g = 0.829):
  sup_{D̄}|Π^P_N f − f| ≤ C·N^κ·[ 12.7·ε^{−4.81}·N^{−4.81}  +  (M_d|Γ_d|/(2πδ_d))·e^{−εNd} ],
Schlitzteil 2·(0.95/(0.933π))·αΓ(5.447α)·c_g^{−5.447α} = 12.7 (exakt nachgerechnet), Konturteil mit M_d = max_{Γ_d}|f_K|, δ_d = dist(Γ_d, K).
Bei d = 0.1 (M_d ≈ 3.4, |Γ_d| ≈ 9.6, δ_d ≈ 0.1) dominiert der Konturteil bis N ≈ 350; mit grösserem d (f_K ist holomorph weit über K hinaus,
nächste Singularitäten bei L + 2πik) ist er schon ab N ≈ 30–40 vernachlässigbar. Dann C' ≈ 2·12.7·ε^{−4.81}·C ≈ 40·C (ε ≈ 0.9) und
**R_N = K₀C'N^{κ−4.81} ≈ 240·C·N^{κ−4.81}** (K₀ = 6). Gemessen ist sup|p_N − f_K| ≈ 1.3·N^{−4.85}: Reserve ≈ 200 bei C ≈ 1, κ ≈ 0.

**(2) Was (A) bei N₀ verlangt:** (i) N₀ ≥ N₀(H) (ab dort gelten (H-i)–(H-iii) mit den Konstanten aus (1)); (ii) ein rigoroser Einschluss von p_{N₀}
**als Zweig** (Übergänge 1 → 2 → … → N₀ entlang der σ-Homotopie zertifiziert); (iii) zertifizierte Kneser-Werte auf D̄ mit Fehler ≪ R_{N₀}
(bei N₀ = 40: R ≈ 5e-6·C, also f_K auf D̄ auf ≈ 1e-6 genau — die Motorwerte sind dafür um 25 Stellen zu genau, aber NICHT zertifiziert).
Dann ist p_{N₀} ∈ B_{N₀} eine endliche Intervallrechnung; ab N₀ übernimmt H2 die Übergänge.

**(3) Zweig-Zertifikat (p18, rigoros, arb):** Krawczyk mit Parameter s ∈ [s_a, s_b] und c-Intervall auf Boxen, die pad(p_N) und den Pfad enthalten,
mit adaptiver Unterteilung von [0,1] — zertifiziert die σ-Homotopie N → N+1 lückenlos für **N = 1 … 45** (Zweig p₁ → p₄₅ rigoros, Box-Radien ·ρ^k ≤ 6e-37; Daten p18_zertifikat.json). Nötige s-Teilschritte je Übergang N = 37…44: 14, 1, 33, 8, 75, 23, 122, 60 (abwechselnd, Wachstum ≈ 1.6–2 je Grad). Grenze der Methode:
Koordinaten-Boxen müssen den ganzen Schritt Δ = −c J⁻¹x^{N+1} enthalten; die Intervall-Jacobi-Matrix über so eine Box hat Breite ~ N·1.62^N·|Δ|
(Komposition verliert den Radius, 2.24/ρ = 1.62, vgl. C/D: ω ~ 1.25^N), also braucht man ~1.62^N·|Δ| Teilschritte — exponentiell (bei N = 40: 33).
Für grössere N₀ braucht der Übergang eine kohärente Verifikation (Taylor-Modell in s längs Δ mit Punkt-Matrizen statt Koordinaten-Boxen) —
oder man beweist (H) mit N₀(H) ≤ 45. Nach (1) ist N₀(H) ≈ 30–40 mit grossem d nicht unplausibel.
**(3b) Kohärentes Zweig-Zertifikat (p24, Runde 3, rigoros, arb) — das 1.62^N-Wachstum ist beseitigt.** Statt einer Koordinaten-Box um den GANZEN
Pfad eine Röhre Z(s) = z̃(s) + B: z̃ = Interpolationspolynom (Grad 6) durch Newton-Punkte des Pfades an 7 Tschebyschow-Knoten in [0,1], B eine winzige
Box (Radien r·ρ^{−k}, r ≈ 1e-55…1e-60). Krawczyk für alle s ∈ [0,1] und c ∈ C zugleich:
K − z̃(s) = −Y G_{s,c}(z̃(s)) + (I − Y J(z̃(s)))B − Y(J(Z) − J(z̃(s)))B ⊂ int B. Dabei werden s ↦ G(z̃(s)) und s ↦ J(z̃(s)) (Polynome in s) NUR an den
Knoten in Punktarithmetik ausgewertet (kein Wrapping) und über die Tschebyschow-Restformel (|ω| ≤ 2^{1−2n}, Hermite–Genocchi) auf [0,1] fortgesetzt;
die n-te Ableitung wird über Majoranten-Reihen auf |s| ≤ R (R ≈ 1/(M‖â₁‖) ~ 1e8…1e11) mit Cauchy beschränkt. Das Wrapping 1.62^M steht nur noch
vor R^{−7} bzw. vor r (Präzision 2M + 320 Bit genügt; Knoten auf 2^{−prec} gerundet — |ω| ≤ 2^{1−2n} und Λ ≤ (2/π)ln n + 1 gelten
dann bis auf relative 1e-100, gedeckt durch die Reserve 0.5; Kopfkoeffizient c als Intervall C über die Box von p_N). Ergebnis: **alle Übergänge N → N+1 für N = 45 … 192 mit EINEM s-Segment**, Reserve 0.5,
linearer Term θ = sup_s‖I − YJ(z̃(s))‖_w ≤ 6e-8 (fällt mit N), Wrapping-Term ≤ 3e-8 (relativ; wächst wie ~1.62^N·r — für N ≳ 230 braucht es
kleineres r, d. h. mehr Präzision und Knoten); Rechenzeit 1 s (N = 46) … ≈ 45 s (N ≈ 190) je Übergang
(p18 brauchte bei N = 43/44 schon 122/60 Segmente). Mit p18 zusammen: **Zweig p₁ → p_193 rigoros als σ-Homotopie-Zweig** (Box-Radien ·ρ^k ≤ 1e-120;
Daten p24_zustand.json, Protokolle p24_lauf*.log). Weiter ist nur Rechenzeit (~M³).
Vorstudien: p22 (eine Kugel in sup über D̄ ∪ P(D̄) trägt NICHT: max_k sup_{P(D̄)}|J⁻¹x^k|/ρ^k = 88 → 1.6e4 für N = 20 → 40, ~1.3^N), p23 (die Richtungs-
Ableitung A = J⁻¹∂_tJ(P₀ + tΔ₁) ist winzig: ‖A‖_w = 6e-6 → 4e-9 für N = 20 → 60 — das Wachstum in p18 war reines Intervall-Wrapping).
**Fazit:** (A) ist kein unabhängiger Test, sondern: [(H) mit expliziten Konstanten ab N₀] + [Zertifikat bis N₀ — **erledigt für jedes N₀ ≤ 193**] +
[zertifizierte f_K-Werte auf D̄ auf ≈ R_{N₀} (N₀ = 100: ≈ 6e-8·C; N₀ = 150: ≈ 8e-9·C) — **erledigt durch Satz P-V (§9.8), Genauigkeit 1e-168**].

### 9.8 Zertifizierte Kneser-Werte — Satz P-V BEWIESEN (rigoros, arb; P, Stand 11.10. 03:57)

**Ergebnis.** Knesers sexp zur Basis e ist auf dem Streifen −½ ≤ Re z ≤ ½ (alle Höhen) mit Fehler **ε ≤ 1.3e-168** durch explizite Daten
(p26_iter6.json) eingeschlossen: sexp_K = T̃ ± ε für |Im z| ≤ 0.45 (T̃ reelles Polynom vom Grad 519 um 0) und sexp_K = Ũ ± ε für Im z ≥ 0.35
(Ũ(z) = P(L^{z+θ̃(z)}), P die Poincaré-Funktion an L, θ̃ trigonometrisches Polynom mit 241 Moden), unten gespiegelt, sonst über exp/Log.
Daraus **zertifiziert: Knesers f_K-Koeffizienten a_0 … a_640 mit |Fehler|·ρⁿ ≤ 2.9e-127** (≥ 126 gewichtete Stellen; a_0 auf 8.8e-143,
a_600 auf 3.8e-211 absolut) und |a_n| ≤ 2.128·1.30^{−n} für alle n — Datei **gemeinsam/kneser_halbexp_koeff_zert_640.json** (Os Anfrage).
Gegenprobe: Strang As 50-stellige Datei stimmt damit auf 1.1e-43 (gewichtet) überein.

**Beweis (Fassung 2 — quasikonforme Glättung; das Cousin-Netz der Fassung 1 hatte am Tripelpunkt einen log-Verlust und ist verworfen):**
1. *Daten (Numerik, p26):* Cousin-Newton auf dem Zylinder (Schleife Im z = 0.25, Naht Re z = −½, Knotenkonstante per Rampe), Start aus
   dem Motor (dps 50; komplexe Werte des Forks versagen ab dps 60), sechs Schritte, Residuen 2e-31 → 2e-59 → 4e-95 → 5e-130 → 5e-169.
2. *P zertifiziert* (kc_rig, p25): p_n per Rekursion, |P − T_{n0}| ≤ C|w|^{n0+1} auf |w| ≤ ½ (Fixpunktargument); Ũ ist ganz mit Ũ(z+1) = e^{Ũ(z)} EXAKT.
3. *Defekte (p27):* e mit Ũ(ζ+e) = T̃(ζ) auf G_B = [−0.55, 0.55]×[0.30, 0.50], d mit T̃(ζ+1+d) = e^{T̃(ζ)} auf G_S = [−0.65, −0.45]×[−0.55, 0.55]:
   Maximumprinzip + Taylor-Modelle (Ordnung 240, Cauchy-Schwanz, 900 Bit) ⇒ sup|Ũ − T̃| ≤ 2.7e-173, sup|T̃(·+1) − e^{T̃}| ≤ 2.1e-173; |Ũ'| ≥ 0.83,
   |T̃'| ≥ 1.14 (Taylor-Modelle auf Kästchen) ⇒ **|e| ≤ 6.5e-173, |d| ≤ 3.7e-173** (skalares Kantorovich, eindeutig im Radius 0.02).
4. *Exakte Lipschitz-Lösung F₀ (p28):* W₀ = χ(Re z)·d(z−1) (χ: 0 → 1 auf [0.4, 0.467]), Y = W₀ + e(·+W₀), V₀ = ψ(Im z)·Y (ψ: 1 → 0 auf
   [0.40, 0.45]); F₀ = T̃∘(id+W₀) für |Im z| ≤ 0.40, F₀ = Ũ∘(id+V₀) für Im z ≥ 0.35 (dort gleich). **F₀(z+1) = e^{F₀(z)} gilt exakt**: an der Naht
   nach Definition von d, oben weil Y periodisch ist — wegen der exakten Identität e(ζ) − e(ζ+1+d(ζ)) = d(ζ) (folgt aus der Exaktheit von Ũ und
   der Eindeutigkeit von e). Damit gibt es KEINE Tripelpunkt-Verträglichkeitsbedingung mehr (die Fassung-1-Rampe wird nicht gebraucht).
5. *Beltrami (p28):* |μ₀| = |∂̄F₀/∂F₀| ≤ 1.3e-171, Träger in |Im z| ≤ 0.45. Messbarer Riemannscher Abbildungssatz auf dem Zylinder (q = e^{2πiz},
   Hauptlösung ψ₀ = q + C[h], h = (I − μS)⁻¹μ in L⁴, ‖S‖₄ ≤ 12) ⇒ qk φ mit φ(z+1) = φ(z)+1, **|φ − id| ≤ 1.3e-169**; Stoilow: F = F₀∘φ⁻¹ ist
   analytisch, F(z+1) = e^{F(z)}, reell auf ℝ, F → L (Im → +∞), L̄ (Im → −∞).
6. *Identifikation (p28):* Prüfung A (F ∉ (−∞, 0] auf dem Streifen) und B (F ∉ (0, 1] für 0 < Re ζ < ½; nahe 0 analytisch über F(0)=1,
   F'(0) ≥ 1.0918, |F''| ≤ 1.71) auf 2958 Scheiben für 0 ≤ Im ≤ 3, Schwanz Im ≥ 3: |F − L| ≤ 0.065 < Im L ⇒ F holomorph auf Re z > −2 mit
   beiden Grenzwerten für jedes x > −2 ⇒ **Paulsen–Cowgill Prop. 2** (Fassung Paper VI Thm 11.11): F(· + τ) = sexp_K, |τ| ≤ 5.2e-169.
7. *Fehler (p29):* ε_T, ε_U ≤ 1.27e-168 (mit sup|T̃'| ≤ 1.95, sup|Ũ'| ≤ 1.69).
*Zitierte Sätze:* Ahlfors–Bers (messbarer RMT, Hauptlösung), Stoilow-Faktorisierung, Beurling-Transformation ‖S‖_p ≤ 4(p*−1)
(Bañuelos–Wang 1995; bei |μ| ~ 1e-171 genügt jede endliche Schranke), Paulsen–Cowgill, Adv. Comput. Math. 43 (2017), Prop. 2.
**f_K-Koeffizienten (p29):** f_K(x) = sexp(slog x + ½) auf |x| = 1.2 (4096 Punkte; Krawczyk 1D für slog auf D(w̃, r₀) mit sexp-Bällen, dann
sexp-Ball auf D + ½; Fehler je Wert ≲ 3e-167), DFT, Aliasfehler ≤ M r'^{−n}(1.2/1.3)^{4096} mit M = sup_{|x|=1.30}|f_K| ≤ 2.128 (adaptive
Bogen-Krawczyk, 1950 Scheiben).
**gp-tetration (nur gelesen):** dort gab es für b = e nur Vorzertifikate (Paper VI setzt den Anker T_e voraus); Satz P-V liefert genau
dieses Anker-Zertifikat (Def. 11.3) — mit 168 Stellen.
**Als Nächstes (P):** (i) f_K und f_K' auf ∂D inkl. Umgebung von L (fK_zert vorhanden) ⇒ Os Pfadbedingungen P1–P4 streng; (ii) Startbedingung (A):
Π_{N₀}f_K aus den zertifizierten Koeffizienten (T_N J_∞f_K mit Schwanzschranke) und sup_{D̄}|p_{N₀} − Π_{N₀}f_K| ≤ ‖·‖_{ℓ¹_ρ} gegen die Box aus p24;
(iii) o22 mit diesen Koeffizienten (Sache von O).

**Prüfstellen für den Gegenleser von Satz P-V (zuerst ansehen):**
1. p28, Schritt 4: exakte Periodizität F₀(z+1) = e^{F₀(z)} an der Naht und oben; die Identität e(ζ) − e(ζ+1+d(ζ)) = d(ζ).
2. p28, Schritt 5: Beltrami-Schranke, Neumann-Reihe in L⁴ (‖S‖₄ ≤ 12), Normierung φ(z+1) = φ(z) + 1, Schranke |φ − id| ≤ 1.3e-169.
3. p28, Schritt 6: Überdeckung der Prüfungen A/B (2958 Scheiben, Schwanz Im ≥ 3) und ob die Voraussetzungen von Paulsen–Cowgill Prop. 2 wörtlich erfüllt sind.
4. p27: Reste der Taylor-Modelle (Cauchy-Schwanz aus grober Kreisschranke) und das Maximumprinzip für die Defekte e, d.
5. kc_rig._restkonstante: Kontraktionsargument für |P − T_{n0}| ≤ C|w|^{n0+1}.
6. p29 / kc_kneser: Fortpflanzung von ε (sup|T̃'|, sup|Ũ'|, τ); der Ableitungszuschlag „ε·100“ setzt voraus, dass die ε-Schranke auf einer 0.01-Umgebung gilt.
7. Log-Zweige in sexp_ball (Hauptzweig für Re z < −½, gestützt auf Prüfung A).

### 9.9 Stand beim Anhalten (11.10. 04:50; P) — Hypothese, Zertifikate, Prüfliste, Wiederaufnahme

**Hauptsatz, Stand.** p_N = p*_N → f_K ist bewiesen unter **(K)**, **(A)** und **einer benannten Hypothese (H₀)**. (H₀) ist die lokale Pol-Antwort
an P = f_K, also Cs (E) = (E_R) ∧ (E_K) (§4.1). Dazu kommt die Familien-Erweiterung (H₀) → (H); sie ist offen.
Ws Reduktion (W-BERICHT §22–§26) führt den Eckteil (E_K) auf eine Aussage auf **einer** Kurve zurück:
- Satz W1'' (bewiesen): (H*) |Y_i| ≤ CN^κ auf γ* genügt.
- Skizze (§23.2): (E_K) ⇐ (I) Wachstumsschranke + (H*)_loc in der Eckzone + Kreisbedingung (KB) auf |t − L| = A·N^{−0.894}.
- Satz W5 (bewiesen): (A1) + (A2) ergeben (H*) nur bis auf e^{o(N)}. Das reicht an der Ecke nicht (§26.3), dort ist die polynomiale Form nötig.

Keine dieser Eingaben ist für alle N bewiesen; numerisch gelten sie bis N = 1000. (K) ist für b = e in den gebrauchten Teilen durch Satz P-V
zertifiziert (§9.8).

**Neu bewiesen (P, p32a, rigoros):** Knesers exakte obere Struktur. Für Im z ≥ 0.36 gilt sexp_K(z) = P(e^{L(z+θ_K(z))}) mit θ_K
1-periodisch und holomorph und |θ_K − θ̃| ≤ **4.4e-167**. Das ergibt |κ_n − c̃_n| ≤ 4.4e-167·e^{2πn·0.36}.
Beweis: Überdeckung des Streifens [0.36, 3] mit 1742 Quadraten (Taylor-Modell von Φ = P(e^{L·}), Krawczyk mit Radius 1e-120); für Im ≥ 3
Rouché + hebbare Singularität von W·e^{−Lz} in q = 0; dann Identitätssatz.
Folgerung nahe L: f_K(x) = P(√L·w·R(q)) mit w = P⁻¹(x), q = e^{2πi slog x} und R(0) = 1. Daraus **(K3): f_K(L) = L und f_K'(L) = √L**
als Randwerte in D̄, durch Grenzwertbildung (|q| ≤ C|w|^{4.447}).
Probe in p32: ∂_θf an θ = 1.3372 ist −1.46149 − 0.67915i = i·L·√L.

**(A), Stand (p31, rigoros, Boxen p30 = Kette p18 + p24 auf aicgolling reproduziert).** Für N = 140…193:

| N | ‖p_N − Π_N f_K‖_{ℓ¹_ρ} | ‖p_N − T₆₄₀f_K‖_{ℓ¹_ρ} |
|---|---|---|
| 140 | 2.2e-21 | 1.52e-10 |
| 150 | 1.3e-21 | 1.10e-10 |
| 193 | 1.3e-22 | 3.1e-11 |

Alle Werte in p31_A.json. Die Kopplung braucht wegen f'∘f = eˣ/f' nur a₀ … a_{N+1}.
Damit ist sup_{D̄}|p_N − f_K| ≤ (obige Zahl) + S₆₄₀ mit S₆₄₀ = Σ_{k>640}|b_k| (numerisch 4.4e-15). **Rigoros fehlt nur S₆₄₀, also (K2).**
(A) verlangt sup ≤ R_{N₀}; bei N₀ ≈ 150 ist R ≈ 1e-3 (X). Die f_K-Werte müssen also nur auf ≈ 1e-4 bis 1e-5 genau sein; vorhanden ist 1e-127,
die Reserve beträgt ≈ 1e7.
Mit beweisbaren Konstanten ist N₀ ≈ 390 (X 04:10). Die Kette bis 400 übernimmt X; Start aus p30_zustand.json (N = 193).

**(K2) für M und O, Stand: OFFEN.**
- p = 4 reicht nicht für Ms Restsummen. Dafür wäre C₄ ≤ 8e-4 nötig; numerisch ist C₄ ≈ 8.4.
- Nötig ist p = 6: |b_k| ≤ C₆k^{−6} mit C₆ = ‖∂_θ⁶f_K‖_{L¹}/2π, numerisch ≈ 124 (As 4000 Koeff., Fejér).
- Mit C₆ ≤ 500 sind Ms Grenzen (Σ_{k>400}|b_k| ≤ 1e-12, Σ_{k>400} k|b_k| ≤ 1e-9) erfüllt.
- p32 zertifiziert ∫|∂_θ⁶f| zellenweise (Mittelwertform mit ∂_θ⁷f). Rigoros fertig sind nur die Stücke [0, 0.3645] (≤ 30.13) und
  [2.811, π] (≤ 4.54). Beide liegen im T̃-Gebiet.
- Im Strukturgebiet (Im slog ≥ 0.39) überschätzt die Kugelauswertung stark: G_k über einer q-Kugel, Faktor ≈ 4e7 je Einheit Radius.
  Abbruch an dieser Stelle.

**Prüfliste.**
- Erledigt in Meilenstein 6:
  - Satz W1'' korrekt (Exponent κ+1 genügt schon).
  - Lemma W3 korrekt (in der Konturdarstellung den Faktor 1/ζ nachrechnen).
  - Satz W5 korrekt (Saff–Totik I.6.8, Jensen). Die „Umkehrung (Blatt–Saff–Simkani)“ braucht deren Zusatzbedingung; sie ist eine Bemerkung, kein ⇔.
  - Maximumprinzip für den Innenkreis: korrekt und trivial; die Last liegt ganz in (KB).
- Offen:
  - O9/O10: Bestätigung (Runde 3 geprüft; λ_L = 0.98988322 nach M).
  - X4a/b, X5, X6, Lemma X7 (mit den Gitterminima c_in = 0.6741, c_out = 0.3456), X8, Kontur d = 0.3 mit der Sprungschranke bis u ≈ 0.46
    (braucht Kneser-Werte, kc_kneser genügt).
  - Text (H-i) „Tangente nur auf 𝒵“.
  - M: mk_kneser_arb, mk_ecke_kneser, m22–m26, M5, M7–M10, M13.
  - Ketten-Prüfung von Xs Verlängerung bis 400.

**Wiederaufnahme (nächster Schritt und Dateien).**
1. **(K2) fertigstellen:**
   - In p32_rand_ableitungen.py die Strukturzellen zentriert auswerten: G_k und H_k am Zellmittelpunkt (Punkt-Jets, scharf). Die Variation
     über die Zelle entweder per Taylor-Modell in z (Punkt-Jets bis Ordnung ≈ 30, Cauchy-Rest aus grober Kreisschranke) oder über
     ∂_θ⁷f mit derselben zentrierten Form.
   - Dann läuft `python p32_rand_ableitungen.py rechts|links 40` je ≈ 10 min auf aicgolling.
   - Danach `python p33_K_datei.py`; es schreibt gemeinsam/kneser_halbexp_b_zert.json im Format von M. Damit laufen Ms m26 und Os O11/o22.
   - Gebraucht: kc_basis.py, kc_rig.py, kc_kneser.py, p26_iter6.json, p27_zert.json, p28_zert.json, p32a_theta_K.json.
   - Server: ~/halbexp/p_pruefung, venv ~/halbexp-venv mit python-flint 0.9.
2. **(A) abschliessen:** Boxen bis N₀ aus Xs Kette, dann `python p31_A_kopplung.py N…` und S₆₄₀ aus (K2) addieren.
3. **Os P1–P4 (Kneser-Pfad):**
   - f_K und f_K' auf ∂D aus (K1) + (K2) als Koeffizientensumme mit Schwanz.
   - P2 braucht f_K' auf f_t(∂D) ausserhalb von D: über h_K = f_K⁻¹ oder slog_krawczyk aus kc_kneser.
4. Prüfliste (oben, offen).

Es laufen keine Prozesse mehr, weder auf dem Laptop noch auf aicgolling. Der Kneser-Motor aus gp-tetration lief in diesem Meilenstein nicht.
