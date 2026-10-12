# Offenes Problem 2: Polynom-Wurzel von exp → Kneser — bedingter Satz auf K_{r₀} und die verbleibende Kohärenzaussage

Stand 12.10.2026, Strang K (Phase 3), Zahlen von K und V. Kennzeichnung: **bewiesen** (V: von einem zweiten Strang geprüft) / **numerisch** / **Vermutung**.
Das ältere Problemblatt (w_grenzschicht/OFFENES_PROBLEM.md, Hypothese (H₀)) bleibt gültig. Dieses Blatt gibt den Stand nach Phase 3.

## 1. Rahmen
- L ≈ 0.31813 + 1.33724i ist der Fixpunkt e^L = L, ρ := |L|, D := {|x| < ρ}, μ := √L.
- f_K ist Knesers Halbiteration (f_K∘f_K = exp, singulär in L, L̄). p_N sind Polynome vom Grad N mit T_N[p_N∘p_N − exp] = 0 (T_N schneidet nach Grad N ab),
  verbunden durch eine σ-Homotopie ab p₁. Der Zweig bis N = 400 ist streng eingeschlossen; sup_{D̄}|p₄₀₀ − f_K| ≤ 1.75e-12 **(bewiesen, V)**.
- Linearisierung: **J δ := δ∘f + m·δ, m := f′∘f** (Komposition + Multiplikation); Abschnitte J_N := T_NJ|_{P_N}, skaliert J̃ = DJD⁻¹, D = diag(ρ^k).
- Testmodell: ein reelles Kubik f mit f(L) = L, f′(L) = μ. Lokal an L stimmt es mit Kneser auf 4e-4, die Zahlen unten auf 2–3 Stellen.
- Köpfe (Rang-1-Ränderung): c_n = J̃_n⁻¹e_n, ℓ_n = J̃_n^{−T}e_n, δ_n(z) = Σc_{n,k}(z/ρ)^k, Λ_n(t) = Σℓ_{n,k}(ρ/t)^k.
- Kern K_n(z,t) = ε_zᵀJ̃_n⁻¹ẽ_t = t·(J_n⁻¹T_n[1/(t − ·)])(z), mit ΔK_n = δ_n(z)Λ_n(t)/ℓ_nn **(bewiesen, V)**.
- Eck-Pol-Antwort u_n := J_n⁻¹T_n[1/(L − ·)], also ΔK_n(z,L) = L·Δu_n(z). Eckskala d_n = ρn^{−a}, a = 1 + iL/(4π).
- Geometrie: Γin = {|x| = ρ, |f(x)| < ρ} (Bogen über −ρ, der „Steg“ bei 180°), Γout der Rest. K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄)).
- Kontur 𝒞: Log-Spiralen Σ₁ (aussen) und Σ₂ (innen) in L (Steigung 1.623°/e-Faktor), eine Linse |t| = ρ(1 + S) über Γin, der Γout-Bogen, konjugiert bei L̄.
  E := Abschluss von int𝒞 ∪ f_K(int𝒞) ist eine Quasischeibe mit Spiralecke (Öffnung ≈ 132°) **(numerisch)**.

## 2. Bedingter Hauptsatz (K_{r₀}-Fassung)
**Satz A (bewiesen bis auf die Annahmen; Logik V).** Unter (EP), (FA), (S_E) gilt sup_{K_{r₀}}|Π_Nf_K − f_K| ≤ C·N^{γ − s_E}, γ = max(0, 1 + κ_G + κ_A).
Dabei ist Π_N die Galerkin-Projektion. Gemessen: κ_G ≈ −3.2, κ_A = 0, also γ = 0 und die Rate N^{−3.96}.
**Hauptsatz (bedingt).** Gelten diese Annahmen gleichmässig für die Familie der Newton-Nachbarn (F1–F4), dann ist p_N für N ≥ 400 die einzige Lösung in einer Kugel um f_K,
und sup_{K_{r₀}}|p_N − f_K| → 0 polynomial. Daraus folgen a_k^{(N)} → a_k für jedes k und p_N(1) → sexp_K(½).
**Nicht behauptet ist gleichmässige Konvergenz auf D̄.** In einem Modell mit schneller Spirale (gleiche Struktur) konvergiert Π_NU − U in einer festen Γout-Zone nahe L
ab N ≈ 1000 nicht mehr (Plateau bis N = 2800; V, Test T1). Grund: Der primale Kopf wächst dort in einem Sektor, der sich mit 1.45° je e-Faktor in n dreht.
Darum r₀ > r_* (Kneser: r_* ~ 1e-6, Schätzung).

## 3. Die offenen Aussagen (präzise)
- **(G-P)** G_n := sup_{z∈K_{r₀}}|Δu_n(z)| ≤ C·n^{κ_G}. Numerisch κ_G ≈ −3.2 bis n = 1800 (kubisch = Kneser auf 1 %).
- **(A-lok)** sup_{t∈Σ₁∪Σ₂, |t−L|≤r₁}|Λ_n(t)/Λ_n(L)| ≤ C·n^{κ_A}. Numerisch 1.008 bis n = 2000.
- **(FA)** lim sup_n max_{𝒞∖U_{r₁}(L,L̄)} (1/n)ln|Λ_n/Λ_n(L)| < 0. Numerisch: Linse (S = 0.03) −0.011 … −0.013, Γout-Bogen −0.037 … −0.053. Folgt aus einem Grenzproblem von W (W12′).
- **(S_E)** dist_E(f_K, P_N) ≤ C·N^{−5.417α_E}, α_E = Aussenwinkel/π. **Bewiesen** mit Andrievskii–Pritsker–Varga, Satz 2, bis auf die numerisch geprüfte Geometrie von E.
- **Familie (F3′/(BU))**, gezeigt mit einem Bootstrap-Gerüst (Logik bewiesen):
  - Entlang der Homotopie gelte CF(P) := min_g sup_γ|P − f_K + x^{N+2}g| ≤ δ₀ auf der Linsenkontur γ.
  - Dann gilt CF ≤ δ₀/2 für alle N, sobald **(BU)** gilt: ‖T_{N+2}(Σ_{n≤N}[w^m](δ_n∘k)·π_n/ℓ_nn)‖ ≤ C·N^{κ_B − s_E}.
  - Dabei ist k: 𝔻 → int γ konform und π_n = ⟨Λ_n, Daten⟩ sind die Kopfpaarungen.
  - Termweise (nur |π_n|) schliesst (BU) höchstens mit Marge 0, und das nur auf einer kritischen Linse. **(BU) mit Reserve ist eine Aussage über die Phasen der π_n.**
- **Gemeinsamer Kern:** Am Steg sind primaler Kopf × Adjungierte auf n-ter-Wurzel-Ebene genau ausgeglichen (Summe der Raten → 0).
  - Alle Betragsabschätzungen enden deshalb an Stellen mit Marge 0. Entscheidend sind die Kohärenz der Daten und die Eck-Wachstumszone nahe L ausserhalb D.
  - Für die Familie gilt: Der beste Repräsentant ist um 0.0024 je Grad kohärenter als ein Kopf (V).

## 4. Exakte Werkzeuge
1. **(CD1) Christoffel-Darboux-Identität (bewiesen, V; kubisch und Kneser auf 1e-176):**
   (z − t)·K_n(z,t) = q_{n+1}(z)·λ_n(t) − t·ν_n(z) + Ξ_n(z,t), Ξ_n(z,t) = (1/2πi)∮K_n(z,s)·(f(s) − s)·K_n(f(s),t)ds/s.
   - q_{n+1} ist der monische primale Kopf der Ordnung n+1, λ_n der adjungierte Kopf, ν_n = J_n⁻¹·1.
   - Herleitung: J S − S J = T_nM_{f−x}C_f − (T_nf^{n+1})e_nᵀ (Verschiebung erster Ordnung des Kompositionskanals).
   - Ohne Komposition ist das die Toeplitz-CD-Formel. Also gilt: Kern = integrabler Kern vom Rang 2 + Dyson-Term mit dem Faktor f(s) − s.
2. **Explizite Eckprofile (bewiesen, V):** Köpfe und Eck-Pol-Antworten nahe L folgen F_β(u) = Γ(β)⁻¹∫₀^∞t^{β−1}e^{−t−ut^{−iα}}dt, mit α = L/(4π), β_P = ½ + iα, β_A = ½ − 2iα.
   Eckprodukt P_n = ΔK_n(L,L) = aCL/(n d_n)·(1 + O(1/n)), C explizit (Herleitung Skizze, numerisch 1e-5).
3. **K19 (bewiesen, V):** sup_{K_{r₀}}|Π_Nf − f| ≤ (1 + c(1+M)Σ_{n≤N}Q_n)·dist_E(f, P_N). Die Datenseite ist klassische Approximation auf E; es gibt keine abgeschnittenen Pole.
4. **Rang-1-Ränderung, LU-Zerlegung der Köpfe (bewiesen, V).** Lemma K20: Die primale n-te-Wurzel-Rate relativ zur Ecke ist ≡ 0 auf D̄ (bewiesen, V; punktweise bis auf polare Mengen).
5. **Steg-Formel (numerisch):** ΔK_n(t,t)·m(t) → κ(t) (lokaler radialer Exponent des Kopfes, 0.87 … 0.96).
   Ausserdem K_n(t,t)·m(t) = t·(log δ_n)′(t) + Δ(t) + o(1) mit n-unabhängigem Δ. Das Kreismittel von ΔK_n(t,t) ist exakt tr J̃_n⁻¹ − tr J̃_{n−1}⁻¹ (→ 0.6215).

## 5. Messwerte (kubisch, falls nicht anders gesagt)
| Grösse | Wert |
|---|---|
| ‖J̃_N‖, ‖J̃_N⁻¹‖ | e^{0.47N} (Kompositionsüberlauf), e^{0.0118N} (Eckmode) |
| \|P_n\| | 0.535·n^{−0.106} beidseitig, bis n = 2000, Kneser gleich |
| R_P(r₀ = 0.05) = sup_{K_{r₀}}\|δ_n/δ_n(L)\| | ~ n^{−3.1} (1000 … 1800) |
| primal auf Σ₁ bei r | Rate ≈ 0.15·r je Grad (Wachstum in der Eckzone ausserhalb D) |
| Steg (Linse S): primal / adjungiert / Summe | ≈ 0.957S − 0.012 / ≈ −(0.957S − 0.012) / +0.0002, fallend (n ≤ 600) |
| Abschneide-Rest (I − T_n)Ju_n auf Γout | Rate ≈ ½·ln(\|f\|/ρ), bis +0.252 je Grad |
| p_N − f_K auf der Linse S = 0.03 / am Steg / auf Γout, Σ₂ | 1.75 bei N = 2000 (e^{0.015N}) / 1e-4 / konvergiert (V) |
| CF-Repräsentant auf der Linse S | Rate 0.90·(S − 0.0147); CF + adjungierte Steg-Rate ≈ −0.0024, unabhängig von S (V, N ≤ 2000) |
| s_E | gemessen 3.9–4.1, Satz K22: 3.96 |
| Dyson-Fixpunkt aus (CD1), linearisiert | Spektralradius > 1 ab n ≈ 58, ≈ e^{0.004n} (normunabhängig) |
| Ξ_n / Generator auf der Eckskala | 0.3 … 3.7, langsam steigend |

## 6. Gescheiterte Wege (Grund)
- Parametrix aus der Bulk-Inversen, Störung um gelöste Modelle, explizite Parametrizes (Taylor, Approximant auf E):
  Normverlust durch Kompositionstransport; primale Grössen an f(Γout) wachsen.
- Absolute Abschätzungen auf einer Kontur, auch mit gemischten Datendarstellungen: Die Adjungierte hat nur ausserhalb D Reserve, die Primale wächst genau dort.
  Am Steg gleichen sich beide exakt aus (Marge 0).
- Rekursionen in n (Volterra, Generator-Schur), absolute Fast-Orthogonalität, Positivität: Auslöschungen, Eck-Kohärenz, logarithmische Phasendrehung.
- IIKS aus dem Verschiebungsrang 2 (G40): Für die Inverse ist der Kompositionsterm nicht vom kleinen Rang (Dyson-Term in CD1).
- Wronski-Paarung auf dem Abel-Torus an der Ecke: Primaler und adjungierter Kopf sind in komplementären Sektoren homogen.
- (FA) über eine superharmonische Vergleichsfunktion: unmöglich (Laplace-Masse der Köpfe).
- Identität über den Abschneide-Rest: tautologisch und exponentiell.
- Dyson-Kontraktion aus CD1: Der Spektralradius ist grösser als 1. Auch beide Kanäle zusammen ergeben nur X = J⁻¹.
- Lokal integrabler Eckkern: Der Dyson-Term ist auf der Eckskala homogen vom Grad 0. Jede lokale Parametrix braucht beide Kanäle.

## 7. Fragen an Fachleute
1. **Riemann–Hilbert mit Komposition.** CD1 = integrabler Kern vom Rang 2 + Dyson-Term ∮K(z,s)(f(s) − s)K(f(s),t)ds/s.
   - Die Kopfbedingungen [x^k](δ∘f + mδ) = 0 (k < n) sind eine Zweipunkt-Hermite-Padé-Aufgabe (Punkte 0 und f(0), gekoppelt durch f).
   - Gibt es einen RH-Rahmen für solche Kerne, etwa 3×3 oder mit Carleman-Verschiebung bzw. konformer Verschweissung? Er müsste eine Ecke enthalten, in der beide Kanäle
     von gleicher Ordnung sind: Spiral-Kegelpunkt, Mellin-Exponent a = 1 + iL/(4π).
2. **Kohärentes Carathéodory–Fejér.** Gibt es Sätze, nach denen die Keime von Galerkin- oder Padé-artigen Fehlern Repräsentanten haben, die auf einem Gebiet jenseits der
   Konvergenzscheibe besser sind, als die Beträge der einzelnen Beiträge erlauben? Konkret (BU): Die Phasen der Paarungen π_n bringen eine Reserve von 0.0024 je Grad.
3. **Christoffel-Funktionen für Kompositionsoperatoren.** ΔK_n(t,t)·m(t) konvergiert ausserhalb D gegen den lokalen Exponenten der Köpfe
   (wie die Christoffel-Funktion gegen die Gleichgewichtsdichte). Gibt es eine Potentialtheorie (Vektor-Gleichgewichtsproblem für die Nullstellen von primalem und
   adjungiertem Kopf), die die exakte Aufhebung h_P + h_A = 0 am Steg erklärt?
4. **Abschnitte gewichteter Kompositionsoperatoren** u ↦ u∘f + m·u mit abstossendem Randfixpunkt auf dem Abschneidekreis.
   - Bekannte Ergebnisse zur gleichmässigen Invertierbarkeit auf Teilmengen fern vom Fixpunkt (statt auf der ganzen Scheibe)?
   - Bekannte Ergebnisse dazu, dass gleichmässige Konvergenz bis an den Fixpunkt scheitert?

Belege: phase3/HAUPTSATZ_KR0.md (Fassung 3), phase3/KERNPROBLEM.md §4, phase3/k_kern/BERICHT.md §26–§43, phase3/v_pruefstand/BERICHT.md.
