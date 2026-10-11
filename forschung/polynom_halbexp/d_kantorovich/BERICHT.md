# Strang D — Kantorovich adaptiv / Induktion in N

Stand: 10.10.2026, 22:40 (Runde 2 ab „# Runde 2“, Runde 3 ab „# Runde 3“, Runde 4 ab „# Runde 4“). Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.

## Ziel

Schranke für alle N: entweder Induktion N → N+1 mit summierbaren Schritten, oder ein Kantorovich-Argument in angepassten
Normen, das gleichmässig in N gilt — und dabei die Frage auf eine **präzise lineare Aussage über J_N** zurückführen.

## Bezeichnungen

- p_N mit Koeffizienten a = (a_0..a_N); F_N(a)_k = [x^k](p∘p) − 1/k!, k ≤ N; J_N = DF_N.
- L = 0.3181 + 1.3372 i (Fixpunkt von exp), |L| = 1.37456. Gewichtete ℓ1-Normen |v|_w = Σ_k |v_k| w_k;
  X_γ: w_k = |L|^k (1+k)^γ; X = X_0 (beschränkt sup_{|x|≤|L|} |v(x)|).
- q_N = Knesers Taylor-Polynom vom Grad N, q_hi = Knesers Schwanz (Koeffizienten k > N).
- b_k(a) = ([x^j] p_a^k)_{j≤N} (Spalte zum Koeffizienten a_k, k > N), B_N = (b_{N+1}, b_{N+2}, …).
- r_{N+1} = [x^{N+1}](p_N∘p_N) − 1/(N+1)!  (Rest der nächsten Gleichung),  w_{N+1} = J_{N+1}(pad p_N)^{-1} e_{N+1},
  s_{N+1} = 1/(w_{N+1})_{N+1} (Schur-Komplement), ŵ = s·w = (−J_N^{-1} b_{N+1}, 1), t_N = a^{(N+1)}_{N+1} (neuer
  Spitzenkoeffizient), c_N = p_{N+1} − pad(p_N).

## 1. Exakte Struktur (bewiesen)

**Lemma D1 (Schwanz wirkt linear).** Sei p = p_lo + p_hi eine Potenzreihe mit p_lo vom Grad ≤ N, p_hi = Σ_{k>N} a_k x^k,
und p analytisch auf einer Kreisscheibe vom Radius R > |a_0| (Konvergenz von Σ a_k p^k nahe 0). Dann gilt exakt

    [p∘p]_{≤N} = [p_lo∘p_lo]_{≤N} + Σ_{k>N} a_k [p_lo^k]_{≤N},   d. h.   Φ_N(a_lo, a_hi) = F_N(a_lo) + B_N(a_lo) a_hi.

*Beweis.* p∘p = p_lo∘p + p_hi∘p. Erstens p_lo∘p − p_lo∘p_lo = Σ_{m≤N} a_m (p^m − p_lo^m), und p^m − p_lo^m ist durch
p − p_lo = O(x^{N+1}) teilbar. Zweitens p_hi∘p = Σ_{k>N} a_k p^k (konvergent nahe 0, da |p(x)| < R für kleine |x|), und
wieder p^k ≡ p_lo^k mod x^{N+1}. ∎

*Folge D1a.* Für Knesers f_K = q_N + q_hi gilt F_N(q_N) = −B_N(q_N) q_hi exakt. q_N ist also die *exakte* Lösung des
N-Systems mit „Schwanzkraft“, p_N die Lösung ohne. Erste Ordnung: p_N − q_N ≈ J_N(q_N)^{-1} B_N(q_N) q_hi.

**Lemma D2 (Rand / Bordering).** Mit p = p_N (Grad N) und pad = (a, 0) ist

    J_{N+1}(pad) = [[J_N(a), b], [cᵀ, d]],  b = ([x^k]p^{N+1})_{k≤N},  c_j = [x^{N+1}](p^j + x^j p'∘p),  d = [x^{N+1}]p^{N+1} + p'(a_0),

und F_{N+1}(pad p_N) = r_{N+1} e_{N+1}. Das System N+1 im Punkt (a, t) (t = a_{N+1}) ist **gleichwertig** mit

    F_N(a) + t·b(a) = 0   (N+1 Gleichungen),        g_0(a) + t·g_1(a) + (N+1) a_0^N t² = 0   (eine Gleichung),

g_0(a) = [x^{N+1}](p∘p) − 1/(N+1)!, g_1(a) = p'(a_0) + [x^{N+1}] p^{N+1}.

*Beweis.* Mit p̃ = p + t x^{N+1}: p̃∘p̃ = p∘p̃ + t p̃^{N+1}. Modulo x^{N+1} ist p∘p̃ ≡ p∘p und p̃^{N+1} ≡ p^{N+1}
(Lemma D1 mit einem einzigen hohen Koeffizienten) — das ergibt die ersten N+1 Gleichungen, linear in t. Koeffizient bei
x^{N+1}: [x^{N+1}] (p + t x^{N+1})^m = [x^{N+1}] p^m + t m a_0^{m−1}, also [x^{N+1}] p∘p̃ = [x^{N+1}] p∘p + t p'(a_0);
[x^{N+1}] p̃^{N+1} = [x^{N+1}] p^{N+1} + (N+1) t a_0^N. Die Jacobi-Matrix ist die Ableitung dieser Darstellung bei t = 0. ∎

*Folge D2a (Newton-Schritt N → N+1).* Der erste Newton-Schritt von pad(p_N) aus ist δ = −r_{N+1} w_{N+1}, mit
w_{N+1} = (−J_N^{-1}b, 1)/s, s = d − cᵀ J_N^{-1} b (Schur-Komplement). Die Spitze bekommt t ≈ −r/s, der Tiefteil
ändert sich um −t·J_N^{-1} b.

## 2. Induktion N → N+1: Messungen (numerisch, d01_induktion.py, N = 1 … 160, 1200 Bit)

| N | r_{N+1} | s_{N+1} | t_N (Spitze) | \|c_N\|_X | \|w\|_X | \|ŵ\|_X = s\|w\|_X | Newton-Fehler rel. |
|---|---|---|---|---|---|---|---|
| 10 | 3.9e-7 | 5.43 | −7.2e-8 | 8.9e-6 | 23 | 1.2e2 | 3e-7 |
| 40 | −5.4e-14 | 6.39e2 | 8.4e-17 | 3.5e-9 | 6.6e4 | 4.2e7 | 2e-10 |
| 80 | −1.0e-19 | 3.86e5 | 2.6e-25 | 1.9e-10 | 1.9e9 | 7.5e14 | 1e-11 |
| 120 | 4.9e-25 | 2.33e8 | −2.1e-33 | 2.3e-11 | 4.7e13 | 1.1e22 | 3e-12 |
| 160 | −3.8e-30 | 1.41e11 | 2.7e-41 | 4.5e-12 | 1.2e18 | 1.7e29 | 7e-13 |

(X = Σ|v_k||L|^k.) Befunde:
- Der Schritt ist **linear**: p_{N+1} − pad(p_N) = −r w bis auf relativ 1e-12 (N = 160). Der Spitzenkoeffizient ist exakt
  t = −r/s bis auf diese Genauigkeit.
- s_{N+1} wächst streng geometrisch: s_{N+1}/s_N → 1.17371 (N = 159; noch langsam steigend, Grenzwert ≈ 1.1738);
  r_{N+1} ~ |L|^{−N}·(polynomial, oszillierend); t_N ~ 1.61^{−N} (vgl. Strang B: |a_N^{(N)}| ≈ 1.59^{−N});
  |w|_X ~ 1.2885^N, |ŵ|_X ~ 1.512^N; c_N ~ (lokal) 0.94^N bzw. polynomial — oszillierend.
- Folge: Der Schritt ist das Produkt einer exponentiell **kleinen** Spitze und einer exponentiell **grossen** Antwort ŵ der
  tiefen Koeffizienten. Jede Induktionsabschätzung, die |r|, s^{-1} und |J_N^{-1} b| getrennt in Normen abschätzt, verliert
  diese Auslöschung (Exponenten 1.512 gegen 1.613 müssten scharf bewiesen werden). **Sackgasse für einen Beweis in dieser Form.**
- Antwortmatrix R_N = J_N^{-1}B_N (d02_antwort.py): ρ_k = |J_N^{-1}b_k|_X/|L|^k hat ein Maximum bei k ≈ 1.45N mit
  Wert ~1.45^N (2.4e15 bei N = 100). Knesers Schwanz wird also nur durch massive Auslöschung (11 Grössenordnungen bei
  N = 100) zu p_N − q_N ~ 1e-10. Normabschätzungen von R_N sind hoffnungslos.

## 3. Affin-invariantes Kantorovich (Deuflhard) in gewichteten ℓ1-Normen (numerisch, d03_affin.py, d03b_gewichte.py)

h = α·ω, α = |J(x0)^{-1}F(x0)|_w, ω = max_{i≤j} |J^{-1}F''(e_i,e_j)|_w/(w_i w_j) (exakte Bilinearnorm auf ℓ1_w).
Startpunkt x0 = pad(p_{N−1}) (Induktionsschritt). log h ist konvex in log w ⇒ das **Optimum über alle Diagonalgewichte**
wurde global berechnet (L-BFGS, geglättetes Maximum).

| N | h (w=\|L\|^k) | h (\|L\|^k(1+k)^2) | h (ρ = 1.3) | h optimal |
|---|---|---|---|---|
| 10 | 1.5e-3 | 7.4e-2 | 1.2e-3 | 3.5e-4 |
| 40 | 1.0e-3 | 9.1e-3 | 5.3e-4 | 3.0e-5 |
| 60 | 3.6e-3 | 3.2e-2 | 1.8e-3 | 4.1e-5 |
| 80 | 0.25 | 2.1 | 0.14 | 1.2e-3 |
| 100 | 5.1 | — | 3.8 | — |
| 160 | 5.2e5 | — | 1.8e6 | — |

- ω wächst in jeder Gewichtsfamilie ≈ ×1.28/Grad (dieselbe Rate wie |w|), α fällt nur polynomial ⇒ h ~ 1.2^N.
  Selbst die optimalen Diagonalgewichte senken h nur um einen festen Faktor (~200), nicht die Rate. Das Maximum sitzt stets
  bei Paaren (e_i, e_N) mit i ≈ N/5 (Strang C fand dasselbe, 19:31/19:42).
- Reduzierte Form (Tiefteil-System F_N(a) + t b(a) = 0 bei festem t, Startpunkt p_N): gleiches ω, gleiches Wachstum.
- **Ergebnis:** Kantorovich mit Kugeln in gewichteten ℓ1-Normen trägt den Induktionsschritt nur bis N ≈ 80–110.
  Die Nichtlinearität muss auf einem Kegel kontrolliert werden (Strang C, Punkt (c)).

## 4. Hauptziel (a): gleichmässige Stabilität sup_N |J_N⁻¹| (Auftrag der Leitung, 20:00)

### 4.1 Exakte Zerlegung (bewiesen)

**Lemma D3.** Mit u = p − a_0 (u(0) = 0) gilt exakt J_N = U_N S_N + M_N, wobei
U_N = ([x^k] u^m)_{k,m≤N} (untere Dreiecksmatrix, Diagonale a_1^k), S_N = (binom(j,m) a_0^{j−m})_{m≤j≤N} (Taylor-
Verschiebung nach a_0, obere Dreiecksmatrix, Einsen auf der Diagonale), M_N = untere Dreiecks-Toeplitz-Matrix von m = p'∘p.
*Beweis.* p∘p = Σ_m c_m u^m mit c = Taylor-Koeffizienten von p bei a_0 (c = S a). Ableiten nach a: δ(p∘p) =
Σ_m (Sδ)_m u^m + Σ_m c_m m u^{m−1} δ = U S δ + (p'∘p) δ; Abschneiden ist exakt, weil U, M untere und S obere
Dreiecksmatrizen sind (u^m = O(x^m), S bildet Polynome vom Grad ≤ N auf ebensolche ab). ∎
Äquivalent: J_N δ = g ⇔ T^{a_0}_N[δ + f'·(δ∘h)] = T^{a_0}_N[g∘h] (h = f⁻¹, T^{a_0}_N = Taylor-Polynom bei a_0) — dieselbe
Struktur wie Strang Bs J_N = H(S + Γ).

### 4.2 Der Grenzoperator G = lim J_N⁻¹ (numerisch, d05_struktur.py, d06_symbol.py)

- **Eintragsweise Konvergenz:** (J_N⁻¹)_{kj} hängt für feste k, j nicht mehr von N ab (N = 60 gegen 100: 7 Stellen für
  k, j ≤ 5; im Inneren j ≈ 100: N = 160 gegen 194 auf ~1e-3).
- **Symbol:** Umskaliert G̃_{kj} = |L|^{k−j} G_{kj} ist G im Inneren eine **Toeplitz-Matrix**, g_{k−j} = Fourier-Koeffizient von

      σ(φ) = f'(x) e^{−x} · 1[ |f(x)| < |L| ],   x = |L| e^{iφ},

  d. h. „Multiplikation mit 1/m = 1/(f'∘f)“ auf dem Bogen Γin = {x ∈ ∂D : f(x) ∈ D} (von L über −|L| nach L̄) und 0 auf
  Γout. Der Rand von Γin liegt genau bei arg L = 76.618°. Abweichung Matrix ↔ Formel ≤ 7e-3 (N = 160, j = 50…90; die
  Koeffizienten selbst sind bis 0.67 gross). Ohne die Abschneidung 1_{Γin} passt nichts (1/m allein ist untere Dreiecksgestalt).
- **Deutung (Skizze):** Für δ = J∞⁻¹g gilt δ∘f + mδ = g auf V₀. Auf Γin ist f(x) ein innerer Punkt, δ∘f also glatt bis
  zum Rand ⇒ der Randsingulärteil von δ ist der von g/m. Auf Γout gilt δ(x) = g(h(x)) − m(h(x)) δ(h(x)) mit h(x) im Inneren
  ⇒ δ ist dort glatt (setzt sich sogar über D hinaus fort — passt zu Strang C 19:43, Konvergenz bis D(a_0, |L−a_0|)).
- ℓ²-Probe: |J_N⁻¹|_{ℓ²,|L|} = 1.366 für N ≤ 100 = ‖σ‖_∞ = |f'(−|L|) e^{|L|}| = 1.3663 (Toeplitz-Normformel ‖T(σ)‖ = ‖σ‖_∞).

### 4.3 Folgerung: (a) ist in ℓ¹ falsch, nur logarithmisch (numerisch + Argument)

*Argument.* Angenommen J_N⁻¹ → G eintragsweise und G̃_{j+d,j} → g_d (j → ∞, d fest) mit g_d = Fourier-Koeffizienten
von σ. Da σ unstetig ist, gilt Σ_d |g_d| = ∞ (eine Funktion mit absolut summierbarer Fourier-Reihe ist stetig). Für jede
ℓ¹-Norm mit w_k = |L|^k ρ(k), ρ(k+d)/ρ(k) → 1 (z. B. (1+k)^s) folgt mit Fatou
sup_N |J_N⁻¹|_w ≥ sup_j liminf_N |J_N⁻¹ e_j|_w / w_j ≥ Σ_{|d|≤D} |g_d| für jedes D ⇒ **sup_N |J_N⁻¹|_w = ∞.**
(Streng, sobald die beiden Voraussetzungen bewiesen sind; diese sind numerisch.)
Gemessen: Σ_{|d|≤j/10} |G̃_{j+d,j}| (max über j) = 1.48, 1.87, 2.12, 2.33, 2.58, 2.74 (N = 20, 40, 60, 100, 160, 194),
Fourier-Theorie Σ_{|d|≤D}|g_d| = 2.23, 2.47, 2.70, 2.94, 3.18, 3.42 (D = 10, 20, …, 320) ⇒ ≈ 0.34·ln D.
Bei s = 4 ist das unsichtbar (Spalte 0 dominiert mit 28.2), wird aber bei N ~ e^{80} überholt. Für Kantorovich ist ein
Faktor log N harmlos (α_N ~ N^{−1.1}); die richtige Zielaussage in ℓ¹ wäre also (a′) |J_N⁻¹|_{ℓ¹,s} ≤ C + c·log N.

### 4.4 Normvergleich ℓ¹ / ℓ² (numerisch, d07_normen.py, N = 20 … 194)

|J_N(p_N)⁻¹| mit gleichen Gewichten w_k = |L|^k (1+k)^s:

| N | ℓ¹ s=0 | ℓ¹ s=2 | ℓ¹ s=4 | ℓ² s=0 | ℓ² s=1 | ℓ² s=2 | ℓ² s=3 | ℓ² s=4 | ℓ² s=5 |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 2.09 | 2.18 | 20.8 | 1.366 | 1.435 | 1.592 | 2.247 | 6.779 | 43.7 |
| 60 | 2.77 | 2.45 | 25.4 | 1.366 | 1.435 | 1.592 | 2.247 | 6.838 | 51.5 |
| 100 | 3.79 | 2.95 | 26.8 | 1.367 | 1.435 | 1.592 | 2.247 | 6.844 | 54.8 |
| 140 | 5.10 | 3.46 | 27.6 | 1.444 | 1.435 | 1.592 | 2.247 | 6.846 | 57.0 |
| 194 | 7.84 | 4.96 | 28.2 | 1.784 | 1.552 | 1.592 | 2.247 | 6.847 | 59.0 |

- In ℓ² mit s ∈ [2, 4] ist |J_N⁻¹| bis N = 194 auf 3–4 Stellen konstant. ℓ²_{|L|,s} ist für s > ½ eine Algebra,
  liegt in ℓ¹_{|L|} (Cauchy–Schwarz) und enthält Kneser für s < 5.95.
- Das Wachstum bei kleinem s (ℓ¹ alle s ≤ 2, ℓ² s ≤ 1 ab N ≈ 120) kommt aus einer **Eckschicht** (Spalten j ≥ N − 30).

### 4.5 Die Eckschicht: Spitze bei L (numerisch, d08_ecke.py, d09_spitze.py, d10_gross.py; p_N bis N = 250 selbst fortgesetzt)

H_N = max_{j ≥ N−40} sup_{|x|=|L|} |(J_N⁻¹ x^j)(x)| / |L|^j (Eckspalten als Polynome auf dem Kreis):

| N | 40–110 | 120 | 140 | 160 | 180 | 200 | 220 | 250 |
|---|---|---|---|---|---|---|---|---|
| H_N | 1.54–1.61 (φ = 180°) | 1.66 | 2.02 | 2.47 | 3.02 | 3.69 | 4.52 | 6.26 |
| φ − arg L | — | 1.03° | 0.81° | 0.68° | 0.55° | 0.46° | 0.42° | 0.33° |

- Ab N ≈ 120 bildet sich in den Eckspalten eine **Spitze knapp innerhalb Γin neben L** (und spiegelbildlich L̄), die
  näher an L rückt (Abstand ~N^{−1.5}) und wächst: H_N − ‖σ‖_∞ ~ N^{3.72}, N^{3.49}, N^{3.37} (Fits ab N = 120, 150, 180;
  Exponent fällt mit dem Fenster — Grenzwert unklar, Kandidaten 3.447 = 4.447 − 1 oder 0.883·3.447 = 3.04, Strang B).
- Dieselbe Spitze mit **identischen** Werten bei J_N(q_N) (d11) ⇒ Eigenschaft der abgeschnittenen Linearisierung, nicht von p_N.
- Ein Fit der Eckspalten an die Störrichtungen δ_{±1} (gemeinsam/stoerrichtungen_delta.json) erklärt sie **nicht**
  (Restanteil 94–97 %). Mitte (j = N/2): stets 1.3664 = ‖σ‖_∞.
- Auswirkung auf die Normen (Eckblock = Spalten j ≥ N−30, ℓ², Gewicht |L|^k(1+k)^s):

| N | s=0 Ecke/Rest | s=2 Ecke/Rest | s=3 Ecke/Rest | s=4 Ecke/Rest |
|---|---|---|---|---|
| 200 | 1.79 / 1.366 | 1.44 / 1.592 | 1.39 / 2.247 | 1.38 / 6.847 |
| 250 | 2.41 / 1.366 | 1.66 / 1.592 | 1.52 / 2.247 | 1.44 / 6.848 |

  ⇒ |J_N⁻¹|_{ℓ²,2} wächst ab N = 240 (1.592 → 1.651 → 1.720), s = 3, 4 bis N = 250 noch konstant (2.247, 6.848).
- **Fortsetzung bis N = 400** (d14_newtoniterat.py, beim ersten Newton-Iterat q⁽¹⁾ = q_N − J(q_N)⁻¹F(q_N), das bis auf
  ≤ 1.5e-9 an p_N liegt; Kneser-Daten bis k = 400 mit rel. Fehler ≤ 1e-24; Ecke identisch mit J(q_N), d11):

| N | 250 | 300 | 350 | 400 |
|---|---|---|---|---|
| H_N | 6.26 | 10.55 | 18.17 | 32.20 |
| H_N − ‖σ‖_∞ | 4.90 | 9.18 | 16.81 | 30.83 |
| Ecke ℓ², s = 4 (Rest) | 1.44 (6.848) | 1.62 (6.848) | 1.98 (6.848) | 2.61 (6.849) |
| ℓ², s = 3 gesamt | 2.247 | 2.247 | 2.486 | 3.483 |
| ℓ², s = 2 gesamt | 1.720 | 2.216 | 3.068 | 4.496 |
| ℓ¹, s = 4 gesamt | 28.68 | 28.99 | 29.23 | 29.43 |
| 2. Newton-Schritt (ℓ¹,4) | 5.1e-13 | 8.0e-11 | 1.2e-9 | 1.5e-9 |

  H_N − ‖σ‖_∞ wächst um den Faktor 1.87, 1.83, 1.83 je 50 Grade; lokale Potenz-Exponenten d ln(H−‖σ‖)/d ln N = 3.34,
  3.45, 3.92, 4.54 (N = 200→250→…→400) steigen ⇒ **überpolynomial**. Fit ab N = 200: e^{0.0128 N} (rms 0.033) besser als
  N^{3.69} (rms 0.044); ab N = 240: e^{0.0124 N} (rms 0.017) gegen N^{3.86} (0.039). Über 140–400 insgesamt passt N^{3.6}
  ebenfalls (Anlaufphase). Asymptotisches Gesetz (exp oder gestreckt-exp, etwa e^{c N^{0.88}} mit Bs Eckwinkel-Exponent)
  offen. Lage der Spitze → arg L (76.77° bei N = 400).
- **Bewertung:** Die Eckschicht ist eine echte, robuste fast-singuläre Richtung von J_N (σ_min ~ 1.012^{−N} in den
  |L|-Normen), lokalisiert am Sprungpunkt L. |J_N⁻¹| wächst damit in jeder getesteten Norm schliesslich exponentiell
  (ℓ²_{|L|,2} ab N ≈ 240, ℓ²_{|L|,3} ab ≈ 330, ℓ²_{|L|,4} extrapoliert ab ≈ 500; ℓ¹ ohnehin). **(a) als gleichmässige
  Operatorschranke ist numerisch widerlegt** (Status: numerisch, N ≤ 400). Die Rate ist klein (1.2 %/Grad); für N ≤ 250
  ist ℓ²_{|L|,3..4} praktisch stabil.

### 4.7 Stabilität gilt bei p_N, nicht bei q_N (numerisch, d11_kneser_ecke.py, d12_qN_vs_pN.py)

| N | 100 | 150 | 194 | 250 |
|---|---|---|---|---|
| \|J_N(q_N)⁻¹\|, ℓ², s = 4 | 6.848 | 40.0 | 1520 | 6.0e6 |
| \|J_N(p_N)⁻¹\|, ℓ², s = 4 | 6.844 | 6.846 | 6.847 | 6.848 |
| \|J_N(q_N)⁻¹\|, ℓ², s = 0 | 1.367 | 1.486 | 1.836 | 3565 |

- Kein Datenfehler: p_N mit relativen Störungen 1e-20 der Koeffizienten gibt weiter 6.847 (Kneser-Daten: rel. Fehler ≤ 2e-31).
- Entlang x_t = q_N + t(p_N − q_N) (N = 194): |J(x_t)⁻¹|_{ℓ²,4} = 1520, 1140, 760, 380, 152, 6.85 für t = 0, ¼, ½, ¾, 0.9, 1 —
  **linear in 1 − t**. Also kein Faltpunkt, sondern J(x_t)⁻¹ ≈ J(p)⁻¹ + (1 − t)·J(p)⁻¹F''(p − q, ·)J(p)⁻¹ mit riesigem
  zweitem Term (Norm ~1500 bei |p − q|_{ℓ¹,4} = 1.2e-2; wächst ≈ 1.1–1.16^N).
- **Ein Newton-Schritt** von q_N aus genügt: q⁽¹⁾ = q_N − J(q_N)⁻¹F(q_N) hat |q⁽¹⁾ − p_N|_{ℓ¹,4} = 1.5e-12 (N = 194),
  5.1e-13 (N = 250) bei |q_N − p_N| = 1.2e-2, 8.9e-3; zweiter Schritt 2e-33. Bei q⁽¹⁾ ist |J⁻¹|_{ℓ²,4} = 6.847/6.848.
- **Folge für das Beweisgerüst:** Ein Kantorovich-Argument mit Zentrum q_N ist hoffnungslos (ω ≥ |J(q)⁻¹F''(δ,·)|/|δ|
  ~ 10^5 bei N = 194, exponentiell wachsend); die Stabilität (a) wird bei Punkten nahe p_N gebraucht (z. B. q⁽¹⁾). Dann
  zerfällt der Beweis in: (b) |J(q_N)⁻¹F(q_N)| klein (Strang C, Cauchy/Pol-Darstellung), (b′) Rest bei q⁽¹⁾ =
  ½F''[δ₁,δ₁] + … klein in J(q⁽¹⁾)⁻¹-Norm (Strang C: ω_eff ~ 1e-8), (a) |J(x)⁻¹| für x nahe p_N, (c) Lipschitz auf dem Kegel.

### 4.8 Starke Konvergenz auf festen Eingaben (numerisch, d13_stark.py)

J_N(p_N)⁻¹P_N v → G v trotz Eckschicht, mit polynomialer Rate, z. B. Fehler gegen N = 250 (ℓ¹, s = 0):
v = e^x: 1.5e-6, 3.6e-8, 1.6e-9 (N = 40, 100, 200); v = (1 − x/L)^{5.447+1.058i} + konj.: 3.9e-4, 4.9e-6, 1.3e-7;
in ℓ¹_{|L|,4} für v = e^x: 8.0, 3.3, 0.93, 0.42 (N = 40, 100, 200, 230 — Gv hat Knesers Singularitätstyp bei L, liegt nur
knapp in ℓ¹_4). Die Eckschicht stört also nicht, solange v genug Abfall hat — das ist die Form, in der (a) wirklich gebraucht
wird (vgl. Strang B: Π_N f_K → f_K ist die lineare Kernaussage).

**Bis N = 400 (d17_stark_gross.py, J bei q⁽²⁾ ≈ p_N, Referenz N = 400):** obwohl |J_N⁻¹| ab N ≈ 250 in allen Normen wächst
(§4.5), konvergiert J_N⁻¹P_N v glatt weiter, und die Lösung zeigt **keine** Eckspitze (sup auf |x| = |L| konstant, auch
in ±3° um arg L):

| N | 150 | 200 | 250 | 300 | 350 |
|---|---|---|---|---|---|
| v = e^x: Fehler ℓ¹_{\|L\|} | 7.0e-9 | 2.3e-9 | 8.9e-10 | 3.5e-10 | 1.2e-10 |
| v = e^x: Fehler ℓ¹_{\|L\|,2} | 1.2e-4 | 6.0e-5 | 3.2e-5 | 1.6e-5 | 7.1e-6 |
| v Kneser-Typ (1−x/L)^{5.447+1.058i}+konj.: ℓ¹_{\|L\|} | 7.2e-7 | 1.8e-7 | 5.8e-8 | 2.2e-8 | 7.1e-9 |

Raten ≈ N^{−4.8} (e^x; passt zu Bs Randfehler N^{−4.81}) bzw. N^{−5.45} (Singularitätsexponent). ⇒ Reguläre Eingaben
regen die exponentiell wachsende Eckrichtung nicht an. **Präzise lineare Zielaussage (Vermutung D-S):** Für
v ∈ {J∞ f_K, e^x, …} mit Knesers Singularitätstyp bei L, L̄ gilt |J_N(p_N)⁻¹P_N v − G v|_{ℓ¹_{|L|}} ≤ C_v N^{−4.8}
für alle N. Mechanismus (gemessen, N = 300/400): die dominante Richtung von J_N⁻¹ (σ₁ = 3.65 bzw. 8.92, gegen sonst
≤ 1.79/3.39) hat 90 % ihrer Masse auf den letzten 30 Eingabe-Koordinaten; ihr Überlapp mit (normiertem) e^x ist 4.8e-9
bzw. 3.7e-9. Die Eckrichtung wird also nur vom **obersten Teil** der Eingabe (k ≥ N − 30) angeregt; für Eingaben mit
Knesers Abfall |v_k||L|^k ~ k^{−6.45} ist das ~N^{−6.45}·1.012^N — harmlos bis N ~ 10³, bei rein exponentieller
Fortsetzung der Eckrate (N ~ 5000) aber nicht mehr (Extrapolation, ungeprüft). **[Runde 2: zurückgezogen — für kohärente Eingaben fällt der Überlapp ∝ 1/σ₁, siehe §9.4/9.5.]** Ein Beweis muss daher (i) den Eckmodus
explizit abspalten (Rang ~1–2) und (ii) zeigen, dass die tatsächlichen Eingaben ihn nur über ihren Kopfteil anregen.

### 4.9 Injektivität der Grenz-Linearisierung (Skizze mit denselben Voraussetzungen wie kn.satz.halbexp_beschraenkt)

**Lemma D4.** Sei δ analytisch auf D mit |δ(x)| ≤ C |x − L|^{−β} |x − L̄|^{−β} für ein β < 3.447 (z. B. jedes δ aus
ℓ¹_{|L|} oder ℓ²_{|L|,s}, s > ½: dort ist δ beschränkt). Gilt δ∘f_K + (f_K'∘f_K)·δ = 0 auf V₀ (gleichwertig: als
Potenzreihe bei 0, J∞δ = 0), so ist δ ≡ 0.

*Beweis (Skizze).* A = slog, T = sexp, Ω = A(D), B = A(V₀) ⊂ Ω ∩ (Ω − ½). Setze Δ(w) = δ(T(w)) / T'(w + ½) auf Ω.
Dafür braucht man T'(w + ½) = f'(T(w)) T'(w) ≠ 0, d. h. f' ≠ 0 auf D (numerisch: p_194' hat auf |x| = |L| die
Windungszahl 0, min |p_194'| = 0.346 auf dem Rand, kleinste Nullstelle |x| = 1.575 > |L|). Mit f(T(w)) = T(w+½),
f'(f(x)) = T'(w+1)/T'(w+½) wird die Gleichung zu T'(w+1)·(Δ(w+½) + Δ(w)) = 0, also Δ(w+½) = −Δ(w) auf B. Da Ω auf
jeder Höhe eine Strecke der Breite ≥ 2.067 enthält (e50, numerisch geprüft), setzt sich Δ durch Verschieben um ½ zu einer
ganzen Funktion mit Δ(w+½) = −Δ(w) fort: Δ = Σ_{n ungerade} c_n e^{2πinw}. Nahe L ist x − L ≈ C e^{Lw}; auf Ω liegt
Re w ≈ −0.238·Im w + O(1), daher |x − L| ≍ e^{−(|L|²/Im L)·y} = e^{−1.4129 y} (y = Im w → +∞), und |T'(w+½)| ≍ |x − L|.
Also |Δ| ≤ C e^{1.4129(β+1) y}. Für n ≤ −1 gibt c_n = ∫_{Strecke der Länge 1 in Ω auf Höhe y} Δ e^{−2πinw}
die Schranke |c_n| ≤ C e^{(1.4129(β+1) − 2π|n|) y} → 0, denn 2π/1.4129 = 4.447 > β + 1. Genauso bei L̄ (y → −∞) für
n ≥ 1. Also Δ ≡ 0, δ ≡ 0. ∎ (Lücken: dieselben geometrischen Eingaben wie im Eindeutigkeitssatz — Breite von Ω, Asymptotik
bei L — plus f' ≠ 0 auf D.)

*Folge.* Ist (J_N⁻¹ v) für ein v konvergent (in ℓ¹_{|L|}) gegen δ, so gilt J∞δ = v (die Zeilen von J∞ sind stetige
Funktionale auf ℓ¹_{|L|}, da [x^k] f^j ≤ C_r f(r)^j / r^k mit f(r) < |L| für r < 0.8) und δ ist durch v eindeutig
bestimmt. G := lim J_N⁻¹ ist also, wo es existiert, **die** Inverse von J∞ (kein Kern in ℓ¹_{|L|}).

### 4.10 Eckschicht = Randschicht = Induktionsschritt (numerisch)

Der Bordering-Vektor u = J_N⁻¹ b_{N+1} (Tiefteil des Induktionsschritts, ŵ = (−u, 1), §2) hat auf |x| = |L| sein Maximum
bei 77.21°, 77.08°, 76.95° (N = 160, 194, 250) — genau dort, wo die Eckspitze von J_N⁻¹ sitzt (77.30°, 77.12°, 76.95°).
Kosinus-Ähnlichkeit Eckspalte ↔ u: 0.32, 0.46, 0.58 (wachsend). Der Induktionsschritt p_{N+1} − p_N = t·ŵ, die Eckschicht
von J_N⁻¹ und Strang Bs Grenzschicht von p_N bei L (Skala N^{−0.88}) sind also dasselbe Objekt (vgl. B 20:20).

## 5. Gesamtbewertung für das Beweisgerüst (Stand 21:05)

1. **Induktion N → N+1 in Normen: Sackgasse.** Der Schritt ist exakt linear (Newton, rel. 1e-12), aber Produkt aus
   exponentiell kleiner Spitze (t ~ 1.61^{−N}) und exponentiell grosser Antwort (|ŵ| ~ 1.51^N); jede getrennte
   Normabschätzung verliert die Auslöschung. Kugel-Kantorovich pro Schritt: h ~ 1.2^N auch bei optimalen Diagonalgewichten.
2. **(a) gleichmässige Stabilität:** in ℓ¹_{|L|,s} falsch (log, Sprung-Symbol). In ℓ²_{|L|,s} bei p_N konstant, bis
   die Eckspitze bei L übernimmt: s = 2 ab N ≈ 240, s = 3 ab ≈ 330, s = 4 (6.848 bis N = 400) extrapoliert ab ≈ 500. Die
   Spitze wächst überpolynomial (≈ e^{0.0124 N} für N = 240 … 400, lokale Exponenten 3.3 → 4.5). ⇒ sup_N |J_N⁻¹| = ∞ in
   jeder getesteten Norm (numerisch).
3. **Stabilität nur bei p_N**, nicht bei q_N (|J(q_N)⁻¹|_{ℓ²,4} = 1520 bei N = 194, 6e6 bei N = 250). Die Abbildung
   x ↦ J(x)⁻¹ hat eine riesige Ableitung in Richtung p_N − q_N. Kantorovich muss bei einem Punkt nahe p_N zentriert werden
   (erster Newton-Iterat q⁽¹⁾ liegt bis auf 1e-12 an p_N, dort |J⁻¹| = 6.848).
4. **Was wirklich gebraucht wird (Rückführung, präzise Formulierung):** Sei A_N := P_N J∞(f_K) P_N = J_N(q_N) + E_N
   (E_N = Multiplikation mit q_hi'∘q_N, abgeschnitten) der Abschnitt der Linearisierung an Kneser und g := J∞ f_K =
   e^x + (f_K'∘f_K) f_K. **Exakt (bewiesen, Lemma D1):** P_N g = A_N P_N f_K + B_N(q_N) q_hi, also
   A_N⁻¹P_N g − P_N f_K = A_N⁻¹ B_N q_hi (Konsistenzfehler = Strang Cs (b)). Starke Konvergenz der Abschnitts-Inversen
   **auf dieser einen Eingabe** ist damit gleichwertig mit (b); die lineare Kernaussage ist also
   „A_N⁻¹P_N g → G g = f_K“
   (= Strang Bs Π_N f_K → f_K), nicht die Operator-Schranke (a). Numerisch gilt starke Konvergenz auf regulären Eingaben mit
   polynomialer Rate (§4.8), bis N = 400, obwohl |J_N⁻¹| dort schon wächst. Die Operator-Schranke (a) ist dafür weder
   nötig noch (wegen §4.7/Strang B, E) hinreichend. **Vorschlag für die Rückführung:** J_N(x)⁻¹ = (Eckmodus vom Rang 1–2,
   σ_min ~ 1.012^{−N}, nur vom Kopf der Eingabe k ≥ N − 30 angeregt) + (Rest, in ℓ²_{|L|,4} gleichmässig ≤ 6.85).
   Zu beweisen wären dann: (L1) gleichmässige Schranke des Rests, (L2) Kopf-Lokalisierung des Eckmodus, (L3) die
   tatsächliche Eingabe P_N g = P_N J∞ f_K hat nur einen kleinen Kopf (~N^{−5.45}, da g Knesers Singularitätstyp hat).
   Achtung: F(q_N) = −B_N q_hi allein hat einen riesigen Kopf (Strang C: |F(q_N)|_{ℓ¹} ~ 1.47^N) — nur die Summe
   A_N P_N f_K + B_N q_hi = P_N g ist regulär; deshalb muss man A_N⁻¹ auf P_N g anwenden, nicht auf F(q_N) allein.
   Bis N ≈ 10³ trägt das numerisch; darüber unklar (§4.8).
5. **Induktionsschritt bis N = 240:** weiterhin exakt linear (rel. Fehler der Newton-Vorhersage 1e-13), Schur-Komplement
   wächst unverändert ≈ 1.1737^N — der Zweig p_N verhält sich bis N = 250 unauffällig (d16).

## 6. Sackgassen (mit Grund)

- Kantorovich um q_N mit Kugeln (README): |J(q_N)⁻¹| wächst exponentiell (§4.7), zusätzlich ω_worst ~ 1.25^N (§3).
- Induktion mit getrennten Normschranken für r_{N+1}, s_{N+1}, J_N⁻¹b (§2).
- Bordering-Rekursion für |J_N⁻¹| mit Dreiecksungleichung: Rang-1-Zuwachs ε_N = |v|_* |ŵ|/|σ| ≈ 1.1–1.7 pro Schritt (fällt
  nicht), obwohl |J_N⁻¹| selbst kaum wächst — die Zuwächse drehen sich in der Eckschicht (d04_bordering.py).
- Störungsschritt p_N → p_{N+1} in Operatornorm: π_N = |J(pad)⁻¹ΔJ|_{ℓ¹,4} = 0.03 … 0.3 (N ≤ 40), für s = 5 > 1 — zu grob.
- (a) in ℓ¹ mit Polynomgewichten (§4.3).

## 7. Dateien (alle in d_kantorovich/)

- d_kern.py — p_N in arb (1200 Bit), F, J, Löser; Zwischenspeicher cache_p_1200.pkl (N = 1 … 250; 195–250 per Newton von
  pad(p_N) aus fortgesetzt).
- d01_induktion.py (+json): r, s, t, c_N, w, Newton-Fehler, N ≤ 160. d02_antwort.py (+json): ρ_k = |J⁻¹b_k|/|L|^k.
- d03_affin.py, d03b_gewichte.py (+d03_daten/, json): affin-invariantes h in Gewichtsfamilien und global optimal.
- d04_bordering.py (+json): Rang-1-Zuwachs und Störung. d05_struktur.py (+npz, png): Wärmebilder von C, M, J⁻¹.
- d06_symbol.py: Toeplitz-Koeffizienten von G gegen das Sprungsymbol. d07_normen.py (+json): ℓ¹/ℓ²-Normen N ≤ 194.
- d08_ecke.py, d09_spitze.py (+json), d10_gross.py (+json, N ≤ 250), d11_kneser_ecke.py (+json, q_N bis N = 400):
  Eckschicht und Spitze bei L. d12_qN_vs_pN.py: Stabilität bei q_N gegen p_N. d13_stark.py: starke Konvergenz.
- d14_newtoniterat.py (+json): Normen/Ecke beim ersten Newton-Iterat q⁽¹⁾ bis N = 400. d15_newton_quad.py: |δ₁|, |δ₂| von
  q_N aus. d16_induktion_gross.py: Induktionsschritt N = 160 … 240 (weiter linear, rel. 1e-13). d17_stark_gross.py:
  starke Konvergenz bis N = 400.

## 8. Nächste Schritte

- Wachstumsgesetz der Eckspitze (≈ e^{0.0124 N}) erklären (Mellin-Analyse an der einspringenden Ecke von D ∪ D(a₀,|L−a₀|) bei L, Bs Symbol
  μ^s + μ); die zugehörige fast-singuläre Richtung explizit machen und zeigen, dass die tatsächlichen Eingaben
  (F(q_N), F''[δ,δ]) sie nur polynomial anregen — das ist die Form, in der (a) im Beweis ersetzt werden muss.
- Beweisidee für das Symbol: J∞ = C_f + M; auf Γout-Frequenzen dominiert C_f (|f| > |L|), auf Γin ist C_f glättend —
  eine Pseudo-Differential-/Toeplitz-Zerlegung modulo kompakter Operatoren in ℓ²_{|L|} würde G = T(σ) + K beweisen.
- Starke Konvergenz J_N⁻¹P_N g → G g für g = J∞ f_K als präzise lineare Zielaussage mit Rate (Strang B, C).


# Runde 2 (ab 21:08): Die Eckmode — Ursache, Wachstum, Gefahr für die Vermutung

Auftrag der Leitung: (1) Theorie der Eckmode (Toeplitz-Abschnitte mit stückweise stetigem Symbol), Wachstumsrate, welchen Anteil
der Eingabe sie sieht; (2) Zerlegung J_N⁻¹ = Eckmode + beschränkter Rest, Bedingung an die Eingabe, Knesers Eingabe für alle N;
(3) Abstimmung mit A (Π_N f_K bis N = 1000) und B (inneres Problem).

## 9.1 Wachstum bis N = 1000 (numerisch, d19_sigma.py, d22_eck_pole.py)

J_N an x₀ = Π_N f_K (Strang A, gemeinsam/lin_PiN_fK_bis_1000.json; = p_N bis auf 2. Ordnung), letzte 100 Spalten der Inversen in
arb (0.85N + 300 Bit; Kontrolle: N = 300/400 gegen volle Inverse 3.616/8.866 statt 3.649/8.915 — 1 % zu tief, da ~10 % der
Masse von v₁ ausserhalb der letzten 100 Spalten liegt), |L|-umskaliert, ℓ²:

| N | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| σ₁ | 1.84 | 3.64 | 8.90 | 24.1 | 68.4 | 199 | 592 | 1779 | 5400 |
| σ₂ | 1.366 | 1.78 | 3.37 | 7.77 | 19.8 | 53.1 | 147 | 416 | 1195 |
| σ₃ | 1.3663 | 1.3663 | 1.3663 | 1.3663 | 1.3663 | 1.384 | 1.847 | 3.25 | 6.70 |

Fit ln σ = A + cN − a ln N (N ≥ 600): σ₁: **c = 0.01197**, a = 0.82 (rms 2e-4); σ₂: c = 0.01206, a = 1.42. Ab N ≥ 400/500
dieselben c auf ±0.0002. Also **σ₁,₂ ≈ C·N^{−a}·e^{0.0120 N}** — exponentiell, beide mit derselben Rate (zwei reelle
Kombinationen derselben komplexen Mode an L, L̄); σ₃ ist die nächste Mode (Anlauf ab N ≈ 700). Alle übrigen Singulärwerte:
1.3663 = ‖σ‖_∞ (Toeplitz-Teil). Die Rate 0.0120 ist **nicht hergeleitet** (Zahlenspiel e^{−Re γ} = 0.0117 — nicht behauptet).

## 9.2 Gestalt: lokal die Konkurrenz-Mode (x − L)^{1−γ} (numerisch, d18_eckmode.py)

- Ausgabe-Vektor u₁ (fast-Kernvektor von J_N), radial bei L, |u₁(L(1−r))|: lokaler Exponent d ln|u₁|/d ln r = −3.66, −3.50,
  −3.00 (r = 0.1 … 0.03), Sättigung für r ≲ 1/N; Phasensteigung d arg/d ln r = −0.98, −0.80 ⇒ **u₁ ≈ (x − L)^{1−γ}**
  (1 − γ = −3.447 − 1.058i, der Konkurrenz-Exponent von A/B), abgeschnitten bei |x − L| ~ 1/N. u₂ hat denselben lokalen Exponenten
  (−3.74 … −3.04, Phase −0.93). Auf dem Kreis fällt u₁ zur Γout-Seite sehr schnell ab (−5°: 1e-2 des Maximums, −20°: 1e-5),
  zur Γin-Seite langsam (+90°: 5 %). Global ist u₁ KEIN abgeschnittenes δ_{±1} (Fit-Rest 94–97 %, §4.5) — nur lokal.
- Eingabe-Vektor v₁ (linker fast-Kernvektor = adjungierte Mode): sitzt auf den Gleichungen k ≈ N − d* mit d* ≈ 10 (N = 200),
  20 (400), 30 (650), 40 (800); die obersten 3 Gleichungen fast unbeteiligt (|v₁| ~ 1e-10 … 1e-27). Als trigonometrisches
  Polynom V(φ) = Σ v₁,k e^{−ikφ} lebt v₁ auf **Γin** (|V| ≈ 0.6–1.5 auf 76.6° … 283.4°, ≤ 1e-7 auf Γout).

## 9.3 Ursache: abstossender Randfixpunkt — nicht das Sprungsymbol (numerisch, d20, d21, d23, d24)

1. **Toeplitz-Modell** T_N(a), a = m auf Γin, a = Λe^{iψ} konstant auf Γout (d20): stabil (σ₁ = 1.366 = ‖1/m‖_∞, N ≤ 400) für
   Λ = 10 … 10⁴, ψ = 0, ±π/2; nur ψ = π (Strecke von m(L) nach −Λ läuft an 0 vorbei) instabil. Die klassische Theorie der
   Abschnitte mit stückweise stetigem Symbol erklärt die Eckmode also **nicht** — die Komposition ist auf Γout kein konstanter Faktor.
2. **Lokaler Exponent irrelevant** (d23, Strang Bs λ-Deformation J_λ = C + λM, Konkurrenz-Exponent s_c(λ) = 1 − γ + 2 log λ/L):
   λ·σ₁ = 3.43, 3.55, 3.64, 3.78, 3.88 (N = 300; λ = 0.5, 0.75, 1, 1.5, 2), Rate 400→500: 0.00985 … 0.01007 — fast
   unabhängig von λ, obwohl Im s_c von −0.08 (fast reell!) bis −2.04 läuft. Die Mode hängt am **Kompositionsteil**.
3. **Möbius-Modelle** mit Fixpunkten L, L̄: Strang Es reelles Möbius-exp (elliptisch, |f'(L)| = 1, d21): **keine** Eckmode bis
   N = 600 (σ₁ = 1.3224 konstant). Komplexes Möbius-Modell f_μ = M⁻¹(μM), M = (x−L)/(x−L̄), f_μ'(L) = μ (d24):

| μ | e^{0.705i} | 1.05e^{0.3i} | 1.05e^{0.669i} | 1.05e^{1.2i} | 1.1e^{0.669i} | 1.172e^{0.3i} | √L = 1.172e^{0.669i} | 1.172e^{1.2i} |
|---|---|---|---|---|---|---|---|---|
| Rate c (N ≈ 250–350) | 0 (stabil) | 0.016 | 0.017 | 0.022 | 0.044 | 0.071 | 0.077 | 0.051 |

   (|μ| = 1.02: Anlauf erst ab N ≈ 300.) ⇒ **Ursache der Eckmode ist der abstossende Fixpunkt auf dem Rand (|f'(L)| > 1)**;
   |μ| = 1 ist genau die Grenze. Die Rate wächst mit |μ| (c/ln|μ| ≈ 0.36 … 0.49), hängt aber auch von der globalen Lage ab: exp hat
   beim gleichen Multiplikator √L nur c = 0.012 statt 0.077 (im Möbius-Modell ist L̄ anziehend, bei exp abstossend).
4. **Reelles kubisches Modell (d25) reproduziert Kneser quantitativ.** f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx) mit
   α + βL = (μ − 1)/(L − L̄) hat Fixpunkte L, L̄ mit f'(L) = μ, f'(L̄) = μ̄. Für μ = √L: α = 0.26463, β = 0.02238, f(0) = 0.5000
   (Kneser: a₀ = 0.4986), max_{|x|=|L|}|f| = 2.232 (Kneser 2.2327). J_N(f) (exakt, f kubisch):

| N | 200 | 300 | 400 | 600 | 800 |
|---|---|---|---|---|---|
| σ₁ kubisch | 1.85 | 3.69 | 9.09 | 70.7 | 619 |
| σ₁ Kneser (Π_N f_K) | 1.84 | 3.64 | 8.90 | 68.4 | 592 |
| σ₂ kubisch / Kneser | 1.365 / 1.366 | 1.79 / 1.78 | 3.42 / 3.37 | 20.4 / 19.8 | 154 / 147 |

   ⇒ Die Eckmode ist durch die Fixpunktdaten (L, L̄, μ = √L) und die grobe globale Form bestimmt, nicht durch Knesers
   Singularität. Für Strang B ist das ein explizites Modell des inneren Problems. Parameter-Scan (lokale Rate d ln σ₁/dN bei
   N ≈ 450; Kneser selbst hat dort 0.0100, asymptotisch 0.0120):

| Variation (sonst μ = √L, Fixpunkt bei arg L = 76.6°) | Rate |
|---|---|
| μ reell 1.1724 (Kreis wird ganz nach innen abgebildet, Γout leer) | 0 (σ = 1.5047 konstant) |
| \|μ\| = 1.05 / 1.1 / 1.1724 / 1.25 (arg μ = 0.669) | Anlauf / 0.0050 / 0.0099 / 0.0147 |
| arg μ = 0.3 / 0.669 / 1.0 (\|μ\| = 1.1724) | 0.0095 / 0.0099 / 0.0118 |
| Fixpunktwinkel 60° / 76.6° / 90° / 105° (Γout-Anteil 0.33 / 0.43 / 0.50 / 0.58) | 0.0075 / 0.0099 / 0.0117 / 0.0133 |

   **Asymptotische Raten** (d26_rate_scan.sh, N = 400 … 800, Fit ln σ₁ = A + cN − a ln N):

| Konfiguration | c | c/ln\|μ\| |
|---|---|---|
| μ = √L, Fixpunkt 76.6° (Kneser-Lage) | 0.0123 (Kneser selbst bis N = 1000: 0.0120) | 0.077 (0.075) |
| \|μ\| = 1.25, 76.6° | 0.0165 | 0.074 |
| \|μ\| = 1.35, 76.6° | 0.0220 | 0.073 |
| μ = √L, Fixpunkt 105° | 0.0151 | 0.095 |
| μ = √L, Fixpunkt 60° | 0.0104 | 0.066 |

   ⇒ **c ≈ 0.074·ln|f'(L)|** bei fester Lage (ln|√L| = Re L/2 = 0.159 ⇒ c = 0.0118 … 0.0123), wachsend mit dem
   Fixpunktwinkel (empirisch ≈ (φ/arg L)^{0.6}; φ ändert hier aber auch die globale Form des Modells). Eine geschlossene Formel für
   c∞ = 0.0120 habe ich nicht (die Probe c = ln|μ|·(½ − φ/π) = 0.0118 passt zwar für exp, wird aber vom φ-Scan widerlegt).
5. **Deutung (Skizze/Vermutung):** J = U(S + K) mit K = f'·C_h, h = f⁻¹ (§4.1). An L ist h anziehend (h'(L) = 1/μ); auf
   Funktionen ε^s (ε = x − L) wirkt K wie μ^{1−s}. Auf H²-artigen Räumen ist Re s = −½ die Randgerade, und {μ^{1−s} : Re s = −½}
   ist für arg μ ≠ 0 eine logarithmische Spirale |μ|^{3/2}e^{t·arg μ}·e^{i(…)} (t ∈ ℝ); für |μ| > 1 liegt −1 „innerhalb“
   (vgl. Cowen–MacCluer: wesentlicher Spektralradius gewichteter Kompositionen mit Rand-Denjoy–Wolff-Punkt, |w(ζ)|·|h'(ζ)|^{−1/2}
   = |μ|^{3/2}; Literaturstelle aus dem Gedächtnis, ungeprüft), für |μ| = 1 genau auf der Grenze. Approximative Eigenvektoren von K zum Wert −1 sind lokal ε^{s} mit
   1 + μ^{1−s} = 0 — genau die Konkurrenz-Exponenten, und u₁ sieht lokal so aus (§9.2). Endliche Abschnitte schneiden sie bei
   Skala 1/N ab ⇒ fast-singuläre Richtungen. Warum die Rate exponentiell in N ist und woher 0.0120 kommt, bleibt offen.

## 9.4 Kohärente Eingaben sehen die Eckmode nur polynomial (bewiesen + numerisch)

**Lemma D5 (exakt, bewiesen).** Sei x₀ ein Polynom vom Grad N, J := J_N(x₀) invertierbar, h analytisch nahe 0 (wie in Lemma D1).
Dann

    J⁻¹ P_N[ h∘x₀ + (x₀'∘x₀) h ] = P_N h + J⁻¹ B_N(x₀) h_hi .

*Beweis.* Lemma D1: [h∘x₀]_{≤N} = [(P_N h)∘x₀]_{≤N} + B_N h_hi; [(x₀'∘x₀)h]_{≤N} = [(x₀'∘x₀)P_N h]_{≤N}, da h_hi = O(x^{N+1}).
Summe = J P_N h + B_N h_hi. ∎

**Folge.** Mit der SVD J⁻¹ = Σ σ_i u_i v_iᵀ (ℓ²_{|L|}) ist ⟨u_i, J⁻¹w⟩ = σ_i⟨v_i, w⟩, und für kohärente w = P_N J∞^{(x₀)} h:

    σ_i⟨v_i, w⟩ = ⟨u_i, P_N h⟩  (natürlicher Anteil, ≤ |P_N h|)  +  σ_i⟨v_i, B_N h_hi⟩  (Eck-Anteil des Konsistenzfehlers).

Insbesondere ist ⟨v_i, w⟩ = O(N^κ/σ_i) genau dann, wenn der Eck-Anteil des Konsistenzfehlers polynomial bleibt.

Gemessen (d22; Pole k_ζ = 1/(ζ − x), ζ = L + u e^{iψ}, ψ = 87.15° = Winkelhalbierende des Aussenkeils, Strang B):

| N | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| ganze Antwort, u = N^{−0.883} | 3.82 | 4.51 | 5.08 | 5.58 | 6.03 | 6.44 | 6.82 | 7.18 | 7.51 |
| Eck-Anteil des Konsistenzfehlers, u = N^{−0.883} | 0.92 | 1.4 | 1.8 | 2.1 | 2.4 | 2.6 | 2.8 | 3.0 | 3.1 |
| dto., u = 4N^{−0.883} | 0.10 | 0.15 | 0.19 | 0.22 | 0.24 | 0.26 | 0.27 | 0.29 | 0.31 |
| dto., u = 0.05 fest | 4.4e-2 | 2.1e-2 | 1.2e-2 | 7.0e-3 | 4.5e-3 | 3.1e-3 | 2.2e-3 | 1.6e-3 | 1.3e-3 |
| dto., ζ = 2.0 | 5.4e-5 | 2.4e-5 | 8.9e-6 | 8.7e-7 | 2.6e-6 | 3.2e-6 | 2.1e-6 | 6.0e-7 | 5.5e-7 |
| Eck-Anteil von Π_N f_K (d19) | 8.1e-7 | 2.1e-7 | 7.8e-8 | 1.6e-7 | 1.3e-7 | 6.7e-8 | 8.0e-9 | — | — |

Während σ₁ um den Faktor 2900 wächst, wächst der Eck-Anteil an der Eckskala nur wie ~N^{0.7} (zuletzt sättigend), bei festem
Abstand fällt er, und Π_N f_K enthält die Eckmode praktisch nicht. Überlapp der Eckmode mit e^x (ebenfalls kohärent,
e^x = J∞(G e^x)): 5e-8, 4.8e-9, 3.7e-9 (N = 250, 300, 400) — fällt ∝ 1/σ₁, nicht wie der Kopf der Eingabe.

## 9.5 Zerlegung und Bedingung an die Eingabe

- **Zerlegung (numerisch):** J_N⁻¹ = Σ_{i ≤ r(N)} σ_i u_i v_iᵀ + R_N mit ‖R_N‖_{ℓ²_{|L|}} = ‖σ‖_∞ = 1.3663 (genau der
  Toeplitz-Wert), r = 1 (N ≈ 120 … 240), r = 2 (bis ≈ 690), r = 3 danach. Beweisen kann ich sie nicht (sie hängt an §9.3).
- **Warum sie allein nicht reicht:** R_N ist beschränkt, aber der Konsistenzvektor B_N h_hi ist in |L|-Gewichten exponentiell gross
  (Vorwärtsseite, Strang B/C); auch der Rest braucht also Kohärenz (der Γout-Inhalt wird vom Symbol σ = 0 getilgt).
- **Bedingung an die Eingabe:** nicht „kleiner Kopf“, sondern **Kohärenz**: w = P_N J∞ h mit h holomorph um Ḡ' bis auf Schlitze an
  L, L̄ mit Sprung ≤ C u^{s}. Mit der Cauchy-Darstellung (Strang B Satz B-Max / Strang C Satz 9.5) ist der Eck-Anteil dann
  ≤ ∫|Sprung|·(Eck-Anteil der Pol-Antwort) ≲ N^{0.7}·N^{−0.883(1+s)} (gemessenes Gesetz von §9.4 eingesetzt). Für Kneser
  (s = 5.447): ≲ N^{−5.0} → 0 **für alle N, falls das N^{0.7}-Gesetz der Pol-Antworten an der Eckskala hält** — und genau das ist
  (H*) eingeschränkt auf den Eck-Anteil.
- **Korrektur zu §4.8:** Die Extrapolation „ab N ~ 5000 Probleme für Kneser“ beruhte auf dem falschen Modell „Überlapp ≈ Kopf der
  Eingabe ~ N^{−6.45}“. Für kohärente Eingaben ist der Überlapp ∝ 1/σ₁ (Lemma D5 + Messung) — **zurückgezogen**.

## 9.6 Bewertung: bedroht die Eckmode die Vermutung?

**Nein, soweit messbar (N ≤ 1000).** Sie ist eine echte, exponentiell wachsende fast-singuläre Richtung von J_N (verursacht durch den
abstossenden Fixpunkt L auf dem Rand), aber die Eingaben, die im Beweis vorkommen (P_N J∞ f_K, Pol-Antworten), sind kohärent und sehen
sie nur polynomial. Sie widerlegt nur jede Beweisführung über Operatornormen (a). Was bleibt, ist dieselbe offene lineare Aussage wie
bei B und C: (H*) — Pol-Antworten im Randband polynomial beschränkt —, und die Eckmode ist der Grund, warum (H*) nicht aus einer
Normschranke folgen kann. Ein Beweis von (H*) muss die Eckmode explizit abspalten (Lemma D5 zeigt, wie: Eck-Anteil =
σ_i⟨v_i, B_N h_hi⟩) und für Pol-Schwänze ⟨v_i, B_N(k_ζ)_hi⟩ = O(N^κ/σ_i) zeigen.

**Für die weitere Arbeit (Abstimmung B):** Das innere Problem samt Eckmode lässt sich am reellen kubischen Modell aus §9.3 (4)
studieren — dort ist J_N(f) exakt, f explizit, und die Eckmode stimmt mit Kneser auf 3–5 % überein (N ≤ 800). Offene Punkte:
(i) die Rate c∞ = 0.0120 aus (μ, Fixpunktwinkel) herleiten; (ii) für Pol-Schwänze zeigen, dass der Eck-Anteil polynomial bleibt
(gemessen ~N^{0.7}, sättigend bis N = 1000); (iii) ob die Eckmode die Fortsetzung p_N selbst je stört (bis N = 250 linear und
glatt, §5 Punkt 5; Strang A bis N ≈ 1000).

### Dateien Runde 2
d18_eckmode.py (+json, npz: Singulärvektoren), d19_sigma.py (+json; σ bis N = 800), d20_toeplitz_modell.py, d21_moebius_ecke.py,
d22_eck_pole.py (+json; N bis 1000), d23_lambda_ecke.py (+d23.log), d24_moebius_mu.py (+d24_scan.log),
d25_kubik_modell.py (reelles kubisches Modell mit Fixpunkten L, L̄; Aufruf: N-Liste, μ-Ausdruck, Fixpunktwinkel).


# Runde 3 (ab 22:03): Eckteil von (H*) — Zerlegung, exakte Eckformel, Anregung nahe L

Auftrag: (1) Zerlegung J_N⁻¹ = Σ_{i≤3}σ_i u_i v_iᵀ + R_N mit ‖R_N‖ ≤ ‖σ‖_∞(1+o(1)); (2) Eckanteil ⟨v_i, B_N(k_ζ)_hi⟩ = O(N^κ/σ_i)
für Pole im Randband, zuerst im kubischen Modell; (3) Rate c. Ergänzung der Leitung: Strang C braucht nach Satz 12.6
(Poisson–Jensen) nur noch die **Anregung der Eckmoden durch Pole η nahe L** (|η − L| ≲ 0.05).

## 10.1 Exakte Eckformel (bewiesen)

Bezeichnungen: x₀ Polynom vom Grad ≤ N (kubisches Modell f, oder x₀ = Π_N f_K), J = J_N(x₀), m = x₀'∘x₀,
k_ζ(x) = 1/(ζ − x), y_ζ := J∞^{(x₀)}k_ζ = 1/(ζ − x₀(x)) + m(x)/(ζ − x), Π_N k_ζ := J⁻¹P_N y_ζ. SVD der |L|-umskalierten Inversen
J̃⁻¹ = Σ σ_i u_i v_iᵀ (reell, orthonormal); ℓ_i(w) := Σ_k ṽ_{i,k}|L|^k w_k und

    V_i(x) := Σ_{k≤N} ṽ_{i,k} (|L|/x)^k      (adjungierte Eckmode als Laurent-Polynom in 1/x).

**Lemma D6 (Eckformel).** Sind die Nullstellen x_j von x₀(x) − ζ einfach und ζ ∉ {0} ∪ {x_j = 0}, so gilt

    ⟨u_i, Π_N k_ζ⟩ = σ_i · [ m(ζ) V_i(ζ)/ζ  +  Σ_j V_i(x_j)/(x_j x₀'(x_j))  +  ℓ_i(P_N Q_ζ) ],   Q_ζ(x) := (m(x) − m(ζ))/(ζ − x).

*Beweis.* (i) ⟨u_i, J̃⁻¹w̃⟩ = σ_i⟨v_i, w̃⟩ = σ_i ℓ_i(w) (SVD). (ii) Für einen einfachen Pol: [x^k] c/(p − x) = c p^{−k−1}, also
ℓ_i(P_N[c/(p − x)]) = (c/p)V_i(p). (iii) Partialbrüche: ζ − x₀(x) = −a_d Π_j (x − x_j) ⇒ 1/(ζ − x₀(x)) = Σ_j 1/(x₀'(x_j)(x_j − x));
m(x)/(ζ − x) = m(ζ)/(ζ − x) + Q_ζ(x) mit einem Polynom Q_ζ. Linearität von ℓ_i. ∎
(Geprüft im kubischen Modell, N = 300 und 500, an vier Polen: relative Abweichung direkte Lösung ↔ Formel 1e-14 … 4e-6, d27.)

**Folgerung D6a (Lokalisierung).** Der Eckanteil einer Pol-Antwort ist σ_i mal Werte der adjungierten Eckmode V_i im Pol ζ und
in seinen Urbildern unter x₀. Da ṽ_i ein Einheitsvektor ist, gilt |V_i(x)| ≤ (1 − |L|²/|x|²)^{−1/2} für |x| > |L|; ausserhalb
der Kreisscheibe fällt V_i wie (|L|/|x|)^{N−d*}.

## 10.2 Zerlegung + Eckanteil ⇒ Pol-Schranke (bewiesen)

**Satz D7.** Sei r ≥ 0, S_N := σ_{r+1}(J̃⁻¹) = ‖J̃⁻¹ − Σ_{i≤r}σ_i u_i v_iᵀ‖ (ℓ²_{|L|}). Liegen alle Urbilder x_j von ζ ausserhalb
D̄ (|x_j| > |L|), so gilt

    |Π_N k_ζ|_{ℓ²_{|L|}} ≤ S_N · ( |m(ζ)|/(|ζ| s(ζ)) + Σ_j 1/(|x_j||x₀'(x_j)| s(x_j)) + |P_N Q_ζ|_{ℓ²_{|L|}} ) + Σ_{i≤r} |⟨u_i, Π_N k_ζ⟩|,
    s(x) := (1 − |L|²/|x|²)^{1/2},

und der letzte Summand ist durch Lemma D6 explizit gegeben.
*Beweis.* J̃⁻¹w = R_N w + Σ_{i≤r}σ_i u_i⟨v_i, w⟩ mit ‖R_N‖ = S_N; Dreiecksungleichung; |P_N[c/(p − x)]|_{ℓ²_{|L|}} =
(|c|/|p|)(Σ_{k≤N}(|L|/|p|)^{2k})^{1/2} ≤ |c|/(|p| s(p)) und die Partialbruchzerlegung aus Lemma D6. ∎
Damit zerfällt (H*) (in ℓ²; sup auf D̄ kostet √N) exakt in **Rest-Teil** (Schranke S_N, Eingabe nur ~dist^{−1/2}) und **Eck-Teil**
(Lemma D6). Mit S_N ≤ S (gleichmässig) ist der Rest-Teil polynomial; der Eck-Teil ist das innere Problem.

## 10.3 Rest-Schranke S_N (numerisch, nicht bewiesen)

σ_{r+1} nach Abspalten von r = 1, 2, 3 Eckmoden: kubisches Modell 1.3645, Kneser 1.3663 (= ‖σ‖_∞ = max_{Γin}|f'(x)e^{−x}| bzw.
das kubische Analogon max_{Γin}|1/f'(f(x))| = 1.36460, bei 150.9°) für alle N bis 1000 (σ₃ wird ab N ≈ 700 selbst Eckmode). Ein Beweis von S_N ≤ ‖σ‖_∞(1+o(1)) ist mir
nicht gelungen: er bräuchte eine Parametrix J_N⁻¹ ≈ P_N T(σ) P_N + (Rang ≤ 3) + klein, d. h. genau die Kontrolle der
Komposition auf Γout-Frequenzen und an der Ecke. Ehrlicher Stand: Rest-Schranke = **numerisch** (Strang B, D übereinstimmend).

## 10.4 Anregung der Eckmode durch Pole nahe L (numerisch; kubisches Modell bis N = 1000, Kneser bis N = 800)

Eckanteil |⟨u₁, Π_N k_η⟩| (ℓ²_{|L|}) für Pole **auf dem D-Bogen** η = |L|e^{i(arg L + τ)} und **auf dem Schlitz**
η = L + u e^{iψ} (ψ = 87.15°), kubisches Modell (d28, d29; Kneser identisch auf 1–3 %, d28_exp):

| Pol | N = 400 | 600 | 800 | 1000 | Rate 800→1000 minus Rate in L |
|---|---|---|---|---|---|
| L (Ecke) | 1.81 | 2.30 | 2.67 | 3.00 | 0 (≈ N^{0.5}) |
| τ = 0.5° | 2.43 | 3.39 | 4.25 | 5.05 | 0.0003 |
| τ = 2° | 2.86 | 3.87 | 5.05 | 6.61 | 0.0008 |
| τ = 4° | 2.82 | 4.45 | 7.11 | 11.8 | 0.0020 |
| τ = 8° | 3.14 | 6.69 | 15.3 | 37.0 | 0.0038 |
| τ = 12° | 3.58 | 9.82 | 29.5 | 92.4 | 0.0051 |
| τ = 103° (−ρ) | 8.1 | 71.6 | — | — | = σ₁-Rate 0.011 |
| Schlitz u = 0.001 | 1.57 | 1.88 | 2.07 | 2.19 | fällt relativ |
| Schlitz u = 0.003 | 1.19 | 1.26 | 1.24 | 1.19 | konstant |
| Schlitz u = 0.01 | 0.459 | 0.336 | 0.234 | 0.162 | fällt |
| Schlitz u = 0.03 | 0.041 | 0.014 | 0.0059 | 0.0031 | fällt |
| σ₁ | 9.09 | 70.7 | 619 | 5713 | 0.0111 |

Befunde:
- **WKB-Form:** Eckanteil(η) ≈ A·N^{κ}·e^{λ(τ)N}·(Pol-Abstandsfaktor) mit λ(0) = 0 (in L rein polynomial, κ ≈ 0.5) und λ(τ)
  stetig wachsend bis zur σ₁-Rate (≈ 0.011) ab τ ≳ 60°. Im Bogenparameter t = τ/103.38° von Strang C:
  λ = 0.0003 (t ≤ 0.01), 0.0008 (0.019), 0.0014 (0.029), 0.0020 (0.039), 0.0030 (0.058), 0.0038 (0.077), 0.0051 (0.116);
  nahe L etwa λ ≈ 0.04·t (Exponent zwischen 1 und 1.3 — mit den Daten nicht schärfer zu trennen).
- Nach Lemma D6 ist das das Profil der **adjungierten** Eckmode: |V₁(|L|e^{i(arg L+τ)})| ≈ e^{−(σ₁-Rate − λ(τ))N}·poly, d. h. V₁ ist
  in einer Umgebung von L exponentiell klein und erst weit auf Γin von der Grösse 1.
- **Auf dem Schlitz** (dort, wo Strang Cs Konsistenzintegral Knesers Sprung ~u^{5.447} integriert) ist der Eckanteil bis N = 1000
  polynomial (u = 0.001: ~N^{0.36}) oder **fallend** (u ≥ 0.003) — obwohl σ₁ um den Faktor 630 wächst.
- **Kneser (J_N bei Π_N f_K, d29_exp, N = 600/800/1000)** stimmt mit dem kubischen Modell auf 1–3 % überein: L 2.29/2.67/3.00,
  τ = 4°: 4.43/7.06/11.7, τ = 12°: 9.69/28.9/90.1, Schlitz u = 0.001: 1.88/2.06/2.19, u = 0.01: 0.336/0.234/0.161, u = 0.03:
  0.0140/0.0059/0.0031 (σ₁ = 68/592/5400).

**Prüfung von Cs Poisson–Jensen-Bedingung (Satz 12.6).** Cs harmonisches Mass am Schlitz nahe L hat (aus seinen Werten 14.65, 2.41, 1.11
für t₁ = 0.02, 0.1, 0.2 zurückgerechnet) die Grenzdichte ν(dt) ≈ 0.199·t^{−2.13}dt relativ zu g_E. Hinreichend für eine
polynomiale Schranke am Schlitz ist ∫ λ dν < 1. Mit dem gemessenen λ(t) (lineare Interpolation, λ ≤ 0.0111 für t ≥ 0.6,
Abschnitt bei t = 1/N): **∫ λ dν ≈ 0.11 (N = 1000)** — Faktor 9 Reserve. Bei λ ≈ 0.04t wächst das Integral nur wie N^{0.13}
(≈ 0.29 bei N = 10⁶, ≈ 0.83 bei N = 10⁹); ist der Exponent nahe L > 1/α = 1.13, konvergiert es gleichmässig in N.

## 10.5 Satzkette und was offen bleibt

1. **Bewiesen:** Lemma D5, D6 (exakte Eckformel), Satz D7 (Zerlegung ⇒ (H*) = Rest-Teil + Eck-Teil), Folgerung D6a.
2. **Bedingung (E_D)** (Eckteil von (H*), in der von C benötigten Form): Für die Pol-Antwort gilt auf ∂W
   log|⟨u_i, Π_N k_η⟩| ≤ κ log N + C + N λ(η) mit λ ≥ 0, λ(η) ≤ A·dist(η, L)^p nahe L (p ≥ 1) und sup λ < 0.068,
   sowie ∫ λ dν < 1. Dann liefert Cs Satz 12.6 mit Satz D7 (Rest-Teil) eine polynomiale Pol-Antwort am Schlitz und Satz 12.3 die
   Konsistenz α_N → 0. **Numerisch erfüllt bis N = 1000** (λ ≤ 0.011, ∫λdν ≈ 0.11, p ≈ 1 … 1.3).
3. **Nicht bewiesen:** (a) S_N ≤ ‖σ‖_∞(1+o(1)) (Rest-Schranke); (b) (E_D) für alle N — dafür müsste man die adjungierte Eckmode
   V_i als WKB-Lösung der transponierten Gleichung m W + 𝓛W ≈ 0 (𝓛 = Transferoperator von f, W = V/ζ, Lemma D6 mit beliebigem ζ)
   beschreiben und zeigen, dass ihr Profil nahe L um ≥ e^{−(c−λ)N} unter dem Γin-Wert liegt. Lokal bei L ergibt die
   transponierte Gleichung μW(ζ) + W(h(ζ))/μ = Rest, also lokale Exponenten W ~ (ζ − L)^{−2−(2k+1)γ} (adjungiert zu den
   Konkurrenz-Exponenten) — eine vollständige Analyse habe ich nicht.
4. **Rate c** (nachrangig): keine geschlossene Form; Daten §9.3: c ≈ 0.074·ln|f'(L)| bei fester Lage.

### Dateien Runde 3
d27_eck_formel.py (Prüfung Lemma D6, Terme längs Kurven; +json, npz), d28_bogen_anregung.py (Pole auf dem D-Bogen, kubisch und
exp; +json), d29_bogen_fein.py (fein nahe L + Schlitz, kubisch und exp bis N = 1000; +d28_bogen_kubisch_fein.json, d28_bogen_exp_fein.json,
d29.log, d29_exp.log).


# Runde 4 (ab 22:27): Rest-Schranke (a) und Bedingung (E_D) (b)

## 11.1 (a) Rückführung auf eine Kontraktionsaussage (bewiesen)

Sei f das kubische Modell oder (für Kneser) A_f = P_N J∞(f) P_N das saubere Abschnittsobjekt von Strang C (§12.6 dort), m = f'∘f
(für Kneser m = e^x/f'), H = H²(|x| < |L|) ≅ ℓ²_{|L|}, P_N die Projektion auf Polynome vom Grad ≤ N, C_f δ = δ∘f.

**Satz D8.** (i) Ist m ohne Nullstelle auf D̄ und 1/m analytisch auf D̄, so ist M_N := P_N M_m P_N invertierbar mit
M_N⁻¹ = P_N M_{1/m} P_N, und

    J_N = M_N (I + K_N),    K_N := P_N K P_N,    Kδ := (δ∘f)/(f'∘f).

(ii) Für jedes r ≥ 0:  S_N = σ_{r+1}(J_N⁻¹) ≤ σ_{r+1}((I + K_N)⁻¹) · max_{|x|=|L|} |1/m(x)|.

*Beweis.* (i) M_m bildet x^{N+1}H in sich ab (Multiplikation mit einer analytischen Funktion senkt den Grad nicht), also
P_N M_m (I − P_N) = 0 und ebenso für M_{1/m}; damit P_N M_m P_N · P_N M_{1/m} P_N = P_N M_m M_{1/m} P_N = P_N. Ferner
J_N = P_N C_f P_N + M_N und M_N⁻¹ P_N C_f P_N = P_N M_{1/m} P_N C_f P_N = P_N M_{1/m} C_f P_N = K_N.
(ii) J_N⁻¹ = (I + K_N)⁻¹ M_N⁻¹, σ_{r+1}(XY) ≤ σ_{r+1}(X)‖Y‖, und ‖M_N⁻¹‖ ≤ ‖M_{1/m}‖_{H→H} = sup_D |1/m| = max_{∂D}|1/m|. ∎

**Lage des Maximums (numerisch):** max_{∂D}|1/m| wird auf Γin angenommen und ist genau ‖σ‖_∞:
kubisch 1.36460 (bei 150.85°), Kneser 1.36630 (bei 180°); auf Γout höchstens 0.8528 bzw. 0.8529. Damit:

    **(a) ⟸ (a′):  σ_{r+1}((I + K_N)⁻¹) ≤ 1 + o(1)**   („I + K_N ist bis auf r Richtungen fast normnichtverkleinernd“).

Achtung für Kneser: die Aussage gilt für A_f (Multiplikator f'∘f), nicht für J_N(q_N) oder J_N(p_N) mit dem abgeschnittenen
Multiplikator x₀'∘x₀ (dort ist 1/m nicht beschränkt — Strang C: |[q'∘q]_{≤N}| bis 3e30).

## 11.2 (a′) numerisch — und warum die naheliegenden Beweiswege scheitern (d30, d31)

Kubisches Modell, ℓ²_{|L|}:

| N | σ(J_N⁻¹) | σ((I + K_N)⁻¹) = σ(Y) | ‖T_N(χ)‖ | σ(Y − T_N(χ)) |
|---|---|---|---|---|
| 100 | 1.366, 1.364, 1.3645 | 1.048, 1.002, **1.0001**, 1.0000 | 1.0000 | 1.049, 1.002, 0.999, 0.952 |
| 200 | 1.852, 1.365, 1.3646 | 1.492, 1.100, **1.0001**, 1.0000 | 1.0000 | 1.497, 1.105, 1.000, 1.000 |
| 300 | 3.698, 1.803, 1.3646 | 2.886, 1.552, **1.0001**, 1.0000 | 1.0000 | 2.890, 1.559, 1.000, 1.000 |
| 400 | 9.100, 3.449, 1.3646 | 6.913, 2.848, **1.0002**, 1.0001 | 1.0000 | 6.915, 2.853, 1.000, 1.000 |

- (a′) gilt numerisch mit σ_{r+1}((I + K_N)⁻¹) − 1 ≈ 1e-4 (r = 2). Gleichwertig: (I+K)ᵀ(I+K) − I hat ausser den Eckrichtungen nur
  Eigenwerte ≥ −2.6e-4 (N = 40) bzw. ≥ −1.1e-4 (N = 60) (hochgenau, mpmath).
- **Akkretivität scheitert massiv:** der hermitesche Teil von K_N hat 11 (N = 40) bzw. 20 (N = 60) negative Eigenwerte bis −6.8e4
  bzw. −2.8e9 (bei positiven bis 3e11). (a′) ist also eine echt nicht-normale Eigenschaft: der riesige Teil |Kδ|² überwiegt
  2Re⟨δ, Kδ⟩ fast überall, nur nicht in den Eckrichtungen.
- **Toeplitz-Parametrix allein ist falsch in der Norm:** Y − T_N(χ) hat viele Singulärwerte ≈ 1 (nicht kompakt, nicht klein).
  Nach Abzug der zwei Eckmoden sitzen ~90 % der Frobenius-Masse von D := Y − Eckmoden − T_N(χ) in den letzten 100 Zeilen/Spalten
  (63 % in den letzten 20); der Eckblock konvergiert aber in Flip-Koordinaten NICHT (Änderung 40 % von N = 300 auf 400), und der
  innere Teil hat Norm 0.55 → 0.62 (Hankel-/Mellin-Anteile an den Sprungpunkten L, L̄). Die eintragsweise Toeplitz-Beschreibung
  von G (§4.2) gilt nur schwach.

**Präzise Lücke (a).** Zu beweisen ist (a′). Eine Parametrix muss drei Teile haben: T_N(χ) im Inneren, lokale Mellin-Operatoren an
L, L̄ (Bs Symbol μ^s + μ bzw. die adjungierten Exponenten −2 − (2k+1)γ) und einen N-abhängigen Eckteil an der Abschneidekante, der
die r ≤ 3 Eckmoden enthält. Die klassische Theorie (Gohberg–Feldman, Böttcher–Silbermann: Stabilität der Abschnitte ⇔
Invertierbarkeit von T(a) und T(ã) und der lokalen Operatoren) setzt ein beschränktes, invertierbares Symbol voraus; hier ist das
Symbol von I + K auf Γout „unendlich“ (Komposition), auf Γin 1 + (glättender Teil). Ich habe keine Version dieser Theorie für
„Toeplitz + Komposition mit abstossendem Randfixpunkt“ gefunden oder hergeleitet.

## 11.3 (b) Bedingung (E_D): Stand

**Exakte Seite (bewiesen):** Nach Lemma D6 gilt für Pole ζ ausserhalb D̄ nahe L mit arg(ζ − L) ∈ (24.9°, 166.6°) (dann liegt auch das
nahe Urbild h(ζ) ausserhalb D̄, da h um arg μ = 38.3° dreht)

    |⟨u_i, Π_N k_ζ⟩| ≤ σ_i ( |m(ζ)||V_i(ζ)|/|ζ| + |V_i(h(ζ))|/(|h(ζ)||f'(h(ζ))|) + Fernterme ),

und |V_i| ist ausserhalb des Kreises subharmonisch in 1/ζ: log|V_i(ζ)| ≤ ∫ P_ζ log|Ṽ_i| (Poisson–Jensen für das Äussere; Ṽ_i =
Randwerte auf |x| = |L|). **(E_D) folgt also aus einer Aussage über das Randprofil der adjungierten Eckmode:**
log(σ_i|Ṽ_i(arg L + τ)|) ≤ κ log N + C + N λ(τ) mit λ(τ) ≤ A|τ|^p (p > 1) nahe L und λ < 0.068 sonst. Für p > 1 ist
∫ λ(τ)τ^{−2}dτ < ∞, und die Poisson-Mittel bleiben auf der Eckskala polynomial. Die Umkehrung (Eckanteile ⇒ Profil) gilt nicht
direkt: die Rückwärtsbahn h^k(ζ) dreht sich um L und läuft in D hinein.

**Warum die transponierte Gleichung allein nicht reicht.** Lemma D6 für allgemeines ζ ergibt die exakte transponierte Gleichung
m(ζ)W(ζ) + 𝓛W(ζ) + ℓ(Q_ζ) = σ⁻¹U(ζ) + T_hi(ζ) (W = V/ζ, 𝓛 Transferoperator, T_hi(ζ) = ℓ(B_N(k_ζ)_hi) = Eckanteil/σ). Eine
WKB-Lösung W ~ e^{NS} müsste ohne T_hi Re S längs der Rückwärtsbahnen von h konstant halten, also bis L — das widerspricht dem
gemessenen Abfall nach L hin. Der Abfall entsteht also aus T_hi, d. h. aus der Abschneidung selbst; die transponierte Gleichung ist
in dieser Hinsicht zirkulär (sie verknüpft Profil und Eckanteil, bestimmt aber keines von beiden).

**Numerik bis N = 1200 (kubisches Modell, d29 mit Lauf _r4; Kneser identisch bis N = 1000):** Zusatzrate λ(τ) = Rate des Eckanteils
bei τ minus Rate in L:

| τ hinter L | 0.25° | 0.5° | 1° | 1.5° | 2° | 3° | 4° | 6° |
|---|---|---|---|---|---|---|---|---|
| λ (800→1000) | 0.00027 | 0.00028 | 0.00024 | 0.00048 | 0.00077 | 0.00136 | 0.00195 | 0.00298 |
| λ (1000→1200) | 0.00024 | 0.00022 | 0.00030 | 0.00061 | 0.00094 | 0.00160 | 0.00222 | 0.00325 |

(σ₁ = 619, 5713, 54379; Eckanteil in L: 2.67, 3.00, 3.29 ~ N^{0.5}.) Unterhalb 1° liegt ein Sockel 0.0002–0.0003, verträglich mit einem
polynomialen Faktor (Δκ/N ≈ 0.25/1000). **Exponent zwischen 1° und 4°: p = 1.52 bzw. 1.44 — also > 1/α = 1.13.** Mit Cs Mass
ν(dt) ≈ 0.199 t^{−2.13}dt und λ ≈ A t^{1.5} unterhalb 1.5°: **∫ λ dν = 0.086, gleichmässig bis t → 0** (0.071 / 0.085 / 0.086 für
Abschnitt t ≥ 1e-3 / 1e-6 / 1e-12). Die Zusatzraten steigen von 1000→1200 gegenüber 800→1000 noch um 10–25 %; selbst eine
Verdopplung liesse ∫λdν < 0.2. Auf dem Schlitz: u = 0.001: 1.88, 2.07, 2.19, 2.28 (N = 600 … 1200, ~N^{0.3}); u = 0.003: 1.26,
1.25, 1.19, 1.11 (fallend).

**Präzise Lücke (b).** Zu beweisen ist: das Randprofil der adjungierten Eckmode erfüllt σ_i|Ṽ_i(arg L + τ)| ≤ C N^κ e^{λ(τ)N} mit
λ(τ) ≤ A|τ|^p, p > 1.13, und λ < 0.068. Dafür bräuchte man eine Beschreibung der Eckmode als Grenzschicht an der Abschneidekante
(die Quelle T_hi), nicht nur die transponierte Transfergleichung. Numerisch: p ≈ 1.45–1.5, sup λ ≈ 0.011, bis N = 1200.

## 11.4 Ergebnis Runde 4 in einem Satz

(a) ist auf die Kontraktionsaussage (a′) für (I + K_N)⁻¹ zurückgeführt (bewiesen, mit scharfer Konstante ‖σ‖_∞), (a′) gilt
numerisch auf 1e-4; (b) ist auf eine Profilaussage über die adjungierte Eckmode zurückgeführt, die mit p ≈ 1.5 > 1.13 und
∫λdν ≈ 0.09 bis N = 1200 numerisch mit Reserve erfüllt ist. Beide Beweise fehlen; die Lücken sind oben benannt.

### Dateien Runde 4
d30_faktor_K.py (Faktorisierung, σ((I+K)⁻¹), Toeplitz-Parametrix), d31_parametrix.py (Lage von Y − T_N(χ)), d29_bogen_fein.py mit
Lauf „_r4“ (d28_bogen_kubisch_fein_r4.json, d29_r4.log; N bis 1200).
