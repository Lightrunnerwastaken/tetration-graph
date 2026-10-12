# Kernproblem — Synthese (Strang K, 12.10.2026 01:22; Auftrag Leitung/Janis)

Kennzeichnung wie überall: bewiesen / numerisch / Vermutung. Quellen: ARCHITEKTUR §4–5, phase3/k_kern/BERICHT §1–15, Berichte W, O, M, X, G, D, B, C, E, V.

## 1. Die eine offene Aussage — gleichwertige Fassungen und ihr Zusammenhang

Alles Übrige ist bewiesen oder endlich zertifiziert ((K)⇒(F), Satz K1, K1″, H1/H2 bedingt, (K2) mit C₆ = 461, Zweig bis N = 400). Offen ist EINE lineare
Aussage über die Abschnitte J_N = T_NJ_∞|_{P_N} (und ihre Familie). Ihre Fassungen, von stark nach schwach (⇒ = bewiesene Implikation):

| Fassung | Inhalt | Status |
|---|---|---|
| (a) | sup_N‖J̃_N⁻¹‖ < ∞ in x | **falsch** (Eckmode σ₁ ~ e^{0.0118N}) |
| (S_γ) | sup_N‖Ĵ_N⁻¹‖ < ∞ in der konformen Koordinate mit Linse | num. bis N = 2000 (1.250337); ⇒ (H1) trivial |
| Cotlar–Stein | A₁, A₂ beschränkt (Paar-Skalarprodukte der Köpfe) | ⇒ (S_γ) (K7); asymptotisch verletzt (Eck-Kohärenz, K §15) |
| VE4 | ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ < ∞ für ein Diagonalgewicht W | ⇒ (S_γ) (K5(5)); num. 1.50–1.61 bis N = 200 (V) |
| **(H1)** | w-Köpfe beschränkt: τ_n = ‖c_n‖·‖r̂_n‖ ≤ C | ⇒ ‖Ĵ_N⁻¹‖ ≤ CN; num. τ_n → 0.015 (N ≤ 1000) |
| G19 / G10 | sup_γ*|𝒟_n| ≤ C_Δ (Kopf auf Cauchy-Daten) | ⇒ (H*_G) mit κ = 1; num. bis 1000 |
| G37 | Ecke: B_n = √n|𝒟_n(L)| ≤ C_A, ‖S_n‖_{W,1.05} ≤ C_S | ⇒ |𝒟_n| ≤ C/√n auf den Eckkreisen; num. |
| (H*_G) | G_N ≤ CN^κ auf γ (Linsenkurve) | ⇒ (H₀) (Satz K1, κ < 2.79; Familie κ < 1.79) |
| (H₀) = (E) | lokale Pol-Antwort polynomial auf Schlitz ∪ Γ_d | ⇒ (LP), mit Familie (H) ⇒ Hauptsatz |
| (E₀) | n-te Wurzel: G_N(L) = e^{o(N)}, glatte Daten Rate 0 ((H-I)) | Teil jeder obigen Fassung; num. (G_N(L) ~ N^{0.44}) |

**Zusammenhang:** (S_γ) ⇒ (H1) ⇒ (H*_G, κ = 1) ⇒ (H₀); (H1) ⇐ (S_γ) bis auf Faktor N. Satz K5 (LU): Ĵ_N⁻¹ = Σ_{n≤N}c_n ⊗ r̂_n, also ist jede Fassung
eine Aussage über die **Köpfe** (Gs Q_nA_nᵀ/h_n). Für den Hauptsatz genügt **(H1)** (κ = 1). Asymptotisch ist (H1) eine reine **Eckaussage**:
die adjungierten Köpfe konzentrieren sich an der Ecke (Wachstumsrate je Grad dort maximal, 0.0656 gegen 0.060 längs Γin), und (H1) heisst
„Eckrate des adjungierten Kopfes = Abklingrate des primalen Kopfes, mit beschränktem Vorfaktor“ (lokale Raten +0.0651 / −0.0653 bei N ≈ 900).
Die **Familie** ist erledigt bis auf (Fam-koh), dieselbe Ecke (Ĵ_N^P hängt nur vom Keim ab; an p_N identisch bis N = 1500).

## 2. Gescheiterte Ansätze: wo sie brechen, und was sie nicht benutzt haben

| Ansatz (Strang) | bricht an | nicht benutzt |
|---|---|---|
| Gleichmässige Stabilität in x, Operatornorm-Induktion N→N+1, Kantorovich mit Normen (C, D, B-Stab) | Eckmode e^{cN} in J̃_N⁻¹ | dass nur Auswertung auf einer Kurve nötig ist (Linse, K1) |
| Akkretivität, Herglotz, (FA)-Kontraktion (O, D) | Mellin-Symbol μ^{s−1} spiralt (|μ| > 1, arg μ ≠ 0) | — (strukturelles Hindernis an jedem drehenden Randfixpunkt) |
| Toeplitz/Szegő, Fisher–Hartwig (C, G1, K7/V2) | kein Ein-Variablen-Symbol, zwei Kanäle (G39) | die Kompositionsstruktur selbst |
| Gonchar–Rakhmanov-Orthogonalität (W §25) | Paarung nicht lokal (Komposition) | — |
| Tropik/Potentialtheorie n-te Wurzel, W12 (W, K3) | Verschiebungsinvarianz: Niveau sitzt in L | die Quelle/Normierung an der Ecke |
| Starke Ebene: geschweisste RH-Aufgabe mit Quelle (W §27, G43) | nicht gebaut (äussere Phase + Eckparametrix + Kleinheit) | — (offenes Programm, nicht widerlegt) |
| Spielzeuge: Möbius, affin, „Pol im Fixpunkt“, radial (V3, K §9.1, M) | Möbius: anziehendes L̄; affin: J nicht injektiv (e^{0.70N}); radial: lösbar, aber ohne Beweis | zwei ABSTOSSENDE Fixpunkte mit Rückkehr über Γin (globale Schleife) |
| Klassische Abschnittstheorie, eigenständiges Eckmodell (K §8, §9.1) | unbeschränkt; lokal invertierbar reicht nicht; Ecke im Gradabschnitt nicht abtrennbar | die globale Schleife |
| Reserve-Weg, „exponentiell ohne Ecke“ (K §6, §11) | Raten 0 in L (L ∈ γ) und bei glatten Daten | — |
| Parametrix aus der Bulk-Inversen (K k07/k08) | Fehlskalierung e^{0.47N} | — |
| Rekursionen in n: Volterra G22, Generator-Schur (G, V1, K §14) | Kern e^{0.14n}, Köpfe durch Auslöschung | eine stabile (Zukunfts-)Darstellung |
| IIKS aus Verschiebungsrang 2 (K §13) | kein gemeinsames Sylvester-Paar (Krylov + Toeplitz) | — |
| Absolute Fast-Orthogonalität (K §15) | Eck-Kohärenz (Laplace-Konzentration der adjungierten Köpfe) | Vorzeichen/Phasen an der Ecke |
| Störung um ein gelöstes Modell (G quadratisch) bzw. um f_K (V v18) | keine Kontraktionsreserve am Randfixpunkt; Lip e^{0.011N} | Kohärenz der Fast-Lösungen (K4) |
| Positivität/Monotonie in N (X §2.1) | logarithmische Phasendrehung | — |
| Doppelkeil (VE2) | Fehlskalierung wandert in die Inverse (e^{0.42N}) | — |

**Muster:** Fast alle Ansätze behandeln f als BELIEBIGE Abbildung mit zwei abstossenden Randfixpunkten L, L̄ (Multiplikator μ) und dem Gewicht m = f'∘f.
Benutzt wurden: μ (Mellin, Kegel κ₁), die Lage von L, die Analytizität von f_K um K_γ, die Kohärenz der Lösung (K4). **Nicht** benutzt wurde, dass f
eine Lösung von f∘f = exp ist — ausser in m = eˣ/f' und in der Lage der Kneser-Singularität (M §3a: das Verfahren funktioniert, weil sie genau in L sitzt).

## 3. Leitfrage: welche Struktur hat noch kein Ansatz ausgenutzt?

Prüfung der Kandidaten der Leitung und eigene Ergänzungen:

1. **Abel/Kneser-Koordinate (K8) — exakt konstante Koeffizienten. [bewiesen, eine Zeile; neu ausgeschrieben]** Mit Knesers Abelfunktion α (α∘f = α + ½) und
   φ := 1/(α'∘f) gilt φ∘f = m·φ (denn m = f'∘f = α'∘f/α'∘exp), also **J_∞(φu) = (1/α'∘exp)·(u∘f + u)**, d. h. in der Koordinate s = α(x):
   **J_∞ ≅ (Multiplikation)·(I + S_{½})·(Multiplikation)** mit der reinen Verschiebung S_{½}U(s) = U(s + ½). Nahe L ist α ≈ log(x − L)/(2 log μ), die
   Abelkoordinate lebt auf dem Zylinder ℂ/(2πi/L)ℤ (denn 2 log μ = L), und die Ecke ist der flache Torus ℂ/(ℤ·log μ + ℤ·2πi) = ℂ/(ℤ·L/2 + ℤ·2πi).
   Die ganze Schwierigkeit steckt dann im Gradabschnitt: Grad n ⇔ Lage s_n ≈ −a·log n/L auf dem Zylinder, a = 1/κ₁ = 1 − L/(4πi) **exakt**.
   **Gradskalierung n → λn ist (auf der Kegelskala) eine TRANSLATION um −a·log λ/L** — Vermutung: der Renormierungsoperator aus R2 ist eine Translation
   auf dem Torus, seine Eigenwerte λ^{−a}, λ^{−1}, λ^{−2a} sind Moden dieses translationsinvarianten Problems. **Nicht ausgenutzt:** keine Rechnung hat die Köpfe in dieser Koordinate angesehen. Für das kubische Modell gilt
   dasselbe lokal (Koenigs-Linearisierung bei L) — und lokal genügt nach der Universalität von D_N.
2. **log L = L und μ² = L — ja, aber als Teil von 1:** sie machen den Torus explizit (log μ = L/2), damit κ₁ = 2πi/(2πi − L/2) und a = 1 − L/(4πi).
   Im kubischen Modell gilt μ² = L nur durch Wahl; die Ecke sieht nur μ.
3. **Realität L ↔ L̄ — schwach:** benutzt nur für die Paarstruktur von D_N; eine erzwingende Symmetrie (Krein, R10) ist nicht in Sicht (V: keine indefinite
   Metrik erkennbar). Verwerfen als Hauptweg.
4. **Exakte Universalität von D_N — ja, als Werkzeug:** die Ecke ist lokal (kubisch = Kneser auf 4e-4); sie darf in JEDEM Modell mit gleichem μ gerechnet
   werden — insbesondere in der Koenigs-Koordinate, wo die Komposition exakt Dilatation ist. Mit Keil-Linse sättigt D_N (VE1) — sauberes Grenzobjekt.
5. **Ergänzung K — (H1) hängt nicht von der Amplitude ab:** τ_n = ‖Q_n‖‖r̂_n‖/h_n ist ein Quotient; Gs unsichere Amplitude B_n kürzt sich. (H1) ist eine
   Aussage über die normierten Kopfgestalten und den Pivot — also über Formen, die numerisch konvergieren (Gs Wiener-Abstand 0.006), nicht über B_∞.
6. **Ergänzung K — Herkunft aus f∘f = exp über das Gewicht:** m = eˣ/f' (K4) ist bisher nur für die Familie benutzt; im Abel-Bild ist es genau das
   Gewicht φ, das J_∞ zur reinen Verschiebung macht. 1, 5 und 6 gehören zusammen.

**Leitfrage (Antwort):** Kein Ansatz hat die Ecke in der Koordinate angesehen, in der die Gleichung konstante Koeffizienten hat und die Gradskalierung eine
Translation ist. Dort wäre (H1) — das exakte Gleichgewicht zweier Raten — die Aussage, dass zwei Lösungen eines translationsinvarianten Problems
(primaler und adjungierter Kopf) eine erhaltene Paarung haben (Wronski-artig): für Lösungen von u(s + ½) = −u(s) und v(s + ½) = −v(s) ist u·v
½-periodisch, also eine Funktion auf dem Torus, deren Mittel erhalten bleibt [elementar]. Sind die Köpfe nahe L genähert solche Eckmoden an den Lagen s_n,
so wäre τ_n ∝ 1/|Torusmittel von u_nv_n| — beschränkt, ohne dass Amplituden konvergieren müssen. **Das ist der Kandidat für die gesuchte Struktur (Vermutung).**


## 4. Nachtrag (Strang K, 12.10.2026 12:44) — der gemeinsame Kern in seiner jetzigen, schärfsten Form

Quellen: k_kern/BERICHT §26–§38, HAUPTSATZ_KR0.md (Fassung 2). Kennzeichnung wie oben.

**Was sich seit 01:22 verschoben hat.** In der K_{r₀}-Fassung ist bei f_K alles auf drei Annahmen zurückgeführt: (EP) (Eck-Pol-Antworten, polynomial),
(FA) (adjungierte Fernreserve, n-te Wurzel, ⇐ W12′) und (S_E) (Satz K22 auf dem gelesenen Satz 2 von Andrievskii–Pritsker–Varga; Geometrie von E bei f_K
von V numerisch geprüft, v55; s_E gemessen 3.9–4.1).
Die Datenseite ist klassisch (Satz K19, bewiesen). Die Primale hat in D̄ keine n-te-Wurzel-Reserve (Lemma K20, bewiesen). Die Familie hängt an (Fam-koh)^{r₀}.

**Die offenen Stellen sind eine.** (EP) über eine Parametrix (§37), die Familie an Linse/Steg und an Σ₁ (§34, §37.1) und F3 brauchen dasselbe:
eine **primale** Grösse (Eck-Pol-Antwort u_n, Kopf δ_n, Fast-Lösung P) **ausserhalb D̄**: Auf der Γin-Seite hat die Adjungierte nur ausserhalb D̄ Reserve
(Linse, Spirale Σ₁), und auf der Γout-Seite holt die Komposition die Primale nach f(Γout), ebenfalls ausserhalb D̄. Dort wächst die Primale exponentiell, und zwar so, dass keine Betragsabschätzung Reserve hat
[numerisch, kubisch, n ≤ 1000]:
1. **[FALSCH, zurückgezogen 13:07, k_kern §39.2: Steigung mit Niveau addiert; gleichartig gemessen ist die Summe ≈ 0] Steg (θ ≈ 180°), Kreise |t| = ρ(1 + s):** primal ≈ ln(1 + s) − 0.012 je Grad (k51c), adjungiert ≈ +0.016 − 0.93·s (k52). Die Steigungen in s heben sich auf
   (+0.94 gegen −0.93). Die Summe ist **≈ +0.003 für alle gemessenen s ∈ [0, 0.03]**: Kein Radius schliesst die Lücke, und nach innen wächst die Adjungierte weiter.
   V bestätigt am echten p_N: sup_Linse|p_N − f_K| = 1.7e-9 (N = 200) → 1.75 (N = 2000) (v_pruefstand §3.19–§3.23).
   (Die primale Rate auf ∂D ist laut K20 vorübergehend; ist sie asymptotisch 0, wird die Lücke grösser, nicht kleiner.)
2. **f(Γout):** die Eck-Pol-Antwort wächst mit ≈ ½·ln(|w|/ρ) je Grad, der halben Bernstein-Walsh-Rate (k56; bis +0.252 bei φ = 0), während die echte
   Fortsetzung s dort O(1) ist. Über die exakte Identität Ju_n = T_nk_L + R_n ist das wörtlich der Abschneide-Rest R_n auf dem Γout-Bogen.
   **Hebel (i) gibt daher keinen unabhängigen Griff** (BERICHT §38).
3. **Nahe L (Σ₁, primaler Wachstumssektor θ ∈ (10.6°, 170.6°) − 1.451°·ln n):** die echte Antwort folgt dem inneren Profil nur bis |ζ| ≈ 2–3. Danach folgt sie
   der äusseren WKB-Lösung mit global bestimmter Amplitude (k55c; Rate längs Σ₁ +0.0074 bei r = 0.05, k51d). Polynomial bleibt R_n nur auf der Eckskala
   |y − L| ≲ C·(ln n)/n.

**Schärfste Form (Skizze der Folgerung, Vermutung im Kern).** Die offene Aussage ist keine Grössen-, sondern eine **Auslöschungsaussage über den Rand von D**:
> **[ZURÜCKGEZOGEN 13:07, siehe Nachtrag 2]** **(KK)** Die Paarung einer primalen Grösse mit der Adjungierten längs einer Kontur, die D̄ auf der Γin-Seite verlässt, ist polynomial,
> obwohl ihr Integrand am Steg im Betrag mit ≈ +0.003 je Grad wächst, unabhängig vom Radius der Kontur.

(EP) ist numerisch gestützt (k49, k51). Für (Fam-koh)^{r₀} spricht nur indirekt, dass die Familie numerisch wie f_K läuft (V); direkt gemessen ist sie nicht.
Jeder Beweis muss diese Auslöschung benutzen. Betragsabschätzungen auf einer geeigneten Kontur, explizite Parametrizes
(Taylor, E-Approximant, innen V/aussen s) und die Identität über R_n scheitern alle an derselben Stelle (§34.3, §37, §38). **(EP) bleibt benannte Hypothese**
(HAUPTSATZ_KR0 §2.1). Vermutung, nicht geprüft: Weil die Summe der Raten im Radius konstant ist, verhält sich das Produkt primal × adjungiert dort wie
eine Funktion ohne radiales Wachstum. Ein Weg müsste also das Produkt als EIN Objekt behandeln, z. B. die Wronski-artige Paarung auf dem Torus aus §3, Leitfrage,
statt beide Faktoren getrennt abzuschätzen.

### Nachtrag 2 (Strang K, 12.10.2026 13:07) — Korrektur: am Steg gibt es keine Lücke; der Kern sitzt woanders

Quelle: k_kern/BERICHT §39 (k57, k57b, k57c) und Vs Carathéodory-Fejér-Test (v_pruefstand §3.24). [numerisch, kubisch, n ≤ 600 bzw. N ≤ 2000]
- **Punkt 1 und (KK) oben sind falsch.** Ich hatte eine Steigung der Primalen mit einem n-te-Wurzel-Niveau der Adjungierten addiert. Gleichartig gemessen
  (Steigung 450 → 600) ist am Steg die Summe h_P + h_A = +0.0002 und fallend, im ganzen Ring s ∈ [−0.01, 0.05] und für θ ∈ [85°, 180°].
- **Das Produkt ist dort sogar explizit:** ΔK_n(t,t)·m(t) → κ(t) (n-unabhängig; κ = lokaler radialer Exponent des Kopfes, 0.87 … 0.96), und
  K_n(t,t)·m(t) = t·(log δ_n)′(t) + Δ(t) + o(1) mit n-unabhängigem Δ. Das ist eine Christoffel-Darboux-artige Form (Vermutung), die gesuchte Paarung, aber am Steg.
  Exakt [bewiesen]: Das Kreismittel von ΔK_n(t,t) ist auf jedem Kreis tr J̃_n⁻¹ − tr J̃_{n−1}⁻¹ (gemessen → 0.6215).
- **Vs Test:** bester Keim-Repräsentant der Familie plus adjungierte Steg-Rate ≈ −0.0024, unabhängig von der Linsenhöhe. Also strikte Reserve für die Datenpaarung.
  Aus Vs zwei Ausgleichsgeraden folgt ein Linsenfenster S ∈ (0.0126, 0.0147), in dem Repräsentant und Adjungierte beide Reserve haben (≈ 0.001 je Seite).

**Der gemeinsame Kern, neue schärfste Form (KK′):** Die exponentiellen Hindernisse sitzen NICHT im Produkt am Steg. Es bleiben:
1. **Komposition der Daten:** Die Familie braucht einen Repräsentanten P̃ punktweise nahe f_K, wo der Integrand lebt, denn W = (f − q)∘P̃ + …. Das Gewicht der Adjungierten
   hilft dort nicht linear. Je Kopf genügt der Keim mod x^{n+1} (aus K13). Damit trifft die Bedingung nur die Köpfe n ≳ 0.95N, und es gibt laut V ein Linsenfenster.
2. **Eck-Wachstumszone:** Ausserhalb D nahe L wächst auch das Produkt (θ = 85°, s = 0.05: +0.015 je Grad; Σ₁ bei r = 0.2: +0.028). Σ₁ bis r = 0.1 ist polynomial.
   Kontur und r₁ müssen also innerhalb dieser Zone bleiben.
3. **Nur für den Parametrix-Weg zu (EP):** Die Primale an f(Γout) wächst (§38). Ein Beweis von (G-P) darf keine Werte der Primalen ausserhalb D brauchen,
   oder er muss die Christoffel-Darboux-Struktur benutzen.

**Nächster Kandidat:** eine exakte Christoffel-Darboux-Identität für die Abschnitte J̃_n, und Vs Test des Linsenfensters S ≈ 0.0137 auf der ganzen Linse.

**Nachtrag 3 (12.10. 13:19):** Die Identität ist gefunden und bewiesen (k_kern §40): (z − t)K_n(z,t) = q_{n+1}(z)λ_n(t) − tν_n(z) + Ξ_n(z,t). Der ganze Kern ist also
durch EIN Kopfpaar, die glatte Antwort ν_n und den Dyson-Term Ξ_n = ∮K_n(z,s)(f(s) − s)K_n(f(s),t)ds/s gegeben. Die offene Grösse (KK′) ist damit genau Ξ_n,
der Kompositionskanal. Er verschwindet in den Fixpunkten (Faktor f(s) − s).

**Nachtrag 4 (12.10. 13:37):** (CD1) als Dyson-Fixpunktgleichung kontrahiert für grosse n in keiner Norm (Spektralradius ab n ≈ 58 über 1; k_kern §41).
Der wachsende Anteil ist glatt, nicht an der Ecke. (CD1) bleibt eine exakte Darstellung, aber kein Kontraktionswerkzeug.

**Nachtrag 5 (12.10. 14:05): Jede lokale Parametrix an L muss BEIDE Kanäle enthalten** (k_kern §42, k60).
- Auf der Eckskala ist f(s) − s ~ d_n, der Kern ~ 1/d_n. Der Dyson-Term von (CD1) ist homogen vom Grad 0, also von derselben Ordnung wie das Kopfpaar.
  Gemessen: |Ξ|/|q_{n+1}λ_n| = 0.3 … 3.7, langsam steigend.
- Ein lokal integrabler Eckkern mit nur einem Generatorpaar (Rang-2-IIKS) existiert nicht. Das passt zu G39, wo sich in L beide Symbolkanäle treffen.
- Ein RH- oder Parametrix-Ansatz an der Ecke braucht mindestens drei Komponenten (Kopfpaar + Kompositionskanal) oder eine Verschiebung.

