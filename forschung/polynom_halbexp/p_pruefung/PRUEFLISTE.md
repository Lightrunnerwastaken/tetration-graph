# Prüfliste Strang P — Status aller Aussagen aus ARCHITEKTUR.md §2–§4

Status: **verifiziert** (Schritt für Schritt nachvollzogen; wo möglich numerisch gegengeprüft) · **Fehler** (mit Gegenbeispiel bzw.
Korrektur) · **Lücke** (genau benannt) · **nicht geprüft**. „verifiziert|X“ = Implikation geprüft, gilt unter Annahme X.
Fundstellen wie in ARCHITEKTUR.md (A/B/C/D/E §x). Prüfskripte und Logs in `p_pruefung/` (p01 … p07). Alle Zeiten 10.10.2026 (aus `date`).

## 0. Gesamturteil (Stand 23:57, ergänzt 11.10. 00:11 und 00:56)

**Ab 11.10. 00:42 ist ARCHITEKTUR.md §9 (von P gepflegt) die massgebliche Statusübersicht**; diese Liste enthält die Einzelnachweise.
Ergänzung 00:56: korrigierte Zielaussage = gemischter Weg (Rest in K' über (a′)_r, Eckteil in W̄ über (V1)–(V3)), §9.3; Reduktionen P-I–P-IV
bewiesen (§9.4); (A) ausgearbeitet (§9.7); rigoroses Zweig-Zertifikat p18 (σ-Homotopie) bis N = 45; Prüfskripte jetzt p01–p19.

Ergänzung Runde 4 (11.10. 02:30–03:57): **Satz P-V bewiesen** — Knesers sexp (Basis e) auf dem Streifen mit Fehler ≤ 1.3e-168 zertifiziert
(Daten aus eigenem Cousin-Newton p26; Defektschranken p27; quasikonforme Glättung + Paulsen–Cowgill-Identifikation p28), daraus f_K-Koeffizienten
a_0…a_640 mit Fehler·ρⁿ ≤ 2.9e-127 (p29; gemeinsam/kneser_halbexp_koeff_zert_640.json). Strang As 50-Stellen-Datei bestätigt auf 1.1e-43.
Das Cousin-Netz mit Rampe (Knotenkonstante) bleibt numerisches Werkzeug; für den Beweis war es wegen des log-Verlusts am Tripelpunkt ungeeignet.

Ergänzung Runde 3 (11.10. ab 01:10): (1) **kohärentes Zweig-Zertifikat p24** — jeder Übergang ab N = 45 mit EINEM s-Segment (Röhre um einen
Interpolations-Prädiktor, G und J nur an Knoten exakt, Rest über Majoranten), Zweig rigoros bis N = 193 (ARCHITEKTUR §9.7 (3b)); (2) **Kneser-Werte
mit Schranken gibt es in gp-tetration nicht** (nur Vorzertifikate) — Weg Satz P-V (§9.8), Baustein K1 (zertifizierte Poincaré-Funktion, p25) erledigt;
(3) Ws Korrektur p = 1.5–1.6 ⇒ K'-Ecke nicht grenzwertig (p17b); (4) Prüfungen O7g–O7j, M2/M3, Xs Schwanzfrage (unten §5). Skripte jetzt p01–p25.

Ergänzung 00:11: (H-ii)/(H-iii) auch an nicht-polynomialen Punkten des λ-Pfads numerisch bestätigt (N = 60, p09); Warnung: M_P nie per
Reihendivision (p09b). (a′) ⇒ (E_R): selbst der Rang r der Eckmoden ist nicht fest (Strang O, o07: 4 Werte > 1.01 bei N = 800–1000) —
(a′) muss mit wachsendem r formuliert werden (noch nicht geprüft, O hat nicht berichtet).

- **Hauptkette (H) ⇒ Hauptsatz ist logisch korrekt.** (K) ⇒ (F) [L8.2–8.4] ⇒ 4.1 ⇒ H1 ⇒ H2 ⇒ p_N → f_K: jeder Schritt nachvollzogen;
  die benutzten exakten Identitäten (Sekanten-Identität Q = λΠ^P f an NICHT-polynomialem P, H2-Identität, B §9.5(ii)) am echten
  Problem (f_K aus 4000 Koeffizienten, N = 22, 70 Stellen) auf 1e-50 bestätigt (Grenze = Genauigkeit der Kneser-Daten).
  Offen bleiben genau die benannten Annahmen (K), (H), (A), b35b′ — plus zwei kleine Präzisierungen (§1: Zweig-Definition,
  Kopplung (A) ↔ Konstanten von (H)).
- **Die Reduktion „(a′) + (E_L) ⇒ (H₀)“ ist NICHT bewiesen** (BRETT 23:46): der A-Bogen ∂D(a₀,R₀)∖D liegt ganz in f_K(D);
  dort sind die Eingaben exponentiell gross (Lücke in (a′) ⇒ (E_R)), und der Urbild-Term der Eckformel D6 ist durch (E_L) nicht
  kontrolliert (Lücke in (E_L) ⇒ (E_K)); auf dem D-Bogen liefert der Urbild-Term nur eine in |τ| lineare Randrate (nach C 13.4(a)
  zu schwach für alle N). Korrekturvorschlag: Kompaktum K' = D̄ ∪ f_K(D̄) (Abschnitt 3).
- Kleinere Fehler: „(E) ≡ (H₀)“ ist nur „(E) ⇒ (H₀)“; Bemerkung in 4.1 („frühere Form ⇒ (H₀) mit ε = 1 − c/d“) falsch auf dem
  Schlitz (braucht Poisson–Jensen, c < 0.068); Faktor 3 in einer Schwanzschranke von satz_k_beweis.py (wirkungslos, 1e-70).

---

## 1. Hauptkette (ARCHITEKTUR §4)

| ID | Aussage | Status | Zeit | Bemerkung / Prüfung |
|---|---|---|---|---|
| 4.1 | (F) + (H₀) ⇒ sup_{D̄}|Π_N f_K − f_K| ≤ C'(ε^{−4.81}N^{κ−4.81} + N^κe^{−εNd}) | **verifiziert|(F),(H₀)** | 23:57 | Cauchy-Darstellung (Π_N hängt nur von N-Jets in 0, a₀ ab; kleine Kreise um L → 0, da E_N(·,z) bei L beschränkt, z ≠ L); Schlitzintegral ∫u^{4.447}e^{−εNc_g u^{1/α}}du = αΓ(5.447α)(εNc_g)^{−5.447α} ✓ (5.447·0.882964 = 4.8095); dist(L(1+u),Ḡ') ≥ 0.93ρu, g(L(1+u)) ≥ 0.829u^{1/α} (u ≤ 0.2, p01). Existenz eines d mit Γ_d im Holomorphiegebiet folgt aus L8.4 (W offen); der Wert d = 0.1 ist nur numerisch. Übertrag auf Π^P_N: wörtlich richtig (nur Cauchy-Darstellung nötig, nicht einmal L7.2). |
| 4.1-Bem. | „Bandschranke mit e^{cN} fern der Ecke ⇒ (H₀) mit ε = 1 − c/d“ | **Fehler (klein)** | 23:57 | Gilt auf Γ_d, nicht auf dem Schlitz: Bernstein–Walsh von {g = 1/N} trägt e^{cN} bis an den Schlitz. Richtig mit Poisson–Jensen (C 12.6): ε = 1 − c·sup_ζ ω_ζ(fern)/g(ζ), braucht c < 1/14.65 = 0.068. Hauptkette unberührt ((H₀) wird direkt angenommen). |
| H1(1) | a-priori: (λ,Q) ∈ 𝒵 ⇒ Q = λΠ^P_N f, P − f = λ(Π^P_N f − f) | **verifiziert** | 23:57 | Algebra ✓; numerisch p03 T1: λ = 0.4 und 1 (P = (1−λ)f + Q nicht polynomial), max|Q − λΠ^P f| = 3.6e-50. |
| H1(2) | lokal: IFT mit D_Q = J_N(P), Nachbarlösungen bleiben in B_N | **verifiziert|(H-i)** | 23:57 | strikte Innenlage (K₀ > 1) + Stetigkeit von (λ,Q) ↦ sup_{D̄}|Q − λf|. |
| H1(3) | 𝒵 kompakt, π: 𝒵 → [0,1] Überlagerung | **verifiziert** | 23:57 | 𝒵 abgeschlossen + beschränkt in [0,1]×P_N; eigentliche lokale Homöomorphie ⇒ endliche Überlagerung des (offen-abgeschlossenen, nichtleeren: (0,0) ∈ 𝒵) Bildes [0,1]. |
| H1(4) | ein Blatt (Faser über 0 = {(0,0)}) | **verifiziert|(H-i)** | 23:57 | braucht T_N𝒥_P invertierbar an jedem (0,Q) ∈ 𝒵 — durch (H-i) gedeckt. Realität (iii) ✓, Eindeutigkeit in B_N ✓. |
| H2 | σ-Homotopie: p_{N₀} ∈ B_{N₀} ⇒ p_N = p*_N ∀N ≥ N₀ | **verifiziert|(H),(A)** | 23:57 | Identität c_M(T_M𝒥_{P₀})⁻¹x^M = P₀ − Π^{P₀}_M f numerisch 5e-50 (p03 T2); Fortsetzungsargument: ODE P' = −c_M J_M(P)⁻¹x^M, J_M⁻¹ gleichm. beschränkt auf der kompakten Menge 𝒮_M ∩ B_M ∩ {λ=1}, Schranke 5C'N^{κ−4.69} < K₀C'M^{κ−4.69} für K₀ = 6, N gross ✓. |
| H2-Def. | „Zweig“ = σ-Homotopie | **Präzisierung** | 23:57 | Newton-Definition (README) ist damit nicht bewiesen. Brücke (neu, trivial): jede Methode, deren Grenzwert eine Lösung in B_M ist, liefert p*_M (Eindeutigkeit aus H1). Für Newton reicht also „Newton ab pad p_N konvergiert und bleibt in B_M“ (numerisch: Schritt ~1e-13 ≪ R_M). |
| (A) | p_{N₀} ∈ B_{N₀} | **Lücke (Kopplung), ausgearbeitet** | 00:54 | Siehe ARCHITEKTUR §9.7: C' ≈ 12.7ε^{−4.81}C (+ Konturteil, mit grossem d vernachlässigbar ab N ≈ 30–40), R_N ≈ 240·C·N^{κ−4.81}; (A) = [(H) ab N₀ mit expliziten Konstanten] + [Zweig-Zertifikat bis N₀ (p18, rigoros)] + [zertifizierte f_K auf D̄ auf ~1e-6 bei N₀ = 40]. |
| 4.4 | Hauptsatz | **verifiziert|(K),(H),(A),b35b′** | 23:57 | Kette vollständig; Rate N^{κ−4.81} (4.69 als sichere Schranke). |

## 2. (K) ⇒ (F)  (B §8.3)

| ID | Aussage | Status | Zeit | Bemerkung |
|---|---|---|---|---|
| (K) | Knesers Struktursatz | **nicht geprüft** (Literatur) | | Strang X; numerisch gestützt (A). |
| L8.2 | f_K = φ(√LΨe^{Lε(q)}), q = Ψ^γ | **verifiziert|(K)** | 23:57 | χ = id+θ injektiv für Im w gross, ε(u) = χ(χ⁻¹(u)+½) − u − ½ 1-periodisch → 0 ⇒ nur positive Moden ✓; erste Ordnungen ε = 2κ₁q + 4πiκ₁²q² + … selbst nachgerechnet (aus ε = ϑ(u) − ϑ(u+½+ε)). μ = f_K'(L) = e^{L/2} auf 12 Stellen aus der Reihe (p01). |
| L8.3 | Fortsetzung auf D(L,r)∖Σ_L, Sprung ≤ C_J|ζ−L|^{5.447} | **verifiziert|(K)** | 23:57 | Ufer: q₋ = q₊e^{±2πiγ}, |e^{2πiγ}| = e^{−2π·1.058}; Schranke C|Ψ|(|q₊|+|q₋|) ✓. r ist durch (K) nicht quantifiziert. |
| L8.4 | f_K holomorph um Ḡ' bis auf Schlitze | **verifiziert bis auf b35b′ und r** | 23:57 | (ii) Halbmond ⊂ {Re z > 0} (selbst gezeigt: Re z ≤ 0, |z| > ρ ⇒ |z−a₀|² > R₀²), Subharmonizität von |Log z|² ✓, b35 (arb) gelesen und plausibel; innerer Bogen analytisch ≤ 0.979 ✓. (iii) b35b′ unabhängig mit dem **Motor** bestätigt (p04 M2: min|f_K(e^z)| = 0.7345, Abstand arg zu ±π ≥ 1.814 rad) — weiterhin nur numerisch. **Kleine Lücke:** Überlapp der Eckstücke braucht r(K) > 0.02 (Radius in b35) — nur numerisch (A: Formel bis |ζ−L| ≈ 0.14). |
| A-Daten | Kneser-Koeffizientendatei (4000) | **verifiziert (numerisch, Motor)** | 23:52 | f_K(x) = sexp(slog x + ½) an 7 Punkten in D, bis |x| = 0.98ρ: Abweichung ≤ 6e-32 (Motor dps 40) (p04 M1). |

## 3. Reduktion von (H₀) auf (a′), (E_L)  (ARCHITEKTUR §0, §3, §5.2)

| ID | Aussage | Status | Zeit | Bemerkung |
|---|---|---|---|---|
| Geometrie | A-Bogen ⊂ f_K(D); D-Bogen ∩ f_K(D̄) = ∅ | **neu, numerisch** | 23:46 | Windungszahl von f_K(∂D): 1 um A-Bogen-Punkte, 0 um D-Bogen-Punkte; Urbild |h(η)|/ρ = 0.990 (|η−L| = 0.058) … 0.826 (Scheitel); f_K(Γout) ausserhalb D̄(a₀,R₀) (p01). Motor bestätigt f_K = exp∘h auf f_K(Γout) (p04 M3, 1e-21). |
| (a′)⇒(E_R) | über D8 | **Lücke** | 23:46 | Auf dem A-Bogen |[1/(η − f_K)]_{≤N}|_{ℓ²_ρ} = 1.4e4, 2.0e8, 3.9e16, 7.8e24 (N = 50…300, Scheitel, Rate 1.2104 = ρ/|h|) (p02) ⇒ σ_{r+1}·|w_η| wertlos. D7 setzt selbst „Urbilder ausserhalb D̄“ voraus. |
| (E_L)⇒(E_K) | D6 + Poisson aussen + C13.2 | **Lücke** | 23:46 | (i) A-Bogen: D6-Term V_i(h(η)) mit h(η) ∈ D (V_i wächst innen wie (ρ/|x|)^N); D selbst nur arg(ζ−L) ∈ (24.9°,166.6°), A-Bogen verlässt L unter 7.7°. (ii) D-Bogen: V_i(h(η)) nur über äusseres Poisson-Integral, das bei positivem Fernprofil einen in dist(h(η),∂D) ∝ |τ| LINEAREN Beitrag gibt ⇒ nach C 13.4(a) Θ_N ~ N^{0.12} (kritisch erst ab N ~ 1e9–1e10, aber nicht „für alle N“). Ds Satz „Poisson-Mittel auf der Eckskala polynomial“ stimmt so nicht: N·dist·I/π mit dist ~ N^{−0.883} ergibt e^{c N^{0.117}}. |
| Korrektur | K' := D̄ ∪ f_K(D̄) statt W̄ | **Vorschlag** | 23:46 | ∂K' = Γin ∪ f_K(Γout); Aussenkeil bei L (24.93°, 166.62°), 141.69°, α' = 0.78717, Schlitzrate N^{−4.288} (Schwelle κ < 4.29); alle nahen Urbilder von η ∈ ∂K' auf ∂D oder aussen ⇒ |w_η|_{ℓ²_ρ} ≤ C√N ⇒ (a′) ⇒ (E_R') wie D7; g_{K'} ≤ g_W ⇒ (H₀') schwächer als (H₀); f_K auf Umgebung von K' bis auf Schlitze holomorph (f_K = exp∘h auf f_K(D̄), p04 M3). (E_L) dann beidseitig (Γin- und Γout-Seite) und vorzeichenbehaftet; Teil (ii) bleibt — sauberstes Ziel ist Cs (E_K^∂) bzw. Θ_N (C 13.4) direkt. |
| Daten (kubisch) | Grösse der D6-Terme, Profil beidseitig | **numerisch, TEILWEISE FALSCH (P-Fehler, korrigiert 00:37)** | 00:04 | p10/p10b: Der A-Bogen-Wert +0.031 (N = 1000) war float-Horner-Rauschen (W 00:40; bei N = 400 stimmt float −0.0152 mit Ws exaktem −0.0153). Richtig bleiben: D-Bogen-Urbildterm −0.004 … −0.014 (Horner aussen stabil). Γout-Profil: Werte sind ein subexponentieller Sockel, asymptotisch λ = 0⁻ (W), nicht −B|τ|^q. |
| K'-Rahmen | Green-Funktion, (a′)_r ⇒ (E_R'), Θ'_N, Eck-Kriterium | **gerechnet/bewiesen** | 00:42 | p13/p13b/kprime (cap 1.6556, κ₀' 5.817, Exponent 1.27018 vs 1.27037), p14/p15 (Θ'_N bis N = 800 ≤ 0.03), p15b (konservativ ≤ 0.110 bis 1000), p16 (Eingaben ≤ 12, Univalenz 0.3456), p17 (K' grenzwertig bei p = 1/α'). Sätze P-I–P-IV in ARCHITEKTUR §9.4. |
| (E)≡(H₀) | §5.2 | **Fehler (klein)** | 23:57 | Bewiesen ist nur (E) ⇒ (H₀) (über |E_N| ≤ |y_ζ| + (ρ/|ζ|)^{N+1}/dist und g ≤ log(|ζ|/ρ)); (H₀) liefert (E_R) auf ∂W nicht. |

## 4. Einzelsätze der Stränge

### B
| ID | Status | Zeit | Bemerkung |
|---|---|---|---|
| B1.1 J_N = H(S+Γ), Γ_ii, det | verifiziert | 23:57 | Γ_ii = p'(a₀)a₁^{−i} (Dreiecksstruktur), = a₁^{−i−1} nur an Lösungen. |
| B1.2 J-Projektion (a)–(d) | verifiziert | 23:57 | (d) = Sekanten-Identität, numerisch p03 T1. |
| B1.3 lineare Eindeutigkeit | verifiziert|(G2) | 23:57 | Fourier-Schranke ✓ (2π/κ = Re γ). Nebenweg. |
| B1.4–1.6 | verifiziert (Algebra/Standard) | 23:57 | nur Nebenweg. |
| B-lin | verifiziert|(G2) | 23:57 | Nebenweg. |
| B-Stab | verifiziert (Implikation), Voraussetzung falsch | 23:57 | leer. |
| L7.1 = C12.1 | verifiziert | 23:57 | p03 T7: Formel ζ^{−N−1}Σz_m(ζ−a₀)^{−m−1} trifft A_f⁻¹B_N[k]_{>N} auf 1e-72; |E_N|·|ζ|^{N+2} → const. |
| L7.2 = C12.2 | verifiziert | 23:57 | Standard (Subharmonizität + Abfall O(ζ^{−N−2})). |
| L7.3 Green explizit | verifiziert | 23:40 | p01: Aussenwinkel 158.9335133°, α = 0.8829639627, cap = 1.516889, g/u^{1/α} = 0.6553 (Halbierende), 0.9193 (radial), u₀(0.1) = 0.1512, e^{−g(2.2)} = 0.8049. |
| B-Max (§7.3) | verifiziert|(H*) | 23:57 | überholt durch 4.1; (H*) global für grosse N falsch (C). |
| B-Sek | verifiziert | 23:57 | a posteriori. |
| L8.1, B-Q | verifiziert|(H*_n^{loc}) | 23:57 | überholt durch H1. |
| B-Id | verifiziert | 23:52 | p03 T4: beide Gleichungen und Formel 9.0.1 auf 1e-70 mit Zufallspolynomen und f = f_K (unabhängig von b40). |
| Kor. B-Ξ₂ / d^{(2)} | verifiziert | 23:52 | p03 T5: direkte d^{(2)} gegen Formel 2e-58 (Kopfanteil 4.6e-20 bei N = 22). |
| L9.1 | verifiziert | 23:57 | von Hand: W∘exp = ga/2 = eˣW, ∂_t p_t = f'W − W∘f = a. |
| B-NF | Skizze (nicht geprüft) | | als Skizze markiert. |
| B9.5(i)–(iii) | verifiziert | 23:52 | p03 T3: x₁(N+1) − pad x₁(N) = −ρ^f_N t_{N+1} auf 3e-74; c_{N+1} = ρ^f_N (−6.87e-10). |
| B-Zweig | verifiziert|(A),(B),(C) | 23:57 | überholt durch H2. |

### C
| ID | Status | Zeit | Bemerkung |
|---|---|---|---|
| C9.1–9.4 | verifiziert | 23:57 | 9.2 (q∘q = f∘q − t∘q), 9.3 (geometrische Reihe, Majorante), 9.4 (Fubini) ✓. |
| C9.5, C12.3 | verifiziert|(F),(H_R) bzw. (F),(S) | 23:57 | Integralformel αΓ(α(σ+1))λ^{−α(σ+1)} ✓, Faktor |L|/π ✓, c_g = 0.829 (u₀ = 0.2) ✓ (p01). |
| C10.1, C10.4 | verifiziert | 23:57 | Terme erster Ordnung in t heben sich ✓. |
| C12.5 | verifiziert|(S^f) | 23:57 | |
| C12.6 | verifiziert | 23:57 | Konstanten (p11, 00:09): min g/ω(D-Bogen) auf ∂(W+0.1) ohne |ζ−L| < 0.25 = 0.0659 (C: ≥ 0.0649 ✓); Schlitz ω(D-Bogen ab |s| ≥ s₁)/g → 1/(κ₀s₁) für u → 0 und fällt in u (14.62 bei s₁ = 0.0109 ✓). |
| C13.0 | verifiziert | 23:57 | |
| C13.1 | verifiziert|(F),(E) | 23:57 | (E_R) wird auf GANZ ∂W gebraucht (Bernstein–Walsh) — siehe Lücke §3. |
| C13.2, C13.3 | verifiziert | 23:57 | 13.3: Gleichheit braucht c(ζ)ζ^{N+2} ↛ 0 ✓. |
| C13.4 (a) Grenzwertformel, κ₀ | verifiziert | 23:57 | κ₀ = 2πIm τ₀/|τ₀|² = 6.2732 (τ₀ = W(∞) = −0.056282+0.998415i, p01); 1/(κ₀|s₁|) = 14.63 für |s₁| = 0.0109 ✓ (C: 14.65). |
| C13.4 Satz (einteilig) | verifiziert|(F) | 23:57 | sauberstes Ziel für (H₀). |
| C13.5 | verifiziert | 23:57 | trivial. |
| C-Trans | nicht geprüft | | Nebenaussage. |

### D
| ID | Status | Zeit | Bemerkung |
|---|---|---|---|
| D1, D2, D3 | verifiziert | 23:57 | |
| D4 | Skizze|(G2) | | Nebenweg. |
| D5 | verifiziert | 23:57 | von Hand (Zerlegung P_N h + h_hi). |
| D6 | verifiziert (polynomiales x₀) | 23:57 | Partialbrüche + ℓ_i(P_N[c/(p−x)]) = (c/p)V_i(p) ✓. Für A_f (nicht-polynomiales f) nicht formuliert (Fernterme als Konturintegral nötig). |
| D7 | verifiziert (unter „alle Urbilder ausserhalb D̄“) | 23:57 | genau diese Voraussetzung fehlt auf dem A-Bogen. |
| D8 | verifiziert | 23:52 | p03 T6: |J_N − M_N(I+K_N)| = 4e-68; σ_{r+1}(XY) ≤ σ_{r+1}(X)‖Y‖, ‖P_N M_{1/m}P_N‖ ≤ sup_D|1/m| ✓. |
| D-(b) Reduktion | **Lücke** | 23:46 | siehe §3. |

### E
| ID | Status | Zeit | Bemerkung |
|---|---|---|---|
| E-Ers, C-Abel | verifiziert | 23:57 | Lδ = (V(α+½)+V(α))/α'(g) nachgerechnet. |
| E-P | verifiziert | 23:57 | Koeffizient von a_k: a₁^k + a₁. |
| T | verifiziert | 23:57 | (4): N^{−σ} setzt einen nicht oszillierenden Anteil von c_n n!h_n voraus — vorhanden (Pole ±2πi). |
| K | verifiziert | 23:57 | Kompaktheit von C auf W_r, ‖(I−T_N)K‖ → 0 ✓. |
| K′ | **verifiziert (rechnergestützt)** | 23:56 | satz_k_beweis.py neu gelaufen (p05_satzK_lauf.log), alle Zusicherungen erfüllt; Logik gelesen. Formaler Schönheitsfehler: Schwanz Σ_{k≥400}|c_k| ≤ w₁.₅·1.5^{−400}·**3** (Faktor 3 fehlt; Grösse 1e-70, folgenlos). |
| M1, M4 | verifiziert | 23:57 | Residuenrechnung r = (1+τ²)/τ² ✓. |
| M2 | verifiziert | 23:58 | Schwelle (1−τ²)² < τ²(1+τ²) ⇔ c > π/3 nachgerechnet (p06). |
| M3 | verifiziert | 23:58 | p06: Green-Funktion direkt, e^{−g(ζ_p)} = Formel auf 10 Stellen (c = 0.4 … π/3−0.01), Keilwinkel = π − c/2. |
| M-lin, M-voll | verifiziert|(H*)_M bzw. (H*_sek) | 23:57 | Brouwer-Argument ✓. |
| E-R3 | **verifiziert für N ≤ 4** | 00:59 | p21 (sympy): Nenner (1−τ²)^{2N+2}, Faktor (1+τ²)^N, Q_N(τ²) irreduzibel vom Grad 4/7/11/16 = (N+1)(N+2)/2+1, Q_N(0) = 2^{N+1}. |

### README
| ID | Status | Bemerkung |
|---|---|---|
| R1, R2 | nicht geprüft | Nebenweg (nur mit (G2)). |
| R3 (Krawczyk N ≤ 120) | verifiziert wie formuliert; **Lücke für (A)** | e51 gelesen: K(X) ⊂ int X komponentenweise, Box x̃ ± 2^{−250}ρ^{−k} — zertifiziert eine isolierte Lösung je N nahe der numerischen Näherung. NICHT zertifiziert: dass der Zweig (Fortsetzung von p₁) diese Lösung ist — die Übergänge N → N+1 fehlen (wie in K′ Teil 3 nötig). Ds Messung: Kugel-Kantorovich für den Übergang versagt ab N ≈ 80–110 (h ~ 1.2^N); Kasten-Krawczyk mit Radius ~ Schrittweite wäre der Weg. (A) braucht das. |

## 5. Phase 2 — neue Aussagen der Stränge W, O, M, X (laufend)

| Eintrag | Aussage | Status | Zeit | Prüfung |
|---|---|---|---|---|
| X 23:53 | Cs „Halbmond-Divergenz ab N ≈ 600–900“ ist Rundungsartefakt | **verifiziert (numerisch)** | 00:05 | p08: x = 1.8: Fehler 7.7e-17, 7.8e-18 (N = 400, 500), ab 600 unter dem Rauschpegel (4.5e-12 < ν = 1.7e-10; 0.06 < 2.2 bei 700; 2.65e29 < 3.2e31 bei 1000); x = 1.4+0.9i fallend bis 700 (8.3e-22), ab 800 im Rauschen. ⇒ ARCHITEKTUR §6.4 und C 23:05(5) „divergiert im Inneren“ sind nicht belegt; nur der Scheitel a₀+R₀ (b42) bleibt. |
| X 23:53 §2.2 | Π_N f_K = T^{a₀}_N f_K − T^{a₀}_N[f_K'·(e_N∘h)] | **verifiziert** | 00:05 | aus J_∞e_N = O(x^{N+1}) bei x = h(y): e_N(y) + f'(y)e_N(h(y)) = O((y−a₀)^{N+1}); T^{a₀}_N Π_N f = Π_N f. |
| X 23:46 §1.1 | Julia–Carathéodory: kein Gebiet mit Tangentialkegel in L ist f- oder h-invariant | **verifiziert** (elementar) | 00:05 | f wirkt bei L als Drehstreckung um arg μ = 38.31°; ein Sektor geht nicht in sich über. Literaturstellen nicht geprüft. |
| X §1.5 (00:03) | Möbius: T_N[P∘f] = L_N((cx+d)^{−N})π_N(g)P; π_N(g)e_j = κ₂^Nμ^j e_j, κ₁/κ₂ = f'(L₁) | **verifiziert** | 00:12 | (a−cL)(x−L) = ax+b−L(cx+d) aus der Fixpunktgleichung; a−cL₁ = cL₂+d (Wurzelsumme); (cL₁+d)(cL₂+d) = ad−bc ⇒ f'(L₁) = κ₁/κ₂. |
| X §2.4 (00:03) | Selektion nur am Rand: f_c = g_c⁻¹∘f_K∘g_c exakte reelle Lösungen ≠ f_K, holomorph um Ū ⋐ D, → f_K | **verifiziert|(K)** | 00:12 | g_c vertauscht mit exp (ψ_c 1-periodisch), f_c ≠ f_K da ψ_c ½-antiperiodisch; Holomorphie um Ū für kleines c (A(Ū) kompakt). |
| O1 (00:14) | J_N = U_N(I+K'_N)S, (I+K'_N)⁻¹ = S(I−Y_N)S⁻¹, Mischformel Π_N g = Y_N T⁽⁰⁾_N g + (I−Y_N)T^{(a₀)}_N g, Pol-Form mit D_ζ | **verifiziert** | 00:25 | Algebra nachgerechnet: P_N[g∘f] = U_N·(a₀-Jet von g), U_N⁻¹ = a₀-Jet von (·)∘h, m∘h = f'; I − Y_N = J_N⁻¹U_N S. |
| O2 (00:14) | Herglotz-Kriterium: Re b ≥ 0 auf D_H ⇔ K akkretiv auf Polynomen ⇒ ‖(I+K_N)⁻¹‖ ≤ 1 | **verifiziert** (Struktur); O2b-Zahlen (Re b ≥ 0.607, 1.204) nicht nachgerechnet | 00:25 | Szegő-Kern-Rechnung 𝔎 = (b(y)+conj b(z))/Δ_H ✓, Herglotz-Kern positiv, Schur-Produkt, Taylor-Matrix eines positiven Kerns psd ✓; Korollar Kodim. r ⇒ σ_{r+1} ≤ 1/(1−ε) per Minimax ✓. |
| O3 (00:14) | lokales Mellin-Modell, keine mellin-diagonale akkretive Norm bei |μ| ≠ 1 | **verifiziert im Modell** | 00:25 | K_μ ist auf H²(H_φ) nur über Mellin-Fortsetzung definiert (με ∉ H_φ); als Modellaussage korrekt, Übertrag auf K heuristisch. |
| O5 (00:14) | Möbius-exp: sup_{D̄}|Π_N f_M − f_M| ≤ CN^κ·0.535^N unbedingt (K'-Kette) | **verifiziert|O2b** | 00:25 | (1) D8 + O2b; (2) für ζ ∉ K' Urbilder ausserhalb D̄ ⇒ Eingabe polynomial; (3) Bernstein–Walsh mit g_{K'} (kleine Ungenauigkeit: auf {g' = 1/N} ist dist ≥ c/N, nicht nur N^{−1/α'} — folgenlos); (4) M1. ζ_p ∉ K', da h(ζ_p) = 1.8894 ∉ D̄. Erste vollständige Durchführung der K'-Kette. |
| X1 (00:58) | (LP) ⇐ (F_{K'}) + (Univ) + (a″) + (Ẽ); Formel c_N = J_N⁻¹Zᵀ(I+ZZᵀ)⁻¹o(g) | **verifiziert** | 01:00 | Beweisschritte (1)–(6) nachvollzogen (Π_N f − q_N = Π_N t_N, Cauchy über ∂Ω', ŷ = R⁺T_{N+d}J_∞[k]_{>N}, P-II mit N+d ≤ 2N Zeilen, Bernstein–Walsh); Kleinste-Quadrate-Formel per Normalgleichung. (Ẽ) in K' grenzwertig an der Ecke (wie §9.3). |
| X 3.1 (00:38) | Rechteck-Abschnitt R = [J_N; B]: J_N⁻¹ = R⁺_top + R⁺_bot S⁻¹BR⁺_top, S = I − BR⁺_bot invertierbar, σ_{d+1}(J_N⁻¹) ≤ 1/σ_min(R) | **verifiziert** | 00:46 | R⁺R = I ⇒ R⁺_top J_N = I − R⁺_bot B; det S = det(R⁺_top)det J_N (Sylvester); Linksinverse nachgerechnet; Weyl für Rang ≤ d. |
| X 3.2 (00:45) | Π_N g = Π̂g + R⁺_bot S⁻¹ρ̂(g), ρ̂ = Überlauf-Residuum | **verifiziert** | 00:46 | R⁺_top w = Π̂g − R⁺_bot·over, BΠ̂g − over = ρ̂(g) eingesetzt. „Π̂ ohne Spitze“ numerisch (nicht nachgerechnet). |
| O6 (00:32) | gemeinsame Eckformel | **verifiziert für polynomiales f; Lücke für A_f** | 00:40 | (b) verlangt 1/(m(η−f)) meromorph auf D_R, R > ρ — bei Kneser falsch (f_K singulär in L auf |x| = ρ), Rest nur ~k^{−6.45} ⇒ Zusatz (E_sing), ARCHITEKTUR §9.4 P-IV. |
| W 00:40 | P's A-Bogen-Wert war float-Rauschen | **richtig** | 00:37 | Bei N = 400 stimmt mein float-Wert (−0.0152) mit Ws arb-Wert (−0.0153); bei N = 1000 versagt Horner innen. |
| O7 / O7a (00:48) | lokales Eckmodell T_L, Fredholm T: E²(Ω) → H²(D), Index 0 per Homotopie | **O7a verifiziert; O7 bis auf zitiertes Lokalprinzip** | 00:56 | Mellin–Plancherel im Sektor nachgerechnet; Kneser-Pfad (Anfrage O) p20: alle Pfadbedingungen erfüllt (Gitter): min|f_t'| ≥ 0.3456, Selbstabstand ≥ 0.030, (S3) ✓, min|m_t| ≥ 0.7319, min|1+μ_t^{s−1}| ≥ 0.999999; Os P1–P5 (o21) für Kneser (p20b, Gitter t = 0…1 in 0.05): min Re f_t' = 0.3456, min Re m_t = 0.7159, Pol-Abstand ≥ 1.45, Vorzeichen/Monotonie ✓, min X_t = 4.447. |
| M 4.2 (00:43) | adjungierte Zwei-Punkt-Form | **verifiziert** | 00:50 | p19: an f_K, N = 14, 1.3e-68 (eigene Implementierung). |
| M 23:55 | radiales Quadrat f = z + ε(z−e^{iθ})(z−e^{−iθ}): f(e^{iφ}) = s(φ)e^{iφ}, μ = 1 + 2iε sin θ; einzige reelle Familie mit dieser Eigenschaft | **verifiziert** | 00:05 | g = f/z reell auf |z| = 1, holomorph bis auf einfachen Pol in 0 ⇒ Spiegelung: g = a/z + b + az ⇒ f = a + bz + az² (2 Parameter = (ε, θ)). σ-Werte numerisch, nicht nachgerechnet. |

| W 01:06 / 01:29 (§14) | Exponent 1.27 zurückgezogen (W1); maßgebliches Γin-Profil p = 1.5–1.6; (E_L)-Weg trägt (Θ_∞ 0.22/0.10, Kneser 0.215/0.103) | **übernommen; Folgerung nachgerechnet** | 01:40 | p17b mit Ws neuem Profil: Ecken-Grenzwert K' 0.191 (p = 1.5) / 0.155 (p = 1.6), W̄ 0.089/0.080 (u → 1e-10, konvergent) ⇒ K' ist NICHT grenzwertig (p·α' = 1.18 > 1). Ws Fit-Methode (ln X = a + κ ln N + λN, Teilfenster) gelesen, Zahlen nicht nachgerechnet. |
| X 01:08 | kanonischer Eckteil in K': Θ ≈ 0.13, Ecken-Grenzwert −30; eckfreier Teil Θ < 0 | **plausibel, konsistent mit p17b** (nicht nachgerechnet) | 01:45 | Vorzeichen auf der f_K(Γout)-Seite (Abklingen) macht den Ecken-Grenzwert negativ; p17b nimmt dort konservativ λ = 0 und bleibt endlich. |
| X 01:22 | Schwanz ‖T_{>3N}J_∞(f_K)δ‖ ≤ CN^{−1}‖δ‖ (Anfrage an P) | **gilt sogar mit N^{−3.95} (Beweisskizze, unter Bedingung (T3))** | 01:50 | Darboux/Transfer an L, L̄: [x^k]f_K^j ≈ jL^{j−1}c·k^{−2−γ}L^{−k}(1 + O(j/k)) gleichmässig für k ≥ 3N ≥ 3j (Glattfaktor (1+μu/L)^{j−1} auf |u| ~ 1/k beschränkt); (f'∘f)δ ~ δ(L)(x−L)^γ ⇒ k^{−1−γ}. Grosskreis |x| = R (Schlitze ausgespart, R < |L ± 2πi|): der Kreisbeitrag ist ≤ max|f_K|^N R^{−3N} und damit ≪ ρ^{−3N}N^{−6.45}, SOFERN max_{|x|=R}|f_K| < (R/ρ)³ für ein solches R (Bedingung (T3); noch nicht nachgerechnet, Motor nötig); Schlitzbeitrag fällt, da (1+|μ|u)/(1+u)³ < 1. ⇒ unter (T3): ℓ¹_ρ-Schwanz ≤ C N^{1−5.45}‖δ‖_{ℓ¹_ρ} ≤ C N^{−3.95}‖δ‖_{H²}. Bei 2N Zeilen bräuchte man max|f_K| < (R/ρ)² — deutlich schärfer; 3N ist die sichere Wahl. |
| O 01:11 O7g | Kommutator [S, c] kompakt auf Sehnen-Bogen-Kurve | **verifiziert** | 01:50 | Lipschitz-Symbol ⇒ Kern ≤ λ/(πκ) beschränkt ⇒ Hilbert–Schmidt; Dichtheit + Abgeschlossenheit der Kompakten. |
| O 01:11 O7h | Spuren auf nichttangentialen Bögen ≤ C‖F‖_{E²}; T beschränkt | **verifiziert bis auf Pommerenke-Eckasymptotik (Standard)** | 01:50 | Mittelwert + Fubini ⇒ ∫|F|²dA/|w−L|; Exponentenbilanz 1−a+2(a−1)−a = −1 nachgerechnet; |ζ−ζ_L|⁻¹dA ist Carleson (Box an ζ_L: πh; fern: ≤ h). |
| O 01:11 O7i | Normstetigkeit der Familie, Index konstant | **verifiziert bis auf CMM/Calderón/Kato (Standard)** | 01:50 | A = ψ_t − ψ_s Lipschitz mit A(L) = A(L̄) = 0; erster Calderón-Kommutator; U_t-Konstruktion Standard. |
| O 01:11 O7j / o22 | T invertierbar für Kneser (qualitativ); σ_min(T_{3N}J_∞|P_N) > 0.7315 für N = 20…100 | **O7j verifiziert|(G2),(S),P1–P4,D4** (D4 = B1.3 nur Skizze; E²(Ω) erfüllt die Wachstumsschranke β < 3.447); o22 **bedingt zertifiziert** | 01:50 | o22 benutzt Knesers Koeffizienten als Bälle mit 10× Strang-A-Fehlerschätzung — nicht zertifiziert (dafür ARCHITEKTUR §9.8). P1–P4 sind von mir nur auf Gittern geprüft (p20/p20b); zertifiziert erst mit §9.8. |
| M 01:04 §8.1 | keine Gårding-Ungleichung (FIO mit imaginärer Phase) | **übernommen** (Korrektur der eigenen Aussage) | 01:50 | — |
| M §8.3 Lemma M3 | σ_min(T_{N+d}J_∞|P_N) ≤ min_Γin|m| + C/N für alle d | **verifiziert bis auf F = 0.8751 < 1 (Gitter)** | 01:50 | Testpolynom zⁿ(1+e^{−iφ₀}z)^{2n}: ‖x‖² = binom(4n,2n), Laplace ⇒ ‖mx‖/‖x‖ = |m(φ₀)| + O(1/n) (Minimum ⇒ kein Term 1. Ordnung), |x∘f| = (s^{1/2}|1+se^{i(φ−φ₀)}|)^{2n} ≤ (2F)^{2n}, Abschneiden verkleinert ℓ². |
| M §1 Lemma M2 | f injektiv auf D̄, (S) exakt, J_N über Legendre-Funktionen | **verifiziert** | 01:50 | f(z₁) = f(z₂) ⇔ (z₁−z₂)(c₀ + ε(z₁+z₂)) = 0 und |ε(z₁+z₂)| ≤ 2ε < c₀; Laplace-Integral für P_k^n mit (c₀ + 2ε cos φ) = √(c₀²−4ε²)(x₀ + √(x₀²−1)cos φ). |
| O 02:09 Satz O9 | ‖Tδ‖_p ≥ λ‖Bδ‖_p ⇒ c_∞ ≥ c(λ) | **verifiziert** | 02:22 | t ≥ λ‖Bδ‖ (Cauchy–Schwarz aus (FA_λ)); a = ‖mδ‖_Γin, b = ‖δ∘f‖_Γout, a² + b² ≤ t²/λ²; auf Γout mδ = J_∞δ − δ∘f ⇒ ‖mδ‖ ≤ t + b; Maximum in β = m_in²/(m_out² − m_in²) = 0.639 ⇒ c(0.98988) = 0.5682 nachgerechnet. |
| O 02:09 Lemma O10 | Eckmodell erfüllt (FA_λ) genau bis λ_L = inf_σ q | **verifiziert (Formel); Zahl: λ_L = 0.98988321** | 02:22 | q = 1 + Re z(1+r)/(1+r|z|²), z = μ^{s−1}, r = e^{−2πσ}: mpmath, Gitter σ ∈ [−60, 60] + findroot: Minimum 0.98988321 bei σ = −5.0384, sup 1.54619. O schreibt 0.9898822, M 0.9898832 — Abweichung in der 7. Stelle (vermutlich Tippfehler bei O). |
| P 03:57 (eigen) | Satz P-V: sexp_K = T̃ ± 1.3e-168 (|Im| ≤ 0.45), = Ũ ± 1.3e-168 (Im ≥ 0.35) | **rigoros** (arb) bis auf zitierte Standardsätze | 03:57 | p27: η_B ≤ 2.7e-173, η_S ≤ 2.1e-173, |e| ≤ 6.5e-173, |d| ≤ 3.7e-173; p28: |μ₀| ≤ 1.3e-171, |φ−id| ≤ 1.3e-169, Prüfungen A/B auf 2958 Scheiben + Schwanz; p29: ε ≤ 1.27e-168. Zitiert: Ahlfors–Bers, Stoilow, ‖S‖₄ ≤ 12, Paulsen–Cowgill Prop. 2. |
| P 03:57 (eigen) | f_K-Koeffizienten a_0…a_640 | **rigoros**: Fehler·ρⁿ ≤ 2.9e-127; |a_n| ≤ 2.128·1.3^{−n} | 03:57 | p29: 4096 Krawczyk-Auswertungen auf |x| = 1.2 + DFT + Alias (M auf |x| = 1.30 über 1950 Bogen-Krawczyks); Gegenprobe A-Datei 1.1e-43. |
| P 02:19 (eigen) | kohärentes Zweig-Zertifikat N > 45 (p24) | **rigoros** (arb) | — | siehe ARCHITEKTUR §9.7 (3b). |
| P (eigen) | zertifizierte Kneser-Werte | **nicht vorhanden** (gp-tetration: nur Vorzertifikate) | — | Weg Satz P-V, ARCHITEKTUR §9.8; Baustein K1 (Poincaré-Funktion an L) zertifiziert: p25, |P − T₈₀| ≤ 3.3e-51·|w|^81 auf |w| ≤ 0.5. |
| W §22.2 | Satz W1'': (H*) auf γ* ⇒ |T_i| ≤ C₁N^{κ+2}e^{−(N+2)g_{K_T}} | **bestätigt** (Beweis korrekt; κ+1 genügt) | 04:50 | Maximumprinzip auf Ext(γ*) ∪ {∞}, Cauchy auf |x| = 2 (braucht ρ(1+s) < 2), D6-Terme, Bernstein–Walsh. Nur für polynomiales f (Kneser: Bem. 3). |
| W §22.1 | Lemma W3: (H-V3) ⇔ |T| ≤ CN^κ auf dem A-Bogen; Λ = 0 auf Γout ⇔ |T| ≤ CN^κ auf f(Γout) | **bestätigt** | 04:50 | D6 nach dem h-Term aufgelöst, |h f'(h)| beidseitig beschränkt. Konturdarstellung: Faktor 1/ζ (k_ζ = ζ⁻¹/(1 − x/ζ)) nachrechnen. |
| W §22.8 | Satz W5: (A1) + (A2) ⇒ limsup max_E[(1/N)ln|Y_N| − G] ≤ 0 | **bestätigt** | 04:50 | Abstiegsprinzip (Saff–Totik I.6.8), Jensen auf |w| = ½; Identität gilt, weil ν = harmonisches Mass von Ω_y in 0. „Umkehrung (BSS)“ braucht deren Zusatzbedingung — nur Bemerkung, kein ⇔. Liefert e^{o(N)}, an der Ecke zu schwach (W §26.3). |
| W §23.2(iii) | Innenscheibe per Maximumprinzip aus (KB) | **bestätigt** (trivial) | 04:50 | Y_i, T_i analytisch in D(L, r_N) (Pole nur 0, a₀). Die Last liegt ganz in (KB) (O). |
| P 04:50 (eigen) | p32a: sexp_K = P(e^{L(z+θ_K)}) auf Im z ≥ 0.36, |θ_K − θ̃| ≤ 4.4e-167; (K3) f_K(L) = L, f_K'(L) = √L | **rigoros** (arb) | 04:50 | 1742 Quadrate (Taylor-Modell von Φ, Krawczyk r = 1e-120), Im ≥ 3 über Rouché + hebbare Singularität; ARCHITEKTUR §9.9. |
| P 04:50 (eigen) | (A): ‖p_N − Π_N f_K‖_{ℓ¹_ρ} ≤ 2.2e-21, ‖p_N − T₆₄₀f_K‖ ≤ 1.5e-10 (N = 140…193) | **rigoros** bis auf den Schwanz S₆₄₀ (= (K2)) | 04:50 | p31_A.json; Boxen p30 (Kette auf aicgolling reproduziert, Radien identisch mit p24). |
| P 04:50 (eigen) | (K2) |b_k| ≤ C₆k^{−6} (numerisch C₆ ≈ 124) | **offen** (Teilintegrale im T̃-Gebiet rigoros) | 04:50 | p32: Strukturgebiet braucht zentrierte Auswertung (§9.9 Wiederaufnahme). |
| offen | O9/O10 (Bestätigung), X4a/b, X5, X6, X7 (c_in = 0.6741, c_out = 0.3456 — vgl. p12: min|f_K'| auf ∂D = 0.3456), X8, Kontur d = 0.3 + Sprungschranke bis u ≈ 0.46, (H-i)-Text, M: mk_kneser_arb, mk_ecke_kneser, m22–m26, Ketten-Check X bis 400 | **offen** | 04:50 | beim Anhalten nicht mehr geprüft |

### Numerische Stützen (Stichproben)
| Grösse | Status | Bemerkung |
|---|---|---|
| (H-ii),(H-iii) am λ-Pfad (nicht-polynomiales P) | **bestätigt N = 60** (00:10) | p09: P(λ) = (1−λ)f + Q(λ) exakt gelöst (Residuum 1e-81), λ = 0, 0.5, 1: dist·max|Π^P k_ζ| unverändert (Ecke 0.2075/0.6152/0.9816, links 0.2975/0.7534/0.9955); (H-iii)-Quotient sup|J(P)⁻¹xᴺ|/sup|(T𝒥_P)⁻¹xᴺ| = 1.0000. Bisher war (H) nur an Polynomen geprüft. |
| **Numerik-Warnung** | M_P = (f∘P − f∘f)/(P − f) per Reihendivision ist INSTABIL | **Fehlerquelle** (00:08) | p09b: P − f hat eine Nullstelle sehr nahe 0 ((P−f)(0) = −5.6e-16 bei sup 2.8e-10, N = 60) ⇒ Division liefert Unsinn (Fehler 9e18, Kopfantwort 2.1e5 statt 4.23e6). Richtig: M_P = Σ_j (f^{(j)}∘f/j!)(P−f)^{j−1} (so bei E r2_07, B b46 ✓). Meine p03-Prüfungen T1/T2 sind davon unberührt (dort geht nur M_P·(P−f) = f∘P − f∘f ein, algebraisch exakt). |
| D8-Konstanten | max_{∂D}|1/m| | **bestätigt** (00:09) | p12: max_Γin |f'e^{−x}| = 1.3663 (bei 180°), max_Γout = 0.8516 (bei 76.5°), min_{|x|=ρ}|f_K'| = 0.3456. |
| (H₀)-Stütze b25 (dist·max|Π_N k_ζ|) | **bestätigt** (00:00) | p07, eigene Matrix: N = 60 Ecke 0.2075/0.6154/0.9826, links 0.2975/0.7534/0.9955; N = 100 Ecke 0.1956/0.5927/0.9774, links 0.2961/0.7526/0.9958 — = B b25 auf 3 Stellen. |
