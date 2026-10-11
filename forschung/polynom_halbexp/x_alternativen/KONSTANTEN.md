# KONSTANTEN — explizite Konstanten und N₀ für den Gesamtbeweis (Strang X)

Stand 11.10.2026, 03:55; Meilenstein 6 ergänzt 04:10 (Zeiten aus `date`). Auftrag der Leitung (Runde 5): für jeden Baustein der Kette die Konstanten, ihre Herkunft,
wie sie in R_N und N₀ eingehen, welches N₀ heute erreichbar wäre, und wo N₀ ≤ 193 knapp wird (P hat den Zweig rigoros bis N = 193).
Kennzeichnung: **[bew.]** bewiesen · **[zert.]** rigoros gerechnet (arb/Intervall) · **[num.]** Gleitkomma/Gitter, nicht zertifiziert ·
**[Skizze]** · **[offen]**. Rechnungen: x21 (Rand- und Bernstein-Konstanten, Schwanz), x22 (Eingangsnormen), x23/x23b (N₀-Bilanz gemessen/beweisbar), x24 (Kontur),
x26 (X7-Minima, A_θ, Schlitzlängen).

---

## 000. Schlussstand (11.10. 04:48, Team angehalten)

| N₀ | Annahmen | Wert |
|---|---|---|
| bewiesen (heute) | Lemma X7 (C_R = 12.7), **Lemma X8′ (A_θ = 13.4)**, Kontur d = 0.3 bzw. Hybrid, Ecke C_E = 1, κ_E = 2 | **≈ 318** (mit X8: 387) |
| repariert | wie oben, A_θ = 9.2 (Zeilen bis 3N; braucht die Sprungschranke bis u ≈ 0.46 von P) | ≈ 262 |
| Optimum A_θ | A*(N) ≈ 0.63 (x27, N = 50–193) | ≈ 93 |
| gemessen | alle Konstanten gemessen | ≈ 47 |

Eck-Konstanten sind die offene Stelle (W); Zertifikat bis 400 (zert/, vorbereitet, nicht gelaufen) würde N₀ ≤ 400 abdecken. Lemma X7/X8′ bewiesen
bis auf [num.]-Randwerte (c_in, c_out, min|f'|, max|m/f'|, M2); Kellogg–Warschawski: Original nicht gelesen, für X4a nicht mehr nötig.
(Fam-Eck): Bogenmodell (Glättung über Bögen fester Länge, ε ≈ 0.06) sagt alle Eckantworten auf 1–2 R voraus; Beweis offen (mit W).

## 00. Stand Meilenstein 6 (11.10. 04:10) — BEWEISBARE Konstanten (Details BERICHT §3.9)

- **Lemma X7 bewiesen** mit expliziter Konstante: |w_η|_{ℓ²_ρ} ≤ 6.3·√M (statt gemessen 1.17) ⇒ C_R = 12.7, κ_R = 1. Offen nur: Zertifikat der
  Gitter-Minima c_in = 0.6741, c_out = 0.3456 und max|m| = 2.3376 [P].
- **Lemma X8 bewiesen** (elementar, Koeffizienten von J_∞δ bis Zeile 5N): **A_θ = 19.7** bei c₀ = 0.54 (17.3 bei 0.7319, 22.2 bei 0.4249; mit Zeilen bis
  3N 15.5). Die gemessene Bernstein-Konstante gab 4.05. Lemma X6 und damit **Kellogg–Warschawski werden nicht mehr gebraucht.**
- **Kontur d = 0.3:** Lage nachprüfbar (x24/x26: |Γ| = 12.871, δ = 0.5095, max|z| = 2.357, M_d = 4.251, Konstante 17.1), Schlitz bis u₀(0.3) = 0.451.
  Fehlt: Sprungschranke |Δf(L(1+u))| bis u ≈ 0.46 (heute bis 0.2) — dieselbe Lücke wie Lemma X5 mit R = 2. Ohne sie: Hybrid-Kontur, gleiches N₀.
- **N₀ mit beweisbaren Konstanten:** 387 (Ecke C_E = 1, κ_E = 2), selbst mit perfekter Ecke ≈ 214; N₀ ∝ A_θ^{0.5} (A_θ = 2/4/8/13/20/30/50 ⇒
  N₀ = 121/172/244/314/387/482/629). **N₀ ≤ 193 ist mit beweisbaren Konstanten nicht erreichbar** — die Reserve kommt vom Zertifikat:
- **Zweig-Zertifikat verlängern** (Hochrechnung aus p24): bis N = 400 ≈ 12 h sequentiell bzw. ≈ ½ h auf 32 Kernen (NC ≈ 12 Knoten, 1120 Bit,
  ≈ 0.2 GB/Prozess); bis 1000 ≈ 32 Tage bzw. 1–1.5 Tage (NC ≈ 26, 2320 Bit, ≈ 2 GB/Prozess). Bis 400 deckt das Szenario C_E ≤ 1, κ_E ≤ 2 ab.
- Ältere Teile unten (§0–§3) rechnen noch mit den gemessenen 1.17 und 4.05; die beweisbare Bilanz steht in BERICHT §3.9.4.

## 0. Ergebnis in Kürze (Stand 03:55, mit gemessenen Konstanten)

1. **Mit den heute bewiesenen Exponenten wäre N₀ ≈ 3500**, also weit über 193. Der Engpass ist die Familien-Bedingung (Satz X4a:
   θ_N = A_θ·N·R_N ≤ ½, mit Bernstein-Faktor N), kombiniert mit zwei groben Stellen der Buchhaltung:
   (a) die Konturwahl d = 0.1 in Satz P-I (Konturteil ≈ 38·N^κ·e^{−0.069N} dominiert R_N bis N ≈ 400) und
   (b) die Eingangsschranke von Satz P-II (|w_η| ≤ e·N^{3/2}/(c_fρ), dadurch κ_Rest = 2).
2. **Zwei Reparaturen bringen N₀ auf ≈ 150 (< 193):** (a) Kontur d = 0.3 (x24: Konturkonstante 17.1, Γ_{0.3} bleibt in |x| ≤ 2.36, Abstand
   1.23 zur nächsten Singularität −0.362 ± πi); (b) **Lemma X7** (unten, Skizze): |w_η|_{ℓ²_ρ} ≤ C_w√(M+1) mit C_w ≈ 1.17 (gemessen, M ≤ 800)
   ⇒ κ_Rest = 1. Dann hängt N₀ nur noch an den **Eck-Konstanten** C_E, κ_E der Pol-Antwort (Ws Teil): N₀ ≤ 193 verlangt etwa
   κ_E ≤ 2 bei C_E ≤ 1, oder κ_E ≤ 2.5 bei C_E ≤ 0.1 (Tabelle §3). Mit κ_E = 2.5, C_E = 1 (Satz W1 mit „+2“) wird es 900 — **das ist die knappe Stelle.**
3. **Mit den gemessenen Werten** (R̂(L) ≈ 0.33N^{1.1}, Rc(L) ≈ 0.01N^{1.6}, ε = 0.87, c₀ = 0.73) wäre N₀ ≈ 47 — die Lücke ist Buchhaltung, nicht Sache.
4. **Startbedingung (A):** bei N₀ ≈ 150 ist R_{N₀} ≈ 1e-3 (Szenario S3) — zertifizierte f_K-Werte auf D̄ auf ≈ 1e-5 genügen (Ps Satz P-V zielt auf ≪ 1e-20).
5. **Wer was explizit liefern muss** (§4): W — C_E, κ_E, ε für (H-I), (H-V1), (H-V3) bzw. (Ẽ) für alle N ≥ N₀ (und den „+2“ im Satz W1 auf „+1“ drücken);
   O/M — c₀ = c_∞ zertifiziert für Kneser (heute 0.540 num.) und die Kanten-Konstante für (a″)_{αN}; P — Kontur d = 0.3 prüfen, Lemma X7 prüfen,
   C_J/c₁/b35b′/Univalenz zertifizieren, f_K-Zertifikat (P-V).

---

## 1. Die Kette und wo die Konstanten eingehen

(K) ⇒ (F) [L8.2–8.4] · (H₀) für J_∞ [Rest in K' + Eckteil in W̄ bzw. X1/X3] · (Fam) ⇒ (H) für alle P ∈ 𝒮_N ∩ B_N [X4a + (Fam-Eck)] ·
H1 [a-priori + Überlagerung] · H2 [σ-Homotopie, (H-iii)] · (A) [Zweig bis N₀ + zertifizierte f_K].

**Ballradius.** H1 braucht sup_{D̄}|Π^P_N f − f| ≤ R_N/K₀ (K₀ = 6 aus H2). Mit Satz P-I (Cauchy über Schlitze und Kontur Γ_d):

  R_N/K₀ = 2·S_{K'}·C_R^P·N^{κ_R − 4.2877} + 2·S_{W̄}·ε_E^{−4.8094}·C_E·N^{κ_E − 4.8094} + (M_d|Γ_d|/(2πδ_d))·C_E·N^{κ_E}·e^{−ε_E N d},

Rest-Teil in K' (volle Green-Rate, Schwelle 5.447·α' = 4.2877), Eckteil in W̄ (Schwelle 4.8094), C_R^P = 2C_R (für 𝒥_P gilt (a″) bei θ ≤ ½ mit c₀/2).

**Bedingungen an N₀** (alle für N ≥ N₀):
- (B1) **Familie, eckfrei (Satz X4a):** θ_N := A_θ·N·R_N ≤ ½, A_θ = K_B(1 + max|m|/c₀) + ½M2/(c₀N) ≈ 4.05 (c₀ = 0.54).
- (B2) Lemma X6 anwendbar: R_N ≤ c/(2N) — schwächer als (B1).
- (B3) Schwanz (nur Variante d = 4N): τ_N ≤ 0.01 — gilt ab N ≈ 50 (x21).
- (B4) (Fam-Eck): C_F·N^a·R_N ≤ ½ — mit den gemessenen C_F schwächer als (B1).
- (B5) Eckaussage (Θ^c ≤ 1 − ε, R^c(L) ≤ C_E N^{κ_E}) und (a″) mit c₀ — für ALLE N ≥ N₀ (das ist (H₀) selbst; heute nur numerisch).
- (B6) (A): Zweig bis N₀ zertifiziert (erledigt bis 193) und f_K auf D̄ zertifiziert auf ≪ R_{N₀}.

---

## 2. Tabelle der Konstanten

| Baustein | Konstante | Wert | Status / Quelle | Rolle |
|---|---|---|---|---|
| (K) ⇒ (F) | min|f_K(e^z)| (b35b′) | 0.7345 | [num.] Motor, P p04 | (F) |
| | Radius der lokalen Formel bei L | > 0.02 | [num.] | (F) |
| Sprung an L | C_J (u ≤ 0.2), σ | ≤ 0.95, 5.447 | [num.] A k25, P | S_{K'}, S_{W̄}, X5 |
| Univalenz | min|f(x)−f(y)|/|x−y| auf ∂D = min|f'| | 0.3456 | [num.] Gitter, P p16/p20; x21 | P-II, X7, X6 |
| Rand ∂D | max|f'|, max|f| = f(ρ) | 1.6912, 2.2327 | [num.] x21 | X6, Tangente |
| | max|m|, min|m| (m = eˣ/f') | 2.3376, 0.7319 | [num.] x21 (min auch P) | A_θ, O9 |
| | M2 = sup_{∂D}|f''∘f| | 0.8183 | [num.] x21 | |M_P − m| ≤ ½M2·R |
| Geometrie W̄ | α, c_g, c_d, κ₀ | 0.882964, 0.829, 0.933, 6.2732 | [num.] P p01 | S_{W̄} |
| Geometrie K' | α', c_g', c_d', κ₀' | 0.787172, 0.754–0.839, 0.787, 5.8170 | [num.] P p13; c_d' x (heute) | S_{K'} |
| P-I Schlitz | S_{W̄} = 2(C_J/(c_dπ))αΓ(5.447α)c_g^{−5.447α} | 12.7 | [bew.] Formel, Zahlen [num.] (P) | R_N Eckteil |
| | S_{K'} (c_g' = 0.754, c_d' = 0.787) | 17.7 | wie oben (X) | R_N Rest |
| P-I Kontur | M_d|Γ_d|/(2πδ_d) für d = 0.1/0.2/0.3/0.5 | 37.9/22.0/17.1/14.2 | [num.] x24 (W̄-Niveaulinien; |Γ_d| = 10.5–15.7, δ_d = 0.15–0.95, M_d = 3.4–5.4) | R_N |
| | Holomorphie von f_K innerhalb Γ_d | bis |x| ≤ 2.8 nötig | [Skizze] slog-Argument (X §3.5c), nächste Sing. |x| = 3.16 | P-I |
| (a″)_∞ | c₀ = c_∞ | ≥ 0.540 (O11), wahr 0.7319 | [num.] O §30; radial streng 0.4249 (M6) [zert. im Modell] | A_θ, C_R |
| | O9: c(λ) | c(0.99) = 0.568 | [bew.] Formel, P verifiziert | c₀ |
| Lemma X5 (d = 4N) | R, β, τ_N | 1.649, ≤ 0.798, τ_N ≤ 1.6e-4/2.9e-9/9e-12/3e-12 (N = 50/100/150/193) | [bew.|(K)], c₁ [Skizze]; x21 | (B3) |
| (a″)_{αN} (kleines d) | d₀(N), c | 0.02N − 2; c = 0.7319 ab d₀ | [num.] x20 (N ≤ 800) | Variante kleines d |
| Lemma X6 / X4a | K_B = max‖δ'∘f‖/(N‖δ∘f‖) | 0.733/0.741/0.751/0.756/0.759 (N = 25/50/100/150/193) | [num.] x21 (Arnoldi auf f(∂D)) | A_θ |
| | A_θ = K_B(1 + max|m|/c₀) | 4.05 (c₀ = 0.54), 4.94 (c₀ = 0.4249) | aus obigen | (B1) |
| P-II (Rest-Eingang) | |w_η|_{ℓ²_ρ} | ≤ 5.8·N^{3/2} [bew.]; gemessen ≤ 1.17√(M+1) | P-II; x22 (Maximum bei η = L) | C_R, κ_R |
| Lemma X7 (neu) | C_X7 = e^{1/2}(1/c_K + max|m|)/ρ, c_in, c_out | **6.28–6.30**; c_in = 0.6741, c_out = 0.3456 | [bew.] bis auf Gitter-Minima [num.] (x26); gemessen 1.17 (x22) | C_R = 12.7, κ_R = 1 |
| Lemma X8 (neu) | A_θ = (1/min|f'|)[5/ρ + (max|m|/ρ + max|m'|/N)/c₀] | **19.7** (c₀ = 0.54); max|m'| = 1.3839 | [bew.] bis auf [num.]-Randwerte (x26) | (B1), ersetzt K_B·(…) |
| Lemma X8′ (Meilenstein 7) | A_θ = (5/ρ)/min|f'| + max|m/f'|/(ρc₀) + M2/(c₀N) | **13.4** (c₀ = 0.54); 9.2 mit Zeilen bis 3N; max|m/f'| = 2.1176 | [bew.] bis auf [num.]-Randwerte (x27) | (B1) |
| Optimum (Vergleich) | A*(N) = max ‖δ'∘f‖/(N‖J_∞δ‖) | 0.610/0.621/0.625/0.627 (N = 50/100/150/193) | [num.] x27 (präzisionsgeprüft) | Ziel für A_θ |
| Rest-Pol-Antwort | C_R, κ_R | 10.7·N² [bew.]; 2.9·N [mit X7]; gemessen R̂(L) ≈ 0.33N^{1.1} | P-II; X7; x12/x13 | R_N |
| Eck-Pol-Antwort | C_E, κ_E, ε_E | gemessen Rc(L) = 2.13/3.43/7.53 (N = 200/400/600, kleines d), ε ≈ 0.87 (Θ^c ≈ 0.13); W1: ε_Γ = 0.694 | [num.]; [offen] für alle N | R_N, (B5) |
| (Fam-Eck) | δσ₁/σ₁ | ≤ 0.34R (kohärent); glatte Störungen: N·(0.34·δr + Re((0.042 − 0.035i)δμ)) + O(R), Radialkoeffizient a(N) ≈ 0.34N − 46 | [num.] x18/x19 (N ≤ 600; Modell auf 2 %, BERICHT §3.8.5) | (B4), (H-iii) |
| H1/H2 | K₀ | 6 | [bew.] (P) | R_N |
| (A) | Zweig-Zertifikat | N ≤ 193 | [zert.] P p18/p24 | (B6) |
| | f_K auf D̄ zertifiziert | fehlt (P-V geplant, Ziel ≪ 1e-20) | [offen] | (B6) |

**Lemma X7 (Skizze).** Für η ∈ ∂K' sei ξ := h(η) (für η ∈ f_K(Γout) ist ξ ∈ Γout ⊂ ∂D; für η ∈ Γin nahe L liegt ξ ausserhalb D̄ nahe L; fern von
f_K(D̄) ist nichts zu tun). Dann ist 1/(η − f(x)) = 1/((ξ − x)·Q_ξ(x)) mit Q_ξ(x) = (f(ξ) − f(x))/(ξ − x) ≠ 0 auf D̄ (Univalenz; Q_ξ(ξ) = f'(ξ) ≠ 0),
und wegen |ξ| ≥ ρ ist |[1/(η − f)]_k|ρ^k ≤ ρ^{−1}·Σ_j |[1/Q_ξ]_j|ρ^j = ‖1/Q_ξ‖_{A_ρ}/ρ, also |T_M[1/(η − f)]|_{ℓ²_ρ} ≤ √(M+1)·‖1/Q_ξ‖_{A_ρ}/ρ.
‖1/Q_ξ‖_{A_ρ} < ∞ gleichmässig in ξ, weil 1/Q_ξ auf ∂D gleichmässig Hölder-stetig mit Exponent > ½ ist (f ∈ C^{1,α} auf ∂D, auch an L: C^{5.4};
Satz von Bernstein über absolut konvergente Fourier-Reihen). Ebenso |T_M[m k_η]| ≤ ‖m‖_{A_ρ}√(M+1)/ρ. Gemessen (x22): max_η |w_η|/√(M+1) =
1.173/1.169/1.167/1.166 (M = 100/200/400/800), Maximum bei η = L; dort ist |T_M[1/(L − f)]| = 6.25/8.81/12.43/17.57 ≈ √M/(ρ|μ|) = 0.621√M — genau der Pol-Anteil 1/((L − x)μ). ⇒ in Satz P-II N^{3/2} → N^{1/2}, also κ_R = 1 statt 2.
Offen für die Strenge: die gleichmässige Hölder-Schranke von 1/Q_ξ für ξ → L (Differenzenquotient über den Schlitz hinweg) — P bitte prüfen.

---

## 3. N₀-Bilanz (x23)

θ_N = A_θ·N·R_N mit R_N aus §1; N₀ = kleinstes N mit θ_N ≤ ½ (dazu: Konturteil ≤ Schlitzteil).

| Szenario | Annahmen | R_193 | θ_193 | **N₀** |
|---|---|---|---|---|
| S1 heute bewiesen | P-II wie bewiesen (C_R = 10.7, κ_R = 2), d = 0.1, Ecke C_E = 1, κ_E = 2, ε_E = 0.69, c₀ = 0.54 | 14 (Kontur) | 1e4 | **3498** |
| S2 | wie S1, aber Lemma X7 (C_R = 2.9, κ_R = 1) | 14 (Kontur) | 1e4 | **396** |
| S3 | wie S2, Kontur d = 0.3 (Konstante 17.1) | 3.8e-4 | 0.30 | **147** |
| S4 | wie S3, Ecke κ_E = 3 | 6.7e-2 | 52 | 59842 |
| S5 | wie S3, c₀ = 0.4249 | 3.8e-4 | 0.37 | **163** |
| S6 gemessen | R̂ ≈ 0.33N^{1.1}, Rc ≈ 0.01N^{1.6}, ε = 0.87, c₀ = 0.73, d = 0.3 | 7.4e-6 | 4.6e-3 | **47** |

**Wo N₀ ≤ 193 knapp wird — Eck-Konstanten in S3:**

| κ_E \ C_E | 0.1 | 1 | 10 |
|---|---|---|---|
| 1.5 | 78 | 91 | 137 |
| 2.0 | 91 | 147 | 492 |
| 2.2 | 99 | 258 | 1058 |
| 2.5 | 166 | 899 | 5208 |
| 3.0 | 3481 | 59842 | — |

Lesart: Ws Satz W1 liefert für den Eckteil C₁N^{κ+2}e^{−(N+2)g} mit κ aus (H-I)/(H-V1)/(H-V3) (gemessen κ ≈ 0.5) ⇒ κ_E ≈ 2.5; dann
braucht N₀ ≤ 193 eine Konstante C_E ≲ 0.12. Drückt man den Polynomfaktor von „+2“ auf „+1“ (κ_E ≈ 1.5), genügt jedes C_E ≤ 10.
Der eckfreie Teil ist mit Lemma X7 unkritisch; ohne X7 (κ_R = 2) liegt N₀ auch mit d = 0.3 bei ≈ 3500 (S1 mit d = 0.3: 3523). Bei κ_E = 2.5 liegt die Grenze für N₀ ≤ 193 bei C_E ≈ 0.12 (C_E = 0.08/0.10/0.12/0.15 ⇒ N₀ = 144/166/187/219).

**Weitere Stellschrauben (nicht ausgereizt):** (i) K₀ = 6 ist nur für H2 nötig; die Bedingung (B1) gilt eigentlich nur auf 𝒵 (sup|P − f| ≤ R_N/K₀),
falls man (H) nur dort verlangt (vgl. X §3.8.3) — Faktor 6 in θ; (ii) θ ≤ ½ ist willkürlich (θ < 1 genügt mit kleinerem c₀); (iii) d = 0.5
statt 0.3 (Konstante 14.2, braucht Holomorphie bis |x| = 2.8).

---

## 4. Was W, O, M, P noch explizit liefern müssen (an das BRETT)

- **W:** für die Eckmoden-Aussagen (H-I), (H-V1), (H-V3) bzw. (Ẽ)/Θ^c **explizite** C, κ, Λ (bzw. ε) **für alle N ≥ N₀** — heute numerisch bis
  N = 1000. Ziel nach §3: Eck-Polynomfaktor κ_E ≤ 2 mit C_E ≤ 1, oder κ_E ≤ 2.5 mit C_E ≤ 0.1. Konkret: der Faktor (2x₀N + 1) und die
  Bernstein–Walsh-ℓ¹-Umrechnung in Satz W1 kosten zusammen N^2 — lässt sich eine Potenz sparen?
- **O / M:** c₀ = c_∞ für Kneser **zertifiziert** (O11 heute 0.540 num.; Ms 𝒯-Weg; braucht Ps zertifizierte Koeffizienten K1/K2); und für die
  Variante kleines d die Kanten-Konstante (a″)_{αN} ≥ c für α ≥ α* (O: Doppel-Skalenlimes; X misst α* ≈ 0.02, d₀ = 0.02N − 2, und
  Θ^c verträgt α bis ≈ 0.1, siehe BRETT).
- **P:** (i) Kontur d = 0.3: Holomorphie von f_K innerhalb Γ_{0.3} (max|z| = 2.36) und die Konstanten von x24 nachrechnen; (ii) Lemma X7 prüfen
  (gleichmässige Hölder-Schranke von 1/Q_ξ nahe L); (iii) die [num.]-Konstanten zertifizieren, die in R_N eingehen: C_J, c_g, c_d, min|f'|,
  max|m|, M2, K_B (letzteres über eine zertifizierte Arnoldi-/Gram-Rechnung auf f(∂D) — Grösse 0.76·N); (iv) f_K-Zertifikat P-V auf ≈ 1e-5
  reicht für N₀ ≈ 150.
- **X (ich):** Lemma X7 ausarbeiten; (Fam-Eck) mit expliziter Konstante, sobald W die Eckasymptotik beweist; K_B rigoros über Lemma X6 mit
  expliziter K–W-Konstante (konforme Abbildung von f(∂D) numerisch mit Fehlerschranke) — sonst bleibt K_B [num.].

---

## 5. Literatur zu Kellogg–Warschawski (Auftrag 3) — **Original nicht gelesen**; für X4a seit Lemma X8 entbehrlich

Original (Pommerenke, *Boundary Behaviour of Conformal Maps*, Grundlehren 299, Springer 1992, §3.3) ist nur über Login erreichbar — **nicht selbst
gelesen**. Sekundär [S]: arXiv 1901.07844 zitiert „[Pom92, Theorems 3.5, 3.6]“ für Kellogg–Warschawski (Dini-glatt ⇒ φ' stetig fortsetzbar, mit
Stetigkeitsmodul von w'; C^{m,α} ⇒ φ^{(m)} Hölder); arXiv 1902.03753 zitiert Thm 3.5 für den Dini-glatten Fall; arXiv 1803.09458 zitiert Thm 3.6 für C^{m,α}.
Die **Nullstellenfreiheit von φ' auf dem Rand** (für |dz| ≍ dθ in Lemma X6 nötig) habe ich in keiner zugänglichen Quelle wörtlich gefunden [U];
sie ist Teil der üblichen Fassung (Pommerenke Thm 3.5: „φ' has a continuous extension to D̄ which has no zeros“ — aus dem Gedächtnis, ungeprüft).
Für uns genügt C^{1,α} (f(∂D) ist C^{5,0.44}). **Bitte (Leitung/Janis):** falls jemand Zugang zum Buch hat, Thm 3.5/3.6 wörtlich nachtragen.
