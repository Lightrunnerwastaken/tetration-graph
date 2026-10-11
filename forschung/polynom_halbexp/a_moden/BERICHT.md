# Strang A — Moden & Kneser-Daten (Bericht)

Stand: 10.10.2026, 23:03 (Zeiten aus `date`). Alle Skripte in diesem Ordner, Logs in `daten/`.
Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.

## 0. Kurzfassung

1. **Daten (numerisch):** Knesers Taylor-Koeffizienten von f_K bis k = 1000 (≈ 47 Stellen, rel. 1e-29 bei k = 300),
   slog-Koeffizienten, Störrichtungen δ_{1,3,5}, Kneser-Werte für andere Basen; alles in `gemeinsam/`.
2. **Ein-Moden-Struktur (numerisch, unabhängig von Strang B mit dem Motor bestätigt):** p_N − f_K ist im Inneren von D bis auf
   relativ 1e-14 (N = 150) gleich c_N δ_{−1} + c̄_N δ_{+1}; damit p_N(1) − T(½) = −r_N(1) + 2T'(½) Re c_N auf 5 Stellen ab N = 55.
   |c_N| ~ N^{−7.9}, arg c_N dreht mit −2.1…−2.2 pro Einheit log N ⇒ c_N ≈ C·N^{1−2γ}, γ = 2πi/L.
3. **Knesers Singularität an L (Rechnung + numerisch):** f_K = φ(√L Ψ e^{Lε(Ψ^γ)}) mit Schröder-Funktion Ψ an L und Knesers
   θ-Koeffizienten κ_m; Leitterm d₁ t^{1+γ}, d₁ = 2L^{3/2}κ₁; Darboux-Reihe trifft die Taylor-Koeffizienten bis 1e-18 (k = 300).
4. **Amplitudengesetz (numerisch):** freier Fit c_N ~ N^{−7.888−2.107i} (Soll 1 − 2γ = −7.894 − 2.116i);
   c_N = K·N·u_N/(δ_{−1})_N mit Knesers Leitkoeffizient, K ≈ 0.002 − 0.051i. **Blindtest:** daraus vorhergesagter
   Vorzeichenwechsel von p_N(1) − T(½) zwischen N = 216 und 218 — gerechnet (p_N bis 232): Pendelzone 215–218, ab 219 positiv.
5. Die fast-singulären Richtungen der Gewichtung (1.3 | 0.5) sind (nach Strang C) ein Normartefakt; Aufgabe 3 wurde auf
   Wunsch der Leitung durch 2.–4. ersetzt. „Anteile des Rests entlang der Moden“ gibt es koeffizientenweise nicht: c_N entsteht
   aus Auslöschung um 1e27 (N = 150); Beweisrichtung in §4d.
6. **Runde 2 (lineares Problem und grosse N, numerisch):** Π_N f_K exakt aus a_0…a_{N+1} bis N = 1000 (Läufe zu 1500 in §R2.7):
   max_D|Π_N f_K| = 2.232651 für alle N (Satz B-lin mit a = 0); c_N ~ N^{1−2γ} bestätigt, K∞ ≈ 0.0017 − 0.0533i;
   |Π_N f_K(L) − L| ~ N^{−4.85}; Eckschicht zweistufig (Kern ~N^{−0.88}, Anschluss ~1/N). **Kein Kippen bis N = 1000:** echte p_N
   (Newton von Π_N f_K, quadratisch konvergent) mit p_N − Π_N f_K ≈ 1e-3·N·|E_N(L)|², Zweig bis N = 300 lückenlos bestätigt,
   kein exponentieller Anteil (Rate < 1e-4 pro Grad). N = 1500: Gesetze gelten weiter, p_1500 quadratisch konvergent, aber nur
   mit ungerundeten Zwischenständen (Rundungsverstärkung ≳ e^{0.16N}, §R2.7). Bedingung (F) für Strang C geprüft.

## 1. Kneser-Taylor-Koeffizienten (Aufgabe 1) — numerisch

- Verfahren `k1_kneser_koeff.py`: Werte A = slog(z), f = sexp(A + ½) auf |z| = r, M = 4096 (obere Hälfte, Rest durch Spiegelung),
  DFT mit arb. Motor: fatou.gp-Fork direkt über `FatouGP` (eigener Speicher in `daten/kreis_*.json`, der 190-MB-Zwischenspeicher
  wird nicht angefasst), dps 57.
- **Motorgrenzen (gemessen):** komplexe sexp gelingt bis dps 57 (Fehler ≈ 10^{−(dps−9)}, also 4e-48), ab dps 59 „exp overflow“
  bzw. bei dps 80 „log domain error“ (auch frisch initialisiert). slog bei dps 80 ist genau (A(e^z) − A(z) − 1 ≈ 3e-71), hängt
  aber bei einzelnen Punkten nahe L (|z − L| ≈ 0.035) minutenlang; dps 100 springt auf einen anderen Zweig. Eine 70-stellige
  Fassung (Newton mit slog bei dps 80, `k1b_hochgenau.py`) ist deshalb noch nicht fertig.
- Ergebnis: Radien 1.34, 1.30, 1.25, 1.20; Prüfrest max|T(A+1) − e^z| ≤ 4.1e-48. Fehler von a_k ≈ 4e-48·1.34^{−k}:
  relativ 1e-46 (k ≤ 3), 1e-35 (k = 100), 1e-29 (300), 1e-24 (400), 7e-21 (600), 3e-9 (1000). Σ_{k≤1200} a_k = T_e(½) auf 40
  Stellen. Hüllkurve |a_k||L|^k(k+1)^{6.447} = O(1…100) (siehe §4: 2|D₁/Γ(−1−γ)| = 71.57).
- Dateien: `gemeinsam/kneser_halbexp_koeff_300.json`, `…_1000.json` (a, fehler_abs, rel_fehler), `…_400.json` (Liste für D),
  `gemeinsam/kneser_slog_koeff_300.json`.
- Andere Basen (Anfrage E, `k4_andere_basen.py`, dps 50): `gemeinsam/kneser_andere_basen.json`. sexp₂(½) = 1.45878181603642170068…,
  sexp₁₀(½) = 2.47700560634496477508…, sexp_{1.5}(½) = 1.28087727794027296901…, sexp_e(1/3) = 1.40302138219710974153…,
  sexp_e(−2/3) = 0.33862804134425504883…; halbe 2- und 10-Funktion k ≤ 40 (Fehler 3e-41 bzw. 1e-32).

## 2. Störrichtungen δ_n (Aufgabe 2) — numerisch

δ_n(z) = T'(A(z)+½)·e^{2πinA(z)}, T'(A+½) = f_K'/A' spektral auf dem Kreis (`k2_stoerrichtungen.py`). δ_{−n}: konjugierte Koeffizienten.
Radien je n wegen Aliasing (r/|L|)^M·M^{4.447n−2}: n = 1: 1.34/1.30, n = 3: 1.30/1.25, n = 5: 1.25/1.20.

| Prüfung | n = ±1 | n = ±3 | n = ±5 |
|---|---|---|---|
| δ(f(z)) + f'(f(z))δ(z), relativ (6 Punkte) | ≤ 7e-35 | ≤ 9e-27 | ≤ 7e-17 |
| koeffizientenweise (k ≤ 60, Gewicht 0.5^k) | 7.5e-36 | 5.6e-30 | 2.9e-22 |
| Wachstumsexponent an der wachsenden Seite (t = 0.995) | −3.466 (Soll −3.447) | −12.39 (−12.34) | −21.31 (−21.24) |
| Exponent an der anderen Seite | 5.454 (Soll 5.447) | 14.38 (14.34) | 23.30 (23.24) |
| |δ_{n,k}||L|^k ~ k^{4.447n−2} | 7.6e5 (k=100) → 1.15e7 (300) | 2.6e18 → 7.2e23 | 5.2e27 → 2.5e37 |

δ_{+n} wächst bei L̄, δ_{−n} bei L. Datei `gemeinsam/stoerrichtungen_delta.json` (k ≤ 1500, Re/Im, Radienabstand je k).
Bemerkung: δ_{+1} hat bei L denselben Exponenten 1+γ wie Knesers eigene Leitsingularität (§4).

## 3. Modengehalt von p_N (statt Aufgabe 3/4) — numerisch

**Abel-Koordinaten (`k3b_abel_moden.py`).** Liegt E_N = p_N − f_K im Kern der Linearisierung, so ist θ_N(w) = E_N(T(w))/T'(w+½)
1-periodisch und ½-antiperiodisch. Gemessen auf w = x + iy, x ∈ [−1.6, −0.6), y ∈ {0, 0.25, 0.5}, T(w) vom Motor, m = 64:

| N | |c_{−1}| | arg c_{−1} | |c_0|/|c_{−1}| | |c_{±2}|/|c_{−1}| | |c_{−3}|/|c_{−1}| | Höhen-Konsistenz | p_N(1) − T(½) | Modell |
|---|---|---|---|---|---|---|---|---|
| 50 | 7.83e-15 | −1.43 | 4.7e-6 | 4.6e-6 | 4.6e-6 | 1.2e-5 | 3.4661e-15 | 3.4907e-15 |
| 100 | 3.27e-17 | −3.00 | 2.3e-14 | 2.3e-14 | 7.3e-13 | 2.3e-12 | −1.0193e-16 | −1.0193e-16 |
| 150 | 1.34e-18 | 2.36 | 1.4e-19 | 4.4e-18 | 3.1e-14 | 2.5e-18 | −2.9922e-18 | −2.9922e-18 |
| 194 | 1.73e-19 | 1.79 | 5.7e-21 | 5.4e-19 | 4.0e-15 | 3.1e-19 | −1.1953e-19 | −1.1953e-19 |

Modell: p_N(1) − T(½) = −r_N(1) + 2 T'(½) Re c_{−1}, r_N(1) = Σ_{k>N} a_k, T'(½) = δ_{±1}(1) = 1.5751579377843008.
- Bis N ≈ 45 dominiert der Kneser-Schwanz −r_N(1) (Rate 1/|L| = 0.7275 pro Grad), danach der Störmodus; der Vorzeichenwechsel
  bei N ≈ 51 ist der Übergang. Die „0.93/Grad“ sind die lokale Steigung von |c_N| ~ N^{−7.9} (bei N ≈ 100: e^{−7.9/100}).
- Exponent: zwischen N = 50, 100, 194 gilt Re β = 7.90 bzw. 7.91, Im β ≈ 2.2 (c ~ N^{−β}); lokale Exponenten schwanken mit
  Periode 2π/arg L = 4.70 in N (Amplitude ~ 1/N). Passt zu c_N = C N^{1−2γ} (1 − 2γ = −7.894 − 2.116i) und zu Strang B
  (c̄_N/N^{1−2γ̄} = 0.18 − 0.08i ± 5 %): aus meinen Daten |C| = 0.20, arg C passend auf ~0.05 rad.
- **Vorhersage:** arg c_{−1} erreicht π/2 bei N ≈ 215 ⇒ nächster Vorzeichenwechsel von p_N(1) − T(½) bei N ≈ 215 (bestätigt, §4c).
- Koeffizientenbild (`k3a_zerlegung.py`): d_N = p_N − q_N ist für k ≲ N/2 gleich c δ_{−1} + c.c., nahe k = N gibt es eine
  Randschicht: der Spitzenkoeffizient von p_N ist fast 0 (a_N^{(N)}/a_N = −4.7e-4 bei N = 100, 1.3e-7 bei 194), d.h. d_{N,N} ≈ −a_N,
  während c δ_{−1,N} ≈ 50·a_N. Deshalb passt eine Kleinste-Quadrate-Zerlegung über alle k ≤ N schlecht (Rest 30–70 %), die
  Abel-Koordinaten-Messung (nur Inneres) dagegen auf 1e-14.

## 4. Knesers Singularitätsentwicklung an L (neuer Schwerpunkt) — Rechnung + numerisch

**Herleitung (Skizze, formal vollständig bis auf Konvergenzfragen der Doppelreihe):** Ψ = Schröder-Funktion an L (Ψ(e^z) = LΨ(z),
Ψ(L) = 0, Ψ'(L) = 1), φ = Ψ^{−1} ganz, φ(Ls) = e^{φ(s)}. Reguläre Abel-Funktion u = log Ψ(z)/L (log L = L). Kneser:
A = u + ϑ(u), ϑ 1-periodisch, beschränkt für Im u → +∞ (das ist z → L aus D), also ϑ = Σ_{m≥0} κ_m q^m, q = e^{2πiu} = Ψ^γ.
f_K = T(A + ½) = φ(e^{LU}), U + ϑ(U) = u + ϑ(u) + ½ ⇒ U = u + ½ + ε, ε = ϑ(u) − ϑ(U) = Σ κ_m q^m (1 − (−1)^m e^{2πimε}) ⇒
ε = 2κ₁q + 4πiκ₁²q² + (2κ₃ − 8πiκ₁κ₂ − 16π²κ₁³)q³ + …, und mit e^{L/2} = √L:

    f_K(z) = φ( √L · Ψ(z) · e^{L ε(Ψ(z)^γ)} ).

Folgen: (i) der analytische Teil von f_K an L ist die reguläre Schröder-Halbiterierte f_reg = φ(√L Ψ) (κ₀ fällt heraus);
(ii) Singularitäten genau mit Exponenten 1 + nγ + j (n ≥ 1, j ≥ 0), keine Logarithmen (γ nicht reell);
(iii) Leitterm d₁ t^{1+γ}, **d₁ = 2L^{3/2}κ₁** (L^{3/2} := L·e^{L/2}).
Zweig: t^α := (−L)^α (1 − z/L)^α mit Hauptzweigen; in D nahe L stimmt das mit dem Zweig von u (Hauptzweig von log Ψ) überein.

**Numerik (`k5_singularitaet.py`, `k6_darboux.py`):**
- φ-Rekursion φ_k = L R_k/(L^k − L) exakt; Prüfung φ(Ls) = e^{φ(s)} auf 2.5e-71. Ψ als Reihe (Newton), f_reg-Koeffizienten
  c₁ = 0.227781627188 − 0.086224117118i, c₂ = 0.013294291655 − 0.007673409723i, c₃ = 0.001851727926 − 0.002159093656i
  (|c_m|^{1/m} → ≈ 0.63–0.7, Konvergenzradius von f_reg an L ≈ 1.4–1.6).
- κ_m aus dem Motor: z = φ(e^{Lu}) auf u ∈ [u₀, u₀+1) + iv, v = 1 und 1.5 (|z − L| ≤ 0.17 bzw. 0.09); ϑ = A(z) − u:
  κ₀ = 0.4884150884437966033074 + 0.9223068629260569360028i,
  κ₁ = −8.352159682325323364e-6 − 2.756799271724180753e-4 i,
  κ₂ = 2.79649047063e-9 − 1.26586334878e-7 i, κ₃ = 7.9031e-12 − 5.7247e-11 i, κ₄ = 7.118e-15 − 2.523e-14 i.
  Negative Moden ≤ 2e-52 (Zweig richtig), Höhenvergleich rel. 3e-49 (κ₀), 1e-41 (κ₁), 3e-34 (κ₂), 3e-26 (κ₃).
- Prüfung der Formel am Motor: |f_K − φ(√L Ψ e^{Lε})| ≤ 5e-37 bei |z − L| = 0.11…0.14 (ε bis q⁵); nur f_reg: Fehler 5e-8.
- d₁ = 8.17118842880927215457e-4 + 3.50075914309583953e-4 i, d₂ = 1.76741354298e-6 + 6.1083780512e-7 i,
  d₃ = 5.1406e-9 + 1.4355e-9 i; e_{1,1} = 0.6631957864 + 1.010512463i, e_{1,2} = −0.5590128626 + 0.7943547907i.
- (1 − z/L)-Form: D₁ = d₁(−L)^{1+γ} = −0.032003287962207982 − 0.011244050502060257i.
- **Darboux-Abgleich** a_k = 2Re Σ_{n≤3, j≤20} d_n e_{n,j} (−L)^α L^{−k}Γ(k−α)/(Γ(−α)k!): relative Abweichung 1.4e-7 (k = 60),
  6.7e-12 (100), 1.6e-16 (200), 1.5e-18 (300), 1.1e-21 (500), 3e-24 (800). Nur Leitterm (n = 1, j = 0): 1–10 %.
  ⇒ a_k = 2Re[D₁ L^{−k} k^{−2−γ}/Γ(−1−γ)](1 + O(1/k)), Hüllkurve 2|D₁/Γ(−1−γ)| k^{−6.447}|L|^{−k} = 71.57 k^{−6.447}|L|^{−k}.
- Datei: `gemeinsam/kneser_singularitaet_L.json` (c_m bis m = 24, κ_m, d_n, e_{n,j} bis j = 24, D₁).

## 4b. Amplitudengesetz c_N und Innenschicht an L — numerisch

**Exponent (`k7_amplitude.py`, alle N = 30…194, Linie y = 0, T(w) vom Motor):** Anpassung log c_N = a + β log N + b/N +
Schwingung(N·arg L)/N über N = 60…194 ergibt **β = −7.8877 − 2.1069i**; Soll 1 − 2γ = −7.8939 − 2.1159i (Abweichung 0.1 % / 0.4 %).
Ohne 1/N-Glied: −7.941 − 2.182i. Also c_N = C·N^{1−2γ}(1 + O(1/N)) mit einer Schwingung der Periode 2π/arg L = 4.70 in N (±3 %).

**Konstante (`k9_konstante.py`):** Strang Bs Regel c_N ≈ K·N·u_N/(δ_{−1})_N mit u_N = [z^N]D₁(1 − z/L)^{1+γ} (Knesers Leitterm)
und (δ_{−1})_N = [z^N]D_δ(1 − z/L)^{1−γ}, D_δ = L^{3/2}e^{−2πiκ₀}(−L)^{1−γ}. Die lokale Form δ_{−1} ≈ L^{3/2}e^{−2πiκ₀}t^{1−γ} ist
am Motor geprüft (`k12_delta_lokal.py`: Quotient 0.863, 0.965, 0.991, 0.998 bei |t| = 0.077, 0.019, 0.0046, 0.0011).
K_N = 0.0029 − 0.0468i (N = 45), 0.0012 − 0.0496i (120), 0.0020 − 0.0519i (194); Mittel N = 166…194: 0.0025 − 0.0507i;
grob in 1/N extrapoliert K_∞ ≈ 0.0016 − 0.053i. (Fast rein imaginär; geschlossene Form unbekannt.)
Damit in Kneser-Daten: **c_N ≈ K · N^{1−2γ} · (D₁/D_δ) · Γ(−1+γ)/Γ(−1−γ)**.

**Innenschicht (`k11_innen.py`):** Φ_N(v) := (p_N − f_K)(L(1 − v/N))·N^{γ} (f_K nahe L aus der exakten lokalen Formel):
Φ_N(8) = −0.00227 (N = 40), −0.00305 (80), −0.00326 (100), −0.00336 (200) — konvergiert; Φ_N(4) ≈ −0.025 stabil;
näher an L schrumpft Φ_N langsam (Φ_N(0): −0.60 → −0.43 → −0.32 bei N = 40, 100, 200). Freier Fit p_N(L) − L ~ N^{−4.845 − 1.090i}
(Strang B: N^{−4.85}, Eckenwinkel 201°). Deutung: in der Anpassungszone |z − L| ≈ (4…8)|L|/N ist p_N − f_K ≈ N^{−γ}Φ(v) —
**N-mal grösser** als Knesers Singulärteil dort (~N^{−1−γ}); Anschluss an c δ_{−1} ~ c (v/N)^{1−γ} erzwingt c ~ N^{−γ}·N^{1−γ}
= N^{1−2γ}. Der Faktor N ist also eine Eigenschaft der linearen Innenaufgabe (Strang B: generisch für jeden Exponenten s,
Amplitude ~ N^{2−γ−s}); im Zentrum (v ≲ 2) kommt eine Eckschicht dazu.

**Blindtest der Vorhersage (`k10_fortsetzung.py`, p_N für N = 195…232, Newton wie polynom_wurzel.py bei 600 Bit,
Abbruch bei Stillstand ~1e-103):** Vorhersage aus k9 (vor der Rechnung, K = Mittel 166…194): p_N(1) − T(½) = −1.03e-19 (196),
−6.99e-20 (200), −1.82e-20 (210), Nulldurchgang zwischen N = 216 und 218. Gerechnet: −9.65e-20 (196), −6.12e-20 (200),
−1.70e-20 (210), … (Ergebnis unten in §4c).

## 4c. Ergebnis des Blindtests — numerisch

| N | 196 | 200 | 204 | 210 | 214 | 215 | 216 | 217 | 218 | 219 | 220 | 224 | 230 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gerechnet | −9.65e-20 | −6.12e-20 | −3.79e-20 | −1.70e-20 | −7.8e-21 | +3.8e-22 | −3.7e-21 | +5.6e-21 | +1.0e-22 | +9.6e-21 | +3.8e-21 | +1.07e-20 | +1.88e-20 |
| Vorhersage | −1.03e-19 | −6.99e-20 | −4.48e-20 | −1.82e-20 | −5.8e-21 | — | −9.4e-22 | — | +3.3e-21 | — | +6.8e-21 | +1.24e-20 | +1.78e-20 |

Der vorhergesagte Vorzeichenwechsel tritt ein (Pendelzone 215–218 wegen der 4.7-Schwingung, ab 219 stabil positiv); die Beträge
treffen auf 5–15 %, was der ignorierten ±3-%-Schwingung von c_N entspricht. Die Polynome p_195 … p_232 (70 Stellen) liegen in
`gemeinsam/polynome_195_bis_232.json`.

**Adjungierte Mode (`k13_adjungiert.py`).** c_N^lin = −⟨m_N, F(q_N)⟩ mit m_N = J_N^{−T}ℓ trifft c_N aus p_N auf 7e-9 (60), 1e-9 (100),
4e-10 (150) — die Linearität (B) gilt also auch für die Modenamplitude bis O(c_N). Aber Σ_k|m_k F_k|/|c_N| = 3.7e10, 4.6e18, 2.9e27
(N = 60, 100, 150), grösste Einzelbeiträge aus den obersten Gleichungen k ≈ N − 5; in Strang Bs a₀-Jet-Basis (c = ⟨Hᵀm, j_{a₀}⟩)
noch 8.8e5, 4.2e10, 2.3e15. Dort ist m̃ = Hᵀm im Dualgewicht |L − a₀|^{−i} gleichmässig beschränkt (max|m̃_i| 1.349^{−i} = 0.0671
für alle N), aber der a₀-Jet des Schwanzes wächst in diesem Radius (0.011, 60, 7.9e5), weil q_N ausserhalb von D gross ist.
Fazit (numerisch): **Komponentenweise „Anteile des Rests entlang der Moden“ gibt es in keiner der beiden Koeffizientenbasen** —
c_N entsteht aus exponentieller Auslöschung; eine Abschätzung muss in Funktionsform (Abel-Koordinaten bzw. Innenschicht an L)
geführt werden.

## 4d. Richtung Beweis (Aufgabe 5) — Formulierung, keine Beweise

Was die Numerik festlegt (alle Stränge zusammen): (M1) p_N = Π_N f_K + O(|c_N|²) mit der linearen J-Projektion
Π_N u := J_N(f_K)⁻¹ T_N[eˣ + (f_K'∘f_K)·u] (B); (M2) Π_N f_K − f_K = c_N δ_{−1} + c̄_N δ_{+1} + ρ_N, ρ_N auf Kompakta des Korridors
um 1e-14…1e-20 relativ kleiner (A, B); (M3) c_N = K·N·u_N/(δ_{−1})_N, K ≈ 0.002 − 0.051i (A, B); (M4) J_N⁻¹ in l¹_{|L|,4} beschränkt
(C, D), aber der Konsistenzfehler wächst dort exponentiell (B) — „stabil + konsistent“ greift in diesen Normen nicht.

**Hinreichende strukturierte Abschätzung (Vermutung S).** Es gibt C, N₀ mit: für alle N ≥ N₀ und jedes u, das in D holomorph ist
und an L, L̄ eine Singulärentwicklung Σ D_s (1 − z/L)^s + konj. mit Re s ≥ 1 + Re γ hat,
  (S1) Π_N u − u = c_N(u) δ_{−1} + c̄_N(u) δ_{+1} + ρ_N(u) auf {|z| ≤ r} (r < |L| fest), |ρ_N(u)| ≤ C (r/|L|)^{N/2} ‖u‖;
  (S2) |c_N(u)| ≤ C N^{2 − Re γ − Re s}·Σ|D_s| (für Knesers u = f_K: N^{1 − 2Re γ} = N^{−7.894}).
Zusammen mit einer quadratischen Kontrolle (Strang C: Nichtlinearität entlang des Zuwachses ~ |c_N|²) folgt p_N → f_K lokal
gleichmässig auf D, sogar mit Rate N^{−7.894}. Für die blosse Konvergenz genügt in (S2) jede Schranke o(1).
Bemerkenswert: (S2) verlangt nur Re s > 2 − Re γ = −2.447 — die lineare Auswahl ist sehr robust; schwierig ist allein, dass
J_N⁻¹ die exponentiell grossen Konsistenzfehler exakt auslöscht (Strang B: Faktor H kürzt sich).

**Wie man (S) beweisen könnte (Skizze, Lücken offen):**
1. *Adjungierte Mode.* c_N(u) = ⟨m_N, r_N(u)⟩ mit dem Mode-Funktional ℓ (Fourier-Mittel von θ auf einer Korridorlinie) und
   m_N = J_N^{−T}ℓ, r_N(u) = T_N[Schwanz_N(u)∘f_K] (B: Π_N u − q_N(u) = J_N⁻¹ r_N(u) exakt). (S2) ist dann eine Aussage über
   die adjungierte endliche Sektion mit glatter rechter Seite — vermutlich leichter als über J_N⁻¹ selbst.
2. *Toeplitz-Modell.* Strang D findet lim J_N⁻¹ = Toeplitz mit stückweise stetigem Symbol (Sprung bei arg L). Für endliche
   Sektionen von Toeplitz-Operatoren mit Sprung-/Fisher–Hartwig-Symbolen gibt es genaue Asymptotik (Böttcher–Silbermann, Kap.
   „finite sections, PC symbols“; Fisher–Hartwig-Asymptotik der Inversen). Erwartung (ungeprüft): der Sprungparameter ist γ
   (mod 1), die Randschicht-Inverse liefert genau die Potenzen N^{±γ} und den Faktor N. → Strang C/E: an einem reinen
   Fisher–Hartwig-Modell prüfen, ob „Amplitude ≈ K·N·u_N/mode_N“ mit einem K wie 0.05 herauskommt.
3. *Abel-Koordinaten (Alternative).* Θ_N := A∘p_N∘T − id − ½ erfüllt Θ_N∘P_N + Θ_N = O((w+1)^{N+1}) exakt. Gelingt eine
   a-priori-Schranke |p_N − f_K| ≤ ε_N|z − L|^{−3.447} bis zum Abstand R/N von L (und L̄), so folgt aus dem Bôcher/Fourier-
   Argument |c_N| ≲ ε_N (R/N)^{3.447} → 0. Lücke: diese a-priori-Schranke ist so schwer wie die Schranke (i) im README.

## 5. Sackgassen / Hinweise

- `kneser()` aus experimente/gemeinsam.py schreibt nach jedem Aufruf den 190-MB-Speicher neu; ich benutze den Motor direkt.
  Bei den ersten Tests (dps 60/80) sind ein paar slog-Werte bei dps 60/80 in den gemeinsamen Speicher gelangt (richtig, geprüft).
- Kleinste-Quadrate-Zerlegung von p_N − q_N in abgeschnittene δ_n über alle Koeffizienten: schlecht (Randschicht bei k ≈ N).

- Gewichtete SVD in |L|-Gewichten (l², `k8_resonanz.py`): 1/σ_min = |L| = 1.366 konstant bis N = 80, danach Anstieg (1.49 bei 150,
  2.71 bei 190, Richtung nähert sich d_N) — der Faktor N in c_N ist **kein** kleiner Singulärwert ~1/N in dieser Norm.
- Die Anfrage von Strang D nach ≥ 60 Stellen ist **nicht** erfüllt (Motorgrenze, §1); vorhanden sind ≈ 47 Stellen.

## 6. Nächste Schritte

- Innenaufgabe an L als lineare Skalenaufgabe formulieren (Laplace-Profil auf [0, 1] für p, Dilatation K₀E(v) = √L E(v/√L),
  Bedingung bei a₀) und die Konstante K ≈ 0.002 − 0.051i bzw. den Faktor N daraus herleiten; Vergleich mit Toeplitz-/
  Fisher–Hartwig-Theorie (Strang D: Grenzoperator mit Sprungsymbol).
- 70-stellige Kneser-Koeffizienten nur bei Bedarf (Newton mit slog dps 80, hängende Punkte nahe L einzeln mit Zeitlimit).

---

# Runde 2 (Auftrag der Leitung ~20:30): das lineare Problem P_N = Π_N f_K bei grossen N

## R2.1 Verfahren (`lin_kern.py`, `k21_lauf.py`) — numerisch

P_N := Π_N f_K = J_N(f_K)⁻¹ T_N[eˣ + (f_K'∘f_K)·f_K] (Strang B). Mit der exakten Identität f'∘f = eˣ/f' (aus f∘f = exp) sind
J_N[k][j] = [f^j]_k + [eˣ/f']_{k−j} und die rechte Seite T_N[eˣ + (eˣ/f')·f] aus a_0 … a_{N+1} allein berechenbar. Koeffizienten:
DFT-Werte bis k = 700, darüber die Darboux-Reihe aus §4 (n ≤ 3, j ≤ 20; trifft die DFT-Werte bei k = 800 auf 3e-24 relativ).
arb-Löser (algorithm "approx") mit 0.72·N + 200 Bit (|J| wächst ≈ 1.63^N). Prüfungen: +150 Bit ändern 50 Stellen nicht
(N = 100 … 600); Rauschen 4e-48·1.34^{−k} in den a_k ändert P_N um ≤ 1e-60; P_194, P_232 stimmen mit p_194, p_232 auf 3e-30
bzw. 3e-33 bei x = 1 überein (Linearität, Strang B). Laufzeit 0.1 s (N = 100), 21 s (600), ~110 s (1000).
Berechnet: N = 100 … 500 (Schritt 25), 550 … 1000 (Schritt 50) und Fünfergruppen 296–300, 496–500, 696–700, 996–1000.

## R2.2 Voraussetzung von Satz B-lin — numerisch

Da |P_N·g^a| mit g = (z−L)(z−L̄) (nullstellenfrei in D) subharmonisch ist, liegt sup_D |P_N|(|z−L||z−L̄|)^a auf |z| = |L|
(`k22_auswertung.py`: FFT-Gitter 16N Punkte + 6001 Punkte um arg L). Für **alle** N = 100 … 1000 und a = 0, 1, 2, 3, 3.4 liegt
das Maximum bei z = |L|: S_a = 2.232651·2.904^a (6.48, 18.8, 54.7, 83.8), unabhängig von N. Also hält die Voraussetzung von
Satz B-lin schon mit **a = 0, B = 2.232651** (max|P_N| auf D). Nahe L (|θ − arg L| ≤ 0.3) ist |P_N| ≤ 1.668 (Fensterrand).

## R2.3 Amplitudengesetz bis N = 1000 — numerisch (`k26_fits.py`)

| Fit-Bereich | ohne 1/N | mit 1/N |
|---|---|---|
| N = 100 … 400 | −7.937 − 2.167i | −7.891 − 2.020i |
| N = 300 … 1000 | −7.906 − 2.129i | −7.873 − 2.110i |
| N = 500 … 1000 | −7.902 − 2.125i | −7.865 − 2.119i |
| Soll 1 − 2γ | −7.894 − 2.116i | |

K_N = c_N(δ_{−1})_N/(N u_N), Mittel über Fünfergruppen (die Einzelwerte schwingen mit Frequenz ≈ 2 arg L, Streuung ±0.0007):
0.002304 − 0.051702i (N = 296–300), 0.001935 − 0.052366i (496–500), 0.001796 − 0.052508i (696–700), 0.001739 − 0.052697i
(996–1000). Anpassung K_N = K∞ + b/N + Schwingung (N ≥ 200): K∞ = 0.00173 − 0.05328i, b = 0.10 + 0.50i; mit N^{−0.117}
bzw. N^{−1/2} statt 1/N ergäbe sich K∞ = −0.0004 − 0.0613i bzw. 0.0015 − 0.0543i bei gleichem Rest. Die Gruppenmittel passen
am besten zu 1/N ⇒ **K∞ ≈ 0.0017 − 0.0533i (Vermutung, ±0.001)**. Moden ≠ ±1 bleiben bei |c₀|/|c| ≤ 2e-23.

## R2.4 Randfehler und Eckschicht — numerisch (`k23_kreisfehler.py`, `k24_profil.py`, `k27_form.py`)

- |P_N(L) − L| = 5.5e-10 (100), 1.9e-11 (200), 2.6e-12 (300), 2.2e-13 (500), 4.3e-14 (700), 7.7e-15 (1000); lokale Exponenten
  4.78 … 4.92 (Mittel ≈ 4.86), Fit N ≥ 300: N^{−4.83−1.08i} (mit 1/N) bzw. N^{−4.85−1.10i}. Vergleich: Eckenexponent
  5.447·0.883 = 4.810 (B, C); Imaginärteil eher Im γ = 1.058 als 0.883·Im γ = 0.934.
- max über den Rand |P_N − f_K| (E_N-Koeffizienten bis k = 4000, Schwanz aus der Darboux-Reihe) ist nur ≈ 1 % grösser als
  |E_N(L)|, Ort (θ − arg L)·N = 1.42, 1.15, 1.00, 0.92, 0.84, 0.72 (N = 100, 200, 400, 500, 700, 1000), also knapp jenseits von L
  auf Γout, schrumpfend etwas schneller als 1/N.
- Profil E_N(L(1 − v/N^α)), v = ρe^{iφ}: (i) α = 1, Φ = E·N^γ: für ρ ≥ 8 stabil (Φ(8) = −0.00327 + 0.00028i (N = 100) …
  −0.00320 + 0.00115i (1000), Φ(16) ≈ −0.00021 + 0.00026i), das ist die Anschlusszone an c_N δ_{−1}; im Kern (ρ ≤ 2) fällt Φ
  weiter (Φ(0): −0.43 → −0.17, ≈ N^{−0.41}). (ii) Formfunktion R_N(v) = E(L(1−v/N^α))/E(L): für α = 0.88 Abstand zu N = 1000
  max_{|v|≤4} = 0.36, 0.19, 0.14, 0.069, 0.033, 0.012 (N = 100, 200, 300, 500, 700, 900) — konvergiert etwa wie N^{−1/2};
  für α = 1: 0.67 … 0.024 (langsamer). Werte R_1000 (α = 0.88): R(1) = 0.273, R(2) = 0.086, R(i·1) ≈ 0.217 − 0.847i,
  R(−i·1) ≈ 0.147 + 0.954i (φ = ±90° ≈ Rand). Deutung (Vermutung): zwei Skalen — Kern auf N^{−0.883} (Ecke von W, Strang B/C)
  mit Amplitude ≈ N^{−4.84}, Anschlusszone auf 1/N mit Amplitude N^{−γ}; bei N ≤ 1000 sind sie nur um den Faktor N^{0.12} ≤ 2.3
  getrennt, daher die langsame Konvergenz. Daten für Strang B: `gemeinsam/eckschicht_profil.json`.

## R2.5 Bedingung (F) für Strang C — numerisch (`k25_bedingungF.py`)

Sprung über den radialen Schlitz ζ = L(1+u) aus der exakten lokalen Formel: |Δf|/u^{5.447} = 0.938 (u = 0.005) … 0.760 (u = 0.2),
Leitterm 0.943 ⇒ C_f ≤ 0.95. Motor-Schnitt bei 89.5° (von L aus), nicht radial (76.6°). max|f_K| auf ∂(W + d) ohne
|ζ − L| < 0.25: 3.216 / 3.320 / 3.535 für d = 0.05 / 0.1 / 0.2 (am rechten Ende auf ℝ). Cauchy-Integrale über Kreise im Ring
|L| < |ζ| < |L| + 0.15: ≤ 1.5e-32 relativ (keine Singularitäten).

**Anschlusszone geprüft:** aus c_N ≈ K N u_N/(δ_{−1})_N und δ_{−1} ≈ D_δ (1 − z/L)^{1−γ} folgt für z = L(1 − v/N), v ≫ 1:
Φ(v) = N^γ E_N → C_out v^{1−γ}, C_out = K·D₁·Γ(−1+γ)/Γ(−1−γ) = 2.20 − 4.46i (mit K∞). Vorhergesagt Φ(16) = −0.000217 + 0.000277i,
Φ(32) = 2.2e-6 + 3.2e-5i; gemessen (N = 1000) −0.000210 + 0.000264i bzw. 1.2e-6 + 3.0e-5i (3–7 %); bei v = 8 noch 15 % Innenanteil.

## R2.6 Kippt es für grosse N? Echte Lösungen p_N bis N = 1000 — numerisch (`k28`, `k29`, `k30`, `k31`)

**Antwort: bis N = 1000 kein Anzeichen.** Vorgehen: Newton für F_N(p) = T_N[p∘p − eˣ] = 0 von P_N = Π_N f_K aus (gleiche
Präzision 0.72N + 200 Bit; Präzisionsprüfung N + 250 Bit bei N = 800, 1000: identisch auf 60 Stellen).

| N | 1. / 2. Newton-Korrektur (max|d_k||L|^k) | max_{|z|=|L|}|p_N − P_N| (in L) | |Δ(L)|/|E_N(L)|² | max_{|z|=1}|p_N − P_N| | rel. Unterschied der Modenamplituden | p_N(1) − T(½) |
|---|---|---|---|---|---|---|
| 100 | 1.5e-21 / 3.8e-41 | 3.0e-20 | 0.097 | 6.1e-24 | 1.6e-10 | −1.02e-16 |
| 300 | 3.5e-26 / 1.4e-49 | 1.8e-24 | 0.26 | 1.3e-29 | 2.0e-12 | +1.08e-20 |
| 500 | 2.4e-28 / 1.8e-53 | 2.0e-26 | 0.41 | 3.0e-32 | 2.7e-13 | +2.93e-22 |
| 700 | 9.0e-30 / 4.5e-56 | 1.0e-27 | 0.55 | 5.4e-34 | 7.0e-14 | +1.22e-23 |
| 1000 | 2.9e-31 / 8.5e-59 | 4.5e-29 | 0.75 | 7.9e-36 | 1.7e-14 | −1.56e-25 |

Die Newton-Iteration konvergiert sauber quadratisch (Lösung nicht entartet). Die nichtlineare Korrektur sitzt in L, fällt wie
≈ N^{−8.8} und ist ≈ 0.75·10^{−3}·N·|E_N(L)|² (Verhältnis zu |E_N(L)|² wächst nur wie ≈ N^{0.9}). Newton **direkt von q_N** aus:
N = 300: 4.3e-14, 4.7e-22, 6.5e-42, 5.1e-81; N = 400: 8.0e-15, 2.7e-21, 3.9e-40; N = 500: 4.0e-15, **5.6e-15**, 2.6e-27, 2.0e-51 —
der zweite Schritt wächst mit N (Strang D), die Iteration landet aber auf derselben Lösung (Abstand 0 auf 60 Stellen). Der
Zuwachs ist also eine Eigenschaft des Startpunkts q_N bzw. der Norm von J⁻¹, nicht der Lösung.
Vorbehalt: Dass diese Lösung der von p_1 aus stetig fortgesetzte Zweig ist, ist direkt nur bis N = 232 geprüft (Fortsetzung bis
300 läuft); dazwischen müsste der Zweig von P_N weglaufen, was bei Abständen ~1e-26 nicht zu erwarten ist.

**Versteckter exponentieller Anteil?** Anpassung log|·| = a + b log N + cN über N = 300…1000: für |P_N(L) − L| b = −4.865 ± 0.011,
c = (2.8 ± 1.8)·10⁻⁵ pro Grad; für |c_N| b = −7.933 ± 0.018, c = (4.5 ± 2.9)·10⁻⁵ (Reststreuung 0.004 bzw. 0.006 in log).
Beide c sind mit 0 verträglich und um mehr als den Faktor 200 kleiner als Ds Eckrate 0.0124. Ein Anteil ~e^{0.0124N} müsste bei
N = 1000 relativ unter ~10⁻³ liegen, um hier unsichtbar zu sein (Vermutung: er fehlt, weil die Eingabe die Eckrichtung kaum trifft — D 21:04).
**Zweigprüfung:** Fortsetzung p_232 → p_300 Schritt für Schritt (`k10_fortsetzung.py`, 450 Bit): p_300 stimmt mit der Lösung
von Π_300 f_K aus auf 6e-61 überein (`gemeinsam/polynome_195_bis_300.json`). Echte p_N für N = 100, 200, …, 1000:
`gemeinsam/pN_N100_bis_1000.json`.

## R2.7 N = 1500 und die Rundungsempfindlichkeit — numerisch (`k32_gross.py`, `k33_praez_newton.py`, `k34_volle_praez.py`)

- Π_1500 f_K (1280 Bit, 544–600 s): |c_N| = 1.613e-26 (Vorhersage aus N = 1000 mit N^{1−2γ}: 1.61e-26), K_N = 0.00180 − 0.05298i
  (passt zu K∞ + b/N mit K∞ = 0.0017 − 0.0533i, b = 0.50i), |Π f(L) − L| = 1.082e-15 (Exponent 4.83 gegen N = 1000),
  max_D|Π_N f_K| = 2.232651.
- Echtes p_1500, Newton in voller Präzision (Iterate mit allen 415 Stellen weitergereicht): Korrekturen 5.73e-33 → 7.05e-62,
  ω_eff = 2.1e3 (N = 1000: 1.0e3). |p − Π_N f|(L) = 1.28e-30, |E(L)|² = 1.17e-30.
- Vergleich mit gerundetem Start: vom 60-stellig gespeicherten Π_1500 f_K aus liefert derselbe (bei 1750 Bit identisch
  wiederholte) Newton-Schritt 1.1e-24, der nächste 8.2e-25 — keine Konvergenz. Bei N = 1000 ändert die 60-Stellen-Rundung
  nichts (2.859e-31 in beiden Fällen). Rundungsfehler ~1e-60 werden also bei N = 1500 um ~1e36 verstärkt (bei N = 1000 um
  < 1e30): Verstärkung in allgemeinen Richtungen ≳ e^{0.16 N}, entlang des tatsächlichen Zuwachses nicht (ω_eff ≈ 1e3).
- Folgerungen: Die 60-stelligen Dateien in `gemeinsam/` reichen bis N ≈ 1000; darüber Zwischenstände ungerundet halten.
  Für einen Beweis sind Kugel-Abschätzungen (Kantorovich) aussichtslos; nötig ist die Struktur des Zuwachses (Strang C/D).
