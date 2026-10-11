# Strang M — Modell mit Eckmode

Stand: 11.10.2026, 04:50 — ANGEHALTEN (Runde 2 ab §8, Runde 3 in §11, Runde 4 in §13, Runde 5 in §15, Stand/Wiederaufnahme in §17; Beginn 10.10. 23:29; Zeiten aus `date`). Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.
Ordner `m_modell_eckmode/`, Python Miniconda 3.13 (python-flint/arb, mpmath, numpy). Kern: `m_kern.py`.

Auftrag: einfachstes Modell mit Eckmode (abstossender, drehender Randfixpunkt) finden, dort die Struktur prüfen und (H) — oder
direkt die Konvergenz — vollständig beweisen; übertragbare Schritte benennen.

## 0. Kurzfassung

1. **Modell (bewiesen + numerisch):** das *radiale Quadrat* f(z) = z + ε(z − e^{iθ})(z − e^{−iθ}) auf der Einheitsscheibe. Auf ∂D gilt exakt
   f(e^{iφ}) = s(φ)e^{iφ}, s(φ) = 1 + 2ε(cos φ − cos θ) > 0 reell. Es ist die *einzige* reelle Familie, die auf D̄ holomorph ist und ∂D
   winkeltreu-radial abbildet (Satz M1). Folge: in ℓ² ist C_f e^{ikφ} = s(φ)^k e^{ikφ}, also
   **J̃_N = P_N[Op(s(φ)^k) + M_m]P_N** — ein Pseudo-Differentialoperator ohne Winkeltransport (Einträge = zugeordnete Legendre-Funktionen).
   Exp-Lage θ = 76.618°, ε = 0.31382 (|μ| = |f'(e^{iθ})| = 1.1717 = |√L|).
2. **Eckmode wie bei Kneser (numerisch):** σ₁ = 1.75, 3.22, 7.38, 18.9, 51.2, 142, 404 (N = 200 … 800, Schritt 100 ab 300), Fit-Rate
   **c = 0.0119** (Kneser 0.0120, kubisch 0.0123). Rate hängt bei festem μ vom Fixpunktwinkel ab (θ = 50/60/76.6/90/105/120°:
   c = 0.0093/0.0102/0.0119/0.0133/0.0145/0.0153) ⇒ global, nicht lokal; c ≈ 0.075·ln|μ| bei θ = 76.6°.
3. **Exakte Umformungen (bewiesen, geprüft 1e-68):** (a) Pol-Antwort an jedem ζ ⇔ Kern eines Rechteck-Abschnitts (A_ζ, §4.1);
   (b) **adjungierte Gleichung in Zwei-Punkt-Form** J^Tv = r ⇔ r_{N−l} = [t^l](m W − PP_{a₀}[W(h)h'(y/h)^{N+1}]) (§4.2) — das Gegenstück
   zur Ersatzform für die Eckmoden; (c) Transfergleichung (m + 𝓛)R = Π_N k − J⁻¹Q mit dem Perron-Frobenius-Operator 𝓛 von f (§4.3).
4. **Adjungiertes Randprofil mit vollen Vektoren (numerisch, N ≤ 1200):** σ₁|V₁(θ)| ~ N^{0.6} (polynomial in der Ecke); Raten ab N ≈ 800
   stabil, Zusatzrate Δλ(τ) = 0.0003 (≤ 1°), 0.0010/0.0016/0.0022/0.0032/0.0040/0.0053/0.0069/0.0099 bei τ = 2/3/4/6/8/12/20/60° — dieselben
   Werte wie Ds Kneser/kubisch. Lokaler Exponent 1.56 (1°→2°), 1.18 (2°→4°), ≈ 1 (3°→8°); Schwelle des Modells 1/α = 1.113.
5. **Der Ort der Singularität entscheidet (numerisch, §3a — neu):** u = (1 − z/P)^{1.5}: P in der Ecke ⇒ Π_N u → u; P = −1 auf dem D-Bogen
   (gleiche Stärke) ⇒ **Divergenz wie e^{cN}N^{−2.5}, Fehler sitzt in L** (Exp-Lage ab N ≈ 600: 3.5e-5 → 9.6e-4 bei N = 1000; ε = 0.5: 1.25 bei
   N = 700). Vorhergesagt aus Ds Eckformel. **Das Verfahren funktioniert bei exp, weil Knesers Singularität genau in der Ecke sitzt.**
6. **Nicht gelungen: ein Beweis von (H) in einem Modell mit Eckmode.** Gründe (§5): keine exakte Lösbarkeit (affine und
   „Pol im Fixpunkt“-Modelle kippen in eine Bulk-Explosion statt einer Ecke), die Rate ist global, und jede Normschranke für J_N⁻¹
   verliert e^{cN} genau an der Stelle, an der Knesers Singularität sitzt. Was ein Beweis leisten muss, ist im Modell präzise: das
   Profil der adjungierten Eckmoden nahe L (= W) und die Bulk-Kontraktion (= O, im elliptischen Fall von O schon bewiesen).
7. **Runde 3 — Injektivität und (a″)_∞ im Modell, rechnergestützt bewiesen (§11; Prüfung durch P steht aus):** der Randoperator
   𝒯x = x + x∘f/m (Γin), x + f'·x∘h (Γout) (nur Werte in D̄, ker 𝒯 = ker T) ist in H²(D) — sogar in L²(∂D) mit Poisson-Fortsetzung —
   **akkretiv: Re⟨𝒯x, x⟩ ≥ 0.8211‖x‖² (Satz M6)**. Folgen: **T injektiv auf E²(Ω)**, **σ_min(J_∞|P_N) ≥ 0.4249 für alle N** (Runde 4 verfeinert: ≥ 0.5596; wahr 0.6905).
   Beweis in Os Streifenkoordinate ζ = log((x−L)/(x−L̄)): Eckmodell (Mellin, arb: ≥ 0.99), Endstück (Schwarz-Lemma, ε ≤ 0.044),
   Kommutator der IMS-Zerlegung (Schur, ≤ 0.125), Mittelstück (Nyström 1664², verifizierte Cholesky, Interpolationsfehler ≤ 2e-4: ≥ 0.990).
   (G2) wird nicht gebraucht; für Kneser fehlen nur zertifizierte f_K-Einschlüsse nahe ∂D und die ϑ-Schranke in L (§11.5).
8. **Runde 4 — Kneser (§13), Vorzertifikat:** derselbe Beweis mit f_K läuft durch: Re⟨𝒯x,x⟩ ≥ 0.8062‖x‖², **σ_min(J_∞|P_N) ≥ 0.5645 ∀N**
   (wahr 0.7319; Os O11 numerisch 0.538). Streng relativ zu den Strang-A-Koeffizienten + angenommenem Schwanz; P muss nur liefern:
   (K1) b_k = a_kρ^k, k ≤ 400…1000, als Kugeln, (K2) |b_k| ≤ Ck^{−p} (p > 3) für k > K, (K3) f_K(L) = L, f_K'(L) = √L. Werte ausserhalb D̄,
   h_K und die Eckschranken leitet die Schnittstelle selbst ab (Funktionalgleichung, Krawczyk, S2). Kein (G2), keine Univalenz.
9. **Runde 5 (§15; Stand beim Anhalten §17):** mit Ps zertifizierten a_0…a_640 (Satz P-V) **σ_min(J_∞|P_N) ≥ 0.6553 für Kneser (T₁ = 7, T₂ = 23), streng bis auf die Schwanzschranke (K2)**
   (Bedarf: Σ_{k>640}|b_k| ≤ 1e-9 usw.; Daten: Reserve 10⁶; Gegenprobe mit C = 500, p = 5 unverändert 0.5645). Kante (a″)_{αN}: perturbativ nur
   für α ≥ 1 (radial, exakt) bzw. α ≥ 2.74 (Kneser, = Xs X5); für α ∈ [0.03, 0.1] ist der Überlauf exponentiell (1.48^N bzw. 1.62^N) — eine
   Eck-Aussage vom Typ (K_d), nicht mit 𝒯 allein. P-V-Gegenprobe: Funktionalgleichung ✓, A-Abgleich 1.1e-43 ✓, Rundungsfehler der Datei gefunden.

## 1. Das radiale Quadrat-Modell

**Definition.** D = Einheitsscheibe (Normierung |L| = 1; zurück zu Kneser-Masstab mit x = ρz, ρ = |L|), θ ∈ (0, π), ε > 0 mit
c₀ := 1 − 2ε cos θ > 2ε. f(z) = ε + c₀ z + ε z² = z + ε q(z), q(z) = (z − e^{iθ})(z − e^{−iθ}). g := f∘f. Dann:
- Fixpunkte e^{±iθ} ∈ ∂D, **μ = f'(e^{iθ}) = 1 + 2iε sin θ** (|μ| > 1, arg μ = arctan(2ε sin θ)); a₀ = f(0) = ε; R₀ = |e^{iθ} − ε|.
- **Radialität:** f(z) = (az + b)(a + bz) mit a² + b² = c₀, ab = ε (a, b > 0); auf |z| = 1 ist |a + bz| = |az + b|, also
  f(e^{iφ}) = e^{iφ}|a + be^{iφ}|² = s(φ)e^{iφ}. Äquivalent: f(1/z) = f(z)/z².
- Normalform: w = εz + c₀/2 konjugiert f zu P_c(w) = w² + c, c = ε² + c₀/2 − c₀²/4 > 1/4 (Exp-Lage c = 0.3432: Cantor-Julia-Menge).
- Γin = {|φ| > θ} (s < 1), Γout = {|φ| < θ} (s > 1); Anteil von Γin = 1 − θ/π (Exp-Lage 0.574, wie Kneser).

**Satz M1 (bewiesen).** Sei f holomorph auf einer Umgebung von D̄ mit f(e^{iφ}) = s(φ)e^{iφ}, s reell. Dann ist f(z) = c̄₁ + c₀z + c₁z²
(c₀ reell); ist f reell, so c₁ reell — also genau die obige Familie (bis auf die Lage der Fixpunkte).
*Beweis.* S := f/z ist auf D̄ meromorph mit höchstens einfachem Pol in 0 und reell auf ∂D. Schwarz-Spiegelung S(z) = conj S(1/z̄) setzt S zu
einer rationalen Funktion auf Ĉ fort, deren Pole in D nur 0 (einfach) und daher ausserhalb nur ∞ (einfach) sind: S = c₀ + c₁z + c̄₁/z. ∎
(Geprüft von Strang P, BRETT 00:11.)

**Lemma M2 (bewiesen).** (i) f ist injektiv auf D̄ (kritischer Punkt −c₀/(2ε) hat Betrag > 1 ⇔ c₀ > 2ε ⇔ s > 0), f(D) = Ω_f :=
{re^{iφ} : r < s(φ)} (sternförmig). (ii) Os globale Struktur (S) gilt: f(Γin) ⊂ D (da s < 1) und h(Γout) ⊂ D (da e^{iφ} ∈ Ω_f für s(φ) > 1),
h = f⁻¹: Ω_f → D. (iii) In ℓ² (Gewicht 1) ist (J_N)_{jk} = [s^k]_{j−k} + m_{j−k}, [s^k]_n = (c₀² − 4ε²)^{k/2}(k!/(k+n)!)P_k^n(x₀),
x₀ = c₀/√(c₀² − 4ε²) > 1 (Laplace-Integral); m = f'∘f = c₀ + 2εf. Also J̃_N = P_N[Op_KN(e^{k log s(φ)}) + M_m]P_N.
*Beweis.* (i) quadratisches Polynom, injektiv auf jeder Scheibe um 0, die den kritischen Punkt nicht enthält; Randbild = Kurve r = s(φ).
(ii) direkt aus (i). (iii) f^k = z^k(c₀ + ε(z + 1/z))^k. ∎

**Rückübersetzung auf exp:** Kneser hat Γin-Anteil 0.574, |μ| = 1.1717, arg μ = 38.3° (Modell: 31.4°, da Re μ = 1 erzwungen). Die Ecke bei L
hängt nach Ds Scans an |μ| und an der Γout-Länge, kaum an arg μ — beides stimmt überein.

## 2. Eckmode im Modell (numerisch; m01, m03, m07)

σ-Werte der Inversen in ℓ² (letzte 150 Spalten bzw. volle Inverse in arb, 0.9N + 300 Bit):

| N | 200 | 300 | 400 | 500 | 600 | 700 | 800 |
|---|---|---|---|---|---|---|---|
| σ₁ radial (Exp-Lage) | 1.753 | 3.216 | 7.382 | 18.93 | 51.19 | 142.5 | 404.0 |
| σ₂ | 1.448 | 1.650 | 2.881 | 6.234 | 15.11 | 38.75 | 102.7 |
| σ₁ kubisch (D) / Kneser (D) | 1.85 / 1.84 | 3.69 / 3.64 | 9.09 / 8.90 | – | 70.7 / 68.4 | – | 619 / 592 |

Rest: 1.4481 (= max|1/m| auf Γin für das Modell; σ₃ setzt bei N ≈ 800 ein (1.7055) wie Os „Paar 2“).
Ratenscan (Fit ln σ₁ = A + cN − a ln N, N = 400…800; m03_raten.jsonl):

| θ | ε | |μ| | c (Fit) | lokale Rate 700→800 | c/ln|μ| |
|---|---|---|---|---|---|
| 76.618° | 0.20 | 1.0730 | 0.0059 | 0.0033 | 0.084 |
| 76.618° | 0.25 | 1.1120 | 0.0090 | 0.0064 | 0.085 |
| 76.618° | 0.31382 | 1.1717 | **0.0119** | 0.0104 | 0.075 |
| 76.618° | 0.40 | 1.2672 | 0.0174 | 0.0164 | 0.073 |
| 76.618° | 0.50 | 1.3951 | — | 0.0237 (600→700) | — |
| 50° | 0.3985 | 1.1717 | 0.0093 | 0.0069 | 0.059 |
| 60° | 0.3525 | 1.1717 | 0.0102 | 0.0083 | 0.065 |
| 90° | 0.3053 | 1.1717 | 0.0133 | 0.0119 | 0.084 |
| 105° | 0.3161 | 1.1717 | 0.0145 | 0.0133 | 0.092 |
| 120° | 0.3525 | 1.1717 | 0.0153 | 0.0142 | 0.097 |

Die Rate ist bei festem μ eine Funktion des Fixpunktwinkels θ (Γout-Länge 2θ) — also eine **globale** Grösse (Bestätigung von D).
Keine geschlossene Formel gefunden (c/ln|μ| ≈ 0.065√θ passt auf ±4 %, ohne Begründung — nicht behauptet).

**Affines Modell zum Vergleich** (f = L + μ(x − L), μ = √L; m01): σ₁ = 4.5e5, 1.1e12, 1.5e18, 1.6e24 (N = 20…80), u₁ bei 237° — Bulk-Explosion,
keine Ecke (Γin-Anteil 0.32; (S) verletzt; vgl. O §4). Ebenso das „Pol im anziehenden Fixpunkt“-Möbius-Modell (rationale Funktionen mit
Polen nur in L̄): dort ist J_N exakt diagonal (Basis w^k, w = (x−L)/(x−L̄), Eigenwerte μ^k + μ), aber E_N wertet Q_ζ auf ∪_n μ^{−n}H aus
(alle Halbebenen durch 0 ⇒ ganz ℂ) — exakt lösbar, aber gescheitert. **Exakte Lösbarkeit und Eckmode schliessen sich in diesen Familien aus:**
exakt lösbar wird J_N genau dann, wenn der Pol des Ansatzraums in einem Fixpunkt liegt, und dann fehlt die Abschneidekante, die (zusammen mit
der Drehung) die Ecke erzeugt.

## 3. Gestalt der Eckmoden (numerisch; m02, m04, m07)

- u₁ (Ausgang): auf dem Kreis bei θ = 76.6° (|u₁| = 7.6, N = 400), zur Γout-Seite sehr schnell abfallend (70°: 1e-3, 60°: 6e-6), zur Γin-Seite
  langsam (90°: 0.72, 180°: 0.15). Koeffizienten bei k/N ∈ [0.5, 0.95] (90 % der Masse in k ≥ 0.84N). Phasenraum: Hauptteil bei (θ, κ ∈ [0.5, 0.85]),
  Nebenteil in der Deckschicht κ ≈ 1 über Γin.
- v₁ (Eingang, adjungierte Mode): Koeffizienten bei k ∈ [N−25, N−10] (N = 400), v_N ≈ 4e-11; auf dem Kreis |V| ≈ 1.3–1.6 auf Γin (120°–240°),
  0.26 in θ, ≤ 1e-5 auf Γout. Wie D §9.2 für kubisch/Kneser.
- **Randprofil mit vollen Vektoren** (m07; m05 mit nur 150 Spalten war nahe L um 30 % falsch — nicht verwenden), log(σ₁|V₁(θ+τ)|):

| τ | −2° | −1° | 0 | 0.5° | 1° | 2° | 4° | 6° | 8° | 12° | 20° | 30° | 60° | 90° |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N = 400 | −3.62 | −0.55 | 0.668 | 0.964 | 1.117 | 1.182 | 1.185 | 1.238 | 1.307 | 1.460 | 1.741 | 2.019 | 2.437 | 2.478 |
| N = 600 | −4.73 | −1.27 | 0.920 | 1.315 | 1.457 | 1.497 | 1.645 | 1.844 | 2.052 | 2.440 | 3.057 | 3.596 | 4.372 | 4.529 |
| N = 800 | −6.05 | −2.55 | 1.078 | 1.554 | 1.675 | 1.771 | 2.113 | 2.499 | 2.868 | 3.510 | 4.466 | 5.264 | 6.397 | 6.666 |
| N = 1000 | −6.91 | −5.11 | 1.196 | 1.737 | 1.847 | 2.046 | 2.619 | 3.205 | 3.737 | 4.630 | 5.917 | 6.971 | 8.457 | 8.839 |
| N = 1200 | −7.48 | −4.20 | 1.291 | 1.885 | 2.007 | 2.333 | 3.158 | 3.943 | 4.636 | 5.775 | 7.391 | 8.698 | 10.537 | 11.031 |

  In der Ecke polynomial (σ₁|V₁(θ)| ≈ N^{0.5…0.6}). Raten (1000→1200) minus Eckrate 0.0005: Δλ = 0.00027 (0.5°), 0.00033 (1°), 0.00097 (2°),
  0.0016 (3°), 0.0022 (4°), 0.0028 (5°), 0.0032 (6°), 0.0040 (8°), 0.0047 (10°), 0.0053 (12°), 0.0060 (15°), 0.0069 (20°), 0.0082 (30°),
  0.0093 (45°), 0.0099 (60°), 0.0105 (90°); (800→1000) auf 2–5 % gleich ⇒ stabil. Lokaler Exponent: 1.56 (1°→2°), 1.18 (2°→4°), 1.0 (3°→6°),
  0.86 (4°→8°); unter 1° ein Sockel 0.0003 (verträglich mit einem τ-abhängigen Polynomfaktor). **Schwelle im Modell:** Winkel der Kreise
  ∂D, ∂D(a₀,R₀) in L ist 18.22°, Aussenwinkel απ mit α = 0.8988, 1/α = **1.113** (K'-Rahmen nach P strenger). Zum Vergleich D (Kneser/kubisch,
  N ≤ 1200): Δλ = 0.0003 (0.5°), 0.0008–0.0009 (2°), 0.0020–0.0022 (4°), 0.0030–0.0033 (6°) — **dieselben Zahlen**: das radiale Modell trifft
  auch das Profil. Ob p asymptotisch > 1/α ist, lässt sich mit N ≤ 1200 nicht entscheiden (der Exponent hängt stark davon ab, welches
  τ-Intervall man nimmt). σ₃,₄ (zweites Paar) ab N ≈ 800: 1.71, 5.24, 25.3 (N = 800/1000/1200), mit flacherem Profil (≈ +1.2 von L bis 45°).

### 3a. Der Ort der Singularität (numerisch, m09) — warum exp funktioniert
Lineares Problem Π_N u = J_N⁻¹T_N[u∘f + m u] mit u = (1 − z/P)^σ, σ = 1.5 (bzw. Summe mit konjugiertem P), Fehler sup_{|z|=1}|Π_N u − u|:

| Modell | P | N = 300 | 400 | 500 | 600 | 700 | 800 | 1000 |
|---|---|---|---|---|---|---|---|---|
| ε = 0.5 (c ≈ 0.024) | e^{±iθ} (Ecke) | 1.5e-4 | 1.0e-4 | 7.8e-5 | 6.1e-5 | 5.0e-5 | | |
| ε = 0.5 | −1 (D-Bogen) | 4.5e-4 | 2.8e-3 | 1.9e-2 | **0.15** | **1.25** | | |
| ε = 0.5 | a₀+R₀ = 1.509 (A-Scheitel) | 1e-15 | | 1e-15 | | 1e-15 | | |
| Exp-Lage (c = 0.0119) | e^{±iθ} | | 6.8e-5 | | 3.8e-5 | | 2.5e-5 | 1.8e-5 |
| Exp-Lage | −1 | | 3.7e-5 | | 3.5e-5 | | 1.7e-4 | **9.6e-4** |
| Exp-Lage | a₀+R₀ = 1.290 | | 2.9e-12 | | 7e-15 | | 1e-15 | 1e-15 |

Bei P = −1 sitzt das Fehlermaximum ab dem Umschlag (ε = 0.5: N ≈ 400; Exp-Lage: N ≈ 600) in L bzw. L̄ — es ist die Eckmode u₁. Die Rate
(Exp-Lage 800→1000: ln(5.7)/200 = 0.0087) passt zu c − 2.5·d(ln N)/dN = 0.0119 − 0.0028. **Vorhersage vor der Rechnung** (aus Lemma D6):
Eckanteil ≈ σ₁·V₁(P)·m(P)/P·Γ(σ+1)N^{−σ−1}·u₁, V₁ ~ 1 auf Γin; ε = 0.5, N = 600: Prognose ≈ 0.1, gemessen 0.15.
**Folgerung (numerisch belegt, Mechanismus = D6):** Eine schwache Singularität an einem Punkt des D-Bogens zerstört das Verfahren, dieselbe
Singularität in der Ecke nicht; auf dem A-Bogen ist sie harmlos (V₁ fällt ausserhalb von D̄ wie |z|^{−0.95N}). Knesers f_K hat auf ∂W nur L, L̄ —
das Verfahren funktioniert bei exp *wegen* des Orts der Singularität. Ein Beweis muss das benutzen (kein Argument, das nur die Stärke der
Singularität sieht, kann reichen; vgl. X §2.4).

## 4. Exakte Umformungen (bewiesen)

### 4.1 Pol-Antwort als Kern eines Rechteck-Abschnitts
Für ζ ∉ {0, a₀} sei A_ζ(z) := −(ζ − z)E_N(ζ, z), E_N = Π_N k_ζ − k_ζ. Dann ist A_ζ ein Polynom vom Grad ≤ N+1 mit A_ζ(ζ) = 1 und
  T_N[(ζ − z)·A_ζ(f(z)) + m(z)(ζ − f(z))·A_ζ(z)] = 0,
und umgekehrt bestimmt diese Bedingung A_ζ eindeutig, wenn J_N invertierbar ist. *Beweis.* E + k_ζ ∈ P_N ⇔ A = 1 − (ζ − z)p mit p ∈ P_N;
J_∞E = O(z^{N+1}) ⇔ multipliziert mit (ζ − z)(ζ − f) die angegebene Bedingung (beide Faktoren ≠ 0 in z = 0 für ζ ≠ 0, a₀). ∎
Für ζ = L vereinfacht sich das (ζ − f = (L − z)(μ − ε(L − z)) im Modell) zu T_N[A∘f + m·q_L·A] = 0, q_L = (f − L)/(z − L):
**die Pol-Antwort in der Ecke ist der eindimensionale Kern eines (N+1)×(N+2)-Abschnitts mit lokalem Multiplikator μ²**.
(H₀) in der Ecke ⇔ |A_L| ≤ C N^κ |L − z| auf D̄. Numerisch (m06, N = 300): |A_L| ≈ 1e-2 im Koeffizientenmittel, auf dem Kreis 1 bei L,
0.02 bei 180°, 2e-3 auf r = 0.95 — polynomial, keine Spur der Eckmode.

### 4.2 Adjungierte Gleichung in Zwei-Punkt-Form (das Gegenstück zur Ersatzform)
Sei f holomorph nahe 0 mit a₀ = f(0), f'(0) ≠ 0, h = f⁻¹ der Zweig mit h(a₀) = 0, m = f'∘f. Für v ∈ ℂ^{N+1} sei W(t) := Σ_i v_{N−i} t^i
(umgekehrte Reihenfolge; D's V(x) = x^{−N}W(x)) und Φ_W(y) := W(h(y))·h'(y)·(y/h(y))^{N+1}. Dann gilt exakt
  **(J_N^T v)_{N−l} = [t^l]( m·W − PP_{a₀}[Φ_W] )   (l = 0, …, N),**
PP_{a₀} = Hauptteil in a₀. Gleichwertig: J^T v = 0 ⇔ T⁽⁰⁾_N[(t − a₀)^{N+1} m W] = T^{(a₀)}_N[W(h)·h'·Λ^{N+1}], Λ(t) := t(t − a₀)/h(t).
*Beweis.* (J^Tv)_{N−l} = [t^l](W g^{N−l}) + [t^l](mW) mit g = f/t (Laurent). Erzeugende Funktion Σ_l t^l[z^l](Wg^{N−l}) =
(1/2πi)∮_{|z|=r} W f^{N+1}/(z^{N+1}(f − t)) dz für kleines |t|; Substitution y = f(z) (f injektiv auf der kleinen Scheibe) und Residuen in y = t
(liefert W(h)h'(t/h)^{N+1} = O(t^{N+1})) und y = a₀ (liefert −PP_{a₀}[Φ_W](t)). ∎ Geprüft (m08, mpmath 80 Stellen, Zufallsvektor):
radial N = 12: 1.6e-68, N = 25: 1.4e-52; kubisch N = 15: 9.3e-68. Gilt wörtlich für Kneser (f Potenzreihe, h = f_K⁻¹ bei a₀).
*Lesart:* Die adjungierte Eckmode ist ein Polynom W, dessen m-faches beim Entwicklungspunkt 0 bis zur Ordnung N mit dem Hauptteil des nach a₀
transportierten Φ_W übereinstimmt; die Gewichte |(t − a₀)|^{N+1} bei 0 und |Λ|^{N+1} bei a₀ sind die adjungierten Gegenstücke der Zwei-Kreis-
Gewichte |z|^N, |z − a₀|^N der direkten Aufgabe. Das könnte die richtige Startgleichung für Ws Grenzschicht sein (keine riesigen Matrixeinträge;
J hat Diagonale ~1.48^N, z. B. J_{400,400} = 7.4e66, die adjungierte Mode löscht das aus).

### 4.3 Transfergleichung für die rohe Antwort
R_ζ := J_N⁻¹T_N k_ζ (Erzeugende der Spalten von J_N⁻¹). Mit J_∞k_ζ = m(ζ)k_ζ + Σ_{f(x_j)=ζ} k_{x_j}/f'(x_j) + Q_ζ (Q_ζ Polynom) gilt
  m(ζ)R_ζ + 𝓛R(ζ) = Π_N k_ζ − J_N⁻¹Q_ζ,   𝓛R(ζ) := Σ_j R_{x_j(ζ)}/f'(x_j(ζ))  (Perron-Frobenius-Operator von f in ζ).
Projiziert auf u_i ist das Ds Eckformel D6. Im Modell ist f eine hyperbolische quadratische Abbildung (c ausserhalb der Mandelbrotmenge) —
für 𝓛 auf holomorphen Funktionen nahe der Julia-Menge gibt es eine Spektraltheorie (Ruelle); ob sie hier greift, ist offen.

## 5. Bewertung: warum (H) im Modell nicht bewiesen ist, und was genau fehlt

- **K'-Kette mit grober Inversen-Schranke** (Os O5 + Ps K' = D̄ ∪ f(D̄)): Mit ‖J_N⁻¹‖ ≤ C e^{c'N} folgt |E_N(ζ,z)| ≤ CN^κe^{(c' − g_{K'}(ζ))N} —
  Konvergenz für Eingaben, die über {g_{K'} ≤ c'} hinaus holomorph sind. Für eine Singularität **in** L versagt das: g_{K'} → 0 in L, das
  Cauchy-Integral über den Schlitz gibt ~e^{c'N}δ^{6.4}. Jede Normschranke verliert genau dort.
- **Folgerung (Skizze, numerisch gestützt):** Eine Kneser-artige schwache Singularität an einem *anderen* Punkt P des D-Bogens würde das
  Verfahren zerstören: der Eckanteil der Pol-Antwort ist σ₁·V₁(P)·m/P ~ e^{cN} (V₁ ~ 1 auf Γin), das Cauchy-Integral ~e^{cN}N^{−(σ+1)}.
  **Das Verfahren funktioniert bei exp, weil Knesers Singularität genau in der Ecke sitzt — dem einzigen Punkt von Γ̄in, an dem die
  adjungierte Eckmode exponentiell klein ist.** (Cs „Spitze“ auf dem linken D-Bogen ist dieselbe Erscheinung.)
- **Lokal reicht nicht:** Os Satz O3 (keine mellin-diagonale Norm), Ds λ-Test und die θ-Abhängigkeit der Rate (§2) zeigen: die Rate c ist eine
  globale Grösse der Kreislage; ein Beweis von (E_L) braucht die globale WKB-Struktur (Strang W).
- **Exakte Lösbarkeit:** in allen geprüften Familien (affin, Möbius mit Pol im Fixpunkt, Trappmann-Typ) schliessen sich exakte Lösbarkeit und
  Eckmode aus (§2). Möbius-Modelle sind laut X exakt Toeplitz × SL(2)-Darstellung (Krawtchouk); im loxodromischen Fall (Ds d24) gibt es
  eine Ecke (Rate 0.077) — der aussichtsreichste Kandidat für eine Riemann–Hilbert-Analyse, aber nicht reell.

- **Bezug zu X (BRETT 00:38):** Eckmoden sind ein Effekt der quadratischen Abschneidung; Rechteck-Abschnitte T_{N+d}J_∞|_{P_N} sind gleichmässig
  injektiv. Im radialen Modell hat das Symbol |s(φ)^{κN} + m(φ)|² überall den Wert ≥ min_{Γin}|m|² (an der Ecke |t + μ| ≥ |μ| für t ≥ 0, da Re μ = 1 > 0)
  — (a″) ist dort genau eine Gårding-Ungleichung für Op(s^k) + M_m ohne Ausgangs-Abschneidung. Das ist der aussichtsreichste *beweisbare* Teil.

## 6. Nächste Schritte

1. (a″) im radialen Modell beweisen: inf_{x ∈ P_N} ‖x∘f + m x‖_{H²}/‖x‖_{H²} ≥ min_{Γin}|m| − o(1) (X-Lemma ⇒ (a′) ohne Eckmoden-Theorie).
   Radiale Struktur: (x∘f)(e^{iφ}) = P_φ(s(φ)) mit dem radialen Polynom P_φ(r) = Σ x_k e^{ikφ}r^k.
2. Mit W: WKB im radialen Modell aus der Zwei-Punkt-Form §4.2 (Phasen |t − a₀|, |Λ|), Ziel Formel für c(ε, θ) und für das Profil;
   der D-Bogen-Test (§3a) als quantitative Probe für jede Profil-Theorie (Umschlagpunkt N ≈ 600 in Exp-Lage).
3. Mit O: „Akkretivität modulo Ecke“ im radialen Modell (Op(s^k) ist ein reiner Pseudo-Differentialoperator; Bulk + Mellin-Teil explizit).
4. Loxodromisches Möbius-Modell (X: Krawtchouk-Struktur) auf Riemann–Hilbert-Tauglichkeit prüfen.

---

# Runde 2 (Auftrag 00:55): (a″) im radialen Modell, Ort-Test als Prüfstein, Zwei-Punkt-Form an W

## 8. (a″) im radialen Modell — Stand und Korrektur

**8.1 Korrektur (bewiesen im Sinn einer Klarstellung):** Mein Hinweis „Gårding-Ungleichung für Op(s^k) + M_m mit Symbol |s^{κN} + m|²“ ist kein
Beweisweg. In semiklassischer Skalierung (h = 1/N, κ = k/N) ist das „Symbol“ e^{κ log s(φ)/h}: eine Komposition = Fourier-Integraloperator mit
imaginärer Phase, kein Pseudo-Differentialoperator (Δ_k a/a = s − 1 = O(1), ∂_φ a/a = O(k)). Die Heuristik trifft den *Wert* (siehe 8.3), mehr nicht.
An der Ecke ist zudem jede Akkretivitäts-/Kontraktionsform falsch (Os O3; auch in der Randdaten-Form 8.4 wird die Auswertung auf der inneren
Kurve nahe L zur Mellin-Dilatation μ^{s}, |·| unbeschränkt in σ).

**8.2 Ausgangszeilen (numerisch, m10; Gewicht 1, Exp-Lage):** σ_min(T_{N+d}J_∞|_{P_N}):

| N | d = 0 | 1 | 2 | 3 | 4 | 6 | 8 | 12 | voll (d = N+2) |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.690730 | 0.690741 | 0.690746 | 0.690748 | | 0.690752 (d=5) | | | 0.690764 |
| 200 | 0.570513 | 0.690573 | 0.690580 | 0.690581 | | 0.690581 (d=5) | | | 0.690586 |
| 400 | 0.135467 | 0.277717 | 0.448748 | 0.606699 | 0.690539 | 0.690539 | 0.690539 | 0.690539 | |
| 600 | 0.019533 | 0.049280 | 0.099676 | 0.173068 | 0.268733 | 0.500909 | 0.690526 | 0.690532 | |

min_{Γin}|m| = min_{∂D}|m| = **0.690525** (bei 143.94°). Jede Zusatzzeile hebt den kleinsten Singulärwert um einen Faktor 2.5 → 1.4 (N = 600:
×2.52, 2.02, 1.74, 1.55, dann 1.37/Zeile), bis er den Bulk-Wert erreicht: d(400) = 4, d(600) = 8 — passt zu X' d ≈ N/60. Die Eckmoden selbst
haben riesige Bilder oberhalb N (m12: |[J_∞x]_{N+1}|/‖x‖ ≈ 3e9 bei N = 400); die Minimierer für kleines d sind *tiefere* Kombinationen, die
weitere Zeilen auslöschen. **Folgerung:** (a″)_d mit d ≈ r(N) ist dieselbe Aufgabe wie (a′)_r — X' Lemma verschiebt die Ecke in die d×d-Matrix
S_N, entfernt sie nicht. Ohne Kantenanalyse ist nur (a″)_∞ (= Invertierbarkeit des unendlichen Operators, Os O7 + Injektivität) erreichbar.

**8.3 Lemma M3 (Schärfe; bewiesen bis auf eine numerisch geprüfte Konstante).** Im radialen Modell gilt für alle N, d ≥ 0:
  σ_min(T_{N+d}J_∞|_{P_N}) ≤ min_{Γin}|m| + C/N.
*Beweis.* φ₀ := Minimalstelle von |m| auf Γin, n := ⌊N/3⌋, x(z) := zⁿ(1 + e^{−iφ₀}z)^{2n} ∈ P_N. Auf ∂D ist |x| = |2cos((φ−φ₀)/2)|^{2n},
‖x‖² = binom(4n,2n) ≍ 2^{4n}n^{−1/2}, und ‖m x‖²/‖x‖² = |m(φ₀)|² + O(1/n) (Laplace-Methode, |m|² glatt mit Minimum in φ₀). Ferner
|x(f(e^{iφ}))| = (s(φ)^{1/2}|1 + s(φ)e^{i(φ−φ₀)}|)^{2n} ≤ (2F)^{2n} mit F := max_φ s^{1/2}|1 + se^{i(φ−φ₀)}|/2 = **0.8751** (numerisch, Gitter 4·10⁵,
bei φ = −56.3°; F < 1 ist mit Intervallarithmetik trivial nachprüfbar). Also ‖x∘f‖/‖x‖ ≤ C n^{1/4}0.8751^{2n} → 0 und
‖T_{N+d}J_∞x‖ ≤ ‖J_∞x‖_{L²(∂D)} ≤ ‖x∘f‖ + ‖mx‖ ≤ (|m(φ₀)| + C/N)‖x‖. ∎ Damit ist 0.690525 die **scharfe** Konstante in (a″) (Daten 8.2: von oben).

**8.4 Randdaten-Form von T (bewiesen, Rechnung):** Ω = D ∪ f(D) ist sternförmig, ∂Ω(φ) = R(φ)e^{iφ} mit R = max(1, s); die innere Kurve
γ(φ) = r(φ)e^{iφ}, r = min(1, s). Für x ∈ E²(Ω) mit Randdaten u(φ) = x(∂Ω(φ)) gilt
  (Tx)(e^{iφ}) = A(φ)·[u(φ) + B(φ)·x(γ(φ))],  A = 1, B = 1/m auf Γin;  A = 1/m, B = m auf Γout,
also T = M_A(I + M_B E_γ) mit der Auswertung E_γ auf der inneren Kurve (strahlweise: γ(φ) und ∂Ω(φ) liegen auf demselben Strahl; das ist die
radiale Struktur). Re B > 0 überall (Re m ≥ 0.68 auf ∂D). E_γ ist kompakt fern von L, L̄ (γ dort im Inneren von Ω); an L, L̄ ist I + M_BE_γ
genau Os lokales Modell T_L (Mellin-Multiplikator 1 + μ^{s−1}, μ = 1 + 2iε sin θ). Fredholm-Eigenschaft im radialen Modell also wie O7, ohne
die Voraussetzung (S) prüfen zu müssen (Lemma M2). **Injektivität bleibt offen** (gleichwertig mit einer Eindeutigkeitsaussage vom Typ D4 für das
Quadrat; die Julia-Menge des Quadrats schneidet D in einer Spirale bei L, darum greift der Böttcher-Abel-Ansatz nicht direkt).

## 9. Ort-Test als quantitativer Prüfstein für (E_L)' (numerisch, m11; Exp-Lage, σ = 1.5)

Eingabe u = (1 − z/P)^σ + (1 − z/P̄)^σ. **Fehler sup_{|z|=1}|Π_N u − u|:**

| P | N = 400 | 600 | 800 | 1000 | 1200 |
|---|---|---|---|---|---|
| L (τ = 0, Ecke) | 6.8e-5 | 3.8e-5 | 2.5e-5 | 1.8e-5 | 1.4e-5 |
| e^{i(θ+1°)} | 5.5e-5 | 2.9e-5 | 1.8e-5 | 1.25e-5 | 9.2e-6 |
| e^{i(θ+2°)} | 4.8e-5 | 2.5e-5 | 1.5e-5 | 1.05e-5 | 7.7e-6 |
| e^{i(θ+4°)} | 4.5e-5 | 2.4e-5 | 1.5e-5 | 1.03e-5 | **1.0e-5** |
| e^{i(θ+8°)} | 4.3e-5 | 2.3e-5 | 1.4e-5 | **2.1e-5** | 3.4e-5 |
| e^{i(θ+15°)} | 4.1e-5 | 2.2e-5 | **2.9e-5** | 7.9e-5 | 1.8e-4 |
| e^{i(θ+30°)} | 3.6e-5 | **2.8e-5** | 8.2e-5 | 2.8e-4 | 1.1e-3 |
| e^{i(θ+60°)} | 3.8e-5 | 3.0e-5 | 2.0e-4 | 1.1e-3 | 3.1e-3 |
| −1 (τ = 103.4°) | 7.4e-5 | 6.9e-5 | 3.4e-4 | 1.9e-3 | 1.2e-2 |
| L(1 + 0.003) (Schlitz) | 1.9e-5 | 6.4e-6 | 2.7e-6 | 1.3e-6 | 6.5e-7 |
| L(1 + 0.01) (Schlitz) | 1.6e-6 | 2.1e-7 | 3.4e-8 | 6.4e-9 | 1.3e-9 |

(fett: Umschlag.) Das Fehlermaximum sitzt nach dem Umschlag in L oder L̄ (= Eckausgang u₁,₂). **Befund:** Singularität in der Ecke und auf dem
Schlitz ausserhalb von L: Konvergenz; auf dem D-Bogen im Winkelabstand τ: Divergenz ab N*(τ) ≈ 600 (τ ≥ 30°), 800 (15°), 1000 (8°), 1200 (4°),
> 1200 (≤ 2°). Der Umschlagpunkt wandert mit N zur Ecke ⇒ λ(τ) > 0 für jedes τ > 0, aber klein nahe L.

**Formel (D6, nur Polterm; numerisch bestätigt):** Eckamplitude a_i := ⟨u_i, Π_N u⟩ ≈ σ_i·(sin πσ/π)·[m(P)·I_i(P) + Konj.],
I_i(P) = ∫₀^∞ t^σ V_i(P(1+t)) dt = Σ_k v_{i,k}P^{−k}B(σ+1, k−σ−1). Bei P = −1: a₁ = 2.455e-6/6.98e-6/2.910e-5/1.477e-4/8.473e-4 gegen Formel
2.465e-6/7.016e-6/2.922e-5/1.482e-4/8.502e-4 (N = 400…1200; Abweichung ≤ 0.4 %). Nahe L kommt der Urbildterm V_i(h(P)) hinzu (Formel dort nur
obere Schranke). Wachstumsrate 1000→1200: ln(8.473e-4/1.477e-4)/200 = 0.00872 = λ(103°) − (σ+1)·ln(1.2)/200 = 0.0110 − 0.0023 ✓ (λ aus m07).

**Prüfstein für W (präzise Form):** Für jede Theorie des adjungierten Profils V_i (Rand und aussen):
(i) sie muss a_i(N, P) = σ_i(sin πσ/π)[m(P)I_i(P) + Urbildterm + Konj.] auf ≤ 1 % reproduzieren (P auf dem D-Bogen, auf dem Schlitz);
(ii) daraus folgt der Umschlag: Π_N u_P divergiert ⇔ λ(τ) > 0, mit N*(τ, σ) bestimmt durch λ(τ)N* ≈ (σ + 1 − κ(τ))ln N* + const;
(iii) gemessene Umschläge (σ = 1.5): N* ≈ 600/800/1000/1200 bei τ = 30/15/8/4°. Für Kneser zählt nur τ = 0 (Ecke) und der Schlitz — dort ist
die Eckamplitude polynomial bzw. fallend (Tabelle). Konsistenz mit dem Profil aus m07 (Raten 1000→1200 minus Eckrate):
λ(4°/8°/15°/30°) = 0.0022/0.0040/0.0060/0.0082 ⇒ N*·λ = 2.6/4.0/4.8/4.9 = σ + 1 − κ(τ), also κ(τ) zwischen −0.1 und −2.4 (zusätzlicher
polynomialer Abfall der Eckamplitude abseits der Ecke). Die Umschläge skalieren nur wie N* ~ τ^{−0.34} (4°–30°), nicht wie τ^{−p}: bei
endlichem N bestimmt der Vorfaktor κ(τ) den Umschlag mit — eine Profil-Theorie muss κ(τ) mitliefern.

## 10. Zwei-Punkt-Form und Ws Laurent-Form
Ws §5 (Erzeugende bei ∞, Pole nur 0, ã₀) und meine §4.2 (Entwicklung bei 0 und a₀ mit Gewichten (t−a₀)^{N+1} und Λ^{N+1}) sind dieselbe Gleichung;
P hat §4.2 an f_K nachgerechnet (1e-68, p19). Mitgeteilt (BRETT 01:04).

## 11. Runde 3 — Injektivität und (a″)_∞ im radialen Modell: rechnergestützter Beweis über den Randoperator 𝒯

**Ergebnis (rechnergestützt bewiesen; Prüfung durch P steht aus).** Exp-Lage θ = 76.618°, ε = 0.31382.
- **Satz M6.** Für alle x ∈ L²(∂D) gilt Re⟨x + K_P x, x⟩ ≥ **0.8211**‖x‖² (K_P siehe 11.1). Insbesondere Re⟨𝒯x, x⟩ ≥ 0.8211‖x‖² auf H²(D).
- **Korollar M11.** T ist auf E²(Ω) injektiv.
- **Korollar M12 ((a″)_∞).** ‖J_∞x‖ ≥ **0.4249**‖x‖ für jedes Polynom x, also σ_min(J_∞|P_N) ≥ 0.4249 für alle N (numerisch 0.6905).
- Mit Os Satz O7 (Fredholm, Index 0) folgt: T: E²(Ω) → H²(D) ist invertierbar. Zusammen mit O7 ist das der erste vollständige
  Bulk-Baustein in einem Modell mit Eckmode; für (a″)_∞ selbst wird O7 nicht gebraucht.
- Der Vorschlag aus Runde 2 (Regularität der Kernelemente mit Exponent 6.29, dann Widerspruch zu einer zertifizierten Schranke für
  σ_min(J_∞|P_{N₀})) ist damit überholt: die Akkretivität liefert Injektivität UND (a″)_∞ ohne Regularitätstheorie.

### 11.1 Der Randoperator 𝒯 (bewiesen)

h := Zweig von f⁻¹, h(w) = (−c₀ + √(c₀² − 4ε(ε − w)))/(2ε) (Hauptzweig; h(L) = L). Es gilt f(Γin) ⊂ D (|f(e^{iφ})| = s(φ) < 1) und
h(Γout) ⊂ D (in Streifenkoordinaten Im H ∈ [2.5926, 2.6344] ⊂ (0, π), m20 mit Einschluss). Für x ∈ L²(∂D) sei
  (K_P x)(φ) := W(φ)·(P_D x)(π(φ)),  π = f, W = 1/m auf Γin;  π = h, W = f' auf Γout  (P_D = Poisson-Fortsetzung),
und 𝒯 := (I + K_P)|_{H²(D)}, also 𝒯x = x + (x∘f)/m auf Γin, x + f'·(x∘h) auf Γout — es werden nur Werte in D̄ benutzt.
**Lemma M5.** (i) Für y ∈ E²(Ω) mit y|_D ∈ H²(D) ist 𝒯(y|_D) = Ty auf Γin und = f'·(Ty)∘h auf Γout; also ker 𝒯 ∩ (E²(Ω)|_D) = ker T.
(ii) Für Polynome x: 𝒯x = J_∞x/m auf Γin und 𝒯x = (J_∞x)∘h auf Γout (wegen m∘h = f'∘f∘h = f').
*Beweis.* (Tx)(v) = x(v) + x(f(v))/m(v); bei v = h(ζ), ζ ∈ Γout: f(v) = ζ und m(v) = f'(ζ). ∎

### 11.2 Streifenkoordinaten (Os ζ = log M) und unitäres Bild

M(x) = (x − L)/(x − L̄), ζ = log M − iθ bildet D konform auf S = {0 < Im ζ < π} ab: Γin → Unterkante, Γout → Oberkante, L ↔ Re ζ = −∞,
L̄ ↔ +∞; Umkehrung ξ(ζ) = L̄ + (L̄ − L)/(e^{iθ}e^ζ − 1). Wegen f(x) − L = (x − L)(c₀ + ε(x + L)) ist exakt
  G(t) := ξ⁻¹(f(ξ(t))) = t + Λ(ξ(t)),  H(t) := ξ⁻¹(h(ξ(t+iπ))) = t + iπ − Λ(h(ξ(t+iπ))),  Λ(x) = log[(c₀+ε(x+L))/(c₀+ε(x+L̄))],
G − t → ℓ = log μ (L-Ende), → −ℓ̄ (L̄-Ende). Kantenmasse dφ/2π = w_e(t)dt mit w₀ = sinθ/(2π(cosh t − cosθ)), w₁ = sinθ/(2π(cosh t + cosθ)).
Unitär Y_e(t) := x(ξ(t + iπe))√w_e(t): K_P wird zum Integraloperator auf L²(ℝ)² mit Kern
  K̃_{ee'}(t,t') = √(w_e(t)/w_{e'}(t'))·W_e(t)·P_{e'}(Φ_e(t), t'),  Φ₀ = G, Φ₁ = H,
  P₀(u+iy, t') = sin y/(2π(cosh(t'−u) − cos y)),  P₁(u+iy, t') = sin y/(2π(cosh(t'−u) + cos y))  (Streifen-Poissonkern).
Alle Innenauswertungen liegen gleichmässig in Abstand ≥ arg μ = 0.548 von der jeweils massgeblichen Kante (Im G ∈ [0.548, 0.923],
Im H ∈ [2.593, 2.634]) — auch in den Ecken. Das ist der Grund, warum alles gleichmässig abschätzbar wird.
Numerisch (m15, Nyström auf [−20, 20]): λ_min(I + Re K̃) = 0.990086 (Eckwert 0.990018).

### 11.3 Beweis von Satz M6

**IMS-Zerlegung.** χ_L = cos(πs/2), χ_m = sin(πs/2) mit s = (t + T₂)/W auf [−T₂, −T₁], gespiegelt rechts (χ_R), χ_m = 1 auf [−T₁, T₁];
T₁ = 5, T₂ = 13, W = 8, Σχ_j² = 1. Für jeden beschränkten Operator gilt K̃ = Σ_j χ_jK̃χ_j + R, R = ½Σ_j[χ_j,[χ_j,K̃]], Kern
R(t,t') = ½Σ_j(χ_j(t) − χ_j(t'))²K̃(t,t'). Also Re⟨(I + K̃)Y, Y⟩ ≥ (min_j c_j − ‖R‖)‖Y‖², c_j := inf über Y mit Träger in supp χ_j.

**(A) Lemma M7 — Eckmodell (bewiesen, arb; m16_mellin_zert.py, m16_zert.log).** Das Modell am L-Ende,
(AY)₀(t) = μ⁻¹∫e^{(t−t')/2}P_{e'}(t + ℓ, t')Y_{e'}(t')dt', (AY)₁(t) = μ∫e^{(t−t')/2}P_{e'}(t + iπ − ℓ, t')Y_{e'}(t')dt', ist eine Faltung;
mit ∫e^{κs}/(cosh s − cos y)ds = 2π sin(κ(π−y))/(sin y·sin κπ) (|Re κ| < 1) ist ihr Symbol (s = −½ + iσ, b = arg μ)
  A(σ) = [[μ⁻¹|μ|^s F₀(b), μ⁻¹|μ|^s F₁(b)], [μ|μ|^{−s} F₀(π−b), μ|μ|^{−s} F₁(π−b)]],  F₀(y) = sin(s(π−y))/sin sπ, F₁(y) = sin(sy)/sin sπ
— dasselbe Symbol wie im Sektor (Runde 3, m13). **c_mod := inf_σ λ_min(I + Re A(σ)) ≥ 0.99** (Intervall-Bisektion auf [−30, 30], Schwanz
‖A(σ)‖_F ≤ 2√2 e^{−|σ|b}(|μ|^{−3} + |μ|³)^{½} ≤ 3.1e-7); Wert 0.9900183 bei σ = −6.50. Für μ = √L: ≥ 0.9898 (Wert 0.9898832 = Os λ_L).

**(B) Lemma M8 — Endstück (bewiesen, arb; m17_ende.py, m17_ende.log).** Für Y mit Träger in (−∞, −T₁] ist |⟨(K̃ − A)Y, Y⟩| ≤ ε_L‖Y‖²,
**ε_L ≤ 0.0439** (T₁ = 5; 0.0159 bei T₁ = 6). Beweis: g₁ = Λ(ξ(M)) − ℓ, g₂ = 1/m − 1/μ, g₃ = Λ(h(ξ(M))) − ℓ, g₄ = f'(ξ(M)) − μ sind auf
|M| ≤ 0.3 analytisch (Kästchenüberdeckung: Log- und Wurzelargumente meiden (−∞, 0], f'∘f ≠ 0) und verschwinden in 0; Schwarz-Lemma
|g_i| ≤ γ_i e^t mit γ = (0.4795, 2.0538, 0.3453, 1.7470). Mittelwertsatz für P (|∇P₀| = 1/(2π(cosh α − cos y))), Gewichtsfaktor
√(w_e(t)/w_{e'}(t')) = e^{(t−t')/2}ρ_e(t)/ρ_{e'}(t') mit ρ = 1/|1 ∓ e^{iθ}e^t|; Schur-Test mit gleichmässigen Majoranten in Δ = t' − t und den
geschlossenen Integralen ∫e^{∓α/2}/(cosh α − cos y)dα = π/sin(y/2). Das L̄-Ende folgt aus der Symmetrie (JY)(t) = conj Y(−t), K̃J = JK̃
(reelle Koeffizienten: ξ(−ζ̄) = conj ξ(ζ), G(−t) = −conj G(t)); also c_R = c_L ≥ c_mod − ε_L ≥ 0.9461.

**(C) Lemma M9 — Kommutator (bewiesen, arb; m18_kommutator.py, m18_T5_13.log).** **‖R‖ ≤ 0.1250.** Schur-Test für |R| auf einem Raster
der Weite 0.1 auf [−32, 32] (Kantengrössen per Unterteilung in 24 Teilkugeln eingeschlossen), S := ½Σ(χ_j(t) − χ_j(t'))² ≤ min(1, π²Δ²/8W²)
und S = 0 im selben flachen Bereich; Paare ausserhalb ⇒ |Δ| ≥ 19, Majorante |K̃| ≤ 1.61·3·(2e/π)e^{−|Δ|/2} (Voraussetzungen |Re Φ − t| ≤ 1,
|W| ≤ 3 auf dem Raster geprüft, ausserhalb aus den Schwarz-Schranken) ⇒ Schwanz ≤ 5.0e-3. Zeilensumme ≤ 0.1349, Spaltensumme ≤ 0.1066.
(Numerisch ist ‖R‖ = 0.048 und λ_min(Re R) = −0.018 — der Schur-Schätzer ist der grösste Verlust im Beweis.)

**(D) Lemma M10 — Mittelstück (bewiesen, arb + verifizierte Gleitkomma-Cholesky; m19_mitte.py, m19_T13.log).**
**λ_min(I + Re K̃) auf L²([−13, 13])² ≥ 0.99004.** 52 Paneele der Länge ½ je Kante, je 16 Gauss-Legendre-Knoten (arb.legendre_p_root).
Endlicher Rang k_N(t,t') = Σ k(t_i,t'_j)ℓ_i(t)ℓ_j(t'); in der Orthonormalbasis ℓ_i/√ω_i (exakt orthogonal) ist K_N die Matrix √ω_i k(t_i,t'_j)√ω_j
(1664 × 1664, Einträge in arb, ‖B − B_float‖ ≤ 3.8e-13). λ_min ≥ s − γ‖R̂‖²_F, wenn die Gleitkomma-Cholesky von B − sI durchläuft
(Higham Satz 10.3, komplexer Faktor 4): ≥ 0.99023452. Interpolationsfehler: je Paneelpaar sup|k − k_N| ≤ C_int(M_{t'} + Λ₁₆M_t) (Hermite-
Restformel auf der Bernstein-Ellipse ρ = 2.8, C_int ≤ 1.63e-6, Λ₁₆ ≤ 6.91), M_{t'} = sup|k| für t' auf der Ellipse, M_t = sup|k^c| für t auf
der Ellipse mit der analytischen Fortsetzung P₀^c = (1/2πi)[φ(Φ(t)−t') − φ(Φ*(t)−t')], φ(z) = 1/(1−e^z), Φ*(t) = conj Φ(conj t)
(Fernpaare mit reeller Intervallrechnung); ‖k − k_N‖_HS ≤ 1.9e-4. Voraussetzung der Restformel — Analytizität im INNERN der Ellipsen — geprüft: in t per Kästchenüberdeckung von {|Im t| ≤ 0.306} (Log-/Wurzelargumente meiden (−∞,0], f'∘f ≠ 0, 0 < Im G < 2π, 0 < Im H < π, ebenso für conj t), in t' analytisch (Pole von P bei |Im t'| = y bzw. π − y ≥ 0.488 > Halbachse 0.305). Analytizitätsbreite des Kerns ≈ 0.4 (der kritische Punkt von G bei
t = −0.24i ist keine Singularität; Im G ≥ 0.31 für |Im t| ≤ 0.4).

**Zusammen:** c ≥ min(c_mod − ε_L, c_m) − ‖R‖ ≥ min(0.9461, 0.9900) − 0.1250 = **0.8211**. ∎

### 11.4 Folgerungen

**Korollar M11 (Injektivität; bewiesen | Satz M6 + Einbettungslemma).** Ist y ∈ E²(Ω) mit Ty = 0, so ist 𝒯(y|_D) = 0 nach Lemma M5,
also 0 = Re⟨𝒯y, y⟩ ≥ 0.8211‖y|_D‖², y|_D = 0, y = 0. *Einbettungslemma (Skizze, Standard):* E²(Ω)|_D ⊂ H²(D) stetig — auf Γin ist ∂D = ∂Ω,
und Γout liegt im Innern von Ω und trifft ∂Ω nur in L, L̄, dort nicht tangential (Γout ist im Eckensektor von Ω, Öffnung π + arg μ, ein innerer
Strahl); Bogenlänge auf Γout ist daher ein Carleson-Mass für E²(Ω) (z. B. nach konformer Abbildung von Ω auf einen Streifen wie in 11.2).
**Korollar M12 ((a″)_∞; bewiesen | Satz M6, Lemma M13).** Für Polynome x: 0.8211‖x‖ ≤ ‖𝒯x‖ ≤ C_E‖J_∞x‖ mit C_E² = 1/m_in² + C_h²
(Lemma M5 (ii)), also **σ_min(J_∞|P_N) ≥ 0.8211/1.9324 = 0.4249** für alle N. **Verfeinert (Runde 4, m20): ≥ 0.5596** — die Kanten getrennt: ‖𝒯x‖² ≤ ‖Z₀‖²/m_in² + R∫|Z(s+iy(s))|²ds, R = r_max/d₀, und ‖Z(·+iy₁)‖² ≤ α‖Z₀‖² + β‖Z₁‖² mit α = (1−θ)s^θ, β = θs^{θ−1}, θ = y₁/π (Young für e^{−2σy₁} ≤ α + βe^{−2πσ}); ⇒ C_E'² = max(1/m_in² + R[(1+ε)α + κ'], R[(1+ε)β + κ']) ≤ 1.4672².
**Lemma M13 (Carleson-Konstante; bewiesen, arb; m20_carleson.py, m20_carleson.log).** ∫_Γout|g∘h|²dφ/2π ≤ C_h²‖g‖²_{H²(D)} mit
**C_h ≤ 1.2794**; m_in = min_Γin|m| ≥ 0.690509 (zertifiziert). Beweis im Streifen-Hardyraum: Z = g∘ξ·(ξ'/2π)^{½}, ∫_Γout|g∘h|² = ∫|Z(H(t))|²r(t)dt,
r = |ξ'(t+iπ)|/|ξ'(H)| ≤ 1.5404, Re H' ≥ 0.9991 (Umparametrisierung s = Re H), Höhe Im H ∈ [2.5926, 2.6344]; |Z(s+iy)|² ≤ (1+ε)|Z(s+iy₁)|² +
(1+1/ε)(y₂−y₁)∫|Z'(s+iη)|²dη, ‖Z(·+iy)‖ ≤ ‖Z‖, ‖Z'(·+iη)‖² ≤ κ‖Z‖², κ ≤ 0.526 (Fourier). Numerisch C_h = 1.052.

### 11.5 Übertragung auf Kneser; kann der Weg (G2) ersetzen?

| Schritt | radiales Modell | Kneser | was fehlt für Kneser |
|---|---|---|---|
| 𝒯, Lemma M5 | bewiesen | gleich, sobald f_K(Γin) ⊂ D, h_K = f_K⁻¹ auf Γout mit h_K(Γout) ⊂ D | Os (S1)–(S3) (bisher numerisch, P p01) |
| Streifen ζ = log M | exakt | exakt (nur L, L̄ ∈ ∂D nötig) | — |
| Eckmodell (A) | ≥ 0.99 | **≥ 0.9898 zertifiziert** (μ = √L, m16) | — |
| Endstück (B) | g_i analytisch in M | f_K ist in L nicht analytisch; im Sektor f_K = reguläre Halbiterierte + O(|z−L|^{1+4.45}) (Knesers ϑ fällt wie e^{2πiu}) ⇒ G_K = t + ℓ + O(e^t) | Schranken für Schröder-Funktion (analytisch, zertifizierbar) und für ϑ (Fourier-Koeffizienten, Strang A) |
| Kommutator (C), Mitte (D), C_h | arb-Programme m18–m20 | dieselben Programme | **zertifizierte Einschlüsse von f_K, f_K', h_K auf den Streifenkanten und auf Ellipsen |Im t| ≤ 0.3** (Motor nur über P) |
| Ergebnis | Injektivität + (a″)_∞ ≥ 0.42 | Injektivität + (a″)_∞ mit expliziter Konstante | — |

**Antwort auf die Frage der Leitung:** Ja — der 𝒯-Weg **ersetzt (G2)**. D4 braucht die globale Abel-Koordinate (slog), die Geometrie von
slog(D) (= (G2), nur numerisch) und eine Wachstumsklasse und liefert nur Injektivität. Der 𝒯-Weg braucht keine Abel-Koordinate und keine
Wachstumsklasse (er arbeitet in H²(D)/E²(Ω)), sondern: (S1)–(S3) (die braucht O7 ohnehin), das Eckverhalten von f_K in L (ϑ-Schranke) und
eine endliche zertifizierte Rechnung mit f_K-Einschlüssen nahe ∂D — und er liefert zusätzlich (a″)_∞ mit Konstante. Numerischer Hinweis:
Os λ_N für Kneser (0.9954 bei N = 100) liegt wie im Modell knapp über dem Eckwert; die Margen sind gleich gross.
**Gegenüber Os (FA):** gleiche Eckkonstante (die Infima stimmen überein, obwohl die σ-Funktionen verschieden sind); 𝒯 braucht keine
Begradigung von f(Γout) (kein Theodorsen), weil alle Daten auf den geraden Kanten von D liegen und nur Innenauswertungen vorkommen.

### 11.6 Prüfaufträge (→ P) und offene Punkte

1. Lemma M5; IMS-Identität und Symmetrie J. 2. m16: Symbol (Fourier-Formel, Vorzeichen, Gewichte μ^{∓1}). 3. m17: Analytizität auf |M| ≤ 0.3,
Schwarz-Lemma, Majoranten (Funktion J), Faktor d_in/d_out. 4. m18: Schur-Summen, S-Schranke, Schwanzmajorante. 5. m19: Orthonormalität,
Hermite-Restformel (C_int, Λ₁₆), analytische Fortsetzung in t (P₀^c, P₁^c mit φ₁ = 1/(1+e^z)), Cholesky-Schranke, Fernpaare. 6. m20.
7. Weiter offen aus Runde 2: Lemma M2/M3 (§1, §8.3).
Offen (nicht nötig für M6/M11/M12): das Einbettungslemma ausgeschrieben; schärfere Konstanten (gewichteter Schur-Test für R ⇒ c ≈ 0.95).

## 13. Runde 4 — der 𝒯-Weg für Kneser: Vorzertifikat, Datenschnittstelle zu P, Vergleich mit O11

**Ergebnis (Vorzertifikat — streng relativ zu den Strang-A-Koeffizienten und einer angenommenen Schwanzschranke; wird mit Ps Daten zum Beweis):**
für f_K gilt **Re⟨𝒯x, x⟩ ≥ 0.8062‖x‖²** auf H²(D) und **σ_min(J_∞|P_N) ≥ 0.5645 für alle N** (wahr 0.7319), ohne (G2), ohne Abel-Koordinate,
ohne Univalenz von f_K. Gesamtlauf: `python m26_kneser_gesamt.py [quelle]` (≈ 20 min; Teilskripte m22–m25 + m16).

### 13.1 Was aus dem Modell übernommen wird
Lemma M5 (für Polynome braucht es nur F(Γin) ⊂ D, h(Γout) ⊂ D und F(h(ζ)) = ζ — beides wird unten mit Einschluss gezeigt; Univalenz
nicht), Streifen ζ = log M (θ = arg L ist im Modell und bei Kneser gleich, also dieselben Gewichte w_e), IMS-Zerlegung, Eckmodell
(m16 'sqrtL': **c_mod ≥ 0.9898**, zertifiziert; Wert 0.9898832 = Os λ_L), Schur- und Nyström-Programme.

### 13.2 Datenschnittstelle (mk_kneser_arb.py): alles aus Koeffizienten, streng
F(z) := f_K(ρz)/ρ = (1/ρ)Σ b_k z^k, b_k = a_kρ^k, z ∈ D̄ (ρ = |L|, L = conj(−W₀(−1)) in arb). Aus
 **(K1)** Kugeln für b_k, k ≤ K, **(K2)** |b_k| ≤ C·k^{−p} für k > K (p > 3), **(K3)** f_K(L) = L, f_K'(L) = √L (Randwerte; Kneser/Strang A)
und der Funktionalgleichung f_K∘f_K = exp (Definition) folgen streng:
- F, F' an Punkten und Kugeln in D̄ und auf Kreisbögen: Horner mit EXAKTEN Mitten bei 400 Bit (die Rechteck-Umwicklung ≤ √2 je Schritt
  wirkt nur auf Rundungsradien), Koeffizientenradien und Restsummen getrennt addiert, Mittelwertsatz mit S1 = sup|F'| ≤ 1.7166,
  S2 = sup|F''| ≤ 1.0459 (aus Σk^j|b_k|);
- h = F⁻¹ auf Γout und auf Kästchen beiderseits von Γout: parametrischer Krawczyk in D̄ (Existenz, Eindeutigkeit im Kästchen, F' ≠ 0);
- Fortsetzung nahe Γin auch AUSSERHALB D̄: F(z) = log(ρF(w))/ρ, F'(z) = F'(w)w/F(w), w = e^{ρz}/ρ ∈ D; W_in = 1/m = F'(z)e^{−ρz};
  W_out = F'(y) = e^{ρh(y)}/F'(h(y)) — **keine Werte von f_K ausserhalb D̄ aus anderen Quellen nötig**;
- Ecken (mk_ecke_kneser.py): f_K ist in L nicht analytisch (Terme t^{1+nγ}), aber Taylor zweiter Ordnung längs der Sehne mit S2 genügt:
  |G − t − ℓ|, |W_in − 1/μ|, |H − (t+iπ−ℓ)|, |W_out − μ| ≤ explizite O(e^t); h nahe L = einziger Fixpunkt einer Kontraktion in einer Kugel ⊂ D̄.

### 13.3 Vorzertifikat (Daten Strang A; Radien 1e-15·|b_k|; Schwanz C = 100, p = 6.4 für k > 4000 ANGENOMMEN)
| Baustein | Skript | radial (streng, M6) | Kneser (Vorzertifikat) | Kneser numerisch (m21) |
|---|---|---|---|---|
| Eckmodell c_mod | m16 | ≥ 0.99 | ≥ 0.9898 (streng) | 0.98988 |
| Endstück ε_L | m22 | 0.0439 (T₁ = 5) | 0.0351 (T₁ = 6) | — |
| Kommutator ‖R‖ | m23 | 0.1250 (W = 8) | 0.1485 (W = 8) | 0.055 (‖Re R‖) |
| Mitte c_m | m24 | ≥ 0.99004 | ≥ 0.98995 | 0.99017 |
| ⇒ Akkretivität c | — | 0.8211 | 0.8062 | 0.98998 (λ_min gesamt) |
| C_E' (Kanten getrennt) | m20/m25 | 1.4672 | 1.4281 | — |
| ⇒ σ_min(J_∞|P_N) ≥ | — | 0.5596 | 0.5645 | wahr 0.6905 / 0.7319 |
Weitere Kneser-Grössen (Einschluss): F(Γin) ⊂ D (|F| ≤ 0.99927), h(Γout) ⊂ D, Im G ∈ [0.669, 0.931], Im H ∈ [2.456, 2.537],
min_Γin|m| ≥ 0.731378 (wahr 0.731905), C_h ≤ 1.357 (numerisch 1.048).
**Vergleich mit O11 (Konturverschiebung, Gleitkomma-Schur, BRETT 03:05):** radial ≥ 0.408, kubisch ≥ 0.520, Kneser ≥ 0.538. Gleiche Eckkonstante
λ_L = 0.98988 in beiden Wegen; der 𝒯-Weg liegt jetzt mit 0.5645 darüber und ist (bis auf die Daten) streng.


**Technische Anmerkungen (für P).** (i) Komplexe Kästchen in m24: Mittelwertform Φ(t) ∈ Φ(t₀) + |t − t₀|·sup|Φ'| (Ableitungen über
dieselben Fortsetzungsformeln, F'' ∈ Kreis(S2) auf D̄) — die naive Kugelrechnung überschätzt um Faktor 10–30. (ii) Matrixeinträge mit
cosh(t' − u) direkt (die Produktform löscht bei |t|, |t'| ≈ 14 aus). (iii) Reelle Teilkugeln nahe der Mitte: Vereinigung über 16 Teile.
(iv) Analytizität im Innern der Ellipsen: Kästchen 0.03, adaptiv bis Tiefe 3 unterteilt, 0 kritische. (v) Konsistenz der Zweige von h:
Krawczyk-Kästchen (|t| ≤ 13) liegen in der Fixpunktkugel der Ecke (|t| ≥ 6) — Eindeutigkeit dort ⇒ derselbe Zweig.
**Gegenprobe numerisch (m21, Gleitkomma):** Kneser λ_min(I + Re K̃) = 0.98998, Mitte 0.99017, ‖Re R‖ = 0.055, C_h = 1.048 — die Vorzertifikats-
zahlen liegen konsistent darunter. **Hinweis an P (Os 7. Stelle):** λ_L = 0.98988322 (m16, Raster 0.01), P: 0.98988321 — O's 0.9898822 ist ein Tippfehler.

### 13.4 Was Strang P liefern muss (→ BRETT)
Datei in gemeinsam/ im Format {"K": int, "b_mitte": [str], "b_radius": [str], "C": str, "p": str}:
- **(K1)** b_k = a_kρ^k für 0 ≤ k ≤ K, K ≥ 400 (besser 1000), als Kugeln mit Σ_{k≤400} Radius ≤ 1e-12 (z. B. jeder Radius ≤ 1e-15);
- **(K2)** C, p mit |b_k| ≤ C·k^{−p} für ALLE k > K, p > 3 (beobachtet p = 6.45, C ≈ 74); Bedarf: Σ_{k>400}|b_k| ≤ 1e-12,
  Σ_{k>400} k|b_k| ≤ 1e-9, Σ k²|b_k| endlich (S2);
- **(K3)** Bestätigung (mit Fundstelle): f_K(L) = L und f_K'(L) = √L als Randwerte in D̄.
Dann `python m26_kneser_gesamt.py gemeinsam/<datei>.json` — keine weitere Arbeit nötig. (Dieselben Daten genügen für Os O11.)

### 13.5 Prüfaufträge (→ P)
Runde 3: Lemma M5, M7–M10, M13 (m16–m20). Runde 4: mk_kneser_arb (Horner/Mittelwertsatz/Krawczyk/Fortsetzung), mk_ecke_kneser
(Fixpunkt- und Sehnenargument), m22–m26; Verfeinerung C_E' (Young-Ungleichung, §11.4).

## 15. Runde 5 — Kneser mit Ps zertifizierten Koeffizienten; die Kante (a″)_{αN}; Gegenprobe zu Satz P-V

### 15.1 (a″)_∞ für Kneser — streng bis auf die Schwanzschranke (K2)
Mit Ps Datei (a_0…a_640, Satz P-V; mk_kneser_arb Quelle "pzert") läuft m26 unverändert durch (m26_kneser_pzert.log):
c_mod ≥ 0.9898, ε_L ≤ 0.0351, ‖R‖ ≤ 0.1485, c_m ≥ 0.98996 ⇒ **Re⟨𝒯x,x⟩ ≥ 0.8062‖x‖², σ_min(J_∞|P_N) ≥ 0.5646 für alle N** (wahr 0.7319).
Einzige offene Voraussetzungen: **(K2)** |b_k| ≤ C k^{−p} für k > 640 (Lauf mit C = 100, p = 6.4) und **(K3)** f_K(L) = L, f_K'(L) = √L. (K3) folgt aus P-V selbst (Skizze): oberhalb Im z ≥ 0.35 ist sexp_K = Ũ = P(L^{z+θ̃(z)}) mit trigonometrischem θ̃ (Moden e^{2πikz}, k ≥ 1, klingen ab) ⇒ sexp_K(w) → L und (sexp_K(w + ½) − L)/(sexp_K(w) − L) → L^{1/2} für Im w → ∞; mit f_K = sexp_K(slog + ½) also f_K(L) = L, f_K'(L) = √L als Randwerte — P bitte bestätigen.
**Was (K2) leisten muss (aus den Skripten):** T₀ := Σ_{k>640}|b_k| ≤ 1e-9 (Genauigkeit von F − L an den Fensterenden |δ| ≈ 1e-6),
T₁ := Σ_{k>640} k|b_k| ≤ 1e-6, T₂ := Σ_{k>640} k²|b_k| ≤ 1e-3 (S2 = sup|F''| ändert sich dann um < 1e-3). Zum Beispiel genügt
p = 5 mit C ≤ 600, p = 6 mit C ≤ 5e5; die Daten zeigen C ≈ 74, p = 6.45 (T₀ ≈ 7e-15) — **Reserve sechs Grössenordnungen**.
Gegenprobe mit absichtlich schlechtem Schwanz (C = 500, p = 5; T₀ = 7.5e-10, T₁ = 6.4e-7): c_m ≥ 0.98666 (statt 0.98996), Ergebnis unverändert c ≥ 0.8062, σ_min ≥ 0.5645 (m26_kneser_pzert_500_5.log).
**Breiterer Übergang (gleiche Skripte, nur T₁, T₂):**
| T₁ / T₂ (W) | ε_L | ‖R‖ | c_m | c | σ_min(J_∞|P_N) ≥ | Log |
|---|---|---|---|---|---|---|
| 6 / 14 (8) | 0.0351 | 0.1485 | 0.98996 | 0.8062 | 0.5646 | m26_kneser_pzert.log |
| 6 / 18 (12) | 0.0351 | 0.0705 | 0.98981 | 0.8842 | 0.6192 | m26_kneser_pzert_T6_18.log |
| **7 / 23 (16)** | 0.0126 | 0.0414 | 0.98740 | **0.9359** | **0.6553** | m26_kneser_pzert_T7_23.log |
(Mitte ≈ 9 min, Gesamtlauf ≈ 11 min.) Der Rest zum wahren Wert 0.7319 steckt fast ganz in C_E' = 1.428 (wirksam ≈ 0.99/0.7319 = 1.35).
**Für X (KONSTANTEN.md):** **c₀ = 0.6553** (Kneser, streng bis (K2)/(K3)); radial **0.6395** (streng; T₁ = 6, T₂ = 22: ε_L 0.0159, ‖R‖ 0.0358, Mitte 0.98983 ⇒ c = 0.9383, C_E' = 1.4672; m18_T6_22.log, m19_T22.log).
**P (04:19):** liefert (K2) mit p = 6, C₆ ≤ 300 (Σ_{k>640}|b_k| ≤ 6e-13, Σ k|b_k| ≤ 5e-10, Σ k²|b_k| ≤ 4e-7 — alle Grenzen erfüllt) und (K3) über p32
im Format {"K","b_mitte","b_radius","C","p"}; mk_kneser_arb liest es (Dezimalrundung der Mitten und Aufrundung der Radien eingebaut).

### 15.2 Die Kante (a″)_{αN}: was der 𝒯-Weg liefert und was nicht
**Lemma M14 (trivial, bewiesen).** σ_min(T_M J_∞|P_N) ≥ σ_min(J_∞|P_N) − ‖Q_M J_∞ P_N‖.
**Wann der Überlauf klein ist.** [x∘f]_k = Σ_j x_j[F^j]_k, Cauchy auf |z| = R (F analytisch bis auf Schlitze an L, L̄):
‖Q_{(1+α)N}(x∘F)‖ ≲ (sup_{|z|=R}|F|·R^{−(1+α)})^N‖x‖. Exponentiell klein nur, wenn sup_{|z|=R}|F| < R^{1+α} für ein zulässiges R:
- radial: f quadratisch ⇒ J_∞P_N ⊂ P_{2N}: **(a″)_{2N} = (a″)_∞, σ_min(T_{2N}J_∞|P_N) ≥ 0.5596 bewiesen** (α = 1, ohne Verlust).
- Kneser: sup_{|z|=R}|F| = F(R) (positive Achse); R = 1.2 (Holomorphie bis |x| = 1.66): F(1.2) = e/ρ = 1.977 ⇒ α > log 1.977/log 1.2 − 1 = **2.74**
  — das ist Xs Lemma X5 ((a″)_{4N} ⇐ (a″)_∞); mit c₀ = 0.5646 also σ_min(T_{4N}J_∞|P_N) ≥ 0.5646 − (Xs Schwanz, O(N^{−4.45})).
**Für α ∈ [0.03, 0.1] ist das unmöglich:** auf dem Rand ist max|F| = 1.624 (Kneser, φ = 0, lokaler Grad zF'/F = 1.041), 1.483 (radial) — der
Überlaufblock ist exponentiell gross. Gemessen radial (m28, arb): ‖Q_{(1+α)N}J_∞P_N‖ (grösste Spaltennorm) = 2.1e9 / 1.2e16 bei N = 60/100
(α = 0.03) ⇒ Rate 1.48^N. Gleichzeitig σ_min(T_{(1+α)N}J_∞|P_N) = 0.6911/0.6907 (N = 60/100, α = 0.03…1) und nach m10 (N = 400/600)
schon ab d = 4/8 Zusatzzeilen der volle Wert 0.6905. ⇒ (a″)_{αN} mit kleinem α ist **keine Störung von (a″)_∞**: die ersten αN Überlauf-
zeilen müssen genau die Eckrichtungen erfassen, während der Rest des Überlaufs exponentiell gross bleibt. Das ist eine Aussage vom Typ (K_d)
(X: Gram-Bedingung der Überlaufprofile; O §32: Innenschicht; W: Kegel) — mit dem 𝒯-Werkzeug allein nicht zu erreichen.
**Was 𝒯 dazu beitragen kann (Vorschlag, nicht ausgeführt):** die Koerzivität zerlegt P_N in Eckrichtungen E (wenige, r ≈ N/100) und den Rest;
auf dem Rest liefert schon der quadratische Abschnitt die Schranke ((a′)-Typ), auf E braucht man, dass die Zeilen N+1…N+d die Eckrichtungen
treffen — genau Xs S_N (d×d) und Os Kantenraten 32.3. Ein Beweis braucht die lokale Parametrix an L (O/W).

### 15.3 Gegenprobe zu Satz P-V (m27_pv_gegenprobe.py, m27_pv_gegenprobe.log)
1. Ps Koeffizienten gegen Strang As unabhängige Datei: max |a_n^P − a_n^A|ρ^n = 1.1e-43 (n = 613) ✓ (Leitung: 1.1e-43).
2. Funktionalgleichung mit Ps Kugeln (arb, 500 Bit; Schwanz n > 640 über Ps Cauchy-Schranke |a_n| ≤ 2.128·1.30^{−n}, die mit
   max_n |a_n|1.30^n = 1.139 konsistent ist): |f_K(f_K(x)) − e^x| ≤ 2.4e-138 (|x| = 0.3), 3.4e-72 (0.5), 3.3e-41 (0.6) ✓.
3. **Befund (an P):** mit Ps Radien e_n schliesst der Einschluss die 0 NICHT ein (Residuum ≈ −4e-141, Radius 3e-142 bei |x| = 0.3) — Ursache:
   die Dezimalmitten haben 140 Stellen (p29: `mid().str(140, radius=False)`), ihre Rundung (bis 5e-141 bei a_0) fehlt in e_n (e_0 = 8.8e-143);
   ausserdem werden die Radien mit 5 Stellen ausgegeben und dabei teils ABgerundet (bis 4.5e-5 relativ, Spalte 3 vs. Spalte 2). Mit e_n + 10^{−139}|a_n|
   schliesst der Einschluss die 0 überall ein. Für alle Anwendungen bedeutungslos (gewichtet ≤ 1e-140 ≪ 2.9e-127), für die Strenge der Datei aber
   zu korrigieren: Mitten mit mehr Stellen schreiben oder die Dezimalrundung addieren; Radien aufrunden. mk_kneser_arb addiert sie bereits.
4. Reihe bis 640 bei L: |f_K(L) − L| ≈ 3e-15, |f_K'(L) − √L| ≈ 2e-12 — verträglich mit (K3) bei angenommenem Schwanz (kein Beweis).
Noch nicht gegengelesen: die Beweisschritte 3–6 (p27/p28: Taylor-Modelle, Beltrami/Stoilow, Paulsen–Cowgill-Prüfungen) — Ps Prüfliste abwarten.

## 17. Stand beim Anhalten (11.10.2026, 04:50) und Wiederaufnahme

**Bewiesen (rechnergestützt, arb; Prüfung durch P noch offen):**
- Radiales Modell: Satz M6 (𝒯 akkretiv), Korollar M11 (T injektiv auf E²(Ω); Einbettungslemma nur Skizze), Korollar M12 ((a″)_∞):
  σ_min(J_∞|P_N) ≥ 0.4249 → 0.5596 (Kanten getrennt) → **0.6395** (breiterer Übergang T₁ = 6, T₂ = 22). Damit auch (a″)_{2N} (J_∞P_N ⊂ P_{2N}).
- Eckkonstante: radial ≥ 0.99, Kneser (μ = √L) ≥ 0.9898 (λ_L = 0.98988322).
- Kneser: **σ_min(J_∞|P_N) ≥ 0.6553 für alle N** (T₁ = 7, T₂ = 23), streng mit Ps zertifizierten a_0…a_640 — bis auf (K2) und (K3).
**Von P geprüft:** Lemma M2, Lemma M3 (bis Gitter F = 0.8751). **Noch nicht geprüft:** M5, M7–M10, M13, mk_kneser_arb, mk_ecke_kneser, m22–m26,
die Kanten-Verfeinerung von C_E' (§11.4) und das Einbettungslemma.
**Numerisch:** Kneser-Budget m21 (λ_min 0.98998, ‖Re R‖ 0.055, C_h 1.048); Kante radial m28/m10 (σ_min(T_{N+d}J_∞|P_N) erreicht 0.6905 ab d ≈ N/100).
**Offen:**
1. **(K2) von P** (p = 6, C₆ ≤ 300; angekündigt als p32, Datei im Format {"K","b_mitte","b_radius","C","p"}) und (K3) (f_K(L) = L, f_K'(L) = √L;
   folgt aus P-V, P bestätigt mit p32). Danach gilt der Kneser-Satz ohne Annahme.
2. **Kante (a″)_{αN}, α ∈ [0.03, 0.1]:** nicht bewiesen. Der Störungsweg trägt nur für α ≥ 1 (radial) bzw. α ≥ 2.74 (Kneser, = Xs X5); für kleines α ist
   der Überlauf exponentiell (1.48^N radial gemessen, max|F| = 1.624 bei Kneser) — eine Eckaussage vom Typ (K_d), braucht die lokale Parametrix an L
   (O/W/X). Beobachtung für X: Lemma 3.1 (σ_{d+1}(J_N⁻¹) ≤ 1/σ_min[J_N; S]) gilt für JEDE Wahl von d Zusatzfunktionalen S; die aus 𝒯 gebildeten
   (B Q_N J_∞) sind dafür zu viele (grob geschätzt ≈ 0.3N grosse Singulärwerte), die Überlaufzeilen bleiben die bessere Wahl.
3. **P-V gegenlesen:** erledigt — Abgleich mit A (1.1e-43), Funktionalgleichung mit Ps Kugeln (bis 2e-138), Cauchy-Schranke (max|a_n|·1.30^n = 1.139).
   Befund: die Rundung der 140-stelligen Dezimalmitten fehlt in den Radien, die Radien sind teils abgerundet (an P gemeldet). Nicht gegengelesen:
   Beweisschritte 3–6 (p27 Taylor-Modelle/Kantorovich, p28 Beltrami/Stoilow/Paulsen–Cowgill-Prüfungen A/B) — Ps Prüfliste stand noch aus.

**Wiederaufnahme (nächster Schritt):**
1. Sobald Ps (K2)-Datei in gemeinsam/ liegt: `cd m_modell_eckmode; python m26_kneser_gesamt.py ../gemeinsam/<datei>.json 7 23` (≈ 11 min) —
   erwartet c ≥ 0.9359, σ_min ≥ 0.6553; Ergebnis in §15.1 und auf dem BRETT eintragen (für X: KONSTANTEN.md).
2. P um Prüfung von §11 (M5, M7–M10, M13) und §13 bitten; Einbettungslemma ausschreiben.
3. P-V: Ps Prüfliste abarbeiten (p27/p28); Gegenprobe mit neuer Datei wiederholen (`python m27_pv_gegenprobe.py korr`).
4. Kante: nur in Abstimmung mit X/O/W (Thema (K_d)); kein eigener Strang ohne neuen Ansatz.
**Wichtige Dateien:** BERICHT.md §11/§13/§15/§17; Module mk_kneser_arb.py (Datenschnittstelle), mk_ecke_kneser.py, mk_kneser.py;
Skripte m16–m20 (radial), m22–m26 (Kneser; Gesamtlauf m26), m27 (P-V-Gegenprobe), m28 (Kante radial); Logs m26_kneser_pzert_T7_23.log,
m18_T6_22.log, m19_T22.log, m20_carleson.log, m27_pv_gegenprobe.log. Es laufen keine Prozesse.

## 18. Dateien (alle Runden)

m_kern.py (J_N für Polynom-f, Eck-SVD, Modelle `quadratisch`, `kubisch`) · m01_eckmode_modelle.py (σ für kubisch/quad/affin) ·
m02_eckvektoren.py (+npz) · m03_raten_scan.py (+m03_raten.jsonl, m03_scan_eps.log, m03_scan_theta.log) · m04_phasenraum.py ·
m05_adjungiert_profil.py (NUR letzte B Spalten — nahe L ungenau, ersetzt durch m07) · m06_pol_ecke.py (A_ζ) ·
m07_voll_profil.py (+m07_profil.log, npz, N = 400…1200) · m08_adjungiert_identitaet.py (Prüfung §4.2) ·
m09_singularitaet_ort.py (+m09_eps0.5.log, m09_eps0.31382.log; Ort-Test §3a) ·
m10_rechteck.py (+m10_rechteck.log; Rechteck-Abschnitte §8.2) · m11_ort_scan.py (+m11_ort_scan.log, m11_ort_scan_radial.log; §9) ·
m12_ueberlauf.py (Überlauf der Eckmoden, §8.2).
Runde 3: m13_randoperator.py / m13b_akkretiv.py / m13c_L2.py / m13d_bloecke.py (𝒯 auf P_N, Akkretivität numerisch, Frequenzblöcke) ·
m14_lokalisierung.py (IMS in der Scheibe, ersetzt durch m15) · m15_streifen.py (Streifen-Nyström, IMS-Budget numerisch) ·
**m16_mellin_zert.py (+m16_zert.log; Lemma M7) · m17_ende.py (+m17_ende.log; Lemma M8) · m18_kommutator.py (+m18_T5_13.log; Lemma M9) ·
m19_mitte.py (+m19_T13.log; Lemma M10) · m20_carleson.py (+m20_carleson.log; Lemma M13)** — Aufrufe: `m16 0.99 30 radial`, `m17 0.3 "4,5,6"`,
`m18 5 13 0.1 32`, `m19 13 0.5 16 2.8 96 8 5`, `m20 6 0.02` (zusammen < 2 min).
Runde 4 (Kneser): mk_kneser.py (Gleitkomma, Daten Strang A) · **mk_kneser_arb.py** (Datenschnittstelle (K1)–(K2), Einschlüsse F/F'/h/Fortsetzung) ·
**mk_ecke_kneser.py** (Eckschranken aus S2) · m21_kneser_streifen.py (+m21_kneser.log, numerisch) · **m22_kneser_ende.py** (+log) ·
**m23_kneser_kommutator.py** (+m23_kneser_T6_14.log) · **m24_kneser_mitte.py** (+m24_kneser_T14.log) · **m25_kneser_carleson.py** (+log) ·
**m26_kneser_gesamt.py** (Gesamtlauf, `python m26_kneser_gesamt.py [quelle] 6 14`, ≈ 10 min; +m26_kneser_vorzertifikat.log).
Runde 5: mk_kneser_arb Quelle "pzert[:C:p]" (Ps Datei + Schwanzparameter) · m26_kneser_pzert.log, m26_kneser_pzert_500_5.log ·
m27_pv_gegenprobe.py (+log; Aufruf mit "korr" für korrigierte Radien) · m28_kante_ueberlauf.py (Überlauf/Rechteck radial, arb).
