# Strang C — Literatur & Analogien

Stand: 10.10.2026, 23:10 (Beginn 19:12; Runde 4 ab 22:10). Kennzeichnung der Aussagen: **bewiesen** / **Skizze** / **numerisch** /
**Vermutung**. Prüfstatus der Quellen:

- **[O]** Originaltext selbst gelesen (für die hier zitierte Aussage),
- **[S]** sekundär geprüft (Abstract, Verlagsseite, Literaturverzeichnis einer selbst gelesenen Arbeit, Suchtreffer mit Daten),
- **[U]** ungeprüft, aus dem Gedächtnis — vor Weiterverwendung nachschlagen.

Eigene Prüfrechnungen: `c01_normtest.py` … `c09_effektiv_s.py` (Logs `*.log` daneben). Kneser-Motor nicht benutzt; Kneser-Daten nur aus
`gemeinsam/`.

---

## 0. Kurzfassung

1. **Es gibt keinen Konvergenzbeweis für irgendein Abschneide-Verfahren, das gegen Kneser laufen soll** — weder für unsere
   Polynom-Wurzel noch für das lineare Gegenstück (Walker/Robbins-Matrix-slog, abgeschnittene Julia-Gleichung, Helms'
   Carleman-Wurzel). Belegt durch Walker 1991 [O] („no proof of convergence is known“), Paulsen 2016 [O] („very difficult to
   prove“), Trappmann 2011 [O] (Beweis nur für f(x) = bx, exp „out of reach“), Forum-TPID 18 [O] (offen seit 2022) und
   Levensteins 2021 angekündigten, bis 04/2026 nicht erschienenen Beweis [O].
2. Numerisch wählen **alle** diese Verfahren Kneser: Walkers c₁^{(N)} (N = 30…80) nähert slog_K'(0) algebraisch (≈ N^{-2.4});
   Paulsens abgeschnittene Julia-Gleichung trifft 1/slog_K'(0) auf 18 Stellen (eigene Nachrechnung, Abschnitt 3).
   Warnbeispiel: für Basis √2 sah „intuitiv = regulär“ bei Grösse 100 richtig aus und ist falsch (Trappmann, Forum tid 153 [O]).
3. **Eigener Fund (numerisch), wichtig für die Beweisstrategie:** Misst man Urbild **und** Bild mit denselben Gewichten
   |L|^k (Radius = Konvergenzradius), so ist J_N gut konditioniert: |J_N⁻¹|_{l¹} = 1.8 … 7.8 (N = 10 … 194), mit
   Polynomgewicht (1+k)⁴ sogar praktisch konstant 25.4 → 28.2 (N = 60 … 194); die kleinsten l²-Singulärwerte häufen sich bei
   0.73 ≈ 1/|L|. Die „fast-singulären Richtungen“ des README (1.48e-12 …) entstehen nur durch die Gewichtspaarung (1.3 | 0.5);
   die Konkurrenz-Tangenten (Koeff. ~ k^{2.45}|L|^{-k}) liegen gar nicht im Raum, Kneser (~ k^{-6.45}) schon (für s < 5.45).
   Damit passt der klassische Rahmen **„diskrete Stabilität + Konsistenz ⇒ Konvergenz“** (Krasnosel'skii–Vainikko, Keller).
   Offen bleibt die **nichtlineare Kontrolle**: Vorwärtsoperator und F'' verlieren Radius (≈ 1.63^N), das affin-kovariante
   Worst-Case-ω (Deuflhard–Heindl) wächst ≈ 1.25^N (h = αω < 1/2 nur bis N ≈ 75); **entlang des tatsächlichen Zuwachses**
   p_N − q_N ist die Konstante aber klein und α·ω_eff fällt (in l¹_{|L|,4}: h_eff = 7.5e-9 → 2.5e-10 für N = 40 → 194) —
   ein Kantorovich/Krawczyk in Kegelform sollte für alle N tragen.
3a. Nebenbefund **[nur bis N ≈ 600–900 richtig, siehe Korrektur in §3.5]**: p_N konvergiert (numerisch) auch auf dem Kreis |x − a₀| < |L − a₀| um a₀ = f_K(0) (z. B. x = 1.8,
   1.4+0.9i), nicht aber bei x = 2 — passt zur Ersatzform, nicht zum Walsh/Hermite-Lemniskatenbild.
4. Transferoperator-Theorie (Bandtlow–Jenkinson u. a.) ist **bei 0 nicht anwendbar**: Sie braucht ein Gebiet U mit
   f⁻¹(U) ⋐ U; wegen f⁻¹∘f⁻¹ = log kann kein solches U den Punkt 0 enthalten (log 0).
5. Rangliste (Abschnitt 5): (1) Stabilität in l¹_{|L|,s} (s ≈ 4) + affin-kovarianter Kantorovich in Kegelform,
   (2) Abel-Koordinaten: Linearisierung = Differenzenoperator u ↦ u(w+½) + u(w), Kneser = Nörlund-Hauptlösung (liefert die
   a-priori-Schranke für (1)), (3) Minimallösungs-/Schwanzkorrektur-Technik (Miller–Olver–Pincherle, Fox-Beschleunigung,
   Flajolet–Odlyzko) für Konsistenz und nichtlineare Kontrolle, erst an linearen Modellen (Paulsens Julia-Gleichung).

### 0b. Kurzfassung Runde 2 (20:00–21:00, Details §9–§11)

- **(b) Konsistenz — Reduktion bewiesen (Satz 9.5):** α_N = |J_N(q_N)⁻¹F_N(q_N)| ≤ (1/2π)∫_Γ|μ|R_N|dζ| mit der skalaren
  Pol-Antwort R_N(ζ) = |J_N⁻¹B_N[(ζ−x)⁻¹]_{>N}| (geschlossene Form über [q^{N+1}/(ζ−q)]_{≤N}, Lemma 9.3). Unter (F) (Holomorphie
  von f_K um W̄ = D̄ ∪ D̄(a₀,|L−a₀|) mit radialen Schlitzen, Sprung ~ u^{5.447}) und (H_R) (R_N(L(1+u)) ≤ C N^a e^{−κN^a u}, exp.
  Zerfall auf einer Kontur ausserhalb W̄) folgt α_N ≤ C N^{−5.447a}; mit a = 0.883 (Aussenwinkel bei L) genau die gemessene Rate
  N^{−4.81}. (H_R) numerisch bis N = 194 (Skalengesetz kollabiert; Schranke nur Faktor ≈ 2 über α_N). **Keine** Schranke für
  |J_N⁻¹| nötig — wichtig, da (a) nach Strang D in keiner Norm gilt.
- **Negativ (numerisch):** |F_N(q_N)| wächst exponentiell (4.5e20 bei N = 160); ebenso die Antwort auf einzelne Schwanzkoeffizienten
  (1e25) ⇒ „(a) + |F| → 0“ und alle Koeffizienten-Majoranten sind aussichtslos; nur kohärente (Cauchy-)Darstellungen tragen.
- **Gesetz R (numerisch, 44 Gitterpunkte):** R_N(ζ)^{1/N} → max(ρ/|ζ|, R₀/|ζ−a₀|); W ist das Konvergenzgebiet von p_N.
  Erklärt Es Möbius-Schwelle exakt: Pol von f_M in D(a₀,R₀) ⇔ c > π/3, Rate sin(c/2)/cos c.
  **[Korrektur Runde 4, §15: ausserhalb W̄ ist das wahre Gesetz e^{−N g_E(ζ)} (Bernstein–Walsh scharf), genauer Poisson–Jensen mit
  den Randdaten auf ∂W (trifft auf 0.0004); die Max-Formel sagt zu schnellen Zerfall voraus. Im Inneren bleibt sie Beobachtung.]**
- **(c) Nichtlinearität:** bewiesen sind die exakte Zerlegung N_q(δ) = [Q_f(δ)] − [Q_t(δ)] (Lemma 10.1) und der Sekantensatz
  p_N − q_N = Ā⁻¹B_N t_N (Satz 10.4). Kantorovich mit Betragsabschätzungen hält um x₁ = q+δ1 bis N ≈ 230 (h = 1.2e-6 bei
  N = 160), wächst aber ≈ 1.19^N; für alle N ist ein kohärentes Argument nötig (§10.3). **Nicht bewiesen.**
- **Korrektur zu Runde 1 (Punkt 3 oben):** „|J_N⁻¹|_{l¹_{|L|,4}} praktisch konstant“ gilt nur bis N ≈ 194; Strang D (BRETT 20:41,
  20:53) zeigt eine exponentiell wachsende Eckspitze bei L (≈ 1.012^N) und die Instabilität am Zentrum q_N für N ≳ 190 —
  (a) im Sinn sup_N |J_N⁻¹| < ∞ ist vermutlich in jeder natürlichen Norm falsch. Die Runde-2-Ergebnisse brauchen (a) nicht.
- **Vermutung Q (§11)** als gemeinsame lineare Aussage (statt (a)): J_N⁻¹ reproduziert reguläre Funktionen auf W bis auf
  poly(N)·(Taylor-Reste bei 0 und bei a₀); impliziert Gesetz R und (b).

### 0c. Kurzfassung Runde 3 (21:05–22:10, Details §12)

- **Bewiesen:** Die Pol-Antwort ist rational in ζ: y_ζ = ζ^{−N−1}Σ_{m≤N} z_m(ζ−a₀)^{−m−1} (Lemma 12.1). Daraus mit dem
  Maximumprinzip in ζ (Bernstein–Walsh, Satz 12.2): R_N(ζ) ≤ S_N e^{−(N+2)g_E(ζ)}, g_E explizit (Möbius + Potenz, Öffnungswinkel
  158.93°, α = 0.882964, Kapazität 1.5169 — identisch mit Strang Bs unabhängigem Lemma 7.2/7.3). Damit folgen (H_R1), (H_R2) und
  (H*) aus der einen Randgrösse S_N = max_{∂W} R_N; Satz 12.3: (F) + S_N = o(N^{5.69}) ⇒ α_N ≤ C S_N N^{−5.69}. Satz 12.5:
  Vermutung Q für auf W̄ holomorphe Funktionen (exponentiell) unter derselben Bedingung. Satz 12.6: Poisson–Jensen-Variante mit
  explizitem harmonischem Mass.
- **Numerisch:** S_N = R_N(L) = 10.4 … 145.8 (N = 40 … 400), ≈ 0.33N^{1.12}, Maximum in der Ecke bis N ≈ 400; Schranke ≈ 10·α_N.
  A = J_N(q_N) und A_f = P_N J_∞(f)P_N geben identische Antworten, obwohl ihr Unterschied (Multiplikator [t'∘q]) bis 3e30 gross ist.
- **Neu und wichtig:** Die Pol-Antwort auf dem fernen D-Bogen wächst ≈ 1.0115^N (Strang Ds Eckspitze) und überholt die Ecke bei
  N ≈ 420 — (S)/(H*) in polynomialer Form ist für grosse N vermutlich falsch. Die Konsistenz δ1 sieht die Spitze bis N = 400 nicht
  (Exponent −4.87). Satz 12.6 rettet das Argument, solange die Spitze die Eckumgebung |η−L| ≲ 0.05 meidet (bis N = 400 erfüllt);
  für alle N braucht man, dass die Spitzenanregung durch Pole nahe L exponentiell klein ist — offene Frage ans innere Problem.

### 0d. Kurzfassung Runde 4 (22:10–23:10, Details §13–§15)

- **Eck-Bedingung (E), so schwach wie möglich — Satz 13.1 bewiesen:** Zerlege die Pol-Antwort mit Strang Ds Eckmoden,
  y_ζ = Q_N y_ζ + (I − Q_N)y_ζ (Q_N = Projektor auf die Eck-Ausgänge u_i). (E_R): Restanteil polynomial auf ∂W (Ds Rest-Schranke);
  (E_K): Eckanteil ≤ C N^κ e^{−ε(N+2)g(ζ)} **nur** auf dem Schlitz L(1+u), u ≤ u₀, und auf einer Niveaulinie {g = d}, ε > 0 beliebig
  klein, κ < 5.693. Dann α_N = O(N^{κ−5.693}) für alle N. Massgeblich sind nur Pole im Abstand ~N^{−0.883} ausserhalb der Ecke.
  Das Eckfunktional muss ⟨u_i, y_ζ⟩ sein (fällt wie ζ^{−N−2}), nicht ⟨u_i, Π_N k_ζ⟩.
- **Randform (Lemma 13.2, Satz 13.3):** (E_K) ⇐ ω_ζ(λ_N) ≤ (1−ε)g mit **vorzeichenbehafteter** Randrate λ_N; ohne Nullstellen auch
  notwendig. Ecke: Grenzwert (1/κ₀)∫λ(s)s^{−2}ds, κ₀ = 6.273; Stufenform: spitzenfreier Eckbogen |η−L| ≤ 0.011 trägt Rate 0.0119,
  ≤ 0.025 trägt 0.031. „Feste Rate × |s|^p“ reicht nicht — die Rate selbst muss zur Ecke hin verschwinden; das Vorzeichen auf dem
  A-Bogen zählt (Schlitz bei arg τ = 78°, neigt zum A-Bogen). Einteilige Form ohne Zerlegung: Satz 13.4 (Θ_N ≤ 1−ε).
- **Gegenprobe bis N = 1000 (numerisch, §14):** α_N aus As Π_N f_K: 9.95e-15 bei N = 1000, Exponent −4.85, Exponentialrate
  (2 ± 2)e-5 ⇒ δ1 spürt die Spitze nicht. Pol-Antworten: Ecke ~N^{1.10}, Schlitzpunkte in festem Verhältnis zur Ecke; die Spitze
  (0.0118/Grad) überholt die Ecke bei N ≈ 450, ihre Front erreicht |η−L| ≈ 0.025 bei N ≈ 900; ihre Rate fällt zur Ecke hin
  (0.0118 → 0.0016 bei |η−L| = 0.012, WKB-artig); A-Bogen fällt. Einteilige Bedingung direkt gemessen:
  Θ_N ≤ 0.105 für N ≤ 1000 (Schlitz negativ, positiv nur auf der Niveaulinie durch die ferne Spitze, dort ≤ 0.23 beschränkt).
- **Gesetz R korrigiert (§15):** ausserhalb W̄ ist es Poisson–Jensen bzw. e^{−Ng} (Streuung 0.0004 bzw. 0.009, Mittel −0.001 bzw. −0.004, statt
  0.021 / +0.032 für die Max-Formel); die Max-Formel bleibt nur im Inneren als Beobachtung. **Korrektur §3.5/3a:** p_N konvergiert im Halbmond nur bis
  N ≈ 600–900, dann exponentielle Divergenz (c26b, wie Strang B).

---

## 1. Problemstruktur, wie sie die Literaturwahl bestimmt

- System: F_N(a)_k = [x^k](p∘p) − 1/k!, k ≤ N; Jacobi-Matrix J = C(p)ᵀ-artig + Toeplitz: J[k][j] = [x^k]p^j + [x^{k−j}](p'∘p).
  Der erste Summand ist die (transponierte) **Carleman-Matrix** von p, voll besetzt, weil p(0) = a₀ ≠ 0.
  Faktorisierung: p = a₀ + g, g(0) = 0 ⇒ δ∘p = (τ_{a₀}δ)∘g (Taylor-Verschiebung um a₀, obere Dreiecksmatrix, mal untere
  Dreiecks-Carleman-Matrix von g); für Polynome ist der Abschnitt davon exakt — das ist die „Ersatzform“ des README.
- Linearisierung an f = f_K in **Abel-Koordinaten** (eigene Rechnung, **bewiesen** formal): mit A = slog, T = sexp,
  δ = T'(A+½)·u(A) gilt δ∘f + (f'∘f)·δ = T'(A+1)·[u(A+½) + u(A)]. Also ist die Linearisierung zu I + S_{1/2}
  (Verschiebung um ½) konjugiert; Kern = ½-antiperiodische u (README). Nahe L: T'(A+½) ~ c(x−L), e^{±2πiA} ~ (x−L)^{±p},
  p = 2πi/L ≈ 4.44695 + 1.05794i ⇒ Konkurrenz-Tangenten ~ (x−L)^{1−p}, Koeffizienten ~ k^{2.45}|L|^{-k}; Kneser selbst
  ~ (x−L)^{1+p} (Koeffizienten ~ k^{-6.45}|L|^{-k}, README). Gleiche Form bei Levenstein, Forum tid 1043 [O] (Abschnitt 2.2).
- Alle Konkurrenzlösungen haben **denselben** Konvergenzradius |L|; sie unterscheiden sich nur in der **Stärke** der
  Randsingularität bei L, L̄. Jede Technik, die „Selektion durch Abschneiden“ erklären will, muss also algebraische (nicht
  exponentielle) Trennung verkraften: Verhältnis der Koeffizienten ~ k^{-8.9}.

---

## 2. Kommentierte Literaturliste

### 2.1 Abschneide-Verfahren für Iterations- und Abel-Gleichungen (Thema 1)

**[A1] P. L. Walker, Infinitely differentiable generalized logarithmic and exponential functions, Math. Comp. 57 (1991),
no. 196, 723–733.** [O] (JSTOR-Scan, über eretrandre.org/rb/files/Walker1991_111.pdf)
- Kern: C^∞-Lösung von G(eˣ) = G(x)+1 über Hilfsfunktion h(eˣ) = e^{h(x)} − 1; §4 „matrix method“: G = Σ cₙxⁿ,
  Σcₙ = 1, (1/k!)Σ nᵏcₙ = c_k; N×N-Abschnitt. Wörtlich: „no proof of convergence is known“. Werte c₁^{(N)} (berechnet von
  F. W. J. Olver, J. M. Smith, D. W. Lozier): 0.9159437391 (N=30), …7811 (40), …52632 (50), …55362 (60), …57310 (70),
  …58315 (80). Walkers eigene Lösung: G'(0) ≈ 0.91536. Schlusssatz: möglicherweise „two or three distinct generalized
  logarithms“.
- Eigene Nachrechnung (**numerisch**): slog_K'(0) = 0.91594605649953 (gemeinsam/kneser_slog_koeff_300.json) ⇒ Fehler
  2.32e-6, 1.28e-6, 7.9e-7, 5.2e-7, 3.3e-7, 2.25e-7: algebraisch ≈ N^{-2.4}, Grenzwert = Kneser, nicht Walkers C^∞-Lösung.
- Übertragbar: lineares Gegenstück unseres Problems (Taylor bei 0 abgeschnitten, Selektion durch Abschneiden). Liefert
  **keinen** Beweisschritt, aber das richtige Modellproblem (→ Strang E) und Beleg für „offen“.

**[A2] P. L. Walker, On the solutions of an Abelian functional equation, J. Math. Anal. Appl. 155 (1991), 93–110,
doi:10.1016/0022-247X(91)90029-Y.** [S] (Literaturverzeichnisse von A1, A3). Nach A1/A3 die Quelle der Matrixmethode; ebenso
**C. W. Clenshaw, D. W. Lozier, F. W. J. Olver, P. R. Turner, Generalized exponential and logarithmic functions, Comput. Math.
Appl. 12B (1986), 1091–1101** [S] (A1 nennt „the author's matrix method [1],[8]“). Nicht gelesen.

**[A3] H. Trappmann, The Intuitive Logarithm, arXiv:1008.1409 (v2, 12.04.2011).** [O]
- Kern: Def. 1 „intuitive solution“ x_m = lim (A↾N⁻¹ b↾N)_m; Theorem 1: für f(x) = bx ist das Abschnittspolynom explizit,
  α̃⁽ⁿ⁾(z) = Σ_{k=1}^n C(n,k)(−1)^{k+1}(1−z^k)/(1−b^k), und α̃⁽ⁿ⁾(x) → log_b x für x, b ∈ (0,1).
- Technik: Zeilenoperationen (alternierende Summen) machen das System dreieckig → Rekursion → Erzeugendenfunktion →
  geschlossene Formel; Konvergenz über Binomialsummen. §3: exp „out of reach with these elementary techniques“; numerisch
  konvergieren die Koeffizienten nicht gleichmässig, aber punktweise.
- Übertragbar: **einziger** vollständiger Konvergenzbeweis eines solchen Verfahrens; Vorlage für Strang E (lineare bzw.
  Möbius-Modelle mit expliziter Lösung). Für exp nicht.

**[A4] W. Paulsen, Finding the natural solution to f(f(x)) = exp(x), Korean J. Math. 24 (2016), no. 1, 81–106,
doi:10.11568/kjm.2016.24.1.81.** [O] (kkms.org, Volltext)
- S. 85–86: Abel-System abgeschnitten → „individual coefficients do not settle down beyond 7 or 8 places“ (seine Werte stimmen
  mit Knesers slog auf ≈ 7 Stellen, eigene Prüfung — passt zu Walkers langsamer algebraischer Konvergenz, nicht zu Divergenz).
- S. 86–87: **abgeschnittene Julia-Gleichung** λ(eˣ) = λ(x)eˣ, Matrix m_{ij} = j^{i−1} − (i−1)!/(i−j−1)!: „seem to converge
  to at least 16 places“, c₀ ≈ 1.091767351258320992; „it would be very difficult to prove that as we let n → ∞, the
  coefficients … converge to a single solution“.
- Eigene Nachrechnung (**numerisch**): 1/slog_K'(0) = 1.09176735125832099180 ⇒ Paulsens c₀ = Kneser auf ≈ 18 Stellen.
- Rest der Arbeit: Störung eˣ + ε(x) mit reellem Fixpunkt; „log ratio metric“ auf Abel-Funktionen; Prop. 1; Abb. 3: ersetzt man
  eˣ durch Taylor-Polynome T_m (Böttcher bei ∞), **wachsen** die periodischen Abweichungen exponentiell in m (Warnung: naive
  Polynom-Ersetzung von exp konvergiert nicht); Endlösung über 2 sinh, ≠ Kneser.
- Übertragbar: Die Julia-Gleichung ist **linear und hat genau den Operatortyp unserer Linearisierung** (Komposition mit der
  Iterationsfunktion plus Multiplikation mit ihrer Ableitung). Bestes lineares Modellproblem (→ E); ein Beweis dort wäre
  Schritt (a) „Stabilität“ in vereinfachter Form.

**[A5] G. Helms (Gottfried), Forum tid 842 (08./09.01.2017).** [O] Halbe Iterierte aus der Wurzel der abgeschnittenen
Carleman-Matrix (Newton N ← (MN⁻¹+N)/2), Grösse 16: 0.498568472273 + 0.876337510066x + 0.247418943917x² + …; Helms selbst:
die Wurzel ist keine Carleman-Matrix mehr, „not yet discussed systematically“. Eigene Prüfung (**numerisch**): Abweichung von
Kneser a₀..a₂ = 5.2e-6, 1.4e-6, −1.3e-4; unser p₁₆: 9.3e-11, 7.3e-10, −1.3e-9. **Anderes Verfahren**, ~10⁵-mal ungenauer.

**[A6] G. Helms, TPID 18 im Thread „open problems survey“ (Forum tid 162, 04.07.2022).** [O] „Is "polynomial" tetration (by
truncated Carleman-matrices) … asymptotic to the Kneser-solution?“ — „never had an idea how to approach at all a proof“.
Tabelle für 4^x, Grössen 4…24, Konvergenz auf ≈ 4 Stellen.

**[A7] J. D. Fox (jaydfox), Forum tid 63 (28.09./10.10./01.12.2007) und tid 101 (13.12.2007).** [O]
- Mechanismus des Abschneidefehlers: Ersetzt man im unendlichen System die Zeilen unter N durch die Einheitsmatrix und die
  rechte Seite dort durch die wahren Koeffizienten, ist das exakt; das gewöhnliche Abschneiden setzt sie 0 und „introduces
  errors in the first n coefficients to compensate“. Beschleunigung: die bekannten log-Singularitäten bei L, L̄ in den Schwanz
  einsetzen. Raten: beschleunigt ≈ 6.8 Bit pro Verdopplung von N (≈ N^{-6.8}), natürlich 2.5–4 Bit (≈ N^{-2.5…-4}, passt zu Walkers N^{-2.4}). „we don't know enough
  to formally prove it“.
- Übertragbar: genau D's Antwort-Matrix J_N⁻¹B_N; die Beschleunigung ist ein **besseres Startpolynom** für Kantorovich
  (Abschnitt 5, Technik 3).

**[A8] H. Trappmann (bo198214), Forum tid 153 (27.04.2008, 21.06.2010).** [O] Basis √2: intuitives slog vs. reguläres slog —
bei Matrixgrösse 100 scheinbar gleich, bei 300 „dilated oscillation“ ⇒ verschieden. **Warnung**: Übereinstimmung auf 1e-6
war trügerisch; „Abschneiden wählt die reguläre Lösung am Fixpunkt“ ist falsch. (Unser Fall: 2e-19 bei N = 194 und
monotone Verbesserung — deutlich stärker, aber kein Beweis.)

**[A9] S. Levenstein (sheldonison), Forum tid 1292 (09.02.2021), tid 1017 (04/2021), tid 1349 (11.09.2021), tid 1017
(27.04.2026).** [O] Kündigt 2021 einen „rigorous proof that one of my pari-gp matrix based slog programs converges to Kneser's
…“ an (für ein vereinfachtes Programm mit eingebauter Schröder-/θ-Struktur von Kneser). Im Forum-Spiegel bis 27.04.2026 nicht
erschienen; letzter Beitrag: er wolle sich nach der Pensionierung evtl. wieder dem Thema widmen. **Status: existiert nicht.**

**[A10] A. Robbins, Solving for the analytic piecewise extension of tetration and the super-logarithm (2005).** [U] Laut A3
Wiederentdeckung der Walker-Matrixmethode; kein Beweis.

**[A11] W. Paulsen, S. Cowgill, Solving F(z+1) = b^{F(z)} in the complex plane, Adv. Comput. Math. 43 (2017), 1261–1282,
doi:10.1007/s10444-017-9524-1.** [S] (Abstract, arch.astate.edu) Natürliche Bedingungen erzwingen Kneser; neue Iteration
(180 Knoten, 10⁻⁵⁰, 1–2 Stellen je Iteration). Bedingungen laut Forum-Zusammenfassung (JmsNxn, tid 162, 23.08.2021 [O]):
holomorph auf ℂ∖(−∞,−2], F(0)=1, F(z̄)=F(z)‾, lim_{y→∞}F(x+iy) = L, F(z+1) = e^{F(z)}.
- Übertragbar: Die Grenzwertbedingung in vertikaler Richtung ist dieselbe Selektionsinformation, die l¹_{|L|} kodiert
  (Konkurrenzen wachsen vertikal exponentiell ⇔ starke Singularität bei L). Paulsen–Cowgills Iteration löst pro Schritt
  genau das linearisierte Kneser-Problem — Kandidat für die a-priori-Schranke in Technik 2.

**Fazit 2.1:** Der Stand ist eindeutig offen; selbst das lineare Gegenstück (Walker 1991) ist seit 35 Jahren unbewiesen.

### 2.2 Iterationswurzeln, Eindeutigkeit, lokale Struktur bei L (Thema 4)

- **[B1] H. Kneser, Reelle analytische Lösungen der Gleichung φ(φ(x)) = eˣ und verwandter Funktionalgleichungen, J. reine
  angew. Math. 187 (1950), 56–67, doi:10.1515/crll.1950.187.56.** [S] (Wikipedia-Literatur, Paulsen 2016, Walker 1991).
- **[B2] G. Szekeres, Regular iteration of real and complex functions, Acta Math. 100 (1958), 203–258.** [S] (A3).
  **G. Szekeres, Fractional iteration of exponentially growing functions, J. Austral. Math. Soc. 2 (1961), 301–320** [S]
  (quellen.yaml); **K. W. Morris, G. Szekeres, Tables of the logarithm of iteration of eˣ−1, ebd. 321–333** [S] (A1).
- **[B3] J. Écalle, Théorie des invariants holomorphes, Publ. Math. Orsay 67-74 09 (1974).** [S] (A3).
- **[B4] M. Kuczma, B. Choczewski, R. Ger, Iterative Functional Equations, CUP 1990 (Encycl. Math. Appl. 32).** [S].
- **[B5] K. Baron, W. Jarczyk, Recent results on functional equations in a single variable, perspectives and open problems,
  Aequationes Math. 61 (2001), 1–48** [S]; **W. Jarczyk, Recent results on iterative roots, ESAIM Proc. Surveys 46 (2014)**
  [S]. Behandeln stetige/monotone Wurzeln, Stabilität, mengenwertige Varianten — **keine** Abschneide-Resultate, nichts
  Spezifisches zu analytischen Wurzeln von exp gefunden.
- **[B6] H. Trappmann, D. Kouznetsov, Uniqueness of holomorphic Abel functions at a complex fixed point pair, Aequationes
  Math. 81 (2011), 65–76** [S]; **D. Kouznetsov, Solution of F(z+1) = exp(F(z)) in complex z-plane, Math. Comp. 78 (2009),
  1647–1670** [S]; **W. Paulsen, Tetration for complex bases, Adv. Comput. Math. 45 (2019), 243–267** [S].
- **[B7] L. J. Crone, A. C. Neuendorffer, Functional powers near a fixed point, J. Math. Anal. Appl. 132 (1988), 520–529,
  doi:10.1016/0022-247X(88)90080-7** [S] (Wikipedia). Wikipedia gibt die Behauptung wieder, es gebe keine halbe e-Funktion,
  die analytisch und reell-auf-reell ist — widerspricht Kneser, muss sich auf einen engeren Rahmen beziehen. **Unsicher**,
  nicht gelesen; nicht verwenden.
- **[B8] S. Levenstein, Forum tid 1043 (23.12.2015–01.01.2016).** [O] Vermutete Form von Knesers halber e-Funktion bei L:
  h_K(z) = h_f(z) + Σ_{n≥1}Σ_m c_{nm}(z−L)^{np+m}, h_f = formale (reguläre) halbe Iterierte S⁻¹(√L·S(z)), p = 2πi/L;
  „first four derivatives … defined and equal to zero“ (für die Differenz). Mit c_{n0} = 0 (Leitterm (z−L)^{p+1}, wie in
  Abschnitt 1 hergeleitet) passt das zur Koeffizienten-Asymptotik k^{-6.45} des README. → Strang A: prüfen, ob c_{10} = 0.
- **[B9] P. Walker, The exponential of iteration of eˣ−1, Proc. Amer. Math. Soc. 110 (1990), 611–620** [S] (A1, A4).

Übertragbar: Die Eindeutigkeit ist im README schon erledigt (Bôcher). Die Literatur liefert hier nur die **lokale Struktur bei
L** (B8, Szekeres/Écalle-Transreihen), die man für präzise Schwanz-Asymptotik (Konsistenzschritt) braucht.

### 2.3 Projektions-/Abschnittsverfahren und Newton–Kantorovich (Themen 2 und 5)

- **[C1] M. A. Krasnosel'skii, G. M. Vainikko, P. P. Zabreiko, Ya. B. Rutitskii, V. Ya. Stetsenko, Approximate Solution of
  Operator Equations, Wolters-Noordhoff 1972, doi:10.1007/978-94-010-2715-1** [S] (Kapitel Projektionsverfahren, glatte
  Operatoren). Kernsatz [U]: isolierte Lösung mit stetig invertierbarer Ableitung + stabile Projektion ⇒ Projektionsgleichungen
  haben nahe Lösungen, quasi-optimale Fehler.
- **[C2] G. M. Vainikko, Galerkin's perturbation method and the general theory of approximate methods for non-linear
  equations, USSR Comput. Math. Math. Phys. 7(4) (1967), 1–41, doi:10.1016/0041-5553(67)90140-1** [S]; **G. M. Vainikko,
  Regular convergence of operators and approximate solution of equations, Itogi Nauki 16 (1979) / J. Soviet Math. 15 (1981)**
  [S]. Begriff „diskrete Stabilität + Konsistenz ⇒ Konvergenz“ für nichtlineare Gleichungen.
- **[C3] H. B. Keller, Approximation methods for nonlinear problems with application to two-point boundary value problems,
  Math. Comp. 29 (1975), 464–474** [U]; **F. Brezzi, J. Rappaz, P.-A. Raviart, Finite dimensional approximation of nonlinear
  problems I, Numer. Math. 36 (1980), 1–25** [U].
- **[C4] P. Deuflhard, G. Heindl, Affine invariant convergence theorems for Newton's method and extensions to related
  methods, SIAM J. Numer. Anal. 16 (1979), no. 1, doi:10.1137/0716001** [S] (Titel/Heft/DOI auf der SIAM-Seite; Autorenliste
  dort nicht angezeigt); **P. Deuflhard, Newton Methods for Nonlinear Problems: Affine Invariance and Adaptive Algorithms,
  Springer SCM 35 (2004, Nachdruck 2011)** [S]. Kantorovichs Satz in affin-kovarianter Form: Bedingung an
  |F'(x₀)⁻¹(F'(x)−F'(y))| statt |F'(x₀)⁻¹|·|F'(x)−F'(y)| [U: Kantorovich–Akilov].
  Übertragbar: **direkt**; ersetzt das h = βγη von Strang D durch h = αω, numerisch um Faktor bis 10⁷ kleiner (Abschnitt 3).
- **[C5] A. Böttcher, B. Silbermann, Introduction to Large Truncated Toeplitz Matrices, Springer 1999** [U]; **M. Lindner,
  Infinite Matrices and their Finite Sections, Birkhäuser 2006 (Frontiers in Mathematics)** [S]; **K. Gröchenig,
  Z. Rzeszotnik, T. Strohmer, Quantitative estimates for the finite section method and Banach algebras of matrices,
  arXiv:math/0610588** [S]. Stabilitätskriterien für Abschnitte bandbeherrschter Operatoren (Grenzoperatoren).
  Übertragbar: nur das Prinzip „Abschnittsverfahren konvergent ⇔ stabil“; unsere J ist in l¹_{|L|} **nicht**
  bandbeherrscht (Spaltennormen ~ 1.63^j), die fertigen Kriterien greifen nicht.
- **[C6] L. V. Kantorovich, V. I. Krylov, Approximate Methods of Higher Analysis (engl. 1958), Kap. I: unendliche Systeme,
  Reduktionsmethode für „reguläre“ Systeme** [U]. Diagonaldominanz liegt bei uns nicht vor.
- **[C7] S. Olver, A. Townsend, A fast and well-conditioned spectral method, SIAM Review 55 (2013), 462–489** [U]: unendliche
  Bandsysteme, adaptive QR, Konvergenz der Abschnitte. Prinzip ja, Struktur (Band) nein.
- **[C8] A. Griewank, On solving nonlinear equations with simple singularities or nearly singular solutions, SIAM Rev. 27
  (1985), 537–563, doi:10.1137/1027141** [S]. Geränderte Systeme bei fast-singulären Lösungen. Nach c04 **vermutlich unnötig**
  (keine fast-singulären Richtungen in l¹_{|L|}).
- **[C10] P. P. Zabrejko, K-metric and K-normed linear spaces: survey, Collect. Math. 48 (1997), 825–859** [S] (Fixpunktsätze
  in K-normierten Räumen = Normen mit Werten in einem Kegel, auch Skalen von Banachräumen) und **S. M. Rump, Verification
  methods: Rigorous results using floating-point arithmetic, Acta Numer. 19 (2010), 287–449, doi:10.1017/S096249291000005X**
  [S] (Krawczyk/Intervall-Newton, komponentenweise Einschliessung). Übertragbar: **Rahmen für die Kegelform** von Platz 1
  (komponentenweise Majoranten statt einer Norm; Kantorovichs Majorantenmethode [U]).
- **[C9] S. Day, J.-P. Lessard, K. Mischaikow, Validated continuation for equilibria of PDEs, SIAM J. Numer. Anal. 45 (2007)**
  [S, nur zitiert]; Radii-Polynome (Galerkin + Schwanzabschätzung + Kontraktion); **A. Haro u. a., The Parameterization Method
  for Invariant Manifolds, Springer 2016** [U]. Übertragbar: für endliche N (dort schon Krawczyk bis 120); für alle N fehlt
  eine Schwanzstruktur (Komposition ist im Schwanz nicht diagonaldominant).

### 2.4 Transfer- und Kompositionsoperatoren (Thema 2)

- **[D1] O. F. Bandtlow, O. Jenkinson, Explicit eigenvalue estimates for transfer operators acting on spaces of holomorphic
  functions, arXiv:0802.1638** [S] (Zeitschriftenangabe Adv. Math. 218 (2008) [U]): |λₙ| ≤ A·exp(−a n^{1/d}).
- **[D2] O. F. Bandtlow, J. Slipantschuk, Lagrange approximation of transfer operators associated with holomorphic data,
  arXiv:2004.03534 (2020)** [S]: Lagrange-Interpolation, exponentielle Konvergenz der Spektraldaten, Rate durch komplexe
  Kontraktionsverhältnisse der Zweige.
- **[D3] C. C. Cowen, B. D. MacCluer, Composition Operators on Spaces of Analytic Functions, CRC 1995; M. Kamowitz (Spektren
  gewichteter Kompositionsoperatoren)** [U].
- Bewertung (**bewiesen**, elementar): Der Rahmen braucht einen Symbol-Zweig φ mit φ(U) ⋐ U. In unserer Linearisierung ist das
  Symbol f (expandierend) bzw. f⁻¹ (halber Logarithmus). Aus f⁻¹∘f⁻¹ = log folgt: ein f⁻¹-vorwärtsinvariantes U mit 0 ∈ U
  müsste log 0 enthalten — unmöglich. (Konkret: f⁻¹(0) = T(−3/2) = log f_K(0) ≈ −0.696, und dort hat f⁻¹ seinen
  Verzweigungspunkt.) Die Theorie greift nur nahe L (attraktiver Fixpunkt von f⁻¹, Multiplikator L^{-1/2}, |·| = 0.853) —
  dort lebt Knesers eigene Konstruktion, nicht unser Taylor-Abschnitt bei 0. In Abel-Koordinaten ist die Linearisierung
  I + S_{1/2}: ein Differenzen-, kein Kontraktionsoperator. **Nicht übertragbar** auf den Kern des Problems.
- **[H1] A. Amini, Q. Zheng, Q. Sun, N. Motee, Carleman linearization … finite-section approximations, arXiv:2207.07755**
  [S] und **K. Glashoff, C. P. Ortlieb, arXiv:1705.00325 (Finite-Section-Methode für Kompositionsoperatoren zwischen
  Flächen)** [S]: andere Struktur (Fixpunkt/Zeitentwicklung bzw. Hilbertraum-Formen), nicht übertragbar.

### 2.5 Schlecht gestellte Probleme und „Abschneiden wählt die schwächste Singularität“ (Thema 3)

- **[E1] F. Natterer, Regularisierung schlecht gestellter Probleme durch Projektionsverfahren, Numer. Math. 28 (1977),
  329–341, doi:10.1007/BF01389972** [S]. **[E2] S. Kindermann, Projection methods for ill-posed problems revisited, Comput.
  Methods Appl. Math. (2016; online 2015, doi:10.1515/cmam-2015-0036, arXiv:1507.03364)** [S]: lokale vs. globale Konvergenz
  diskreter Lösungen bei exakten Daten; Beispiel einer beschränkten, nicht konvergenten Folge. **[E3] P. C. Hansen, The
  discrete Picard condition for discrete ill-posed problems, BIT 30 (1990), 658–672** [S]. **[E4] H. W. Engl, M. Hanke,
  A. Neubauer, Regularization of Inverse Problems, Kluwer 1996** [U].
  Bewertung: Passte, solange man „fast-singuläre Richtungen“ annahm (Rest fast senkrecht = diskrete Picard-Bedingung). Nach
  c04 ist das Problem in l¹_{|L|} **gut gestellt** (beschränkte Inverse) — der Regularisierungsrahmen wird für die lineare
  Seite nicht gebraucht. Kindermanns Gegenbeispiel zeigt, warum die Montel-+-Eindeutigkeits-Reduktion des README wertvoll ist:
  dort genügt Beschränktheit.
- **[E5] W. Gautschi, Computational aspects of three-term recurrence relations, SIAM Rev. 9 (1967), 24–82** [O für S. 24–28,
  Scan von cs.purdue.edu]: Definition „minimal at infinity“ (f_n/g_n → 0; Fussnote: Pincherles „soluzione distinta“, bei
  linearen Dgl. „principal solution“); S. 27 Gl. (0.6): die Minimallösung als Lösung des **abgeschnittenen** tridiagonalen
  Systems mit y_N = f_N ≈ 0; S. 27–28: für allgemeine Funktionalgleichungen mit Zerlegung S = S₁ ⊕ S₂ (S₁ dominiert) sei das
  Berechnen der Minimallösungen „not been thoroughly studied“ (Hinweis auf Clenshaw, Schäfke). Die genaue Fehlerformel des
  Miller-Algorithmus (§3) habe ich nicht gelesen [U];
  **F. W. J. Olver, Numerical solution of second-order linear difference equations, J. Res. Nat. Bur. Standards 71B (1967),
  111–129** [S]; **J. Wimp, Computation with Recurrence Relations, Pitman 1984** [S]; Pincherle 1894 [S, über Suchtreffer].
  Kern: Abschneiden (Miller rückwärts, Olver als Randwertproblem mit Abschnitt N) konvergiert genau gegen die **Minimallösung**,
  wenn sie existiert; Trennung darf algebraisch sein. Übertragbar: als **Mechanismus-Modell** der Selektion (Kneser ~ k^{-6.45},
  Konkurrenz ~ k^{2.45}, Verhältnis k^{-8.9}); als Beweiswerkzeug nur für lineare Modelle mit endlicher Rekursionsordnung
  (Strang E), nicht für unser volles System. Bemerkenswert: Olver selbst hat Walkers c₁^{(N)} gerechnet (A1).
- **[E6] J. P. Boyd, Chebyshev and Fourier Spectral Methods (2. Aufl. Dover 2001), „behavioral boundary conditions“** [U für
  die genaue Stelle]: Spektralabschnitte wählen an singulären Randpunkten automatisch die reguläre Lösung. Gleiche Analogie.
- **[E7] C. Lanczos (τ-Methode, 1938), E. L. Ortiz, The tau method, SIAM J. Numer. Anal. 6 (1969)** [U]: Abschneiden = exakte
  Lösung eines gestörten Problems (τ-Terme im höchsten Grad). Übertragbar als Sichtweise: p_N ist **exakte** halbe Iterierte von
  E_N := p_N∘p_N = exp + O(x^{N+1}); Fragestellung wird „Stabilität der Kneser-Wurzel unter Störung von exp bei hohen Graden“.
- **[E8] J. L. Walsh (1935), Bernstein–Walsh-Satz, maximale Konvergenz, Überkonvergenz** [U]; **H.-P. Blatt, E. B. Saff,
  M. Simkani, Jentzsch–Szegő type theorems for the zeros of best approximants, J. London Math. Soc. (2) 38 (1988), 307–316,
  doi:10.1112/jlms/s2-38.2.307** [S]. → Strang B (b05 Nullstellen, b02 Überkonvergenz): Für maximal konvergente
  Polynomfolgen häufen sich Nullstellen auf dem Rand des Konvergenzkreises (|x| = |L|); Überkonvergenz von **Teilsummen** gibt
  es nur bei Lücken (Ostrowski [U]). p_N sind keine Teilsummen: ihre Konvergenz ist algebraisch (B: N^{-7.9}), also weit von
  maximaler Konvergenz (0.727^N auf |x| ≤ 1) entfernt — die Abweichung bestimmt der Abschneideeffekt (Ein-Moden-Störung, B),
  nicht der Kneser-Schwanz; ihr Konvergenzgebiet D ∪ {|x−a₀| < |L−a₀|} (3.5) folgt aus der Ersatzform.
- **[E9] P. Flajolet, A. Odlyzko, Singularity analysis of generating functions, SIAM J. Discrete Math. 3 (1990), 216–240**
  [S für Inhalt, Seiten nicht bestätigt]; Flajolet–Sedgewick, Analytic Combinatorics, Kap. VI [S]. Übertragungssätze
  (1−z/L)^α ↦ L^{-n}n^{-α-1}/Γ(−α) mit Fehlerterm auf Δ-Gebieten. Übertragbar: **Konsistenzschritt** — präzise Asymptotik
  des Kneser-Schwanzes und damit von F_N(q_N), α_N, J_N⁻¹B_N q_hi, sobald die Singularitätsentwicklung bei L (B8) steht.
- **[E10] M. D. Kruskal, H. Segur, Asymptotics beyond all orders in a model of crystal growth, Stud. Appl. Math. 85 (1991)**
  [U]: Selektion einer Lösung aus einer Familie durch exponentiell kleine Effekte. Nur Analogie (unser Selektionsrest
  ~ (a₀/|L|)^N), kein Werkzeug.

### 2.6 Newton in Skalen (Thema 5)

- **[F1] E. Zehnder, Generalized implicit function theorems with applications to some small divisor problems I, Comm. Pure
  Appl. Math. 28 (1975), 91–140** [S]; **M. Berti, P. Bolle, M. Procesi, An abstract Nash–Moser theorem with parameters and
  applications to PDEs, Ann. IHP (C) 27 (2010), 377–399, doi:10.1016/j.anihpc.2009.11.010** [S] (Galerkin-Nash–Moser: Lösungen
  abgeschnittener Gleichungen, zahme Abschätzungen der abgeschnittenen inversen Linearisierung); Moser 1966, Nirenberg 1972,
  Nishida 1977 (abstrakter Cauchy–Kowalewski) [U].
- Bewertung nach c01–c03: Bei uns ist die **Inverse nicht** das Problem (beschränkt bei Radius |L|), sondern der Vorwärtsteil
  und F'' verlieren Radius um einen **festen** Faktor (Komposition: r ↦ |p|_r ≈ 1.63 r bei r = |L|). Nash–Moser-Schemata
  sind für kleine Verluste je Schritt gebaut (Kosten δ^{-τ}); ein fester Verlust kostet auf Grad-N-Polynomen c^N. **Nicht
  direkt übertragbar**; höchstens die Buchhaltungs-Idee (Normen-Skala r_k) für Technik 1(c).

### 2.6a Nächstes Analogon mit Komposition: Feigenbaum-Gleichung

- **[F2] O. E. Lanford III, A computer-assisted proof of the Feigenbaum conjectures, Bull. Amer. Math. Soc. (N.S.) 6 (1982),
  427–434, doi:10.1090/S0273-0979-1982-15008-X** [S] (AMS-Seite, Einleitung); Eckmann–Koch–Wittwer (Mem. AMS 1984) und
  Koch–Schenkel–Wittwer (SIAM Rev. 1996) [U]. Die Funktionalgleichung g(x) = −α g(g(−x/α)) enthält wie unsere g∘g; gelöst
  mit abgeschnittenen Potenzreihen + Newton + Schwanzschranken in Räumen analytischer Funktionen auf Kreisen.
- Unterschied (entscheidend): Dort sorgt die Reskalierung x ↦ x/α dafür, dass die innere Funktion den Kreis **kompakt in sich**
  abbildet — der Renormierungsoperator glättet (kompakt), die Abschnittsfehler sind exponentiell klein. Bei uns gibt es keine
  solche Reskalierung (Abschnitt 2.4: kein f⁻¹-invariantes Gebiet um 0), die Komposition **verliert** Radius. Übertragbar ist
  nur die Technik der Fehlerbuchhaltung (gewichtete l¹-Normen, Schwanzoperator), nicht die Kontraktion.

### 2.7 Differenzengleichungen in Abel-Koordinaten

- **[G1] N. E. Nörlund, Vorlesungen über Differenzenrechnung, Springer 1924** [S für Existenz/Begriff „Hauptlösung“; Inhalt
  der „Wechselsumme“ u(x+ω)+u(x) = φ [U]]. Hauptlösung = durch Summierbarkeit bzw. minimales Wachstum ausgezeichnete Lösung.
  Übertragbar: unsere Linearisierung ist u(w+½)+u(w) = g (Abschnitt 1); Kneser entspricht der Lösung ohne exponentielles
  vertikales Wachstum, Konkurrenzen addieren ½-antiperiodische ganze Funktionen. Liefert die **a-priori-Schranke** der
  unendlichen Linearisierung auf beschränkten Funktionen (Technik 2). Carlson-Satz/Bôcher (README) ist die Eindeutigkeitsseite
  desselben Arguments.
- **[G2] F. Loray, F. Touzet, S. M. Voronin, Two dimensional neighborhoods of elliptic curves: analytic classification in the
  torsion case, arXiv:2001.03033 (v2 24.09.2020), Theorem 9.5 mit Beweis, S. 37–38.** [O] Inhalt (sinngemäss):
  Δ holomorph auf dem **senkrechten** Streifen Σ_{A,B} = {A ≤ Re ζ ≤ B}, A < A + λ₁ < B ≤ 1 (λ₁ vermutlich Re λ), mit
  |Δ|_m := sup|Δ(ζ)||ζ|^m < ∞ für ein m ≥ 3 ⇒ es gibt eine **beschränkte** holomorphe Lösung von φ(ζ+λ) − φ(ζ) = Δ(ζ) auf
  Σ_{A,B+λ₁}, eindeutig bis auf Konstanten, und sup|φ| ≤ C(B−A)·|Δ|_m·∫_{Re t = B}|dt|/|t|^m. Beweis: Cauchy-Zerlegung
  Δ = F₀⁺ − F₀⁻ über die beiden Randgeraden, F_n^± = F₀^±(ζ ∓ nλ), φ = Σ_{n≥1}F_n⁺ + Σ_{n≥0}F_n⁻ (Teleskop).
  Übertragbar (**Skizze**, Lücken benannt): Das ist genau die Schablone für Technik 2 — Verschiebung quer über einen
  senkrechten Streifen, der breiter ist als die Verschiebung (bei uns ½ gegen Breite von Ω ≥ 2.067). Anpassungen: (i) Plus-
  statt Minuszeichen → alternierende Reihe Σ(−1)^n…, Eindeutigkeit dann exakt (beschränkte ½-antiperiodische ganze Funktion
  = 0); (ii) Ω ist kein gerader Streifen (linker Rand ≥ −2 − 0.24|y|) → Cauchy-Integral über ∂Ω; (iii) die Abfallbedingung
  |ζ|^{-m} in w entspricht in x nur einem **schwachen** Verschwinden bei L, L̄ (vertikal ∞ in w ↔ x → L, |x−L| ~ e^{-c|Im w|}),
  ist also für Reste, die bei L wie |x−L|^ε verschwinden, mehr als erfüllt.

---

## 3. Eigene Prüfrechnungen (alle **numerisch**)

### 3.1 Normen der Inversen (c01, J an p_N, l¹ mit w_k = r^k(1+k)^s in Urbild und Bild)

| N | r = 0.5 | r = 0.8 | r = 1.0 | r = 1.3 | r = |L| | r = |L|, s = 4 |
|---|---|---|---|---|---|---|
| 10 | 3.3e1 | 3.5 | 2.2 | 1.79 | 1.79 | 16.3 |
| 20 | 5.6e4 | 6.5e1 | 6.5 | 2.18 | 2.09 | 20.8 |
| 30 | 3.0e8 | 3.1e3 | 3.5e1 | 2.58 | 2.27 | 22.7 |
| 40 | 2.0e12 | 2.2e5 | 3.3e2 | 3.16 | 2.48 | 24.0 |
| 50 | 1.5e16 | 2.2e7 | 3.4e3 | 3.95 | 2.58 | 24.8 |
| 60 | 1.8e20 | 2.3e9 | 4.8e4 | 4.87 | 2.77 | 25.4 |

|J_N| selbst wächst bei r = |L| wie 6.8e1, 7.6e3, 9.0e5, 1.1e8, 1.3e10, 1.6e12 (≈ 1.63^N = (Σ|a_k||L|^k / |L|)^N).

Grössere N (c07, r = |L|, s = 0): |J_N⁻¹| = 2.77 (60), 3.24 (80), 3.79 (100), 4.30 (120), 5.63 (150), 7.84 (194). Das wächst
langsam, aber mit zunehmender log-log-Steigung (0.6 → 1.3) — verträglich mit ≈ 1.008^N. Gegen α ~ N^{-6} unschädlich bis
N ~ 10⁴, für einen Satz „für alle N“ aber in genau dieser Norm **nicht** gleichmässig. Gewichts-Scan dazu (c08):

| N | r=|L|,s=0 | r=|L|,s=2 | **r=|L|,s=4** | r=1.40,s=0 | r=1.40,s=2 | r=1.45,s=0 | r=1.50,s=4 |
|---|---|---|---|---|---|---|---|
| 60 | 2.77 | 2.45 | **25.38** | 2.55 | 2.43 | 2.45 | 1.6e2 |
| 100 | 3.79 | 2.95 | **26.80** | 3.12 | 2.64 | 2.60 | 1.2e3 |
| 150 | 5.63 | 3.89 | **27.69** | 3.80 | 3.05 | 2.90 | 2.3e4 |
| 194 | 7.84 | 4.96 | **28.20** | 4.34 | 3.39 | 3.30 | 3.5e5 |

Mit Polynomgewicht (1+k)⁴ bei Radius |L| ist |J_N⁻¹| praktisch konstant (Zuwächse 1.42, 0.89, 0.51) — **gleichmässige
Stabilität (numerisch)**. Kneser liegt in l¹_{|L|,s} genau für s < 5.45 (Koeffizienten ~ k^{-6.45}|L|^{-k}), die
Konkurrenz-Tangenten (~ k^{2.45}) für kein s ≥ 0. Arbeitsfenster: **4 ≤ s < 5.45**. Radien r > |L| verlieren schnell.

Lineare Modelle (c08, r = |L|, s = 0, Paulsens Zeilenwahl i = 0…N−1): Walker-Abel-System 3.54, 4.09, 4.27, 4.49, 4.74
(N = 20, 40, 60, 80, 120) — ebenso stabil; Paulsens Julia-System 4.34, 4.95, 5.64, 6.74, 10.67 (bei r = 1.45: 4.6 → 6.2).

### 3.2 Singulärwerte (c04, N = 30, l² mit D = diag r^k)

| (r_X, r_Y) | kleinste Singulärwerte |
|---|---|
| (1.3, 0.5) | 1.48e-12, 5.40e-12, 1.81e-11, 6.07e-11, 1.99e-10 (= README) |
| (1.3, 1.3) | 0.744, 0.745, 0.746, 0.748, 0.752 |
| (|L|, |L|) | 0.732, 0.732, 0.733, 0.735, 0.737 |
| (1, 1) | 0.076, 0.146, 0.802, 0.804 |
| (0.5, 0.5) | 6.4e-9, 2.4e-8, 1.2e-2, 2.1e-2 |

Fast-singuläre Richtungen erscheinen paarweise erst für r < |L| (L, L̄) — genau dann, wenn die Konkurrenz-Tangenten in der
Norm „klein“ werden. Bei r = |L| verschwinden sie.

### 3.3 Kantorovich an q_N (c02/c03): Produktform vs. affin-kovariant (ω = Worst-Case an x₀ über e_i, e_j)

| N | r | d = |p_N−q_N| | α | ω | h_aff = αω | h_prod = αβγ |
|---|---|---|---|---|---|---|
| 10 | |L| | 2.1e-5 | 2.1e-5 | 1.2e1 | 2.5e-4 | 1.5e-2 |
| 30 | |L| | 2.0e-7 | 2.0e-7 | 2.9e3 | 5.7e-4 | 6.8 |
| 50 | |L| | 1.9e-8 | 1.9e-8 | 4.0e5 | 7.7e-3 | 1.8e4 |
| 60 | |L| | 8.5e-9 | 8.5e-9 | 4.9e6 | 4.1e-2 | 1.2e6 |
| 70 | |L| | 3.9e-9 | 3.9e-9 | 5.2e7 | 0.20 | 8.9e7 |

- α = d auf allen Stellen (erster Newton-Schritt trifft p_N) — wie bei D.
- ω wird bei (i, j) ≈ (N/5, N) angenommen (höchster Grad); r = 1.45, 1.55, 1.7 und Polynomgewichte (1+k)^s verschlechtern h.
- Extrapolation: h_aff = 1/2 bei N ≈ 75 (Krawczyk-Zertifikat reicht schon bis 120). **Das Hindernis für alle N ist allein die
  nichtlineare Konstante.**

### 3.4 Konstante entlang des tatsächlichen Zuwachses (c05)

δ = p_N − q_N, l¹_{|L|}; ω_eff = |J(q_N)⁻¹F''(q_N)[δ,δ]|/|δ|², κ = |q_N − J⁻¹F(q_N) − p_N|/|δ|² (beobachtete Newton-Konstante).

| N | |δ| = α | ω_eff | h_eff = αω_eff | κ |
|---|---|---|---|---|
| 20 | 1.27e-6 | 0.12 | 1.6e-7 | 0.06 |
| 40 | 5.06e-8 | 0.19 | 9.5e-9 | 0.09 |
| 60 | 8.52e-9 | 0.56 | 4.8e-9 | 0.28 |
| 80 | 1.94e-9 | 1.17 | 2.3e-9 | 0.58 |
| 100 | 6.80e-10 | 2.74 | 1.9e-9 | 1.37 |
| 120 | 2.81e-10 | 3.57 | 1.0e-9 | 1.79 |

Mit Polynomgewichten (c09, w_k = |L|^k(1+k)^s, N bis 194):

| N | s=0: α / ω_eff / h_eff | s=3: α / ω_eff / h_eff | s=4: α / ω_eff / h_eff |
|---|---|---|---|
| 40 | 5.1e-8 / 0.19 / 9.5e-9 | 2.0e-3 / 4.0e-6 / 7.9e-9 | 7.1e-2 / 1.1e-7 / 7.5e-9 |
| 80 | 1.9e-9 / 1.17 / 2.3e-9 | 4.9e-4 / 4.0e-6 / 2.0e-9 | 3.3e-2 / 5.7e-8 / 1.9e-9 |
| 120 | 2.8e-10 / 3.57 / 1.0e-9 | 2.1e-4 / 3.9e-6 / 8.3e-10 | 2.1e-2 / 3.7e-8 / 7.8e-10 |
| 160 | 7.0e-11 / 7.52 / 5.3e-10 | 1.2e-4 / 3.7e-6 / 4.4e-10 | 1.5e-2 / 2.8e-8 / 4.2e-10 |
| 194 | 2.8e-11 / 11.1 / 3.1e-10 | 7.9e-5 / 3.3e-6 / 2.6e-10 | 1.2e-2 / 2.1e-8 / 2.5e-10 |

Bei s = 4 (gleichmässig stabile Norm, 3.1) gilt also numerisch alles zugleich: |J_N⁻¹| ≤ 28.2, α_N → 0 (≈ N^{-1.1}),
ω_eff ≤ 1.1e-7 und fallend, h_eff ≤ 7.5e-9. Entlang der tatsächlichen Richtung wächst die nichtlineare Konstante (bei s = 0)
nur langsam (≈ N^{1.9}), h_eff **fällt** in jeder Gewichtung. Der
exponentielle Worst-Case (3.3) kommt aus Richtungen, die der Zuwachs nicht besucht. Konsequenz: Kantorovich/Krawczyk in
**Kegel-/Majorantenform** (Kugel ersetzt durch {|δ_k| ≤ ε m_k |L|^{-k}} mit passendem Profil m_k) statt in der vollen
l¹-Kugel — die Schranke für J⁻¹F'' wird dann nur auf dem Kegel gebraucht. (Achtung: δ hat **grosse** relative Anteile in den
höchsten Graden — p_N's Spitzenkoeffizienten sind ~1.59^{-N} (B), also δ_k ≈ −(Kneser-Koeffizient) nahe k = N; das Profil
muss das abbilden.)

### 3.5 Konvergenzgebiet ausserhalb von D (c06)

Cauchy-Test |p_N(x) − p_194(x)| für N = 40…180 (a₀ = p(0) = 0.49856, |L − a₀| = 1.34935, |L||L − a₀| = 1.85476):

| x | |x − a₀| | |x(x−a₀)| | in D? | Verlauf | Urteil |
|---|---|---|---|---|---|
| 1.6 | 1.10 | 1.76 | nein | 5.9e-12 … 6.0e-19 | konvergent |
| 1.8 | 1.30 | 2.34 | nein | 3.6e-9 … 9.4e-14 | konvergent |
| 1.4+0.9i | 1.27 | 2.12 | nein | 2.2e-9 … 2.9e-15 | konvergent |
| 0.5±1.3i | 1.30 | 1.81 | nein | 2.0e-8 … 1.8e-13 | konvergent |
| 2.0 | 1.50 | 3.00 | nein | 1.0e-2 … 1.1e-2 | nicht konvergent |

Befund: konvergent auf (mindestens) D ∪ {|x − a₀| < |L − a₀|} — der Vereinigung der Konvergenzkreise von Knesers Reihe um 0
und um a₀ = f_K(0); die Hermite/Walsh-Lemniskate |x(x−a₀)| < |L||L−a₀| ist widerlegt (1.8, 1.4+0.9i). Erklärung über die
Ersatzform: der N-Jet von p_N bei a₀ ist an exp∘g⁻¹ gebunden. Raten überall algebraisch (≈ N^{-7.9}, wie B).
Literatur dazu: keine direkte; Walsh (1935) wäre das Vergleichsbild für Interpolation [U].
**Korrektur (Runde 4, c26b; Strang B, BRETT 22:22):** Das gilt nur bis N ≈ 600–900. Mit Strang As Π_N f_K (= p_N bis 1e-24) divergiert die Folge im Halbmond danach exponentiell: x = 1.84 ab N ≈ 600, x = 1.8 ab N ≈ 650 (Π_1000 f_K(1.8) = 2.7e29), x = 1.4+0.9i ab N ≈ 900; x = 1.6 bis N = 1000 noch konvergent. Ursache ist vermutlich die Eckspitze (§12.9, §14.2), deren Anteil ausserhalb |x| = ρ mit (|x|/ρ)^N verstärkt wird. Das Konvergenzgebiet **für alle N** ist also nur D̄ (bzw. ein Teil davon); W bleibt das richtige Gebiet für die **Polvariable** ζ der Pol-Antwort (Lemma 12.1), nicht für die Auswertungsvariable x.

### 3.6 Literaturwerte gegen Kneser

- Walker c₁^{(N)} → slog_K'(0): Fehler ≈ N^{-2.4} (A1). Paulsens Julia-c₀ = 1/slog_K'(0) auf 18 Stellen (A4). Paulsens
  Abel-Koeffizienten = Kneser auf ≈ 7 Stellen (A4). Helms' Carleman-Wurzel (16) ≈ 1e-4 genau, unser p₁₆ ≈ 1e-9 (A5).

---

## 4. Welcher Beweisschritt wird wovon getragen?

| Schritt | Inhalt | Technik / Quelle | Stand |
|---|---|---|---|
| (a) Stabilität | sup_N |J_N⁻¹|_{l¹_{|L|,s}} < ∞ (s ≈ 4) | Vainikko/Keller-Rahmen [C1–C3]; a-priori-Schranke über Abel-Koordinaten + Nörlund-Hauptlösung [G1], Paulsen–Cowgill-Iteration [A11]; Finite-Section-Prinzip [C5] | numerisch ≤ 28.2 (s = 4, N ≤ 194); kein fertiger Satz anwendbar |
| (b) Konsistenz | α_N = |J⁻¹F(q_N)| → 0 (s = 0: ≈ N^{-6}; s = 4: ≈ N^{-1.1}) | Singularitätsanalyse [E9] mit lokaler Entwicklung bei L [B8]; Fox' Schwanzbild [A7] = D's J_N⁻¹B_N | numerisch klar |
| (c) Nichtlinearität | Kegel-Version von ω_N·α_N < 1/2 für alle N | affin-kovarianter Kantorovich [C4]; besseres Startpolynom (Schwanzkorrektur [A7], τ-Sicht [E7]); Carleman-Multiplikativität | **Engpass**: Worst-Case ω ~ 1.25^N, aber ω_eff·α ≤ 7.5e-9 (c09) |
| (d) Grenzübergang | p_N → f_K | Montel + Eindeutigkeit (README, bewiesen) oder direkt aus (a)–(c) | erledigt, sobald (a)–(c) |
| Selektion (Erklärung) | warum Kneser | Minimallösung [E5], behavioral BC [E6], Konkurrenzen ∉ l¹_{|L|} | numerisch belegt (c04) |

---

## 5. Rangliste der drei vielversprechendsten Techniken mit Angriffsplan

### Platz 1 — Diskrete Stabilität in l¹_{|L|,s} + affin-kovarianter Kantorovich in Kegelform (Vainikko/Keller, Deuflhard–Heindl, Krawczyk)

Warum: c01/c04/c08 zeigen, dass die abgeschnittene Linearisierung bei Radius |L| gut konditioniert ist (s = 0: |J_N⁻¹| ≤ 7.8
bis N = 194, σ_min ≈ 0.73; s = 4: |J_N⁻¹| ≤ 28.2, praktisch konstant) und Kneser dort konsistent ist (α → 0, c09); die
Konkurrenzen sind aus dem Raum ausgeschlossen. Das ist exakt das Schema „stabil + konsistent ⇒ konvergent“ — die affin-kovariante Form verhindert, dass der
Radiusverlust des Vorwärtsoperators (1.63^N) in β·γ explodiert.
Angriffsplan:
1. (D) Kantorovich auf affin-kovariante Form umstellen (h = αω), Norm l¹ mit |L|^k(1+k)^s, s ∈ [4, 5.45), in X **und** Y.
   Mit Worst-Case-ω trägt das numerisch bis N ≈ 70; mit arb als Satz für einzelne N (ergänzt das Krawczyk-Zertifikat).
2. (B/D) Engpass ω: Worst-Case-Paare (N/5, N) wachsen ≈ 1.25^N, entlang des tatsächlichen Zuwachses aber nur ω_eff ≈ N^{1.9}
   (c05, h_eff fällt). Also Kantorovich/Krawczyk in **Kegel-/Majorantenform** formulieren: statt der l¹-Kugel die Menge
   K_N(ε) = {δ : |δ_k| ≤ ε m_k |L|^{-k}} (komponentenweise Majorante, Rahmen [C10]) mit einem Profil m_k, das die grossen Spitzenanteile von p_N − q_N (≈ −Kneser-
   Koeffizienten bei k ≈ N) enthält, und zeigen: (i) J⁻¹F(q_N) ∈ K_N(ε/2), (ii) J⁻¹(F(q+δ) − F(q) − Jδ) ∈ K_N(ε/2) für
   δ ∈ K_N(ε) (Selbstabbildung des Kegels, Brouwer/Banach). Das ist die Struktur des vorhandenen Krawczyk-Zertifikats
   (Boxen |a_k − ã_k| ≤ 2^{-250}|L|^{-k}), nur mit analytisch statt numerisch kontrollierten Konstanten.
3. Strukturidentität nutzen: ∂_{a_i}C(p) = M_{x^i}·C(p)·D (D: e_j ↦ j e_{j−1}) und C(p∘q) = C(q)C(p); damit
   J⁻¹F''[e_i,·] ≈ (Multiplikation mit (p⁻¹)^i)∘Ableitung — sollte eine Schranke mit |f⁻¹(0)|^i = 0.696^i statt 1.25^N
   liefern, sofern der Toeplitz-Anteil kontrolliert wird. (**Vermutung**, nicht geprüft.)
4. Lemma (a) als eigenständige lineare Aussage formulieren und an E übergeben: „sup_N |J_N(q_N)⁻¹|_{l¹_{|L|,4}} < ∞“
   (numerisch ≤ 28.2). Für die linearen Modelle (Walker-Abel: 3.5 → 4.7 bei s = 0) dieselbe Aussage zuerst.

### Platz 2 — Abel-Koordinaten: Linearisierung = I + S_{1/2}, Kneser = Hauptlösung (Nörlund; Paulsen–Cowgill)

Warum: liefert die **unendlichdimensionale** a-priori-Schranke, ohne die Schritt (a) nicht zu beweisen ist. In w = slog(x)
ist J ≅ u ↦ u(w+½)+u(w) (Abschnitt 1). Beschränkte Lösungen auf Ω = slog(D) sind eindeutig (Bôcher, README); Existenz und
Schranke |u|_∞ ≤ C|g|_∞ wären der „Kneser-Teil“ (Kouznetsovs Cauchy-Integral bzw. Paulsen–Cowgills Iteration lösen genau
dieses lineare Problem).
Angriffsplan:
1. (B) Beweisen: für g beschränkt-holomorph auf Ω' = Ω ∩ (Ω−½) existiert genau eine auf Ω beschränkte Lösung, mit
   |u| ≤ C|g| (Riemann–Hilbert/Cauchy-Integral über den Rand von Ω, Phragmén–Lindelöf in vertikaler Richtung; Breite von
   Ω ≥ 2.067 > ½ ist die entscheidende Geometrie — README).
2. Zurück nach x: T'(A+½) ~ (x−L) nahe L; Normäquivalenz zwischen „beschränkt auf Ω“ und gewichteten Räumen auf D
   (H^∞ mit Gewicht |x−L|) klären; daraus J_∞⁻¹ beschränkt in einem l¹_{|L|}-nahen Raum.
3. Abschnitts-Stabilität aus der unendlichen Schranke + Kleinheit der Kopplung an den Schwanz (Stärke (a₀/|L|)^N ≈ 0.363^N,
   aus der Taylor-Verschiebung um a₀) gewinnen — klassisches Störargument (C5-Prinzip).

### Platz 3 — Minimallösung/Schwanzkorrektur (Miller–Olver–Pincherle, Fox-Beschleunigung, Flajolet–Odlyzko), erst linear

Warum: erklärt die Selektion (Kneser ist die „Minimallösung“, Trennung k^{-8.9}) und liefert die Werkzeuge für Konsistenz und
für ein besseres Startpolynom, das den Engpass (c) entschärft. **Neu durch Strang B (BRETT 19:53):** im Korridor ist p_N bis
auf O(c_N²) eine exakte halbe e-Funktion f_{c_N} = Kneser + **ein** Störmodus c_N e^{2πiw} + c̄_N e^{−2πiw}, mit
c_N ≈ (0.18 − 0.08i)·N^{1−2γ̄}. Das ist genau die Fehlerform von Millers Algorithmus (Gautschi 1967; Formel aus dem Gedächtnis [U]):
berechnete Lösung = Minimallösung − (Verhältnis minimal/dominant am Abschneideindex) × dominante Lösung. Damit wird ein
**Lyapunov–Schmidt-Ansatz** konkret: p = T_N(f_c) + r mit komplexem Parameter c und Rest r im Komplement des Modus; die
abgeschnittene Gleichung zerfällt in ein stabiles Restsystem (Platz 1/2) und **eine** skalare Gleichung G_N(c) = 0, deren
Lösung c_N ~ N^{1−2γ̄} die Rate liefert (G_N(0) = Kopplung des Kneser-Schwanzes, G_N'(0) = Kopplung der Konkurrenz-Tangente).
Literatur: Lyapunov–Schmidt-Reduktion z. B. Golubitsky–Schaeffer, Singularities and Groups in Bifurcation Theory I (1985) [U];
Gautschi/Olver für die Fehlerform (Quellen [S], Formel [U]).
Angriffsplan:
1. (E) Lineares Modell zuerst: Paulsens abgeschnittene Julia-Gleichung λ(eˣ) = eˣλ(x) (konvergiert numerisch schnell, gleicher
   Operatortyp) und Walkers Abel-System (langsam, N^{-2.4}). Ziel: Satz „Abschnitt → Kneser“ für das lineare Modell via
   Stabilität in l¹_{|L|,s} + Konsistenz (Abel-System: Kneser-slog ~ k^{-1}|L|^{-k} liegt NICHT in l¹_{|L|} — dort
   Konsistenz nur in schwächerer Norm, daher langsam); Trappmanns Technik (A3) für explizit lösbare Spielzeuge (bx, Möbius).
2. (A) Lokale Entwicklung bei L (B8) mit c_{10} = 0 bestätigen; daraus mit Singularitätsanalyse (E9) die Asymptotik des
   Kneser-Schwanzes und von F_N(q_N) (Kopplung über a₀: Grösse ~ N^{-6.45}((a₀+a₁s)/|L|)^N in l¹_s).
3. (D) Startpolynom q̃_N = q_N + (Schwanzantwort aus der Singularitätsentwicklung, Fox' Idee) — verkleinert α in den hohen
   Graden, wo ω gross ist; dann Kantorovich wie Platz 1.

---

## 6. Warnungen und Sackgassen

- Numerische Übereinstimmung von Abschneide-Verfahren mit einem Kandidaten kann täuschen (A8: √2, 1e-6; Forum „Bummer!“:
  reguläre Iterationen an zwei Fixpunkten von √2^x unterscheiden sich um 1e-24).
- Walkers c₁^{(N)} wirkte 1991 wie „vielleicht eine andere Lösung“ (Diskrepanz zu seiner C^∞-Lösung in der 4. Stelle) — beide
  Lösungen sind verschieden, der Matrix-Grenzwert ist Kneser (3.6).
- Paulsens „7–8 Stellen, dann nicht weiter“ beim Abel-System (A4) ist mit langsamer algebraischer Konvergenz verträglich —
  kein Beleg für Divergenz.
- Transferoperator-Rahmen (2.4) und Nash–Moser (2.6) sind am Kern vorbei; Regularisierung/Picard (2.5) wird durch die Wahl
  r = |L| überflüssig.

## 7. Offene Literaturfragen / Anfragen

- Walker JMAA 155 (1991) und Clenshaw–Lozier–Olver–Turner (1986) im Volltext: enthalten sie Konvergenzaussagen für den
  Matrix-Abschnitt? (bisher nur [S]).
- Kantorovich–Akilov: genaue affin-kovariante Formulierung zitierfähig nachschlagen.
- Gibt es Finite-Section-Sätze für „Komposition + Multiplikation“ in gewichteten l¹-Räumen (Carleman-artig, nicht
  bandbeherrscht)? Suche bisher erfolglos.

## 8. Gezielte Antworten an die anderen Stränge

- **A:** Levensteins Form der Kneser-Entwicklung bei L (B8, Forum tid 1043) prüfen, insbesondere ob der Koeffizient von
  (z−L)^p verschwindet (Leitterm (z−L)^{p+1} ⇔ Koeffizienten ~ k^{-6.45}). Für Schwanz-Asymptotik: Flajolet–Odlyzko (E9).
- **B:** (i) Die Linearisierung ist in Abel-Koordinaten u ↦ u(w+½)+u(w) (Abschnitt 1, formal bewiesen); Schablone für die
  beschränkte Lösung: Loray–Touzet–Voronin Thm 9.5 (G2, gelesen). (ii) Ein-Moden-Befund (BRETT 19:53) = Fehlerform des
  Miller-Algorithmus (E5) ⇒ Lyapunov–Schmidt mit Parameter c (Platz 3). (iii) Konvergenzgebiet D ∪ D(a₀, |L−a₀|) (3.5);
  Raten an der einspringenden Ecke bei L: konstruktive Funktionentheorie (Dzjadyk, Andrievskii–Belyi, Gaier) [U].
- **D:** Normen mit **gleichen** Gewichten in Urbild und Bild, Radius |L|, Polynomgewicht s ≈ 4 (3.1); Kantorovich in
  affin-kovarianter Form (C4) und Kegelform (C10, Rump); worst-case ω ist exponentiell, ω entlang des Zuwachses nicht (3.4).
- **E:** Lineare Modelle: Walker-Abel-System (in l¹_{|L|} stabil, aber Kneser-slog nicht im Raum ⇒ langsam) und Paulsens
  Julia-System (schnell, Kneser auf 18 Stellen); Trappmanns explizite Technik (A3) für f = bx; für Möbius-Modelle Trappmann,
  „The iterational paradise of fractional linear functions“, Forum tid 1608 (quellen.yaml: geprüft original).

---

# Runde 2 (Auftrag der Leitung 20:00): Beweise für (b) Konsistenz und (c) Nichtlinearität

Bezeichnungen: ρ = |L|, a₀ = f_K(0), R₀ = |L − a₀| = 1.34935, q = q_N = [f_K]_{≤N}, t = t_N = f_K − q_N (Schwanz),
A = A_N = J_N(q_N), B_N(p) = ([p^k]_{≤N})_{k>N} (Strang D), |v|_s = Σ_j |v_j| ρ^j (1+j)^s, |·| = |·|_0.
δ1 := −A⁻¹F_N(q_N) (erster Newton-Schritt), α_N := |δ1|.  W := D(0,ρ) ∪ D(a₀,R₀).

## 9. (b) Konsistenz

### 9.1 Negativbefunde (numerisch, entscheiden über den Beweisweg)

| N | |F_N(q_N)|_0 | |F_N(q_N)|_4 | α_N (s=0) | α_N (s=4) | worst-case Antwort sup_k |A⁻¹b_k|/ρ^k | α_N/|t_N|_0 | α_N/|t_N|_4 |
|---|---|---|---|---|---|---|---|
| 20 | 1.7e-3 | 2.4e2 | 1.27e-6 | 0.127 | 3.1e2 (k=29) | 1.43 | 0.22 |
| 40 | 0.34 | 7.5e5 | 5.06e-8 | 0.071 | 5.2e5 (k=58) | 2.77 | 0.40 |
| 80 | 7.7e4 | 2.5e12 | 1.94e-9 | 0.033 | 1.4e12 (k=116) | 5.07 | 0.56 |
| 120 | 7.4e12 | 1.4e21 | 2.81e-10 | 0.021 | 4.2e18 (k=174) | 6.84 | 0.68 |
| 160 | 4.5e20 | 2.7e29 | 7.04e-11 | 0.015 | 1.2e25 (k=232) | 8.36 | 0.77 |

(c10, c11.) Das Residuum ist in jeder Gewichtung mit Radius ρ **exponentiell gross** (Masse in den obersten Graden), die
Antwort auf einzelne Schwanzkoeffizienten ebenso (Spitze bei k ≈ 1.45N) — trotzdem ist α_N klein und ≈ |t_N|·O(1).
**Folgerung:** (b) folgt **nicht** aus (a) + Koeffizientenschranken; jede Abschätzung der Form |A⁻¹|·|F| oder
Σ_k |A⁻¹b_k||f_k| ist um Faktoren bis 1e30 zu grob. Die Kleinheit von α_N ist eine **Kohärenz**-Eigenschaft des
Kneser-Schwanzes (Schwanz einer Funktion mit bekannter analytischer Fortsetzung). Strang E (BRETT 20:08) sieht dasselbe
im Möbius-Modell. Zusätzlich (D, 20:15): (a) gilt in l¹ nur bis auf log N; in l²_{ρ,s} gleichmässig — für die folgenden
Schranken unerheblich, da alle hier in l¹ gemessen sind und |·|_{l²} ≤ |·|_{l¹}.

### 9.2 Exakte Identitäten (bewiesen)

**Lemma 9.1 (Schwanzkopplung; = Strang D, Lemma D1).** Für p ∈ P_N und h = Σ_{k>N} h_k x^k mit
Σ_k |h_k| r^k < ∞ für ein r > |p(0)| gilt [(p+h)∘(p+h)]_{≤N} = [p∘p]_{≤N} + B_N(p)h, B_N(p)h = [h∘p]_{≤N}.

**Lemma 9.2 (Konsistenz-Residuum).** F_N(q_N) = −[t_N∘q_N]_{≤N} = −B_N t_N, also δ1 = A⁻¹B_N t_N.
*Beweis.* q∘q = f∘q − t∘q; f∘q − f∘f = f∘(f−t) − f∘f = O(x^{N+1}) (f holomorph bei a₀, t = O(x^{N+1})); f∘f = exp. ∎

**Lemma 9.3 (Pol-Kopplung in geschlossener Form).** Für |ζ| > |a₀|:
B_N[(ζ−x)⁻¹]_{>N} = ζ^{−(N+1)}·[q^{N+1}/(ζ − q)]_{≤N} =: ζ^{−(N+1)} v_ζ.
*Beweis.* (ζ−x)⁻¹_{>N} = Σ_{k>N} ζ^{−k−1}x^k, also B_N(·) = Σ_{k>N} ζ^{−k−1}[q^k]_{≤N}. Mit der Majorante Q (Beträge der
Koeffizienten von q) ist Σ_k |ζ|^{−k−1}[Q^k]_j = [Q^{N+1}/(|ζ|^{N+1}(|ζ|−Q))]_j < ∞, weil Q(0) = |a₀| < |ζ|; also darf man
summieren: Σ_{k>N} q^k ζ^{−k−1} = q^{N+1}ζ^{−N−1}/(ζ − q) als Potenzreihe in x. ∎

**Definition (Pol-Antwort).** R_N(ζ) := |A⁻¹B_N[(ζ−x)⁻¹]_{>N}| = |A⁻¹v_ζ| / |ζ|^{N+1}  (|ζ| > |a₀|). A reell ⇒ R_N(ζ̄) = R_N(ζ).

**Lemma 9.4 (Cauchy-Darstellung der Antwort).** Sei Γ ein stückweise glatter Weg in {|ζ| ≥ r₁}, r₁ > |a₀|, μ messbar mit
∫_Γ |μ(ζ)||ζ|^{−N−2}|dζ| < ∞, und h_k = (1/2πi)∫_Γ μ(ζ)ζ^{−k−1}dζ für alle k > N. Dann
A⁻¹B_N h_{>N} = (1/2πi)∫_Γ μ(ζ)·A⁻¹B_N[(ζ−x)⁻¹]_{>N} dζ  und  |A⁻¹B_N h_{>N}| ≤ (1/2π)∫_Γ |μ(ζ)| R_N(ζ) |dζ|.
*Beweis.* Koordinate j von B_N h_{>N} ist Σ_{k>N} h_k [q^k]_j. Einsetzen und Vertauschen von Σ und ∫ (Fubini mit der Majorante
aus 9.3, gleichmässig in |ζ| ≥ r₁) gibt (1/2πi)∫μ(ζ)Σ_k ζ^{−k−1}[q^k]_j dζ = (1/2πi)∫μ(ζ)(B_N[(ζ−x)⁻¹]_{>N})_j dζ. A⁻¹ ist eine
feste Matrix; Norm des Integrals ≤ Integral der Norm. ∎

### 9.3 Satz 9.5 (Konsistenz-Reduktion) — bewiesen unter den genannten Voraussetzungen

**Voraussetzungen.**
- (F) [analytisch, über f_K] Es gibt u₀ > 0, einen geschlossenen Weg Γ' ⊂ ℂ ∖ W̄ von L(1+u₀) nach L̄(1+u₀) und ein Gebiet Ω,
  berandet von den beiden Ufern der radialen Schlitze S_L = {L(1+u): 0 ≤ u ≤ u₀}, S_{L̄} und Γ', mit D̄∖{L, L̄} ⊂ Ω, so dass f_K
  holomorph auf Ω ist, stetig bis an beide Schlitzufer, |f_K| ≤ M auf Γ', und für den Sprung über S_L gilt
  |Δf(L(1+u))| ≤ C_f u^{σ}, σ = Re(1+γ) = 5.447 (0 ≤ u ≤ u₀).
- (H_R1) [linear] R_N(L(1+u)) ≤ C₁ N^{a} e^{−κ N^{a} u} für 0 ≤ u ≤ u₀ (a, κ > 0).
- (H_R2) [linear] sup_{ζ∈Γ'} R_N(ζ) ≤ C₂ e^{−κ' N}.

**Satz 9.5.** Unter (F), (H_R1), (H_R2) gilt für alle N
  α_N = |A_N⁻¹F_N(q_N)| ≤ (|L|/π)·C_f C₁ Γ(σ+1) κ^{−σ−1} · N^{−aσ} + (M·|Γ'|/2π)·C₂ e^{−κ'N}.
Insbesondere α_N → 0 polynomiell; mit a = 0.883 (Strang B: Aussenwinkel der Ecke von W bei L) ist der Exponent
aσ = 0.883·5.447 = 4.81 — genau die gemessene Rate von α_N (c14: N^{−4.8}) und von |p_N(L) − L| (B: N^{−4.83}).

*Beweis.* (i) Für r < ρ ist f_k = (1/2πi)∮_{|ζ|=r} f(ζ)ζ^{−k−1}dζ. Nach (F) ist f auf Ω holomorph; Ω enthält den Kreis
|ζ| = r und das Gebiet zwischen ihm und ∂Ω. Cauchy: ∮_{|ζ|=r} = ∮_{∂Ω}; die kleinen Kreisbögen um L, L̄ verschwinden, da f bis
L stetig (beschränkt) ist und ihre Länge → 0. Auf den Schlitzen bleibt das Integral über den Sprung Δf, auf Γ' das über f.
Für jedes k > N ist das Integral absolut konvergent (|ζ| ≥ min(ρ, dist) > |a₀| auf ∂Ω). (ii) Lemma 9.4 mit μ = Δf auf
S_L ∪ S_{L̄} und μ = f auf Γ' und h = t_N (dessen Koeffizienten genau die f_k, k > N, sind): α_N = |A⁻¹B_N t_N| (Lemma 9.2)
≤ (1/2π)[2∫_0^{u₀}|Δf(L(1+u))|R_N(L(1+u))|L|du + ∫_{Γ'}|f|R_N|dζ|] (Faktor 2: konjugierter Schlitz, R_N(ζ̄) = R_N(ζ),
|Δf| symmetrisch). (iii) Einsetzen von (H_R1): ∫_0^∞ u^σ C₁N^a e^{−κN^a u}du = C₁Γ(σ+1)κ^{−σ−1}N^{a−a(σ+1)} = …N^{−aσ};
(H_R2) für Γ'. ∎

**Was damit bewiesen ist:** (b) ist **vollständig** auf zwei Arten von Aussagen zurückgeführt: (F) reine Eigenschaften von
f_K (Holomorphie auf einem festen Gebiet um W̄ mit radialen Schlitzen, Sprung ~ u^{5.447} — letzteres folgt aus Strang As
lokaler Entwicklung f_K − L = √L t(…) + Σ d_n t^{1+nγ}(…), deren Sprung über den Schlitz der Leitterm d₁ t^{1+γ} bestimmt),
und (H_R1/2) eine **skalare** Aussage über die abgeschnittene Linearisierung: die Pol-Antwort R_N(ζ). (H_R) ersetzt die
falsche Hoffnung „|F| → 0“ und ist das genaue lineare Gegenstück zu (a) für die Konsistenz.

### 9.4 Numerische Prüfung der Voraussetzungen von Satz 9.5 (c14, c15)

**(H_R1) — Skalengesetz am Schlitz.** R_N(L(1+v/N^{0.883}))/N^{0.883}:

| N | v=0.01 | v=0.1 | v=0.3 | v=1 | v=2 | v=4 | v=8 | v=16 |
|---|---|---|---|---|---|---|---|---|
| 40 | 0.394 | 0.347 | 0.262 | 0.0989 | 0.0255 | 0.00189 | 1.5e-5 | 3.4e-9 |
| 80 | 0.484 | 0.424 | 0.318 | 0.117 | 0.0286 | 0.00184 | 9.5e-6 | 5.2e-10 |
| 120 | 0.540 | 0.473 | 0.352 | 0.127 | 0.0301 | 0.00180 | 7.6e-6 | 2.2e-10 |
| 160 | 0.581 | 0.508 | 0.377 | 0.134 | 0.0310 | 0.00177 | 6.6e-6 | 1.3e-10 |
| 194 | 0.610 | 0.532 | 0.393 | 0.138 | 0.0316 | 0.00175 | 6.1e-6 | 9.5e-11 |
| 250 | 0.648 | — | 0.416 | 0.144 | — | 0.00172 | 5.5e-6 | — |
| 300 | 0.677 | — | 0.432 | 0.148 | — | 0.00170 | 5.1e-6 | — |

Die Kurven fallen zusammen (bei kleinem v schwacher Anstieg ≈ N^{0.27}, für das Integral unerheblich); c18 bei N = 250, 300:
unverändert, obwohl Strang D dort die exponentiell wachsende Eckspitze von J_N⁻¹ sieht — die Pol-Antwort (kohärente Eingabe)
sieht sie nicht. Konturpunkte (Abstand 0.1 von W): ζ = −0.01+1.47i: 8.8e-6 / 5.8e-7, ζ = 1.95: 2.6e-10 / 3.2e-12, ζ = −1.47:
1.2e-6 / 6.5e-8 (N = 250 / 300). α_N = 8.07e-12 (250), 3.39e-12 (300), lokaler Exponent −4.90 bzw. −4.76. Hüllkurve z. B.
R_N(L(1+u)) ≤ 0.7·N^{0.883}·e^{−1.15·N^{0.883}u} für alle Messpunkte. Der Exponent 0.883 ist Strang Bs Aussenwinkel der Ecke von
W bei L (π·0.883 = 158.9°); mit einem Exponenten 1 kollabieren die Daten deutlich schlechter.

**(H_R2) — Kontur ausserhalb W̄** (Rand von W um d nach aussen, ohne |ζ−L| < 0.25): max R_N

| d | N=40 | N=80 | N=120 | N=160 | N=194 |
|---|---|---|---|---|---|
| 0.05 | 1.61 | 0.94 | 0.41 | 0.165 | 0.073 |
| 0.10 | 0.43 | 0.074 | 0.0097 | 0.0012 | 1.9e-4 |
| 0.20 | 0.042 | 7.8e-4 | 1.2e-5 | 1.6e-7 | 4.1e-9 |

**Wachstumsgesetz (Vermutung R, numerisch):** R_N(ζ)^{1/N} → max(ρ/|ζ|, R₀/|ζ−a₀|) =: 1/r_W(ζ) für ζ ∉ W̄ (Zerfall) und
für ζ ∈ W (Wachstum); Beispiele: ζ = 2.2 vorhergesagt 0.793, beobachtet Faktor 9.2e-5 je 40 Grade (≈ 0.793^40); ζ = −2.2
vorh. 0.625, beob. ≈ 0.65; ζ = 1.5 vorh. Wachstum 1.35/Grad, beob. ≈ 1.3; Vorfaktoren ~N^{3…5}. W ist also genau das Gebiet,
auf dem p_N konvergiert (3.5) — das Abschneiden „sieht“ f auf W. **[Runde 4: in x nur bis N ≈ 600–900, §3.5; als Gebiet der Polvariable ζ richtig.]**

**Gittertest von Vermutung R (c17, N = 80 → 160, 44 Punkte in −2.2 ≤ Re ζ ≤ 2.4, 0 ≤ Im ζ ≤ 2.4):**
log(beobachtete Rate) − log(1/r_W): Mittel +0.016, Streuung 0.039, min −0.12, max +0.07. Ausserhalb W stimmt der Zerfall überall
(z. B. ζ = 1.9: 0.963 vorh./0.945 beob.; ζ = −1.6: 0.859/0.874; ζ = 2.4+2.4i: 0.441/0.464), innerhalb W das Wachstum
(ζ = −0.4: 3.436/3.433; ζ = 0.2+0.6i: 2.173/2.173; ζ = −1+0.6i: 1.179/1.185). Grösste Abweichung im Inneren beider Kreise auf ℝ⁺
(ζ = 0.8: 4.48/3.96, ζ = 1.4: 1.50/1.42). Der systematische Vorlauf +0.016 entspricht einem Vorfaktor ≈ N^{1.9}.
**Korrektur (Runde 4, §15):** Der Vorlauf ist kein Vorfaktor, sondern genau log(e^{−g_E}/Max): ausserhalb W̄ gilt die Green-Funktions-
Rate e^{−g_E(ζ)} (Abweichung −0.0035 ± 0.0091 an den 32 Aussenpunkten statt +0.032 ± 0.021 für die Max-Formel), mit den gemessenen
Randdaten über Poisson–Jensen sogar −0.0008 ± 0.0004 (c29). Die Max-Formel ist für ζ ∉ W̄ falsch (exponentiell, Strang E).

**Schranke gegen Wahrheit** (Leitterm D₁ des Sprungs, Schlitzintegral bis |ζ−L| = 0.25): I_N/α_N = 1.34, 2.28, 2.20, 2.12,
2.06 (N = 40, 80, 120, 160, 194) — die Dreiecksungleichung im Integral verliert nur einen Faktor ≈ 2. Und die Antwort auf
D₁(1−x/L)^{1+γ} + konj allein trifft δ1 bis auf 0.37, 0.18, 0.12, 0.089 (relativ, ~1/N: nächste Terme d₂, e_{1,j}).

**Möbius-Gegenprobe (Strang E):** f_M = tan(arctan x + c/2), Fixpunkte ±i, a₀ = τ = tan(c/2), R₀ = √(1+τ²), Pol 1/τ.
Pol ∈ D(a₀,R₀) ⇔ (1−τ²)² < τ²(1+τ²) ⇔ τ > 1/√3 ⇔ c > π/3 — genau Es Schwelle (f³(0) = ∞, keine Lösung nahe f_M); und
R₀/|Pol − a₀| = τ√(1+τ²)/(1−τ²) = sin(c/2)/cos c ist genau Es gemessene Rate. Satz 9.5 + Vermutung R sagen also Es Modell
vollständig vorher (bewiesen ist dort die Geometrie, die Rate ist numerisch).

**(F) numerisch bestätigt (Strang A, BRETT 20:59, a_moden/k25):** Sprung über den radialen Schlitz |Δf|/u^{5.447} = 0.938 …
0.760 für u = 0.005 … 0.2 (Leitterm 0.943, also C_f ≤ 0.95); max|f_K| auf ∂(W + 0.1) ohne |ζ−L| < 0.25: 3.32 (d = 0.2: 3.54);
keine weiteren Singularitäten nahe dem Kreisrand, Hauptzweig von slog nur mit L, L̄ (nächste bei L + 2πik). Damit sind alle
Konstanten von Satz 9.5 bis auf (H_R) belegt. **Explizite Schranke:** Schlitzteil I_N (oben) plus Konturteil ≈ (max|f|·Länge/2π)
·max_{Γ'}R_N ≈ 6·max R_N; für d = 0.2: 0.25, 4.7e-3, 7e-5, 1e-6, 2.5e-8 (N = 40 … 194). Der Konturteil zerfällt exponentiell
(≈ e^{−0.11N}), dominiert aber die Schranke bis N ≈ 270; danach bestimmt der Schlitzteil die Rate N^{−4.81}. (Die Wahrheit α_N
liegt schon ab N = 40 beim Schlitzteil: die Konturbeiträge löschen sich offenbar weitgehend aus.)

**Zusammenfassung (b):** bewiesen ist die Reduktion „(F) + (H_R) ⇒ α_N ≤ C N^{−4.81}“ (Satz 9.5); (F) betrifft nur f_K (von
Strang A mit dem Motor prüfbar: Holomorphie auf einem Gebiet um W̄ mit radialen Schlitzen, Schranke M auf Γ'); (H_R) ist eine
skalare lineare Aussage, numerisch für N ≤ 160 mit Skalengesetz und exponentiellem Zerfall belegt. Offen: Beweis von (H_R) —
vermutlich über Ds Grenzoperator (Toeplitz mit Sprungsymbol) + Bs Mellin-Symbol an der Ecke (lokales Problem bei L).

## 10. (c) Nichtlinearität

### 10.1 Exakte Zerlegung des nichtlinearen Rests (bewiesen)

**Lemma 10.1.** Für δ ∈ P_N gilt mit J_∞(f)δ := δ∘f + (f'∘f)δ:
  F_N(q+δ) = [(f+δ)∘(f+δ) − exp]_{≤N} − [t∘(q+δ)]_{≤N},
  A δ = [J_∞(f)δ]_{≤N} − [(t'∘q)δ]_{≤N},
  N_q(δ) := F_N(q+δ) − F_N(q) − Aδ = [Q_f(δ)]_{≤N} − [Q_t(δ)]_{≤N},
mit Q_f(δ) = (f+δ)∘(f+δ) − f∘f − J_∞(f)δ (Rest zweiter Ordnung des **exakten** Problems an seiner Lösung) und
Q_t(δ) = t∘(q+δ) − t∘q − (t'∘q)δ (Taylor-Rest zweiter Ordnung des Schwanzes längs q).
*Beweis.* (q+δ)∘(q+δ) = f∘(q+δ) − t∘(q+δ) + δ∘(q+δ); modulo x^{N+1} ist f∘(q+δ) ≡ f∘(f+δ) und δ∘(q+δ) ≡ δ∘(f+δ), weil
(q+δ) − (f+δ) = −t = O(x^{N+1}) und f, δ bei a₀ + δ₀ holomorph sind. Zweitens: (f'∘f) − (q'∘q) = (f'∘f − f'∘q) + t'∘q ≡ t'∘q
mod x^{N+1}. Einsetzen; die Terme erster Ordnung in t heben sich (t∘(q+δ) − t∘q − (t'∘q)δ). ∎
Folge: Die Tail-Terme erster Ordnung fallen exakt weg; der nichtlineare Rest besteht aus dem Rest des wahren Problems und einem
Rest zweiter Ordnung des Schwanzes (über Lemma 9.3 wieder durch Pol-Kerne darstellbar:
[Q_t(δ)]_{≤N} = (1/2πi)∫ μ(ζ){B_N(q+δ) − B_N(q) − ∂B_N(q)[δ]}[(ζ−x)⁻¹]_{>N} dζ).

### 10.2 Vorzeichen-unempfindliche Kantorovich-Formen scheitern für grosse N (numerisch)

| Ansatz | N = 20 | 40 | 60 | 80 | 100 | Wachstum |
|---|---|---|---|---|---|---|
| affin, ganze l¹_ρ-Kugel um q (c02/c03), h = αω | 2.5e-4 | 1.8e-3 | 4.1e-2 | — (N=70: 0.20) | — | ω ~ 1.25^N |
| Kastennorm um q, Profil |δ1| (c12, P1) | 3.0e-5 | 6.9e-5 | 6.1e-4 | 9.3e-3 | 0.33 | ≈ 1.19^N |
| Kastennorm um x₁ = q+δ1, Profil |δ2| (c13, P1) | 2.7e-12 | 3.6e-13 | 2.3e-12 | 1.3e-11 | 2.5e-10 | ≈ 1.13^N |
| entlang des echten Zuwachses (c05/c09), α·ω_eff | 1.6e-7 | 9.5e-9 | 4.8e-9 | 2.3e-9 | 1.9e-9 | fällt |

Diagnose (c16, N = 60, Zentrum x₁): die Kastensumme wird von Paaren (i, j) ≈ (N/4, N) getragen, d. h. von den obersten
Koeffizienten; jede Newton-Korrektur (δ1, δ2, …) hat einen Anteil vergleichbarer Grösse in den obersten Graden. Weitere
Newton-Schritte als Zentrum verschieben die Schwelle nur: am Zentrum x_k ist α_k ≈ (κα)^{2^k}/κ polynomiell klein, ω_Kasten
wächst aber exponentiell (x₁: 2.5e-10, 2.7e-9, 8.0e-8, 1.2e-6 für N = 100, 120, 140, 160, zuletzt ≈ 1.19^N; Schwelle h = ½
bei N ≈ 230, c13b); für alle N
bräuchte man k ~ log₂ N Schritte, jeder einzelne analytisch (kohärent) kontrolliert. **Folgerung:** Mit
Abschätzungen, die nur Beträge der Koeffizienten benutzen (Kugel, Kasten, Kegel aus Majoranten), ist (c) für alle N **nicht**
zu haben; im Zertifikatsbereich (N ≤ 100 numerisch, Strang D/Krawczyk bis 120) genügt es.

### 10.3 Was für (c) nötig ist (Programm, nicht bewiesen)

Der Rest muss **kohärent** behandelt werden wie in (b): Man braucht eine Klasse C_N von Zuwächsen, die (i) δ1 enthält (Satz 9.5:
δ1 = (1/2πi)∫ μ(ζ) A⁻¹B_N[(ζ−x)⁻¹]_{>N} dζ), (ii) unter δ ↦ −A⁻¹N_q(δ) in sich abgebildet wird und (iii) eine Norm trägt,
in der Lemma 10.1 quadratisch abgeschätzt werden kann. Kandidat: C_N = Bild der Pol-Antworten mit Dichten μ auf Γ (Norm
∫|μ| R_N |dζ|) — dann ist [Q_t(δ)] wegen 10.1 wieder ein Pol-Integral (zweiter Ordnung), und [Q_f(δ)] muss über die Lösung V
des **wahren** linearen Problems J_∞(f)V = Q_f(δ) auf Ω zurückgeführt werden (A⁻¹[J_∞V]_{≤N} = V_{≤N} + A⁻¹B_N V_{>N} + klein,
exakte Identität aus 9.2). Dazu fehlen: (α) die Lösbarkeit von J_∞(f)V = H mit Schranken auf Ω (Abel-Koordinaten, Platz 2,
Loray–Touzet–Voronin-Schablone G2), (β) eine bilineare Pol-Antwort-Schranke (H_R^{(2)}). Strang Bs Ein-Moden-Befund
(p_N ≈ f_{c_N} im Inneren) ist die numerische Bestätigung, dass die Newton-Folge in einer solchen Klasse bleibt.

### 10.4 Satz 10.4 (Sekanten-Darstellung des Fehlers) — bewiesen, a posteriori

**Satz 10.4.** Sei p ∈ P_N eine Lösung von F_N(p) = 0 (z. B. p_N) und Ā := ∫_0^1 J_N(q + s(p − q)) ds invertierbar. Dann gilt
exakt p − q = Ā⁻¹B_N t_N, und für jede Cauchy-Darstellung des Schwanzes wie in Satz 9.5
  |p − q| ≤ (1/2π)∫_Γ |μ(ζ)| R̄_N(ζ) |dζ|,   R̄_N(ζ) := |Ā⁻¹B_N[(ζ−x)⁻¹]_{>N}|.
*Beweis.* F_N ist ein Polynom in den Koeffizienten, also F(p) − F(q) = Ā(p − q) (Hauptsatz entlang der Strecke). Mit F(p) = 0
und F(q) = −B_N t (Lemma 9.2) folgt p − q = Ā⁻¹B_N t; Rest wie Lemma 9.4. ∎
**Bedeutung.** Gibt es den Zweig p_N (numerisch bis N = 1200 bei Strang D, mit Zertifikat bis 120), so ist die Vermutung
**äquivalent** zu einer linearen Antwort-Schranke für den Sekantenoperator Ā_N — die Nichtlinearität (c) ist dann vollständig
in Ā versteckt, kein Kantorovich-Schritt nötig. Offen bleibt, Ā⁻¹ ohne Kenntnis von p zu kontrollieren; Strang D (BRETT 20:41)
zeigt, dass J(x_s)⁻¹ längs der Strecke linear von J(p)⁻¹ (stabil) zu J(q)⁻¹ (für N ≳ 190 instabil) übergeht — also ist das
**Zentrum q_N für grosse N ungeeignet**, und auch Ā erbt einen Teil davon. Folgerung für die Beweisarchitektur: Konsistenz
(Satz 9.5) und Nichtlinearität sollten an einem kohärenten Zentrum nahe p_N formuliert werden (erstes Newton-Iterat
x₁ = q + δ1 oder die Ein-Moden-Lösung f_{c_N} von Strang B), nicht an q_N.

### 10.5 Stand (c)

- **Bewiesen:** Lemma 10.1 (exakte Zerlegung des Rests), Satz 10.4 (Sekanten-Darstellung), das Kantorovich-Kriterium in
  beliebiger Norm (klassisch).
- **Numerisch:** am Zentrum x₁ = q + δ1 hält die Kantorovich-Bedingung in der Kastennorm mit h ≤ 1.2e-6 für N ≤ 160 (c13/c13b),
  Schwelle bei N ≈ 230; entlang des echten Weges ist die Nichtlinearität winzig (h_eff ≤ 1e-8, c05/c09).
- **Nicht erreichbar mit Betragsabschätzungen:** ein Beweis für alle N (Abschnitt 10.2). Für (c) „für alle N“ wird ein
  kohärentes Argument gebraucht (10.3): Newton- oder Fortsetzungsweg in einer Klasse von Zuwächsen, die als Pol-Antworten mit
  kontrollierter Dichte darstellbar sind, oder ein Ljapunow–Schmidt-Ansatz um die exakte Familie f_c (Strang B), bei dem die
  Nichtlinearität längs der Familie exakt ist.

## 11. Gemeinsame lineare Aussage: Vermutung Q (starke Konvergenz auf regulären Eingaben)

Strang D (BRETT 20:53) findet, dass sup_N |J_N⁻¹| in keiner getesteten Norm beschränkt ist (Eckspitze bei L wächst ≈ 1.012^N).
Satz 9.5 ist davon **nicht** betroffen: er benutzt keine Operatornorm von J_N⁻¹, nur die skalare Pol-Antwort R_N(ζ) auf einer
festen Kontur. Die passende lineare Aussage ist daher eine **starke** (eingabeabhängige) Konvergenz:

**Vermutung Q.** Mit J_∞(f)v := v∘f + (f'∘f)v und y_N(v) := J_N(q_N)⁻¹[J_∞(f)v]_{≤N} − [v]_{≤N} gilt für jedes v, das auf
einer Umgebung von W̄ holomorph ist (oder dort stetig mit Hölder-Verhalten bei L, L̄):
  |y_N(v)| ≤ C N^p (‖v − T⁰_N v‖_{D} + ‖v − T^{a₀}_N v‖_{D(a₀,R₀)}),
T⁰_N, T^{a₀}_N = Taylor-Polynome vom Grad N bei 0 bzw. a₀.

Folgerungen (bewiesen, falls Q gilt): (i) Für v = (ζ−x)⁻¹ ist y_N(v) = A⁻¹B_N v_{>N} + A⁻¹[(t'∘q)v_{≤N}] (exakte Identität
[J_∞v]_{≤N} = Av_{≤N} + B_N v_{>N} + [(t'∘q)v_{≤N}]_{≤N}, aus Lemma 9.1 und 10.1), und die Taylor-Reste sind ~ (ρ/|ζ|)^N bzw.
(R₀/|ζ−a₀|)^N ⇒ Gesetz R (§9.4). (ii) Für v = t_N (Kopf 0) ist y_N(t_N) = δ1; der Taylor-Rest bei 0 ist t_N selbst, der bei a₀
ist f − T^{a₀}_N f, weil q_N als Polynom vom Grad N von T^{a₀}_N exakt reproduziert wird — **das erklärt, warum die riesigen
Werte von q_N ausserhalb von D die Konsistenz nicht stören.** Beide Reste sind O(N^{−5.45}) (Singularität L auf beiden
Kreisrändern) ⇒ α_N ≤ C N^{p−5.45}; gemessen p ≈ 0.64. (iii) Möbius-Modell (Strang E, 20:52): Max-Gesetz in vier Lagen bestätigt.

Q ist die natürliche Form von (a) für dieses Problem: kein „sup |J_N⁻¹| < ∞“, sondern „J_N⁻¹ reproduziert reguläre Funktionen
auf dem Zwei-Kreis-Gebiet W mit Taylor-Genauigkeit“. Beweisidee (offen): Ds Grenzoperator (Toeplitz mit Sprungsymbol auf D)
plus Bs Mellin-Analyse der Ecke bei L; die Ecke (Aussenwinkel 0.883π) liefert den Vorfaktor N^p und die Skala N^{−0.883}.

---

# Runde 3 (Auftrag 21:05): Beweis von (H_R) über Gesetz R

## 12. (H_R) folgt aus einer einzigen Randgrösse: Rationalität in ζ + Bernstein–Walsh

### 12.1 Lemma 12.1 (Pol-Antwort ist rational in ζ) — bewiesen

Mit g := q − a₀ (g(0) = 0) und z_m := A⁻¹[q^{N+1}g^m]_{≤N} (m = 0 … N, feste Vektoren) gilt für alle ζ ∉ {0, a₀}
  y_ζ := A⁻¹B_N[(ζ−x)⁻¹]_{>N} = ζ^{−N−1} Σ_{m=0}^{N} z_m (ζ − a₀)^{−m−1}.
*Beweis.* Lemma 9.3 gibt B_N[(ζ−x)⁻¹]_{>N} = ζ^{−N−1}[q^{N+1}/(ζ − q)]_{≤N}. Als formale Potenzreihe in x ist
1/(ζ − a₀ − g) = Σ_{m≥0} g^m (ζ − a₀)^{−m−1} (jeder x-Koeffizient ist eine endliche Summe, da g = O(x)); wegen g^m = O(x^m)
tragen zu [·]_{≤N} nur m ≤ N bei. A⁻¹ ist eine feste Matrix. Die Formel setzt die für |ζ| > |a₀| bewiesene Identität
analytisch fort. ∎
**Folgerung.** Jede Koordinate von y_ζ ist eine rationale Funktion von ζ mit Polen nur in 0 und a₀ (je Ordnung ≤ N+1),
und y_ζ = O(|ζ|^{−N−2}) für ζ → ∞. Also ist ζ ↦ log R_N(ζ) = log Σ_j ρ^j |(y_ζ)_j| subharmonisch auf ℂ ∖ {0, a₀}.

### 12.2 Satz 12.2 (Bernstein–Walsh für die Pol-Antwort) — bewiesen

Sei K ⊂ ℂ kompakt, regulär für das Dirichlet-Problem, mit 0, a₀ ∈ int K, E := Ĉ ∖ K zusammenhängend, g_E die Greensche
Funktion von E mit Pol in ∞. Dann gilt für alle ζ ∈ E:
  R_N(ζ) ≤ S_N(K) · e^{−(N+2) g_E(ζ)},   S_N(K) := max_{ζ∈∂K} R_N(ζ).
*Beweis.* u(ζ) := log R_N(ζ) + (N+2) g_E(ζ) ist auf E ∖ {∞} subharmonisch (Lemma 12.1; g_E harmonisch). Bei ∞ ist
R_N(ζ)|ζ|^{N+2} beschränkt und g_E(ζ) − log|ζ| beschränkt, also u nach oben beschränkt und daher (hebbar) subharmonisch auf E.
Auf ∂E = ∂K ist g_E = 0 und u = log R_N stetig (die Pole 0, a₀ liegen im Inneren von K). Maximumprinzip: u ≤ log S_N(K). ∎

**Wahl K = W̄ = D̄(0,ρ) ∪ D̄(a₀,R₀).** g_E ist explizit: M(ζ) = (ζ − L)/(ζ − L̄) bildet beide Randkreise auf Strahlen durch 0
ab und E auf einen Sektor vom Öffnungswinkel πα; τ = (e^{−iθ}M)^{1/α} bildet ihn auf die obere Halbebene ab, und
g_E(ζ) = log|(τ(ζ) − τ̄₀)/(τ(ζ) − τ₀)|, τ₀ = τ(∞). Numerisch (c19): **Öffnungswinkel 158.9335°, α = 0.882964** (= Strang Bs
Aussenwinkel), Robin-Konstante −0.41669, und bei L: g_E(L(1+u)) / u^{1/α} = 0.9188, 0.9143, 0.8949, 0.8718, 0.8290 für
u = 0.001, 0.01, 0.05, 0.1, 0.2 (1/α = 1.1325). Im Vergleich: φ(ζ) = min(log|ζ|/ρ, log|ζ−a₀|/R₀) ≥ g_E (E ist kleiner als das Äussere jedes einzelnen Kreises), d. h.
die Max-Formel aus §9.4 sagt schnelleren Zerfall voraus als Bernstein–Walsh. **[Korrigiert in Runde 4, §15:]** Bernstein–Walsh ist die
Wahrheit (scharf), die Max-Formel ausserhalb W̄ nicht — e^{−g_E} trifft die gemessenen Raten auf 0.009, Poisson–Jensen mit Randdaten auf 0.0004.

**Was das bedeutet.** Die Skala N^{−α} = N^{−0.883} am Schlitz (Strang B, §9.4) ist jetzt *erklärt*: g_E ~ c u^{1/α} ⇒
e^{−N g_E} = e^{−c (N^α u)^{1/α}}. Und (H_R1), (H_R2) folgen beide aus der einen Grösse S_N := max_{∂W} R_N.

### 12.3 Satz 12.3 (Konsistenz unter einer Randbedingung) — bewiesen

**Voraussetzungen:** (F) wie in Satz 9.5 (von Strang A numerisch bestätigt: C_f ≤ 0.95 für u ≤ 0.2, max|f| ≤ 3.6 auf Γ');
**(S)** S_N = max_{ζ∈∂W} R_N(ζ) ≤ C_S N^b für ein b < α(σ+1) = 0.883·6.447 = 5.693.
**Aussage:** α_N ≤ (|L|/π) C_f C_S · α Γ(α(σ+1)) · c_g^{−α(σ+1)} · N^{b} (N+2)^{−α(σ+1)} + (M|Γ'|/2π) C_S N^b e^{−(N+2)γ'},
mit c_g := min_{0≤u≤u₀} g_E(L(1+u))/u^{1/α} (≈ 0.829 für u₀ = 0.2) und γ' := min_{Γ'} g_E > 0. Insbesondere α_N → 0 wie
N^{b − 5.693}.
*Beweis.* Satz 9.5 mit (H_R1) ersetzt durch Satz 12.2: R_N(L(1+u)) ≤ S_N e^{−(N+2)c_g u^{1/α}}, und
∫_0^∞ u^σ e^{−λu^{1/α}} du = α Γ(α(σ+1)) λ^{−α(σ+1)}; auf Γ' ist R_N ≤ S_N e^{−(N+2)γ'}. ∎

Damit ist **(H_R) bewiesen bis auf (S)** — eine Bedingung an die Pol-Antwort auf der festen Kurve ∂W, und zwar nur in der
schwachen Form „langsamer als N^{5.69}“. Numerisch liegt das Maximum auf ∂W genau in der Ecke ζ = L (bzw. L̄), und
S_N = R_N(L) wächst wie ≈ N^{1.1…1.2} (Tabelle 12.4).

### 12.4 Die Randgrösse S_N numerisch (c19, c20) — Maximum genau in der Ecke

| N | 40 | 80 | 120 | 160 | 200 | 250 | 300 |
|---|---|---|---|---|---|---|---|
| S_N = max_{∂W} R_N | 10.38 | 23.50 | 37.53 | 52.12 | 67.11 | 86.28 | 105.8 |
| Ort des Maximums | L̄/L | L̄/L | L̄/L | L | L | L | L |
| (N+1)/ρ = |[(L−x)⁻¹]_{≤N}| | 29.8 | 58.9 | 88.0 | 117.1 | 146.2 | 182.6 | 219.0 |
| Schranke 12.3 / α_N | 4.5 | 9.2 | 9.8 | 10.2 | 10.6 | — | — |

Weitere Randpunkte (N = 250 / 300): ρe^{i(arg L+0.01)}: 80.6 / 98.5; a₀+R₀e^{i(arg(L−a₀)−0.01)}: 46.4 / 50.5; −ρ: 28.1 / 44.6;
ρe^{2.4i}: 32.5 / 50.4; rechter Scheitel a₀+R₀: 0.0043 / 0.0014. Lokale Exponenten von S_N: 1.18, 1.15, 1.14, 1.13, 1.13, 1.12
⇒ S_N ≈ 0.33·N^{1.12}, also κ ≈ 1.1 ≪ 4.69 (Bs Schwelle) bzw. ≪ 5.69 (Satz 12.3). Bemerkenswert: S_N ist durchweg ein fester
Bruchteil (0.35 → 0.48) der Norm des **Kopfes** des Cauchy-Kerns in der Ecke, (N+1)/ρ — der „Faktor N“ (Strang B) ist die
Kopfnorm des Pols bei L. Probe von Satz 12.2: R_N/(S_N e^{−(N+2)g_E}) ≤ 0.98 an allen 9 Testpunkten (Schlitz, Gitter), N ≤ 200.

### 12.5 Abgleich mit Strang B (Satz B-Max) und (H*)

Strang B hat unabhängig dieselbe Struktur gefunden (b_reformulierung/BERICHT.md §7: Lemma 7.2 Maximumprinzip in ζ, Lemma 7.3
Green-Funktion). Gegenprobe: meine Robin-Konstante −0.41669 ⇔ Kapazität e^{0.41669} = 1.5169 = Bs Wert; α = 0.882964 = Bs
Exponent. **(H*) folgt aus (S):** Bs Π_N k_ζ ist (bis auf den Unterschied A ↔ A_f, §12.6) y_ζ + [k_ζ]_{≤N}; für ζ ∉ Ḡ' im
Randband ist |y_ζ| ≤ S_N (Satz 12.2, sogar ohne 1/dist) und sup_{D̄}|[k_ζ]_{≤N}| ≤ 1/(|ζ|−ρ) ≤ 1/dist(ζ,Ḡ'). Also
|Π_N k_ζ(z)| ≤ S_N + 1/dist(ζ,Ḡ') ≤ (S_N·diam + 1)/dist(ζ,Ḡ') — (H*) mit κ = b. **(H*) ⇔ (S) bis auf Konstanten:** (H*) ist
also genau die Aussage „Pol-Antwort auf ∂Ḡ' polynomial“, und das Maximum liegt in der Ecke.

### 12.6 A = J_N(q_N) gegen A_f = P_N J_∞(f) P_N (c21) — dieselben Antworten, obwohl astronomisch verschieden

A_f δ := [δ∘q + (f'∘f)δ]_{≤N} mit f'∘f = eˣ/f'(x) ist der exakte Abschnitt der wahren Linearisierung ([δ∘q] = [δ∘f], da q ≡ f).
A − A_f = −(Multiplikation mit [t'∘q]_{≤N}). Gemessen (l¹_ρ): |[q'∘q]_{≤N}| = 9.8, 2.6e6, 4.7e14, 3.7e22, 3.0e30 (N = 40 … 200),
dagegen |[f'∘f]_{≤N}| = 2.339 für alle N. **Trotzdem** stimmen α_N und R_N(L) für A und A_f auf alle gedruckten Stellen überein
(Differenz relativ ≈ |δ1|, d. h. A⁻¹M_{[t'∘q]}y ist für kohärente y von zweiter Ordnung). Deutung: Die Instabilität von J_N(q_N)
(Strang D, 20:41) sitzt im riesigen, aber kohärenten Multiplikator [q'∘q] — die a₀-Jets von q_N sind „Partialsummen-Müll“.
Für die lineare Theorie ist **A_f** das saubere Objekt: dort gilt exakt A_f w_{≤N} + B_N w_{>N} = [J_∞(f)w]_{≤N}, also
y^f(w) := A_f⁻¹B_N w_{>N} = A_f⁻¹[J_∞(f)w]_{≤N} − w_{≤N} ohne t'-Korrektur.

### 12.7 Satz 12.5 (Vermutung Q für auf W̄ holomorphe Funktionen) — bewiesen unter (S^f)

Sei S^f_N := max_{∂W}|A_f⁻¹B_N[(ζ−x)⁻¹]_{>N}| (numerisch = S_N). Ist v holomorph auf einer Umgebung U von W̄ und Γ_d = {g_E = d} ⊂ U
(d > 0), so gilt
  |A_f⁻¹[J_∞(f)v]_{≤N} − v_{≤N}| ≤ (|Γ_d|/2π)·max_{Γ_d}|v|·S^f_N·e^{−(N+2)d}.
*Beweis.* Exakte Identität 12.6: linke Seite = |A_f⁻¹B_N v_{>N}|. Cauchy: v_{>N} = (1/2πi)∮_{Γ_d} v(ζ)[(ζ−x)⁻¹]_{>N}dζ (|ζ| > ρ auf Γ_d),
Lemma 9.4, Lemma 12.1 und Satz 12.2 (mit A_f statt A; der Beweis braucht nur, dass A_f fest und invertierbar ist) auf Γ_d, wo
g_E = d. ∎ Für v, die auf W̄ holomorph sind, ist der Fehler also **exponentiell** klein (Rate e^{−d}, d bis zum Abstand der
nächsten Singularität von v in Green-Masseinheiten) — stärker als die Form „poly(N)·Taylor-Reste“ der ursprünglichen
Vermutung Q; für v mit Randsingularität in L (wie f_K) gibt Satz 12.3 die polynomiale Rate.

### 12.8 Was für einen vollständigen Beweis fehlt — und warum es nicht „weich“ geht

Alles ist auf **(S): S_N = max_{∂W} R_N = O(N^κ), κ < 4.69** zurückgeführt (gemessen κ ≈ 1.12, Maximum in der Ecke L).
Weiche Argumente reichen dafür nicht: y_ζ ist rational in ζ mit Polen nur in 0, a₀ — innerhalb von W wächst R_N wie r_W^{−N}
(Gesetz R), ausserhalb fällt es; jedes Maximumprinzip braucht eine Schranke auf einer Kurve um die Pole, und die einzige
„billige“ Information (|ζ| → ∞: Leitkoeffizient A⁻¹b_{N+1}) führt nur nach aussen. (S) ist daher genau das innere Eckproblem
(Strang B: Mellin-Symbol, Strang D: Toeplitz mit Sprungsymbol). Hinweise für dieses Problem aus meinen Rechnungen:
- J_∞(f) wirkt auf den Pol in der Ecke **nicht-resonant**: J_∞(f)(L−x)⁻¹ = c_L·(L−x)⁻¹ + h, c_L = (1+L)/√L (Bs Symbol μ^s + μ bei
  s = −1), h stetig auf D̄ (Singularität (x−L)^{γ−1}, Re = 3.45). Der Pol ist also für das innere Problem eine reguläre Eingabe.
- S_N/[(N+1)/ρ] = 0.35 … 0.48 (N = 40 … 300): die Antwort ist von der Grösse des Kopfes des Cauchy-Kerns; ein Beweis müsste
  zeigen, dass A_f⁻¹ auf den Kopf v_{≤N} = [(L−x)⁻¹]_{≤N} wie Multiplikation mit 1/c_L plus O(N^κ) wirkt.

### 12.9 Die Eckspitze von Strang D macht (S) für grosse N falsch — Poisson–Jensen-Verfeinerung (c23–c25)

**Befund (numerisch, c23).** Pol-Antwort längs ∂W, D-Bogen von L (t = 0) über −ρ (t = 1), A-Bogen = Rand von D(a₀,R₀):

| N | t=0 (L) | 0.02 | 0.05 | 0.1 | 0.2 | 0.35 | 0.5 | 0.7 | 0.9 | 1 (−ρ) | A-Scheitel a₀+R₀ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 67.1 | 55.7 | 47.2 | 39.9 | 32.4 | 26.4 | 22.7 | 19.7 | 18.1 | 17.9 | 1.3e-2 |
| 300 | 105.8 | 86.9 | 75.9 | 68.4 | 62.2 | 57.0 | 52.6 | 48.0 | 45.1 | 44.6 | 1.4e-3 |
| 400 | 145.8 | 120.0 | 109.8 | 108.6 | 116.8 | 126.8 | 129.6 | 127.2 | 123.2 | 122.2 | 1.6e-4 |

In der Ecke und bei t = 0.02 wächst R_N wie N^{1.12} (Faktor 1.38 von 300 auf 400), auf dem fernen D-Bogen schneller (t = 1:
Faktor 2.74), auf dem A-Bogen **fällt** es exponentiell. Die Antwort auf den Pol bei −ρ hat im Ausgang ihr Maximum bei
arg x = 76–77° mit 6.0, 18.6, 62.7 (N = 200, 300, 400) ≈ 1.0115^N — das ist Strang Ds Eckspitze. Zerlegt man R_N = P + S
(P ~ N^{1.12}, S ~ 1.0115^N), ist der Spitzenanteil bei N = 300: 0.3 %, 4 %, 12 %, 29 %, 63 %, 79 % (t = 0.02, 0.05, 0.1, 0.2,
0.5, 1), im Randparameter τ (Ecke ↦ 0) also ~|τ|^{1.2…2.5}. **Folge:** (S) bzw. (H*) mit polynomialem κ wird ab N ≈ 420 auf dem
fernen D-Bogen verletzt (dort überholt die Spitze die Ecke).

**Die Konsistenz selbst sieht davon nichts (c25):** α_N = 3.39e-12 (N = 300), 8.35e-13 (N = 400), lokaler Exponent −4.87; das
Maximum von |δ1| auf |x| = ρ liegt bei arg x = 76.6–77° (Randschicht an L, wird mit N schmaler: Wert bei 75° relativ 0.25 bzw.
0.10), keine exponentielle Komponente (N = 500 wegen der 20-min-Grenze nicht fertig). Das passt zu Strang A (Zweig bis 300 =
Lösung bei Π_N f_K) und B (Überlapp 8.7e-9).

**Satz 12.6 (Poisson–Jensen statt Maximum) — bewiesen.** Mit dem harmonischen Mass ω_ζ von E gilt für ζ ∈ E
  log R_N(ζ) ≤ ∫_{∂W} log R_N dω_ζ − (N+2) g_E(ζ).
(*Beweis:* u = log R_N + (N+2)g_E ist subharmonisch auf E, nach oben beschränkt, stetig bis ∂W; also u ≤ Poisson-Integral
seiner Randwerte.) Ist log R_N(η) ≤ κ log N + log C + N Λ(η) auf ∂W, folgt
  R_N(ζ) ≤ C N^κ exp(N[ω_ζ(Λ) − g_E(ζ)]).
**Explizite Konstanten (c24):** auf der Kontur ∂(W + 0.1) ist g_E/ω(ganzer D-Bogen) ≥ 0.0649 — eine Wachstumsrate λ < 0.065 auf
dem ganzen D-Bogen ist dort unschädlich (Spitze: λ ≈ 0.0114). Am Schlitz nahe L ist ω_ζ(D-Bogen ab t ≥ t₁)/g_E(ζ) ≤ 14.65, 2.41,
1.11 für t₁ = 0.02, 0.1, 0.2 (u → 0). Also: **Wenn die exponentielle Rate auf dem Teil des D-Bogens mit |τ| < 0.011 (|η − L| ≲ 0.05)
verschwindet und sonst λ < 1/14.65 = 0.068 ist, bleibt Satz 12.3 mit der polynomialen Eckgrösse gültig.** Bis N = 400 ist das
erfüllt (t = 0.02 rein polynomial). Für „alle N“ reicht das nicht, wenn der Spitzenanteil nahe L nur polynomial klein in |τ| ist
(dann wandert die Grenze t₁(N) ~ e^{−λN/p} gegen L und ω/g_E ~ 1/t₁ explodiert). Gebraucht wird dann eine Aussage über die
**Anregung der Spitze** durch Pole nahe der Ecke: sie muss exponentiell klein in N sein (WKB-artig, c_N(η) ~ e^{−N h(τ)} mit h > 0
nahe 0) — oder die Spitze wird für den Kneser-Eingang strukturell nicht angeregt (was die Messung zeigt).

### 12.10 Stand von (H_R)/(H*) nach Runde 3

- **Bewiesen:** Lemma 12.1 (Rationalität), Satz 12.2 (Bernstein–Walsh), Satz 12.3 ((F) + (S) ⇒ α_N ≤ C S_N N^{−5.69}), Satz 12.5
  (Q für holomorphe v unter (S)), Satz 12.6 (Poisson–Jensen), explizite Green-Funktion und harmonisches Mass (Möbius + Potenz).
- **Numerisch:** S_N = R_N(L) ≈ 0.33 N^{1.12} bis N = 400 (Maximum bisher in der Ecke); Schranke/α_N ≈ 10.
- **Nicht bewiesen, und in der Form „polynomial auf ganz ∂W“ für grosse N wohl falsch:** (S)/(H*). Richtige Form: polynomial in
  einer festen Eckumgebung (|τ| < 0.011) und Rate < 0.068 sonst (Satz 12.6) — numerisch bis N = 400; für alle N hängt es an der
  Spitzenanregung nahe L (Strang D/B: inneres Problem).

---

# Runde 4 (Auftrag 22:10): schwächste Eck-Bedingung, Gegenprobe bis N = 1000, Gesetz R = Green-Funktion

Alles für A = A_f = P_N J_∞(f_K) P_N (das saubere Objekt, §12.6; für J_N(q_N) numerisch identisch, c21). Zusätzliche Bezeichnungen:
- g = g_E, ω_ζ = harmonisches Mass von E = Ĉ∖W̄; τ: E → obere Halbebene (L ↦ 0, L̄ ↦ ∞, −ρ ↦ −1, a₀+R₀ ↦ +1, ∞ ↦ τ₀ =
  −0.05628+0.99842i); Randparameter s = τ(η) ∈ ℝ für η ∈ ∂W, D-Bogen ↦ s < 0, A-Bogen ↦ s > 0, Konjugation ↔ s ↦ |τ₀|²/s.
  κ₀ := 2π Im τ₀/|τ₀|² = 6.2732. Umrechnung auf dem D-Bogen (c28): |η−L| = 0.005 / 0.0107 / 0.025 / 0.05 / 0.1 ↔ |s| = 0.00081 /
  0.0019 / 0.005 / 0.011 / 0.024 (Winkel hinter L: 0.21° / 0.45° / 1.04° / 2.08° / 4.17°).
- Γ_d := {g = d} (Niveaulinie; für d = 0.1 radial 0.15…0.20 ausserhalb ∂W, also im von Strang A für (F) geprüften Bereich bis
  ∂(W+0.2)); S_L(d) := {L(1+u): 0 < u ≤ u₀(d)} mit g(L(1+u₀)) = d (u₀(0.05) = 0.079, u₀(0.1) = 0.151 < 0.2 = Prüfbereich von C_f).
  Der Schlitz liegt im τ-Bild unter arg τ = 78.07° (u → 0; 79.1° bei u = 0.15) — er neigt zum A-Bogen (Winkelhalbierende 90°).
- y_ζ := A_f⁻¹B_N[k_ζ]_{>N}, k_ζ = (ζ−x)⁻¹, R_N(ζ) = |y_ζ| (ℓ¹_ρ); w_ζ := [J_∞(f)k_ζ]_{≤N}, Π_N k_ζ = A_f⁻¹w_ζ.
- **Eckzerlegung (Strang D):** Singulärwertzerlegung in ℓ²_ρ, A_f⁻¹ = Σ_i σ_i u_i⟨v_i,·⟩; r = Zahl der Eckmoden (D: r ≤ 3 bis N = 800,
  die übrigen σ_i ≤ 1.3663). Q_N := Σ_{i≤r} u_i⟨u_i,·⟩ (Orthogonalprojektor auf die Eck-Ausgänge), 𝓟_N := Q_N A_f⁻¹ = Σ_{i≤r}σ_i u_i⟨v_i,·⟩,
  𝓡_N := A_f⁻¹ − 𝓟_N (‖𝓡_N‖_{ℓ²_ρ} = σ_{r+1}). **Eckanteil** K_N(ζ) := |Q_N y_ζ|, **Restanteil** R°_N(ζ) := |(I − Q_N)y_ζ|.

## 13. Die Eck-Bedingung (E) und Satz „(F) + (E) ⇒ α_N → 0 für alle N“

### 13.1 Lemma 13.0 (die Zerlegung erhält Rationalität und Abfall) — bewiesen

(i) Q_N y_ζ = 𝓟_N B_N[k_ζ]_{>N} und (I − Q_N)y_ζ = 𝓡_N w_ζ − (I − Q_N)[k_ζ]_{≤N}.
(ii) ⟨u_i, y_ζ⟩ = ⟨u_i, Π_N k_ζ⟩ − ⟨u_i, [k_ζ]_{≤N}⟩ (Strang Ds Eckfunktional, Lemma D6, minus Kopf-Überlapp).
(iii) Q_N y_ζ und (I − Q_N)y_ζ sind koordinatenweise rational in ζ mit Polen nur in 0 und a₀ und O(ζ^{−N−2}) für ζ → ∞. Also sind
log K_N und log R°_N subharmonisch auf E∖{∞}, und log K_N + (N+2)g, log R°_N + (N+2)g sind bei ∞ nach oben beschränkt.
*Beweis.* (i) Nach Definition der Singulärwertzerlegung ist ⟨u_i, A_f⁻¹x⟩ = σ_i⟨v_i, x⟩, also Q_N A_f⁻¹ = 𝓟_N und (I−Q_N)A_f⁻¹ = 𝓡_N.
Die exakte Identität A_f[k]_{≤N} + B_N[k]_{>N} = [J_∞(f)k]_{≤N} (§12.6) gibt y_ζ = A_f⁻¹w_ζ − [k_ζ]_{≤N}; anwenden von I − Q_N.
(ii) folgt aus y_ζ = Π_N k_ζ − [k_ζ]_{≤N}. (iii) Q_N ist eine feste Matrix; Lemma 12.1 für y_ζ; der Logarithmus der ℓ¹-Norm eines
holomorphen Vektors ist subharmonisch. ∎
**Wichtig für die Messung:** Ds Funktional ⟨u_i, Π_N k_ζ⟩ selbst fällt bei ∞ nur wie ζ^{−1}; der Kopf-Überlapp ⟨u_i, [k_ζ]_{≤N}⟩ ist
auf ∂W durch √(N+1)/ρ beschränkt, zerstört aber den Faktor e^{−(N+2)g}. Für Maximumprinzip und Poisson–Jensen ist ⟨u_i, y_ζ⟩ das
richtige Eckfunktional (am Schlitz sind beide fast gleich, sobald u_i in den oberen Koeffizientenbändern sitzt:
|⟨u_i,[k_ζ]⟩| ≤ Σ_j |u_{i,j}|ρ^j|ζ|^{−j−1} ≈ e^{−0.8Nu} bei Träger j/N ≥ 0.8).

### 13.2 Die Eck-Bedingung (E)

- **(E_R) [Rest — Strang D, Teil (1)]:** max_{η∈∂W} R°_N(η) ≤ C_R N^{κ_R}.
- **(E_K) [Ecke, Schlitzform]:** Es gibt d > 0 und ε ∈ (0,1] mit K_N(ζ) ≤ C_K N^{κ_K} e^{−ε(N+2)g(ζ)} für alle ζ ∈ S_L(d) ∪ Γ_d.
- **(E) := (E_R) ∧ (E_K)** mit κ := max(κ_R, κ_K) < α(σ+1) = 0.882964·6.447 = 5.693.

**Gebiet und Rate.** (E_K) betrifft nur Pole **ausserhalb** W̄ nahe der Ecke (Schlitz) und auf einer festen Niveaulinie, nicht den
Rand ∂W selbst. Wegen des Gewichts u^σ (σ = 5.447) und des Faktors e^{−εNg} zählt am Schlitz nur u ≲ (C log N/(εN))^α, also Pole im
Abstand ~N^{−0.883} von L; dort ist Ng = O(1) und (E_K) verlangt nur „polynomial“. Weiter aussen verlangt (E_K) Zerfall mit einem
beliebig kleinen Bruchteil ε der Green-Rate. ε geht nur als Faktor ε^{−5.693} in die Konstante ein; ε = 0 (bloss polynomial ohne
Zerfall) reicht nicht (dann ist das Schlitzintegral O(1)).

### 13.3 Satz 13.1 ((F) + (E) ⇒ α_N → 0 für alle N) — bewiesen

**Voraussetzungen:** (F) wie in Satz 9.5 mit Γ' = Γ_d und u₀ = u₀(d) (Konstanten C_f, M := max_{Γ_d}|f_K|); (E).
**Aussage:** Für alle N ≥ 1 ist
  α_N ≤ (|L|/π)·C_f·αΓ(α(σ+1))·c_g^{−α(σ+1)}·[C_R N^{κ_R}(N+2)^{−α(σ+1)} + C_K N^{κ_K}(ε(N+2))^{−α(σ+1)}]
      + (M|Γ_d|/2π)·[C_R N^{κ_R}e^{−(N+2)d} + C_K N^{κ_K}e^{−ε(N+2)d}],
c_g := min_{0<u≤u₀} g(L(1+u))/u^{1/α} (≥ 0.829 für u₀ ≤ 0.2, c19). Insbesondere α_N = O(N^{κ−5.693}) → 0.
*Beweis.* (1) Lemma 9.2: α_N = |A_f⁻¹B_N t_N|. (2) Wie im Beweis von Satz 9.5: Ω := {g < d} ∖ (S_L ∪ S_L̄); nach (F) ist f_K holomorph
auf Ω und stetig bis an beide Ufer, also f_k = (1/2πi)∮_{∂Ω} μ(ζ)ζ^{−k−1}dζ (k > N) mit μ = Δf auf den Schlitzen, μ = f_K auf Γ_d.
Lemma 9.4 (gilt für jede feste Matrix, also auch für A_f) gibt
  α_N ≤ (1/2π)[2|L|∫_0^{u₀}|Δf(L(1+u))|·R_N(L(1+u)) du + ∫_{Γ_d}|f_K|·R_N|dζ|]   (R_N(ζ̄) = R_N(ζ), da A_f reell).
(3) R_N ≤ R°_N + K_N (Dreiecksungleichung). (4) Restanteil: (I − Q_N)y_ζ erfüllt die Voraussetzungen von Satz 12.2 (Lemma 13.0 (iii)),
also R°_N(ζ) ≤ C_R N^{κ_R}e^{−(N+2)g(ζ)} auf ganz E. (5) Eckanteil: (E_K) auf S_L(d) ∪ Γ_d. (6) |Δf| ≤ C_f u^σ, g(L(1+u)) ≥ c_g u^{1/α}
und ∫_0^∞ u^σ e^{−λc_g u^{1/α}}du = αΓ(α(σ+1))(λc_g)^{−α(σ+1)} mit λ = N+2 bzw. ε(N+2); auf Γ_d ist g = d. ∎
(Mit den gemessenen κ ≈ 1.1 gäbe das N^{−4.6}; gemessen ist N^{−4.85}, §14.1.)

### 13.4 Randform von (E_K): Rate gegen Gebiet — Lemma 13.2 und Satz 13.3 bewiesen, Konstanten numerisch (c28)

**Lemma 13.2 (Poisson–Jensen, vorzeichenbehaftet).** Sei h ein Vektor wie in Lemma 13.0 (iii) (z. B. h = Q_N y_ζ), H := |h|, und auf ∂W
gelte log H(η) ≤ log(C N^κ) + N λ_N(η) mit einer messbaren, nach oben beschränkten Funktion λ_N: ∂W → [−∞, ∞) — **das Vorzeichen ist
frei**. Gilt
  **(E_K^∂)**  ω_ζ(λ_N) ≤ (1 − ε) g(ζ) für alle ζ ∈ S_L(d) ∪ S_L̄(d) ∪ Γ_d,
so ist H(ζ) ≤ C N^κ e^{−(εN+2)g(ζ)} dort, insbesondere (E_K).
*Beweis.* Satz 12.6 für log H: log H(ζ) ≤ ω_ζ(log H) − (N+2)g(ζ) ≤ log(CN^κ) + Nω_ζ(λ_N) − (N+2)g(ζ). ∎

**Satz 13.3 (Schärfe).** Ist der Eckanteil vom Rang 1, Q_N y_ζ = c(ζ)·u, und hat c in E keine Nullstelle (auch nicht in ∞, d. h.
c(ζ)ζ^{N+2} ↛ 0), so gilt Gleichheit: log|c(ζ)| = ω_ζ(log|c|) − (N+2)g(ζ). Dann ist (E_K^∂) mit λ_N := (log|c| − log CN^κ)/N auch
**notwendig** für (E_K); Nullstellen können nur helfen.
*Beweis.* log|c| + (N+2)g ist harmonisch auf E einschliesslich ∞ und stetig bis ∂W, wo g = 0; das Dirichlet-Problem ist eindeutig. ∎
Die Randform verliert also nichts: die Randdaten bestimmen die Werte am Schlitz. (Für die volle Norm R_N bestätigt der Gittertest §15,
dass Poisson–Jensen bis auf einen Faktor e^{−(0.03…0.55)} eine Gleichung ist.)

**Explizite Konstanten (c28, c24).**
(a) *Grenzwert in der Ecke.* Ist λ_N(s) = O(|s|^{1+δ}) bei s = 0, so gilt lim_{u→0} ω_{L(1+u)}(λ_N)/g(L(1+u)) = (1/κ₀)∫_ℝ λ_N(s)s^{−2}ds
(Probe c28: Formel 1.502 / 0.3004 gegen exakt 1.500 / 0.3003 für zwei glatte Profile). Für ein in s **lineares, vorzeichenwechselndes**
Profil λ = −a·s nahe 0 (Wachstum auf der D-Seite, Zerfall auf der A-Seite — so sieht der Realteil einer in τ holomorphen Phase aus) ist
der lokale Beitrag ω(λ) ≈ −a Re τ, im Verhältnis −aπ cot(78.07°)/κ₀ = −0.106a < 0: harmlos, gleich wie gross a ist. Ein **einseitig**
lineares Profil (λ = a|s| nur auf der D-Seite) gibt (a/κ₀)·ln(1/|τ|) → ∞ (logarithmisch), ein in der Bogenlänge t lineares (λ = a·t ∝
|s|^{0.883} nur auf der D-Seite — so liest Strang D sein Kubik-Modell, BRETT 22:20) ~|τ|^{−0.117} ~ N^{0.12} am massgeblichen Ort
|τ| ~ 1/N: bei N = 1000 klein (Ds Wert 0.11), formal aber nicht „für alle N“. **Für alle N entscheidet das Vorzeichen auf dem A-Bogen
nahe L** (oder ein Einsatz schneller als linear in s). **Nachtrag (Strang D, BRETT 22:40):** genauer gemessen setzt Ds Eckrate wie
t^{1.44…1.52} ein (Kubik-Modell bis N = 1200), also wie |s|^{1.27…1.34} — schneller als linear in s. Dann ist der Grenzwert endlich, Ds
Mass ν(dt) ≈ 0.199t^{−2.13}dt ist genau (1/κ₀)ds/s² (s ≈ 0.91t^{1.1325}), und ∫λdν = 0.086 gleichmässig in t → 0 heisst: (E_K^∂) gilt
dort schon einseitig (ohne A-Bogen) mit 1 − ε ≈ 0.09.
(b) *Stufenform.* λ_N ≤ 0 auf den Eckbögen |s| < s₁ (bei L) und |s| > |τ₀|²/s₁ (bei L̄) und auf dem A-Bogen, λ_N ≤ λ₀ auf dem übrigen
D-Bogen ⇒ (E_K^∂) mit 1 − ε = λ₀/λ*(s₁), λ*(s₁) := 1/max(sup_Schlitz ω/g, sup_{Γ_d} ω/d) (für das Indikatorprofil):

| spitzenfreier Eckbogen |η − L| ≤ | 0.006 | 0.0106 | 0.025 | 0.050 | 0.12 | 0.25 |
|---|---|---|---|---|---|---|
| s₁ | 0.001 | 0.0019 | 0.005 | 0.0109 | 0.03 | 0.066 |
| λ*(s₁), d = 0.1 | 0.0063 | 0.0119 | 0.0314 | 0.0684 | 0.103 | 0.104 |
| λ*(s₁), d = 0.05 | 0.0063 | 0.0119 | 0.0314 | 0.0508 | 0.0508 | 0.0509 |

Für kleine s₁ ist λ* = κ₀s₁ (Schlitz-Grenzwert, exakt bestätigt), für grosse begrenzt die Niveaulinie: λ* ≈ d/sup_{Γ_d}ω(D-Bogen) ≈ 1.03d.
Mit der Eckrate 0.0117…0.0120 (D; c27: 0.0118 bei −ρ) genügt ein spitzenfreier Eckbogen |η−L| ≤ 0.011 knapp, |η−L| ≤ 0.025 mit Reserve 2.6.
(c) *Was nicht reicht* (Folgerung aus Satz 13.3): ein Eckanteil, der auf ∂W beiderseits nahe L ≍ e^{cN}|s|^p·poly(N) ist (feste Rate,
nur polynomiale Unterdrückung zur Ecke hin) und keine Nullstellen hat. Dann ist am Schlitz bei u ~ N^{−α} (|τ| ~ 1/N) K_N ≥ e^{cN}N^{−p−O(1)}:
exponentiell, (E_K) falsch, und die Betragsabschätzung von α_N bricht für grosse N. (E) verlangt also, dass die **Rate selbst** zur Ecke
hin verschwindet (WKB-artig, K_N ≈ e^{Nλ(η)} mit λ → 0 bei L im Sinn von (a)), nicht nur die Amplitude.

### 13.5 Einteilige Form ohne Zerlegung; kohärente Minimalform

**Satz 13.4 — bewiesen.** Sei λ_N(η) := (1/N)log(R_N(η)/R_N(L)) (vorzeichenbehaftet, = 0 in der Ecke). Gelten (F), R_N(L) ≤ C N^κ mit
κ < 5.693 und Θ_N := sup_{ζ∈S_L(d)∪Γ_d} ω_ζ(λ_N)/g(ζ) ≤ 1 − ε, so ist α_N = O(N^{κ−5.693}).
*Beweis.* Lemma 13.2 mit h = y_ζ (Lemma 12.1 statt 13.0) gibt R_N ≤ R_N(L)e^{−(εN+2)g} auf S_L ∪ Γ_d; weiter wie Satz 13.1. ∎
Diese Form braucht Ds Zerlegung nicht und ist direkt aus Randdaten messbar (Θ_N, §14.3).

**Satz 13.5 (kohärente Minimalform) — bewiesen (trivial).** (F) + (E_R) + (E_koh) ⇒ α_N → 0, mit (E_koh): |Q_N δ1| = |𝓟_N B_N t_N| → 0
(Eckanteil der Konsistenz selbst). Kette: (E_K^∂) ⇒ (E_K) ⇒ [∫u^σK_N(L(1+u))du + ∫_{Γ_d}K_N|dζ| → 0] ⇒ (E_koh). (E_koh) ist die
schwächste Form, aber keine Aussage über die Pol-Anregung k_η mehr. Numerisch |Q_Nδ1|_{ℓ¹} ≤ √(N+1)·α_N ≤ 3.2e-13 bei N = 1000 (§14.1).

## 14. Gegenprobe jenseits N ≈ 420 (c26, c27, c27b, c29) — numerisch

### 14.1 Die Konsistenz selbst bis N = 1000 (c26, aus Strang As Π_N f_K, ohne eigene Matrix)

Exakt (§12.6): δ1^f = A_f⁻¹B_N t_N = Π_N f_K − q_N, also α_N = |Π_N f_K − q_N|_{ℓ¹_ρ} aus `gemeinsam/lin_PiN_fK_bis_1000.json` (60 Stellen) und
denselben a_k, die A benutzt hat (`a_moden/daten/kneser_koeff_lang.json`). Probe: N = 300 / 400 gibt 3.389e-12 / 8.350e-13 = meine c25-Werte mit
A = J_N(q_N).

| N | 100 | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|
| α_N (s = 0) | 6.80e-10 | 2.39e-11 | 3.39e-12 | 8.35e-13 | 2.86e-13 | 1.18e-13 | 5.57e-14 | 2.94e-14 | 1.66e-14 | 9.95e-15 |
| Max von |δ1| auf |x| = ρ bei | 77.09° | 76.82° | 76.75° | 76.71° | 76.69° | 76.67° | 76.65° | 76.65° | 76.64° | 76.64° |

Lokale Exponenten −4.72 … −4.98 ohne Trend; Fit N ≥ 300: log α_N = 1.25 − 4.851 log N + (2.1 ± 1.7)e-5·N (N ≥ 500: (3.6 ± 4.5)e-5·N) — mit 0
verträglich und 500-mal kleiner als die Eckrate 0.0118 (ein Spitzenanteil wäre von N = 300 bis 1000 um e^{8.3} ≈ 4000 relativ gewachsen).
Das Maximum von |δ1| auf dem Kreis wandert gegen arg L = 76.62° (Randschicht an L wird schmaler; Werte 1° daneben: 0.04 bzw. 0.15·Max bei N = 1000),
|δ1(−ρ)| ≈ 1e-18. **δ1 spürt die Spitze bis N = 1000 nicht.** (Gewicht (1+j)⁴: α_N = 2.6e-2 → 1.75e-3, Exponent −1.19, ebenfalls ohne
Exponentialanteil.)

### 14.2 Pol-Antwort längs ∂W und am Schlitz bis N = 1000 (c27, A_f, Präzision 1.4N+300 Bit; N = 400 und 800 mit +440/+480 Bit identisch)

D-Bogen: Anteil t des Wegs von L (t = 0) nach −ρ (t = 1), |η−L| = 0.0124 / 0.025 / 0.050 für t = 0.005 / 0.01 / 0.02. A-Bogen: Anteil t von L zum
Scheitel a₀+R₀. Schlitz S_v := L(1 + v/N^{0.883}).

| N | L | t=0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.5 | −ρ | A 0.02 | A-Scheitel | S_0.3 | S_1 | S_4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 67.1 | 63.4 | 60.4 | 55.7 | 47.2 | 39.9 | 22.7 | 17.9 | 14.2 | 1.3e-2 | 42.6 | 14.9 | 0.188 |
| 400 | 145.8 | 135.7 | 128.8 | 120.0 | 109.8 | 108.6 | 129.6 | 122.2 | 8.95 | 1.6e-4 | 90.8 | 30.5 | 0.331 |
| 500 | 186.7 | 173.2 | 164.7 | 155.5 | 151.4 | 169.4 | 344 | 365 | 6.89 | 1.8e-5 | 115.4 | 38.2 | 0.398 |
| 600 | 228.4 | 211.4 | 201.9 | 193.8 | 204 | 266.6 | 961 | 1140 | 5.37 | 2.0e-6 | 140.4 | 45.9 | 0.462 |
| 700 | 270.6 | 250.3 | 240.5 | 235.6 | 272.6 | 427.6 | 2761 | 3648 | 4.23 | 2.3e-7 | 165.5 | 53.6 | 0.524 |
| 800 | 313.3 | 289.9 | 280.5 | 281.5 | 364.2 | 701.3 | 8052 | 1.18e4 | 3.35 | 2.6e-8 | 190.7 | 61.2 | 0.585 |
| 900 | 356.4 | 330.3 | 322.2 | 332.4 | 489.2 | 1174 | 2.36e4 | 3.83e4 | 2.66 | 3.0e-9 | 216.1 | 68.8 | 0.644 |
| 1000 | 399.9 | 371.4 | 365.8 | 389.3 | 662.7 | 2000 | 6.96e4 | 1.25e5 | 2.12 | 3.5e-10 | 241.6 | 76.3 | 0.702 |

- **Ecke und Schlitz polynomial bis N = 1000:** R_N(L) ≈ N^{1.10} (400 → 1000); die Schlitzwerte stehen in festem Verhältnis zur Ecke
  (S_0.3/L = 0.623 → 0.604, S_1/L = 0.209 → 0.191, S_4/L = 0.00227 → 0.00176 für N = 400 → 1000, leicht fallend) — **die für α_N
  massgeblichen Pole (Abstand ~N^{−0.883} ausserhalb der Ecke) sehen keine Spitze**; R_N(L(1+N^{−0.883}))/N^{0.883} = 0.148 (300) → 0.171 (1000).
- Die Spitze überholt die Ecke bei N ≈ 450 (−ρ: 365 > 187 bei N = 500); Rate bei −ρ: 0.0118 je Grad (800 → 1000). Ihre Front wandert
  langsam zur Ecke: sichtbarer Überschuss (≥ 2 % über dem kleinsten Eckverhältnis) bei t = 0.05 (|η−L| = 0.12) ab N ≈ 350, t = 0.02 (0.05)
  ab ≈ 550, t = 0.01 (0.025) ab ≈ 900; t = 0.005 (0.012) bis N = 1000 nur +0.4 % (Verhältnis zu R_N(L): 0.931, 0.928, 0.926, 0.925, 0.925,
  0.927, 0.929 für N = 400 … 1000).
- Der A-Bogen fällt überall exponentiell (Scheitel: Rate −0.021 je Grad; t = 0.02: −0.0024).

### 14.3 Zerlegung in Eck- und Spitzenteil (c27b) — die Spitzenrate geht zur Ecke hin zurück

Fit je Randpunkt R_N(t) = a(t)·R_N(L) + b(t)e^{r(t)N} über N = 400 … 1000 (log-Residuen ≤ 0.009):

| t (D-Bogen) | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | 0.35 | 0.5 | 0.7 | 1 |
|---|---|---|---|---|---|---|---|---|---|---|
| Spitzenrate r(t) | 0.0016 | 0.0023 | 0.0032 | 0.0035 | 0.0057 | 0.0080 | 0.0098 | 0.0108 | 0.0115 | 0.0118 |
| Spitzenanteil bei N = 1000 | 4.5 % | 11 % | 26 % | 81 % | 94 % | 99 % | > 99 % | > 99 % | > 99 % | > 99 % |

(Zum Vergleich: der Eckteil wächst wie N^{1.1}, d. h. mit „Rate“ ≈ 0.0012 in diesem N-Bereich.) Am Schlitz und auf dem A-Bogen findet der
Fit keinen Spitzenteil. **Deutung (numerisch):** Die Anregung der Spitze durch Pole auf dem D-Bogen hat ein Ratenprofil, das zur Ecke hin
abfällt — WKB-artig, nicht „feste Rate mal |s|^p“ (sonst wäre r(t) ≈ 0.0118 überall). Überschuss über die Eckrate: 0.0004 / 0.0011 / 0.0020 /
0.0045 bei t = 0.005 / 0.01 / 0.02 / 0.1 (lokale Exponenten in t zwischen 0.9 und 1.5 — bei so kleinen Spitzenanteilen unsicher). Strang D
misst das Eckfunktional im Kubik-Modell direkt: Einsatz ∝ t^{1.44…1.52} bis N = 1200 (BRETT 22:40) — verträglich, und nach §13.4 (a)
ausreichend. Auf dem A-Bogen ist die Antwort fallend (zusätzliche Reserve).

### 14.4 Die einteilige Bedingung direkt gemessen: Θ_N (c29)

λ_N(η) = (1/N)log(R_N(η)/R_N(L)) auf 135 Randpunkten (beide Bögen, Ecke geometrisch verdichtet bis |η−L| ~ 1e-9, Spiegelung für die untere
Hälfte), Poisson-Integral exakt für stückweise lineare Randdaten. Θ_Schlitz := sup_{0<u≤0.15} ω(λ_N)/g am Schlitz, Θ_Γ := sup_{g=d} ω(λ_N)/d.
Satz 13.4 braucht Θ ≤ 1 − ε.

| N | 100 | 200 | 300 | 400 | 600 | 800 | 1000 |
|---|---|---|---|---|---|---|---|
| Θ_Schlitz (sup über 0 < u ≤ 0.15, erreicht bei u = 0.15) | −0.168 | −0.126 | −0.106 | −0.094 | −0.079 | −0.069 | −0.063 |
| ω/g am Schlitz, u = 1e-6 | −4.70 | −4.43 | −4.29 | −4.20 | −4.08 | −3.99 | −3.92 |
| ω/g am Schlitz, u = 0.01 | −0.627 | −0.539 | −0.490 | −0.455 | −0.406 | −0.372 | −0.347 |
| Θ_Γ, d = 0.05 | −0.165 | −0.099 | −0.057 | −0.016 | +0.045 | +0.082 | +0.105 |
| Θ_Γ, d = 0.1 | −0.110 | −0.062 | −0.035 | −0.012 | +0.018 | +0.036 | +0.048 |

**Befund:** Θ_N ≤ 0.105 (d = 0.05) bzw. ≤ 0.048 (d = 0.1) für alle gemessenen N — die Bedingung von Satz 13.4 hält bis N = 1000 mit
ε ≥ 0.89, also mit grosser Reserve. Am Schlitz ist das Verhältnis **negativ** (die Randschicht in der Ecke — R_N fällt von L aus nach beiden
Seiten — und der fallende A-Bogen überwiegen); positiv wird es nur auf der Niveaulinie durch die ferne Spitze, und dort ist es durch deren
Sättigungsrate begrenzt: 0.0118/λ* ≈ 0.23 (d = 0.05) bzw. 0.11 (d = 0.1), Tabelle §13.4 (b) — für alle N unschädlich. Die Schlitzwerte
steigen langsam (u = 1e-6: +0.27, +0.23, +0.21 je Verdopplung von N, abnehmend); mit Ds Einsatz t^{1.5} ist der Spitzenbeitrag zum
Eck-Grenzwert gleichmässig beschränkt (≤ 0.09), während der Beitrag der Randschicht skaleninvariant ist (Profil in N^α|η−L|) — daher
**Vermutung:** Θ_N bleibt für alle N unter 0.25.

### 14.5 Eckprofil von y_L (c22-Frage von Strang B; c22 selbst nicht gelaufen, ersetzt durch c27 mit python-flint, gleiche Grösse)

y_L/v_lo mit v_lo = [(L−x)⁻¹]_{≤N} in Bändern j/N (Mittel der komplexen Quotienten, N = 400; N = 200 … 1000 gleich bis auf langsame Drift):
[0, 0.5): −0.11+0.05i; [0.5, 0.7): −0.67; [0.7, 0.8): −0.99; [0.8, 0.853): −1.006; [0.853, 1]: **−1.000** (auf 3 Stellen). Also
Π_N k_L = v_lo + y_L ≈ (1 − β(j/N))·L^{−j−1}x^j mit einem festen **Fenster** β: ≈ 0.1 unten, Übergang bei j/N ≈ 0.5…0.75, **= 1 im ganzen
Streckband [1/|√L|, 1]** — die abgeschnittene Projektion löscht die oberen ~25 % der Koeffizienten des Eck-Pols, wie bei Π_N f_K (b_N ≈ 1e-40·a_N
bei N = 1000, c26). Daraus R_N(L) = |y_L| ≈ (N/ρ)∫β (der „Faktor N“, §12.4). Antwort auf Bs Frage: ja, im Band [0.853, 1] ist y_L = β·v_lo
mit festem β = −1. Das Fenster wandert langsam nach unten (Band [0.5, 0.7): |β| = 0.56, 0.68, 0.74, 0.79, 0.82 für N = 200, 400, 600, 800,
1000). (Numerisch.)

## 15. Gesetz R = Green-Funktion bzw. Poisson–Jensen (Abgleich mit Strang E) — Korrektur von §0b, §9.4, §12.2

**Aussage.** Für ζ ∉ W̄ gilt
  log R_N(ζ) = ω_ζ(log R_N|_{∂W}) − (N+2)·g(ζ) − Δ_N(ζ),   Δ_N(ζ) ≥ 0.
Δ_N ≥ 0 ist Satz 12.6 (**bewiesen**); dass Δ_N beschränkt ist (0.03 … 0.55 an 32 Gitterpunkten, von N = 80 auf 160 nur um ≤ 0.11 wachsend),
ist **numerisch** (c29). Bei polynomialen Randdaten (bis N ≈ 400: max_{∂W}R_N = R_N(L) ~ N^{1.1}, auf dem A-Bogen fallend) heisst das: Rate
**e^{−g(ζ)} je Grad** (Bernstein–Walsh, scharf), **nicht** max(ρ/|ζ|, R₀/|ζ−a₀|) = e^{−φ(ζ)}. Wegen φ ≥ g sagt die Max-Formel zu schnellen
Zerfall voraus, und der Fehler wächst exponentiell in N (Strang E im Möbius-Modell: Faktor 411 bei N = 160; BRETT 22:01).

**Gittertest (c29; dieselben 32 Aussenpunkte wie c17, beobachtete Rate (R_160/R_80)^{1/80}):**

| Vorhersage | log(beob./Vorhersage): Mittel | Streuung | min | max |
|---|---|---|---|---|
| Max-Formel e^{−φ} (alt, §9.4) | +0.0324 | 0.0214 | −0.018 | +0.070 |
| Green-Funktion e^{−g} | −0.0035 | 0.0091 | −0.021 | +0.009 |
| Poisson–Jensen mit gemessenen Randdaten | −0.0008 | 0.0004 | −0.0014 | −0.0001 |

Beispiele (e^{−g} / Max / beob. / PJ): ζ = −2.2: 0.6403 / 0.6248 / 0.6442 / 0.6451; ζ = 0.2+2.4i: 0.6130 / 0.5708 / 0.6123 / 0.6127; ζ = 1.9:
0.9652 / 0.9628 / 0.9454 / 0.9454 — dort zerfällt die Antwort schneller als e^{−g}, weil die Randdaten auf dem nahen A-Bogen selbst fallen
(Rate ≈ −0.022 im Scheitel a₀+R₀, §14.2); Poisson–Jensen erfasst das exakt. Der „systematische Vorlauf +0.016“ von §9.4 war genau
log(e^{−g}/Max) (Strang E: ζ = 2.2: 0.805/0.793).

**Im Inneren von W** (Wachstum) bleibt die Max-Formel eine Beobachtung (ζ = −0.4: 3.436/3.433; Abweichung in der Linse beider Kreise,
ζ = 0.8: 4.48/3.96) — **Vermutung**; naheliegende Erklärung über die Partialbruchzerlegung y_ζ = Σ_k α_kζ^{−k} + Σ_m β_m(ζ−a₀)^{−m} (je ein Teil
pro Pol), nicht geprüft. Die Beweise (Satz 12.2, 12.6, 13.1–13.4) brauchen nur die Aussenform.

**Konsequenz:** Gesetz R ist kein eigenes Postulat mehr, sondern Satz 12.6 plus die Randdaten log R_N|_{∂W}; die einzige offene Grösse sind
diese Randdaten — polynomial in der Ecke, exponentiell (Spitze, Rate 0.0118) auf dem fernen D-Bogen, fallend auf dem A-Bogen — und genau
darüber spricht (E).

### 15.1 Stand nach Runde 4

- **Bewiesen:** Lemma 13.0, Satz 13.1 ((F) + (E_R) + (E_K) ⇒ α_N = O(N^{κ−5.693}) für alle N), Lemma 13.2 (vorzeichenbehaftete Randform),
  Satz 13.3 (Schärfe ohne Nullstellen), Satz 13.4 (einteilige Form), Satz 13.5 (kohärente Minimalform); Gesetz R als Ungleichung (Satz 12.6).
- **Numerisch:** α_N ~ N^{−4.85} bis N = 1000 ohne Exponentialanteil (c = (2 ± 2)e-5); Ecke und Schlitz polynomial bis N = 1000, die Spitze
  (0.0118/Grad) beherrscht ab N ≈ 450 den fernen D-Bogen, ihr Ratenprofil fällt zur Ecke hin ab (WKB-artig); die einteilige Bedingung von Satz 13.4 hält bis N = 1000 mit Θ_N ≤ 0.105 (ε ≥ 0.89);
  Poisson–Jensen ist bis auf O(1) eine Gleichung.
- **Offen:** (E_R) und (E_K) für alle N (Strang D: Rest über Satz D8, Ecke über das Einsatzprofil); entscheidend ist das Ratenprofil
  der Eckanregung nahe L **mit Vorzeichen**: einseitig linear in t wäre nur bis N ~ 10^{10} unschädlich; Ds gemessener Einsatz t^{1.5}
  (Kubik-Modell bis N = 1200) und der fallende A-Bogen (Kneser, bis N = 1000) reichen — bewiesen ist beides nicht.
