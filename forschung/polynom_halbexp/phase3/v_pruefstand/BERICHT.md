# Strang V (Prüfstand) — Bericht, Phase 3

Auftrag (Leitung, 11.10.2026 19:50): Prüfschlange abarbeiten, Experimente für K, endliche Bausteine ((K2), Zweig-Zertifikat bis 400, (A)),
Beweis-Buchhaltung (`ABHAENGIGKEIT.md`), kritische Gegenstimme. Kennzeichnung: bewiesen / Skizze / numerisch / Vermutung; Statuswörter der
Buchhaltung V/S/B (siehe ABHAENGIGKEIT.md). Zeiten aus `date`.

## 0.0 Abschlussnotiz Phase 3 (12.10. 14:15; Team auf Entscheidung des Autors angehalten — keine Läufe mehr, lokal und auf aicgolling geprüft)
**Verifiziert (Logik gegengelesen, mit eigenem Code nachgerechnet):** Satz P-V Schritte 3–6 (§1.1); Gs G8–G44 (§1.2); Ks K1, K1′/K1″, §18, K10/K11, K19, §32, alle Sätze von
HAUPTSATZ_KR0 §7 (§1.3, §1.6, §1.9–§1.12); Ps p32a/rest_schranke (§1.7); Ms 𝒯-Weg (a″)_∞, σ_min ≥ 0.641 (§1.8); Ks (CD1) exakt, kubisch + Kneser, 1e-89 (§1.13).
Streng (Status S/V): (K2) C₆ = 461.36 (§2.2); Zweig-Zertifikat 193 → 400 und Kopplung (A) sup_{D̄}|p₄₀₀ − f_K| ≤ 1.75e-12 (§2.1); Sprungschranke C_J = 1.0378 bis u = 0.48 (§3.15).
**Gemessen (numerisch):** ‖Ĵ_N⁻¹‖ = 1.251155 für N = 100 … 2000 (§3.1); (S_γ) am echten p_N bis 2000 (§3.4, §3.9b); T1: Mechanismus real, M4-Plateau 1.2–2.3e-6
in der Γout-Zone für N = 1000 … 3200 (§3.14); (R_L) |K_n|/n^{0.894} → 0.599, C auf 5 Stellen (§3.17); E3 B_n = 0.626·Q_n·n^{−0.053} (§3.17); E2 η_N/r_N = 5e-5 (§3.18);
s_E = 3.9–4.1 (80°/−90°) (§3.19); F3 mit P selbst falsch (§3.20); E4′ Q_n/|P_n| ≈ 1.008·R_P, R_P ~ n^{−3.1} (§3.22); Primalrand: Exponentialphase endet bei n ≈ 2200 (§3.25);
CF ungewichtet Steigung ≈ 0.90(S − 0.0147) (§3.24); Kopf×Kopf am Steg Summe +0.00006 (fallend), |m·ΔK_n(t,t)| → 0.956 (§3.26a); **gewichtetes CF ≈ 2.95e-16 bei N = 2000,
~N^{−4.7}, linsenunabhängig (S = 0.0126 … 0.03, τ₁ = 4° … 12°)** (§3.26d); Adjungierte auf der K-Form-Linse S = 0.0137: Reserve am Steg (−0.0011), Wachstum nahe L (≤ 8°) (§3.26c).
**Offen:** Prüfstellen HAUPTSATZ_KR0 §5: Nr. 1 (Buchhaltung (G-P)/Budget), Nr. 3 (H2 ohne sup_{D̄} im Wortlaut), Nr. 5 (Startbedingung s′ < s_E) nicht gegengelesen;
Nr. 2 mit Zusatzbedingung P̃(𝒞) ⊂ E′ (§1.14); Ks (BU)-Messungen (R_P/R_A auf S₀ = 0.0125 bis 3000, gewichtetes CF mit π_n) nicht begonnen; V2-Keil zurückgestellt (§3.16).
**Daten:** Skripte v01–v60 und Logs hier (sg/, cf/, cfw/, cd1/, e4/, sE/, t1/, st/); Server aicgolling ~/v_pruefstand/sg/a/b (v33_pN_1500/2000.json = p_N, v48/v51/v56/v58/v59-Logs,
v59_ell_*.json = adjungierte Köpfe ℓ_n als exakte Vektoren bis n = 2000, v59_w_*.npz = Gewichte auf γ), ~/v_pruefstand/t1, zert400, k2. Fallen: Präzision der letzten
Eliminationsschritte (n = 2000 braucht ≥ 1700 Bit), acb-Horner bläht Radien um ≤ √2 je Schritt (+n/2 Bit), Prozesse nur über PID/PGID beenden.

## 0. Stand in Kürze
- 19:55 Zweig-Zertifikat 193 → 400 auf aicgolling gestartet (`zert400/`, Kopie von Xs zert/, venv ~/halbexp-venv, 12 Arbeiter, SID 2892380,
  Server-Ordner ~/v_pruefstand/zert400).
- 20:01 Landkarte Teil B geschrieben (`../LANDKARTE.md`), Abhängigkeitsbild `ABHAENGIGKEIT.md`.
- 20:01 Test: reicht für Ms (a″)_∞-Satz eine (K2)-Schranke C₆ ≤ 5e5? (m26 mit „pzert:5e5:6“, `m_kopie/v01_m26_C5e5_p6.log`).
- 20:20 Satz P-V Schritte 3–6 verifiziert (§1.1); 20:25 Gs Sätze G8–G44 (relevante) verifiziert (§1.2); 20:27 **C₆ ≤ 5e5 genügt** (Fenster 6/14,
  σ_min ≥ 0.5695, §2.2); 20:51 **Ks Satz K1 verifiziert** (§1.3); 20:52 Ideen bewertet, Plan bestätigt (../PLAN.md, A5 neu).
- 20:55 (K2)-Lauf auf aicgolling gestartet: `pv_kopie/v06_k2_rand.py` (Reparatur von Ps p32), 19 Abschnitte (`k2_jobs.txt`), 8 Arbeiter,
  Server-Ordner ~/v_pruefstand/k2, SID 2892826, Logs log/k2_<seite>_<i>.log, Ergebnisse k2_<seite>_<i>.json. Erwartung ≈ 1 h (Abschnitte
  nahe L ≈ 30 min). Danach: Zweigtreue an den Abschnittsgrenzen prüfen (erste/letzte z-Scheibe in den JSON-Dateien; Punkt-Krawczyk am
  Grenzpunkt wie in v06), C₆ = 2·Σ(∫ + Rest)/(2π) bilden, Ziel ≤ 5e5; dann Datei im Format {"K","b_mitte","b_radius","C","p"} für M/O.
- 21:13 **(K2) zertifiziert: C₆ = 461.36** (§2.2), Datei `kneser_halbexp_b_zert_V.json`; 21:15 Möbius-Modell (A2) negativ (../REFLEXION.md);
  21:17 Ks Lemma K1″ verifiziert (§1.6).
- **Laufend auf aicgolling (Stand 21:18):** (i) Zertifikat 193 → 400 (60 von 207 Übergängen fertig, alle ok; Kette am Ende automatisch,
  log/kette.log); (ii) Ms m26 mit der V-Datei, nacheinander Fenster 7/23 und 6/14, 4 Threads (~/v_pruefstand/mrun/x/y/m_kopie/v10_*.log, v11_*.log);
  (iii) Gegenprobe der (K2)-Abschnitte bis δ_min = 1e-40 (rechts_8, links_9). Ein erster m26-Versuch lief mit numpy auf allen Kernen
  (Last 41) — abgebrochen und mit OMP/OPENBLAS_NUM_THREADS = 4 neu gestartet. Lehre: auf aicgolling numpy-Threads immer begrenzen.
- 21:20 Auftrag der Leitung: (S_γ) zuerst. 21:24 Ms m26 mit V-(K2)-Datei: σ_min(J_∞|P_N) ≥ 0.641 (7/23) bzw. 0.5696 (6/14) (§1.8);
  21:29 p32a und rest_schranke verifiziert (§1.7); 21:33 X4a-Lücken eingeordnet (§1.5); **21:45 (S_γ) unabhängig bestätigt** (kubisch bis N = 1500,
  Kneser bis 639), Schwellengesetz κ = 1.00, Gegenbeispiel-Profil (§3.1).
- **Laufend (21:47):** Zertifikat 193 → 400 (120/207 Übergänge fertig, alle ok; Kette danach automatisch, log/kette.log, Erwartung ≈ 22:05);
  Gegenprobe (K2)-Abschnitte rechts_8/links_9 bis δ_min = 1e-40 (redundant, enden spätestens 23:55).
- 22:02 Reserve-Fenster U/Λ ≥ 14 (§3.2, LANDKARTE C10); 22:05 Kette 193 → 400 unabhängig geprüft, (A) bei N₀ ∈ [318, 400]: sup ≤ 3.6e-12 (§2.1);
  22:15 (K2)-Neulauf mit Ks F1/F2: C₆ ≤ 461.354 unverändert ⇒ (K2) **V** (§2.2); 22:31 Ms Lemmata M7–M10, M13 verifiziert ⇒ (a″)_∞ Kneser **V** (§1.8);
  22:23–22:36 Familie in w-Form, Neumann-Grössen (§3.3): naive Neumann-Reihe scheitert (e^{0.18N}), wahre Änderung der Inversen linear in ε, aber
  bei Störungen hoher Grade wächst Δ_inv/ε ≈ e^{0.0127N} (Eckraten-Grösse).
- 23:14 Block 3 (Auftrag Leitung 22:52): (S_γ) am echten p_N bis N = 1000 (§3.4), gewichtete Normen (§3.5), Ks A3 (i)–(iv) (§3.6, §3.7),
  PRUEFLISTE-Notiz (§5), „Stand des Gesamtbeweises“ in ABHAENGIGKEIT.md.
- 01:59 Block 4 (Aufträge 23:20/23:45/00:40): Familie längs H1/H2 bis N = 1000 (§3.8), Satz V2 mit B_max = 30.8 / 76.7 (§3.9), (S_γ) am echten p_N bei
  N = 1500/2000 für drei Linsen (§3.9b), Datenkorrektur (§3.10); Runde 2 Ecke: 6 Ideen, VE1/VE2 abgebrochen, VE4 positiv (IDEEN, PLAN).
- 10:52 Block 5 (Aufträge 08:40/08:58/09:30): Sprungschranke streng C_J = 1.0378 bis u = 0.48 (§3.15); V2-Keil zurückgestellt (§3.16); T1 entschieden (§3.14: Mechanismus real, Plateau in einer Γout-Zone im Modell M4); (R_L)/E3 bis 2000 (§3.17); E2 (EB_N) bis 2000 (§3.18). Laufend: M4 N = 3200 (sg/a/b/t1/E_M4_s1.5_f.log).
- 13:01 Block 6/7 (Leitung 11:40/12:00): §7-Sätze gegengelesen (§1.12), s_E/F3/E4′/f_K-Geometrie (§3.19–§3.23), CF-Test (Fam-koh) (§3.24), Primalrand bis 3000 (§3.25). Keine Serverläufe offen.
- 13:35 Block 8 (Leitung 13:00/13:08/13:20): Übergabenotiz §0.0; (CD1) gegengelesen und mit eigenem Code bestätigt (§1.13). Laufend auf aicgolling seit 13:11:
  `run_v59.sh` (PGID 2928683, ~/v_pruefstand/sg/a/b): v59c-Koeffizienten (fertig), v59_elim kubisch/Kneser N = 2001 (ℓ_n als Vektoren), danach v59_eval je Linse.
- 14:13 Block 8 abgeschlossen: Kopf×Kopf am Steg, CF im Fenster, Adjungierte auf der ganzen Linse, gewichtetes CF (§3.26), Prüfstellen 2/4 (§1.14). Team angehalten;
  keine Prozesse mehr (aicgolling und lokal geprüft).
- **Abholen:** `aic 'cd ~/v_pruefstand/k2 && tar -cf - log k2_*.json' | tar -C pv_kopie/lauf -xf -` bzw. für das Zertifikat
  `aic 'cd ~/v_pruefstand/zert400 && tar -cf - log boxen erg kette_193_400.json' | tar -C zert400/lauf -xf -` (Kette: log/kette.log).

## 1. Prüfprotokolle

### 1.1 Satz P-V, Schritte 3–6 (p27, p28, p29, kc_rig, kc_kneser) — Ergebnis: **verifiziert (Logik), mit drei formalen Lücken ohne Wirkung**
Gelesen: p27_kneser_zert.py, p28_qk_identifikation.py, p29_fK_koeffizienten.py, kc_rig.py, kc_kneser.py, die Logs; Paper VI Thm 11.11
(gp-tetration/papers/paper6, wörtlich: F₁, F₂ holomorph auf H_{−2}, F(z+1) = e^{F(z)}, F(0) = 1, für jedes x > −2 beide Grenzwerte L, L̄ ⇒ F₁ = F₂).
Ps Prüfstellen (ARCH §9.8) einzeln:
1. **Periodizität von F₀ (p28 Schritt 4): richtig.** An der Naht Re z = ½ ist χ ≡ 1 auf [0.467, 0.5], also F₀ = T̃∘(id + d(·−1)) auf BEIDEN
   Seiten (rechts über F₀(z) = e^{T̃(z−1)} = T̃(z + d(z−1))) — holomorph über die Naht hinweg. Oben: Ũ(z+1) = e^{Ũ(z)} exakt (P(Lw) = e^{P(w)},
   θ̃ periodisch); Y ist 1-periodisch wegen e(ζ) − e(ζ+1+d(ζ)) = d(ζ), und diese Identität folgt aus Ũ(ζ+1+e(ζ)) = e^{T̃(ζ)} = Ũ(ζ+1+d(ζ)+e(ζ+1+d(ζ)))
   und der Injektivität von Ũ auf Scheiben vom Radius < m_U/M_U ≈ 0.47 (|Ũ'| ≥ 0.83, |Ũ''| ≤ 1.78). Auf 0.35 ≤ Im z ≤ 0.40 stimmen beide
   Definitionen überein (ψ = 1). μ₀ = ∂̄W/(1+∂W) bzw. ∂̄V/(1+∂V) braucht T̃' ≠ 0 nicht. ✓
2. **Beltrami (p28 Schritt 5): richtig.** k aus Kettenregel (∂̄(χd) = χ'd/2, ∂̄(e∘g) = e'(g)∂̄g, ∂̄ψ(y) = iψ'/2) nachgerechnet; |e'|, |d'| per Cauchy
   mit Rand 0.05 liegen in G_B bzw. G_S. Zylinder: q = e^{2πiz} erhält |μ|; Hauptlösung ψ₀ = q + C[h], und die Normierung ψ₀(q) − ψ₀(0) erzeugt
   genau den Faktor 1/(|ζ||ζ−q|) in ε (Hölder 4/(4/3), Umordnung auf die flächengleiche Scheibe); |φ − id| ≤ ε/(2π(1−ε)) aus |log(1+δ)| ≤ ε/(1−ε). ✓
3. **Identifikation (p28 Schritt 6): richtig, Voraussetzungen von Thm 11.11 erfüllt.** Linksfortsetzung auf H_{−2} braucht genau A auf
   [−½, ½] und B auf (0, ½] (nachgerechnet: zweite Log-Stufe ⇔ F ∉ (0,1]); die Verklebung an Re z = −½ braucht zusätzlich |Im F| < π dort —
   folgt aus A bei Re = ½ (e^{F} ∉ (−∞,0)) plus Stetigkeit ab Im F = 0 auf ℝ und dem Schwanz Im ≥ 3. Überdeckung: Gitter 0.035, Radius 0.025 ≥ halbe
   Diagonale; x-Bereich [−0.525, 0.505], y [−0.025, 3.03], unten per Symmetrie; Grenzwerte für x < −½ über Log L = L. ✓
4. **p27 (Defekte): richtig.** Maximumprinzip auf den Rechteckrändern (Abstand der Scheibenmitten 1.8·RT < 2·RT), Taylor-Modelle mit
   Cauchy-Schwanz, Ableitungsschranken mit Cauchy auf 2r_z; skalarer Kantorovich (h = ηM₂/m² ≤ ¼) mit Radius 2η/m ≤ s₀/2. ✓
5. **kc_rig._restkonstante: richtig.** Rekursion p_n = Lq/(L^n − L) aus P(Lw) = e^{P(w)}; Fixpunktgleichung R(w) = δ(w) + e^{T(w/L)}(e^{R(w/L)} − 1)
   in der Norm sup|R|/|w|^{n₀+1} auf |w| ≤ ½, δ-Schwanz per Cauchy auf |w| = 2 (ganze Funktion), Lipschitz < 1; Gleichheit mit P über formale
   Eindeutigkeit (L^n ≠ L). U_taylor: Koeffizienten von F₁ = exp^k∘T_{n₀}∘(w/L^k) mit ctx.cap = m+1 (kein 10-Term-Fehler), Schwanz per Cauchy,
   Wertfehler P − T_{n₀} mit Lipschitz von exp über konvexe Bälle. ✓
6. **p29 / kc_kneser (Fortpflanzung von ε): richtig** (ε_T, ε_U über |φ⁻¹ − id|, |τ| ≤ 2(|T̃(0)−1| + lip·|φ−id|), d_max, η_B; der Zuschlag ε·100
   für Ableitungen ist gedeckt, weil ε auf einer 0.01-Umgebung jeder benutzten Scheibe gilt: Band bis Im 0.45, Re bis 0.6). Aliasformel
   M r'^{−n} q/(1−q) richtig; Bogenüberdeckung auf |x| = 1.30 lückenlos (Abstand ≤ 0.99(h+h_n)). ✓
7. **Log-Zweige in sexp_ball:** Hauptzweig, konsistent mit der Fortsetzung über Prüfung A; flints Log liefert über dem Schnitt nur weite Bälle (konservativ). ✓

**Formale Lücken (ohne numerische Wirkung, Ränder ≥ 1e-3 gegen Fehler ≤ 1e-166):**
- (L1) p28, Prüfung B nahe 0: das analytische Argument benutzt G(0) = 1 für G = F(·+τ), das Gitter prüft aber F; der Ableitungsfehler wird mit
  „10·(d_max + |φ−id|)“ angesetzt, korrekt wäre per Cauchy ≈ 20·(Wertfehler) plus |τ|·sup|F''|. Behebung: Fehlerbudget um |τ|·sup|F'| und den
  Cauchy-Faktor erweitern (eine Zeile).
- (L2) p28: „Lipschitz-Schranke 10 für |∇F₀|“ wird nicht in p28 hergeleitet; p29 liefert sup|T̃'| ≤ 1.95, sup|Ũ'| ≤ 1.69 auf den benutzten Gebieten
  (plus |∇W₀|, |∇V₀| ≤ 1e-170) — damit gedeckt, Verweis fehlt.
- (L3) p29: die Krawczyk-Nullstelle s_j von sexp_K(s) = x_j an den 2049 DFT-Knoten ist eindeutig in ihrer Scheibe, aber dass sie der **Hauptzweig**
  slog(x_j) ist (Fortsetzung längs |x| = 1.2 ab x = 1.2), wird nur durch Newton-Fortsetzung mit Schritt 0.0018 gestützt, nicht zertifiziert
  (ebenso im Bogenteil). Behebung: x-Scheiben über die Bögen zwischen den Knoten (wie im Bogenteil, xrad = halber Abstand) und Prüfung, dass
  aufeinanderfolgende w-Scheiben sich schneiden ⇒ Fortsetzung bleibt in der Scheibenkette. Numerisch gestützt durch Ms Abgleich mit Strang As
  unabhängiger Datei (1.1e-43 gewichtet) — ein Zweigsprung gäbe grobe Abweichungen.
- (L4) bekannt (M): Radien der Koeffizientendatei ohne Dezimalrundung der Mitten (≤ 1e-139·|a_n|).
**Folgerung:** Satz P-V gilt als **V** (verifiziert) mit den Anmerkungen L1–L4; die Behebungen sind Routine. p32a ((K3)) ist noch ungelesen.

### 1.2 Strang G — die für den Eckweg relevanten Sätze (11.10. 20:25) — Ergebnis: **alle geprüften richtig; sie schliessen keine Lücke**
- **G9 (G15) Rang-1-Zuwachs:** Blockelimination am geränderten System nachgerechnet (x_n = (g_n − rᵀA⁻¹g')/s_n, x' = A⁻¹g' − A⁻¹b·x_n). Voraussetzung
  „linker oberer Block von J̃_n = J̃_{n−1}“ gilt, weil f, m fest sind. ✓ **G16/G17/G18:** Zuwachs = Q_n·(e_nᵀJ̃_n⁻¹g) = Q_n(A_nᵀg)/h_n ✓, Teleskop ✓.
- **G10 Summenkriterium:** Dreiecksungleichung + Σ n^{−1/2} ≤ ∫ t^{−1/2} ✓; Scheibenfassung per Maximumprinzip (𝒟_n analytisch in x ≠ 0,
  Scheiben konzentrisch und schachtelnd) ✓. Trivial, aber korrekt; die Last liegt ganz in der All-n-Hypothese (G19).
- **G11 Kreisabbildung:** |ρ/w − L| = r ⇔ |w − c|² = R² mit c = ρL̄/(ρ²−r²), R = ρr/(ρ²−r²) nachgerechnet ✓; Majorante = Dreiecksungleichung ✓.
- **G22 funktionale Schur-Rekursion:** aus J_{n−1}⁻¹ = Σ_{j<n} ιQ_j(ιA_j)ᵀ/h_j (iteriertes G16) folgen z_n, w_n, h_n und (G22b) ✓.
- **G23 Pullbacks:** Scheibenschachtelung unter w = ρ/x, ‖(a+bζ)^k‖_W ≤ (|a|+|b|)^k ✓.
- **G8:** Eigenvektoren (λ², λ, 1, 0) aus (t−a)λ³ = bλ² + cλ + d ✓ (Monodromie-Zahlen nicht nachgerechnet, für den Beweis unerheblich).
- **G39/G40:** p_{j+2} = (F+z)p_{j+1} − zFp_j für beide Summanden ✓; Spalten von R_n: j ≤ n−2 null, j = n−1: −T_np_{n+1}, j = n: T_n(zFp_n) ✓
  (T_n vertauscht mit Multiplikation durch Potenzreihen als untere Dreiecksmatrix ✓).
- **G44 (F = 2z + ηz², M = F'∘F):** Rouché (17/250) ✓, Subordination ✓, ‖𝕋|_V‖ ≤ (43/20)(r/R)⁴ = 0.174 ✓, Kopfmatrix: Einträge 2η, 3η/4, −η²/8, η/4 aus
  H = z/2 − ηz²/8 + η²z³/16 + … nachgerechnet ✓, Frobenius 0.0260 < 1/16 ✓, Blockinversion √3881/16 < 4 ✓, |[zⁿ]u_k| ≤ 4(2/3)ⁿ(4/5)^k ✓.
  **Bewertung:** F(0) = 0 macht J_n unteres Dreieck — es gibt KEINEN Abschneide-Effekt (J_n⁻¹ = Abschnitt von J_∞⁻¹). G44 ist richtig, aber
  strukturell das triviale Modell (Ks A3.5: exakte Lösbarkeit und Ecke schliessen sich aus); als Schablone für die Ecke ungeeignet.
- **G43:** Neumann-Reihe ✓ (abstrakte Implikation).
**Fazit G:** Gs Sätze sind korrekte, meist elementare Umformungen; der Beweisinhalt steckt in den offenen All-n-Hypothesen G19/G24/G38
(B_n ≤ 0.75, ‖S_n‖ ≤ 4) und in der Numerik. Für die Buchhaltung: G9–G11, G22, G23, G39, G40, G43, G44 → **V**.

### 1.3 Ks Satz K1 (Kurvenkriterium, k_kern/BERICHT §1) — Gegenlesen 11.10. 20:51 — Ergebnis: **verifiziert (Logik)**
- **K1.0:** [x^k]W = (1/2πi)∮_γW(t)t^{−k−1}dt für W ∈ A(Int γ) ∩ C(Int γ‾), also T_NW = (1/2πi)∮W(t)T_Nk_t dt; (T_Nk_t)_kρ^k = ẽ_t/t ⇒ ‖r_t‖ = G_N(t)/|t| ✓.
- **(i):** W = J_∞k_η = 1/(η − f) + m/(η − x) ist analytisch auf Int γ, sobald η ∉ K_γ = Int γ‾ ∪ f(Int γ‾) und m auf Int γ analytisch ist
  (f(Int γ‾) trifft den Schlitz von f' nicht — Ks Geometrie 1.4); |W| ≤ (1 + M_γ)/dist(η, K_γ) auf γ ✓; sup_{D̄} ≤ √(N+1)·ℓ²_ρ ✓.
- **(ii):** y_η = Π_Nk_η − T_Nk_η: Koordinaten rational in η, Pole nur 0 und a₀ ([x^k]1/(η−f) hat nur den Pol a₀, [x^k]m/(η−x) nur 0), und
  y_η = Σ_{j>N}Π_N(x^j)η^{−j−1} = O(η^{−N−2}) ✓ ⇒ Bernstein–Walsh auf Ĉ∖K_γ^{1/N} ✓; g_{1/N} ≥ g_γ − sup_{∂K^{1/N}}g_γ (Maximumprinzip) ✓;
  Lipschitz von g_γ: an L, L̄ ist der Aussenwinkel der von K' (141.69°, α' = 0.787 < 1 ⇒ g ~ r^{1.27}), sonst glatter Rand ✓.
- **Schlitz und Linse:** S_L zeigt in Richtung der äusseren Normalen von ∂D; lokal ist f(D̄) die um arg μ = 38.3° gedrehte Halbebene, und
  Re(uL·L̄/μ) = uρ²Re μ/|μ|² > 0 ⇒ S_L ∩ f(D̄) = ∅ nahe L ✓; die Linse hat bei θ_L die Dicke 0 und trifft den Schlitz (fester Winkel θ_L) nicht ✓.
- **Bewertung:** Satz K1 ist richtig und ersetzt W1'' + D6 + (E_sing) + Q_N durch EINE Aussage (H*_G) sup_γ G_N ≤ C_GN^κ, gültig für f_K.
  Damit ist Landkarte B6.2 (Kneser ≠ kubisch) **erledigt**. Numerisch (k01, Kneser = kubisch) gilt (H*_G) mit κ ≈ 0.45.
- **Nötige Stärke (Korrektur zu Ks „κ < 2.79“):** für (LP) κ + 3/2 < 4.29 ⇒ κ < 2.79 ✓; **für die Familie über X4a (θ_N = A_θNR_N → 0)
  aber κ + 3/2 + 1 < 4.29 ⇒ κ < 1.79**; gemessen 0.45 ✓. Für N₀ ≤ 400: Eck-Polynomfaktor κ_E = κ + 3/2 ≈ 2, C_E ≈ A_γ(1 + M_γ)C_G ≈ 3.3·0.53 ≈ 1.8
  (C_G aus ln G_1200(L) = 2.59) — also etwas über Xs Szenario C_E = 1; die 3/2 sind Buchhaltung und sollten gedrückt werden.
- **Offen (endlich, V):** Geometrie für f_K zertifizieren (Linse ∩ Schlitz = ∅, f(ℬ) ⊂ D, f(D̄ ∪ ℬ̄) ∩ Schlitz von f' = ∅, Lipschitz-Konstante
  von g_γ); (F_{K_γ}) braucht Holomorphie von f_K auf einer Umgebung der Linse (dünn genug ⇒ in Bs L8.4-Gebiet).
- **Familie:** dieselbe Aussage (H*_G) für die Sekantenoperatoren T_N𝒥_P (G^P_N auf γ) — Störung r^P_t − r_t = −(T_N𝒥_P)⁻¹T_N(𝒥_P − J_∞)r_t; braucht
  r_t auf f(Γout) ⊄ D̄ (Ks Beobachtung r_t → s_t auf K_γ wäre genau die nötige Überkonvergenz).



### 1.4 Zweig-Zertifikat: zert_kern.uebergang (= Ps p24 + Xs NC/Newton-Änderungen), zert_job, zert_kette — Ergebnis: **verifiziert (Logik)**
- Röhren-Krawczyk für G_{s,c}(a) = T_M[a∘a − eˣ] − (1−s)c·x^M, für alle s ∈ [0,1] und c ∈ C zugleich: K − z̃(s) = −YG(z̃(s)) + (I − YJ(z̃(s)))B
  − Y(J(Z) − J(z̃(s)))B ⊂ int B ⇒ in jeder Röhrenscheibe genau eine Nullstelle, J invertierbar ⇒ Pfad stetig in der Röhre ✓.
- s ↦ G(z̃(s)), J(z̃(s)) sind Polynome in s; Werte an NC Tschebyschow-Knoten in Punktarithmetik, Fortsetzung über Lebesgue-Konstante
  Λ ≤ (2/π)ln NC + 1 und Restformel |ω| ≤ 2^{1−2NC}·M_R/(R−1)^{NC} (Cauchy auf |s| ≤ R, M_R über Majoranten-Reihen) ✓.
- J(z̃ + b) − J(z̃) mit Majoranten: Einträge [x^k]p^j + [x^{k−j}]p'∘p sind Polynome mit nichtnegativen Koeffizienten in den a_k ⇒
  |Φ(z+b) − Φ(z)| ≤ Φ(|z|+|b|) − Φ(|z|) ✓. Kopfkoeffizient c als Intervall über die Box von p_N (C.rad() im Rest) ✓.
- Start: pad(p_N^exakt) ist Nullstelle von G_{0,c_exakt} und liegt in der Röhre bei s = 0 (start_ok für die ganze Box B_N) ✓; Ende:
  Krawczyk-Verengung bei s = 1 auf Z_b ⊂ Röhre ⇒ die Pfad-Endstelle ist die eindeutige Nullstelle in Z_b ✓. Kette: K_{N+1} ⊂ B_{N+1}
  (Newton-Mitten nur Vorschlag, die Inklusion entscheidet) und P-Box p_193 ⊂ B_193 ✓.
- **Formal (wie P selbst vermerkt):** die Knoten sind auf 2^{−prec} gerundet; |ω| ≤ 2^{1−2n} und die Λ-Schranke gelten für exakte
  Tschebyschow-Knoten, für die gerundeten nur bis auf relative ~1e-100 — gedeckt durch die Reserve ≈ 0.5, aber nicht ausgerechnet.
  Behebung: |ω| und Λ für die tatsächlichen Knoten direkt in Intervallarithmetik bilden (eine Zeile).

### 1.5 X8′ und X4a (Familie, eckfreier Teil) — erstes Lesen 11.10. 21:02 — Ergebnis: **X8′ Hauptterm richtig; zwei Lücken in X4a/X8**
- **X8′ Hauptterm richtig:** δ'∘f = [(J_∞δ)' − m'δ − mδ']/f' (aus (J_∞δ)' = f'·δ'∘f + m'δ + mδ'), m'/f' = f''∘f, ‖δ'‖ ≤ (N/ρ)‖δ‖,
  ‖(J_∞δ)'‖ ≤ (KN/ρ)‖J_∞δ‖ + Schwanz (X5), ‖δ‖ ≤ ‖J_∞δ‖/c₀ ⇒ A_θ = (K/ρ)/min|f'| + max|m/f'|/(ρc₀) + M2/(c₀N) ✓ (13.4 bei c₀ = 0.54; 13.2 bei 0.5695).
- **Lücke 1 (Skizze in X8):** der Verschiebungsterm δ'(f + s(P−f)) = Σ_j δ^{(j+1)}(f)t^j/j! braucht ‖δ^{(j)}∘f‖ ≤ (A'N)^j‖J_∞δ‖ für alle j mit
  explizitem A'; „rekursiv über Faà di Bruno, A' ≈ A_θ“ ist nicht ausgeführt. Wegen N·R_N ~ N^{−3.8} harmlos, aber zu schreiben.
- **Lücke 2 (Schlitz):** für P ∈ B_N kann P(D) den Schlitz S_L in einem Stück der Länge ~R_N/(0.785ρ) treffen (f(D̄) hat bei L zum Schlitz nur den
  Abstand 0.785ρu), und L kann in P(D) liegen. Dann ist f∘P (also M_P) nicht holomorph auf D, und die L²(∂D)-Abschätzung von (M_P − m)δ ist
  so nicht gerechtfertigt. Richtig dürfte sein: nur Jets bei 0 zählen; die Taylor-Koeffizienten von M_P bis 5N sind durch den Sprung
  C_J·R^{5.447} und (1+O(R))^{5N} ≈ 1 kontrolliert (R ≪ 1/N). Muss in X4a ausgeschrieben werden, bevor X4a als Fundament zählt.

**Einordnung der zwei Lücken (Auftrag Leitung 21:20) — schliessbar? umgeht K1 sie? muss die Familie neu gebaut werden?**
- **Lücke 1 ist Routine und schliessbar.** δ ist ein Polynom, also ist δ'(f + t) = Σ_{j<N}δ^{(j+1)}(f)t^j/j! exakt. Aus F = δ∘f = J_∞δ − mδ folgt
  durch wiederholtes Ableiten (Faà di Bruno): δ^{(j)}∘f = f'^{−j}[(J_∞δ)^{(j)} − Σ(niedrigere δ^{(i)}∘f mit f''…f^{(j)}) − (mδ)^{(j)}]. Induktiv mit
  ‖(J_∞δ)^{(i)}‖ ≤ (KN/ρ)^i‖J_∞δ‖ + Schwanz (X5) und ‖δ^{(i)}‖ ≤ (N/ρ)^i‖δ‖ ergibt sich ‖δ^{(j)}∘f‖ ≤ (A'N)^j‖J_∞δ‖ mit A' = A_θ(1 + O(1/N)).
  Wegen N·R_N ~ N^{−3.8} tragen nur j ≤ 2 merklich bei. Rechenaufwand: eine Seite, keine neue Idee.
- **Lücke 2 ist im x-Rahmen schliessbar, im K1-Rahmen braucht sie eine kleine neue Überlegung.** Für T_{5N}𝒥_P zählen nur Jets bei 0;
  liegt der Verzweigungspunkt P⁻¹(L) von f∘P im Abstand ~R von L (möglicherweise in D), so wachsen die Koeffizienten von M_P erst ab k ~ 1/R ≫ 5N
  merklich ((1 − R/ρ)^{−5N} ≈ 1 + 5NR/ρ), und der Verzweigungsanteil hat Grösse ~C_J R^{5.447}. X4a muss dafür von L²(∂D)-Funktionswerten auf
  Koeffizientenschranken umgestellt werden. Im K1-/K1″-Rahmen dagegen braucht die Cauchy-Darstellung für 𝒥_P, dass M_P auf ganz Int γ analytisch
  ist — liegt P⁻¹(L) in Int γ (die Linse ist bei L tangential, Int γ ist dort wie D), versagt das formal. Ausweg: eine von P abhängige Kurve γ_P
  durch P⁻¹(L) statt durch L (Verschiebung O(R) ≪ Eckskala N^{−0.894}) und der Nachweis, dass (S_γ) bzw. (H*_G) unter O(R)-Störungen der Kurve
  stabil ist; oder der Nachweis L ∉ P(Int γ‾) für P ∈ 𝒮_N ∩ B_N.
- **K1 umgeht die Familienfrage NICHT.** X4a liefert nur den eckfreien (rechteckigen) Teil: ‖(𝒥_P − J_∞)δ‖ ≤ θ‖J_∞δ‖ ist eine Störung relativ zum
  VOLLEN J_∞δ, nicht zu seinem Abschnitt T_NJ_∞δ; für quadratische Abschnitte (und genau die braucht K1″(a) über ‖(Ĵ^P_N)⁻¹‖) folgt daraus nichts,
  weil ‖(J_∞δ)∘k‖_{H²(𝔻)} ≫ ‖Ĵ_Nδ‖ sein kann (J_∞δ auf f(Γout) ⊄ D̄). In der Sprache von K1 heisst die Familien-Aussage: **(S_γ) bzw. (H*_G) für
  T_N𝒥_P gleichmässig über P ∈ 𝒮_N ∩ B_N** — das ist (Fam-Eck) in neuer Gestalt.
- **Was neu gebaut werden muss:** die Familien-Erweiterung im K1″-Rahmen. Zwei Wege: (i) (S_γ) als strukturelle Aussage beweisen, die f nur über
  stabile Grössen benutzt (Fixpunkt nahe dem Kreis, μ, Univalenz, Randkurve) — dann gilt sie für alle Sekantenoperatoren (Idee V5); (ii) Störung
  über Xs exakte Formel Π^P_Ng − Π_Ng = −(T_N𝒥_P)⁻¹T_N[(𝒥_P − J_∞)(Π_Ng − g)] mit der kohärenten Eingabe E_N = Π_Nk_η − k_η, die auf K_γ klein ist —
  braucht Kontrolle von E_N auf f(Int γ) ⊄ D̄ (Überkonvergenz r_t → s_t, Ks Beobachtung k03). Numerischer Prüfstein für beides: ‖(Ĵ^P_N)⁻¹‖ für
  P = p_N (echte Lösung) und für zufällige Störungen gleicher Grösse (Plan A3).

### 1.6 Ks Lemma K1′/K1″ (konforme Gleichungskoordinate, k_kern/BERICHT §3, §5.1) — Gegenlesen 11.10. 21:17 — **verifiziert (Logik)**
- Ĵ_Nx = T_N^w[(J_∞x)∘k]: C_k unteres Dreieck wegen k(0) = 0, also T_N^wC_k = T_N^wC_kT_N ✓.
- (a) Π_Nk_η = Ĵ_N⁻¹T_N^w[W∘k] (W − T_NW = O(x^{N+1}) ⇒ ∘k = O(w^{N+1})), ‖T_N^w[W∘k]‖ ≤ ‖W∘k‖_{H²(𝔻)} ≤ sup_γ|W| (k stetig bis ∂𝔻) ✓.
- (b) ⟨u, Ĵ_N⁻¹c⟩ = ⟨y_u, T_N[c∘k⁻¹]⟩ = (1/2πi)∮_{|x|=ρ}(c∘k⁻¹)Y_u dx/x (Terme x^{j>N} paaren mit Y_u zu 0); Verschiebung |x| = ρ → γ ist
  zulässig (Zwischengebiet = Linse, auf Γout fallen die Konturen zusammen; c∘k⁻¹ stetig bis γ, Y_u analytisch auf ℂ∖{0}); Cauchy–Schwarz mit
  |dt| = |k'||dw| ⇒ c_k = (sup|k'|/2π)^{1/2}/min_γ|t| ✓.
- (c) 1/(k(w_t) − k(w)) = q(w)/(w_t − w), q analytisch und |q| ≤ Q_k, ‖L_N(q)‖ ≤ ‖q‖_∞, ‖T_N^w[1/(w_t − ·)]‖ = √(N+1) ✓.
- Voraussetzung: γ in L, L̄ C¹ mit Dini-stetiger Tangente (Linse tangential, s ≍ |α − arg L|^p, p > 1) ⇒ Kellogg–Warschawski: k' stetig, ≠ 0 ✓
  (Satz zitiert, nicht nachgelesen — X §5).
- **Bewertung:** richtig. Damit ist die Zielaussage auf **(S_γ) sup_N‖Ĵ_N⁻¹‖ < ∞** verschärfbar (Eck-Polynomfaktor 3/2). (S_γ) ist eine Aussage über
  die Abschnittsmethode einer FESTEN Operatorfamilie Ĵ_∞ = C_kJ̃_∞ mit Standardprojektionen — dafür gibt es Literatur (Gohberg–Feldman,
  Böttcher–Silbermann, Roch–Santos–Silbermann: Stabilität ⇔ Invertierbarkeit gewisser Grenzoperatoren). Ks Messung (kubisch, N ≤ 800):
  ‖Ĵ_N⁻¹‖ = 1.251 konstant. **Unabhängige Nachrechnung durch V steht aus** (eigene konforme Abbildung, andere Linse, Kneser).

### 1.7 Ps p32a (θ_K, (K3)) und p32.rest_schranke — Gegenlesen 11.10. 21:29 — **verifiziert (Logik)**
- p32a (a): je Quadrat skalarer Krawczyk für g(τ) = Φ(ζ̃(z)+τ) − sexp_K(z), Φ = P(e^{L·}): |g(0)| ≤ ε_U, |g'| ≥ m₁ auf D(ζ_c, s₁) (ζ̃(Z_b) ⊂ D(ζ_c, s₁)
  per Mittelwertsatz), |g''| ≤ M₂ auf D(ζ_c, s₁ + r_τ) ✓; Eindeutigkeit in D(θ̃, r_τ) ⇒ Überlappungen und Periodizität (Φ(z+1+θ) = e^{Φ(z+θ)},
  Quadrate decken Re ∈ [−0.52, 0.52]) ✓. (b) W = P⁻¹(sexp_K) per Rouché eindeutig, W(z+1) = LW(z) (|LW| < 0.1 geprüft) ✓, E = We^{−Lz} periodisch
  mit |E| ≤ C|q|^{−0.213} ⇒ hebbar ✓, Maximumprinzip in q ✓, Zweig von log E passt bei Im = 3 zu (a) (Abweichung 4e-167 ≪ r_τ = 1e-120) ✓.
- rest_schranke: |w| ≤ |x − L|/(1 − κ₂) ✓; |q| = |w|^{Re γ}e^{−Im γ·arg w}e^{2π Im θ} ⇒ C_q mit der richtigen Seite von Im u ✓; Schwarz für G_k
  (G_k(0) = 0, k ≥ 2) ✓; |H_k| ≤ √L(M_k/q_t)C_q|w|^{Re γ+1−k}/k! am Rand w_t maximal für k ≤ 5 (Exponent ≥ 0.447) ✓, H_6 separat
  (|w|^{−0.553}, c₆ = |iP/P'|⁶·6!·|P'(H)|) ✓, Integral mit |w| ≥ c_w|τ| ✓. Kleine Unschärfe: der Bereich von Im u im Rest wird aus der letzten
  Zelle ±0.01 genommen (Variation dort O(δ_min + κ₂) ≪ 0.01) — als Abschätzung begründet, nicht gerechnet.
⇒ Die von P stammenden Glieder der (K2)-Kette (P-V, p32a, rest_schranke) sind von V geprüft; v06/v07 hat K gegengelesen (BRETT 21:57) ⇒ (K2) **V**.

### 1.8 Ms 𝒯-Weg ((a″)_∞) — Struktur gelesen 11.10. 21:33, Lemmata 22:31 — **verifiziert: (a″)_∞ Kneser σ_min ≥ 0.641**
- Lemma M5(ii): 𝒯x = x + (x∘f)/m = J_∞x/m auf Γin; auf Γout 𝒯x = x + f'·(x∘h) = (J_∞x)∘h wegen f∘h = id, m∘h = f'∘f∘h = f' ✓.
- IMS: Σχ_j² = 1 ⇒ K − Σχ_jKχ_j = ½Σ[χ_j,[χ_j,K]], Kern ½Σ(χ_j(t) − χ_j(t'))²K(t,t') ✓; Re⟨(I+K)Y,Y⟩ ≥ (min_j c_j − ‖R‖)‖Y‖² ✓.
- Korollar M12: ‖𝒯x‖² ≤ ‖J_∞x‖²_{Γin}/m_in² + C_h²‖J_∞x‖²_{H²(D)} (Carleson-Schranke für g∘h auf Γout, h(Γout) ⊂ D) ⇒
  ‖J_∞x‖_{H²} ≥ c‖x‖/C_E ✓ — das ist (a″)_∞ für den VOLLEN Operator (ohne Abschneiden) ✓.
- m26 setzt c = min(c_mod − ε_L, c_m) − ‖R‖ und σ_min ≥ c/C_E' korrekt zusammen ✓.
- **Teilprüfung 22:07:** mk_ecke_kneser (Eckschranken für M8-K): gG aus log[(F−L)/(μδ)] − log[(F−L̄)/(z−L̄)] mit |R| ≤ S2|δ|²/2 längs der Sehne ✓,
  |ξ(t) − L| ≤ 2 sinθ·q/(1−q) für eᵗ ≤ q ✓, gWin aus 1/m = F'(z)e^{−ρz} ✓, h(y) als Fixpunkt von Φ(v) = L + (y − L − R(v))/μ (Selbstabbildung und
  Kontraktion nachgerechnet) ✓, gWout = S2·s richtig, weil W_out = f' am Randpunkt selbst ausgewertet wird (M5) ✓. In m22 die Gewichts-/Verschiebungs-
  faktoren der Schur-Majorante nur auf Plausibilität gelesen. M7: Ms c_mod = 0.98988322 stimmt mit Os λ_L (Lemma O10, von P per mpmath bestätigt)
  überein — zwei verschiedene Herleitungen, gleicher Wert ⇒ M7 numerisch doppelt gestützt.
- **M10-K (m24, Mitte) Logik geprüft 22:07:** Lagrange-Basis ℓ_i/√ω_i auf Gauss-Legendre-Paneelen orthonormal (Gauss exakt bis Grad 2n−1) ⇒
  Matrix √ω_i k(t_i,t'_j)√ω_j und λ_min(I + Re K_N) = min(1, λ_min(I + Re A)) ✓; verifizierte Cholesky (Higham 10.3): λ_min(B) ≥ s − γ‖R̂‖_F² ✓,
  ‖B − B_wahr‖ ≤ ‖Brad‖_F ✓; ‖K − K_N‖ ≤ ‖k − k_N‖_HS mit Tensor-Interpolationsfehler E_t + Λ_nE_{t'} (Hermite auf Bernstein-Ellipsen, Analytizität
  per Kästchen geprüft) ✓. Nicht nachgerechnet: die Majoranten M_t, M_t', M_fern im Einzelnen.
- **Gegenlesen der Lemmata M7–M10, M13 (Kneser-Fassungen m16, m22–m25, mk_kneser_arb) 11.10. 22:31 — verifiziert (Logik):**
  - **M7 (m16):** Symbol F0 = sin(s(π−y))/sin(sπ), F1 = sin(sy)/sin(sπ) aus ∫e^{aα}/(cosh α − cos y)dα = 2π sin(a(π−y))/(sin y sin aπ) ✓, Verschiebung
    |μ|^{±s} ✓, λ_min der 2×2-Hermite-Matrix ✓, Schwanz ‖A‖_F ≤ 2√2e^{−|σ|b}(|μ|^{−3}+|μ|³)^{1/2} ✓ (Faktor (1−e^{−2|σ|π})⁻¹ fehlt formal, Wirkung 1e-50).
    Selbst gerechnet: 1175 Intervalle, **0 gescheitert**, c_mod ≥ 0.9898000672 ✓.
  - **M8-K (m22):** Majorante J(·) nachgeleitet: cosh(α−v) ≥ e^{−|v|}cosh α ⇒ Faktor e^{du} und c_eff = e^{du}c ✓, |∇P0| = 1/(2π(cosh α − cos y))
    exakt ✓ (Zähler = (cosh α − cos y)²), Gewichtsabweichung √(w(t)/w(t'))/e^{−Δ/2} ∈ [(1−q1)/(1+q1), (1+q1)/(1−q1)] ⇒ dr = 2q1/(1−q1) ✓,
    Zerlegung K̃ − A = [√r W − e^{−Δ/2}/μ]P(Φ) + e^{−Δ/2}μ⁻¹[P(Φ) − P(Φ_mod)] ✓, Schur ✓.
  - **M9 (m23):** S ≤ min(1, (π²/8W²)(t−t')²) über Bogenlänge der Kurve χ auf der Einheitssphäre ✓, Riemann-Obersummen mit Kugel-Einschlüssen
    über ganze Kästchen ✓, Raster deckt [−Tb, Tb] (2Tb/d ganzzahlig in m26) ✓, Strip-Poissonkerne ✓, Schwanz: C_w = √(2/(1−cos θ)) nachgerechnet ✓,
    P ≤ (4e/2π)e^{−|Δ|} ✓. Kommentar „Schwanzpaare nur, wenn der andere Punkt in [−T2, T2]“ ist falsch (gegenüberliegende flache Bereiche: S = 1),
    die Schranke integriert aber über alle |Δ| ≥ Tb − T2 und deckt diese Paare ✓ — ohne Wirkung.
  - **M10 (m24):** Hermite-Rest mit C_int = a_e·max|ω|/(min_E|ω|·(a_e − 1)) ✓ (Umfang ≤ 2πa_e; Abstand [−1,1]→∂E = a_e − 1, weil Randabstand
    konkav und Brennpunktabstand a − c), Tensorfehler C_int(M_t' + Λ_n M_t) ✓, analytische Fortsetzung P0 = (1/2πi)[φ(Φ−t') − φ(Φ*−t')] mit
    φ = 1/(1−e^z) nachgerechnet ✓ (P1 mit φ1 = 1/(1+e^z) ✓), Pole ⇔ Im Φ ∈ {0, π} ✓, Fernpaare |cosh(α+iβ)| ≥ sinh|α| ✓, φ-Schranken ✓,
    Ellipsenkästchen decken ∂E ✓. Eigener Lauf mit Vs (K2)-Datei: Analytizität „geprüft (0 kritische Kästchen)“, Halbachse 0.305 < min(y, π−y) ✓.
  - **M13 (m25):** κ aus σ²e^{−2σy} ≤ (1/(e·y))² bzw. (1/(e(π−y)))² ✓, |Z(s+iy)|² ≤ (1+ε)|Z(s+iy1)|² + (1+1/ε)(y2−y1)∫|∂Z|² ✓, Young-Schritt
    x^θ ≤ (1−θ)s^θ + θs^{θ−1}x ✓, Massfaktor r = |ξ'(t+iπ)|/|ξ'(H)| ✓, m_in = min e^{ρ cos φ}/|F'| ✓. Eigener Lauf: Im H ∈ [2.456, 2.537] ⊂ (0, π) ✓,
    C_E' ≤ 1.42810 ✓.
  - **mk_kneser_arb:** Horner mit exakten Mitten + Koeffizientenradien + Rest aus (K2) (∫-Schranke, p − j − 1 > 0) ✓, Mittelwertsatz auf Kugeln/Bögen ✓,
    komplexer Krawczyk (konvexe Hülle von F'(X) genügt) ✓, Fortsetzung F(z) = log(ρF(w))/ρ aus f∘exp = exp∘f ✓. Formal: Rundung von k·mid(b_k)
    (120 Bit) fehlt in Rc1 (< 1e-30, ohne Wirkung).
  - **Schwachstelle im Zusammensetzer m26:** liest nur die letzte Zahl je Teillauf, prüft weder „gescheitert: 0“ (m16) noch „Analytizität geprüft“ (m24)
    noch die Höhen in (0, π) (m25). Für die Läufe mit Vs Datei alle drei hier einzeln nachgesehen ✓. Empfehlung an M: Abbruch bei fehlender Prüfung.
  ⇒ **(a″)_∞ für Kneser: σ_min(J_∞|P_N) ≥ 0.641 ∀N — V** (Fenster 7/23; bedingt nur auf (K1) P-V und (K2), beide V).
- **Ergebnis mit der zertifizierten (K2)-Datei (aicgolling, 4 Threads):** Fenster 7/23: Re⟨𝒯x,x⟩ ≥ 0.9153‖x‖², **σ_min(J_∞|P_N) ≥ 0.641 ∀N**;
  Fenster 6/14: ≥ 0.8134, σ_min ≥ 0.5696 (Logs ~/v_pruefstand/mrun/x/y/m_kopie/v10_*, v11_*). Bedingt nur noch auf die Prüfung der Lemmata oben.

### 1.9 Ks §18.1/§18.2 (explizites Grenzprofil des adjungierten Kopfes) — Gegenlesen und Gegenrechnen 12.10. 03:00 — **richtig (Skizze wie angegeben)**
- Faktorisierung μ(1 + μ^{s−2}) = K_LK_R: über die Spiegelungsformel Γ(½ − z)Γ(½ + z) = π/cos πz, z = iα(s−2), wird K_LK_R = 2μ^{s/2}cosh(πα(s−2)) und
  μ + μ^{s−1} = 2μ^{s/2}cosh((s−2)log μ/2) — gleich wegen α = log μ/(2π) ✓; numerisch 2e-31 ✓. α = L/(4π) ✓ (log μ = L/2, da log L = L).
- Nullstellen s_m = 2 + (2m+1)·2πi/L: Re s_m = 2 + 4.447(2m+1) ⇒ m ≥ 0 rechts, m ≤ −1 links ✓; s₀ = 6.44695 + 1.05794i ✓ (Exponent im Keil).
- Mellin-Schritt: μ⁻¹y(u/μ) ↔ μ^{s−1}Φ_y ✓; Wiener–Hopf: Ξ/K_L rechts, K_RΨ links analytisch ⇒ ganzes E ✓. **Annahme, kein Satz:** E = const·λ^{−s}
  („kleinste Wahl“) — ohne Wachstumsschranke für E ist die Lösung nicht eindeutig; das ist Ks angekündigter nächster Schritt (Satzform).
- Ordnung a = 1 + iα: |a| = 0.894 ✓ (das bekannte n^{−0.894}).
- Gegenrechnung mit eigenem mpmath-Code: Taylor-Koeffizienten von Y(λζ) für k = 1…5 und q_2…q_6 stimmen mit Ks Tabellen auf alle angegebenen Stellen;
  k28 und k29 neu ausgeführt: Ausgaben identisch mit k28.log/k29.log (Wiener-Abstand je nach Extrapolation 2.5e-7 … 2.7e-5 — K nennt den besten Wert;
  die Streuung über die vier Varianten ist die ehrliche Unsicherheit).

### 1.10 Ks §19.1–19.3 (Lemma K10.1, Satz K10, Satz K11/K11′) — Gegenlesen 12.10. 08:27 — **richtig (mit Ks eigenen Vorbehalten)**
- **K10.1:** Ergänzungssatz mit z = β_A + iαs = ½ + iα(s − 2): sin πz = cosh((s−2)log μ/2) ⇒ Γ(z)Γ(1−z) = 2πμ^{(s−2)/2}/(1 + μ^{s−2}) ✓; primal mit
  z = ½ + iα(s + 1) ⇒ μ(1 + μ^{−s−1}) = 2πμ^{(1−s)/2}/(ΓΓ) ✓. Pole von Γ(β_A + iαs): s = 2 + (2m+1)·2πi/L ✓ (rechts für m ≥ 0).
- **K10 (ii):** Taylor-Koeffizient von E_AF bei u^k = μ(1 + μ^{−k−2})c_k; mit K10.1 bei s = −k folgt Ks Formel (1 − β_A + iαk = ½ + iα(k+2)) ✓.
  (i) Ordnung 1/Re a = 1/0.8936 = 1.119 ✓. (iii) Barnes-Asymptotik zitiert, nicht nachgeprüft — wie K schreibt.
- **K11:** (1) Mellin längs Spiralstrahlen konvergiert genau für 0 < Re(as) < δRe a ✓ (|u| = e^{(Re a)t − (Im a)θ}). Dilatation: −log μ = a·2πi(1 − κ₁)
  (wegen 2πi(a − 1) = −2πα = −log μ) ⇒ Winkelversatz 2π(Re κ₁ − 1) = 42.55° ✓, die Verschiebung in t ist harmlos. (2) linke Pole ⊂ −ℕ₀ ✓. (3) rechts
  nur die Nullstellen von 1 + μ^{s−2} mit Re(as) > 0, das sind genau m ≥ 0 (Re(as₀) = 5.73, Re(as₋₁) = −2.16) ✓; dafür braucht man E_Ay = O(|u|^{−M}) im
  Sektor und Analytizität in 0 ✓. (4) E ganz ✓. (5) Liouville unter (E2) ✓, E(0) = 1/Γ(β_A) ✓, Π(s) ↔ Π(−u d/du) ✓. K11′ analog ✓.
- **Offen (Ks eigene Liste):** Barnes-Asymptotik (Braaksma), Lemma K12 (wann folgt (E2) aus dem Verhalten von y), und dass die echten Köpfe die
  Voraussetzungen erfüllen (nur numerisch). Bemerkung (a) (Nicht-Eindeutigkeit ohne (E2)) ist richtig und wichtig.

### 1.11 Ks Satz K19 (§30.2), Hypothesenliste §32 und HAUPTSATZ_KR0 Fassung 2 (§1–§5) — Gegenlesen 12.10. 11:26 — **richtig, mit vier Anmerkungen**
**K19.** (1) Identität: für q ∈ P_N ist T_NJ_∞q = J_Nq, also Π_Nq = q und Π_Nf − q = J_N⁻¹T_NJ(f − q) — ✓ (exakt, linear). (2) K13 mit fester Kontur 𝒞 durch L:
⟨ℓ_n, T_NW̃⟩ = (1/2πi)∮Λ_n(t)W(t)dt/t gilt für jede rektifizierbare Jordankurve um 0, auf deren Abschluss des Innern W stetig und im Innern holomorph ist (Cauchy für
Jordangebiete). Die Spiralen haben endliche Länge (ds = √(1 + T²)dr) ✓; W = (f − q)∘f + m(f − q) ist im Innern holomorph, wenn f das Innere nach E abbildet ✓ (Annahme
an E); in L, L̄ reicht Stetigkeit von f_K (Hölder 5.447) ✓. (3) Zweige bei sehr kleinem r: das Innere von 𝒞 bei L ist der Spiralsektor Σ₁ → Σ₂ über die D-Seite; er dreht
mit 1.623°/e-Faktor und überstreicht für r → 0 auch die Richtung von S_L — f_K ist dort die Fortsetzung über den Spiralschnitt. Nachgerechnet mit Ps oberer Struktur
(v44): pro e-Faktor in r ändert sich Im z um −0.0048 (Drehung) gegenüber +0.708 (|W|), die singulären Glieder bleiben O(|t − L|^{5.447}) mit beschränktem Faktor ✓ (Ks
„Faktor ≤ n^{0.027}“ aus §22.1 passt). (4) Abschätzung |W| ≤ (1 + M)dist_E auf 𝒞 und |q − f| ≤ dist_E auf K_{r₀} ⊂ E ✓; daraus die Formel ✓.
**Anmerkung A (E ⊃ D̄ falsch, E ⊃ D̄ ∖ U_{r₀} richtig)** — Ks Korrektur 11:20 bestätigt: mit D̄ ⊂ E schiebt sich die feste Γout-Tangente (346.6°) in die drehende Lücke.
**Anmerkung B (Aussenwinkel):** nahe L (r < r₁) ist E der Spiralsektor Σ₁ → f(Σ₂), Aussenwinkel u_A1 − u_A2 − 38.3° = 131.7° für (80°, −90°), nicht 99° (die „25°“-Grenze war
die gerade Kante von f(D̄), die E nicht enthalten muss). Jenseits r₁ verengt der Verbindungsbogen zur Linse (Linsenpunkt bei 111.6° von L aus) die Lücke auf ≈ 97°.
**Anmerkung C (Budget):** in §32/HAUPTSATZ ist (P) neben (G-P) überflüssig, wenn (G-P) als Schranke für |ΔK_n(z, L)| = |L·Δu_n(z)| statt relativ zu |Δu_n(L)| formuliert
wird (dann braucht man keine Unterschranke für |P_n|); mit der relativen Form ist eine Unterschranke |P_n| ≥ cn^{−b} stillschweigend nötig, sobald man R_P aus
Δu_n gewinnt. Gemessen ohnehin |P_n| = 0.535·n^{−0.106} (beidseitig).
**Anmerkung D (FA):** gebraucht wird die gleichmässige Form lim sup_n max_t (wie formuliert) ✓; dann |Λ_n/Λ_n(L)| ≤ e^{−ηn/2} für n ≥ n₀ auf 𝒞 ∖ U_{r₁} ✓.
**§32-Logik:** (EP) + (FA) + (S_E) ⇒ Σ_{n≤N}Q_n ≤ CN^γ ⇒ sup_{K_{r₀}}|Π_Nf − f| ≤ CN^{γ−s_E} ✓ (Rang 1: Q_n = |ΔK_n(z, L)|·max_t|Λ_n(t)/Λ_n(L)|).
**HAUPTSATZ Fassung 2:** §1 enthält r₀ ∈ (r_*, 0.1] mit r_* ~ 1e-6 (Kneser, Schätzung K) und ≈ 0.13 (M4) — der Auftrag der Leitung (r₀ > r_* ausdrücklich) ist damit erfüllt.
Prüfstellen §5: 1 ✓ (oben); 2 mit Anmerkungen A, B; 3 s_E gemessen (§3.19); 4 Kompaktheit ✓ (P_N endlichdimensional, sup_{K_{r₀}} ist eine Norm, 𝒮_N als
Nullstellenmenge abgeschlossen in jeder Norm), H2 ohne sup_{D̄} noch nicht im Wortlaut von ARCH §4.3 geprüft; 5 F3 gemessen (§3.20).
Fassung 1, Prüfstelle 1 (Faktor 2 in der Schlitzkonstante): **die „2·“ in ARCH §9.7 (1) ist ein Schreibfehler** — der Wert 12.7 ist (1/π)(C_J/c_d)αΓ(5.447α)c_g^{−5.447α}
und enthält beide Schlitze schon (2·1/(2π) = 1/π); nachgerechnet 12.75 (W̄). Für K′ also **8.2** (C_J = 1.0378, c_g′ = 0.754, c_d = 0.933, α′ = 0.787), nicht 16.3. Der weitere
Faktor 2 in „C′ ≈ 2·12.7“ ist der Familienfaktor (X4a), keine Doppelzählung. In v31/v41 steht dieselbe Struktur (pre = 2ρ/2π für zwei Schlitze, Faktor 2 Familie) ✓.

### 1.12 Ks restliche Sätze aus HAUPTSATZ_KR0 §7 — Gegenlesen 12.10. 11:58 — **alle richtig; Anmerkungen bei K2, K5, K20**
- **Lemma K2** (§0 Nr. 3): T_N^wC_kF = C_k^{(N)}T_NF für k(0) = 0 (C_k untere Dreiecksmatrix, Diagonale k′(0)^m ≠ 0) ✓; daher Ĵ_N⁻¹T_N^wC_kJ_∞ = J_N⁻¹T_NJ_∞ = Π_N ✓.
  *Anmerkung:* „wirkt auf Grad-n-Anteile wie das Phasenraumgewicht e^{nε Re φ}“ ist eine WKB-Lesart, kein Teil des Lemmas.
- **Satz K3** (§9.2): bei q < max a_i maximieren ≥ 2 A-Terme, bei q = max a_i wegen Σ = ∅ ebenfalls ≥ 2 A-Terme, q > max a_i widerspricht der Bilanz; Verschiebung um δ erhält
  das ✓.
- **Satz K4** (§9.3): Ableiten von P∘P = eˣ + R_P ✓; für Fast-Lösungen R_P ≡ cx^N mod x^{N+1} ⇒ R_P′ ≡ cNx^{N−1} mod x^N, der Term T_N[M_{R_P′/P′}·] hat nur in den Zeilen N−1, N
  Einträge, und zwar aus δ₀, δ₁ — Rang ≤ 2 ✓. (Gilt für die Tangente J(P); die Sekanten 𝒥_P haben den Überlauf nicht, wie K schreibt.)
- **Satz K5** (§11.3): 1, 3, 4 ✓ (in v48/v51 benutzt, Ks k46/k49-Werte reproduziert); 2: C_kL ist untere Dreiecksmatrix mit Diagonale k′(0)^m, also LU mit nicht normiertem L —
  „dasselbe U“ stimmt; 5: (U⁻¹)_N = (U_N)⁻¹ für obere Dreiecksmatrizen, Kompression ⇒ ‖Ĵ_N⁻¹‖ ≤ ‖U⁻¹W⁻¹‖‖WL̂⁻¹‖ ✓ (Grundlage meiner VE4-Rechnung).
- **Lemma K6** (§11.2): |ρ/t|^k = e^{−ks}, Dreiecksungleichung ✓; Folgerung mit (N+1)(e^{εN} + e^{(c′−ξ₀s)N}) ✓.
- **Lemma K13** (§22.1): Λ_n Laurent-Polynom in 1/t vom Grad n ⇒ (1/2πi)∮Λ_nW dt/t = Σ_{k≤n}ℓ_{n,k}ρ^kW_k auf jeder Kontur um 0 in einem einfach zusammenhängenden Gebiet,
  in dem der Zweig von W mit dem Keim bei 0 analytisch ist ✓; kopfweise Konturen zulässig, da I_n nur Koeffizienten ≤ n sieht ✓.
- **Rang-1-Zuwachs, P_n = K_n − K_{n−1}** (§24.1) ✓ (aus K5; in v48 auf alle Stellen reproduziert); **sup_{D̄}|δ_n| ≤ √(n+1)‖c_n‖** (§23.1) ✓ (Cauchy–Schwarz auf |z| ≤ ρ).
- **Reduktionssatz, lokale Fassung** (§23.2, gegen kn.satz.polynom_halbexp_reduktion gelesen): Schritt 1 (Montel) braucht nur lokale Schranken; Schritt 2 braucht |p_N| ≤ B′ auf
  einer festen Scheibe |y| ≤ ρ₀′ < |L| (dort Cauchy für p_N′ und |p_N∘p_N| ≤ B′ für |x| ≤ s) ✓; Schritte 3, 4 unverändert ✓.
- **Satz K16** (§27.1): Induktion nachgerechnet — Bn₀/n ≤ Bn₀^b|p_n|/c_p folgt aus (n₀/n)^{1−b} ≤ 1 und |p_n|n^b ≥ c_p; ¾|P_n| ≤ (3/2)|p_n| ⇒ |P_n| ≤ 2|p_n|; beide Schranken ✓.
- **Integraldarstellung F_β, Watson-Asymptotik** (§27.3): unabhängig geprüft — V(ζ) = Γ(β)⁻¹∫t^{β−1}e^{−t}(1 − e^{−w})/w dt (w = λ̃ζt^{−iα}) = Reihe auf 10 Stellen; die Mellin-Residuen
  von Γ(s)Γ(β + iαs) sitzen bei s′_m = γ(β + m), γ = i/α — dieselben Exponenten wie Ks Watson-Reihe, mit derselben Konstante (γ/Γ(β))(−1)^mΓ(s′_m)/m!; Asymptotik =
  Integral auf 1e-12 im Potenzsektor (v54, §3.21) ✓.
- **Spiralfenster** (§28.2): arg u = c + T ln|u| mit T = Im a/Re a (Rand {Re(C u^{1/a}) = 0}); mit u = λ(t − L)/d_n wird θ = c″ + T ln r + (T Re a − Im a) ln n, Klammer = 0 ✓.
- **Satz K19, Hypothesenliste §32** ✓ (§1.11).
- **Lemma K20** (§31.1): richtig für die oberhalbstetige Regularisierung u* (subharmonisch); für den punktweisen limsup gilt u ≥ 0 auf ∂U nur bis auf eine polare Menge (Brelot–Cartan) —
  für die Aussage „keine exponentielle Reserve auf K_{r₀}“ genügt das. Die Voraussetzung u = 0 im Innern ist numerisch (Innengesetz n^{−3.05}).
- **O1** (§33.2): (1/n)ln|Λ_n| = (1/n)ln|Polynom_n(t)| − ln|t| hat Laplace-Masse 2π(ν_n − δ₀), Gesamtmasse 0 auf Ĉ; superharmonisch auf Ĉ ∖ {0} mit −ln|t| in 0 gäbe Masse ≤ −2π ✓.
- **Satz K22, Glattheit von f_K auf E und Verlust δ_T** (§35.2 Glied 3): auf dem Spiralsektor ist arg W = θ₀ + T ln|W| + O(1), also |W^{s₀}| = |W|^{5.447}e^{−1.058 arg W} ≍ |W|^{5.447−1.058T}
  = |W|^{5.417} ✓ — dieselbe Zahl wie meine unabhängige Rechnung in v44 (Im z verliert 0.0048 je e-Faktor, 2π·0.0048 = 0.030). f_K^{(5)} ~ W^{s₀−5} ⇒ Hölder 0.417 ✓; APV Satz 2
  (von K gelesen, nicht von mir) mit r = 5, k = 1 ✓; Geometrie E ∈ H* für f_K numerisch (§3.23) ✓; Rate gemessen (§3.19) ✓.

### 1.13 Ks Christoffel–Darboux-Identität (CD1) (k_kern/BERICHT §40.1) — Gegenlesen + eigener Code 12.10. 13:35 — **richtig (bewiesen), numerisch bestätigt, kubisch und Kneser**
**Herleitung Schritt für Schritt nachgerechnet:**
- Spalten von J_n: p_j = T_n[f^j + m x^j]. (JS)e_j = T_np_{j+1} (j < n), (SJ)e_j = T_n(x·T_np_j) = T_n(x p_j) (Abschneiden vertauscht mit x·, weil x·(Grad > n) Grad > n+1 hat).
  p_{j+1} − x p_j = (f − x)f^j ✓; für j = n: −T_n(x p_n) = T_n((f − x)f^n) − T_nf^{n+1} (m·x^{n+1} fällt weg) ✓ ⇒ E = T_nM_{f−x}C_f − (T_nf^{n+1})e_nᵀ ✓.
- SX − XS = X(JS − SJ)X ✓ (nur J_n invertierbar nötig, also (F1)). e(z)ᵀS = z e(z)ᵀ − z^{n+1}e_nᵀ ✓, Sê(t) = tê(t) − t e_0 ✓.
- Einsetzen: (z − t)K_n − z^{n+1}λ_n(t) + tν_n(z) = e(z)ᵀX T_nM_{f−x}C_f Xê(t) − (J_n⁻¹T_nf^{n+1})(z)λ_n(t) ⇒ (CD1) ✓.
  C_f Xê(t) = K_n(f(x), t) als Reihe ⇒ Ξ_n = J_n⁻¹T_n[(f − x)K_n(f(·), t)] ✓; T_nJq_{n+1} = T_nf^{n+1} − J_nJ_n⁻¹T_nf^{n+1} = 0 ✓, q_{n+1} monisch vom Grad n+1 ✓.
- Konturform: J_n⁻¹T_ng(z) = (1/2πi)∮K_n(z,s)g(s)ds/s ✓, **aber nur auf Kreisen im Konvergenzkreis von g**, für f_K also |s| < ρ (nicht „jeder Kreis“; beim Kubik jeder).
- Normierung (kosmetisch): in Ks Konvention (c_nn = ℓ_nn) ist q_{n+1}(z) = ρ^{n+1}·δ_{n+1}(z)/ℓ_{n+1,n+1} (Faktor ρ^{n+1} fehlt in §40.1); λ_n = ρ^{−n}Λ_n ✓.
**Eigener Code** (`v60_cd1.py`, Logs `cd1/`): ungeskalierter Rahmen, J_n direkt aus den Spalten, J_n⁻¹ einmal per Mittelpunkt-Lösung, Ξ ZWEIMAL (Matrixform und
Konturform mit Trapezregel auf |s| = 0.6ρ, 6n + 64 Punkte), Wiederholung mit +256 Bit. Fünf Punktpaare: beide in D; Steg ausserhalb D (z = 1.03ρ∠175°, t = −1.03ρ);
beide ausserhalb (1.1ρ∠115°, 1.05ρ∠−143°); Diagonale z = t = −ρ; Γout-Seite (z in D, t = 1.2ρ∠17°). n = 60/150/300, kubisch und Kneser:
- |LHS − RHS|/(grösster Term) ≤ 2e-89 überall (Arbeitsgenauigkeit 390–750 Bit), Ξ Matrixform/Konturform − 1 ≤ 3e-79. **(CD1) gilt, auch ausserhalb D, für beide Modelle.**
- Ks Tabellenwerte (kubisch, t = −ρ) unabhängig reproduziert: m·ν_n(t) = 0.470077 (n = 60/150/300), B(t) = 1.03970/1.03968/1.03968, m·Ξ_n(t,t)/t = −0.56962/−0.56960/−0.56960 ✓.
- Kneser wie kubisch: am Steg |K_n| = 0.898/1.662 (n = 150/300, 1.03ρ), Ξ(−ρ,−ρ) = 1.074 für alle n. Ohne den Faktor m(t) (f_K′∘f_K bei −ρ liegt auf dem
  Konvergenzrand, nicht ausgewertet): ν_n(−ρ) = 0.639470, q_{n+1}λ_n/t = 1.42049/1.42046, Ξ_n(t,t)/t = −0.78102/−0.78099 (n = 60/150) — dieselbe n-Unabhängigkeit wie kubisch
  (`cd1/v60_kneser_steg.log`).

### 1.14 HAUPTSATZ_KR0 Fassung 3, Prüfstellen §5 (teilweise) — 12.10. 14:05
- **Prüfstelle 4 ((KK′) mit eigenem Code):** bestätigt bis n = 2000, kubisch und Kneser (§3.26 (a)): h_P + h_A = +0.00025 → +0.00006 (fallend wie 1/n), |ΔK_n(t,t)| → 1.304,
  |m·ΔK_n| → 0.956 = κ(−ρ), für jedes S ∈ [0.0126, 0.03] gleich.
- **Prüfstelle 2 ((F3′) über einen Keim-Repräsentanten) — richtig mit einer Zusatzbedingung.** Die Ersetzung ist zulässig: T_N𝒥_P hängt nur vom Keim ab, und K13 braucht
  nur, dass W̃ = (f − q_N)∘P̃ + M_{P̃}(f − q_N) in int 𝒞 analytisch ist (dann ist ∮ΔK_nW̃ = Kopfanteil von T_NW̃ = T_NW). Für M_{P̃} (geteilte Differenz von f) genügt,
  dass die Strecken [f(t), P̃(t)] im Holomorphiegebiet von f liegen. **Aber:** Schritt 3 von Satz A (|W| ≤ (1 + M)·sup_E|f − q_N| auf 𝒞) benutzt f(𝒞) ⊂ E. Für P̃ braucht man
  P̃(𝒞) ⊂ E bzw. sup_{P̃(𝒞)}|f − q_N|. f(𝒞) liegt auf dem Rand von E (E = abgeschlossene Hülle von int 𝒞 ∪ f(int 𝒞)), also verlässt P̃(𝒞) E schon bei kleiner Abweichung,
  und nahe L ist eine feste δ₀-Umgebung von E nicht mehr im Holomorphiegebiet von f_K (Spiralschnitt). (F3′) muss deshalb lauten: P̃(𝒞) ⊂ E′ mit einem festen Kompaktum E′ ⊃ E,
  auf dem (S_E) gilt (oder |P̃ − f| ≤ δ(t) mit δ(t) → 0 bei L, L̄). Ausserdem ist das Supremum über GANZ 𝒞 zu nehmen, nicht nur über 𝒞 ∖ K_{r₀}: Auf K_{r₀} ist P̃ ≠ P.
  (Mein CF misst sup über ganz Int γ ✓.)

## 2. Endliche Bausteine
### 2.1 Zweig-Zertifikat 193 → 400
Skripte: `zert400/zert_kern.py`, `zert_job.py`, `zert_kette.py` (unverändert von Strang X übernommen, Daten-Prüfsummen = Originale in gemeinsam/
und p_pruefung/p24_zustand.json), `lauf.sh` (eigener Pfad, halbexp-venv). Start 11.10. 19:55 (Server 17:55 UTC).

**Ergebnis (11.10. 21:59): Kette 193 → 400 vollständig, alle 207 Übergänge ok** (zert_kette.py: „KETTE 193 → 400: ok = True; Fehler: []“).
Unabhängig nachgeprüft mit eigenem Code `zert400/v16_kette_pruefen.py` (Ps Box p_193 ⊂ B_193; für jeden Übergang ok/streng/start_ok/innen/
in_roehre; K_{N+1} ⊂ B_{N+1} koordinatenweise strikt in arb, 1600 Bit): **keine Fehler**, kleinste Reserve 0.49986, grösster Wrapping-Term t3 = 1.4e-4,
grösstes θ = 3.4e-6. Die Strenge jedes Übergangs liegt in zert_kern.uebergang (§1.4 verifiziert). Daten: `zert400/lauf/` (boxen, erg, log, 59 MB).
Damit ist der σ-Homotopie-Zweig p₁ → p₄₀₀ rigoros (p18 + p24 bis 193 von P, 193 → 400 von V).

**Kopplung (A) bei N₀ ∈ [318, 400] (`v15_A_kopplung.py`, `v15_A.json`, rigoros):** sup_{D̄}|p_N − f_K| ≤ d3 + S₆₄₀ mit
d3 = ‖p_N − T₆₄₀f_K‖_{ℓ¹_ρ} (Boxen B_N bzw. K_400, Ps Koeffizienten mit Dezimalrundung) und S₆₄₀ ≤ 461.36·640^{−5}/5 = 8.6e-13 ((K2)):
| N | 318 | 330 | 350 | 370 | 387 | 400 |
|---|---|---|---|---|---|---|
| sup_{D̄}|p_N − f_K| ≤ | 3.58e-12 | 3.13e-12 | 2.58e-12 | 2.16e-12 | 1.90e-12 | 1.75e-12 |
Verlangt ist sup ≤ R_{N₀} (Xs Abschätzung R ~ 1e-3…1e-4 bei N₀ ≈ 300–400) ⇒ **Reserve ≈ 1e8**. (A) ist damit für jedes N₀ ∈ [193, 400] als endliche
Aussage erledigt; offen bleibt nur, dass (H) mit seinen Konstanten ab diesem N₀ gilt.

### 2.2 (K2)
**Wie scharf muss (K2) sein? (gemessen, 20:01–20:27, Kopien von Ms Skripten in m_kopie/, Daten = Ps zertifizierte a_0…a_640)**
| Lauf | Schwanz | Fenster T₁/T₂ | Ergebnis | Log |
|---|---|---|---|---|
| v01 | C = 5e5, p = 6 | 7/23 | m24 scheitert („Streifenzweig“) | v01_m26_C5e5_p6.log |
| v02 | C = 5e4, p = 6 | 7/23 | m24 scheitert | v02_m26_C5e4_p6.log |
| v03 | C = 100, p = 6.4 (Ms Annahme) | 7/23 | m24: c_m ≥ 0.987400 (= Ms Wert, Kopie reproduziert) | v03_m24_default.log |
| **v04** | **C = 5e5, p = 6** | **6/14** | **c ≥ 0.8134, σ_min(J_∞|P_N) ≥ 0.5695 ∀N** | v04_m26_C5e5_p6_T6_14.log |
Ursache des Scheiterns bei T₂ = 23: das Fensterende liegt bei |z − L| ≈ e^{−23}|L − L̄| ≈ 3e-10; der Schwanz R₀ = Σ_{k>640}|b_k| ≈ 9e-11 (C = 5e4)
ist dann vergleichbar mit |F − L|, und der Log-Zweig in xi_inv wird mehrdeutig. **Folgerung: Ziel für (K2) ist C₆ ≤ 5e5 mit p = 6**; damit
gilt Ms (a″)_∞ mit c₀ ≥ 0.5695 (besser als die 0.54 in Xs Bilanz). Ps Ziel 300 war unnötig scharf.
**Diagnose von p32 (warum P im Strukturgebiet abbrach, v05_dbg_struktur.py):** bei θ = 0.5 und Zellbreite h = 0.004 liefert
P_jet(H₀) für den Zellenball **inf** (Taylor-Modell-Radius R_k = 16·r_y ≈ 1.5 um |H₀| ≈ 1.87 ⇒ grobe Kreisschranke von P über eine
11-fache exp-Kette), und p32s Abbruchtest `hb·s7 > 0.2·c6 + 1` ist für nan FALSCH — die Zelle wurde mit nan angenommen (Integral = nan).
Mit h = 5e-4 / 1e-4: sup|∂⁷f| über die Zelle 1.1e6 / 1.8e5 (Punktwert 200), endlich. Mit der Reserve aus v04 genügt eine Zellbreite
h ≈ 8e-4 (Annahme, wenn hb·sup|∂⁷f| ≤ 1e3) — Gesamtüberschätzung ≲ Faktor 10, weit unter dem Bedarf.

**(K2) ZERTIFIZIERT (11.10. 21:13; rechnergestützt, arb):** ‖∂_θ⁶f_K(ρe^{iθ})‖_{L¹(0,2π)} ≤ **2898.77** ⇒ **|b_k| ≤ C₆k^{−6} für alle k ≥ 1,
C₆ = 461.36** (numerisch ≈ 124, Überschätzung ×3.7). Damit Σ_{k>640}|b_k| ≤ 8.6e-13, Σ k|b_k| ≤ 6.9e-10, Σ k²|b_k| ≤ 5.9e-7.
- Rechnung: `pv_kopie/v06_k2_rand.py` (repariertes p32: Variable d = |θ − θ_L| statt float-θ, nan/inf ⇒ Zelle halbieren, grobe Annahmeregel,
  Zweigtreue je Zellgrenze per Punkt-Krawczyk), 19 Abschnitte auf aicgolling (Logs `pv_kopie/lauf/log/`), Rest |θ − θ_L| < 1e-12 nach Ps
  rest_schranke (rechts 0.0021, links 0.057); Zusammensetzen `v07_k2_zusammen.py rechts_8 links_9` (Abschnittsgrenzen, Lückenlosigkeit,
  Zweigtreue an den Grenzen — alles ok): rechts ∫ ≤ 388.96, links ∫ ≤ 1060.36 (je obere Halbkreishälfte, Faktor 2 für die Spiegelung).
- Datei für M/O: **`v_pruefstand/kneser_halbexp_b_zert_V.json`** (Format {"K","b_mitte","b_radius","C","p"}, C = 461.359, p = 6; Radien mit
  Dezimalrundung der Ps-Mitten, aufgerundet) — `v09_K_datei.py`.
- Benutzt (und von V gelesen): Satz P-V (§1.1), Ps Strukturformel p32a (θ_K) und rest_schranke (**p32a und rest_schranke noch nicht
  gegengelesen** — Prüfschlange). Die 6-fache partielle Integration ist zulässig: ∂_θ⁵f stetig auch in θ_L (|θ−θ_L|^{0.447}-Verhalten),
  ∂_θ⁶f ∈ L¹, b_k = Fourier-Koeffizienten der Randwerte (f ∈ A(D)).
- Gegenprobe des Rests (alte Läufe ohne DMIN, rechts_8/links_9, nach 6660 s bei d ≈ 1.5e-18 beendet, 22:47): ∫ bis d = 1.8e-18 = 14.1105 ± 3e-5
  bzw. 180.620 ± 4e-4 gegen Neulauf ∫ bis 1e-12 = 14.11044 (+ Rest ≤ 0.00206) bzw. 180.6195 (+ Rest ≤ 0.0571) — der Beitrag von [1e-18, 1e-12]
  ist ≈ 6e-5 bzw. 5e-4, also weit unter der Restschranke ✓ (Restschranke konservativ).
- **Zweiter Leser K (BRETT 21:57): verifiziert**, zwei formale Punkte: (F1) Mittelpunkt-Scheibe ⊂ Zellscheibe prüfen, (F2) Summe in arb statt float.
  Beide umgesetzt (v06 Zeile „F1“, v07 arb-Summe mit Faktor 1 + 1e-12), alle 19 Abschnitte neu gerechnet (~/v_pruefstand/k2b, `pv_kopie/lauf_F1/`):
  **identisches Ergebnis C₆ ≤ 461.354, alle Prüfungen bestanden** (F1 hat nie angeschlagen). (K2) zählt damit als **V**.

## 3. Experimente für K

### 3.1 (S_γ) unabhängig nachgerechnet (Plan A2/A1-Vertiefung, Auftrag Leitung 21:20) — `v12_konform_S.py`, Logs `sg/`, Server ~/v_pruefstand/sg
**Aussage (wie K):** Ĵ_N = T_N^wC_kJ_∞|_{P_N} : ℓ²_ρ → ℓ², k: 𝔻 → Int γ (Γout ∪ Linse), Linse s(τ) = S(1 − e^{−(τ/τ₁)^p}), S = 0.0166, τ₁ = 12°, p = 1.3.
**Unabhängig:** eigenes Kubik, eigene Theodorsen-Iteration (M = 2^16) mit Prüfung der Randbedingung, Ĵ_N direkt als Spalten
T_N^w[(f∘k/ρ)^j + (m∘k)(wE)^j] (K bildet C_k·J̃_N), nur arb_poly mit expliziter Abschneidung (kein arb_series, keine ctx.cap-Falle), Präzision
1.2N + 400 Bit, keine float-Auswertung in D (F_n als exakte Zahlen; Randprüfung nur auf ∂𝔻).
**Fallen-Überlegung „Abbildungsfehler täuscht Stabilität vor“:** Ks K1″(a) braucht NUR eine analytische Abbildung k mit k(0) = 0, k'(0) ≠ 0 und
k(𝔻̄) im Analytizitätsgebiet von f, m — keine exakte Konformität. Ĵ_N hängt nur von F_0…F_N ab und wird für genau diese (exakten) Zahlen exakt
gerechnet; getestet wird also die Abbildung, die wirklich benutzt wird. Für Kneser zusätzlich nötig: L ∉ k(𝔻) (sonst hat f_K∘k einen
Verzweigungspunkt in 𝔻); die numerische Randkurve geht bis auf ~1e-5 an L heran — für die Theorie eine Abbildung mit Rand genau durch L nehmen.
Robustheit gegen Abbildungsfehler wird direkt getestet (gröbere Theodorsen-Auflösung, zufällig gestörte F_n, anderes Linsenprofil).
**Ergebnisse (Stand 21:45; Logs sg/*.log):**
| Modell, Linse | N | ‖Ĵ_N⁻¹‖ | Bemerkung |
|---|---|---|---|
| kubisch, Ks Linse S = 0.0166 | 100 / 200 / 400 / 800 / 1000 / 1200 / 1500 / **2000** | **1.251155** (alle gleich, 7 Stellen) | 2./3. Singulärwert 1.2022 / 1.1553 = Ks Werte |
| Kneser (Ps zertifizierte Mitten), S = 0.0166 | 300 / 400 / 500 / 600 / 639 | **1.250337** (alle gleich) | Datengrenze a_640 |
| Falle: Theodorsen M = 4096 statt 65536 | 400 / 800 | 1.251155 | Abbildungsgenauigkeit unerheblich |
| Falle: F_n zufällig gestört 1e-6 / 1e-4 | 400 / 800 | 1.251147 / 1.250381 | stetig, keine Verstärkung |
| Falle: +250 Bit | 800 | 1.251155 | Präzision ausreichend |
| kontrolle S = 0 (k = ρw, Ĵ = J̃) | 200 / 400 / 600 | 1.852 / 9.100 / 70.71 | = ‖J̃_N⁻¹‖ (Ks Werte) — direkte Konstruktion konsistent |
| S = 0.0003 / 0.001 / 0.003 | 200 / 400 / 600 | 1.76/8.14/59.6 · 1.56/6.28/40.1 · 1.33/3.04/12.9 | wächst |
| S = 0.006 | 200 / 400 / 600 | 1.3135 / 1.3135 / 2.4836 | Einsatz ≈ 550 |
| S = 0.008 | 600 / 800 / 1000 / 1200 | 1.3005 / 1.622 / 2.984 / 5.935 | **blinde Vorhersage (Einsatz ≈ 800, Rate ≈ 0.0033) bestätigt**: Rate 0.0031 / 0.0034 |
| S = 0.010 | 800 / 1100 / 1400 | 1.2882 / 1.2882 / 1.2927 | Einsatz um 1400 (vorhergesagt „> 1400“) |
| **anderes Profil („sin“: S·sin(πτ/τ_end)^{1.3}, gleiches S = 0.0166)** | 400 / 800 | 1.264 / **3.681** | **wächst** — Linse zwischen 5° und 40° zu dünn |
| Bezugslinse S = 0.05 (Ks Profilform) | 200 / 400 / 800 / 1500 / **2000** | **1.110071** (alle gleich) | Nachtrag 22:48 |
**Gesetz (gemessen):** ‖Ĵ_N⁻¹‖ ≈ max(B_S, A·e^{(c(N) − κS)N}) mit **κ = 1.00 ± 0.02** (Raten bei N ≈ 500: S = 0/0.0003/0.001/0.003 ⇒ 0.01025/0.00996/
0.00927/0.00724), c(N) = Eckrate von J̃_N⁻¹ (0.0103 bei N ≈ 500, 0.0114 bei 1100, → 0.0118). **Schwelle für Ks Profilform: S* = c_∞ ≈ 0.0118.**
Das „sin“-Profil zeigt, dass es nicht auf s(π) allein ankommt: die Linse muss überall über Ws Γin-Profil Λ(τ) (W §18: 0.0009/0.0025/0.0049/0.0069/
0.0091/0.0116 bei 2°/4°/8°/15°/30°/103°) liegen. Ks Linse liegt überall um den Faktor 1.4–1.7 darüber (bei τ → 0 wegen 1.3 < p_Λ ≈ 1.5 sogar
beliebig weit), das sin-Profil liegt bei 15° mit 0.0024 deutlich darunter.
**Bewertung:** (S_γ) für Ks Linse ist unabhängig **bestätigt** (numerisch, kubisch bis N = 1500, Kneser bis 639, robust gegen Abbildungs- und
Präzisionsfehler). Präzise Fassung, die die Daten tragen: **(S_γ) gilt für jede Linse ℬ, die mit Reserve über dem n-te-Wurzel-Profil Λ der
Eckmoden liegt, und versagt darunter.** **Einschätzung (Vermutung, kein Satz):** In der w-Koordinate sind die Exponentialraten der Eckmoden NEGATIV (c − s < 0). Ein Beweis von (S_γ)
könnte deshalb an der Ecke mit n-te-Wurzel-Genauigkeit plus Reserve auskommen (W11/W12-Niveau), sofern der Bulk gleichmässig kontrolliert wird
(Parametrix mit exponentiell kleinem statt polynomial genauem Fehler). Das wäre ein Weg um die Sackgasse „e^{o(N)} reicht an der Ecke nicht“
(W §26.3) herum — dort war die Rate 0, hier ist sie c − s < 0. Offen ist genau, wie man Bulk und Ecke trennt, ohne die alte Lücke wieder einzuführen.
**Korrektur nach Ks Einwand (BRETT 21:55, akzeptiert):** Die Reserve gilt nur für den Steg-Anteil fern von L. Die Linse ist bei L tangential
(s(0) = 0 = Λ(0)); in der Eckzone (Skala N^{−0.894}) sind alle Raten O(1/N) — dort bleibt polynomiale Genauigkeit nötig. Ausserdem enthält (S_γ)
die (H₀)-artige Aussage für alle glatten Eingaben g (‖J_N⁻¹T_Ng‖ ≤ ‖Ĵ_N⁻¹‖‖g∘k‖), deren Lösungen J_∞⁻¹g bei L die Singularität (x−L)^{s₀} tragen.
Die Vermutung „n-te Wurzel + Reserve genügt“ ist daher auf den Steg zu beschränken; die Ecke L selbst bleibt der Engpass, jetzt aber in einer
Form ohne exponentiell wachsende Moden.

### 3.2 Reserve-Fenster der Linse (Auftrag Leitung 21:50) — `v13_linsen_fenster.py`, `v14_lambda.py`, Logs `sg/v14_*`, `sg/v12_dick_*`
**(a) Untergrenze = Wachstumsprofil Λ(τ) der Eckmoden auf Γin** (G_N auf dem Kreis, ln G = a + κ ln N + ΛN, kubisch N = 600…1500;
Kneser bis 639 zu kurz für eine Anpassung). Eigene Werte (Λ_V) neben Ws Werten (Λ_W, N ≤ 1000):
| τ | 2° | 4° | 6° | 8° | 15° | 30° | 60° | 103.4° |
|---|---|---|---|---|---|---|---|---|
| Λ_W | 0.0009 | 0.0025 | — | 0.0049 | 0.0069 | 0.0091 | 0.0109 | 0.0116 |
| Λ_V (Anpassung) | 0.0006 | 0.0034 | 0.0060 | 0.0074 | 0.0087 | 0.0099 | 0.0115 | 0.0121 |
| Rohrate 1200→1500 | 0.0008 | 0.0022 | 0.0036 | 0.0046 | 0.0067 | 0.0090 | 0.0109 | 0.0116 |
Fern von L stimmen alle überein (Λ(π) = c_∞ ≈ 0.012). Bei 4°–15° steigen die Raten mit N noch deutlich (Anpassung mit grossem negativem κ);
der Grenzwert liegt zwischen Rohrate und Anpassung — **Unsicherheit bis Faktor 1.6 bei τ ≈ 6°–8°**. Konservativ: Λ_max = max(Λ_W, Λ_V).
Probe des Gesetzes „Rate = max_τ(Λ − s)“: Ks Profilform skaliert gibt κ = 1.00 (§3.1); das sin-Profil (bei 15° um 0.0045 unter Λ) wächst.
**(b) Obergrenze = K1-Geometrie** (für Kneser bindend; für das kubische f gibt es keine, siehe unten). Bedingung f(ℬ) ⊂ D (dann trifft f(ℬ) die
Schlitze nicht, K_γ = K' ∪ ℬ̄ behält die K'-Ecke, m = eˣ/f' ist auf ℬ analytisch): U(τ) = sup{s : |f(ρe^{s}e^{i(θ_L+τ)})| < ρ}, kubisch und Kneser gleich:
U = 0.0070 / 0.0139 / 0.0279 / 0.0565 / 0.1153 / 0.224 / ≥ 0.5 bei τ = 0.5° / 1° / 2° / 4° / 8° / 15° / ≥ 30°; bei L exakt U ≈ τ·tan(arg μ) = 0.790τ
(aus f(x) ≈ L + μ(x − L): |f| < ρ ⇔ s·Re μ < τ·Im μ). Weitere Grenzen weiter aussen: Holomorphie von f_K bis |x| < 3.16 (s < 0.83 bei τ ≈ 20°),
f_K' ≠ 0 auf ℬ (automatisch, da f_K'(x) = f_K'(eˣ)eˣ/f_K(eˣ) mit eˣ ∈ D für s < 4.2τ). Tangentialität bei L verlangt s = o(τ) — jedes Profil ~τ^p, p > 1.
**Operatorseite:** dicke Linsen sind sogar BESSER konditioniert (kubisch, N = 200/400/800 jeweils konstant): S = 0.05: 1.110; 0.1: 0.969; 0.2: 0.801;
0.3: 0.712 — S = 0.3 verletzt U bei 8° und bleibt stabil ⇒ die Obergrenze kommt allein aus Knesers Singularitäten, nicht aus (S_γ).
**(c) Reserve U(τ)/Λ_max(τ):**
| τ | 2° | 3° | 4° | 6° | 8° | 11° | 15° | 20° | 30° | 60° | 103° |
|---|---|---|---|---|---|---|---|---|---|---|---|
| U/Λ_max | 31 | 23 | 17 | **14** | 16 | 20 | 26 | 35 | ≥ 50 | ≥ 44 | ≥ 41 |
**Reserve: kleinster Wert ≈ 14 (bei τ ≈ 6°, mit konservativem Λ); ≈ 20–23 mit Ws Λ. Bei L → ∞ (U linear, Λ ~ τ^{1.5}), fern von L ≥ 40.**
Ein Beweis über Wachstumsraten muss also das Profil der Eckmoden auf Γin bis auf einen Faktor ≈ 14 genau nach oben abschätzen (und bei L den
Exponenten p > 1 zeigen) — nicht polynomial genau.
**Befund zu Ks Linse (Gegenstimme):** mit dem konservativen Λ liegt Ks Linse (S = 0.0166) bei 4°–8° nur bei s/Λ ≈ 0.92–1.04 — möglicherweise auf
oder unter der asymptotischen Untergrenze. Bis N = 1500 hält sie, weil die Raten dort noch klein sind; ob sie für N → ∞ hält, ist offen.
**Empfehlung:** Bezugslinse S = 0.05 (Ks Profilform): s/Λ ≈ 3–5 und U/s ≈ 5–10 überall — mittig im Fenster; ‖Ĵ_N⁻¹‖ = 1.110071 konstant
(N ≤ 2000, Server 22:48).
**Nachtrag 22:44 — Ks Linse bei N = 2000: ‖Ĵ⁻¹‖ = 1.251155, unverändert seit N = 100** (= Ks Bulk-Wert 1.2512). Läge sie bei 4°–8° wirklich
unter Λ (Defizit ≈ 8 %), müsste die Norm bis N = 2000 sichtbar wachsen. ⇒ Das konservative Λ (eigene Anpassung) ist dort zu gross, Ws Λ passt;
Ks Linse ist zulässig. Die Gegenstimme ist damit entkräftet; S = 0.05 bleibt die Linse mit mehr Abstand.

### 3.3 Familie in w-Form (Plan A3, `v17_familie_w.py`, `v18_familie_neumann.py`, Logs sg/v17_*, sg/v18_*): (S_γ) für T_N𝒥_P, P = f + η
Kubisches f, M_P = c₁ + c₂(P + f) + c₃(P² + Pf + f²) exakt (keine Reihendivision). η = ε·g(x/ρ), sup_{D̄}|η| ≈ ε; g = „glatt“ (w²), „zufällig“
(alle Grade, Σr_k w^k/√(N+1)), „hoch“ (nur Grade ≥ 0.9N). Auf der Linse (|x| > ρ) ist η um bis zu e^{SN} grösser (gemessen e^{0.05N}/√N bei S = 0.05).
**(a) v17, ‖(Ĵ^P)⁻¹‖ (Bezugswert 1.110071):** unverändert bis ε = 1e-3 bei allen Arten und N = 400/800; bei ε = 1e-2: glatt 1.110532,
zufällig 1.110080 (N = 400) bzw. 1.110666 (N = 800), hoch 1.110071.
**(b) v18, Neumann-Grössen (22:23–22:36, Server):** X = Ĵ⁻¹ΔĴ, t_k = ‖X^kĴ⁻¹‖ (Glieder der Reihe (Ĵ^P)⁻¹ = Σ(−X)^kĴ⁻¹, Produkte in arb),
Δ_inv = ‖(Ĵ^P)⁻¹ − Ĵ⁻¹‖, q_R = ‖ΔĴĴ⁻¹‖:
| Art (Linse) | N | Δ_inv/ε | t₂/(t₁ε) | q_R/ε |
|---|---|---|---|---|
| zufällig (0.05) | 300 / 400 / 600 / 800 | 4.34 / 4.88 / 6.92 / **23.9** | 6.6 / 7.7 / 8.9 / 14 | 2e11 / 2e15 / 1e23 / 8e30 |
| zufällig (K, 0.0166) | 300 / 400 / 600 | 5.17 / 5.76 / 8.33 | 6.4 / 7.5 / 9.4 | 1e10 / 1e13 / 3e19 |
| glatt (0.05) | 300 / 400 / 600 / 800 | 15.5 / 20.3 / 30.0 / 39.7 | 64 / 84 / 124 / 164 | 4e5 / 2e7 / 3e10 / 5e13 |
| hoch (0.05) | 300 / 400 / 600 / 800 / 1000 / 1200 | 0.160 / 0.323 / 1.84 / **23.4 / 163 / 1107** | ≤ 1e-5 | 1e5 / 9e6 / 2e11 / 4e15 / – |
| hoch (K, 0.0166) | 400 / 600 / 800 | 0.405 / 2.78 / 37.3 | – | – |
| hoch (0.1) | 400 / 600 / 800 | 0.265 / 1.20 / 16.3 | – | – |
| L-gerichtet hoch (0.05) | 400 / 600 / 800 | 0.072 / 0.133 / 0.587 | – | – |
| L-gerichtet alle Grade, k^{−1.5} (0.05) | 400 / 600 / 800 | 17.8 / 26.6 / 35.5 | – | – |
(„hoch“: sup_{D̄}|g| ≈ 0.6–0.7, „zufällig“ ≈ 2.7–2.9, L-gerichtete auf sup = 1 normiert; „–“ = Kurzlauf nur Δ_inv.)
**Gradbänder (22:49, Zufallskoeffizienten nur für ξ = k/N im Band, sup = 1, je eine Realisierung), Δ_inv/ε bei N = 400 / 600 / 800:**
| Band ξ | [0.1, 0.2) | [0.3, 0.4) | [0.5, 0.6) | [0.6, 0.7) | [0.7, 0.8) | [0.8, 0.9) | [0.9, 1.0) |
|---|---|---|---|---|---|---|---|
| Δ_inv/ε | 5.5 / 6.9 / 17.9 | 4.5 / 1.9 / 2.6 | 0.49 / 2.4 / 1.5 | 1.37 / 1.44 / 1.55 | 1.03 / 0.92 / 0.98 | 0.39 / 1.03 / 4.5 | 0.50 / 2.8 / **34.2** |
⇒ Die exponentielle Lipschitz-Konstante sitzt **nur in den obersten Graden ξ ≥ 0.8** (bei ξ ≥ 0.9 Rate ≈ 0.0125); die Mittelbänder
0.3 ≤ ξ < 0.8 sind beschränkt (≈ 1–2), das untere Band wächst langsam (Ableitungsverlust; Einzelrealisierung, Faktor ~2 Streuung).
**Befund (numerisch):**
1. Die naive Neumann-Reihe in ℓ² scheitert: q_R/ε und q_L/ε wachsen wie e^{≈0.18N}, weil Ĵ riesige Einträge hat (‖Ĵ‖ ≈ 1.6^N: auf Γout ist
   |f∘k| bis 1.6ρ, die Spalten (f∘k/ρ)^j wachsen entsprechend) und ΔĴ dort mitwächst.
2. Trotzdem ist die wahre Änderung der Inversen linear in ε: Δ_inv = t₁ (erste Ordnung) bis auf < 1 %, die Reihe Σ X^kĴ⁻¹ konvergiert mit
   Quotient t₂/t₁ ≈ 10ε (zufällig), ≈ 0.2Nε (glatt), ≤ 1e-5ε (hoch). Die Dämpfung sitzt im rechten Faktor Ĵ⁻¹. ⇒ Ein Beweis der Familien-
   Erweiterung braucht eine **gewichtete Norm** (etwa δ beschränkt auf f(Int γ)), in der X klein ist; ℓ² und die Graphennorm von Ĵ reichen nicht.
3. **Lipschitz-Konstante Lip_N = Δ_inv/ε von P ↦ (Ĵ^P)⁻¹ (bezüglich sup_{D̄}|P − f|):** glatte und auf L gerichtete Störungen mit allen Graden:
   **linear in N** (≈ 0.045N — Ableitungsverlust δ∘P − δ∘f ≈ δ'·η, Bernstein-Faktor N). Störungen hoher Grade (≥ 0.9N): **exponentiell,
   Lip_N ≈ e^{0.010…0.013·N}** (Raten 400→600→800→1000→1200: 0.0087, 0.0127, 0.0097, 0.0096 — keine Sättigung bis 1200), **unabhängig von der Linsendicke** (S = 0.0166 / 0.05 / 0.1
   bei N = 800: 37 / 23 / 16) — die Grösse der Eckrate c ≈ 0.012 taucht in der P-Ableitung wieder auf, obwohl ‖Ĵ⁻¹‖ selbst beschränkt ist.
   Auf L gerichtete hohe Grade: kleiner (Faktor 5–40) und langsamer (Raten 0.003–0.007), wachsen aber auch. Ab N ≈ 800 dominiert der hohe
   Anteil auch die Zufallsstörung (zufällig 6.9 → 23.9). Gradbänder: nur ξ = k/N ≥ 0.8 wächst exponentiell, 0.3 ≤ ξ < 0.8 bleibt beschränkt.
4. **Realistische Grösse:** Profil von p_N − f_K aus den zertifizierten Boxen (`v19_profil_pN.py`): das ℓ¹_ρ-Gewicht liegt bei k ∈ [0.5N, 0.9N]
   (6.2e-13 bei N = 400), das Band k ∈ [0.9N, N] trägt 4.5e-14 (≈ 7 %), Schwanz k > N 5.3e-14; Abfall ≈ N^{−5.5} (passt zu |b_k| ~ k^{−6.45}).
   ε·Lip_N ist bei N = 400 ≈ 1e-14. Hochgerechnet (Lip ~ e^{0.0104N}, ε_hoch ~ N^{−5.45}) wird ε·Lip ≈ 1 erst bei N ≈ 5000. ⇒ Für alle praktisch
   gerechneten N ist die Familie harmlos; **für einen Beweis für alle N reicht reine Störungsrechnung um f_K nicht**, solange Lip_N exponentiell
   wächst — Lip_N sättigt bis N = 1200 nicht; also muss die Familien-Erweiterung die Struktur von p_N − f_K (L-gerichtet,
   Schwanzgestalt) oder eine Fortsetzung in N benutzen.
5. Ks dünne Linse verhält sich wie S = 0.05 (Werte 20–60 % grösser).

### 3.4 (S_γ) für das echte p_N (Auftrag Leitung 22:52) — `v20_pN_S.py`, Logs sg/v20_* — **hält bis N = 1000, Konstante wie bei f_K**
Ĵ_N^{p_N} = T_N^w[(δ∘p_N + (p_N'∘p_N)δ)∘k] (Tangente des Newton-Verfahrens, p_N'∘p_N exakt als Polynomkomposition), p_N aus
gemeinsam/pN_N100_bis_1000.json (60 Stellen), Vergleich Ĵ^{f_K} (kneser_halbexp_koeff_4000.json):
| N | 100 | 200 | 400 | 600 | 800 | 1000 |
|---|---|---|---|---|---|---|
| ‖(Ĵ^{p_N})⁻¹‖, Ks Linse | 1.250337 | 1.250337 | 1.250337 | 1.250337 | 1.250337 | **1.250337** |
| ‖(Ĵ^{p_N})⁻¹‖, S = 0.05 | 1.107257 | 1.107257 | 1.107257 | 1.107257 | 1.107257 | **1.107257** |
| ε_N = Σ_{k≤N}|p_k − a_k|ρ^k | 6.8e-10 | 2.4e-11 | 8.4e-13 | 1.2e-13 | 2.9e-14 | 1.0e-14 |
| Δ_inv/ε_N (Ks Linse) | 2.4 | 3.9 | 6.7 | 9.5 | 12.0 | 14.6 |
**Befund (numerisch):** (S_γ) gilt für das echte p_N mit derselben Konstante wie für f_K (6 Stellen gleich). Die Änderung gegenüber f_K wächst
nur linear, Δ_inv ≈ 0.015·N·ε_N, und ε_N fällt wie N^{−4.8} (v22: Profil von p_N − f_K, Gewicht bei k ∈ [0.5N, 0.8N], oberste Grade
[0.99N, N] ≈ N^{−5.6}). ⇒ Die Homotopie innerhalb der Lösungsmenge (ARCHITEKTUR §4.2 H1, §5.1(ii): (H) nur an Fast-Lösungen) ist numerisch
gedeckt; dort sind die Störungen kohärent.
**Geometrischer Faktor (v27_geometrie.py, numerisch):** in K1″(c) G_N(t) ≤ |t|·Q_k·√(N+1)·‖Ĵ_N⁻¹‖ ist |t|·Q_k ≈ 1.11 (Ks Linse; min_∂𝔻|k'| = 1.258,
max_γ|t| = 1.398) bzw. 1.36 (S = 0.05; min|k'| = 1.061, max|t| = 1.445); c_k = 0.347 bzw. 0.361. Q_k ist hier über min|k'| und 2e5 Zufallspaare
geschätzt (Untergrenze des Supremums), nicht rigoros eingeschlossen.

### 3.5 Gewichtete Norm für δP (Auftrag Leitung 22:52, Punkt 2) — `v21_grad_antwort.py`, `v23_kohaerent.py`, `v25_m_differenz.py`
**(a) Einzelgrade (Kneser, Tangente, P = f_K + ε(x/ρ)^k):** r_k(N) = ‖∂_ε(Ĵ^P)⁻¹‖:
| ξ = k/N | ≤ 0.6 | 0.7 | 0.8 | 0.9 | 0.95 | 1 |
|---|---|---|---|---|---|---|
| N = 100 | 6.5–15 | 19 | 46 | 1.4e3 | 8.4e4 | 1.2e7 |
| N = 200 | 12–32 | 47 | 353 | 5.5e5 | 3.9e9 | 1.2e14 |
| N = 300 | 17–49 | 92 | 2.2e3 | 2.1e8 | 1.9e14 | 1.1e21 |
| N = 400 | 23–67 | 170 | 1.3e4 | 7.8e10 | 9.0e18 | 1.1e28 |
r_k ≈ e^{N·g(k/N)} mit g(1) = 0.161 (sehr genau), g(0.95) ≈ 0.11, g(0.9) ≈ 0.063, g(0.8) ≈ 0.024; für ξ ≤ 0.6 linear in N (≈ 0.11N).
**Unabhängig von der Linse** (S = 0.05 auf 1 % gleich). Ursache: der Tangenten-Multiplikator P'∘P enthält η'∘f, und |f| reicht auf Γout bis 1.6ρ.
(Im kubischen Sekantenmodell §3.3 fehlt dieser Term; dort ist das Wachstum nur e^{0.010N}.)
**(b) Folgerung — keine diagonale Gewichtsnorm leistet beides:** (a) verlangt ω_k ≳ r_k(N)/poly(N); dann ist für das echte p_N − f_K
Σ_k d_k ω_k ≥ d_N r_N/poly = 7e-5 (N = 100), 12 (200), 8e6 (300) — wächst ⇒ (b) verletzt. Auch die Graphennorm ‖T_N J_∞η‖ scheitert an (b):
für η = p_N − T_N f_K ist ‖J̃_Nξ‖ = 3e9 (N = 100), weil T_N J_∞ den Schwanz von f_K (Grade > N) nicht sieht.
**(c) Was die Kleinheit wirklich trägt (Kohärenz):** Die linearisierte Multiplikator-Änderung T_N(η'∘f) für η = p_N − T_N f_K ist riesig
(1.6e11 bei N = 100, 3e30 bei N = 200), die exakte Änderung T_N(p_N'∘p_N) − T_N(f_K'∘f_K) dagegen klein (2.6e-6 bzw. 3.8e-5): der Schwanz von f_K,
mit f komponiert, hebt den grossen Teil auf. Richtungen der Form η = J_N⁻¹(glatt) (Newton-/Homotopie-Richtungen, z. B. die Kopf-Antwort
J⁻¹x^N aus (H-iii)) sind harmlos: Lip pro ℓ¹_ρ-Einheit von η = 3.0 / 6.2 / 12.8 / 19.5 / 25.9 (Kopf, N = 100/200/400/600/800), 5.3 / 9.9 / 19.5 /
30.5 / 43.0 (g = 1), 1.1 / 2.0 / 4.0 / 6.2 / 8.7 (glatte Zufallsdaten) — überall **linear in N**.
Ecken der Graphennorm-Kugel (η = J̃_N⁻¹(x/ρ)^k): pro ℓ¹_ρ-Einheit von η ebenfalls linear (max 4.8 / 11.5 / 18.2 bei N = 200/400/600), pro Einheit
Residuum aber bei k ≈ 0.95N wachsend (41 / 539 / 4356 ≈ e^{0.0105N}), weil J̃_N⁻¹ in der x-Koordinate die Eckmode anregt (‖η‖ = 8.5 / 47 / 240).
⇒ **Empfehlung für Ks Kopplungssatz A2c:** Die Familie nicht über eine Norm auf beliebigen Koeffizientenstörungen, sondern über die
Lösungsmenge führen (H1-Homotopie, Fast-Lösungen). Nötig ist dann nur: (S_γ) gleichmässig an Fast-Lösungen, plus Lipschitz längs kohärenter
Richtungen J(P)⁻¹(glatt) — numerisch linear in N. Die schädlichen Richtungen (Einzelgrade k ≥ 0.8N, oder Koeffizientenstörung ohne den
Schwanz von f_K) kommen in 𝒮_N ∩ B_N nicht vor.

### 3.6 Ks Plan A3 (i)–(iii): Eck-Korrektur D_N = Ĵ_N⁻¹ − P_NĴ_{3N}⁻¹P_N (`v26_widom.py`, Logs sg/v26_*, 23:05)
| Modell / Linse | N | ‖D_N‖ und σ₂ | σ₃ / σ₄ | führende Vektoren (x / w) | ‖χD_Nχ‖ fern von L (δ = 10° / 40°) |
|---|---|---|---|---|---|
| kubisch, Ks Linse | 100 / 150 / 200 | 0.4774, 0.4602 / 0.4922, 0.4842 / **0.5123, 0.5057** | 0.281/0.269 · 0.220/0.209 · 0.235/0.223 | ξ-Spitze 0.86–0.92, ≥ 93 % in ξ > 0.6; Winkel ±(77–80)° | 0.270, 0.185 / 0.125, 0.104 / 0.060, 0.050 |
| Kneser, Ks Linse | 100 / 150 / 200 | 0.4774, 0.4603 / 0.4921, 0.4841 / 0.5122, 0.5055 | gleich (3–4 Stellen) | gleich | gleich |
| kubisch N = 150: τ₁ = 5°, p = 1.15 (dicker an L) | 150 | 0.329 | 0.142/0.126 | ξ 0.86, 77° | 0.123 / 0.104 |
| τ₁ = 25° (dünner an L) · p = 2.0 | 150 | 0.637 · 0.573 | | | 0.198 · 0.123 |
| S = 0.05 · S = 0.1 | 150 | 0.319 · 0.241 | 0.070/0.065 · 0.041/0.040 | ξ 0.85–0.89, 77° | 4.4e-3, 9.8e-4 · 1.7e-3, 3.2e-4 |
**Befunde (numerisch):** (i) **Universalität:** kubisch und Kneser geben dieselbe Eck-Korrektur auf 3–4 Stellen (Singulärwerte, Paare L/L̄,
Lage). Ks Werte (0.512/0.506 bei N = 200) unabhängig reproduziert. D_N hängt also nur von den Daten an L ab (f(L) = L, f'(L) = √L) und von der Linse.
(ii) **Linsenform:** ‖D_N‖ fällt, wenn die Linse an L dicker wird (τ₁ = 5°: 0.33 statt 0.49) oder insgesamt höher (S = 0.05/0.1: 0.32/0.24);
bis S = 0.1 keine Sättigung, aber abnehmend. ‖D_N‖ wächst mit N langsam (0.477 → 0.512 von 100 bis 200).
(iii) **„Kante fern von L trivial“ bestätigt:** mit Winkel-Abschneidung um L, L̄ fällt der Rest exponentiell, etwa wie e^{−S·N} (Ks Linse: Halbierung
je ΔN = 50 ⇒ Rate 0.015 ≈ S = 0.0166; S = 0.05: schon bei N = 150 bei 4e-3) und hängt kaum von δ (10°–40°) ab. D_N sitzt also ganz an der Ecke.

### 3.7 Ks Plan A3 (iv): Schärfe der Schleifenbedingung profilweise — blinder Kerben-Test (`v12_konform_S.py` Form „kerbe“, Logs sg/v12k_*)
Ks Linsenform mit S = 0.05 (sonst überall weit über Λ), dazu eine Gauss-Kerbe, die s bei τ₀ auf s_min unter Λ(τ₀) drückt (beidseitig L und L̄).
- Vorhersage vor dem Lauf (23:05): (a) τ₀ = 15°, Breite w = 15°, s_min = 0.003 ⇒ Rate 0.0030–0.0035, Einsatz N ≈ 700–900; (b) wie (a), aber w = 4° ⇒
  kleinere Rate, Einsatz > 1200 oder unsichtbar; (c) τ₀ = 60°, w = 15°, s_min = 0.005 ⇒ Rate ≈ 0.004, Einsatz ≈ 600.
| Kerbe | N = 400 | 600 | 800 | 1000 | 1200 |
|---|---|---|---|---|---|
| (a) 15°, w 15°, 0.003 | 1.1154 | 1.1154 | 1.1308 | **2.2729** | **4.8274** (Rate 0.0038) |
| (b) 15°, w 4°, 0.003 | 1.1115 | 1.1115 | 1.1115 | 1.1115 (σ₂ steigt 1.0085 → 1.0188) | **2.0807** (Einsatz zwischen 1000 und 1200) |
| (c) 60°, w 15°, 0.005 | 1.1907 | 1.1907 | **2.5074** | **7.8596** | **25.4886** (Rate 0.0059) |
**Befund (numerisch):** Die Bedingung wirkt **profilweise**: eine Kerbe unter Λ(τ₀) in einem Fenster von ±15° löst Wachstum aus, auch wenn die Linse
sonst überall weit darüber liegt. (a) trifft die Vorhersage (Einsatz zwischen 800 und 1000, Rate 0.0038 statt 0.0030–0.0035). (c) wächst schneller
als vorhergesagt (Rate 0.0058 statt 0.004) ⇒ Λ(60°) ≈ 0.011–0.012, flacher als die lineare Interpolation zwischen Ws Werten bei 30° und 103°.
(b) Die schmale Kerbe (±4°) setzt später ein (zwischen 1000 und 1200; Vorhersage „> 1200 oder unsichtbar“ knapp verfehlt): die Eckmoden haben bei
endlichem N eine Winkelbreite, die mit N schrumpft — sie lösen schmalere Kerben erst bei grösserem N auf. ⇒ **Die Schleifenbedingung ist
asymptotisch punktweise scharf** (jede Stelle mit s(τ₀) < Λ(τ₀) führt irgendwann zu Wachstum mit Rate ≈ Λ(τ₀) − s(τ₀)); bei endlichem N wirkt ein
Winkelmittel. Für Vermutung S passt die Form s ≥ (1 + η)Λ überall — sie ist hinreichend und (bis auf η) nötig.

### 3.8 Familie längs der Homotopie (Auftrag Leitung 23:20, Punkt 2) — `v28_familie_homotopie.py` (Fassung b), Logs sg/v28b_*
Kneser, Linse wie v12. Alle Kompositionen mit f_K über die Taylor-Entwicklung um f(x) (f∘P = eˣ + Σ g_nη^n/n!, g_n = f^{(n)}∘f aus m = eˣ/f');
Fassung a (Σ_k a_k P^k) war ab N = 400 unbrauchbar, weil das Umentwickeln um P(0) ≈ 0.5 die Datenfehler der a_k (k > 500, relativ 1e-21) um bis zu
1e144 verstärkt.
| N | Sekante 𝒥_{p_N} | λ = 0.25 / 0.5 / 0.75: Tangente J_N(P(λ)) und Sekante 𝒥_{P(λ)} | σ-Homotopie N → N+1 (s = 0 / 0.5 / 1) |
|---|---|---|---|
| 100 | 1.250337 | 1.250337 (λ = 0.5) | 1.250337 / 1.250337 / 1.250337 |
| 200 | 1.250337 | alle 1.250337 | alle 1.250337 |
| 400 | 1.250337 (S = 0.05: 1.107257) | alle 1.250337 (S = 0.05, λ = 0.5: 1.107257) | alle 1.250337 |
| 600 | 1.250337 | alle 1.250337 | alle 1.250337 (mit nachgeschärftem p_N, v28c) |
| 800 | 1.250337 | alle 1.250337 | alle 1.250337 (v28c) |
| 1000 | 1.250337 (S = 0.05: 1.107257) | λ = 0.5: 1.250337 (S = 0.05: 1.107257) | alle 1.250337 (v28c) |
Q(λ) − λp_N ist von zweiter Ordnung (1e-20 bei N = 100 … 1.5e-29 bei N = 1000): die Lösungsmenge ist zwischen f_K und p_N fast gerade.
**Befund (numerisch):** Für Sekanten und Tangenten an allen geprüften Fast-Lösungen ist ‖Ĵ⁻¹‖ auf 6–7 Stellen gleich dem Wert bei f_K.
σ-Homotopie: mit den 60-stelligen p_N divergierte Newton bei N ≥ 600 für s < 1 — der Kopfrest c_{N+1} von pad(p_N) war Rundungsrauschen (die
Komposition P∘P verstärkt Koeffizientenfehler um bis zu 1e(0.36N)). Nach Nachschärfen von p_N auf ≈ 0.45N Stellen (v28c) konvergiert Newton bei allen s
in 1–3 Schritten, und alle Normen sind 1.250337.

### 3.9 Konstanten von K1 explizit (Auftrag Leitung 23:20, Punkt 1) — `v30_geometrie_Kgamma.py`, `v31_N0_bilanz.py`, `v32_geometrie_zert.py`
**Geometrie (rigoros, v32, arb mit Vs (K2)-Datei, Ks Linse S = 0.0166, τ₁ = 12°, p = 1.3):**
(a) f(ℬ) ⊂ D und f(Γin) ⊂ D: Intervalle für τ ≥ 0.2° (80 bzw. 12 Teilintervalle), für τ < 0.2° lokal über F = e^{iθ} + μΔ + R,
|R| ≤ ½C₂|Δ|², C₂ ≤ 7.99 (F'' über die Darstellung F = log(ρF(w))/ρ) — bestanden. (b) f(Γout) trifft den Schlitz S_L nicht (arg f < θ; lokal
Im(μ(e^{−iτ} − 1)) ≈ −τRe μ) ⇒ S_L ∩ f(D̄ ∪ ℬ̄) = ∅ (S_L liegt in der unbeschränkten Komponente von ℂ∖f(∂D)); Linse ∩ Schlitz = ∅ trivial.
(c) M_γ = max_{D̄∪ℬ̄}|m| ≤ 2.399 (Intervalle; nahe L |m| ≤ e^{Re x}/(|μ| − C₂|Δ|) ≤ 1.21).
**Green-Funktion g_γ von Ĉ∖K_γ (numerisch, v30):** K_γ sternförmig; äussere Theodorsen-Abbildung (M = 2^15, Randfehler 3e-14), Kapazität 1.212ρ;
Lip(g_γ) = 0.648 pro x-Einheit (Maximum am linken Linsenscheitel); längs S_L: g_γ(L(1+u)) ≥ 0.70·u^{1/α'} für u ∈ [0.005, 0.3] (Werte 0.73–0.77,
α' = 0.787172 wie K'); dist(L(1+u), K_γ) ≥ 0.786·ρu; Niveaulinie Γ_{0.09}: Länge 11.49, trifft den Schlitz bei u₀ = 0.19 (< 0.2, dort gilt
Ks/Ps Sprungschranke C_J = 0.95), max|f_K| darauf 4.31.
**Kette (v31):** K1″(a) ⇒ |E_N| ≤ [√(N+1)·N·B(1+M_γ) + (N+1)/ρ]·e^{(N+2)Lip/N}·e^{−(N+2)g_γ} + Schwanz (Bernstein–Walsh in Ĉ∖K_γ^{1/N}),
nahe L (dist < 1/N) die 1/dist-Form; Satz P-I mit Schlitz bis u₀ und Kontur Γ_{0.09}; Familie: Faktor 2 auf den B-Anteil (wie C_R^P = 2C_R in X4a);
R_N = 6·sup|Π^P_Nf − f|; Bedingung X4a θ_N = A_θ·N·R_N ≤ ½ mit A_θ = 12.93 (Lemma X8′, c₀ = 0.641 aus (a″)_∞).
| N | 300 | 318 | 350 | 400 | 500 | 1000 | 2000 |
|---|---|---|---|---|---|---|---|
| θ_N / B | 0.033 | 0.016 | 0.0095 | 0.0065 | 0.0037 | 6.4e-4 | 1.1e-4 |
| B_max(N) | 15.0 | **30.8** | 52.7 | **76.7** | 136 | 787 | 4500 |
(Bis N ≈ 350 dominiert der Konturteil, danach der Schlitz; B_max wächst etwa wie N^{2.6}.)
**Satz V2 (bedingt; Konstanten der Green-Funktion numerisch).** Es gelte (S_γ) mit Ks Linse in der Form sup_{N≥N₀}‖Ĵ_N⁻¹‖ ≤ B an f_K, dazu die
Familie mit Faktor 2 (Annahme (Fam), siehe unten). Ist **B ≤ 30.8 (N₀ = 318) bzw. B ≤ 76.7 (N₀ = 400)**, so gilt θ_N ≤ ½ für alle N ≥ N₀,
der Ballradius ist R_318 = 6(6.6e-7·B + 4e-9) ≈ 5e-6 bei B = 1.25 (der B-Anteil enthält schon den Familienfaktor 2) — weit über dem zertifizierten sup|p_{N₀} − f_K| ≤ 3.6e-12 —, und mit H1/H2
und dem Zweig-Zertifikat bis 400 folgt die Vermutung. **Gemessen ist B = 1.2503 (bis N = 2000 an f_K, bis N = 1000 an p_N und allen
Fast-Lösungen) — Reserve Faktor 25 (N₀ = 318) bzw. 61 (N₀ = 400).**
**Was in Satz V2 noch steckt (ehrlich):** (i) die Green-Konstanten Lip, c_g, c_d, Γ_d, M_d sind numerisch (wie Ps K'-Werte), C_J = 0.95 [num.],
A_θ-Eingänge min|f'|, max|m/f'|, M2 [num.]; (ii) Lemmata X5, X7, X8′, Satz X4a sind B (kein zweiter Leser); (iii) **(Fam)**: dass sich die
Eckantwort für P ∈ 𝒮_N ∩ B_N höchstens um Faktor 2 ändert. Der naheliegende Beweis über K1″(a) mit 𝒥_P braucht W_P = 1/(η − P) + M_P/(η − x)
analytisch auf Int γ, also P(ℬ) nahe f(ℬ). Das gilt nur bis N_krit(S): gemessen |p_N − f_K| auf Ks Linse ≤ 2.4e-12 bei N = 1000, auf S = 0.05
schon 0.24 bei N = 1000 (v29). Für beliebig grosse N muss (Fam) anders gezeigt werden (Ks Kopplungssatz über die Lösungsmenge).

### 3.9b (S_γ) am echten p_N bei N = 1500 und 2000 (Auftrag Leitung 23:45, Vorrang) — `v33_pN_sprung.py`, `v34_norm_pN_gross.py`
p_1500 per Newton-Sprung von pad(p_1000) (2650 Bit; Schritte ‖δ‖ = 9.1e-15, 1.1e-27, 8.7e-54, 5.3e-106 — quadratisch; Residuum ‖T_N[P∘P − eˣ]‖ = 1.4e-110),
Gegenprobe der Methode: Sprung 100 → 200 trifft Strang As p_200 auf 4e-61. p_2000 per Sprung 1000 → 2000 (3400 Bit, ‖δ‖ = 1.0e-14, 5.0e-27, 3.3e-52, …;
Residuum 1.9e-106), ε_2000 = 3.5e-16; zweiter, unabhängiger Sprung 1500 → 2000 (Residuum 1.5e-226) stimmt damit auf 2.3e-206 überein (max_k|Δp_k|ρ^k).
| N | Linse | Tangente an p_N | Sekante an p_N | f_K | Δ_inv/ε_N |
|---|---|---|---|---|---|
| 1500 | Ks (S = 0.0166) | **1.250337** | **1.250337** | 1.250337 | 20.8 |
| 1500 | S = 0.05 | **1.107257** | **1.107257** | 1.107257 | 19.3 |
| 1500 | S = 0.1 | **0.965153** | **0.965153** | 0.965153 | 17.6 |
| 2000 | Ks / 0.05 / 0.1 | **1.250337 / 1.107257 / 0.965153** | **1.250337 / 1.107257 / 0.965153** | 1.250337 / 1.107257 / 0.965153 | 26.7 (Ks) |
**Befund:** Die Konstante hält bei N = 1500 und 2000 für alle drei Linsen, auch für die dicke Linse S = 0.1, wo das Polynom auf der Linse weit von
f_K abweicht. Δ_inv/ε_N wächst weiter linear (≈ 0.0134·N; 20.8 bei 1500, 26.7 bei 2000), ε_N = 1.4e-15 bzw. 3.5e-16. **Kein Anstieg — die
Linsen-Abweichung ist für (S_γ) an der Familie harmlos** (bestätigt Ks Keimabhängigkeit §10.1).

### 3.10 Datenkorrektur der Familie für h = k_η (Ks Anfrage 23:59) — `v35_datenkorrektur.py`, Logs sg/v35_*
Gemessen ‖Π^{p_N}_Nk_η − Π_Nk_η‖_ℓ¹ρ (Sekante wie in H1, alles als Keim mod x^{N+1}, x-Koordinate, arb) im Verhältnis zu ε_N = ‖p_N − T_Nf_K‖_ℓ¹ρ;
‖Π_Nk_η‖ selbst ist 70 (u = 0.01), 15 (u = 0.05), 4 (u = 0.19), 1–6 auf Γ_{0.09}.
| N | 200 | 400 | 600 | 800 | 1000 | 1500 | 2000 |
|---|---|---|---|---|---|---|---|
| Schlitz u = 0.01 | 1.86 | 1.69 | 1.02 | 0.51 | 0.23 | 0.027 | 0.0027 |
| Schlitz u = 0.05 | 4.8e-3 | 2.0e-5 | 6.0e-8 | 1.6e-10 | 3.7e-13 | 9e-20 | 2e-26 |
| Schlitz u = 0.19 / Γ_{0.09} | ≤ 5e-12 | ≤ 5e-23 | ≤ 4e-34 | ≤ 3e-45 | ≤ 2e-56 | ≤ 2e-84 | ≤ 2e-86 |
**Befund (numerisch):** Die Datenkorrektur ist höchstens 2·ε_N und wird mit N relativ kleiner; fern von L ist sie exponentiell klein (≈ e^{−N g_γ(η)}).
Die „Linsen-Abweichung“ der Polynome taucht hier nicht auf, weil nur Keime eingehen (Ks §10.1). Für (Fam) in Ks Keimform ist das der benötigte Beleg.

### 3.11 Green-Konstanten von Satz V2 streng (Auftrag Leitung 02:05, Punkt 2) — `v39_lune_zert.py`, `v40_lipschitz_zert.py`, `v41_N0_streng.py`
**Vergleichsgebiet (v39, arb):** K* = D₁ ∪ D₂ mit Kreisen durch L, L̄: D₁ = D(−0.09, 1.39813), D₂ = D(1.0, 1.50105). Geprüft: Linsenkurve ⊂ D₁
(102 Intervalle + Ableitung an L), f(Γout) ⊂ D₂ (202 Intervalle + Richtung an L), Γin ⊂ D₁ und Γout ⊂ D₂ (Kreisbüschel), f(Γin), f(ℬ) ⊂ D (v32)
⇒ K_γ ⊂ K*. Unter T(z) = (z − L)/(z − L̄) ist Ĉ∖K* ein Sektor der Öffnung απ, α = 0.75561 (K': 0.787) ⇒ g_{K*} explizit.
⇒ **g_γ ≥ g_{K*} (streng)**: längs des Schlitzes 0.89–0.91 des numerischen g_γ (Tabelle mit 241 Werten, monoton); Kontur Γ*_d = {g_{K*} = d},
d = 0.078: dort g_γ ≥ 0.078, Länge ≤ 12.18, Schlitz-Schnitt u₀ = 0.196 (< 0.2, Sprungschranke C_J gilt); dist(L(1+u), K_γ) ≥ dist(L(1+u), K*) explizit.
(Nebenbefund: bei u = 0.001 lag der numerische Wert von v30 unter der strengen Unterschranke — die Theodorsen-Abbildung löst die Ecke dort nicht auf.)
**Lipschitz (v40):** innere Scheiben (rollende Scheibe): Krümmung der Linsenkurve ≤ 0.859, von f(Γout) ≤ 0.861 (arb, f'' mit Koeffizienten- und
Schwanzradien) ⇒ r₀ ≥ 1.162 ⇒ **Lip(g_γ) ≤ 0.861** (numerisch war 0.648).
**M_d (v41):** |f_K| auf Γ*_d über f = exp∘f∘log (log x ∈ D̄), f = log∘f∘exp (eˣ ∈ D̄) und nahe dem Schlitz-Schnitt f = exp∘exp∘f∘log∘log (zweimal log
zieht zu L hin), 6001 Bogenkugeln (Sehne ≤ 0.071, Radius Sehne/2 + κ·Sehne²/8, κ ≤ 0.7) ⇒ **M_d ≤ 5.03**.
**Ergebnis (v41):** mit diesen strengen Konstanten und M_γ ≤ 2.40:
| N₀ | 318 | 350 | 400 | 500 | 1000 | 2000 |
|---|---|---|---|---|---|---|
| B_max(N) | 1.36 | 7.0 | **17.9** | 31.8 | 176 | 940 |
(Mit den strengen A_θ-Eingängen neu gerechnet, 08:36: 1.355 / 17.88 — unverändert.) Xs Lemma X8/X8′ habe ich gelesen: Logik richtig; formal ungenau ist nur
„A' ≈ A_θ“ in der Faà-di-Bruno-Verschiebung — wegen N·R ≈ 1e-9 ohne Wirkung (jedes A' ≤ 10·A_θ genügt).
**Satz V2 (streng bis auf (Fam), C_J und Xs Lemmata X5/X7/X8′/X4a; die A_θ-Eingänge sind seit 08:35 streng, v43: min|f'| ≥ 0.34558, max|m/f'| ≤ 2.11804, M2 ≤ 0.81835 ⇒ A_θ ≤ 12.934):** Gilt sup_{N≥400}‖Ĵ_N⁻¹‖ ≤ B mit **B ≤ 17.9**, so folgt die
Vermutung (mit dem Zweig-Zertifikat bis 400). Gemessen B = 1.2503 — Reserve Faktor 14. Ab N₀ = 318 bliebe nur B ≤ 1.36 (knapp über 1.25; der
Konturteil e^{−Nd} mit d ≤ 0.08 dominiert, d ist durch den Gültigkeitsbereich u ≤ 0.2 von C_J begrenzt) — das Zertifikat bis 400 ist also wichtig.
Länge von Γ*_d: Sehnensumme ×1.001 (Fehler ≤ (κ·Sehne)²/24 ≈ 1e-4); die Krümmungsschranke κ ≤ 0.7 der expliziten Kurve ist numerisch (0.66 gemessen).

### 3.12 VE4 — faktorisierter Schur-Test mit Diagonalgewicht (Hauptlinie, Auftrag Leitung 02:05) — `v38_ul_gewicht.py`, `v42_gewicht_opt.py`
**Verfahren:** X = Ĵ_N⁻¹ (arb), UL-Zerlegung X = U'L' per Rückwärts-Elimination auf X (Satz K5: X − Kopf_N = Ĵ_{N−1}⁻¹ ⊕ 0; Spalte k von U' = Ĵ_k⁻¹e_k,
Zeile k von L' = e_kᵀĴ_k⁻¹/(Ĵ_k⁻¹)_kk), Mittelpunkt-Arithmetik mit P Bit; dann min über Diagonalgewichte W von ‖U'W⁻¹‖·‖WL'‖ (konvex in log W,
Subgradient mit führenden Singulärvektoren). **Präzisionsfalle:** P = 0.25N + 300 reicht ab N ≈ 900 nicht (N = 1000: 2.78 statt 1.49) —
die Elimination löscht stärker aus als die Kopf-Dynamik e^{±0.07N} vermuten lässt; P = N + 200 ist geprüft stabil (N = 1000: 900 und 1300 Bit gleich).
| Linse / Modell | N = 400 | 600 | 800 | 900 | 1000 | 1500 | 2000 |
|---|---|---|---|---|---|---|---|
| tangential (Ks Linse), kubisch | 1.334 | 1.335 | 1.338 | 1.392 | 1.490 | VE4_T1500 | VE4_T2000 |
| tangential, Kneser | | | | | 1.482 | | VE4_TK2000 |
| Keil (Steigung 0.5), kubisch | | | | | 1.332 | VE4_K1500 | VE4_K2000 |
| Keil, Kneser | | | | | 1.339 | | VE4_KK2000 |
(Untergrenze ist ‖Ĵ_N⁻¹‖ = 1.2512 bzw. 1.2505; alle Werte mit P ≥ N + 200 bzw. bei N ≤ 600 mit P ≥ 450 nachgeprüft. Ohne Gewicht: ‖U'‖·‖L'‖ = 50 / 790 / 2.5e4 / 7.9e5 bei N = 60/100/150/200.)
**Gestalt des Gewichts:** log W_k = ½(log‖u_k‖ − log‖ℓ_k‖) + x_k: der Hauptteil ist der **Kopf-Ausgleich** W_k ≈ (‖c_k‖/‖r̂_k‖)^{1/2}, Steigung
−0.071…−0.073 pro Index (kubisch und Kneser gleich, beide Linsen gleich: log‖c_k‖ ≈ −0.072k, log‖r̂_k‖ ≈ +0.070k); der Feinanteil x_k ist beschränkt
(Spanne 0.6–0.7, ohne Trend). **Universell:** kubisch und Kneser stimmen bei N = 1000 auf 1 % überein (1.490/1.482 tangential, 1.327/1.339 Keil).
**Struktur:** Bei der tangentialen Linse sitzt der führende Singulärvektor von U'W⁻¹ bis N = 800 bei kleinen Indizes (k/N ≈ 0.02, glatter Teil),
bei grösserem N wandert Masse ins Eckband (k/N ≈ 0.8); bei der Keil-Linse bleibt er bei kleinen Indizes.
**Welche Aussage man beweisen müsste (Vermutung VE4):** Für die unendliche LU-Zerlegung Ĵ_∞ = L̂U (existiert, weil alle Abschnitte invertierbar sind)
und W_k := (‖c_k‖/‖r̂_k‖)^{1/2} sind die beiden UNENDLICHEN Dreiecksoperatoren U⁻¹W⁻¹ und WL̂⁻¹ auf ℓ² beschränkt. Dann gilt für jedes N
‖Ĵ_N⁻¹‖ = ‖(U⁻¹W⁻¹)_N(WL̂⁻¹)_N‖ ≤ ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ — (S_γ) ohne Abschnitte. Gleichwertig: die normierten Primal-Köpfe c_k/W_k und die normierten
adjungierten Köpfe W_k r̂_k bilden jeweils eine Bessel-Folge. **Verhältnis zu (H1):** VE4 ⇒ (H1) (jeder Kopfterm ‖c_k⊗r̂_k‖ ≤ ‖U⁻¹W⁻¹‖‖WL̂⁻¹‖) und
⇒ (S_γ). Die exakte Auslöschung der beiden Raten (Ks §11.3) steckt ganz im Gewicht; danach sind beide Teilaussagen ratenfrei und LINEAR (je eine
Familie), statt der bilinearen Fast-Orthogonalität der Kopfprodukte (Cotlar–Stein). Leichter ist VE4 nur, wenn man die n-te-Wurzel-Asymptotik der Köpfe
genau genug kennt, um W festzulegen — das ist Ks (E₀); der polynomiale Teil (E₁) wird zu zwei getrennten Bessel-Schranken.

### 3.13 τ_N-Test (Ks Anfrage 03:05, Vorrang 12.10. 08:25): trägt eine feste Linse (H1)/(S_γ) für alle N? — aus den UL-Läufen v38 (Köpfe aller k ≤ N)
τ_k := ‖c_k‖·‖ℓ_kᵀC⁻¹‖ = ‖Ĵ_k⁻¹ − (Ĵ_{k−1}⁻¹ ⊕ 0)‖ (Kopfterm in w); die Rückwärts-Elimination auf Ĵ_N⁻¹ liefert alle k ≤ N exakt (geschachtelt), P = N + 200 Bit.
Bis N = 1000 stimmen Ks Werte (k12) auf alle Stellen (0.0448 / 0.0287 / 0.0212 / 0.0177 / 0.0152 / 0.0149).
| Linse | k = 300 | 500 | 800 | 1000 | 1200 | 1500 | 1800 | 2000 |
|---|---|---|---|---|---|---|---|---|
| Ks Linse (τ₁ = 12°, p = 1.3), kubisch | 0.0448 | 0.0212 | 0.0152 | **0.0149** | 0.0153 | 0.0166 | 0.0187 | **0.0204** |
| Ks Linse, Kneser | 0.0443 | 0.0210 | 0.0151 | 0.0148 | 0.0152 | 0.0164 | 0.0184 | 0.0201 |
| Vs Linse (τ₁ = 5°, p = 1.15), kubisch | 0.0380 | 0.0127 | 0.0039 | 0.0026 | 0.0020 | (läuft) | | |
| Keil c = 0.5 (τ₁ = 1.9°, p = 1), kubisch | 0.0372 | 0.0120 | 0.0029 | 0.0013 | 0.0008 | 0.0005 | 0.0004 | **0.0004** |
| Keil c = 0.5, Kneser | 0.0367 | 0.0117 | 0.0028 | 0.0013 | 0.0007 | 0.0005 | 0.0004 | 0.0004 |
Zusammen mit VE4 (§3.12, min_W ‖U'W⁻¹‖·‖WL'‖):
| Linse | N = 1000 | 1200 | 1500 | 2000 |
|---|---|---|---|---|
| Ks Linse kubisch / Kneser | 1.490 / 1.482 | | 2.084 | **2.975 / 2.925** |
| Vs Linse kubisch | | 1.329 | | (läuft) |
| Keil c = 0.5 kubisch / Kneser | 1.332 / 1.339 | | 1.333 | **1.332 / 1.338** |
**Befund (numerisch, entschieden für Ks Linse und den Keil):** Ks Vorhersage trifft. Mit der tangentialen Linse hat τ_N ein Minimum bei N ≈ 1000 und
steigt dann wieder (0.0149 → 0.0204 bis 2000, etwa wie N^{0.45}); zugleich wächst das VE4-Produkt (1.33 → 2.98) — die Eck-Köpfe werden kohärent,
der führende Singulärvektor wandert ins Eckband. ‖Ĵ_N⁻¹‖ selbst bleibt bis 2000 bei 1.2512 (das Wachstum der Köpfe ist noch zu klein, um die Norm zu
heben), aber eine gleichmässige Schranke für alle N ist mit dieser Linse nicht zu erwarten. **Mit dem Keil (Ecke an L, Steigung 0.5) fällt τ_N
monoton und sättigt bei ≈ 4e-4, und das VE4-Produkt ist von N = 1000 bis 2000 konstant 1.332 (kubisch = Kneser).** Vs Linse (tangential, aber steil,
p = 1.15) liegt dazwischen: τ fällt bis 1200 noch (0.0020), langsamer als beim Keil — vermutlich späteres Minimum.
**Antwort auf die Frage der Leitung:** Eine feste Linse kann (H1)/(S_γ) für alle N tragen, wenn sie an L eine **Ecke** hat (Keil, s(τ) ≈ cτ mit
0 < c < 0.79 wegen f(ℬ) ⊂ D); tangentiale Linsen (s ~ τ^p, p > 1) verlieren gegen die Spirale. Preis: Für Satz V2 wird das Vergleichsgebiet mit Keil
grösser (α = 0.61 statt 0.76 bei c = 0.5, Schlitz-Schnitt der Kontur jenseits u = 0.2) — die strengen Green-Konstanten müssen für den Keil neu
gerechnet werden; ein flacherer Keil (c = 0.25) läuft als Kompromiss.

### 3.14 Test T1 (Auftrag Leitung 08:58, Ks Warnung §22.6): konvergiert Π_N U auf ∂D nahe L bei schneller Spirale? — `v45_t1_modell.py`
**Aufbau.** Kubisches Modell f = x + (x² − 2Re L·x + ρ²)(α + βx) mit f(L) = L, f'(L) = μ frei (|μ| = 2, arg μ = 38.3°: Drehung Im a = ln|μ|/2π = 6.3°/e).
Statt Daten g^σ (Ks k44: deren Lösung J_∞⁻¹g^σ ist nicht die gesuchte Funktion; dort e^{0.18N} überall) **exakte Lösung U = g^s**, g = (1 − x/L)(1 − x/L̄),
Daten W = J_∞U = U∘f + (f'∘f)U als Taylor-Reihe (Miller-Rekursion für (g∘f)^s). Gemessen wird E_N(z) = |Π_N U − U|(z), Π_N = J_N⁻¹T_NJ_∞ — genau die
Abschnittsprojektion von Ks Kette; die Linse ändert die Lösung nicht (nur die Norm), sie spielt für E_N keine Rolle. Vorlauf (nur Aufbau, N ≤ 200):
|μ| = 2 konvergiert mit exakter Lösung überall (0.9L: 1.5e-7 → 7.3e-9; ∂D bei τ = −0.03: 1.1e-5 → 1.5e-6), die globale Schleife ist also hier zahm.
**Vorhersage (blind, 12.10. 09:05, vor allen Läufen mit N > 200):**
- Ks Kette sagt: sobald die primale Γout-Hälfte im Wachstumssektor von F_P liegt (bei |μ| = 2 nach k43 ab n ≈ 450 bei A = 16, +1.3 Dekaden bei n = 3000),
  hört E_N auf der Γout-Seite (τ < 0) auf zu fallen und wächst schliesslich, während die Γin-Seite und Innenpunkte weiter konvergieren.
- **V sagt voraus:** im erreichbaren Bereich N ≤ 1600 fällt E_N an ALLEN Punkten weiter; für s = 3.5 um ≥ Faktor 10 von N = 400 bis 1600, für s = 1.5 um
  ≥ Faktor 2; die Γout-Seite fällt langsamer als die Γin-Seite (Verhältnis Γout/Γin bei |τ| = 0.03–0.08 wächst, aber um < Faktor 10 von 400 bis 1600).
  Grund: bei festem z wächst A = |z − L|/r_n wie n^{0.89}; Ks Kopfverhältnis steigt bis n = 3000 nur um ≈ 1.4 Dekaden (≈ n^{1.6}), der Eckfaktor
  n^{−1/2}r_n^{s+1} fällt wie n^{−0.5−0.89(s+1)} (s = 3.5: n^{−4.5}; s = 1.5: n^{−2.7}). Erst bei s ≲ 0.5 (n^{−1.8}) wäre im erreichbaren Bereich ein
  Abflachen oder Anstieg zu erwarten — deshalb wird s = 0.5 als empfindliche Variante mitgerechnet; dort sage ich ein Abflachen der Γout-Seite voraus
  (lokaler Abfallexponent → 0), aber keinen Anstieg um mehr als Faktor 2 bis N = 1600.
- Entscheidung: steigt E_N(τ < 0) von N = 800 auf 1600 um mehr als Faktor 2 (bei fallender Γin-Seite), ist Ks Warnung im Modell bestätigt.

- **Nachtrag vor den Läufen (09:34):** das kubische |μ| = 2-Modell (M0, arg μ = 38.3°) ist pathologisch: f(Γin) verlässt D (max|f(Γin)|/ρ = 1.48), die J̃_N
  sind für einzelne N fast singulär (Spitzen bei N = 710, 730, 800, 820, 1100 um Faktor 10–300 an ALLEN Punkten; Präzision geprüft) — Ks „globale Schleife“.
  Gut konditionierte Modelle mit f(Γin) ⊂ D: M2 = kubisch |μ| = 2, arg μ = 90°; **M3 = quintisch f = x + q(x)(α + βx − 0.12x² − 0.01x³)**, |μ| = 2,
  arg μ = 38.3° (gleiche Eckgeometrie wie M0), max|f(Γin)|/ρ = 0.988, max|f|/ρ = 1.34, übrige Urbilder von L bei ≥ 1.46ρ; Scan N = 300…500: glatt.
  Kopf-Test (v47) in M3: Γout-Hälfte driftet ins Wachstum (log10 max|δ_n/δ_n(L)| bei A ∈ [64, 256): −6.6 / −3.1 / −0.8 für n = 200/400/800).
  **Die Vorhersage oben gilt unverändert für M3 (s = 0.5, 1.5) bis N = 1600; zusätzlich für N = 2400: kein Anstieg von E_N auf der Γout-Seite gegenüber
  N = 1200 um mehr als Faktor 2.** Ks Kette sagt in M3 dagegen einen Anstieg im Band τ ∈ [−0.6, −0.1] (dort liegt A ∈ [64, 256)) ab N ≈ 1600 voraus.

- **Ergebnisse (10:05; Läufe sg/a/b/t1/E_*.log, kopf_*.log; v47 = primaler Kopf längs ∂D):**
  **M3** (gut konditioniert, Ecke wie M0): Kopf-Verhältnis log10 max|δ_n/δ_n(L)| auf der Γout-Hälfte bei A ∈ [16, 64) / [64, 256): −0.04/−6.6 (200),
  +0.06/−3.1 (400), +0.12/−0.8 (800), +0.53/+0.32 (1200), +0.42/+0.33 (1600) — die Γout-Hälfte IST im Wachstum, aber beschränkt (Faktor ≈ 3, kein
  F_P-artiges Weiterwachsen), log10‖c_n‖/|δ_n(L)| = −0.60/−0.28/−0.40 (800/1200/1600). E_N (s = 1.5), Bandmaxima auf ∂D:
  | N | 200 | 400 | 800 | 1200 | 1600 | 2000 | 2400 |
  |---|---|---|---|---|---|---|---|
  | Γout τ ∈ [−0.6, −0.1] | 3.0e-5 | 5.5e-6 | 9.8e-7 | 3.4e-7 | 1.8e-7 | 1.1e-7 | 7.6e-8 |
  | Γout τ ∈ [−0.1, −0.02] | 1.2e-4 | 2.5e-5 | 5.3e-6 | 2.1e-6 | 1.2e-6 | 7.4e-7 | 5.3e-7 |
  | sup_∂D (bei L) | 3.3e-4 | 1.1e-4 | 4.0e-5 | 2.1e-5 | 1.4e-5 | 9.6e-6 | 7.2e-6 |
  | Γin τ ∈ [0.1, 0.6] | 1.7e-5 | 2.0e-6 | 1.5e-7 | 2.4e-8 | 5.0e-9 | 1.2e-9 | 3.2e-10 |
  Γout-Band fällt gleichmässig ~N^{−2.1} (kein Abflachen), s = 0.5 ebenso fallend (sup 3.3e-2 → 9.7e-3, Band [−0.6, −0.1] 6.5e-4 → 5.3e-5).
  **Vorhersage: Abfall ✓ (s = 1.5: Faktor 8–30 von 400 bis 1600; s = 0.5: Faktor 2.5–7; kein Anstieg bis 2400 ✓); Verhältnis Γout/Γin < 10 ✗** (wächst auf
  ≈ 240, weil die Γin-Seite viel schneller konvergiert).
  **M2** (kubisch, arg μ = 90°): Kopf auf Γout bleibt fallend (−0.9 bei A ∈ [64, 256), sättigt), E_N fällt monoton (s = 1.5: sup 1.4e-3 → 2.6e-4, N = 400 … 2000).
  **M0** (pathologisch): Kopf auf Γout +0.86/+1.10 (A ∈ [4, 64), n = 1600/2400); E_N auf der Γout-Seite stagniert ab N ≈ 1000 bei 5e-5 … 1e-4 (2000: 9.6e-5,
  2200: 8.9e-5; 2400: 8e-4 bei τ = −0.015; **2800: 5.5e-4**, Maximum bei τ ≈ −0.014, |x − L| ≈ 0.02), Γin konvergiert weiter (2800: ≤ 1.9e-5 nahe L, 1e-8 weiter
  weg) — bei 2400/2800 sind Γin und Innenpunkte NICHT mit hochgegangen, es ist also kein globaler Resonanz-Ausreisser, sondern Wachstum auf der Γout-Seite.
  **M4** (quintisch |μ| = 2.6, arg μ = 38.3°, Drehung 8.7°/e, g = −0.16; leichte Resonanz bei N = 400): Kopf auf Γout +0.27 (400) → +1.24 (800) → +0.72 (1200) →
  +0.95/+1.02 (1600), ‖c_n‖/|δ_n(L)| schwankt um 1; E_N (s = 1.5) Band [−0.6, −0.1]: 2.1e-6 (1000), 1.65e-6 (1200), 1.16e-6 (1600), **1.26e-6 (2000)**,
  Γin-Band fällt 3.9e-8 → 7e-10. Weiter (s = 1.5): Γout-Band 1.56e-6 (1800), 1.36e-6 (2200), **2.58e-6 (2400)**; Profil längs ∂D: das Maximum sitzt in einer
  festen Zone τ ≈ −0.09 … −0.12, |x − L| ≈ 0.13–0.17 (A ≈ 100–130 bei N = 2400, dort Kopf-Verhältnis +1): Wert 1.8e-6 (1000), 1.4e-6 (1200), 1.2e-6 (1600),
  1.6e-6 (1800), 1.3e-6 (2000), 1.4e-6 (2200), 2.3e-6 (2400) — **keine Konvergenz in dieser Zone seit N ≈ 1000**, während näher an L (τ = −0.02: 1.4e-6 → 4.0e-7),
  in L (sup 1.3e-5 → 7.0e-6) und auf Γin (5e-9 → 5e-10) weiter konvergiert. s = 0.5: Γout-Band 6.8e-4 (1000) → 5.4e-4 (1600), fällt kaum. M3-Kopf bei
  n = 2400: +0.51/+0.50 (A ∈ [16, 64)/[64, 256)), also leicht wieder steigend. Läufe M4 N = 2800, 3200 und Kopf n = 2400 laufen (t1/E_M4_s1.5_e/f.log, kopf_M4b.log).
  **Zwischenstand T1:** (1) Das Glied (a) der Kette ist in sauberen Modellen REAL: die primale Γout-Hälfte dreht ins Wachstum (M3 ab n ≈ 1200 um ≈ +0.5 Dekaden,
  M4 ab n ≈ 800 um ≈ +1), aber bis n = 2400 beschränkt/schwankend, nicht F_P-artig super-polynomial. (2) Wirkung auf Π_N U: in M3 keine (alle Bänder fallen bis
  2400, Γout ~N^{−2.1}); in M4 Stillstand in einer Γout-Zone fester Entfernung von L seit N ≈ 1000, erster Anstieg (×1.9) bei 2400 — das ist das von K
  vorhergesagte Muster im Ansatz; ob es Wachstum wird, entscheiden N = 2800/3200. (3) M0 (Ks Modell) ist global fast-singulär und taugt nicht als Beleg.
  **Nachtrag 10:52:** M4, N = 2800 (s = 1.5): Zone wieder 1.55e-6 (τ = −0.097) — die Zonenwerte 1.8/1.4/1.2/1.6/1.3/1.4/2.3/1.55 (×1e-6, N = 1000…2800)
  bilden ein **Plateau ohne Konvergenz** (2400 ist ein Ausreisser nach oben), L und Γin konvergieren weiter (sup 5.5e-6, Γin-Band 4.6e-11). M4 s = 0.5:
  Γout-Band 6.8e-4 (1000) → 5.4e-4 (1600) → **3.6e-3 (2400)**, Γin 1.1e-5 → 5.5e-7. M4-Kopf n = 2400: +1.12/+1.48 (A ∈ [16, 64)/[64, 256)), ‖c_n‖/|δ_n(L)| = 10^{+0.65}.
  **Entscheidung T1 (numerisch, Modelle):** Ks Mechanismus ist real — in sauberen Modellen dreht die primale Γout-Hälfte ins Wachstum, und sobald das
  Kopf-Verhältnis dort ≈ 1 Dekade erreicht (M4), hört Π_N U − U auf, in einer FESTEN Zone der Γout-Seite (|x − L| ≈ 0.13–0.15 in M4, ≈ 0.02 in M0) zu fallen
  (Plateau bzw. Anstieg), während L selbst, die Γin-Seite und das Innere weiter konvergieren. Bei schwächerem Kopfwachstum (M3, +0.5) bis N = 2400 keine Wirkung.
  Meine blinde Vorhersage (kein Anstieg bis 1600; M3 bis 2400) hat gehalten — der Effekt setzt erst danach bzw. erst bei schnellerer Spirale ein.
  Folgerung: „gleichmässig auf D̄“ ist im Modell bei grossem N fraglich; die Zone liegt in festem Abstand r_* von L (modellabhängig), d. h. Ks K_{r₀}-Fassung
  braucht r₀ > r_* — für Kneser nach Ks Profilrechnung r_* ~ 1e-6 (§23.3), hier nicht prüfbar.
  **Nachtrag 11:53 — M4, N = 3200:** Zone 1.59e-6 bei τ = −0.096 (Plateau 1.2–2.3e-6 für N = 1000 … 3200, kein weiteres Wachstum, keine Konvergenz); in L 4.4e-6,
  näher an L (τ = −0.02) 4.5e-7, Γin-Band 1.4e-11 — weiter fallend.



### 3.15 Sprungschranke streng aus Ps oberer Struktur (Auftrag Leitung 08:40, für V2 jenseits u = 0.2 und Xs u ≤ 0.451) — `v44_sprung_zert.py`, `v44_fein.log`, `v44_sprung.json`
**Satz (streng bis auf die Zweigwahl, die über Ps slog-Krawczyk festgelegt ist):** |Δf_K(L(1+u))| ≤ C_J·u^{5.4469} für 0 < u ≤ 0.48 mit **C_J = 1.0378**
(5.4469 < σ = 1 + Re(2πi/L) = 5.446951, also auch für u → 0 richtig; die bisherige Form u^{5.447} ist für u → 0 um u^{−5e-5} zu optimistisch).
Herleitung: sexp_K = P(e^{L(z+θ_K)}) für Im z ≥ 0.36 (p32a), W = P⁻¹(x) ⇒ f_K(x) = P(√L·W·e^{Lδ(z)}), δ = θ_K(z+½) − θ_K(z) = −2Σ_{ungerade}κ_nqⁿ;
die Ufer von S_L sind arg W = θ + 4π (Γin-Seite) bzw. θ + 6π (Γout-Seite); Zweig k₀ = 2 an x₁ = 1.2e^{iθ} mit Ps slog-Krawczyk (Radius 1e-30) eindeutig
(k₀ = 1, 3, 4: ausserhalb). Δf = P(√LWe^{Lδ_in}) − P(√LWe^{Lδ_out}), Mittelwertsatz mit gemeinsamem W. Teil 1 (u ≤ 0.01) analytisch: Rouché auf D(0, 0.02)
(κ₂ ≤ 0.0067), |q_in| ≤ 3.2e-9, Halbebene Im z ≥ 0.45 invariant, Kontraktion 0.036 ⇒ K₀ = 1.0378; Teil 2 (0.01 ≤ u ≤ 0.48) 7900 Intervalle (Faktor
1.0005): Quotient ≤ 0.97 (u = 0.016), fällt auf 0.74 (u = 0.3); Im z_in ≥ 0.515 bei u = 0.48. Punktwerte (Kontrolle): Quotient 0.941 (u = 0.001),
0.875 (0.067), 0.760 (0.2), 0.618 (0.4) — As numerisches C_J ≤ 0.95 ist richtig. Gilt bis u ≈ 0.48 (dann Im z → y₁), deckt Xs u ≤ 0.451 ab.
Mit C_J = 1.0378 (statt 0.95 [num.]) gibt v41 für Ks Linse **B ≤ 16.4 ab N₀ = 400** — damit ist Satz V2 für Ks Linse streng bis auf Xs Lemmata und (Fam).

### 3.16 Satz V2 mit Keil c = 0.25 — Stand beim Zurückstellen (Leitung 09:30; K 09:19: K15 braucht V2/Keil nicht)
Keil s = S(1 − e^{−τ/τ₁}), S = 0.0166, τ₁ = 3.8° (φ = 14.0°, mit K abgestimmt). Streng: Lune c₁ = −0.40 (α = 0.693, Linse ⊂ D₁ inkl. Ordnung 2 an L,
v39), Geometrie (f(ℬ) ⊂ D, Schlitz frei, M_γ ≤ 2.399, v32), Krümmung der Linse ≤ 3.22 (v40, am Knick bei L) ⇒ Lip ≤ 3.22 (innere Scheiben).
Bilanz v41 (C_J = 1.0378): B_max(400) = 0.25 mit Lip 3.22, aber **2.69 mit Lip 0.86** bzw. 2.34 mit Lip 1.0. Der Knick ist nur für die grobe Schranke
schädlich: Vergleich mit D̄ ∪ (innere Scheibe am Knick) (Lune-Green-Funktion explizit) gibt numerisch Lip ≈ 0.875 — strenge Fassung nicht ausgeführt.
VE4/τ_N mit diesem Keil bei N = 2000: min_W-Produkt 1.3323 (Bessel-Faktoren ‖U'W⁻¹‖ = 1.735, ‖WL'‖ = 0.768), τ_N = 0.0017/0.0012/0.00086/0.00071/0.00064
(N = 1000/1200/1500/1800/2000) — fällt; Vs Linse (5°, 1.15): 1.3327, τ_N 0.0026 → 0.0014. Bessel-Schranken (= Quadrate der Faktoren) bei optimalem W: kubisch N = 500/1000/2000: ‖U'W⁻¹‖ = 1.747/1.752/1.735, ‖WL'‖ = 0.763/0.760/0.768 (Produkt 1.3325/1.3323/1.3323); Kneser N = 500/1000: 1.732/1.751 bzw. 0.773/0.765 (1.3384/1.3389). Läufe N = 1500 und Kneser 2000 abgebrochen (Bessel nach K 09:19 nicht mehr gebraucht).

### 3.18 Ks E2: Eckbedingung (EB_N) an den echten p_N (k_kern §24.5) — `v49_E2_EB.py`, `v49_E2.log`, `v49_E2.json`
f_K nahe L auch ausserhalb D̄ aus Ps oberer Struktur (Zweig k₀ = 2, Kontrolle gegen die P-V-Taylorreihe in D: alle 15 Stellen gleich); p_N von A (N ≤ 1000) und
v33 (1500, 2000). Messfalle: acb-Horner nahe |x| = ρ verliert je Schritt den Faktor (|Re x| + |Im x|)/|x| ≈ 1.2 an Radius (bei N = 2000: 1e158) — Auswertung
am Mittelpunkt mit 900 Bit. Kegelscheibe |x − L| ≤ 16r_N (r_N = ρN^{−0.894}), Polargitter 14 × 120:
| N | r_N | Teil in D̄ | aussen, ≤ 4r_N | aussen, 4–16r_N (= η_N) | η_N/r_N |
|---|---|---|---|---|---|
| 200 | 1.2e-2 | 1.9e-11 | 3.3e-9 | 1.9e-3 | 0.15 |
| 400 | 6.5e-3 | 6.5e-13 | 1.1e-10 | 1.1e-4 | 1.7e-2 |
| 800 | 3.5e-3 | 2.3e-14 | 3.8e-12 | 5.1e-6 | 1.5e-3 |
| 1000 | 2.9e-3 | 7.7e-15 | 1.2e-12 | 1.9e-6 | 6.5e-4 |
| 1500 | 2.0e-3 | 1.1e-15 | 1.6e-13 | 2.8e-7 | 1.4e-4 |
| 2000 | 1.5e-3 | 2.7e-16 | 3.9e-14 | 7.4e-8 | 4.8e-5 |
Teil in D̄ ≈ 2·N^{−4.8} (wie vorhergesagt); aussen am Rand der Kegelscheibe η_N ~ N^{−4.65}; Linsen- und Γout-Seite gleich gross (z. B. N = 2000: 7.4e-8 / 6.5e-8),
der Schlitzsektor (±3°) nicht grösser. **(EB_N) gilt mit η_N/r_N ≤ 1.5e-3 ab N = 800, fallend; η_N·ln N ≤ 6e-7 (N = 2000).**
Ruhefenster-Strahlen (θ_phys = u_A − arg Λ₀ − 1.451°·ln N, arg Λ₀ = −71.46°): Fenster 2 (u_A = −115/−95/−75°, θ_phys ≈ 306–347°, in D̄) überall ≈ N^{−4.8} und
nach aussen fallend. Fenster 1 (u_A = 25/45/65/85°, θ_phys ≈ 85–147°, ausserhalb D zwischen S_L und Γin): auf [r_N, 16r_N] ≤ 9e-7 (N = 1000) bzw. 4.8e-8
(N = 2000); WEITER aussen divergiert p_N auf den schlitznahen Strahlen (u_A = 25°, 9° neben S_L: r = 0.05 ⇒ 1.2e-5 / 6e-2 / 7e2 für N = 1000/1500/2000;
u_A = 85°: 2e-10 bei N = 2000) — dort liegt der Rand des Überkonvergenzgebiets. Für (EB_N) also ε ≲ 0.02 auf dem schlitznahen Teil von Fenster 1 wählen.

### 3.17 (R_L), Eckprodukt und Ks E3 bis N = 2000 (Auftrag Leitung 09:30, k_kern §24) — `v48_RL.py`, `v50_E3_norm.py`, Logs sg/a/b/v48_*, v50_*
**Methode (v48):** EINE Vorwärts-Elimination von J̃_N (LU ohne Pivot, Mittelpunkt, Rang-1-Schritte in C) liefert für ALLE n ≤ N zugleich δ_n(L) = w_n
(w = ε_LᵀU⁻¹ per Rekursion), Λ_n(L) = y_n (y = L⁻¹ẽ_L), P_n = w_ny_n, K_n(L,L) = Σ_{k≤n}w_ky_k. Kontrollen: Ks k46-Werte (n = 200/300/400) auf alle Stellen;
δ_n(L) = ε_Lᵀ(J̃_n⁻¹e_n) aus eigener Lösung (v50) = w_n; N = 1000 mit 1100 Bit und N = 2000 mit 1500 Bit stimmen überein (x-Rahmen gut konditioniert,
schon 300 Bit reichen bei N = 400). Kubisch, N = 2000 in 920 s.
| n | K_n(L,L) (kubisch) | |K_n|/n^{Re a} | P_n | C_n = P_n·n·d_n/(aL) |
|---|---|---|---|---|
| 100 | 36.873 − 4.894i | 0.60720 | 0.32620 − 0.03465i | −0.01029 − 0.59894i |
| 400 | 126.465 − 12.487i | 0.60105 | 0.28212 − 0.01993i | −0.01016 − 0.59845i |
| 1000 | 286.664 − 21.680i | 0.59959 | 0.25622 − 0.01213i | −0.01013 − 0.59836i |
| 1500 | 411.888 − 26.922i | 0.59923 | 0.24550 − 0.00909i | −0.01013 − 0.59834i |
| 2000 | 532.698 − 30.937i | 0.59904 | 0.23815 − 0.00708i | −0.01013 − 0.59833i |
⇒ **(R_L) numerisch bis 2000: |ẽ_L*J̃_n⁻¹ẽ_L| = 0.599·n^{0.894}·(1 + O(1/n))**; C konvergiert (5 Stellen ab n = 1000), K_n − CL/d_n = K₀ ≈ 0.70 + 0.05i bei n = 2000
(passt zu Ks K₀∞ = 0.70 ± 0.05). **Kneser (N = 2000, 1500 Bit, 1509 s):** K_2000 = 532.759 − 31.084i, |K_n|/n^{Re a} = 0.61446 (50) → 0.59912 (2000), C_n = −0.01066 − 0.59967i (50) → −0.01029 − 0.59840i (2000) — gleiches Gesetz, |C| = 0.5985 wie kubisch (Re C weicht in der 4. Stelle ab).
**E3:** Q_n = ‖c_n‖/(|δ_n(L)||d_n|^{1/2}) (v50, eigene Lösungen) kubisch / Kneser: 1.694/1.697 (250), 1.482/1.484 (500), 1.425/1.426 (750), 1.397/1.398 (1000),
1.380/1.380 (1250), 1.368/1.369 (1500), 1.3602/1.3605 (1750), 1.3541/1.3544 (2000). B_n = √n‖c_n‖|Λ_n(L)| = √n·Q_n|d_n|^{1/2}|P_n| (Identität, exakt): kubisch 0.792 (250), 0.668 (500),
0.607 (1000), 0.582 (1500), 0.567 (2000), Kneser 0.794 … 0.567; **B_n/(Q_n·n^{−0.053}) = 0.6267 → 0.6262 (kubisch), 0.6268 → 0.6263 (Kneser), n = 250 … 2000** — Ks Gesetz B_n ∝ Q_n·n^{−0.053} hält auf 4 Stellen.
Q_n fällt wie Q_∞ + c·n^{−0.89} (Differenzenquotienten passen zu Exponent Re a), Q_∞ ≈ 1.307 (Fit 1750/2000, c = 42) ⇒ B_n ≈ 0.818·n^{−0.053} → 0 langsam. Kneser = kubisch auf 0.3 %.

### 3.19 Rate s_E von dist_E(f_K, P_N) an der Spiralecke (Auftrag Leitung 10:55, Ks Satz K22) — `v52_sE.py`, Logs sE/
E = (D̄ ∖ U_{r₀}(L, L̄)) ∪ int𝒞 ∪ f_K(int𝒞) mit 𝒞 aus Γout-Bogen, Σ₂, Σ₁, Verbindungsbögen bei r₁ und Linse 1.03ρ (shapely-Vereinigung); f_K auf ∂E: Taylor (A), Ps
obere Struktur nahe L (Übergangskontrolle zur Formel Log f(eᵗ) auf der Linse: alle 10 Stellen gleich), auf f(𝒞) exakt eᵗ. Diskrete Bestapproximation (Arnoldi-
Vandermonde, Lawson 40 Schritte) auf ≈ 2300 Randpunkten, N = 20 … 300; Fehlerpegel 4e-6 → 8e-11.
| (u_A1, u_A2), r₁ | Aussenwinkel bei L (r < r₁) | Vorhersage 5.447α_E/(1+T²) | lokale Exponenten N = 60 → 80 → … → 300 |
|---|---|---|---|
| (80°, −90°), 0.2 | 131.7° | 3.98 | 3.94 / 4.02 / 3.87 / 4.10 / 3.68 / 4.00 |
| (60°, −90°), 0.2 | 111.7° | 3.38 | 3.50 / 3.31 / 3.37 / 3.26 / 3.32 / 3.27 |
| (80°, −90°), 0.05 | 131.7°, jenseits r₁ ≈ 97° | — | 3.31 / 3.28 / 3.02 / 3.07 / 3.19 / 2.90 |
| (60°, −90°), 0.05 | 111.7°, jenseits r₁ ≈ 97° | — | 3.32 / 3.28 / 3.02 / 3.06 / 3.18 / 2.61 |
⇒ **Ks Satz K22 (s_E = 5.417·α_E) bestätigt: 3.96 bzw. 3.36 gemessen auf ±0.1.** Mit r₁ = 0.05 (Ks Kontur) bestimmt bei N ≤ 300 der Verbindungsbogen zur Linse
(Linsenpunkt 111.6° von L aus, Lücke dort ≈ 97° ⇒ 2.94) die Rate — asymptotisch zählt nur der Winkel unterhalb r₁, für die Konstanten bei N₀ = 400 aber die ganze
Geometrie: **Empfehlung r₁ ≈ 0.2 oder Verbindung im selben Winkel.** Budget γ ≈ 0.89 < s_E in allen Fällen mit Reserve ≥ 2.

### 3.20 F3: sup_𝒞 |p_N − f_K| an den echten p_N bis 2000 (HAUPTSATZ Fassung 2 §2.4, K 11:05/11:16) — `v53_F3.py`, `v53_F3.log`
| N | Γout-Bogen | Σ₂ (in D) | Σ₁ (aussen, max bei r = 0.05) | Linse 0.03 (max bei arg ≈ 78°, neben S_L) | Steg (Linse 150–180°) | Verbindung r₁ |
|---|---|---|---|---|---|---|
| 200 | 5.6e-12 | 1.9e-11 | 9.3e-11 | 1.7e-9 | 4.4e-12 | 1.7e-9 |
| 600 | 9.4e-16 | 9.1e-14 | 3.9e-12 | 5.3e-9 | 2.0e-11 | 5.3e-9 |
| 1000 | 4.6e-16 | 7.8e-15 | 3.8e-12 | 4.3e-7 | 5.0e-10 | 4.3e-7 |
| 1500 | 4.6e-16 | 1.1e-15 | 1.5e-11 | 5.7e-4 | 2.2e-7 | 4.6e-4 |
| 2000 | 4.6e-16 | 4.5e-16 | 1.2e-10 | **1.75** | 1.1e-4 | 1.07 |
⇒ **F3 ist auf der Linse verletzt:** p_N − f_K wächst dort ab N ≈ 600 wie e^{+0.015N} (1000 → 2000: Faktor 4e6), am stärksten am schlitznahen Linsenanfang, am Steg
ebenfalls (≈ e^{+0.012N}); Γout-Bogen und Σ₂ konvergieren bis an die Rundungsgrenze, Σ₁ (aussen) wächst langsam ab N ≈ 1000. Das bestätigt Ks Messung (primale Köpfe
auf der Linse s = 0.03: +0.016/Schritt) am echten p_N: **eine Datenbedingung sup_𝒞|P − f| ≤ δ₀ ist mit dieser Linse für N ≳ 1500 falsch.**

### 3.21 Vortest (b), angefangen und an K abgegeben (Leitung 11:40) — `v54_vortest_b.py`, vb/vb_n200.log
Nur n = 200 gerechnet. Bausteine geprüft: A = c_s/C = −0.11151 − 1.03750i stimmt mit dem Mellin-Residuum Γ(β+iα)/(Γ(β)λ̃) auf alle 10 Stellen (Ks C bestätigt).
V(ζ) im Potenzsektor (arg ζ ∈ [195°, 345°], |ζ| ≥ 30) = A/ζ − Σ_m Γ(s′_m)(−1)^m(λ̃ζ)^{−s′_m}/(m!·iα(1 − s′_m)Γ(β)), s′_m = −(β+m)/(iα) (Re s′₀ ≈ 3.45), Abweichung zum
Integral ≤ 1e-12; bei arg ζ ≈ 180° ein subdominanter Exponentialteil (2e-3 bei |ζ| = 30), unter 175° Wachstum (|V| = 1.5e12 bei |ζ| = 100, 165°). Auf |x| = 0.99ρ liegt
arg ζ ∈ [182°, 348°] — dort kein Wachstum, |(C/d_n)V| ≈ |c_s/(x − L)|.
Ergebnis n = 200 (Gewicht w_n = |Λ_n/Λ_n(L)|): ũ = 0: max R_n = 220 (Linse), 150 (Σ₁, Σ₂); **ũ = T_n[s] (s aus u_800): R_n ≈ 5e41, Komposit ≈ 6e40** — das x-Rahmen-Residuum
b − J̃_nũ wird von den Einträgen von J̃ (~1.6^k) aufgebläht, sobald ũ in den hohen Koeffizienten (k ≈ n) von u_n abweicht; ρ_{m,n} entsprechend bis 2e22 (Maximum bei
m ≈ 0.8n). Die naive Komposit-Parametrix trägt also nicht; ein brauchbarer Kandidat muss die Koeffizienten k ≈ n bis auf e^{−cn} treffen (oder das Residuum in
einem Rahmen messen, in dem J̃ nicht verstärkt). Weiterführung bei K.

### 3.22 E4′: Q_n = max_{z ∈ K_{r₀}, t ∈ 𝒞}|ΔK_n(z, t)| bis N = 2000 (k_kern §28.5, PLAN „Anpassung K 10:48“) — `v51_E4.py`, Logs sg/a/b/v51_*
Eine LU von J̃_N mit mitgeführten Zeilen ε_zᵀU⁻¹ (z ∈ ∂K_{r₀}, 240 Kreis- + 3×40 Bogenpunkte, obere Hälfte) und Spalten L⁻¹ẽ_t (t ∈ 𝒞: Σ₁ 80°, Σ₂ −90° bis r₁ = 0.05, Verbindungsbögen,
Linse 0.03, Γout-Bogen) liefert δ_n(z), Λ_n(t) für alle n ≤ N. Kontrollen: |P_n| = v48 auf 5 Stellen (n ≤ 1800); n = 200/400 wie Ks k49 (R_P 0.82/0.51, R_A Σ = 1.008).
Kubisch, N = 2000 (1300 Bit; **n = 2000 selbst ist unzuverlässig** — |P_2000| = 0.2308 statt 0.2383 aus v48, der letzte Schritt verliert Stellen; Werte bis 1800 glatt):
| n | |P_n| | R_P r₀ = 0.1 / 0.05 / 0.02 | R_A Σ | R_A fern (Lage) | R_A Verbindung | Q_n(0.05) | Q_n/|P_n| | ln|δ_n/δ_n(L)| bei 120°/150°/180° |
|---|---|---|---|---|---|---|---|---|
| 200 | 0.3046 | 0.546 / 0.769 / 0.995 | 1.0079 | 0.178 (Linse 180°) | 0.507 | 0.236 | 0.775 | −2.15 / −2.41 / −2.47 |
| 600 | 0.2709 | 0.118 / 0.279 / 0.633 | 1.0082 | 1.3e-3 (Γout 74°) | 0.046 | 0.076 | 0.281 | −6.44 / −7.16 / −7.34 |
| 1000 | 0.2565 | 0.0272 / 0.109 / 0.391 | 1.0082 | 9.3e-6 | 3.5e-3 | 0.028 | 0.110 | −10.66 / −11.84 / −12.14 |
| 1400 | 0.2475 | 0.0064 / 0.043 / 0.248 | 1.0081 | 1.8e-7 (Linse 78.7°) | 2.5e-4 | 0.011 | 0.044 | −14.86 / −16.51 / −16.93 |
| 1800 | 0.2409 | 0.0015 / 0.017 / 0.159 | 1.0083 | 4.2e-8 | 1.8e-5 | 0.0042 | 0.017 | −19.09 / −21.14 / −21.81 |
| 2000 (1700 Bit) | 0.2383 | 0.00074 / 0.011 / 0.127 | 1.0082 | 2.3e-8 | 4.8e-6 | 0.0026 | 0.011 | −21.48 / −22.98 / −23.78 |
⇒ **Mit z ∈ K_{r₀} OHNE z = L fällt Q_n weit schneller als |P_n|:** Q_n/|P_n| = 1.008·R_P, R_P(0.05) ~ n^{−3.1} (1000 → 1800; Ks Potenzgesetz n^{−3.05}), R_P(0.1) noch schneller,
R_P(0.02) ~ n^{−1.5} (vorasymptotisch, näher an L). R_A = 1.008 kommt immer von den Spiralen nahe L ((A-lok), κ_A = 0); die fernen Teile fallen exponentiell ((FA)), die Verbindungsbögen
ebenso. Ks Vorhersage Q_n ≈ 0.54·n^{−0.106} galt für 𝒦 = K_{r₀} ∪ {L}; ohne L ist Q_n bei 1800 um den Faktor 60 kleiner. Σ_n Q_n bleibt damit beschränkt (γ = 0 statt 0.89), sofern R_P weiter wie eine Potenz < −1 fällt.
**Kneser (N = 2000, 1300 Bit):** bis n = 1800 auf < 1 % gleich (z. B. n = 1800: |P| = 0.24097, R_P = 0.00154/0.0175/0.160, R_A Σ = 1.0083, fern 4.2e-8, Q_n(0.05) = 0.0042;
ln|δ/δ(L)| = −19.00/−21.03/−21.69); n = 2000 wieder unbrauchbar (letzter Schritt). Für n = 2000 bräuchte es P ≈ 1600 Bit.
Primal am Rand (Ks Punkt 4): ln|δ_n/δ_n(L)| bei 120°/150°/180° fällt bis n = 1800 linear (Rate −0.0105/−0.0118/−0.0121) — **kein Kippen in die Potenz n^{−3.05} bis 1800**
(Ks Vorhersage 1500–2500). **Nachtrag 12:25 (Wiederholung mit 1700 Bit, |P_2000| = 0.23825 = v48):** von 1800 auf 2000 wird der Abfall bei 150°/180° flacher (−0.0092 bzw. −0.0099 je Schritt statt −0.0118/−0.0121), bei 120° nicht (−0.012) — möglicher Beginn des Übergangs; v56 bis N = 3000 läuft.

### 3.23 f_K-Geometrie von E (Ks k54 mit f_K statt kubisch; Auftrag Leitung 11:40) — `v55_geometrie_fK.py` (Kopie von k_kern/k54, nur f, f′ ersetzt), `v55_geometrie_fK.log`
f_K: Taylor (A, 4000) in D̄ fern von L, Ps obere Struktur (Zweig k₀ = 2) für |x − L| < 0.35 (konjugiert bei L̄), auf dem Linsenstreifen f(t) = Log f_K(eᵗ) (Zweig an f_K(ρe^{iφ})
angeschlossen); f′ per Taylor-Ableitung bzw. zentraler Differenz. Parameter wie Ks Standard: u₁ = 80°, u₂ = −90°, Linse s = 0.03, r₁ = 0.05, r₀ = 0.02.
| Prüfung | kubisch (K, k54) | f_K (V) |
|---|---|---|
| (1) Γin = {|f| < ρ} auf |x| = ρ | Wechsel bei arg L, arg L̄; max|f|/ρ auf Γout 1.624 | gleich; max|f_K|/ρ 1.624 |
| (2) Univalenz auf int 𝒞 | Windungszahl 1.000 (600 Proben), min|f′| = 0.301 | **1.000, min|f_K′| = 0.336** |
| (3) E zusammenhängend, ohne Löcher | 1 / 1 Komponente | **1 / 1** (Rasteranteile 61413 / 60550 / 59371) |
| (4) Schnittspirale in der Lücke meidet E | ja (r ∈ [1e-6, 1]) | **ja** |
| (6) Aussenöffnung bei L | 132.3° … 131.9° (r = 3e-2 … 1e-6) | **132.3° / 132.1° / 132.0° / 131.9° / 131.9°**, Lücke dreht von (13°, 145°) nach (357°, 129°) |
| (5) Dreipunkt-Verhältnis | 1.07 | **1.07** |
⇒ **E ist auch für f_K eine Quasischeibe mit Spiralecke Θ ≈ 132° (numerisch, Raster 0.01 und Stichproben);** damit greift Ks Satz K22 mit s_E = 5.417·0.732 = 3.96, wie in §3.19
gemessen. Nicht streng: Raster/Stichproben, f_K ausserhalb D über die Funktionalgleichung (nicht zertifiziert eingeschlossen).

### 3.24 (Fam-koh)-Test: Carathéodory–Fejér-Repräsentant von p_N − f_K auf Int γ (Auftrag Leitung 12:00; K §34.2) — `v57_cf.py`, `v58_adj_steg.py`, Logs cf/, st/, sg/a/b/v58_*
CF := min_{P̃ ≡ p_N mod x^{N+2}} sup_{Int γ}|P̃ − f_K| = ‖T_{N+2}(c)‖₂, c_m = [w^m](p_N − f_K)∘k (k: 𝔻 → Int γ, k = ρwe^{F}, Theodorsen; Linse in Ks Form
S(1 − e^{−(τ/12°)^{1.3}}); F₀ = F(0) = harmonisches Mittel der Linsenhöhe). Koeffizienten von p_N − f_K mit 500 Bit (P-V bis 640, A darüber; Fehler ≤ 1e-40).
| S (F₀) | N = 800 | 1000 | 1500 | 2000 | Rate 1500 → 2000 | max |c_m| bei m ≈ |
|---|---|---|---|---|---|---|
| 0.0166 (0.0084) | 1.9e-14 | — | 2.3e-14 | 5.3e-14 | +0.0017 | 0.93N |
| 0.022 (0.0111) | 7.6e-13 | — | 2.9e-11 | 7.3e-10 | +0.0064 | 0.95N |
| 0.026 (0.0131) | 1.3e-11 | — | 6.4e-9 | 9.9e-7 | +0.0101 | 0.95N |
| 0.03 (0.0151) | 2.2e-10 | 2.3e-9 | 1.5e-6 | 1.5e-3 | +0.0137 | 0.95N |
(S = 0.0166, N = 400: 5.1e-13; S = 0.03, N = 200/400/600: 2.4e-11 / 1.5e-11 / 5.4e-11.) Inkohärenter Vergleich |δ̃_N|e^{NF₀} bei N = 2000: 7.4e-13 / 1.7e-10 / 8.9e-9 / 4.7e-7.
⇒ **Der beste Keim-Repräsentant wächst auf der Linse exponentiell, sobald S ≳ 0.015**: Rate ≈ 0.90·(S − 0.0147), kleiner als F₀, das Maximum der w-Koeffizienten im Stegband
m ≈ 0.95N. Polynomial (bis 2000 praktisch flach, unter der inkohärenten Schätzung) ist er nur bei der Standardlinse S = 0.0166. **(Fam-koh)^{r₀} in Ks Form (Repräsentant auf der
Linse beschränkt) ist für S = 0.03 bei N ≈ 2300–2500 verletzt** (CF ~ 0.1–1); Ks inkohärente Schätzung (bis N ≈ 6300) war zu optimistisch, weil das Stegband m ≈ 0.95N grössere
Kopfabweichungen trägt als der oberste Koeffizient.
**Gegenstück: adjungierte Rate am Steg** (v58, kubisch, ln|Λ_n(−ρ(1+S))/Λ_n(L)|, Rate 1600 → 2000): S = 0 / 0.0166 / 0.022 / 0.026 / 0.03: +0.0120 / −0.0038 / −0.0088 / −0.0126 / −0.0163
(Nullstelle bei S ≈ 0.0126 = Ks Schwelle c/0.95). **Summe CF-Rate + adjungierte Steg-Rate: −0.0021 / −0.0024 / −0.0025 / −0.0026 — unabhängig von S, streng negativ.**
Lesart (numerisch, Vermutung): die Linse verschiebt nur Gewicht zwischen Repräsentant und adjungiertem Kopf; das Produkt, das in der Datenkorrektur der Familie auftritt, hat eine
linsenunabhängige Reserve ≈ 0.0024 — Ks e^{ΛN} (§10.2, inkohärent) wird durch die Kohärenz der echten p_N zu e^{−0.0024N}. Vorbehalt: CF misst die Gleichverteilung über ganz γ
(der Extremal-Repräsentant hat konstanten Betrag auf γ, auch an den Spiralen nahe L, wo das adjungierte Gewicht ≈ 1 ist); der genaue Test ist das **gewichtete** CF-Problem mit dem
Aussengewicht |Λ_n/Λ_n(L)| (über eine äussere Funktion Φ mit |Φ| = Gewicht auf γ: ‖T(c ∗ Φ)‖) — nächster Schritt, falls gewünscht.

### 3.25 Primale Randrate bis N = 3000: der Übergang in die Potenz kommt bei n ≈ 2200 (Ks Punkt 4, PLAN 10:48) — `v56_rand.py`, sg/a/b/v56_kubisch_N3000.log
Eine LU von J̃_3000 (kubisch, 2000 Bit) mit den Zeilen ε_zᵀU⁻¹ für z = ρe^{i·120°/150°/180°}; ln|δ_n(z)/δ_n(L)|:
| n | 1000 | 1600 | 2000 | 2200 | 2400 | 2600 | 2800 | 3000 |
|---|---|---|---|---|---|---|---|---|
| 120° | −10.66 | −16.97 | −21.48 | −23.07 | −22.81 | −22.99 | −23.20 | −23.42 |
| 150° | −11.84 | −18.85 | −22.98 | −24.05 | −24.41 | −24.71 | −24.82 | −25.18 |
| 180° | −12.14 | −19.33 | −23.78 | −24.62 | −25.34 | −25.50 | −25.55 | (−27.85) |
⇒ **Exponentialphase (Rate −0.0105 … −0.0121) endet bei n ≈ 2100–2300, danach Potenzabfall** (120°: −0.61 von 2400 auf 3000, Exponent ≈ −2.7; Ks n^{−3.05}); Ks Vorhersage
1500–2500 ✓, Lemma K20 (keine exponentielle Reserve am Rand) damit auch am Rand gestützt. Der letzte Wert bei 180° (n = 3000) ist wie in §3.22 vom letzten Eliminationsschritt
verfälscht. v56 reproduziert v48 (K_n, P_n bis 2000 auf alle Stellen) und v51 (Randwerte bis 2000).

### 3.26 Linsenfenster S ≈ 0.0137 und Kopf×Kopf am Steg bis N = 2000 (Leitung 13:00/13:08; K §39) — `v59_elim.py`, `v59c_koeff.py`, `v59_eval.py`, Server sg/a/b/v59_*
**(a) Kopf×Kopf am Steg, gleich gemessen (v59_elim, kubisch, 3001 Bit, z = t = −ρ(1+S), S = 0 … 0.03):** Raten je Schritt (Steigungen), primal + adjungiert:
| n | 300→600 | 600→1000 | 1000→1500 | 1500→2000 |
|---|---|---|---|---|
| Summe (jedes S ∈ [0.0126, 0.03]) | +0.00025 | +0.00014 | +0.00009 | +0.00006 |
| S = 0.0137: primal / adjungiert | +0.00089/−0.00064 | +0.00100/−0.00086 | +0.00105/−0.00096 | +0.00107/−0.00101 |
- Die Summe fällt wie 1/n: ln|ΔK_n(t,t)/P_n| = 1.497 / 1.572 / 1.626 / 1.669 / 1.700 (n = 300/600/1000/1500/2000), mit |P_n| ~ n^{−0.106} also **|ΔK_n(t,t)| → 1.304**
  (n = 1000/1500/2000: 1.3042/1.3042/1.3042), für alle S ∈ [0.0126, 0.03] gleich auf 0.5 %. Mit m(−1.0137ρ) = 0.7334: |m·ΔK_n(t,t)| = 0.956 = Ks κ(−ρ) ✓ (K §39, KK′).
- Ks frühere +0.003 reproduziert: Ks k52-Grösse h = (ln R)/n ist ein Mittel, keine Steigung (bei n = 600: +0.01484/−0.00090/−0.01341 für S = 0/0.0166/0.03, genau Ks Werte).
  Steigung primal (+0.0162) + Mittel adjungiert (−0.0134) = +0.003. **Gleich gemessen ist die Summe ≈ 0 (fallend) — kein Widerspruch zu meinem −0.0024.** Das −0.0024
  ist CF-Repräsentant (Steigung) + adjungiert (Steigung): Der Repräsentant ist um ≈ 0.0025 je Grad kohärenter als ein Kopf.
- Auf ∂D (S = 0) ist die Summe bei 1500→2000 +0.0007: Die Primale verlässt dort ihre Exponentialphase (§3.25), die Adjungierte noch nicht.
**(b) Ungewichtetes CF im Fenster (v59c, Kneser p_N − f_K, K-Form-Linse τ₁ = 12°, p = 1.3):**
| S | N = 600 | 1000 | 1500 | 2000 | Steigung 1500→2000 |
|---|---|---|---|---|---|
| 0.0126 | 7.1e-14 | 6.1e-15 | 8.8e-16 | 2.2e-16 | −0.0028 |
| 0.0137 | 7.2e-14 | 6.2e-15 | 8.8e-16 | 3.7e-16 | −0.0017 (Potenz −3.0) |
| 0.0147 | 7.2e-14 | 6.2e-15 | 2.0e-15 | 2.0e-15 | 0.0000 |
| 0.0166 | 7.3e-14 | 1.6e-14 | 2.3e-14 | 5.3e-14 | +0.0017 |
⇒ Im Fenster fällt der beste Repräsentant bis 2000 (bei 0.0137 noch wie N^{−3}); der Umschlag liegt bei S ≈ 0.0147, wie aus der Geraden 0.90(S − 0.0147).

**(c) Adjungierte auf der GANZEN Linse (v59_eval: ℓ_n als exakter Vektor, Λ_n(t)/Λ_n(L) an 16384 Punkten von γ in Ballarithmetik, rel. Radius ≤ 1e-80).**
Steigung (ln w_2000 − ln w_1500)/500, Maximum je Band des Winkelabstands zu L bzw. L̄, S = 0.0137, K-Form τ₁ = 12°, p = 1.3 (Kneser = kubisch auf ±0.00006):
| Abstand zu L | 0.5–1° | 2–3° | 3–4° | 5–6° | 7–8° | 8–10° | 10–15° | 15–30° | 30–60° | Steg-Band |
|---|---|---|---|---|---|---|---|---|---|---|
| \|t − L\| | 0.024 | 0.063 | 0.072 | 0.121 | 0.169 | 0.193 | 0.241 | 0.361 | 1.38 | 2.17 |
| Linsenhöhe s | 0.0005 | 0.0018 | 0.0021 | 0.0038 | 0.0054 | 0.0061 | 0.0075 | 0.0101 | 0.0137 | 0.0137 |
| Steigung | +0.00085 | +0.00114 | +0.00113 | +0.00082 | +0.00030 | +0.00001 | −0.00057 | −0.00174 | −0.00163 | −0.00110 |
- **Am Steg und auf dem Hauptteil der Linse hat die Adjungierte Reserve (≤ −0.0006 für |t − L| ≥ 0.24), aber nahe L, wo die K-Form-Linse dünn ist (s ≤ 0.006),
  wächst sie (bis +0.0011, Zone |t − L| ≤ 0.19 ≈ 8°).** Die nötige Höhe steigt von 0 in L an: Steigung 0 bei s ≈ 0.0028 (2°), 0.0038 (3.4°), 0.0050 (5.5°), 0.0061 (8°)
  (aus S = 0.03/0.022/0.0166/0.0137). Bei S = 0.0126 ist auch der Steg nur noch ±0.00005.
- Steilerer Anstieg an L (S = 0.0137, τ₁ = 8°/6°/4°): Wachstumszone |t − L| ≤ 0.083/0.045/0.019 (also innerhalb r₁ = 0.1), aber das ungewichtete CF fällt dann nicht mehr:
  CF(1000/1500/2000) = 6.6e-15/4.3e-15/4.9e-15 (τ₁ = 8°, Steigung +0.0003), 1.9e-14/1.7e-14/2.9e-14 (6°), 7.1e-14/1.0e-13/2.8e-13 (4°, +0.0028).
  ⇒ **Getrennt haben Repräsentant und Adjungierte nirgends in dieser Linsenfamilie beide strikte Reserve auf ganz γ ∖ U_{0.1}(L, L̄);** am nächsten kommt τ₁ ≈ 8° (CF flach).

**(d) Gewichtetes CF (Auftrag Leitung 13:00) — `v59_gewichtet.py`, Logs cfw/v59_cfw_*.log:** CF_w(N) = min_{P̃ ≡ p_N mod x^{N+2}} sup_γ |Λ_N(t)/Λ_N(L)|·|P̃ − f_K|
= ‖T_{N+2}((Φ·c) mod w^{N+2})‖₂, Φ äussere Funktion mit |Φ| = w_N∘k (Φ, 1/Φ analytisch ⇒ Nebenklasse erhalten). Gewicht mit kubischem und mit Kneser-J (gleich auf 0.1 %).
| Linse | N = 200 | 600 | 1000 | 1500 | 2000 | Potenz 1500→2000 |
|---|---|---|---|---|---|---|
| S = 0.0126 … 0.022 (τ₁ = 12°) | 9.02e-12 | 7.53e-14 … 7.62e-14 | 7.29e-15 … 7.39e-15 | 1.13e-15 | **2.93e-16 … 2.97e-16** | −4.68 |
| S = 0.03 (τ₁ = 12°) | 9.09e-12 | 7.70e-14 | 7.46e-15 | 1.15e-15 | 2.99e-16 (Gitter 131072, exakte Faltung) | −4.69 |
| S = 0.0137, τ₁ = 8°/6°/4° | — | — | 7.4e-15/7.4e-15/7.6e-15 | — | 2.96e-16/2.99e-16/3.05e-16 | −4.68 |
- **Das gewichtete CF fällt für jede Linse wie N^{−4.7} (Potenz −4.34 bei 200 → −4.68 bei 2000), mit linsenunabhängigem Wert** — auch bei S = 0.03, wo der ungewichtete
  Repräsentant auf 1.5e-3 wächst. Die Familienannahme in Produktform (sup_γ|Λ_N/Λ_N(L)|·|P̃ − f_K| polynomial fallend) hat damit bis N = 2000 strikte Reserve, numerisch.
- Prüfungen: Re h = ln w auf 4e-14; halbes Randgitter gleich auf 4 Stellen (Ausnahme S = 0.03, N = 2000: mit 16384 Punkten instabil, 1.2e-13 gegen 1.05e-12 —
  mit 131072 Punkten und exakter Faltung (flint) 2.993e-16 kubisch / 2.996e-16 Kneser, halbes Gitter 3.09e-16); Gewicht in Ballarithmetik (rel. Radius ≤ 1e-80).
- Lesart: Die Gewichtung tötet das Steg-Wachstum vollständig (Produkt am Steg linsenunabhängig, Teil (a)); übrig bleibt der gemeinsame Anteil nahe L.
- Nicht gemessen (Team angehalten 14:12): Ks (BU)-Messungen (R_P, R_A auf S₀ = 0.0125 bis 3000; gewichtetes CF mit den Paarungen π_n).

## 4. Dateien
- `ABHAENGIGKEIT.md` — Bausteintabelle, Status, Konstanten, Abhängigkeiten.
- `zert400/` — Zweig-Zertifikat 193 → 400 (Skripte, Daten, später Ergebnisse).
- `m_kopie/` — Kopien von Ms Kneser-Skripten (Pfad zu gemeinsam/ angepasst) für Prüfläufe.

## 5. PRUEFLISTE-Notiz (formale Punkte für p_pruefung/PRUEFLISTE.md; Ps und Ms Ordner ändert V nicht)
- **m26 (M, Zusammensetzer von (a″)_∞):** liest je Teillauf nur die letzte Zahl und prüft keine Prüfzeile — weder „gescheitert: 0“ (m16), noch
  „Analytizität … geprüft“ und „Halbachse … geprüft“ (m24), noch Im H ∈ (0, π) (m25). Ein gescheiterter Teil würde unbemerkt eine falsche Konstante
  liefern. Für die Läufe mit Vs (K2)-Datei hat V alle Prüfzeilen einzeln nachgesehen (§1.8): in Ordnung. Empfehlung an M: Abbruch, wenn eine
  Prüfzeile fehlt oder negativ ist.
- m23: Kommentar zu den Schwanzpaaren falsch (gegenüberliegende flache Bereiche haben S = 1), die Schranke deckt sie trotzdem — ohne Wirkung.
- mk_kneser_arb: Rundung von k·mid(b_k) bei 120 Bit fehlt in Rc1 (< 1e-30) — ohne Wirkung.
- m16: Faktor (1 − e^{−2|σ|π})⁻¹ in der Schwanzschranke fehlt (Wirkung 1e-50) — ohne Wirkung.
- P-V: formale Lücken L1–L4 (§1.1) — ohne Wirkung.
