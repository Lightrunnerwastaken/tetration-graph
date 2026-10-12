# Forschungsteam: Polynom-Wurzel von exp → Kneser (Beweis der Schranke)

## PHASE 3 (11.10. 19:50 – 12.10. 14:15, auf Janis' Wunsch beendet): Stränge K und V, siehe `phase3/README.md`; Ergebnis: `phase3/HAUPTSATZ_KR0.md` (Fassung 3), offener Kern `phase3/KERNPROBLEM.md`, Problemblatt `phase3/OPEN_PROBLEM_2_en.md`; im Graphen Strang kn.halbe_exp

## PAUSE seit 11.10.2026 04:51 (Janis: „langsam stoppen“)

Alle fünf Stränge sind sauber angehalten, es laufen keine Prozesse (lokal und auf aicgolling). Neu starten nur auf Janis' Auftrag.

**Stand:** Hauptsatz (K) + (H₀) + (A) ⇒ p_N → f_K gleichmässig auf D̄ (ARCHITEKTUR §9.9). Bewiesen bzw. zertifiziert:
- Satz P-V: Knesers sexp zur Basis e streng zertifiziert (±1.3e-168), a₀…a₆₄₀ in `gemeinsam/kneser_halbexp_koeff_zert_640.json`; (K3) f_K(L) = L, f_K′(L) = √L (p32a).
- (A) rigoros für N = 140–193 bis auf die Tail-Schranke (K2) jenseits k = 640 (numerisch 4e-15, Reserve ~1e7).
- Bulk (a″)_∞: Ms 𝒯-Weg mit zertifizierten Koeffizienten σ_min(J_∞|P_N) ≥ 0.6553 (wahr 0.7319), fehlt nur (K2).
- Zweig bis N = 193 zertifiziert; N₀ mit beweisbaren Konstanten ≈ 318 (X, KONSTANTEN.md) → Zertifikat bis 400 vorbereitet, nicht gelaufen.

**Offen:** (H₀) — globale starke Asymptotik der Eckmoden bei L (W: GR in der Ebene scheitert, n-te Wurzel reicht nicht, Vermutung W12/W12′,
Schweiss-RH-Aufgabe mit Parametrix in L); Kante (a″)_{αN} für kleines α (ist eine Eckaussage, gehört zu (H₀)); (K2) mit Exponent 6;
Prüfschlange P (Rest), Gegenlesen von P-V Schritte 3–6.

**Wiederaufnahme:** ARCHITEKTUR §9.9 (Dateiliste), BRETT-Abschlusseinträge 04:48–04:51, `w_grenzschicht/WIEDERAUFNAHME.md` und
`OFFENES_PROBLEM.md` / `OPEN_PROBLEM_en.md`, `x_alternativen/zert/STARTANLEITUNG.md` (≈ 2.5 h auf aicgolling), `m_modell_eckmode/BERICHT.md` §17,
`o_operator/BERICHT.md` §34. Experten nur mit Janis' Zustimmung anfragen.

## PHASE 2 (ab 10.10. 23:30) — Dauerbetrieb bis zum vollständigen Beweis oder bis Janis stoppt

**Ausgangspunkt ist ARCHITEKTUR.md** (Phase 1, Stränge A–E): Hauptsatz (K) + (H) + (A) ⇒ p_N → f_K; offen ist nur die lineare Annahme
(H) über die Eckmoden bei L, Kern (a′) Rest-Kontraktion und (E_L) Randprofil der adjungierten Eckmoden. Die Phase-1-Ordner a–e sind
abgeschlossen (nur noch lesen).

Neue Stränge (je ein lange laufender Agent, Ordner in diesem Verzeichnis):
- **P — Prüfung** (`p_pruefung/`): prüft jeden Beweis aus ARCHITEKTUR.md Schritt für Schritt; führt `PRUEFLISTE.md`. Nur von P als
  **verifiziert** markierte Sätze gelten als Fundament. Einziger Strang mit Kneser-Motor (für unabhängige Nachrechnungen).
- **W — Grenzschicht/WKB** (`w_grenzschicht/`): asymptotische Theorie der Eckmoden an der Abschneidekante → (E_L), Rate c ≈ 0.012.
- **O — Operatortheorie** (`o_operator/`): Parametrix für I + K_N, Pseudospektren endlicher Abschnitte → (a′); Riemann–Hilbert-Test
  am Möbius-Modell.
- **M — Modell mit Eckmode** (`m_modell_eckmode/`): einfachstes Modell mit abstossendem drehendem Fixpunkt, dort (H) vollständig beweisen.
- **X — Alternativen & Literatur** (`x_alternativen/`): Literatur zu den neuen Werkzeugen; Beweiswege, die (H) ganz umgehen.

Regeln wie unten (Abschnitt 5), dazu: neue BRETT-Einträge mit `[Strang P|W|O|M|X]`; Ergebnisse, auf die andere bauen, erst nach
P-Prüfung als „bewiesen“ führen; Präzisionswarnung von Strang A beachten (Zufallsfehler werden ≳ 1.17^N verstärkt — Präzision
durchreichen, 60-stellige Dateien nur bis N ≈ 1000).


Stand: 10.10.2026, 19:10. Leitung: Hauptsitzung (Claude). Fünf Stränge A–E, je ein Ordner.

## 1. Das Problem

Für jedes N sei p_N(x) = Σ_{k≤N} a_k^(N) x^k das reelle Polynom mit

    p_N(p_N(x)) = e^x + O(x^{N+1})      (N+1 Gleichungen [x^k] p_N∘p_N = 1/k!, k = 0..N)

auf dem Zweig, der stetig von p_1 = x + 1/2 aus fortgesetzt wird (Newton, Start: p_{N−1} mit a_N = 0).
Code: `experimente/polynom_wurzel.py` → `loesungen(nmax, bits)` gibt {N: [a_0..a_N]} (mpmath), `horner(c, z)`.
Jacobi-Matrix: J[k][j] = [x^k] p^j + ([x^{k−j}] p'∘p falls j ≤ k).

**Vermutung (Janis):** p_N → f_K := T_e(slog_e(x) + 1/2) (Knesers halbe e-Funktion) lokal gleichmässig auf
D = {|x| < |L|}, L = 0.31813150520476413 + 1.3372357014306895 i (Fixpunkt von exp), |L| = 1.3745570107437075.

## 2. Was bewiesen ist (Graph-Knoten in `daten/aussagen/05_kneser_slog.yaml`)

1. **Rückführung** (kn.satz.polynom_halbexp_reduktion): Gilt
   (i) sup_N sup_D |p_N| < ∞, (ii) |p_N| ≤ |L| auf jedem Kompaktum des Korridors V₀ = {z ∈ D : f_K(z) ∈ D} für grosse N,
   (iii) p_N(0) ≤ ρ₀ < |L|, dann p_N → f_K. (Montel + Maximumprinzip + Eindeutigkeit.)
2. **Eindeutigkeit** (kn.satz.halbexp_beschraenkt): Jede auf D holomorphe, reelle Lösung von f∘f = exp (nahe 0), die V₀ in D
   abbildet, ist f_K. Beweis vollständig (Bôcher-Argument für die periodische Störung θ in Kneser-Koordinaten), nur die Form
   von Ω = slog(D) ist numerisch geprüft (waagrecht konvex, Breite ≥ 2.067, Re ≤ 0.3122, linker Rand ≥ −2 − 0.24|y|).
3. **Endlich viele N** (kn.satz.polynom_halbexp_endlich): Für N ≤ 120 streng (Krawczyk, arb): eindeutige Lösung in der Box
   |a_k − ã_k| ≤ 2^−250 |L|^−k, und Σ|a_k||L|^k ≤ 2.2405.

**Offen = Ziel:** die Schranke (i)+(ii) für **alle** N — oder direkt ||p_N − q_N|| → 0 mit q_N = Knesers Taylor-Polynom
(das beweist die Vermutung sofort, ganz ohne die Eindeutigkeit).

## 3. Messwerte (alles in `gemeinsam/` bzw. in den Experimenten e49–e51)

- p_N(1) − T_e(1/2): 1.1e-7 (N=5), −3.3e-11 (15), 3.5e-15 (50), −1.0e-16 (100), −2e-19 (194). Rate ≈ 0.73/Grad bis N≈45,
  dann ≈ 0.93/Grad (Vorzeichenwechsel bei N≈51 — Übergang zwischen zwei Abklingraten).
- Koeffizienten konvergieren einzeln; p_N(L) → L etwa wie N^−4.8 (|p_N(L) − L| = 5.5e-10 bei N=100, 2.2e-11 bei 194).
- max_{|x|=|L|} |p_N| = 2.2327 für alle N (bei x = |L|); im Korridor max|p_N|/|L| = 0.9999012.
- Kondition von J (ungewichtet) wächst ≈ ×1.9/Grad.
- **Kantorovich um q_N scheitert** (Normen X: Σ|δ_k|1.3^k, Y: Σ|r_k|0.5^k): η_N = |J⁻¹F(q_N)|_X = wahrer Abstand |p_N − q_N|_X
  (1.3e-5, 5.3e-7, 5.4e-8, 9.4e-9, 2.5e-9, 8.0e-10 für N = 10..60), aber β_N = |J⁻¹|_{Y→X} wächst ×2.57/Grad; h_N = βγη
  von 0.49 (N=10) bis 9e15 (N=60).
- **Spektrum** der gewichteten Linearisierung (N=30): kleinste Singulärwerte 1.48e-12, 5.4e-12, 1.81e-11, 6.07e-11, … (Leiter
  ×≈3.4); kleinster Singulärvektor als Polynom auf |x| = 1.3 am grössten bei arg x ≈ ±80° (arg L = 76.6°).
- **Mechanismus (vermutet):** fast-singuläre Richtungen = abgeschnittene Tangentenvektoren anderer halber e-Funktionen,
  δ_θ = T'(A + 1/2)·θ(A), θ ½-antiperiodisch (θ = e^{±2πi(2m+1)w}), A = slog. Nahe L: A(x) ≈ log(x − L)/L + c,
  |e^{−2πiA}| ≈ |x − L|^{−4.447} (4.447 = Re(2πi/L)). Kneser selbst hat bei L nur eine schwache Singularität
  (Koeffizienten ~ k^{−6.45}|L|^{−k}). Der Rest F(q_N) steht fast senkrecht auf den fast-singulären Richtungen.
- Ersatzform: p ist genau dann Lösung, wenn das Taylor-Polynom von p an der Stelle a_0 = p(0) mit dem Jet von
  exp∘(p − a_0)^{−1} bis Ordnung N übereinstimmt (wahre Lösung: f = exp∘f⁻¹).

### Korrektur der Leitung (10.10. ~20:00, nach Strang C)
Die „fast-singulären Richtungen“ und das Wachstum ×2.57/Grad oben entstehen **nur durch die Gewichtspaarung** X ρ = 1.3,
Y s = 0.5. Mit **gleichen** Gewichten in Urbild und Bild ist J_N gut konditioniert: l¹ mit w_k = |L|^k gibt |J_N⁻¹| ≈ 1.8 → 7.8
(N = 10 … 194, langsam wachsend), mit w_k = |L|^k (1+k)⁴ ist |J_N⁻¹| ≤ 28.2 bis N = 194 (praktisch konstant). Kneser liegt in
diesem Raum genau für Exponenten < 5.45, die Konkurrenz-Tangenten (Koeff. ~ k^{2.45}|L|^{−k}) nicht. Arbeitsrahmen jetzt:
**Stabilität** (a) sup_N |J_N⁻¹|_{|L|,4} < ∞ + **Konsistenz** (b) α_N = |J⁻¹F(q_N)| → 0 + **Nichtlinearität** (c) Lipschitz-
Schranke auf einem Koeffizienten-Kegel (Worst-Case über die ganze Kugel wächst ≈ 1.25^N, entlang des echten Zuwachses winzig).
Ausserdem (Strang B): p_N ist im Inneren bis auf O(c_N²) die exakte halbe e-Funktion „Kneser + ein Störmodus e^{±2πiA}“ mit
Amplitude c_N ≈ (0.18 − 0.08i)·N^{1−2γ̄}, |c_N| ~ N^{−7.9} — die Konvergenz ist **polynomial**, nicht geometrisch
(die „Rate 0.93/Grad“ oben ist lokal dasselbe). Details: BRETT.md, c_literatur/BERICHT.md.

### Beweisgerüst der Leitung (10.10. 21:23)
1. **Linearer Kern** (Strang B, b_reformulierung §1, §7): p_N = Π_N f_K + O(2. Ordnung), Π_N = J_N⁻¹T_N J_∞.
   **Satz B-Max** (bewiesen bis Lemma 7.4/7.5 = Eigenschaften von f_K): (H*) |Π_N k_ζ(z)| ≤ C N^κ/dist(ζ, Ḡ') für Pole im
   Randband {g ≤ 1/N}, κ < 4.69 ⇒ sup_{D̄}|Π_N f_K − f_K| = O(N^{κ+1−4.69}). (H*) numerisch bis N = 300 mit C = 1, κ = 0.
   → **Strang C (Runde 3): (H*) beweisen.**
2. **Quadratische Kontrolle** p_N − Π_N f_K → 0: Kantorovich mit Betragsabschätzungen bis N ≈ 230, darüber offen;
   zweiter Newton-Schritt von q_N aus wächst ab N ≈ 250 in der Eckrichtung (D). → **Strang B (Runde 3)**, dazu Lemma 7.4/7.5.
3. **Gleichmässige Stabilität ist falsch** (D, bestätigt von B): Eckmode bei L, |J_N⁻¹| ~ e^{0.0124N} in ℓ²; Überlapp mit
   kohärenten Eingaben ~1e-8, Anteil an der Pol-Antwort ~N. → **Strang D (Runde 2): Theorie der Eckmode, bedroht sie die
   Vermutung für grosse N?** → **Strang A (Runde 2, Priorität): Π_N f_K bis N ≈ 1000 und echte p_N bis ≈ 500 rechnen.**
4. **Modelle** (Strang E): Satz K′ (rechnergestützt, Kontraktionsfall), Satz T (Translationsmodell: Auswahl kleinsten Typs,
   Fehler N^{1−Lücke} — erklärt den Faktor N), Möbius-Drehung mit scharfer Schwelle (Pol der Wurzel in W ⇔ Scheitern).
   → **Strang E (Runde 2): das ganze Gerüst (H*) + quadratisch im Möbius-Modell vollständig beweisen (Schablone, Prüfstein).**

### Stand nach Runde 1–5 (10.10. abends) — Leitung
Alle fünf Stränge abgeschlossen. **Konsolidiert in ARCHITEKTUR.md** (Strang B): Hauptsatz (K) + (H) + (A) ⇒ p_N → f_K gleichmässig auf D̄
mit Rate N^{κ−4.69}; offen ist nur die lineare Annahme (H) über die Eckmoden bei L (Kern: (a′) Rest-Kontraktion und (E_L) Randprofil der
adjungierten Eckmoden, p > 1.13). Numerisch bis N = 1500 bestätigt (Strang A). In den Graphen übernommen: kn.satz.polynom_halbexp_hauptsatz,
kn.verm.polynom_halbexp_eckmoden, kn.beob.polynom_halbexp_fehlergesetz, kn.beob.polynom_eckmode, kn.satz.halbexp_singularitaet_L,
kn.satz.polynom_modelle, kn.beob.polynom_moebius_schwelle; Prüfdatei experimente/e52_polynom_team.py (5/5).
Nächste Schritte siehe ARCHITEKTUR.md §7.

## 4. Dateien

- `gemeinsam/polynome_bis_194.json` — {"N": [a_0..a_N] als Strings, 60 Stellen}, N = 1..194.
- `gemeinsam/kneser_halbexp_koeff_80.json` — Taylor-Koeffizienten von f_K um 0, k = 0..80 (≈ 35 Stellen).
- `gemeinsam/omega_rand.json` — Rand von Ω = slog(D) (obere Hälfte), siehe e50.
- `gemeinsam/janis_notiz.md` — Janis' ursprüngliche Notiz.
- Experimente: `experimente/e49_polynom_halbexp.py`, `e50_halbexp_korridor.py`, `e51_polynom_schranke.py`.
- Kneser-Motor: `experimente/gemeinsam.py`, `kneser(art, basis, punkte, dps=40)` (art 'sexp'/'slog', basis 'exp(1)',
  komplexe Punkte erlaubt). **Nur Strang A darf ihn benutzen** (gemeinsamer 190-MB-Zwischenspeicher, gleichzeitige Nutzung
  zerstört ihn). Andere Stränge: Daten bei A anfragen (Eintrag im BRETT) oder aus `gemeinsam/` lesen.

## 5. Regeln für alle Stränge

- Python: `C:/Users/janis/miniconda3/python.exe` (mpmath 1.3, numpy 2.2, scipy 1.17, python-flint 0.9 mit arb/acb).
  Skripte mit dem Write-Werkzeug anlegen (Bash-Heredocs verschlucken Backslashes). `sys.stdout.reconfigure(encoding="utf-8")`.
- Nur im eigenen Ordner schreiben (und neue Dateien in `gemeinsam/`). **Nicht** ändern: `daten/`, `experimente/`, `baue.py`,
  `ansicht/`. Kein `baue.py`. Keine Git-Befehle.
- Rechenlast: einzelne Läufe ≤ ~20 min, höchstens ein schwerer Lauf gleichzeitig pro Strang, keine Mehrkern-Orgien.
- Ehrlich kennzeichnen: **bewiesen** / **Skizze (Lücken benennen)** / **numerisch** / **Vermutung**. Zahlen nur aus eigenen
  Läufen oder Dateien, Zeitstempel aus `date`.
- **BRETT.md** (in diesem Ordner) ist das gemeinsame schwarze Brett: zu Beginn lesen, nach jedem wichtigen Fund oben einen
  kurzen Eintrag anhängen (`## [Strang X] dd.mm. HH:MM — Titel` + 2–6 Zeilen + Dateipfad). Regelmässig neu lesen.
- Jeder Strang pflegt `<ordner>/BERICHT.md` (deutsch): Ziel, Versuche, Ergebnisse mit Zahlen, Beweise/Skizzen, Sackgassen,
  nächste Schritte. Am Ende eine kurze Zusammenfassung an die Leitung zurückgeben.
