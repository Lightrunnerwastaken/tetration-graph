# Runde 25.10.2026 — 07 Riemannsche Fläche von Knesers slog (Basis e)

Beginn 13:23, Ende 13:45 (Uhr aus `date`). Motor: eigene Kopie `motor.py`, Zwischenspeicher `daten/motor_cache.json`.

## Kurz in einfacher Sprache

- Knesers slog ist um 0 eine Potenzreihe, die genau bis |L| = 1.37456 reicht; am Rand sitzen nur die beiden Fixpunkte L und L̄, beide als Logarithmus mit Vorfaktor 1/L.
- Daraus folgt ein ganzes „periodisches Blatt“: links der Senkrechten Re w = Re L ist slog(w) = slog(eʷ) − 1, also 2πi-periodisch, und auf dieser Senkrechten sitzt ein ganzes **Gitter** von Logarithmus-Punkten L + 2πik und L̄ + 2πik.
- Rechts davon liegt Knesers gewohntes Blatt. Dort sind dieselben Gitterpunkte **harmlos**: slog(L + 2πi) = 2.03397945… + 0.62603889…i. Damit hat Knesers sexp konkrete Stellen, an denen es den Fixpunkt L trifft: sexp(3.03397945… + 0.62603889…i) = L.
- Der Übergang zwischen beiden Blättern ist ein Umlauf um L. Er hat eine exakte Formel (Verschiebung um γ = 2πi/L in der „regulären Höhe“ u). Die Formel stimmt mit dem Motor auf 15 Stellen.
- Gesamtbild: Die Verzweigungspunkte aller Blätter bilden eine Menge, die unter exp in beide Richtungen erhalten bleibt. Sie enthält L, L̄, 0 und alles, was exp oder ein log-Zweig daraus macht (z. B. 1, e, eᵉ, 2πik, L + 2πik). Diese Menge liegt **dicht in ℂ**.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| R1 | Taylorreihe A(w) = Σ Aₙwⁿ von slog um 0: Aₙ = −2 Re(1/(n L^{n+1}))·(1 + O(n^{−Re γ})); Radius genau \|L\|, Leitglied log(1 − w/L)/L + konj., auf dem Kreis \|w\| = \|L\| keine weitere Singularität mit Exponent Re α ≲ 4 | **numerisch** | Quotient = 1 ± 2·10⁻¹¹ bei n = 250, 300; Rest fällt mit Steigung −4.5 … −5.4 (Theorie −4.447); `s1_radius.py` | Radius \|L\| ist Folklore (Walker/Paulsen); Exponent des Rests neu formuliert |
| R2 | Periodisches Blatt P: slog(w) = A(eʷ) − 1 ist holomorph auf {Re w < Re L}, 2πi-periodisch, und hat auf der Randgeraden Re w = Re L genau die Log-Punkte L + 2πiℤ und L̄ + 2πiℤ, jeweils slog ≈ log(w − w₀)/L + c | **bewiesen** aus R1 (\|e^w\| = e^{Re w} < e^{Re L} = \|L\|); „genau“ hängt an R1 | Motor = P an 7 Punkten (Abstand ≤ 10⁻³⁰, z. B. −0.5+5i, −3+7.9i); Log-Test an L, L+2πi, L+4πi von links: Differenz r=10⁻³ minus r=10⁻² = log(0.1)/L ± 0.0022 (Grösse des (w−w₀)-Glieds) | vermutlich neu (das Gitter auf Re w = Re L) |
| R3 | Knesers Blatt K (rechts der Senkrechten, wie im Motor) ist an L + 2πik, k ≠ 0, **regulär**: s*₁ = slog(L+2πi) = 2.0339794527621306885320165722 + 0.6260388922275826519647882805i, s*₂ = 2.1711479776960512760791065186 + 0.4626241915670184355004755120i, s*₃ = 2.2440413113722070178255433948 + 0.3910177242666462632551004951i, s*₋₁ = 1.9164381432450097009166882644 − 0.8024346897198475543050758983i | **numerisch** | dps 30 und dps 40 gleich auf 28 Stellen; sexp(s*ₖ) = L + 2πik und sexp(s*ₖ+1) = L auf 10⁻³⁰ (= Versatz des Messpunkts); `s4_schnitt_und_sstern.py`, `daten/s4_dps*.txt` | neu |
| R4 | Folgerung: Knesers sexp trifft den Fixpunkt L, z. B. bei z_L = s*₁ + 1 = 3.03397945276213068853 + 0.62603889222758265196i, ebenso bei s*ₖ + 1 + n (n ≥ 0). Über L hat die Fläche also auch reguläre Punkte | **numerisch** (Werte), **bewiesen** (sexp(z_L+n) = L aus sexp(z_L) = L) | wie R3 | neu (konkrete Lösungen von sexp(z) = L) |
| R5 | Exakte Monodromie bei L: mit u = regulärer slog bei L (log σ/L) und slog = u + ϑ(u), ϑ 1-periodisch, gilt beim Umlauf (gegen den Uhrzeiger) u ↦ u + γ, also M_L(s) = (u+γ) + ϑ(u+γ), u = (id+ϑ)⁻¹(s). Alle Zweige nahe L: s_k = (u+kγ) + ϑ(u+kγ), k ∈ ℤ, paarweise verschieden; M_L vertauscht mit s ↦ s+1 | **bewiesen** (σ hat bei L eine einfache Nullstelle, Knesers Darstellung) + **numerisch** | Sprung des Motor-Blatts über die Senkrechte Re w = Re L, y = 1.8 und 2.0: Vorhersage M_L(links) = rechts auf 2·10⁻¹⁵ (= Versatz 10⁻¹⁵ der Messpunkte), dps 30 und 40 gleich; bei y = 2.5, 3 Abweichung 3·10⁻¹¹, 6·10⁻⁸ durch Abbruch der ϑ-Reihe nach 6 Gliedern; `s5b_monodromie_genau.py` | Formel neu ausgeschrieben; Umkehrung von kn.kor.slog_bei_L |
| R6 | Der Motor klebt P (links) und K (rechts) längs der Halbgeraden {Re L + iy : y ≥ Im L} (und gespiegelt) zusammen; der Sprung dort ist genau M_L | **numerisch** | Sprung ≈ γ bei y = 2 (4.4471+1.0554i gegen γ = 4.4470+1.0579i), exakt M_L wie R5; FG-Rest k(w) − k(eʷ) + 1 = \|γ\| = 4.57 direkt über L | Hinweis für Motor-Nutzer |
| R7 | Ψ = σ⁻¹ (Schröder-Umkehrung bei L) ist ganz, ohne kritische Punkte: Ψ′(ζ) = exp(Σ_{k≥1}(Ψ(ζ/Lᵏ) − L)); Ψ lässt 0 aus | **bewiesen** (aus Ψ(Lζ) = e^{Ψ(ζ)}: Ψ′(ζ) = e^{Ψ(ζ/L)−L} Ψ′(ζ/L), e^L = L, Produkt konvergiert wegen \|L\| > 1) | — | Produktformel vermutlich neu |
| R8 | 0 ist asymptotischer Wert von Knesers sexp: z(X) = 2 + slog_K(X + iπ), X → +∞ reell, gibt sexp(z(X)) = exp(−e^X) → 0. Über w = 0 hat die Fläche daher einen logarithmischen Punkt: slog_(0),k(w) = 2 + slog_K(log(−log w − 2πik) + iπ), k ∈ ℤ, alle → ∞ für w → 0 | **bewiesen** (Lücke nur: K-Blatt enthält {Re ζ > R, π/2 < Im ζ < 3π/2}, durch slog_K(ζ) = slog_K(log ζ) + 1 klar) | — | vermutlich neu |
| R9 | Die Menge B aller Projektionen von Verzweigungspunkten der slog-Fläche ist voll invariant: exp(B) ⊂ B und exp⁻¹(B) ⊂ B (alle log-Zweige). Mit L, L̄, 0 ∈ B enthält B die vollen Bahnen G(L) ∪ G(L̄) ∪ G(0), z. B. 1, e, eᵉ, 2πik, log(2πik)+2πij, L+2πik | **bewiesen** (Hochheben von Wegen durch die Überlagerung expⁿ: ℂ∖expⁿ⁻¹(Oₙ) → ℂ∖Oₙ, Oₙ = {0,1,…,exp^{n−1}(0)}; Einzelheiten unten) | — | vermutlich neu |
| R10 | B liegt dicht in ℂ | **bewiesen** aus R9 und J(exp) = ℂ (Misiurewicz 1981): Rückwärtsbahnen nicht-ausnahmehafter Punkte liegen dicht in J | — | beantwortet für slog die offene Frage 1 der Runde 18 (dort für sexp gestellt) |
| R11 | Jeder Punkt w₀ ∈ expⁿ⁻¹(L) ist auf einem passenden Blatt ein Log-Punkt mit slog ≈ log(w − w₀)/L + c, Monodromie = M_L (bis auf Verschiebung um n, die mit M_L vertauscht) | **bewiesen** (slog = A∘expⁿ − n lokal, expⁿ lokal biholomorph) | R2 ist der Fall n = 1 | neu |
| R12 | Die slog-Fläche hat keine algebraischen Verzweigungspunkte (alle Singularitäten transzendent/logarithmisch) | **Vermutung** | über dem Bild des Hauptblatts von sexp bewiesen (sexp′ ≠ 0, kn.satz.ableitung_nullstellenfrei); auf Nebenblättern sexp′(z) = sexp(z)·sexp′(z−1) und sexp(z) ≠ 0 ausser bei „Log-Punkt + 1“ | — |
| R13 | B = G(L) ∪ G(L̄) ∪ G(0) genau | **Vermutung** | Bild: slog = u + ϑ(u), u verzweigt nur über {L} ∪ O (it.satz.regulaerer_erzeuger_fortsetzung + R7); keine Gegenbeispielsuche gemacht | — |
| R14 | Monodromie-Formel R5 mit 6 Fourier-Gliedern von ϑ nicht nutzbar, sobald Im u < −1 (bei y = 5: grösstes Glied 10¹³) | **negativ** | `s5_monodromie.py` | — |

## Beweise und Einzelheiten

**R2.** Für Re w < Re L ist \|eʷ\| = e^{Re w} < e^{Re L} = \|e^L\| = \|L\|. Nach R1 ist A auf der offenen Kreisscheibe \|ζ\| < \|L\| holomorph, also ist A(eʷ) − 1 auf der Halbebene holomorph und 2πi-periodisch. Nahe 0 ist A(eʷ) − 1 = A(w) (Funktionalgleichung, Probe 10⁻⁵⁰), also ist es dieselbe Fortsetzung. Für w → L + 2πik von links geht eʷ von innen gegen L, und A(eʷ) − 1 ≈ log(eʷ − L)/L ≈ log(L(w − L − 2πik))/L: ein Log-Punkt. An allen anderen Randpunkten ist eʷ ein regulärer Randpunkt von A (falls R1 „keine weitere Singularität“ stimmt).

**R5.** Kneser: slog = u + ϑ(u) mit u(w) = log(σ(w))/L, σ(eʷ) = Lσ(w), σ′(L) ≠ 0, ϑ 1-periodisch (Koeffizienten κ_m aus `polynom_halbexp/gemeinsam/kneser_singularitaet_L.json`, nur gelesen). Ein Umlauf um L ändert log σ um 2πi, also u um γ. Da u ↦ u + γ mit u ↦ u + 1 vertauscht und ϑ 1-periodisch ist, vertauscht M_L mit s ↦ s + 1. Die Zweige s_k sind verschieden, weil u + kγ verschieden sind und (id+ϑ) dort umkehrbar ist (Im u gross).

**R7.** Aus Ψ(Lζ) = exp(Ψ(ζ)) folgt LΨ′(Lζ) = exp(Ψ(ζ))Ψ′(ζ), also Ψ′(ζ) = (e^{Ψ(ζ/L)}/L)·Ψ′(ζ/L) = e^{Ψ(ζ/L) − L}Ψ′(ζ/L) wegen e^L = L. Iteriert und Ψ′(0) = 1: Ψ′(ζ) = exp(Σ_{k≥1}(Ψ(ζL^{−k}) − L)); die Summe konvergiert lokal gleichmässig, weil Ψ(x) − L = O(x) und \|L\| > 1. Ψ(ζ) = exp(Ψ(ζ/L)) ≠ 0.

**R8.** Für X > 0 liegt X + iπ rechts von L und über Im L; dort ist slog_K holomorph und sexp(slog_K(X + iπ)) = X + iπ. Dann sexp(z(X) − 1) = e^{X+iπ} = −e^X und sexp(z(X)) = exp(−e^X) → 0, während z(X) → ∞ (slog_K(X+iπ) ≈ slog(log X) + 1). Die Umkehrung längs der Bildkurve ist ein Zweig von slog; ersetzt man log w durch log w + 2πik, so bleibt das Argument von −log w − 2πik in (−π/2, π/2), also Im(log(−log w − 2πik) + iπ) ∈ (π/2, 3π/2), noch über Im L = 1.337: alle Zweige leben auf K und sind verschieden.

**R9.** (a) exp(B) ⊂ B: Als Keime gilt längs jedes Weges slog = (slog längs exp∘Weg) − 1; wäre slog längs exp∘γ am Endpunkt exp(w₀) holomorph, dann auch längs γ bei w₀ (exp lokal biholomorph). (b) exp⁻¹(B) ⊂ B: Sei b ∈ B, erreicht längs δ ab 1 (vorher [0,1]). δ lässt sich (in kleinen Scheiben, in denen der benutzte Zweig holomorph ist) von 0 wegschieben und unter exp von 0 aus hochheben; das Ende ist log b + 2πik₀. Vorgeschaltet: das Stück 0 → 2πik auf der imaginären Achse; dort ist slog = A(e^{it}) − 1 (A holomorph auf \|ζ\| < \|L\| ⊃ Einheitskreis), danach endet der gehobene Weg bei log b + 2πi(k + k₀). Längs dieses Wegs ist slog = (slog längs δ) − 1, singulär am Ende. Also liegt jeder Punkt von exp⁻¹(b) in B. (c) Für L direkt (R11): Weg δ₀ von Eₙ = expⁿ(0) knapp über der reellen Achse bis nahe 0 und in der Kreisscheibe zu L, unter expⁿ hochgehoben; alle Punkte der Faser erreicht man, weil die Schleifen um O_n in einer einfach zusammenhängenden Umgebung der reellen Strecke liegen können, in der A eindeutig ist.

**R10.** J(exp) = ℂ; für eine ganze transzendente Funktion liegt die Rückwärtsbahn jedes nicht-ausnahmehaften Punkts dicht in J. L (und 1) sind nicht ausnahmehaft (einziger Ausnahmewert von exp ist 0). Also ist ⋃ₙ expⁿ⁻¹(L) ⊂ B dicht.

## Gesamtbild der Fläche (zusammengefasst)

1. **Kern:** Kreisscheibe \|w\| < \|L\| (R1).
2. **Periodisches Blatt P:** Halbebene Re w < Re L, 2πi-periodisch, Rand-Gitter L + 2πiℤ, L̄ + 2πiℤ von Log-Punkten mit Vorfaktor 1/L (R2). Das ist die Halbebene, die der Motor links benutzt.
3. **Kneser-Blatt K:** rechts davon; dort sind die Gitterpunkte (k ≠ 0) regulär (R3), die Werte s*ₖ wandern langsam nach rechts (≈ slog(log 2πk + iπ/2) + 1).
4. **Naht:** Halbgeraden von L nach oben und von L̄ nach unten auf Re w = Re L, Sprung = exakte Monodromie M_L (R5, R6).
5. **Über L:** auf dem Hauptzweig Log-Punkt unendlicher Ordnung mit Zweigen s_k (R5); auf anderen Blättern regulär, Werte z mit sexp(z) = L (R4).
6. **Über 0, 1, e, eᵉ, …:** auf dem Hauptblatt regulär, auf Nebenblättern Log-Punkte (R8, R9).
7. **Global:** Verzweigungsmenge B voll exp-invariant und dicht (R9, R10); vermutlich genau die drei vollen Bahnen von L, L̄, 0 (R13) und nur logarithmische Punkte (R12).

## Skripte (alle in diesem Ordner)

- `s1_radius.py` — Darboux-Prüfung der Taylorkoeffizienten (liest `polynom_halbexp/gemeinsam/kneser_slog_koeff_300.json`).
- `s2_blaetter.py 30` — Motor-Blatt gegen A(eʷ) − 1 links der Senkrechten.
- `s3_gitter_L.py 30` — Log-Test an L, L+2πi, L+4πi aus vier Richtungen.
- `s4_schnitt_und_sstern.py 30|40` — Sprung über die Senkrechte, s*ₖ; Ausgaben `daten/s4_dps30.txt`, `daten/s4_dps40.txt`.
- `s5_monodromie.py`, `s5b_monodromie_genau.py 30|40` — Monodromie-Formel gegen den Motor (liest κ_m aus `kneser_singularitaet_L.json`).
- Vorsicht (eigener Fehler, behoben): mpmath `-lambertw(-1)` ist L̄, nicht L.

## Vorschläge für Grapheinträge

- kn.satz.slog_periodische_halbebene — „slog = slog(eʷ) − 1 ist auf Re w < Re L holomorph und 2πi-periodisch; Log-Gitter L + 2πiℤ, L̄ + 2πiℤ auf dem Rand“ — satz — bewiesen (aus numerischem Radius)
- kn.beob.slog_taylor_radius_L — „Taylorreihe von slog um 0: Radius |L|, Aₙ ≈ −2Re(1/(nL^{n+1})), Rest n^{−Re γ}“ — beobachtung — numerisch
- kn.beob.sexp_trifft_L — „sexp(z) = L bei z = 3.0339794527621306885 + 0.6260388922275826520i und s*ₖ + 1 + n“ — beobachtung — numerisch
- kn.satz.slog_monodromie_L — „Umlauf um L: u ↦ u + γ, M_L(s) = (u+γ) + ϑ(u+γ)“ — satz — bewiesen + numerisch (15 Stellen)
- it.satz.schroeder_umkehr_ohne_kritische_punkte — „Ψ′(ζ) = exp Σ(Ψ(ζ/Lᵏ) − L), Ψ ohne kritische Punkte, lässt 0 aus“ — satz — bewiesen
- kn.satz.sexp_asymptotischer_wert_null — „0 ist asymptotischer Wert von sexp längs 2 + slog_K(X + iπ)“ — satz — bewiesen
- kn.satz.slog_verzweigung_dicht — „Verzweigungsmenge der slog-Fläche ist exp-invariant, enthält die Bahnen von L, L̄, 0 und ist dicht in ℂ“ — satz — bewiesen
- kn.verm.slog_verzweigung_genau — „Verzweigungsmenge = volle Bahnen von L, L̄, 0; nur logarithmische Punkte“ — vermutung
- Kante: kn.satz.slog_verzweigung_dicht → offene Frage 1 aus runde_2026_10_18/06 (dort für sexp)

## Offene Fragen

1. Strenger Beweis von R1 („nur L, L̄ auf dem Kreis \|w\| = \|L\|“), etwa über das Bild der Kreisscheibe unter slog.
2. R13: Ist B genau die Vereinigung der drei Bahnen? Gegenbeispielsuche: Blätter um 0 (R8) numerisch abtasten.
3. Gilt die Dichte auch für die Verzweigungspunkte von sexp (Runde 18, Frage 1)? Aus R9/R10 folgt sie nicht direkt.
4. Wie sieht die Naht (R6) oberhalb von L + 2πi aus? Die ϑ-Reihe divergiert dort (R14); man bräuchte ϑ auf ganz ℂ∖(Schnitte).
5. Monodromie über 0 auf dem Blatt aus R8 als geschlossene Formel: s ↦ slog_b(sexp(s − 1) + 2πi) + 1 — welche Gruppe erzeugen M_L, M_L̄ und M_0 zusammen?
