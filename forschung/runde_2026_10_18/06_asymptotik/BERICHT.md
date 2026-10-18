# Runde 18.10.2026 — 06 Asymptotik, Singularitäten und Wachstum von Knesers sexp/slog (Basis e)

## Kurz in einfacher Sprache

- **Bei −2** sieht Knesers Tetration genau aus wie „log(z+2) plus eine glatte Funktion“. Diese glatte Funktion ist viel braver als gedacht:
  Sie ist bis −3 hin glatt (holomorph auf ℂ ohne (−∞,−3]). Die Entwicklung um −2 konvergiert daher genau im Radius 1.
- **Taylor-Koeffizienten um 0:** Sie sind *nicht* alle positiv (offene Frage aus gr.beob.kneser_ableitung). Bis n = 21 sind sie positiv,
  ab n = 21 wechseln sie das Vorzeichen wie (−1)^{n+1} — erzwungen durch den Logarithmus bei −2. Der Rest nach dem Logarithmus hat
  ein überraschendes Verhalten in **drei Stufen**: bis n ≈ 500 bestimmt ein Sattelpunkt auf der *positiven* Achse den Rest (die Funktion
  wird bei z ≈ 3 riesig, sexp(3) ≈ 3.8·10⁶), erst danach übernimmt die Singularität bei −3. Der Kipppunkt n = 506 wurde vorhergesagt und
  dann gemessen.
- **Nebenblätter:** Setzt man sexp über den Schnitt fort, gibt es **unendlich viele weitere logarithmische Verzweigungspunkte**, z. B. bei
  1.00431248 + 0.70941354i. Grund: Knesers sexp nimmt jeden Wert ausser 0 unendlich oft an (die 0 genau einmal, bei −1).
- **Oben in der komplexen Ebene** ist sexp eine konvergente Doppelreihe in e^{Lz} und e^{2πiz}; die reinen e^{kLz}-Koeffizienten sind
  genau p_k κ^k (p_k = Koeffizienten der Schröder-Umkehrfunktion Ψ = σ⁻¹). slog nahe L ist die Umkehrung davon (Hochzahlen j + nγ, ein Logarithmus).
- **slog auf der reellen Achse:** Es gibt *keine* Konstante mit slog(x) = log* x + c + o(1). Exakt gilt slog = (linear interpolierter log*)
  + eine periodische Welle P, scharf eingeschlossen in [−0.0078653, +0.0091705], Mittelwert 1/2 − ∫₋₁⁰ sexp = 0.00078182.
  Folge: Knesers halbe e-Funktion liegt in „Höhe“ gemessen höchstens 0.01704 neben der halben Stufe des linearen log*.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | g(z) := sexp(z) − log(z+2) ist holomorph auf ℂ∖(−∞,−3]; Sprung von sexp über (−3,−2) ist genau 2πi | **bewiesen** (aus Knesers Holomorphie + Funktionalgleichung) | Reihe um −2 trifft Motor bei −2.5+0.3i auf 3·10⁻³⁶, bei −2.8+0.05i auf 2·10⁻¹⁶ (Abbruch) | vermutlich neu als Aussage (die Schnitt-Lage ist bekannt, die Hebbarkeit über (−3,−2) nicht gefunden) |
| B2 | sexp(−2+ε) = log ε + Σ g_n εⁿ, Radius genau 1; g₀ = log sexp′(−1) = 0.08779780632645985179…, g₁ = −0.29721969736603008485…, g₂ = 0.27626137746255679474… | **bewiesen** (Radius); Koeffizienten **numerisch** | zwei Einstellungen (r=0.6/N=400/dps60, r=0.8/N=300/dps45), Abstand ≤ 2·10⁻³² | Entwicklung folklore, Radius neu |
| B3 | (−1)ⁿ g_n · n log n → 1, genauer 1 + (g₀−γ)/log n + O(log⁻²n) | **hergeleitet** (Flajolet–Odlyzko auf log(Λ−g₀) bei −3), Lücke: Restglied nicht abgeschätzt | Messung 0.852 bei n=140; Fit a = −0.461 gegen Theorie −0.489 | neu |
| B4 | Taylor um 0: c_n = (−1)^{n+1}/(n2ⁿ) + d_n, d_n = O(ρ⁻ⁿ) für jedes ρ < 3 | **bewiesen** (aus B1) | Abweichung von der Darboux-Formel relativ 7·10⁻¹² bei n=100, 1·10⁻¹⁸ bei n=140 | Darboux-Leitterm vermutlich bekannt, die exponentiell kleine Fehlerordnung 3⁻ⁿ neu |
| B5 | c_n > 0 für 0 ≤ n ≤ 21, c₂₂ < 0 erste negative Zahl; sign c_n = (−1)^{n+1} für alle n ≥ 21 | **bewiesen bis auf eine numerische Schranke** (max|g| = 177.6 auf |z|=2.5 aus einem Gitter; Cauchy für n ≥ 44, Rechnung für n < 44; numerisch bis 640) | zwei Einstellungen | beantwortet gr.beob.kneser_ableitung f5 („alle positiv?“: nein) |
| B6 | Drei Stufen von d_n: (i) Hayman-Sattel auf der positiven Achse, d_n ≈ g(R)R⁻ⁿ/√(2πb(R)), R g′/g = n; d_n > 0; (ii) Kipppunkt n* ≈ 506; (iii) d_n ≈ (−1/3)ⁿ/(n log n)·(1+(B−γ)/log n+…), B = log 3 + g₀ | **hergeleitet** (Sattelpunktmethode, Lücke: Beweis, dass der Sattel dominiert) + **numerisch** | S_n/d_n = 1.0133 (n=100), 1.0105 (n=140), 1.0042 (n=400); (d_n−S_n)/H_n = 1.129 stabil für n = 500…640; letzte positive Abweichung vom Muster n = 505 | neu |
| B7 | Knesers sexp hat genau eine Nullstelle (z = −1) und nimmt jeden Wert c ≠ 0 unendlich oft an | **bewiesen** (Denjoy–Wolff für log auf der oberen Halbebene + Kneser-Darstellung oben) | Familie w_j ≈ w₀ + jγ für c = 2πi nachgerechnet (Restfehler 3·10⁻⁴⁷) | vermutlich neu im Graph; Literatur nicht gefunden |
| B8 | Die Riemannsche Fläche von sexp hat ausser −2, −3, … unendlich viele weitere logarithmische Verzweigungspunkte: über w−1 für jedes w mit sexp(w) ∈ 2πiℤ∖{0}; Beispiel w−1 = 1.00431248259468800394 + 0.70941354086940362668i (sexp(w) = 2πi), w−1 = 1.15587108746 + 0.49064172727i (sexp(w) = 4πi) | **bewiesen** (Skizze; Blattwechsel über (−4,−3)) | w bei dps 60 und 45: Abstand 5·10⁻⁵⁴; sexp(w−1) = log(2πi) auf 30 Stellen | vermutlich neu |
| B9 | Oben: sexp(z) = L + Σ_{k≥1,m≥0} p_k κ^k e_{k,m} e^{(kL+2πim)z}, e_{k,m} = [q^m] exp(kL(θ(q)−a₀)), absolut konvergent auf ganz Im z > 0; Reihenfolge der Grössen: e^{Lz}, …, e^{5Lz}, dann e^{(L+2πi)z} | **bewiesen** (aus kn.konstr.kneser) | Motor gegen Ψ(κe^{Lz}): Fehler 2·10⁻¹⁴ (y=4) = Grösse des ersten Fourier-Glieds; mit ihm 1·10⁻²⁵ | Struktur klassisch (Kneser), explizite Koeffizientenformel C_{k,0} = p_kκ^k folklore |
| B10 | slog nahe L: slog(w) = ℓ/L − Σ_{n≥1} a_n e^{2πin·slog w}, ℓ = log(σ(w)/κ); Leitform ℓ/L − a₁(σ(w)/κ)^γ; Hochzahlen m + nγ, ein Logarithmus | **bewiesen** (Umkehrung von B9) | 4 Punkte |w−L| = 0.05…0.7: volle Umkehr trifft Motor auf 10⁻⁴⁸, Leitform mit a₁ auf 5·10⁻¹⁹ (y=3) | analog zu kn.satz.halbexp_singularitaet_L; für slog neu formuliert |
| B11 | slog(x) − slog_lin(x) = P(slog x − ⌈slog x⌉), P(u) = u+1−sexp(u) auf [−1,0]; kein Grenzwert; min P = −0.0078653491479050970569 (u = −0.80655771470683940698), max P = 0.0091704769530226506576 (u = −0.22232190393351360972), Nullstelle −0.53308891039947058146, Mittel 1/2 − ∫₋₁⁰sexp = 0.00078181717000866027139 | **bewiesen** (Form, kein Grenzwert); Zahlen **numerisch** | zwei Taylor-Sätze, 20 Stellen gleich | vermutlich neu (Schranken) |
| B12 | Halbe e-Funktion: |slog_lin(f(x)) − slog_lin(x) − ½| ≤ max P − min P = 0.0170358261009277477 für alle x | **bewiesen** aus B11 | Stichprobe 400 Punkte x ∈ [0.05, 5·10⁴]: max 0.016521 | neu |

## 1. Singularität bei −2 und die Riemannsche Fläche

### 1.1 Satz (B1). 
Voraussetzung: Knesers sexp ist holomorph auf ℂ∖(−∞,−2], reell und wachsend auf (−2,∞), sexp(z+1) = e^{sexp z} dort, Schwarz-Spiegelung.

*Beweis.* (a) Nullstellen: Ist sexp(w) = 0 und w−1 nicht auf dem Schnitt, dann e^{sexp(w−1)} = 0, unmöglich. Also liegt w−1 ∈ (−∞,−2], w ∈ (−∞,−1];
auf (−2,−1] ist sexp nur bei −1 null. Einzige Nullstelle: −1.
(b) Für x ∈ (−3,−2) gilt exp(sexp(x ± i0)) = sexp(x+1) < 0, also sexp(x±i0) = log|sexp(x+1)| ± iπ(2k+1) mit einem ganzen k, das wegen
Stetigkeit konstant ist; nahe −2 ist sexp(z) = log(sexp(z+1)/(z+2)) + log(z+2) mit reellem erstem Teil auf (−2,−1), also k = 0. Der Sprung
sexp(x+i0) − sexp(x−i0) = 2πi ist derselbe wie bei log(z+2). Also setzt sich g = sexp − log(z+2) stetig, nach Morera holomorph, über (−3,−2) fort.
Bei −2 selbst: g = log(sexp(z+1)/(z+2)), und sexp(z+1)/(z+2) ist nahe −2 holomorph mit Wert sexp′(−1) > 0. ∎

Folgerung (B2): Die Entwicklung um −2 hat Radius genau 1 (bei −3 ist g singulär: sexp(z) = log(log(z+3) + g(z+1)) ~ log log). Koeffizienten aus
g(−2+ε) = log(sexp(−1+ε)/ε) mit der Taylorreihe um −1 (`a4_minus2.py`). Ersten Werte (40 Stellen in `daten/g_minus2.json`):
g₀ = 0.0877978063264598517937, g₁ = −0.297219697366030084849, g₂ = 0.276261377462556794738, g₃ = −0.124963836571985612803,
g₄ = 0.0965781618001037110934, g₅ = −0.0692202637177039506677.

(B3) Bei ε = −1: mit η = 1+ε, Λ = log(1/η) ist der singuläre Teil log(Λ − g₀). Flajolet–Odlyzko ([ζⁿ]Λ^α ~ (log n)^α/n · Σ binom(α,k) (1/Γ)^{(k)}(0)/log^k n)
gibt [ζⁿ] log(Λ−c) = 1/(n log n) · [1 + (c−γ)/log n + (γ² − π²/6 − 2γc + c²)/log² n + …].

### 1.2 Satz (B7): Wertevorrat.
*Beweis.* Für c in der oberen Halbebene H bildet der Hauptlogarithmus H in den Streifen 0 < Im < π ⊂ H ab, ist kein Automorphismus und hat den Fixpunkt L ∈ H;
nach Denjoy–Wolff gilt log^{[k]}c → L. Nahe L gilt sexp(z) = Ψ(κ e^{L(z+θ(z)−a₀)}) für Im z > 0 (Kneser), mit θ − a₀ → 0 für Im z → ∞.
Für v = log^{[k]}c nahe L ist σ(v) klein, und die Gleichung z + θ(z) − a₀ = (log(σ(v)/κ) + 2πij)/L hat für jedes grosse j genau eine Lösung z_j nahe dem
rechten Seitenwert (Rouché; z+θ(z) ist für grosse Im z fast die Identität; Im(2πij/L) = j·Im γ, γ = 2πi/L = 4.44695 + 1.05794i). Dann sexp(z_j) = v und
sexp(z_j + k) = exp^{[k]}(v) = c. Die w_j = z_j + k sind verschieden. Für c in der unteren Halbebene Spiegelung; für c < 0 liegt log c auf Im = π und
der nächste Logarithmus in H; für c > 0 läuft die reelle Log-Kette nach endlich vielen Schritten ins Negative oder trifft genau 0 (c = sexp(n));
im letzten Fall nimmt man im ersten Schritt log c + 2πi ∈ H. ∎
Numerisch (`e8b_familie.py`): c = 2πi, k = 12: w_j = −5.5487+1.7663i +12, −1.1017+2.8242i +12, 3.3452+3.8822i +12, … mit Abstand genau γ.

### 1.3 Satz (B8): weitere Verzweigungspunkte.
Überquert z+1 den Schnitt (−3,−2) von unten nach oben, so wird sexp(z+1) wegen B1 zu sexp(z+1) − 2πi; das ist ein Zweig auf ganz ℂ∖(−∞,−2].
Der Zweig von sexp(z), der dabei entsteht (z quert (−4,−3)), ist log(sexp(z+1) − 2πi) + 2πil. Hat sexp(w) = 2πi, so ist dieser Zweig bei z = w−1 von der Form
log(sexp′(w)(z−w+1)) + holomorph — ein logarithmischer Verzweigungspunkt (sexp′ ≠ 0 nach kn.satz.ableitung_nullstellenfrei). Mit B7 gibt es unendlich viele.
Allgemein: ein Zweig b über w verschwindet nur, wenn exp^{[n]}(0) = sexp(w+n) für ein n; die Projektionen aller Verzweigungspunkte liegen in
{w − m : m ≥ 1, ∃n: sexp(w+n) = sexp(n−1)}. **Lücke:** dass jeder solche Turm von Logarithmen auch als Zweig erreicht wird, ist nur für n ≤ 2 ausgeführt.
Nicht entschieden: ob die Verzweigungspunkte dicht liegen.

Beispiel (dps 60 und 45, Newton mit Motor): sexp(w) = 2πi für w = 2.004312482594688003944312780444 + 0.709413540869403626679070125811i,
sexp(w−1) = log(2πi) auf 30 Stellen, sexp′(w) ≈ −11.6001 + 6.4945i.

## 2. Taylor-Koeffizienten um 0

### 2.1 Darboux mit exponentiell kleinem Fehler (B4).
Aus B1: c_n − [zⁿ]log(z+2) = [zⁿ]g, g holomorph auf |z| < 3. Also c_n = (−1)^{n+1}/(n2ⁿ) + O(ρ⁻ⁿ), ρ < 3.
Koeffizienten: Cauchy-DFT aus Motorwerten, zwei Einstellungen (r=1, N=512, dps 60 / r=1.5, N=400, dps 50; Übereinstimmung ≤ 10⁻⁴⁰ relativ bis n = 70),
für n bis 640 g direkt (Logarithmus vor der DFT abgezogen: r=1.9/N=1400/dps140 und r=1.8/N=1300/dps150; gleich auf 30 Stellen bis n = 520, 10⁻⁷ bei n = 640).

### 2.2 Drei Stufen (B6).
- n ≲ 500: Sattel auf der positiven Achse, R(n) von 2.49 (n=20) bis 3.23 (n=500). Hayman: S_n/d_n = 1.045, 1.026, 1.013, 1.0105, 1.0042 für n = 20, 40, 100, 140, 400 (Fehler O(1/n)).
- Kipppunkt: vorhergesagt aus log₁₀S_n ≈ log₁₀|H_n| bei n ≈ 500 (`a6_wechsel.py`), gemessen: d_n > 0 bis 505, ab 506 Vorzeichen (−1)ⁿ (`a7_wechsel_pruef.py`).
- n ≳ 520: (d_n − S_n)/H_n = 1.1285 … 1.1282 für n = 520…640, mit H_n = (−1/3)ⁿ/(n log n)·(1+(g₀−γ)/log n). Mit dem richtigen B = log 3 + g₀ und dem
  dritten Glied sagt die Theorie 1.064 · (1/(n log n)) voraus, gemessen 1.043 (n=600) — verträglich (Asymptotik in 1/log n).

### 2.3 Vorzeichen (B5).
c_n > 0 für n ≤ 21, c₂₂ = −1.5492416655467695218·10⁻⁹. Für n ≥ 21: c_n(−1)^{n+1} > 0 numerisch bis 640. Für einen Beweis für alle n: |d_n| ≤ M·2.5⁻ⁿ
mit M = max_{|z|=2.5}|g|; mit M ≤ 180 folgt |d_n| < 1/(2n2ⁿ) ab n ≈ 50. M: `a12_schranke.py` (Gitter 400 und 300 Punkte auf |z| = 2.5, dps 40): max|g| = 177.41 bzw. 177.25 am Gitterpunkt neben arg = 0, g(2.5) = 177.6114 (Maximum auf der positiven Achse, wie beim Sattel erwartet), nahe arg = π nur 0.33. Mit M = 180 gilt |d_n| < 1/(2n2ⁿ) für alle n ≥ 44; zusammen mit der Rechnung bis n = 640 folgt das Vorzeichenmuster für alle n ≥ 21. Die Schranke M ist gitterbasiert, kein Intervallbeweis.

## 3. Oben, slog bei L

(B9) Aus sexp = Ψ(κ e^{Lz} E(q)), E(q) = exp(L Σ_{n≥1} a_n qⁿ), q = e^{2πiz}, Ψ ganz: sexp = Σ_k p_k κ^k e^{kLz} E(q)^k. Absolute Konvergenz auf ganz H:
Σ|p_k||κ e^{Lz}|^k Ẽ(|q|)^k < ∞, Ẽ(t) = exp(|L|Σ|a_n|tⁿ), endlich für t < 1. Also C_{k,0} = p_k κ^k (die „σ⁻¹-Struktur“ ist exakt), C_{1,1} = κ L a₁.
Grössenordnung: |e^{kLz}| = e^{−1.337ky}, |e^{(L+2πi)z}| = e^{−7.62y}: die ersten fünf Potenzen von e^{Lz} kommen vor der ersten Fourier-Korrektur.
p₁…p₆: 1, −0.151315−0.296749i, −0.036976+0.098731i, 0.025812−0.017387i, −0.007944+0.000579i, 0.001972+0.000838i.

(B10) Umkehrung; Prüfung `a10_slogL.py` (z vorgeben, w = sexp(z), zurückrechnen): j = 0 auf dem Hauptblatt, Restfehler 10⁻⁴⁸.

## 4. slog reell und die halbe e-Funktion (B11, B12)
slog_lin(x) = n − 1 + log^{[n]}x mit log^{[n]}x ∈ (0,1]. Beide erfüllen A(eˣ) = A(x)+1, also ist die Differenz 1-periodisch in der Höhe; auf [0,1]
in u = slog x ∈ [−1,0] ist sie u − (sexp(u) − 1) = P(u). Kein Grenzwert: P ≡ const hiesse sexp linear auf [−1,0], dann überall (Analytizität) — Widerspruch zu sexp(1) = e.
Mittel (in u und in s gleich): ∫₋₁⁰P = 1/2 − ∫₋₁⁰sexp (gr.beob.integrale_kneser). Halbe e-Funktion: slog f(x) = slog x + 1/2, also folgt B12 aus zweimal B11.
**Antwort auf Frage (3):** „slog(x) = log* x + const + o(1)“ ist falsch; richtig ist die Form B11 mit einer periodischen Welle der Breite 0.017.
**Frage (5):** Eine Asymptotikreihe aus elementaren Funktionen gibt es nicht (kn.satz.halbe_exp_wachstum); B12 ist die genaueste elementare Beschreibung, die ich fand.

## Hinweis zum Motor
Der Motor (`fatou_fork.gp`) wurde während der Runde von einer anderen Sitzung geändert (Zeitstempel 13:51). Danach hing er bei dps 30 an Punkten mit Im z ≈ 2.5 (z. B. 2.5i), bei dps 40 nicht; ausserdem hing er an Punkten genau auf dem Schnitt (−∞,−2]. Alle Zahlen oben stammen aus je zwei Einstellungen, die übereinstimmen.

## Skripte (alle in diesem Ordner; Motor-Zwischenspeicher `cache/`)
`motor.py` (Kopie von etappe2/kern.py mit eigenem Speicher), `lib3_kopie.py` (Kopie etappe3/lib3.py, θ-Daten nur lesend),
`e1_taylor.py` (Taylor um 0/−1), `e7_gkoeff.py` (d_n bis 650), `a1_darboux.py`, `a1b.py`, `a2_sattel.py`, `a3_vorzeichen.py`, `a4_minus2.py`, `a5_fit.py`,
`a6_wechsel.py`, `a7_wechsel_pruef.py`, `e8_zweigpunkte.py`, `e8b_familie.py`, `a9_oben.py`, `a10_slogL.py`, `a11_slog_reell.py`, `a12_schranke.py`, `a13_halbexp_band.py`.
Reihenfolge: e1 (0 1 512 60; 0 1.5 400 50; −1 0.6 400 60; −1 0.8 300 45), e7 (1.9 1400 140; 1.8 1300 150), dann die a-Skripte.

## Vorschläge für Grapheinträge
- kn.satz.minus2_hebbar — „sexp − log(z+2) ist holomorph bis −3“ — satz — bewiesen
- kn.kor.minus2_entwicklung — „Entwicklung bei −2 mit Radius genau 1, Koeffizienten g_n ~ (−1)ⁿ/(n log n)“ — korollar — bewiesen (Radius) / hergeleitet (Asymptotik)
- kn.satz.taylor_darboux — „Taylor-Koeffizienten um 0: (−1)^{n+1}/(n2ⁿ) + O(ρ⁻ⁿ), ρ<3; Vorzeichen ab n = 21 alternierend, c₂₂ erste negative Zahl“ — satz — bewiesen / numerisch
- kn.beob.taylor_drei_stufen — „Rest der Taylor-Koeffizienten: Hayman-Sattel bis n ≈ 505, dann Singularität bei −3“ — beobachtung — hergeleitet + numerisch
- kn.satz.wertevorrat — „sexp nimmt jeden Wert ≠ 0 unendlich oft an, die 0 nur bei −1“ — satz — bewiesen
- kn.satz.nebenblatt_zweigpunkte — „Unendlich viele weitere log-Verzweigungspunkte auf Nebenblättern, z. B. 1.0043125+0.7094135i“ — satz — bewiesen (Skizze)
- kn.kor.doppelreihe_oben — „sexp = L + Σ p_kκ^k e^{kLz}E(q)^k auf ganz Im z>0“ — korollar — bewiesen
- kn.kor.slog_bei_L — „slog nahe L: log(σ/κ)/L − a₁(σ/κ)^γ − …“ — korollar — bewiesen
- gr.satz.slog_logstern_welle — „slog = linearer log* + Welle in [−0.0078653, 0.0091705], Mittel 0.00078182; kein Grenzwert“ — satz — bewiesen / numerisch
- kn.kor.halbexp_band — „halbe e-Funktion: Höhe um ½ ± 0.0170358 im linearen log*“ — korollar — bewiesen
- Kante: kn.satz.taylor_darboux beantwortet gr.beob.kneser_ableitung (f5).

## Offene Fragen
1. Liegen die Verzweigungspunkte aller Nebenblätter dicht in ℂ? (Menge {w−m : sexp(w+n) = sexp(n−1)}.)
2. Strenger Beweis der Sattel-Dominanz (B6) und Intervall-Schranke für M (B5) — dann ist das Vorzeichenmuster für alle n bewiesen.
3. Gibt es im Streifen −2 < Re ≤ −1 Lösungen von sexp(w) ∈ 2πiℤ∖{0}? (Dann lägen Nebenblatt-Verzweigungspunkte nahe am Hauptschnitt.) Die Gittersuche fand nur w mit Re w ≈ 2.
4. Das stabile Verhältnis 1.128 in B6: drittes und viertes Glied der 1/log n-Reihe ausrechnen (mit B = log 3 + g₀).
5. Andere Basen: Kipppunkt n*(b) als Funktion der Basis.
