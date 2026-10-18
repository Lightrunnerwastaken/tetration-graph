# 07 Ränge — Rangleiter, Grenzprofil und Rangverschiebung

Forschungsrunde 18.10.2026, Gebiet „Rangleiter und gebrochene Ränge“. Alle Zahlen sind mit den Skripten in diesem Ordner reproduzierbar (siehe Ende).

## Kurz in einfacher Sprache

Für kleine Basen (zum Beispiel b = √2, 1.2, 1.5) bleiben alle höheren Rechenarten ab einem Rang beschränkt. Dann ist jede neue Stufe eindeutig festgelegt (Koenigs-Konstruktion), und man kann die Leiter billig bis Rang 16 rechnen. Dabei kam Folgendes heraus:

1. **Grenzprofil.** Mit wachsendem Rang wird H_r(b,x) zu einer Rampe mit Deckel, min(x+1, b). Sie steigt wie „plus 1“ und bleibt dann bei b stehen. Bei b = 2 ergibt das genau das schon beobachtete Bild H_r(2,x) → x+1 auf [0,1].
2. **Gesetz für die Annäherung.** Der Fixpunkt p_r jeder Stufe sinkt gegen die Basis b. Die Geschwindigkeit hängt nur von der Zahl k = (2−b)/(b−1) ab: Die relative Abnahme pro Rang ist etwa k/r. Ein einfaches Zwei-Größen-Modell gibt die echte Leiter Schritt für Schritt auf 0.02–2 % wieder.
3. **Rangverschiebung.** Ändert man die Tetration mit einer θ-Welle ab, die die ganzen Höhen nicht verändert, dann ist die ganze spätere Leiter nur **um einen festen Bruchteil eines Rangs verschoben** (bei ε = 0.05 um −0.3102 Ränge). Alle anderen Abweichungen sterben schneller als jede geometrische Folge. Die θ-Freiheit wird hoch oben in der Leiter also zu einer einzigen Zahl, einer Art „Rangkoordinate“ einer Funktion. Das ist das funktionale Gegenstück zur Fatou-Koordinate.
4. **Negativ:** Eine Konvexitätsannahme würde η_N → 2 beweisen. Sie ist an den echten kritischen Basen aber ab N = 5 verletzt. Für die η/e-Folge fand ich keine saubere neue Konstante.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | Für 1<b<e^{1/e} (und für b<2 ab dem Rang, an dem die Leiter beschränkt wird): jede Koenigs-Stufe H_{r+1}=σ_r^{-1}(λ_r^x σ_r(1)) ist streng wachsend mit H_{r+1}↑p_r. Außerdem gilt b < p_{r+1} < p_r und q_r:=H_{r+1}(b,2) ist streng fallend | **bewiesen** (Voraussetzung: der kleinste Fixpunkt ≥1 ist jeweils nicht entartet, λ_r<1) | Beweis §1 | vermutlich neu (einfach) |
| B2 | Grenzprofil: lim_{r→∞} H_r(b,x) = min(x+1, b) für x ≥ 0, 1<b<2. Insbesondere p_r → b | **Vermutung** | b=√2, 1.2, 1.5 bis r=16: monotone Annäherung, Plateau-Formel trifft das Profil auf 0.01–1 %; Mechanismus hergeleitet (§2) | neu |
| B3 | Plateau-Modell: δ_{r+1}=δ_r(1−μ_r), λ_{r+1}=δ_r L_r μ_r, μ_r=λ_r^{b−1+δ_{r+1}} (δ_r=p_r−b, L_r=−ln λ_r) | **hergeleitet mit Lücke** (Fehler O(L_r δ_r) aus der Linearisierung von σ_r im Plateau, Abschätzung fehlt) | Einschrittfehler bei b=1.2: δ 1.6e−4, λ 4.3e−4 (r=14); bei √2: 2.4e−3 bzw. 1.6e−2, beide fallend | neu |
| B4 | Asymptotik: r μ_r → k, Δ(1/μ_r) = 1/(k+1/M_r) + …, M_r=−ln μ_r; δ_r ~ (b−1)k^k / (r^k ln r); L_r ~ ln r/(b−1); k=(2−b)/(b−1) | **hergeleitet** für das Modell (formale Asymptotik), **numerisch** für das Modell bis r=10⁷; für die echte Leiter **Vermutung** | Modell bei r=10⁷: Δ(1/μ)=0.676678 gegen Formel 0.6767 (√2), 0.245823 gegen 0.2458 (b=1.2); Konstanten konvergieren logarithmisch langsam | neu |
| B5 | Rangverschiebung: H_4 → H_4∘ψ, ψ=x+ε sin2πx (oder x+ε(1−cos2πx)). Für die gestörte Leiter gilt H'_r = H_{r+c}+o(1) mit einer einzigen Konstante c(ε), gleich für p_r, λ_r und das Profil bei 0.25 und 0.5 | **numerisch** (4–5 Stellen) / **Vermutung** (allgemein) | √2: c(0.05)=−0.31022, c(−0.05)=+0.3067, c(0.01)=−0.06134, c_cos(0.01)=−0.07739; Gradient (dc/dε_sin, dc/dε_cos) ≈ (−6.13, −7.74) | neu |
| B6 | Die Querabweichungen (alles außer der Verschiebung) fallen pro Rang etwa um den Faktor μ_r ~ k/r, also überexponentiell (~k^r/r!) | **numerisch** (grob) | c_r−c_∞ (aus λ): 1.6e−2, 4.8e−3, 1.4e−3, 4.0e−4, 1.0e−4 (r=6…10); 1/μ_r = 3.5, 4.1, 4.8, 5.4 | neu |
| B7 | Daraus: Auf beschränkten Leitern gibt es ein Abel-Funktional Φ mit Φ(T f)=Φ(f)+1 (Rangkoordinate). Die θ-Freiheit einer Stufe wirkt asymptotisch nur über Φ. Bedingung für einen echten halben Rang: Φ(H_{3.5})=Φ(H_3)+½. Diese Bedingung legt **eine** Richtung der θ-Freiheit fest, nicht alle | **Vermutung** (Existenz von Φ) / **bewiesen** (die Zählung, falls Φ existiert und glatt ist) | B5, B6 | neu |
| B8 | Bedingter Satz: Ist H_s(b,·) für alle s≥4 konvex auf [1,2], dann gilt H_r(b,2) ≤ b+b(b−1)^{r−2}. Jede Basis b<2 hat dann ab einem Rang einen Fixpunkt, also η_r → 2 (mit kn.satz.schwellen_monoton) | **bewiesen** (bedingt) | §4 | vermutlich neu |
| B9 | Die Voraussetzung von B8 ist an den kritischen Basen ab N=5 **falsch**: F_N(1.5) liegt über der Sehne (N=5…8: +0.0048, +0.0122, +0.0171, +0.0206). F_N ist auf [0,2] konkav und erst danach konvex (S-Form) | **numerisch** (Profildaten rang_gamma, ~20 Stellen) | §4 | neu (negativ) |
| B10 | (4−e_r)/(2−η_r) = 2.308, 2.499, 2.660, 2.777, 2.849, 2.877 (r=3…8). Die Zuwächse fallen fast linear (0.191 → 0.029); vermutlich gibt es ein Maximum und keine Konstante 3 | **numerisch** / negativ | Werte aus eta_values.csv | — |
| B11 | Leiter auf negativen Höhen [−1,2] mit Tschebyschow Grad 60: zu langsam (nach 10 min erst Rang 7) | Fehlschlag | neg_lauf.py | — |

## 1. Die kanonische Leiter im beschränkten Bereich (B1)

Sei f=H_r(b,·) stetig, streng wachsend, f(1)=b>1, und habe f einen Fixpunkt ≥1. p_r sei der kleinste solche Fixpunkt, mit 0<λ_r=f′(p_r)<1. Auf [1,p_r) gilt f(y)>y, die Bahn von 1 steigt gegen p_r, und die Schröder-Funktion σ_r(y)=lim λ^{−n}(f^n(y)−p_r) ist reell-analytisch, streng wachsend und negativ. H_{r+1}(x)=σ_r^{−1}(λ_r^x σ_r(1)) erfüllt H_{r+1}(0)=1, H_{r+1}(x+1)=f(H_{r+1}(x)) und ist streng wachsend mit Grenzwert p_r. Das ist die Koenigs-Lösung, eindeutig unter den Lösungen mit H−p_r ~ C λ^x ohne periodischen Anteil; für b<e^{1/e} stimmt sie bei r=3 mit Kneser überein.

**Beweis von b<p_{r+1}<p_r.** Es gilt H_{r+1}<p_r, also auch p_{r+1}<p_r. Wegen H_{r+1}(b)>H_{r+1}(1)=b ist die Bahn 1, b, H_{r+1}(b), … ab dem zweiten Glied >b. Sie konvergiert gegen p_{r+1}, also p_{r+1}>b. Monotonie von q_r=H_{r+1}(2)=H_r(b): Es ist q_{r+1}=H_{r+2}(2)=H_{r+1}(b)=H_r(H_{r+1}(b−1)). Wegen b−1<1 gilt H_{r+1}(b−1)<b, also q_{r+1}<H_r(b)=q_r. Damit existiert p_∞=lim p_r ≥ b. ∎ **Offen:** p_∞=b (B2).

Numerisch (dps 30, Grad 40/50; Kontrolle der Funktionalgleichung ≤ 2·10⁻¹³; eine zweite Einstellung mit dps 20, Grad 30 stimmt auf 12 Stellen überein):

| r | p_r (b=√2) | λ_r | H_{r+1}(√2,0.5) | p_r (b=1.2) | p_r (b=1.5, Kneser-Saat) |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 0.693147 | 1.243622 | 1.257735 | — |
| 4 | 1.549217 | 0.208628 | 1.269410 | 1.215796 | 1.766710 |
| 6 | 1.459146 | 0.057083 | 1.296375 | 1.203615 | 1.589464 |
| 8 | 1.438372 | 0.025378 | 1.311410 | 1.201355 | 1.550993 |
| 10 | 1.429760 | 0.014014 | 1.321471 | 1.200636 | 1.534694 |
| 12 | 1.425236 | 0.008754 | 1.328846 | 1.200342 | 1.525854 |
| 15 | 1.421547 | 0.004982 | 1.336976 | 1.200158 | 1.518360 |

## 2. Grenzprofil min(x+1,b) — Herleitung (B2, B3)

Im Plateau [b,p_r] ist der Abstand zu p_r höchstens δ_r. Dort ist σ_r bis auf den relativen Fehler O(c₂δ_r)=O(L_rδ_r) die Identität (c₂=f″/(2λ(1−λ)), f″/λ ≈ −L_{r−1}). Daraus folgt H_{r+1}(1+x) ≈ p_r − δ_r λ_r^x. Der Fixpunkt und seine Ableitung ergeben die Modellabbildung B3. Setzt man das Plateau von H_r in die Funktionalgleichung H_{r+1}(x+1)=H_r(H_{r+1}(x)) ein, erhält man für 0≤x≤1:

  1 − λ_{r−1}^{H_{r+1}(x)−1} = (δ_r/δ_{r−1})·(1 − λ_r^x).

Mit μ_{r−1}=1−δ_r/δ_{r−1} ≈ λ_{r−1}^{b−1} wird daraus

  H_{r+1}(x) − 1 = −ln(λ_{r−1}^{b−1} + λ_r^x − …)/L_{r−1} → min(b−1, x),

sobald λ_r→0 und L_r/L_{r−1}→1. Gegenprobe der Profilformel (auswertung_leiter.py): bei b=1.2, Rang 15 gilt H(0.25)=1.146017 gegen die Formel 1.145963 und H(0.75)=1.198815 gegen 1.198720. Bei √2, Rang 15: 1.20277 gegen 1.20187. Die Formel liegt immer näher an den Daten. Die Daten selbst bewegen sich langsam (logarithmisch) auf min(x+1,b) zu.

Verträglich mit der bekannten Lage: Bei b=2 ist min(x+1,2)=x+1 auf [0,1] (kn.beob.zwischenwerte_basis2, H_r(2,½)→1.5). Und es passt zu kn.satz.eta_kriterium (F_N(q)→q+1 an den kritischen Basen).

## 3. Asymptotik der Modellabbildung (B4)

In Log-Variablen D=−ln δ, M=−ln μ (dann L=M/(b−1)) lautet das Modell M′=a(M+D−ln L), D′=D−ln(1−μ) mit a=b−1. M wird schnell an die langsame Kurve M=(D−ln(M/a))/k gebunden (Kontraktion a<1). Ableiten ergibt dM/dr=μ/(k+1/M), also d(1/μ)/dr = 1/(k+1/M). Daraus folgen r μ_r → k, δ_r ~ a k^k r^{−k}/ln r und L_r ~ ln r/a. Die Iteration des Modells bis r=10⁷ (modellabbildung.py, Start mit den echten Werten bei r=15) bestätigt die Zuwachsformel auf 4–5 Stellen. Die Konstanten r μ_r und δ_r r^k ln r nähern sich ihren Grenzwerten (k bzw. a k^k) nur mit Korrekturen der Größe ln ln r/ln r. Für die echte Leiter bis r=16 lassen sich die Exponenten deshalb nicht direkt ablesen (lokale Exponenten von δ: 1.82 bei √2, 3.46 bei 1.2, Ziel k+1/ln r). Belegt ist vor allem die Schritt-für-Schritt-Übereinstimmung mit dem Modell.

## 4. Konvexität und η_N → 2 (B8, B9)

**Satz (bedingt, bewiesen).** Sei 1<b<2, und für alle s≥4 sei H_s(b,·) stetig und konvex auf [1,2]. Dann gilt u_r:=H_r(b,2)=H_{r−1}(b,b) ≤ (2−b)b+(b−1)u_{r−1}, also u_r−b ≤ (b−1)^{r−3}(b²−b). Sobald b(b−1)^{r−2} ≤ 2−b, ist H_r(b,2) ≤ 2, also H_r(b,1)−1 > 0 ≥ H_r(b,2)−2, und der Zwischenwertsatz liefert einen Fixpunkt in [1,2]. Mit kn.satz.schwellen_monoton folgt η_r→2, mit der expliziten Schranke η_r ≥ Lösung von b+b(b−1)^{r−2}=2 (r=4: 1.545; r=8: 1.732). ∎

**Aber:** An den kritischen Basen ist F_N(1.5) − (F_N(1)+F_N(2))/2 = −0.0072 (N=4), +0.0048, +0.0122, +0.0171, +0.0206 (N=5…8). Ab N=5 ist F_N in der Mitte von [1,2] konkav. Auch auf [0,1] liegt F_N(0.5) über der Sehne. Die echten Stufen sind S-förmig, das passt zum Grenzprofil min(x+1,b). Der Konvexitätsweg zu TPID 10 ist damit **für die gewählten Fortsetzungen versperrt**. Die Schranke η_r ≥ η_r^* gilt nur für konvexe Fortsetzungsketten.

## 5. Rangverschiebung und Rangkoordinate (B5–B7)

Getestet wurde H_4(√2,·) ↦ H_4∘ψ mit ψ(x+1)=ψ(x)+1, ψ(k)=k, also dieselben ganzen Werte. Ab H_5 folgt wieder die kanonische Koenigs-Konstruktion. rangverschiebung.py bestimmt für jede Größe Q das c_r mit Q′(r)=Q(r+c_r) (Interpolation von ln|Q−Q_∞| in r).

| r | c_r aus p_r (ε=0.05) | aus λ_r | aus H(0.25) | aus H(0.5) |
|---:|---:|---:|---:|---:|
| 6 | −0.29226 | −0.22572 | −0.4391 | −0.3853 |
| 8 | −0.30950 | −0.30289 | −0.3275 | −0.3192 |
| 10 | −0.31033 | −0.30964 | −0.3124 | −0.3113 |
| 12 | −0.31028 | −0.31019 | −0.3105 | −0.3103 |
| 14 | −0.31022 | −0.31021 | −0.3103 | −0.3102 |

Alle vier Größen laufen auf **dieselbe** Zahl zu. Die gestörte Leiter ist also asymptotisch die ungestörte, nur um c Ränge verschoben. Gegenproben: ε=−0.05 ergibt +0.3067, ε=0.01 ergibt −0.06134. Damit ist c fast linear in ε, c ≈ −6.15ε − 0.7ε². Die Kosinusrichtung mit ε=0.01 ergibt −0.07739. Die Abweichungen der vierten Stelle bei r≥13 kommen aus der Randinterpolation.

**Deutung (Vermutung).** Im Funktionenraum besteht die beschränkte Leiter aus zwei Teilen. Quer zu einer eindimensionalen langsamen Kurve W wird sie überexponentiell stark angezogen (B6). Längs W bewegt sie sich parabolisch, mit Schritten ~k/r. Die Phase längs W ist ein Abel-Funktional Φ(Tf)=Φ(f)+1. Es ist das Funktionenraum-Gegenstück zur Fatou-Koordinate und die „Rangkoordinate“ einer Funktion. Für die θ-Frage folgt: Eine Fortsetzung, die „halb“ zwischen zwei Rängen liegt, muss Φ=Φ(H_3)+½ erfüllen. Das ist **eine** Gleichung. Sie fixiert eine Linearkombination der θ-Parameter (hier mit Gewichten ≈ (−6.13, −7.74) für (sin, 1−cos)), aber nicht die ganze θ-Freiheit. Für einen vollständig kanonischen halben Rang bräuchte man zusätzlich eine kanonische Interpolation der Querkontraktion. Die Faktoren μ_r ~ k/r legen hier eine Gamma-Interpolation ∏μ ↦ k^s/Γ(r+s)-artig nahe. Das ist nicht ausgearbeitet. Bei b=2 (nie beschränkt) gilt diese Analyse nicht direkt; der Rangabstieg R bleibt dort instabil (Etappe 6). Das passt dazu, dass R gerade die Querrichtungen überexponentiell aufbläht.

## Vorschläge für Grapheinträge

- „Beschränkte Rangleiter: Koenigs-Stufen fallen monoton, b<p_{r+1}<p_r“ — satz — bewiesen (Nichtentartung vorausgesetzt)
- „Grenzprofil der Rangleiter: H_r(b,x) → min(x+1,b)“ — vermutung — Vermutung mit Belegen (b=√2, 1.2, 1.5; r≤16)
- „Plateau-Modell der Rangleiter (δ, λ)-Abbildung“ — satz/heuristik — hergeleitet mit Lücke (Fehler O(Lδ))
- „Rate der Rangleiter: r·μ_r → k=(2−b)/(b−1), δ_r ~ (b−1)k^k/(r^k ln r)“ — vermutung — formal hergeleitet, Modell numerisch bis 10⁷
- „θ-Störung einer Stufe = Rangverschiebung der Leiter (Rangkoordinate Φ)“ — beobachtung/vermutung — numerisch, 4–5 Stellen
- „Querabweichungen der Leiter fallen überexponentiell (Faktor μ_r pro Rang)“ — beobachtung — numerisch grob
- „Konvexe Fortsetzungsketten haben η_r → 2, mit Schranke b+b(b−1)^{r−2}=2“ — satz — bewiesen (bedingt)
- „Die kritischen Stufen F_N sind ab N=5 auf [1,2] nicht konvex (S-Form)“ — beobachtung — numerisch (negativ für den Konvexitätsweg)

## Offene Fragen

1. Beweis von p_∞=b. Dafür bräuchte man eine quantitative Plateau-Abschätzung; die Lücke liegt in B3.
2. Existenz und Glattheit von Φ, und ob Φ unabhängig von der Messgröße exakt ist. Hier wäre eine Rechnung mit höherer Genauigkeit und mehr Rängen nötig (Grad und dps erhöhen; die Kosten wachsen, weil das Profil einen Knick bei x=b−1 entwickelt).
3. Grenzprofil auf negativen Höhen, x∈(b−2,0): Die Herleitung legt es dort nicht fest, und die Rechnung war zu langsam (B11).
4. Gibt es eine Rangverschiebung auch in der unbeschränkten Leiter bei b=2? Wenn ja, würde H_r(2,½) ≈ 1.5−0.62/r² eine θ-Abhängigkeit erst in der Ordnung r⁻³ tragen.
5. Kanonische Querinterpolation (Gamma-artig), um aus Φ einen vollen halben Rang zu machen.

## Reproduktion (in diesem Ordner, Python = C:/Users/janis/miniconda3/python.exe)

```
python koenigs_leiter.py "sqrt(2)" 16 30 40          # ~45 s
python koenigs_leiter.py "sqrt(2)" 9 20 30           # Gegenprobe mit anderer Einstellung
python koenigs_leiter.py 1.2 16 30 50                # ~40 s
python kneser_saat.py 1.5 16 30 40                   # Kneser-H4 vom gp-Motor (eigener cache/), ~60 s
python koenigs_leiter.py "sqrt(2)" 16 30 40 0.05     # theta-Störung (auch -0.05, 0.01, c0.01)
python auswertung_leiter.py daten/<datei>.json       # Modellfehler, Profilformel
python modellabbildung.py daten/<datei>.json         # Modell bis r=1e7 (~40 s)
python rangverschiebung.py daten/koenigs_bsqrt2_R16_d30_g40.json daten/koenigs_bsqrt2_R16_d30_g40_th0.05.json
```
Die Rohdaten liegen in `daten/`. Es wurden keine fremden Zwischenspeicher benutzt, und es gab keine Commits.
