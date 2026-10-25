# Runde 25.10.2026 — 12 Höhere Ränge als Funktion der Basis

## Kurz in einfacher Sprache

Jede höhere Rechenart wird an einem Fixpunkt der Rechenart darunter gebaut. Die ungeraden Ränge (Pentation, Heptation, …) sitzen an einem
**negativen** Fixpunkt, den es für jede Basis b > 1 gibt; die geraden Ränge (Hexation, Oktation, …) an einem **positiven**, der bei einer
Grenzbasis verschwindet. Daraus folgt eine einfache Regel: Nur die Grenzbasen mit ungeradem Index (η₃, η₅, η₇, …) sind Stellen, an denen
sich die Bauweise ändert; an η₄ passiert mit der Pentation nichts (geprüft: sie ist dort auf 19 Stellen analytisch in b).
An η₅ dagegen passen die Hexation von links (regulär) und von rechts (Kneser) **nicht** zusammen: sie unterscheiden sich um eine kleine,
genau periodische Höhenverschiebung θ(x) mit Grösse etwa 10⁻⁵ — bei der Tetration an η₃ gibt es keinen solchen Sprung.
Dazu eine Landkarte der negativen Fixpunkte und Singularitäten (Ränge 4–6) für b von 1.2 bis 3.

Notation: H₃ = bˣ, H₄ = Tetration, H₅ = Pentation, H₆ = Hexation; H_{r+1}(b,0) = 1, H_{r+1}(b,x+1) = H_r(b,H_{r+1}(b,x)); ηᵣ = Basis, an der Hᵣ die Diagonale berührt
(η₃ = e^{1/e}, η₄ = 1.6353245, η₅ = 1.7373559). Werte aus dem gp-Motor (fatou_fork, `hyper(r,x)`), Zugang `motor.py` mit eigenem Zwischenspeicher.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | Ganzzahlige Werte links: $H_r(b,-k)=1-k$ für $r\ge4$, $1\le k\le r-3$, jede Basis | bewiesen (Injektivität von $H_{r-1}$ vorausgesetzt) | Induktion unten; Motor: $H_6(-3)+2=0$ bei allen 6 Basen | wohl Folklore |
| B2 | **Lagegesetz der negativen Struktur.** Gerade $r\ge4$: $H_r(b,\cdot)$ lebt auf $(\sigma_r,\infty)$ mit $H_r\to-\infty$ bei $\sigma_r=H_r^{-1}(p_{r-2})-1\in(-(r-2),-(r-3))$ ($\sigma_4=-2$, $p_2:=0$) und hat einen negativen Fixpunkt $p_r\in(\sigma_r,-(r-3))$. Ungerade $r\ge5$: $H_r:\mathbb R\to(p_{r-1},\infty)$ und $H_r(x)>x$ für alle $x\le0$ — kein negativer Fixpunkt | bewiesen unter Monotonie (Lücke: Monotonie der gebauten Stufen nicht bewiesen) | Beweis unten; Motor bestätigt bei b = 1.2, 1.3, 1.5, 1.7, 2, 3 ($\sigma_6<p_6<-3$; $\min_{[-12,0]}(H_5-x)\ge0.90$) | Lage von σ_r vermutlich neu |
| B3 | **Gerade–ungerade-Gesetz der singulären Basen:** Die Bauweise wechselt nur an $\eta_3,\eta_5,\eta_7,\dots$ (gerade Ränge verlieren dort ihren anziehenden Bau-Fixpunkt); die ungeraden Ränge sind an $p_{r-1}$ gebaut, der für alle $b>1$ existiert und einfach ist. Also ist $H_r(b,x)$ höchstens an den $\eta_{2j+1}<\eta_{r-1}$ singulär in $b$; $\eta_4,\eta_6,\dots$ sind Schwellen der Beschränktheit, aber keine singulären Basen | hergeleitet mit Lücke | Lücke: Einfachheit von $p_r$ für alle b nur numerisch (B6), Analytizität der Kneser-Stufen in b angenommen | vermutlich neu |
| B4 | **Pentation ist an η₄ analytisch in b:** Cauchy-Koeffizienten $c_k$ von $b\mapsto H_5(b,x)$ um η₄ auf Kreisen mit Radius 0.03 und 0.06 (16 Punkte, untere Hälfte per Spiegelung) stimmen überein | numerisch | x = 0.5, 1.5, −0.5: $c_0\ldots c_4$ gleich auf $10^{-19}$–$10^{-20}$, $c_7$ auf $10^{-15}$ (Faltungsrest); Gegenprobe Einseiten-Fit: Sprung $3\cdot10^{-12}$ = Fit-Rest | neu (als Messung) |
| B5 | **Sprung der Hexation an η₅:** Grenzwerte $L(x)=\lim_{b\uparrow\eta_5}H_6$ (regulär) und $R(x)=\lim_{b\downarrow\eta_5}H_6$ (Kneser, Zwei-Bogen-Auswahl) verschieden: $R(x)=L(x+\theta(x))$ mit **1-periodischem** θ, θ(k) = 0, $\theta(\tfrac14)=1.247\cdot10^{-5}$, $\theta(\tfrac12)=8.282\cdot10^{-6}$, $\theta(\tfrac34)=-4.185\cdot10^{-6}$; fast reine erste Harmonische $\theta\approx8.33\cdot10^{-6}\sin2\pi x+4.14\cdot10^{-6}(1-\cos2\pi x)$ | numerisch (Deutung offen) | Sprung bei x = 0.5: $-6.10\cdot10^{-6}$ mit Abständen 0.01–0.04 **und** 0.0025–0.01 (Fit-Abbruch würde um Faktor 256 fallen); Kontrollen ohne Bauwechsel: $3\cdot10^{-9}$ (H₆ um 1.66), $3\cdot10^{-12}$ (H₅ um η₄); Periodizität: θ(−½), θ(½), θ(3/2) gleich auf 4 Stellen, θ(¼)=θ(5/4), θ(¾)=θ(7/4) | neu |
| B6 | Landkarte der Bau-Fixpunkte (Tabelle unten): $p_4,\sigma_6,p_6$ fallen monoton in b; $\lambda_4(b)=T_b'(p_4)$ hat ein inneres Minimum bei $b^*=1.5327874$ (±4·10⁻⁷), $\lambda_4(b^*)=4.7798524$; $\lambda_6$ fällt von 12.7 (b = 1.2) auf 9.49 (b = 2) und steigt wieder (11.46 bei b = 3) | numerisch | $p_6(2)=-3.69940781692482$, $\lambda_6(2)=9.4930871143$ = Werte der Heptations-Konstruktion im Graphen (14 bzw. 11 Stellen); λ₄-Minimum nur eine Präzision (dps 24) | neu |
| B7 | Die Basisabhängigkeit des Bau-Fixpunkts der ungeraden Ränge entartet nie: $\inf_{b>1}\lambda_4(b)\approx4.78>1$ (auf dem Raster 1.2–3; Ränder: λ₄ wächst zu b → 1 und b → ∞) | Vermutung | Raster 1.2–3 + Minimumsuche; Ränder nicht bewiesen | — |
| B8 | Der Test „Tetration einseitig an η₃“ mit demselben Fit lief nicht: der Motor wirft für $b<\eta_3$ in diesem Aufruf einen gp-Fehler | negativ | e4c, Rang 4 | — |

### Landkarte (dps 24, `daten/struktur.json`)

| b | p₄ (Fixpunkt der Tetration) | λ₄ | σ₆ (Pol der Hexation) | p₆ (Fixpunkt der Hexation) | λ₆ | Bau der Hexation |
|---:|---:|---:|---:|---:|---:|---|
| 1.2 | −1.41090711782106 | 5.35928758 | −3.46432484986 | −3.41709586535 | 12.6942156 | regulär |
| 1.3 | −1.48119252125389 | 5.00162503 | −3.53074409403 | −3.47523583313 | 11.2856765 | regulär |
| 1.5 | −1.58672927938995 | 4.78300576 | −3.62649163343 | −3.56178274236 | 10.0075471 | regulär |
| 1.7 | −1.66313572515115 | 4.84481697 | −3.69421687256 | −3.62660329249 | 9.5445897 | regulär |
| 2 | −1.74390917613274 | 5.17179726 | −3.76502205053 | −3.69940781692 | 9.4930871 | Kneser |
| 3 | −1.87450725581019 | 7.07353609 | −3.87995816466 | −3.83415305722 | 11.4559004 | Kneser |

Abstand $p_6-\sigma_6$: 0.047, 0.055, 0.065, 0.068, 0.066, 0.046 — nicht monoton. Die Werte gehen ohne Knick über η₅ = 1.737 (zwischen 1.7 und 2).

## Beweise und Details

**B1.** $H_4(-1)=\log_b1=0$. Gilt $H_r(-k)=1-k$ für $k\le r-3$, so $H_{r+1}(-1)=H_r^{-1}(1)=0$ und $H_{r+1}(-k-1)=H_r^{-1}(H_{r+1}(-k))=H_r^{-1}(1-k)=-k$ für $k\le r-3$, also bis $k+1=r-2$. ∎

**B2.** Induktion, Anfang $r=4$: $H_3=b^x$ bildet $\mathbb R$ auf $(0,\infty)$ ab, also lebt $H_4$ auf $(-2,\infty)$; $T(x)-x\to-\infty$ bei −2, $T(-1)+1=1>0$ gibt $p_4\in(-2,-1)$.
*Ungerades r* (gebaut regulär am abstossenden $p_{r-1}\in(-(r-3),-(r-4))$, $H_r\downarrow p_{r-1}$ für $x\to-\infty$): Für $x\le-(r-3)$ ist $H_r(x)>p_{r-1}>-(r-3)\ge x$;
für $-k-1<x\le-k$ mit $k+1\le r-3$ ist $H_r(x)>H_r(-k-1)=-k\ge x$ (B1, Monotonie). Also $H_r(x)>x$ auf $(-\infty,0]$.
*Gerades r:* $H_{r-1}^{-1}$ ist nur auf $(p_{r-2},\infty)$ erklärt und strebt dort am linken Rand gegen −∞. Mit B1 ist $H_r(-(r-4))=-(r-5)>p_{r-2}>-(r-4)=H_r(-(r-3))$,
also gibt es genau ein $y\in(-(r-3),-(r-4))$ mit $H_r(y)=p_{r-2}$, und $\sigma_r=y-1\in(-(r-2),-(r-3))$ mit $H_r(\sigma_r^+)=-\infty$. Da $H_r(-(r-3))+(r-3)=1>0$, liegt ein Fixpunkt in $(\sigma_r,-(r-3))\subset(-(r-2),-(r-3))$ — genau das Intervall, das der nächste ungerade Schritt braucht. ∎
Lücke: Monotonie der gebauten Stufen in x; Eindeutigkeit des negativen Fixpunkts ist nur numerisch (r = 6, sechs Basen, Abstoss λ₆ ≥ 9.49).

**B4.** `e2_cauchy_h5.py 1.6353244967152763993 16 20 0.03 0.06`. Der Motor scheitert bei Im b < 0 in `sexpinit`; die untere Hälfte kommt aus $H_5(\bar b)=\overline{H_5(b)}$.
Stimmen die Koeffizienten zweier Kreise überein, sind die Kreiswerte (gerechnet + gespiegelt) Randwerte einer einzigen analytischen Funktion — das prüft die Spiegelung mit.
$c_k$ bei x = 0.5: 1.36021836220, 0.36612344974, −0.17305719728, 0.10786705950, −0.07798323078, … (Datei `daten/cauchy_h5_*.json`).

**B5.** `e4_h6_ueber_eta5.py`, `e4c_stencil.py`, `e5_theta_sprung.py`: kubische Einseiten-Interpolation in $t=b-\eta_5$ aus vier Basen je Seite; θ = (R−L)/L′ mit L′ aus zentraler Differenz der linken Grenzwerte.
Der Motor meldet 33.6 (links) bzw. 34.2 Stellen (rechts) für H₆ bei η₅ ∓ 0.01; dps 20 und 32 liefern dieselben 22 Stellen. Beide Seiten sind für sich glatt in b (die zwei Stencils geben dieselben einseitigen Grenzwerte auf 2·10⁻⁸).
**Deutung offen:** Zwei Lösungen derselben Abel-Gleichung unterscheiden sich immer um ein periodisches θ; der Befund zeigt also, dass die Kneser-Auswahl des Motors für b ↓ η₅ **nicht** zur regulären Hexation an η₅ strebt. Ob das für die wahre Kneser-Hexation gilt oder an der Zwei-Bogen-Auswahl nahe dem fast parabolischen Fixpunkt liegt, ist nicht entschieden. Bei der Tetration an η₃ fand Runde 18/03 keinen Sprung (Kneser = Fortsetzung der regulären Familie auf 10⁻³³).

**B6.** `e1_struktur.py 24 <basen>` (Sekante bzw. Illinois mit Klammer; die erste Fassung mit Sekante sprang bei b = 1.3 auf einen Wert links von σ₆ — Ast jenseits des Pols, verworfen), `e3_lambda4_min.py 24` (Raster 1.45–1.85, Goldener Schnitt).

## Vorschläge für Grapheinträge

- kn.satz.negative_struktur — Lage von σ_r und p_r, kein negativer Fixpunkt ungerader Ränge — satz — bewiesen unter Monotonie
- kn.satz.ganze_negative_hoehen — H_r(b,−k) = 1−k für k ≤ r−3 — satz — bewiesen
- kn.verm.singulaere_basen_ungerade — nur η₃, η₅, η₇, … sind singuläre Basen der Hierarchie — vermutung/hergeleitet mit Lücke
- kn.beob.pentation_analytisch_eta4 — H₅(b,x) analytisch in b an η₄ (19 Stellen) — beobachtung — numerisch
- kn.beob.hexation_sprung_eta5 — Kneser- und reguläre Hexation an η₅ um periodisches θ ~ 10⁻⁵ verschieden — beobachtung — numerisch, Deutung offen
- kn.beob.landkarte_p4_p6 — p₄, σ₆, p₆, λ₄, λ₆ für b = 1.2…3, λ₄-Minimum bei b* = 1.53279 — beobachtung — numerisch

## Offene Fragen

1. Ist der Sprung an η₅ echt (wahre Kneser-Hexation) oder eine Eigenschaft der Zwei-Bogen-Auswahl? Prüfweg: Kneser-Hexation bei b = η₅ + 0.0025 mit einem zweiten, unabhängigen Verfahren; Verhalten der Amplitude von θ für b ↓ η₅ (die Fits zeigen keinen Trend, θ scheint einen festen Grenzwert zu haben).
2. Die positiven Fixpunkte der Pentation nahe η₅ (Faltengesetz $a_5,c_5\approx e_5\mp\sqrt{\alpha(\eta_5-b)/\beta}$) wurden hier nicht nachgerechnet.
3. Hat $\lambda_6(b)$ ein Minimum nahe b ≈ 2, und ist $b^*_4=1.5327874$ mit etwas Bekanntem verbunden? (Zweite Präzision fehlt.)
4. Grenzverhalten für b → 1⁺ und b → ∞: $p_4\to-1$ bzw. $-2$, $\sigma_6\to-3$ bzw. $-4$ (aus B1/B2 naheliegend, nicht gerechnet).
