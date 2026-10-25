# Runde 25.10.2026 — 23 Basisableitung: Vorzeichen und Formeln

## Kurz in einfacher Sprache

Wir haben untersucht, wie sich die Tetrationskurve T_b(x) verändert, wenn man die Basis b ein wenig grösser macht.
Ergebnis: Das Muster ist immer dasselbe, von b = 1.000001 bis b = 10 000 (regulär unter η, Kneser darüber). Zwischen −2 und −1 steigt die Kurve,
zwischen −1 und 0 sinkt sie, ab 0 steigt sie. Zwei Kurven verschiedener Basen schneiden sich darum nur bei x = −1 und x = 0.
Wir haben bewiesen, dass das ganze Muster nur an einem Stück hängt, nämlich an einer „Elastizität“ E auf dem Stück (−1, 0). Diese Elastizität
gehorcht einer sehr einfachen Stufenregel. An den Rändern hängt alles an einer Zahl ε₀(b): Wie stark wächst die Steigung T′(0), wenn ln b wächst?
ε₀ liegt immer zwischen 0 und 1. Sie fällt von 1 (bei b → 1, mit einer Formel) bis etwa 0.41 (bei b = 10⁴).
Ausserdem gilt bei b knapp über η ein einfaches Kanalgesetz für die Vorwärtssumme der Basisableitung: Sie wächst wie (b − η)^{−3/2}, mit einer Konstante,
die wir hergeleitet und auf 6·10⁻⁵ bestätigt haben.

## Befunde

Bezeichnungen: $D(x)=\partial_bT_b(x)$, $u=D/T'$, **Elastizität** $E(x)=\dfrac{b\ln b\,D(x)}{T(x)}=\dfrac{\partial\ln|T_b(x)|}{\partial\ln\ln b}$,
$\varepsilon_0(b)=\dfrac{b\ln b\,\partial_bT_b'(0)}{T_b'(0)}=\dfrac{\partial\ln T_b'(0)}{\partial\ln\ln b}$, $\delta=b-\eta$.

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | **Elastizitäts-Stufenregel:** $E(x+1)=\bigl(1+E(x)\bigr)\ln T(x+1)$, $E(0)=0$, $E(1)=\ln b$ (jede in b differenzierbare Familie mit $T_b(0)=1$) | bewiesen | Beweis unten; an Gitterdaten geprüft (b = 2, x = −½ → ½: 0.29086 = 0.29087) | vermutlich neu (Graph hat die lineare Gleichung für D, gr.satz.basis_gleichung, nicht diese Form) |
| B2 | **Reduktion:** Für eine reelle, in x steigende Familie mit $T(-1)=0$, $T(0)=1$ gilt das Muster „D > 0 auf (−2,−1) und (0,∞), D < 0 auf (−1,0)“ genau dann, wenn $\max(-1,\ln T(y))<E(y)<0$ für alle $y\in(-1,0)$ | bewiesen | Beweis unten | vermutlich neu |
| B3 | **Randwerte:** $E(-1^+)=\varepsilon_0-1$, $E'(0)=\varepsilon_0\,T'(0)$; das Muster mit transversalen Schnitten bei −1 und 0 verlangt $0<\varepsilon_0<1$, d. h. $T_b'(0)$ steigt in b und $T_b'(-1)=T_b'(0)/\ln b$ fällt in b | bewiesen | Beweis unten | vermutlich neu |
| B4 | Das Muster (B2) gilt für alle gemessenen Basen: 1.01, 1.1, 1.3, 1.44, 1.444 (regulär), η+10⁻³, 1.45, 2, e, 10, 1000, 10⁴ (Kneser); Gitter x = −1.95…2 (Schritt 0.05) | numerisch | Ränder: min(E+1) ≥ 0.42, min((E−ln T)/\|ln T\|) ≥ 0.23; D auf 20 Stellen (dps 34, Schritt 10⁻¹⁰); Gegenprobe dps 30 (Schritt 10⁻¹⁰, eigene Basisliste) gleiches Muster | Muster vermutlich neu für Kneser |
| B5 | **Monotonie in b:** $b\mapsto T_b(x)$ ist steigend für $x\in(-2,-1)\cup(0,\infty)$ und fallend für $x\in(-1,0)$, über 39 Basen von 1.01 bis 1000 (inkl. 1.444 / 1.445 beidseits η) | numerisch | alle 380 Nachbardifferenzen (10 Höhen) mit dem richtigen Vorzeichen | vermutlich neu |
| B6 | **Schnittpunkte:** zwei Kurven $T_{b_1},T_{b_2}$ ($b_1<b_2$) treffen sich auf (−2,∞) nur bei x = −1 und x = 0; $T_{b_2}<T_{b_1}$ genau auf (−1,0) | numerisch (folgt aus B4/B5) | wie B4/B5 | vermutlich neu |
| B7 | Folgerung slog: $\partial_b\operatorname{slog}_b(y)=-u(\operatorname{slog}_by)$ ist > 0 für $y\in(0,1)$ und < 0 für $y<0$ und $y>1$ | bewiesen aus B4 (bedingt) | Umrechnung $D+T'\,G\circ T=0$ | — |
| B8 | E ist auf (−1,0) streng steigend (von ε₀−1 nach 0) und bei fester Höhe in b fallend; z. B. E(−½) = −0.0454 (1.01), −0.190 (η), −0.230 (2), −0.311 (10), −0.405 (10⁴) | numerisch | 13 Basen, Gitter 0.05 | — |
| B9 | **ε₀(b)** (Tabelle unten) liegt für 1 + 10⁻⁶ ≤ b ≤ 10⁴ in (0,1), stetig und fallend durch η (ε₀(η) ≈ 0.61296) | numerisch | 12–17 Stellen (Spanne zweier Einstellungen 10⁻¹⁸ … 5·10⁻¹⁴) | neu |
| B10 | **b → 1⁺:** $\varepsilon_0=1-\dfrac1{\ln(1/\ln b)}+O(\ln b)$, $T_b'(0)\approx(a-1)\ln(1/\lambda)$ ($a$ Fixpunkt, $\lambda=a\ln b$) | hergeleitet mit Lücke | b = 1+10⁻⁶: 0.927619011 gegen 0.927617589, Rest ≈ 1.42·ln b; b = 1+10⁻⁴: Rest 1.38·ln b | vermutlich neu |
| B11 | Grosse Basen: ε₀ fällt fast linear in ln ln b, Steigung −0.063, −0.060, −0.058 (10→100→10³→10⁴), langsam flacher | numerisch | B9-Werte | offen, ob ε₀ → 0 oder → Grenzwert > 0 |
| B12 | **Vorwärtssumme + Wackeln:** $u=u_0+w$, $u_0(x)=\sum_{k\ge0}[r(k)-r(x+k)]$, $r(x)=\dfrac{T(x)}{b\ln b\,T'(x)}=\dfrac{1}{b\ln^2b\,T'(x-1)}$, $w$ 1-periodisch mit $w(0)=0$ (b > η) | bewiesen (Struktur, = Kobis Gleichung in u-Sprache) | Periodizität von w am Gitter: 6·10⁻¹⁴ … 7·10⁻²⁴ | Struktur bekannt (gr.satz.turm_pde); Messung neu |
| B13 | Kneser-Wackeln über die Basen: Spanne von w (Gitter 0.05) 1.512·10⁻³ (η+10⁻³), 1.500·10⁻³ (1.45), 5.66·10⁻⁴ (2), 2.85·10⁻⁴ (e), 3.23·10⁻⁵ (10), 1.0·10⁻⁷ (10³), 7.4·10⁻⁹ (10⁴); relativ zu \|u(−½)\|: 0.64 %, 0.64 %, 0.60 %, 0.59 %, 0.55 %, 0.51 %, 0.50 %. Für b → η⁺ bleibt w endlich | numerisch | b = 2 und e stimmen mit gr.beob.kneser_pde_wackeln (5.7·10⁻⁴, 2.9·10⁻⁴) | Basisabhängigkeit neu |
| B14 | Die Formel $D\approx T'\,u_0$ (ohne Kneser-spezifisches Wackeln) gibt die Kneser-Basisableitung überall auf < 0.7 % relativ und trifft das Vorzeichenmuster | numerisch | B13 | — |
| B15 | **Kanalgesetz:** $S_b(0)=\sum_{k\ge0}r(k)\sim\dfrac\pi2\sqrt{\dfrac{2\eta}{e}}\;\delta^{-3/2}$ ($=1.61946527\,\delta^{-3/2}$), Korrektur $1-c\,\delta$ mit $c\approx0.644$ | hergeleitet mit Lücke (Leitterm); Korrektur numerisch | Verhältnis 0.94432, 0.98151, 0.99365, 0.99808, 0.99936, 0.99981, 0.999936 bei δ = 10⁻¹, 3·10⁻², …, 10⁻⁴ | vermutlich neu |
| B16 | $D(x)\to+\infty$ für $x\to-2^+$ | bewiesen | $D(x)=[D(x+1)/T(x+1)-T(x)/b]/\ln b$, erster Term → $D'(-1)/T'(-1)$ endlich, $T(x)\to-\infty$ | — |
| N1 | Motor bei η+10⁻⁴: T′(0) = 0.6111659 plausibel, aber $\partial_bT'(0)$ falsch (ε₀ = 0.139); b = 10⁵ scheitert (invz0rfunc) | negativ | Trend der Nachbarbasen | Motorgrenze |
| N2 | Kein Beweis der inneren Ungleichungen $\max(-1,\ln T)<E<0$ für Kneser; Versuch über die Vorwärtssumme führt nur auf die gleichwertige Bedingung $R(x+1)<R(0)<R(x)$ (unten) | negativ | — | — |

## Beweise und Details

### B1 Elastizitäts-Stufenregel
$\ln T(x+1)=\ln b\cdot T(x)$. Wende $b\ln b\,\partial_b$ an. Links steht $b\ln b\,D(x+1)/T(x+1)=E(x+1)$. Rechts steht $b\ln b\,[T(x)/b+\ln b\,D(x)]=\ln b\,T(x)\,[1+b\ln b\,D(x)/T(x)]=\ln T(x+1)\,(1+E(x))$.
$E(0)=0$ wegen $D(0)=0$. $E(1)=\ln b\cdot D(1)=\ln b$ wegen $T_b(1)=b$. ∎
Bei ganzen Höhen ergibt das $E(n)=\sum_{k=1}^{n}\prod_{j=k}^{n}\ln T(j)$, die logarithmische Form von gr.satz.basis_gleichung.

### B2 Reduktion
Voraussetzung: T stetig, streng steigend auf (−2,∞), $T(-1)=0$, $T(0)=1$. Dann ist $T<0$ auf (−2,−1), $T\in(0,1)$ auf (−1,0), $T>1$ auf (0,∞).
* (−1,0): $\operatorname{sign}D=\operatorname{sign}E$, weil $T>0$.
* (0,1): Für $y\in(-1,0)$ ist $E(y+1)=(1+E(y))\ln T(y+1)$ mit $\ln T(y+1)>0$. Also $D(y+1)>0\iff E(y)>-1$.
* (1,∞): Ist $E(x)>0$, so ist $E(x+1)=(1+E(x))\ln T(x+1)>0$. Induktion über die Streifen.
* (−2,−1): Für $x=y-1$ ist $E(y)=(1+E(x))\ln T(y)$ mit $\ln T(y)<0$. Also $E(x)=E(y)/\ln T(y)-1$. Wegen $T(x)<0$ ist $\operatorname{sign}D(x)=-\operatorname{sign}E(x)$.
  Damit gilt $D(x)>0\iff E(y)/\ln T(y)<1\iff E(y)>\ln T(y)$.
* Ganze Stellen: $D(-1)=D(0)=0$, $D(n)>0$ für $n\ge1$ (Summenformel). ∎

### B3 Randwerte
$T'(0)=\ln b\,T'(-1)$ (Kettenregel bei x = −1, $T(0)=1$). Bei $y\to-1$ gilt nach l'Hôpital $D/T\to D'(-1)/T'(-1)$. Also
$E(-1^+)=b\ln b\,\partial_b\ln T'(-1)=b\ln b\,\partial_b[\ln T'(0)-\ln\ln b]=\varepsilon_0-1$.
$E'(0)=b\ln b\,[D'(0)T(0)-D(0)T'(0)]/T(0)^2=b\ln b\,D'(0)=\varepsilon_0T'(0)$.
Notwendigkeit: Nahe $0^-$ ist $E(y)\approx\varepsilon_0T'(0)\,y$ und $\ln T(y)\approx T'(0)\,y$. Aus $E<0$ folgt $\varepsilon_0\ge0$, aus $E>\ln T$ folgt $\varepsilon_0\le1$.
Nahe $-1^+$ gibt $-1\le E\le0$ dasselbe. ∎

### B10 b → 1⁺ (Lücke)
Für b → 1 rückt der anziehende Fixpunkt a an 1, und der Multiplikator $\lambda=a\ln b$ geht gegen 0. Die Iteration ist dann bis auf $O(\ln b)$ linear:
$T(x)\approx a-(a-1)\lambda^x$. Daraus folgt $T'(0)\approx(a-1)\ln(1/\lambda)$. Mit $a-1\approx\lambda\approx\ln b$ ergibt sich
$\ln T'(0)\approx\ln\ln b+\ln\ln(1/\ln b)$, also $\varepsilon_0\approx1-1/\ln(1/\ln b)$.
Gemessen: T′(0) bei 1+10⁻⁶ ist 1.38155238735·10⁻⁵, die Formel gibt 1.38155·10⁻⁵. **Lücke:** Der Fehler der Linearisierung (Koenigs-Koordinate) ist nur geschätzt.
Der Rest $(\varepsilon_0-\text{Formel})/\ln b$ wandert noch: 1.27, 1.34, 1.38, 1.42. Vermutlich steckt dort ein Faktor mit ln ln.
Nebenbei: Bei festem $y\in(-1,0)$ gilt $T_b(y)\to1$ (Stufe bei −1). Das passt zu Runde 18 B4 ($T_b(1/2)\approx1+\varepsilon-\varepsilon^{3/2}$).

### B12/B13 Vorwärtssumme
Aus gr.satz.basis_gleichung: $u(x+1)-u(x)=r(x)$. Für b > η wächst $T'(x+k)$ übergeometrisch. Darum konvergiert $u_0$, und es gilt $u_0(x+1)-u_0(x)=r(x)$.
Also ist $w=u-u_0$ 1-periodisch, und $w(0)=0$ wegen $u(0)=u_0(0)=0$.
Das ist Kobis Gleichung (gr.satz.turm_pde): $w(x)=p(0)-p(x)$. Für die Kobi-Familie wäre $w\equiv0$.
$S_b(0)=\sum_kr(k)$ braucht nur T′(0) und die ganzen Türme: $S_b(0)=\frac{1}{b\ln^2b}\bigl[\frac{\ln b}{T'(0)}+\sum_{j\ge0}\frac1{T'(j)}\bigr]$ mit $T'(j)=T'(0)\prod_{i=1}^{j}\ln b\,(b\uparrow\uparrow i)$.
Für die Kobi-Familie ist das Vorzeichenmuster gleichwertig mit $R(x+1)<R(0)<R(x)$ für $x\in(-1,0)$, wobei $R(x)=\sum_{k\ge-1}1/T'(x+k)$ und $R(-1)=R(0)$.
Mit dieser Form kam ich nicht weiter (N2).

### B15 Kanalgesetz (Herleitung, Lücke)
Knapp über η läuft die Bahn lange durch den Engpass bei e. Dort ist $q(y)=b^y-y\approx\alpha+\beta(y-e)^2$ mit $\alpha=e^2\delta/\eta$ (aus $\ln b=1/e+\delta/\eta+O(\delta^2)$) und $\beta=1/(2e)$.
Im Kanal gilt $T'\approx q(T)$. Damit wird
$\sum_j 1/T'(j)\approx\int dx/T'=\int dy/q(y)^2=\pi/(2\alpha^{3/2}\beta^{1/2})$. Mit $1/(b\ln^2b)\to e^2/\eta$ ergibt sich $S\sim\frac\pi2\sqrt{2\eta/e}\;\delta^{-3/2}$.
**Lücke:** Wir ersetzen T′ durch q, die Summe durch das Integral und lassen die Glieder ausserhalb des Kanals weg. Die Messung zeigt, dass die Korrektur von der Ordnung δ ist, nicht δ^{1/2}.
$(1-\text{Verhältnis})/\delta$: 0.6432 (10⁻³), 0.6436 (3·10⁻⁴), 0.6436 (10⁻⁴). Die Konstante c ≈ 0.644 ist nicht identifiziert.

### Tabelle ε₀(b) und T′_b(0) (Auszug; alle Werte in `a4_eps0_tabelle.py`)
| b | ε₀ | T′(0) |
|---|---|---|
| 1+10⁻⁶ | 0.927619011350656 | 1.38155238735013·10⁻⁵ |
| 1+10⁻⁴ | 0.891564964316508 | 9.21121157816745·10⁻⁴ |
| 1.01 | 0.795746922025197 | 0.046470364663981 |
| 1.1 | 0.694300150872462 | 0.252411832369148 |
| 1.3 | 0.634409655317119 | 0.494931191358990 |
| 1.444 | 0.613038361675502 | 0.610624490849637 |
| η+10⁻³ | 0.612836763472541 | 0.611799684319319 |
| 1.45 | 0.612316944671937 | 0.614837428629744 |
| 2 | 0.571578377725955 | 0.889364954620976 |
| e | 0.547308383994694 | 1.09176735125832 |
| 10 | 0.492576113765722 | 1.68428902795601 |
| 100 | 0.448962383770609 | 2.33385048951909 |
| 1000 | 0.424696207332254 | 2.78599024341557 |
| 10⁴ | 0.408139790749556 | 3.14051125641320 |

T′(0) bei e trifft den bekannten Kneser-Wert 1.0917673512583. Bei 2, 100, 1000 und 10⁴ trifft er die Werte aus Runde 18 (B10).

## Reproduktion (Python: miniconda, `PYTHONIOENCODING=utf-8`, im Ordner ausführen; Motor `motor.py` mit eigenem Speicher `cache/`)
- Gitterdaten: `d2_elastizitaet.py 34 10 10 fein 101/100 11/10 13/10 36/25 29/20 2 "exp(1)" 10 1000 10^4`, `... fein 361/250`, `... fein "exp(exp(-1))+1/1000"`.
- ε₀ (Rand): `d2_elastizitaet.py 40 12 12 rand ...` (Basislisten siehe Dateinamen in `cache/`), Gegenprobe `d2_elastizitaet.py 30 9 9 rand 6/5 2 10`; b → 1: `d2_elastizitaet.py 40 14 12 rand "1+1/10^4" "1+1/10^6"`.
- B4, B8, B12–B14: `a2_muster.py cache/d2_fein_....json`. B9–B11: `a4_eps0_tabelle.py`. B15: `a3_kanal.py cache/d2_rand_40_12_12_expexp-1+1o10_*.json`.
- B5: `d3_dicht.py 30`. Erste Gegenprobe (dps 30, Dezimalbasen): `d1_vorzeichen.py 30 1e-10 ...` und `a1_auswertung.py`.
- Jeder Lauf dauert unter 1 Minute.

## Vorschläge für Grapheinträge
- gr.satz.elastizitaet_stufe — „Die Basis-Elastizität E = ∂ln T/∂ln ln b gehorcht E(x+1) = (1+E(x)) ln T(x+1)“ — satz — bewiesen
- gr.satz.vorzeichen_reduktion — „Das Vorzeichenmuster der Basisableitung hängt nur an max(−1, ln T) < E < 0 auf (−1,0)“ — satz — bewiesen
- gr.satz.elastizitaet_rand — „E(−1⁺) = ε₀ − 1, E′(0) = ε₀T′(0); das Muster verlangt 0 ≤ ε₀ ≤ 1“ — satz — bewiesen
- kn.beob.basis_monoton — „T_b(x) ist in b steigend auf (−2,−1) und (0,∞), fallend auf (−1,0); Kurven verschiedener Basen schneiden sich nur bei −1 und 0“ — beobachtung — numerisch
- kn.beob.eps0_tabelle — „ε₀(b) fällt von 1 (b → 1) über 0.613 (η) bis 0.408 (10⁴)“ — beobachtung — numerisch
- kl.satz.eps0_b_gegen_1 — „ε₀ = 1 − 1/ln(1/ln b) + O(ln b) für b → 1⁺“ — satz — skizze (Lücke: Linearisierung)
- gr.satz.kanalgesetz_vorwaertssumme — „Σ r(k) ~ (π/2)√(2η/e)(b−η)^{−3/2}“ — satz — skizze (Lücke: Kanal-Näherung), Korrektur numerisch
- gr.beob.wackeln_basen — „Das Kneser-Wackeln der Basisableitung beträgt 0.5–0.64 % von u, bleibt für b → η⁺ endlich“ — beobachtung — numerisch

## Offene Fragen
1. Beweis von $\max(-1,\ln T)<E<0$ auf (−1,0) für Kneser, etwa über die Fatou-Koordinate oder für die Kobi-Familie über $R(x+1)<R(0)<R(x)$.
2. Wie verhält sich ε₀ für b → ∞? Bleibt es positiv? Fiele es linear in ln ln b weiter, würde es nahe ln ln b ≈ 9 null, und das Muster bräche bei x = 0.
   Dafür bräuchte der Motor Basen über 10⁵.
3. Konstante c ≈ 0.644 in der Korrektur des Kanalgesetzes; endlicher Grenzwert des Wackelns für b → η⁺ (≈ 1.51·10⁻³ Spanne); hängt er mit der Écalle-Invariante von $e^u-1$ zusammen?
4. Gilt ein Mittelwertsatz für $\partial_b\theta_b$ (Mittel über eine Periode = $\partial_b\theta_\infty$, wie bei den Mittelwert-Gesetzen aus Runde 18/09)? Damit hätte man eine geschlossene Formel für das Mittel des Wackelns. Nicht gerechnet.
5. Das Muster für komplexe Höhen und die zweite Basisableitung $\partial_b^2T$ (Krümmung in b) wurden nicht untersucht.
