# Runde 5 — Etappe 1: Pivot-Grenzwert und Modellfamilien

**11.10.2026, Strang G. Nur (H₀); K2 bleibt bei P. Zwischenstand nach Etappe 1.**

## 1. Ergebnis und Reichweite

**Bewiesen:** vollständige Determinanten- und Schurformeln für reine Komposition; vollständige starke Determinantenasymptotik in einem Dreiecksmodell; affine Komposition mit konstantem M. Der einfachste nichttriviale Randfixpunktfall lässt auch Q, A und 𝒟 exakt bestimmen.

**Numerisch:** neun reale Familien, jeder Grad 25…400, wachsende Präzision 1200+n Bit. Eine Störung vom Grad 4 lässt L und μ exakt unverändert und verändert den extrapolierten Grenzwert deutlich. Vier direkt ausprobierte geometrische Mittel/Kapazitätsformeln passen nicht. Die kubische Grenzzahl wird auf etwa 1,1737972269, Kneser auf etwa 1,1737312358 extrapoliert. Dies sind **Fits, keine bewiesenen oder zertifizierten Grenzwerte**.

**Offen:** eine geschlossene, bewiesene Formel für λ des ursprünglichen kubischen/Kneser-Operators. Ein Quotient konformer Radien auf der geklebten Fläche wird unten nur als **bedingter Ansatz**, nicht als Ergebnis ausgegeben.

## 2. Normierung und Determinanten

Wir arbeiten in z=x/ρ mit F(z)=f(ρz)/ρ, M(z)=m(ρz) und
J_n=([z^k](F(z)^j+M(z)z^j))_{0≤k,j≤n}. D_n=det J_n, D_{−1}=1, h_n=D_n/D_{n−1}. Radiusgewichtung ist eine diagonale Ähnlichkeit und ändert die Determinante nicht. Ein Kreisradius hat daher keine eigenständige Determinantenwirkung; die Abbildung F und die Abschneidestelle 0 bleiben entscheidend.

**Lemma G25 [bewiesen].** Sind alle h_n positiv und h_n/h_{n−1}→λ>0, dann log h_n=n log λ+o(n) und log D_n=(n²/2)log λ+o(n²). **Die Umkehrung gilt nicht.**

*Beweis:* Die erste Aussage folgt durch Summation von log(h_j/h_{j−1}) und Cesàro, die zweite durch erneute Summation. Für die Umkehrung setze h_n=λ^n exp(√n sin n); die zusätzlichen Summen sind jedenfalls O(n^{3/2})=o(n²), während die Quotienten nicht konvergieren. ∎

## 3. Reine Komposition: ein vollständiger Satz

**Satz G26 [bewiesen].** Sei F holomorph nahe 0, a=F(0), b=F′(0)≠0. Bei M=0 gilt für alle n

\[
D_n=b^{n(n+1)/2},\qquad h_n=b^n,\qquad h_n/h_{n-1}=b,\qquad Q_n(z)=(z-a)^n.
\]

Ist G(z)=F(z)−a und H dessen lokaler Inversenzweig bei 0, dann
\[
A_{nk}=b^n[z^n]H(z)^k,\qquad
\mathcal D_n(w)=\|Q_n\|_2[z^n]\frac1{1-wH(z)},\quad w=\rho/x.
\]

Alle Aussagen betreffen endliche Jets; F muss die Einheitskreisscheibe nicht in sich abbilden. Für komplexes b ist der **komplexe** Pivotquotient b; die entsprechende Modulrate ist |b|.

*Beweis:* Entwickle u(a+G(z)). Die Translation u(s)↦u(a+s) ist auf P_n dreieckig mit Diagonale 1; die Jetmatrix der Komposition mit G ist unterdreieckig mit Diagonale b^j. Damit folgen Determinante und Pivot. (F−a)^n=b^nz^n+O(z^{n+1}) beweist Q. Die lokale inverse Jetkomposition hat Matrix ([z^n]H(z)^k); die abschließende Rücktranslation ändert den höchsten Koeffizienten nicht. Somit ist die letzte Zeile der Inversen genau diese Zeile. Summation der geometrischen Reihe beweist die 𝒟-Formel. ∎

**Wichtige Konsequenz:** In diesem Modell ist λ vollständig lokal an **0**, nicht am Randfixpunkt. Eine Formel nur aus einer Bildkapazität kann daher nicht ohne weitere Hypothesen für die ganze Operatorklasse gelten.

## 4. Dreiecksmodell mit M: starker Szegő-Typ vollständig bewiesen

**Satz G27 [bewiesen].** Sei F(0)=0, b=F′(0), m=M(0). Dann ist J_n unterdreieckig,
D_n=∏_{j=0}^n(b^j+m), h_n=b^n+m, sofern die nötigen Hauptminoren nicht null sind. Insbesondere für b>1 und m≥0 gilt
\[
D_n=b^{n(n+1)/2}\,C(b,m)\,\exp(-R_n),\quad
C(b,m)=\prod_{j=0}^\infty(1+mb^{-j}),\quad
0\le R_n\le \frac{m b^{-n-1}}{1-b^{-1}}.
\]
\[
\frac{h_n}{h_{n-1}}=b-\frac{(b-1)m b^{-(n-1)}}{1+m b^{-(n-1)}}.
\]

*Beweis:* Jede Spalte F^j+Mz^j beginnt bei Ordnung j. Die Diagonale ist b^j+m. Das unendliche Produkt konvergiert, da ∑mb^{−j}<∞; log(1+t)≤t kontrolliert den Schwanz. Die Quotientenformel folgt durch Division. ∎

Für 0<b<1 und m>0 gilt ebenfalls vollständig explizit D_n=m^{n+1}C_<(b,m)exp(−R_n), mit C_<=∏_{j≥0}(1+b^j/m), 0≤R_n≤b^{n+1}/[m(1−b)], h_n=b^n+m und λ=1. Ein **quadratisches Modell mit repulsivem Randfixpunkt** ist F(z)=(z+z²)/2, M=F′∘F=1/2+F: F(1)=1 und F′(1)=3/2, trotzdem λ=1. Hier sind h_n=2^{−n}+1/2 und der obige Produktschwanz exakt kontrolliert. Dies zeigt im einfachsten Quadratikmodell, dass der Randmultiplikator allein λ nicht bestimmt.

Das liefert explizite λ=b und **exponentielle** nächste Ordnung. Im vollständig diagonalen Teilmodell F=bz, M=m sind Q=A=e_n und 𝒟_n(w)=w^n/(b^n+m). Auf |w|≤1 fällt die Funktion exponentiell. Ein nichtverschwindender √n-Grenzvorfaktor ist also kein allgemeines Merkmal der Klasse.

## 5. Affiner Randfixpunkt und eine exakte Warnung zum Vorfaktor

**Satz G28 [bewiesen].** Für F(z)=a+bz, b≠1, M=m konstant, p=a/(1−b), sind (z−p)^j Eigenpolynome mit Eigenwerten b^j+m. Daher D_n=∏(b^j+m), h_n=b^n+m und A_n=e_n. Genauer
\[
Q_n(z)=(b^n+m)\sum_{j=0}^n\binom nj\frac{p^{n-j}(z-p)^j}{b^j+m}.
\]

*Beweis:* F(z)−p=b(z−p); entwickle z^n in der verschobenen Eigenbasis und löse J_nQ_n=h_ne_n. Die letzte Zeile der oberen Dreiecksinversen ist e_n^T/(b^n+m). ∎

**Korollar G29 [bewiesen, explizites einfachstes Randmodell].** Für F(z)=2z−1, M=0 liegt der Fixpunkt p=1 auf dem Abschneidekreis. Es gilt h_n=2^n, Q_n=(z+1)^n, A_n=e_n und
\[
\mathcal D_n(w)=2^{-n}\sqrt{\binom{2n}{n}}\,w^n.
\]

Insbesondere bei w=1:
\[
\log|\mathcal D_n(1)|=-\tfrac14\log(\pi n)-\frac1{16n}+\varepsilon_n,\qquad
|\varepsilon_n|\le\frac1{360n^3}\quad(n\ge1).
\]

*Beweis:* Die Koeffizientennorm folgt aus ∑binom(n,k)²=binom(2n,n). Verwende die Stirling-Abschätzung log(n!)=(n+½)log n−n+½log(2π)+1/(12n)+r_n mit −1/(360n³)≤r_n≤0. Die kombinierte Restgröße ist (r_{2n}−2r_n)/2 und erfüllt die angegebene Schranke. ∎

Der Vorfaktor ist hier n^{−1/4}, **nicht n^{−1/2}**. G19 mit C/√n scheitert sogar am Randfixpunkt in diesem reinen Kompositionsmodell. Das widerlegt G19 für das ursprüngliche Modell nicht; es zeigt aber, dass Pivotasymptotik allein den √n-Faktor nicht liefert und dass die Störung zur echten Geometrie zusätzliche Struktur braucht.

## 6. Modellfamilien und hochpräzise Numerik

Sei l=L/ρ, P(z)=(z−l)(z−l̄)=z²−2Re(l)z+1. Ausgangspunkt ist F₀ des realen Kubikmodells. Wir verwenden:

- `jet_ε`: F_ε=F₀+εP². L, L̄ und μ, μ̄ sind **exakt** fest; Grad 4. Ein reales Kubikpolynom mit diesen vier reellen Bedingungen hat keinen verbleibenden Parameter.
- `mu_s`: μ_s=(1+s)μ, L fest, das passende reale Kubikpolynom.
- `quad_κ`: F=z+κP, ein reales Quadrat mit denselben Randfixpunkten und F′(l)=1+2iκ Im(l).
- `tau_t`: F=F₀, M=tF′∘F. t=0 ist der exakt gelöste reine Fall; t>0 prüft den Effekt des Multiplikationsterms.

Standard ist M=F′∘F. Jede Familie berechnet alle Grade 25…400 in einer O(n²)-Bordering-Rekursion bei 1200+n Bit. Die folgende Extrapolation benutzt log(h_n/h_{n−1})=log λ+c₂/n²+c₃/n³+c₄/n⁴ auf n=250…400; andere Fenster und ein freier 1/n-Term stehen vollständig im JSON. Die Zahl c₂ ist **numerisch**, keine bewiesene Konstante.

| Modell | extrapoliertes λ | c₂ | Sitzung in Sekunden |
|---|---:|---:|---:|
| base | 1.1737972269 | -0.427268 | 43.0 |
| jet_0.01 | 1.1786504203 | -0.427266 | 45.6 |
| jet_-0.01 | 1.1689495448 | -0.427270 | 44.9 |
| mu_0.01 | 1.1741204529 | -0.422181 | 49.1 |
| mu_-0.01 | 1.1734724195 | -0.432158 | 49.9 |
| quad_0.2 | 1.0860217742 | -1.150486 | 39.7 |
| quad_0.31382 | 1.1359010295 | -0.609173 | 41.0 |
| tau_0.1 | 1.1737972269 | -1.621696 | 50.3 |
| tau_10 | 1.1737972270 | -1.621761 | 48.6 |

Die Differenz zwischen jet_−0,01 und jet_+0,01 beträgt ungefähr 0,00970. Bei identischen Eckdaten ist dies erheblich größer als Fensteränderungen des Fits. Für positive t=0,1 und t=10 bleibt der extrapolierte λ-Wert nahe dem Ausgangswert, während c₂ stark wechselt. Ein universeller Übergang bei t=0 ist ausgeschlossen: dort ist λ=F₀′(0)≈0,87390617 exakt. **Konjektur:** im geometrisch stabilen Bereich hängt die N²-Hauptordnung von F, der logarithmische nächste Term aber zusätzlich von M ab; keine Uniformität bis t=0.

## 7. Direkt geprüfte Kandidaten

**Numerisch verworfen, nur in der unten genau bezeichneten Form:**

- exp(mean_θ log max(1,|F(e^{iθ})|)): kubisch 1,14770091 statt etwa 1,17379723.
- exp(2 mean_θ log max(1,|F(e^{iθ})|)): kubisch 1,31721738.
- cap(W̄)/ρ: kubisch 1,10392992.
- (cap(W̄)/ρ)²: kubisch 1,21866128.

Die Kreisquadratur benutzt 131072 Mittelpunkte in gewöhnlicher Arithmetik, nur für diese groben Ausschlüsse. Die W̄-Kapazität wurde unabhängig aus der expliziten Außenabbildung der beiden Kreise berechnet: mit c=F(0), R=|l−c|, α₋=2arg(c+R−l), α₊=arg l und p=π/(α₊−α₋) ist cap(W̄)/ρ=p Im(l)/sin(−pα₋). Die Möbiusabbildung (z−l)/(z−l̄) schickt den Außenraum auf den Sektor (α₋,α₊); seine Potenzabbildung und der Halbebenen-Green-Kern liefern die Formel. Die genannten Zahlen sind keine Intervallzertifikate.

Das verwirft **keine** passend gewichtete Symbolformel auf ℳ und keinen korrekt formulierten Gleichgewichtswert. Ω_y und die geklebte Fläche sind bisher selbst nur durch ein numerisches freies Randproblem bestimmt. Ihre Verwendung als bekannte, bewiesene Eingabe würde die offene Mathematik verschieben.

## 8. Ein exakter Determinantenansatz statt einer fertigen Kapazitätsformel

**Lemma G30 [bewiesen].** Mit m=n+1, p_j(z)=F(z)^j+M(z)z^j, Δ(s)=∏_{i<j}(s_j−s_i) gilt auf einem kleinen Kreis Γ um 0
\[
D_n=\frac1{m!}\oint_{\Gamma^m}\Delta(z_1^{-1},\ldots,z_m^{-1})
\sum_{\sigma\in\{0,1\}^m}\left(\prod_{\sigma_i=0}M(z_i)\right)
\Delta(s_1,\ldots,s_m)\prod_i\frac{dz_i}{2\pi iz_i},\quad
s_i=\begin{cases}z_i&\sigma_i=0,\\F(z_i)&\sigma_i=1.\end{cases}
\]

*Beweis:* Entwickle beide Determinanten im Integral von det(z_i^{−k})det(p_j(z_i)). Fubini und Koeffizientenextraktion ergeben m!det J_n. Entwickle die zweite Determinante zeilenweise in die beiden Zweige. ∎

Das ist eine exakte Zwei-Zweig-Integralform. **Kein positives logarithmisches Gas:** Die Integrale und Gewichte sind komplex, und der führende Beitrag kann sich auslöschen. Ein Maximum des Betrags darf nicht als Determinantenasymptotik ausgegeben werden. Diese Form umgeht Ws Ausschluss der falschen ebenen Orthogonalität, beweist aber noch kein Gleichgewichtsprinzip.

**Bedingter Ansatz G31 [bewiesene Implikation; Hypothesen offen].** Angenommen, auf einer festen Kontur Γ liefern die monischen Paare starke asymptotische Profile
JQ_n(z)=κ_∞^{−n}ζ(z)^n B(z)(1+ε_{Q,n}(z)),
A_n(1/z)=κ_0^nζ(z)^{−n}C(z)(1+ε_{A,n}(z)),
mit gleichmäßig ε→0 und H=(1/2πi)∮BC dz/z≠0. Dann
h_n=(κ_0/κ_∞)^n(H+o(1)), also λ=κ_0/κ_∞.

*Beweis:* Setze die Profile in die exakte Bilinearform W §25.1(c) ein und majorisiere nur die bereits multiplizierten, normierten Fehler. ∎

Wenn κ₀ und κ∞ tatsächlich die Ableitungen einer Uniformisierenden von ℳ an 0 bzw. ∞ sind, wäre λ der skalierungsinvariante Quotient cap(Eᶜ)/rad(Ω_y,0). **Das ist derzeit ein geometrischer Kandidat**, nicht eine für unseren Operator abgeleitete Formel: Zuordnung der beiden Schurprofile, gemeinsame Kontur, nichtverschwindendes H und quantitative Klebungsfehler sind alle offen. G31 enthält starke, quellenverträgliche Profilhypothesen und ist kein Ersatz durch n-te-Wurzel-Asymptotik.

## 9. Literatur und Prüfstatus

- [Böttcher–Widom, *Szegő via Jacobi*, 2006](https://arxiv.org/pdf/math/0604009): **gelesen** wurden die Definitionen, Sätze 1.1–1.2, Jacobi 2.1 und Korollar 2.2 mit Beweis. Die dortige Szegő-Form besitzt ein festes Toeplitz-Symbol und einen kontrollierten Spurklassenschwanz. Unser J hat keine konstante Toeplitzdiagonale und eine N²-Hauptordnung; eine direkte Anwendung wurde nicht behauptet.
- [Cowen–MacCluer, *Composition Operators on Spaces of Analytic Functions*, 1995](https://www.math.purdue.edu/~cowen/Book.html): **nur bibliographisch zitiert**; Autoren-Webseite und Vorwort gelesen, keine Kapitel oder Abschnitte über Determinanten gelesen. Die übliche Selbstabbildungsgeometrie löst unsere über den Abschneidekreis reichende Geometrie nicht automatisch.
- Widom/Szegő–Widom und Böttcher–Silbermanns Buch: **nur zitiert**, kein vollständiger Satz für unsere Klasse gelesen oder benutzt.
- Determinantenasymptotik über Riemann–Hilbert: **Strategiebezeichnung**, hier kein Literaturtheorem als gültig übertragen.

`exact_models.json` enthält exakte rationale Checks der Determinanten/Paare bis n=6 und der Zwei-Zweig-Integralform bis n=2. Diese endlichen Checks sind ergänzende Fehlerkontrolle, nicht die Beweise. +300-Bit-Gegenproben für zwei nichttriviale Familien werden in `family_analysis.json` aufgeführt. Alle Ausgangs-, Fit- und Kontrollfenster bleiben erhalten.

## 10. Verbleibende Lücke und nächste Etappe

Neu bewiesen sind G25–G30 und die **bedingte Implikation** G31. λ und die starke Profilhypothese für das ursprüngliche Modell sind offen. Nächster Schritt ist die getrennte Untersuchung der nächsten Ordnung und der **komplexen**, nicht bloß betragsmäßigen 𝒟-Profile. Ein guter Pivot-Fit allein reicht dafür nicht.
