# Runde 2: (H₀) als exaktes Cauchy-Problem auf festen Konturen

Bearbeitung vom 11.10.2026. Auf Janis' Rückmeldung hin behandelt diese Runde allein die ursprünglichen quadratischen Abschnitte J_N und (H₀). Keine neue Grad-5N-Konstruktion, keine K2-Arbeit. Unabhängige Prüfung durch Strang P steht aus.

**Ergebnis:** Die feste Cauchy-Operator-Fassung ist vollständig ausgeschrieben, einschließlich des endlichen Projektors und der gemeinsamen Korrektur I_f − B_y. Im kubischen Modell hat diese Korrektur auf einer großen festen Kontur exakt Rang höchstens 6. Dazu kommen eine exakte Viererrekursion der transportierten Testkerne und neue SVD-freie Rechnungen. **Eine polynomiale Schranke für alle N und damit (H₀) ist noch nicht bewiesen.**

Gelesen: W-BERICHT §15–§27, W/OFFENES_PROBLEM.md, O-BERICHT §32–§34, X/KONSTANTEN.md. Die dort ausgeschlossenen lokalen, rein ℓ²-variationellen und rein logarithmischen Wege werden nicht als Beweis verwendet.

## 1. Konventionen

Hier wird die physische Variable x benutzt, damit die Zahlen unmittelbar mit O38 vergleichbar sind. Es sei ρ = |L| und

\[
\widetilde J_N=D J_ND^{-1},\quad D=\operatorname{diag}(\rho^k),\qquad
\widetilde J_N^T y=u,
\]
\[
Y(x)=\sum_{k=0}^N y_k(\rho/x)^k,\qquad
\mathcal U(x)=\sum_{k=0}^N u_k\rho^k x^{-k-1},\qquad m=f'\circ f.
\]

Für Kneser gilt die letzte Formel als Keim; m = exp(x)/f_K′(x) wird nur für die exakte Halbiteration verwendet. In normierter Variable t = x/ρ ist die rationale Funktion T um den Faktor ρ größer als das hier benutzte physische T. G ist dimensionslos; G^T hat die physische Normierung von O38.

Die Transposition in den Keimidentitäten ist bilinear. In Normabschätzungen wird das hermitesche ℓ² benutzt. Alle Matrizen der hier gerechneten beiden Modelle sind reell.

## 2. Eine gemeinsame Cauchy-Operation statt getrennter Restabschätzungen

**Satz G3 — bewiesen.** Sei Ω ein beschränktes Jordan-Gebiet mit rektifizierbarem Rand Γ, 0 ∈ Ω. f und m seien holomorph in einer Umgebung von Ω̄. Setze a = f(0) ≠ 0 und

\[
\mathcal C_\Gamma Y(t)=\frac1{2\pi i}\int_\Gamma
 \frac{Y(z)}z\left\{\frac1{t-f(z)}+\frac{m(z)}{t-z}\right\}\,dz.
\tag{G3}
\]

Es seien t ∉ Γ ∪ f(Γ), t ≠ 0,a. Mit den rationalen Funktionen aus Runde 1,

\[
R_y(t)=\sum_{k=0}^N y_k\rho^k[z^k]\frac1{t-f(z)},\quad
M_y(t)=\sum_{j=0}^N\left(\sum_{k=j}^N y_k\rho^k m_{k-j}\right)t^{-j-1},
\]

gilt:

1. Für t im Außengebiet von Ω und f(Ω̄) ist \(\mathcal C_\Gamma Y=R_y+M_y\).
2. Für t ∈ Ω ist \(\mathcal C_\Gamma Y=I_f-B_y\), wobei

\[
I_f(t)=\frac1{2\pi i}\int_\Gamma\frac{Y(z)}{z(t-f(z))}\,dz,
\quad B_y(t)=\frac{m(t)Y(t)}t-M_y(t).
\tag{G4}
\]

3. Für t ∈ Ω und einfache Urbilder in Ω gilt exakt

\[
T_y(t)=\frac{m(t)Y(t)}t+
 \sum_{\substack{f(z)=t\\z\in\Omega}}\frac{Y(z)}{zf'(z)}
 -\mathcal U(t)+\mathcal C_\Gamma Y(t),
\quad T_y=R_y+M_y-\mathcal U.
\tag{G5}
\]

Bei mehrfachen Urbildern ersetzt man die Summe durch die entsprechenden Residuen.

**Beweis.** Im Außengebiet liefern die Residuen in 0 für die beiden Integranden R_y und M_y. Im Inneren hat der zweite Integrand zusätzlich in z = t das Residuum −m(t)Y(t)/t; sein Integral ist also −B_y. Für den ersten Integranden ist das Residuum an einem einfachen Urbild −Y(z)/(zf′(z)); daher I_f = R_y − Summe. Einsetzen gibt (G5). ∎

Dies ist die korrigierte Gleichung (10) der ersten Runde. I_f und B_y werden in (G3) gemeinsam transportiert. Die Operation besitzt je nach Lage des Auswertungspunkts verschiedene analytische Zweige; ihre Innen- und Außenwerte dürfen nicht gleichgesetzt werden.

## 3. Der endliche Projektor gehört zur Gleichung

**Satz G4 — bewiesen.** Wähle eine feste zweite Kontur Ξ um Ω̄ ∪ f(Ω̄), auf der die Außenform von (G3) gilt. Für eine dort durch eine Laurent-Reihe bei ∞ gegebene Funktion h definiere

\[
(\Pi_N^-h)(t)=\frac1{2\pi i}\int_\Xi h(s)
 \sum_{j=0}^N\frac{s^j}{t^{j+1}}\,ds.
\tag{G6}
\]

Dann ist das ursprüngliche System mit N+1 Unbekannten genau äquivalent zu

\[
\Pi_N^-\mathcal C_\Gamma Y=\mathcal U,
\qquad Y\in\mathcal L_N:=\operatorname{span}\{1,(\rho/z),\ldots,(\rho/z)^N\}.
\tag{G7}
\]

**Beweis.** Der j-te Laurent-Koeffizient des Außenwerts von \(\mathcal C_\Gamma Y\) ist
\(\sum_k y_k\rho^k[z^k]\{f(z)^j+m(z)z^j\}\).
Dies ist \(\rho^j(\widetilde J_N^T y)_j\). Der Projektor entnimmt genau j = 0,…,N. ∎

Die Konturen bleiben fest, während N wächst. Damit ist die Cauchy-Fassung konkret; sie ersetzt das quadratische Abschnittsproblem nicht durch ein System mit zusätzlichen Freiheitsgraden.

### Gewichtete Fassung

**Bewiesen als exakte Äquivalenz.** Für feste positive, stetige Gewichte w_N auf Γ und v_N auf Ξ, beispielsweise w_N = exp(−N H_Γ), setze \(\widehat Y=w_NY\). Die Gleichung wird

\[
\Pi_N^{-,v}\left[v_N\mathcal C_\Gamma(w_N^{-1}\widehat Y)\right]=v_N\mathcal U,
\quad \Pi_N^{-,v}=v_N\Pi_N^-v_N^{-1},
\quad\widehat Y\in w_N\mathcal L_N.
\tag{G8}
\]

Der gewichtete Kern enthält den Faktor v_N(t)/w_N(z). Sowohl der konjugierte Projektor als auch der gewichtete endliche Funktionenraum sind wesentlich. Das bloße Ersetzen von J durch einen unendlichen gewichteten Kompositionsoperator ist nicht äquivalent zu (G8).

Die Datennorm bleibt ausdrücklich ‖u‖₂. Äquivalent ist
\(\|u\|_2^2=\rho^2(2\pi)^{-1}\int_0^{2\pi}|\mathcal U(\rho e^{i\theta})|^2d\theta\).
Die gewöhnliche L²-Norm von U auf der großen Außenkontur Ξ ist hierzu nicht N-unabhängig äquivalent und wird nicht an ihre Stelle gesetzt. Eine neue Radialgewichtung der ursprünglichen Koeffizientennorm wird damit nicht als Ausweg eingeführt.

Damit ist **noch keine Inversennorm** gewonnen. (G8) ist die Ausgangsgleichung für eine globale gewichtete Infsup-Abschätzung. Sie hält Quelle, Phasen, Abschneidung und die Korrektur gleichzeitig fest.

## 4. Kubisch: die Korrektur auf einer großen Kontur ist endlichdimensional

**Satz G5 — bewiesen.** Sei f ein Polynom vom Grad d ≥ 1 und m ein Polynom vom Grad s. Wähle einen festen Kreis Γ_R, der für alle betrachteten t sämtliche Urbilder von t und t selbst enthält. Dann

\[
I_f(t)=0,\qquad \mathcal C_{\Gamma_R}Y(t)=-B_y(t),
\]
\[
B_y(t)=\sum_{j=0}^{s-1}\left(\sum_{k=0}^{\min(N,s-j-1)}y_k\rho^k m_{j+k+1}\right)t^j.
\tag{G9}
\]

Die Korrekturoperation hat Rang höchstens s. Für das kubische Modell ist s = 6: sie hängt nur von y₀,…,y₅ ab, für jedes N ≥ 5.

**Beweis.** Der Integrand von I_f fällt bei ∞ mindestens wie z^(−d−1), hat also dort Residuum 0. Sein Integral um sämtliche endlichen Pole ist 0. Der holomorphe Teil von mY/t enthält genau die in (G9) aufgeführten Potenzen. ∎

Für das reelle Modell

\[
f(x)=x+(x^2-2\Re L\,x+|L|^2)(\alpha+\beta x),\quad f'(L)=\sqrt L,
\]

genügt Γ_R = {|z| = 16ρ} für alle |t| ≤ 1.5ρ. **Rigorose geometrische Gegenprüfung**, 512 Bit und L = −W₋₁(−1):

\[
|f_3|(16\rho)^3-\sum_{j=0}^2|f_j|(16\rho)^j-1.5\rho
\in[95.1508465794872708000714188333\pm3.31\cdot10^{-29}]>0.
\]

Rouché zählt damit alle drei Urbilder innerhalb dieses Kreises. Die Aussage I_f = 0 ist unabhängig von einer Inversenschranke für J_N.

### Explizite Majorante für alle Quellen

**Bewiesen, abhängig von den sechs Antwortnormen.** Schreibe \(\ell_{N,k}=\widetilde J_N^{-1}e_k\). Aus y_k = u^Tℓ_{N,k} folgt für ‖u‖₂ = 1:

\[
|B_y(t)|\le\sum_{k=0}^5\|\ell_{N,k}\|_2\rho^k
 \sum_{j=0}^{5-k}|m_{j+k+1}|\,|t|^j.
\tag{G10}
\]

Eine einheitliche Schranke für diese sechs Spalten würde also die gemeinsame Korrektur auf Γ_R mit einer expliziten, N-unabhängigen Konstante kontrollieren. **Diese Spaltenschranke für alle N ist hier nicht bewiesen.**

## 5. Ein exakter Transport mit vier Zuständen

**Satz G6 — bewiesen.** Für f(z) = a+bz+cz²+dz³ definiere die dualen Testkerne

\[
\psi_k(t)=[z^k]\frac1{t-f(z)}+m(t)t^{-k-1}.
\]

Ihre Erzeugende ist

\[
\sum_{k\ge0}\psi_k(t)z^k=\frac1{t-f(z)}+\frac{m(t)}{t-z}.
\tag{G11}
\]

Für jedes k ≥ 0 gilt die exakte Rekursion

\[
t(t-a)\psi_{k+4}-((t-a)+bt)\psi_{k+3}
 +(b-ct)\psi_{k+2}+(c-dt)\psi_{k+1}+d\psi_k=0.
\tag{G12}
\]

**Beweis.** Multiplikation von (G11) mit
\((t-f(z))(t-z)=t(t-a)-((t-a)+bt)z+(b-ct)z²+(c-dt)z³+dz⁴\)
gibt einen Zähler vom Grad höchstens 3 in z. Alle Koeffizienten ab Grad 4 verschwinden. ∎

Schreibt man \(q_k=[z^k](t-f(z))^{-1}\), ist eine äquivalente Vierzustandsform

\[
\begin{pmatrix}q_{k+1}\\q_k\\q_{k-1}\\t^{-k-2}\end{pmatrix}
=\begin{pmatrix}
b/(t-a)&c/(t-a)&d/(t-a)&0\\
1&0&0&0\\0&1&0&0\\0&0&0&1/t
\end{pmatrix}
\begin{pmatrix}q_k\\q_{k-1}\\q_{k-2}\\t^{-k-1}\end{pmatrix}.
\]

Die vier charakteristischen Zweige sind 1/t und die drei Kehrwerte der Urbilder von t. Damit lässt sich der globale Transport der Testkerne ohne SVD und ohne angenommene Grenzfläche verfolgen. W §25.4 hatte eine Struktur der Größe d+1 bereits vorgeschlagen; hier ist die zugehörige Rekursion ausgeschrieben. **Eine klassische 4×4-RH-Aufgabe oder deren starke Asymptotik ist damit noch nicht konstruiert.** Insbesondere wird keine Gonchar–Rakhmanov-Orthogonalität behauptet.

Beachte die beiden verschiedenen Kerne: (G3) enthält m(z); die duale Momentenerzeugende (G11) enthält m(t). Das ist kein austauschbarer Faktor.

## 6. Neue Rechnungen am kubischen Modell

**Status: hochpräzise numerisch, keine Intervallzertifikate der inversen Matrizen.** `experimente/h0_cauchy.py` löst die Systeme mit 720 bis 2000 Bit und `algorithm="approx"`. Es benutzt keine SVD. Pro N werden 194 reelle rechte Seiten gelöst: Kreisproben für G/G^T, 20 feste γ*-Proben für G und die ersten sechs Koordinatenquellen.

Die Kreisproben verwenden genau O38s Winkelraster von 10° und dessen gerundeten Exponenten β = 1/1.1182. Die γ*-Proben liegen auf dem durch das kubische Modell definierten h(A-Bogen) und der W-Referenzlinse mit x₀ = 0.7. Genau dieselben festen Punkte werden anschließend für Kneser verwendet. Es sind **Stichpunkte, keine zertifizierte Überdeckung der Kurve**. Die Bucht wird durch numerische Fortsetzung des zentralen Urbilds konstruiert; die Linse durch das Poisson-Integral des vorhandenen Referenzmodells. Die Eingabe dieser Geometrie ist allein das kubische Modell, kein Kneser-Halbiterationswert.

| N | gemischte KB / √N | max G auf ganzem Kreis / √N | max G^T auf ganzem Kreis / √N | max G an γ*-Stichpunkten / √N |
|---:|---:|---:|---:|---:|
| 200 | 0.842994 | 0.842994 | 1.972115 | 0.536181 |
| 300 | 0.820292 | 0.820292 | 2.123209 | 0.511939 |
| 400 | 0.805286 | 0.805286 | 2.239391 | 0.493950 |
| 600 | 0.785414 | 0.785414 | 2.417277 | 0.480025 |
| 800 | 0.772095 | 0.772095 | 2.554679 | 0.471537 |
| 1000 | 0.762158 | 0.762158 | 2.668383 | 0.468121 |

Die relative Abweichung der neuen gemischten KB-Werte von O38 ist höchstens 1.7·10⁻¹¹; dies ist ein numerischer Konsistenzcheck derselben Modellparameter und desselben Rasters. Der wichtige neue Prüfgegenstand ist die **volle Majorante G an festen γ*-Punkten**, nicht nur einzelne Eckmoden dort. Sie deckt sämtliche Einheitsquellen ab.

**Präzisierung:** O38s KB mischt innen G und außen G^T. Eine Schranke dieser gemischten Werte allein ist keine Randwertschranke für beide Funktionen auf dem ganzen Kreis. Hier wurde der ganze G-Kreis zusätzlich geprüft. Der ganze G^T-Kreis liegt deutlich darüber; für ihn wird 0.84√N nicht behauptet. Die Zahlen prüfen β ≈ 0.8943 auf dem vorgegebenen Raster, beweisen aber weder diese asymptotische Skala noch den Vorfaktor für alle N.

### Korrektur und niedrige Spalten

Die sechs gemessenen Spaltennormen bei N = 1000 sind

\[
(0.5232435530,\ 0.6917165418,\ 0.8141014016,\
 0.8733889548,\ 0.8999132513,\ 0.9140741476).
\]

Für N = 200,…,1000 stimmen sie auf rund 14 Dezimalstellen überein. Einsetzen in (G10) liefert bei |t| ≤ 1.5ρ in allen sechs Läufen die **numerische** Majorante 1.016000964. Das ist starke Evidenz für eine harmlose Korrektur auf der großen Kontur, kein All-N-Beweis der Spaltennormen.

### Rate c

Das vorhandene, vom kubischen Modell her aufgebaute Referenzpotential liefert H(−1) = **0.01176285525**, passend zu W's c ≈ 0.0118. Das ist die Auswertung des vorhandenen Referenzmodells, keine neue rigorose Lösung des freien Randwertproblems.

Ein freier Dreiparameterfit ln σ₁ = cN+b ln N+a an den vorhandenen O38-Werten N = 200,…,1000 liefert c = 0.0127976 und b = −1.3473. Er zeigt die Vorasymptotik und identifiziert c nicht verlässlich mit einem All-N-Grenzwert. Dieser Fit wird ausdrücklich **nicht** als neue Herleitung von c = 0.0118 ausgegeben. Die Rate bleibt hier ein Test für das äußere Gewicht; die starke Schranke wird daraus nicht abgeleitet.

## 7. Gemeinsames Integral: echte Auslöschung, positive Majoranten verlieren sie

**Status: numerisch.** `experimente/h0_kernel_probe.py` benutzt die gespeicherten approximativen W-Eckvektoren y = σ₁Vt[0]. Ihre Float-Koeffizienten werden als dyadische Eingangsdaten behandelt; die Auswertung erfolgt vollständig mit **2400-Bit-acb** auf dem festen Kreis |z| = ρ. Der Zielpunkt ist der reelle Punkt h(a+|L−a|), etwa 0.827ρ.

Für die kleine Kontur enthält Ω nur das zentrale Urbild. Daher ist (G3) exakt die Summe der zwei fernen Urbildterme minus B_y. Die trapezförmige Rechnung mit 2048 und 4096 Punkten wurde mit dieser Residuenform verglichen.

| N | Betrag des gemeinsamen E = I_f−B_y | Integral der Beträge | gewichtete Cauchy–Schwarz-Majorante |
|---:|---:|---:|---:|
| 200 | 5.64·10⁻⁸ | 1.415 | 2.377 |
| 600 | 1.14·10⁻⁸ | 46.19 | 192.22 |
| 1000 | 2.94·10⁻⁹ | 3358.89 | 19365.07 |

Für N = 1000 ist das Verhältnis „Integral der Beträge / Betrag des Integrals“ etwa **1.14·10¹²**. Die Residuenabweichung ist in allen drei Rechnungen < 10⁻²⁷¹, die Differenz nach Verdoppelung des Gitters < 10⁻¹⁰². Diese Zahlen betreffen die gewählten gespeicherten Polynome; die Genauigkeit ihrer Eckmoden-Eigenschaft wird dadurch nicht auf 2400 Bit verbessert. Eine rigorose Quadraturfehlerschranke wurde nicht hergeleitet; daraus folgt kein H₀-Zertifikat.

Das getestete Gewicht ist w_N = exp(−N Λ_ref) auf Γ. Der gemessene gewichtete L²-Spurnormwert von Y bleibt etwa 0.35. Die duale gewichtete Kernnorm wächst aber von 6.85 auf 54359.4. **Befund:** Ein bloßes positives Gewicht plus Schur/Cauchy–Schwarz liefert keine polynomiale Kontrolle. Der nächste Schritt muss die komplexe Auslöschung und den endlichen Projektor behandeln; die Normierung der Quelle darf nicht entfernt werden.

**Bewiesene Grenze dieses Abschätzungsschritts.** Für jedes positive Gewicht w auf dem festen Kreis gilt mit normiertem Bogenmaß
\[
\|wY\|_2\,\|K/w\|_2\ \ge\ \int_\Gamma |YK|\,\frac{|dz|}{2\pi\rho}.
\]
Das Infimum über positive Gewichte ist genau die rechte Seite: Cauchy–Schwarz beweist die Ungleichung; das regularisierte Gewicht
\(w_\varepsilon^2=(|K|+\varepsilon)/(|Y|+\varepsilon)\) erreicht sie im Grenzwert. Somit verliert auch ein für dieses Y optimiertes positives Gewicht die Auslöschung des Integrals. Die gemessenen großen Verluste sind nicht allein eine schlechte Wahl des getesteten Gewichts. Ein exponentieller All-N-Grenzwert dieses Betragsintegrals wird aus den drei Zahlen nicht behauptet.

Technischer Fehlversuch: 320 Bit gaben bei langen komplexen Horner-Auswertungen zu breite Intervalle (Wrapping). Die erfolgreiche Rechnung wurde deshalb mit 2400 Bit vollständig neu ausgeführt. Die Zahlen oben stammen nur aus diesem Lauf.

## 8. Das stärkere, SVD-freie Ziel als präzise Formungleichung

**Satz G7 — bewiesen als Äquivalenz.** Für invertierbares reelles \(\widetilde J_N\) und \(e_x=((\rho/x)^k)_{k=0}^N\) gilt

\[
G_N(x)=\|\widetilde J_N^{-1}e_x\|_2
=\sup_{\|u\|_2=1}|Y_u(x)|.
\]

Das Supremum läuft über komplexe Einheitsquellen; G majorisiert insbesondere sämtliche reellen Eckmoden. Die Schranke \(G_N(x)\le C N^\kappa\) ist äquivalent zur hermiteschen Formungleichung

\[
\overline{e_x}e_x^T\ \preceq\ C^2N^{2\kappa}
 \widetilde J_N\widetilde J_N^T.
\tag{G13}
\]

**Beweis.** Y_u(x) = u^T J̃_N⁻¹e_x; das Supremum ist die euklidische Norm. Für y = J̃_N⁻ᵀu ist (G13) genau |y^Te_x|² ≤ C²N^(2κ)‖J̃_N^Ty‖². ∎

Damit wäre eine globale Schranke G_N|γ* ≤ C√N ein **stärkerer Ersatz für (H₀)**: sie kontrolliert alle Einheitsquellen, ohne Eckmoden wählen oder zählen zu müssen. Die neue Stichpunktrechnung unterstützt dieses stärkere Ziel. **(G13) ist eine überprüfbare Zielungleichung, kein hier bewiesener All-N-Baustein.**

Eine globale gewichtete Randfassung benutzt den festen Gramoperator

\[
W_{\gamma,N}=\int_{\gamma^*}\overline{e_x}e_x^T\,d\nu(x),
\quad \int_{\gamma^*}|Y(x)|^2d\nu=y^*W_{\gamma,N}y.
\]

Gesucht ist eine polynomiale Infsup-Abschätzung von (G8) in diesem Randnormraum, mit anschließender expliziter Auswertungsnorm. Dies ist eine Aussage über den gesamten Cauchy-Transport, nicht eine Rayleigh-Extremalität in der ursprünglichen ℓ²-Koeffizientennorm. Die rationale Momentprojektion bleibt dabei erhalten.

## 9. Übertragung auf Kneser

### Algebraische Übertragung

**Bewiesen:** (G3), (G7) und (G8) gelten auf jeder festen zulässigen Keimkontur ebenfalls für nichtpolynomiales f_K und m_K. Das System benötigt nur f_K-Koeffizienten bis N+1. B_y bleibt der holomorphe Teil mit

\[
[t^j]B_y=\sum_{k=0}^N y_k\rho^k m_{j+k+1}.
\]

**Nicht übertragen:** Der große kubische Kreis Γ_R und Rang 6 gelten nicht für Kneser. Ein Kreis, der die Singularitäten L und L̄ ins Innere nimmt, ist kein zulässiges holomorphes Gebiet. Ebenso gibt es für nichtpolynomiales f_K keine feste Viererrekursion wie (G12).

Wenn f_K und m_K auf D holomorph und stetig auf D̄ sind, darf man die Cauchy-Operation auf |z| = ρ als Grenzwert innerer Kreise verwenden, für Ziele mit positivem Abstand zu ∂D und f_K(∂D). Für Ziele an L benötigt man Spurgrenzwerte und eine Kontrolle des singulären Grenzoperators. Diese Kontrolle ist ein **offener Teil des starken Beweises**, nicht durch formales Fortsetzen der kubischen großen Kontur erledigt. Außerhalb D bleibt die rationale T-Fassung aus den Jets verfügbar, ohne globale inverse Äste von Kneser vorauszusetzen.

### Numerischer Vergleich

**Numerisch, zertifizierte Mittelpunkte als Eingabe, keine Intervallmatrix-Zertifizierung.** Die eigenen Kneser-Läufe lesen a₀,…,a₆₄₀ und stoppen bei N = 600. Ein Lauf N > 639 wird vom Skript ausdrücklich zurückgewiesen.

| N | kubisch: KB / √N | Kneser: KB / √N | Kneser: G an denselben festen γ*-Proben / √N |
|---:|---:|---:|---:|
| 200 | 0.842994 | 0.843190 | 0.536299 |
| 300 | 0.820292 | 0.820468 | 0.511953 |
| 400 | 0.805286 | 0.805450 | 0.493953 |
| 600 | 0.785414 | 0.785563 | 0.479998 |

Das ist die erwartete numerische Nähe auf deutlich besser als 0.1 %. Auch die Majorante an denselben festen γ*-Stichpunkten stimmt eng überein. Die Kurve wurde dafür bewusst vom kubischen Modell übernommen; ihre Zulässigkeit für Knesers Geometrie und eine vollständige kontinuierliche Kurvenüberdeckung sind nicht zertifiziert.

### Was eine quantitative Übertragung verlangen würde

**Bewiesen als Differenzidentität:** Auf gemeinsamen zulässigen Konturen ist

\[
K_K(z,t)-K_c(z,t)
=\frac{f_K(z)-f_c(z)}{(t-f_K(z))(t-f_c(z))}
 +\frac{m_K(z)-m_c(z)}{t-z}.
\tag{G14}
\]

Für einen gemeinsamen gewichteten endlichen Raum sei A_{c,N} die Operation in (G8) und Δ_N ihre Differenz zur Kneser-Operation. Wenn
\(\|A_{c,N}^{-1}\|\le M N^a\) und \(\|A_{c,N}^{-1}\Delta_N\|\le\eta<1\), folgt
\(\|A_{K,N}^{-1}\|\le M N^a/(1-\eta)\) durch die Neumann-Reihe. **Die zweite Hypothese ist nicht bewiesen.** Nahe Übereinstimmung von f und f′ in L oder der Messwerte reicht dafür nicht; die Konturen, Gewichte und Projektoren müssen gemeinsam kontrolliert werden.

## 10. Konstanten für die Gesamtbilanz

Die für X relevante Ausgabe bleibt eine All-N-Schranke mit expliziten C_E, κ_E und Green-Reserve ε_E; Stichpunkttabellen ersetzen sie nicht.

Auf W1''s Weg gilt, **bedingt auf die dortigen Geometrievoraussetzungen und deren Konstanten**,

\[
\sup_{\partial K_T}|T_u|\le C(f)\sup_{\gamma^*}G_N
 +\frac{\sqrt{N+1}}{\min_{\partial K_T}|\zeta|}.
\tag{G15}
\]

W §23.6 gibt für das kubische Modell die numerischen Werte C(f) = 4.95 und 1/min|ζ| = 0.728. Würde die Stichpunkthülle G_N ≤ 0.54√N als All-N-Satz auf der ganzen Kurve bewiesen und diese Geometrie zertifiziert, gäbe (G15) für N ≥ 200 z. B. **C_E = 3.5, κ_E = 1/2**. Diese Zahlen sind eine **bedingte Zielbilanz**, keine bewiesenen Eckkonstanten. Ein neuer N₀ wird deshalb nicht als zertifiziert ausgegeben.

## 11. Was neu ist, welche Lücke bleibt, nächster Schritt

1. **Neu ausgeschrieben und bewiesen:** die gemeinsame Cauchy-Operation (G3), die exakte endliche und gewichtete Fassung (G7)/(G8), die Rang-6-Korrektur samt expliziter Majorante (G9)/(G10), die Vierzustandsrekursion (G12), sowie die SVD-freie Formulierung (G13). Der große kubische Kreis ist unabhängig bei 512 Bit geprüft. Keine Priorität gegenüber anderer Literatur wird beansprucht.
2. **Präzise Lücke:** Die globale gewichtete Infsup-Abschätzung mit dem endlichen Projektor ist offen. Auch die sechs niedrigen Spaltennormen sind nur für sechs Grade gemessen. Es fehlt eine polynomiale All-N-Kontrolle des komplexen Cauchy-Transports; ein Betrags-Schurargument schließt sie ausdrücklich nicht. Für Kneser kommen die Spuren an L/L̄ und die gewichtete Operatorstörung hinzu.
3. **Ein nächster Schritt:** Am kubischen Modell eine globale Parametrix für den Vierzustands-Transport und den Momentprojektor (G8) auf den festen Kurven konstruieren, mit der Quelle im System, und daraus (G13) auf γ* mit einem expliziten polynomialen Rest beweisen. Die Rang-6-Korrektur wird dabei als endliche Störung separat geführt; die äußere Phase wird durch W's Referenzgeometrie geprüft. Erst eine solche Operatorabschätzung erlaubt die Übertragung mittels (G14), nicht die lokalen Vergleichsdaten allein.

## Dateien und Reproduktion

- `experimente/h0_cauchy.py`: SVD-freie Antwortnormen, feste γ*-Stichpunkte, niedrige Spalten.
- `experimente/h0_kernel_probe.py`: gemeinsamer Integralkern, gewichtete Spurnorm, Auslöschung.
- `experimente/h0_exact_checks.py`: exakte Symbolprüfungen und rigorose Rouché-Marge.
- `experimente/h0_plot.py`: Prüfbild aus den gespeicherten Daten, ohne Neuberechnung.
- `ergebnisse/h0_cauchy_kubisch.json`, `h0_cauchy_kneser_z.json`, `h0_kernel_probe.json`, `h0_exact_checks.json`.

Aus dem Forschungsordner `D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\g_gpt_hilfssaetze` in PowerShell mit dem vorhandenen Interpreter:

```powershell
& 'C:\Users\janis\miniconda3\python.exe' -B .\experimente\h0_cauchy.py --model kubisch --degrees 200,300,400,600,800,1000
& 'C:\Users\janis\miniconda3\python.exe' -B .\experimente\h0_cauchy.py --model kneser_z --degrees 200,300,400,600
& 'C:\Users\janis\miniconda3\python.exe' -B .\experimente\h0_kernel_probe.py
& 'C:\Users\janis\miniconda3\python.exe' -B .\experimente\h0_exact_checks.py
```

Alle neuen Skripte schreiben nur in den eigenen Ergebnisordner. `--project-root` und `--output-dir` erlauben ausdrücklich andere Ein- und Ausgabeorte. Keine fremden Skripte wurden geändert. Es wurden keine externen Sätze als Beweisbausteine verwendet; die ausgeschriebenen Identitäten folgen aus Koeffizientenentnahme, Residuensatz und elementarer Matrixrechnung.
