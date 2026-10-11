# Polynom-Halbiteration: erste Ergebnisse

Bearbeitung vom 11. Oktober 2026. Ablage: `g_gpt_hilfssaetze/`. **Unabhängige Prüfung durch Strang P steht noch aus.** Die Bezeichnungen „bewiesen“ unten beziehen sich auf die hier ausgeschriebenen Herleitungen; sie bedeuten keine bereits erfolgte P-Freigabe. Ausgangspunkt: `PROMPT_GPT.md` und die dort genannten Forschungsberichte. Die ursprüngliche Vermutung mit **Grad N und Übereinstimmung bis Ordnung N** ist weiterhin offen. Ebenso ist (H₀) nicht bewiesen.

Die beiden folgenden Sätze haben vollständige, elementare Beweise. Der erste liefert eine schwächere Konvergenzaussage mit zusätzlichen Polynomkoeffizienten. Der zweite präzisiert die exakte Ausgangsgleichung für Knesers nichtpolynomiale Funktion. „Neu“ bedeutet hier: in dieser Bearbeitung hergeleitet; eine Priorität gegenüber der Literatur wird nicht behauptet.

## 1. Bewiesen: Konvergenz bei Grad höchstens 5N

### Satz 1 — Allgemeine Konstruktion

**Status: bewiesen.** Sei ρ > 0. Sei f holomorph auf D = {|z| < ρ}, mit reellen Taylor-Koeffizienten, und schreibe

\[
F(t)=\frac{f(\rho t)}\rho=\sum_{k\ge0}b_k t^k.
\]

Es gelte:

1. \(\sum_{k\ge0}|b_k|<\infty\), sodass f auf dem abgeschlossenen Kreis durch seine absolut konvergente Taylor-Reihe gegeben ist;
2. \(\alpha=F(0)\in(0,0.37]\);
3. \(f(f(z))=e^z\) in einer Umgebung von 0. Die äußere Funktion ist dort definiert, da f(0) im Inneren von D liegt.

Dann gibt es für jedes N ≥ 1 ein reelles Polynom P_N vom Grad höchstens 5N mit

\[
P_N(P_N(z))=e^z+O(z^{N+1}),
\qquad
P_N\longrightarrow f\quad\text{gleichmäßig auf }\overline D.
\]

Es gilt sogar, mit \(B=\sup_k|b_k|\) und \(q=0.675622435625<1\),

\[
\|P_N-f\|_{\infty,\overline D}
\le \rho\sum_{k>5N}|b_k|
 +\frac{\rho B}{0.47}(N+1)q^N.
\tag{1}
\]

Wenn zusätzlich \(|b_k|\le Ck^{-s}\) für k ≥ 1 und s > 1 gilt, folgt

\[
\|P_N-f\|_{\infty,\overline D}=O(N^{1-s}).
\tag{2}
\]

### Vollständiger Beweis

Setze M = 5N, \(Q_M(t)=\sum_{k=0}^M b_k t^k\) und \(E_M=F-Q_M\). Weil E_M in 0 mindestens Ordnung M+1 hat, ist

\[
G_{M,N}(t)=\frac{E_M(t)}{t^{N+1}}
 =\sum_{k>M}b_k t^{k-N-1}
\]

auf der Einheitskreisscheibe holomorph, auch in 0. Bezeichne mit \(T_{N,\alpha}\) den Taylor-Abschnitt bis Grad N um α. Definiere

\[
H_N(t)=t^{N+1}T_{N,\alpha}G_{M,N}(t),
\qquad
\widehat P_N(t)=Q_M(t)+H_N(t),
\qquad
P_N(z)=\rho\widehat P_N(z/\rho).
\tag{3}
\]

H_N hat Grad höchstens 2N+1. Damit hat P_N Grad höchstens 5N. Alle Koeffizienten sind reell, weil F und α reell sind.

In 0 verschwindet H_N mindestens zur Ordnung N+1. In α verschwindet

\[
F-\widehat P_N=t^{N+1}(G_{M,N}-T_{N,\alpha}G_{M,N})
\]

mindestens zur Ordnung N+1. Folglich stimmen die Jets von \(\widehat P_N\) und F bis Ordnung N sowohl in 0 als auch in α überein. Da \(\widehat P_N(0)=F(0)=\alpha\), ergibt Einsetzen der beiden Jets

\[
\widehat P_N(\widehat P_N(t))
 =F(F(t))+O(t^{N+1})
 =\rho^{-1}e^{\rho t}+O(t^{N+1}).
\]

Nach Rückskalierung gilt die behauptete Kompositionsgleichung.

Für die Normabschätzung wähle r = 0.16 und R = α+r ≤ 0.53. Der Kreis \(|t-\alpha|=r\) liegt in der Einheitskreisscheibe. Auf ihm ist

\[
|G_{M,N}(t)|\le B\sum_{k>M}R^{k-N-1}
 =\frac{B R^{M-N}}{1-R}.
\]

Schreibt man \(T_{N,\alpha}G=\sum_{j=0}^N c_j(t-\alpha)^j\), liefert die Cauchy-Formel

\[
|c_j|\le\frac{B R^{M-N}}{(1-R)r^j}.
\]

Für |t| ≤ 1 gilt daher

\[
|H_N(t)|
\le\frac{B R^{M-N}}{1-R}
 \sum_{j=0}^N\left(\frac{1+\alpha}{r}\right)^j
\le\frac B{0.47}(N+1)
 \left(0.53^4\frac{1.37}{0.16}\right)^N.
\]

Die Klammer ist exakt 0.675622435625. Zusammen mit
\(\|F-Q_M\|_{\infty,|t|\le1}\le\sum_{k>M}|b_k|\) folgt (1), einschließlich des Randes. Die rechte Seite geht wegen der absoluten Summierbarkeit gegen 0. Schließlich gibt
\(\sum_{k>M}k^{-s}\le M^{1-s}/(s-1)\) die Rate (2). ∎

### Folgerung für Kneser

**Status: bewiesen unter den ausdrücklich zugelassenen Strukturannahmen aus §2(a) des Prompts.** Für f = f_K ist α ≈ 0.3627083373365253 < 0.37. Die dort angegebene Koeffizientenasymptotik impliziert \(|b_k|=O(k^{-6.446950720\ldots})\) und absolute Summierbarkeit. Daher liefert Satz 1

\[
\deg P_N\le5N,
\quad P_N\circ P_N=e^z+O(z^{N+1}),
\quad \|P_N-f_K\|_{\infty,\overline D}
 =O(N^{-5.446950720\ldots}).
\]

Der Exponent der Rate ist exakt \(1+\operatorname{Re}\gamma\); die Dezimaldarstellung oben ist abgekürzt. Keine Aussage (H₀), (H), Bulk-Koerzivität oder Korridor-Eindeutigkeit wird dafür benötigt. Die Kneser-Struktur selbst wurde hier nicht unabhängig neu zertifiziert. Die Voraussetzung α < 0.37 wurde zusätzlich mit 512-Bit-Intervallarithmetik aus dem vorhandenen a₀-Zertifikat geprüft; eine zusätzliche Fehlertoleranz 10⁻¹³⁹ deckt dessen Dezimalrundung ab (`ergebnisse/alpha_check.json`).

**Grenze des Ergebnisses:** P_N ist eine neue Folge. Sie ist nicht mit dem ursprünglichen p_N identifiziert. Insbesondere erfüllt ein Polynom vom Grad 5N hier nur die Gleichungen bis Ordnung N; die Gleichungen der Ordnungen N+1 bis 5N werden nicht behauptet.

### Endliche Berechnung der Konstruktion

**Status: bewiesen.** Formel (3) scheint die ganze Funktion F zu benötigen. Tatsächlich genügen die Koeffizienten b_0,…,b_M.

Aus der Kompositionsgleichung folgt F′(0) ≠ 0. Sei χ(v) der inverse Funktionskeim zu F(t)−α in 0. Dann

\[
F(\alpha+v)=\rho^{-1}\exp(\rho\chi(v)).
\tag{4}
\]

Den N-Jet von χ bestimmt man aus b_1,…,b_N durch formale Reiheninversion. Mit (4) erhält man den N-Jet von F in α. Subtraktion des Jets von Q_M und Division durch \((\alpha+v)^{N+1}\) geben genau die Koeffizienten c_0,…,c_N in (3). Alle Operationen sind endlich.

### Numerische Gegenprüfung

**Status: numerisch, keine Intervallzertifizierung.** Mit 160 Dezimalstellen und den Mittelpunkten der vorhandenen 640 zertifizierten Koeffizienten wurden folgende Beispiele berechnet. Der Jetdefekt betrifft die normierte Gleichung \(\widehat P_N\circ\widehat P_N=\rho^{-1}e^{\rho t}+O(t^{N+1})\).

| N | Grad | maximaler Jetdefekt | \(\|P_N-Q_{5N}\|_{\ell^1_\rho}\) | bekannte Taylorschwanzsumme bis 640 |
|---:|---:|---:|---:|---:|
| 5 | 25 | 2.67·10⁻¹⁶¹ | 5.42·10⁻¹⁰ | 2.41·10⁻⁷ |
| 10 | 50 | 1.78·10⁻¹⁶¹ | 4.61·10⁻¹³ | 5.06·10⁻⁹ |
| 20 | 100 | 1.22·10⁻¹⁶⁰ | 5.16·10⁻¹⁷ | 1.12·10⁻¹⁰ |
| 40 | 200 | 2.10·10⁻¹⁵⁸ | 1.70·10⁻²³ | 2.47·10⁻¹² |

Der winzige Jetdefekt prüft die endliche Hermite-Konstruktion. Er bedeutet nicht, dass die zugrunde gelegten Kneser-Koeffizienten 160 Stellen zertifiziert wären. Die letzte Spalte enthält keine Schranke für die unbekannten Koeffizienten oberhalb von 640.

## 2. Bewiesen: rationale adjungierte Identität für holomorphe Keime

### Satz 2 — Keine Polynomannahme an F oder m nötig

**Status: bewiesen.** Seien F und m holomorphe Keime in 0, \(a=F(0)\ne0\), und

\[
A_N[k,j]=[z^k]\{F(z)^j+m(z)z^j\},\qquad 0\le k,j\le N.
\]

Für beliebige komplexe Vektoren y,u definiere

\[
c_\ell=\sum_{k=0}^N y_k[z^k](F(z)-a)^\ell,
\quad
d_j=\sum_{k=j}^N y_k[z^{k-j}]m(z),
\]

\[
R_y(t)=\sum_{\ell=0}^N\frac{c_\ell}{(t-a)^{\ell+1}},
\quad
M_y(t)=\sum_{j=0}^N\frac{d_j}{t^{j+1}},
\quad
U(t)=\sum_{j=0}^N\frac{u_j}{t^{j+1}}.
\]

Dann gilt die exakte Äquivalenz

\[
A_N^{\mathsf T}y=u
\quad\Longleftrightarrow\quad
T_y(t):=R_y(t)+M_y(t)-U(t)=O(t^{-N-2})\quad(t\to\infty).
\tag{5}
\]

Im äquivalenten Fall besitzt T_y die Form

\[
T_y(t)=\frac{q_N(t)}{t^{N+1}(t-a)^{N+1}},
\qquad \deg q_N\le N.
\tag{6}
\]

Die Transposition ist bilinear; hier wird kein komplex konjugiertes Skalarprodukt verwendet.

### Vollständiger Beweis

Für großes |t| gilt als Identität von holomorphen Keimen in z

\[
\frac1{t-F(z)}
 =\sum_{\ell\ge0}\frac{(F(z)-a)^\ell}{(t-a)^{\ell+1}}.
\]

Weil F(z)−a in 0 verschwindet, tragen bei der Entnahme der Koeffizienten z⁰,…,zᴺ nur ℓ ≤ N bei. Also

\[
R_y(t)=\sum_{k=0}^N y_k[z^k]\frac1{t-F(z)}
 =\sum_{j\ge0}\left(\sum_{k=0}^N y_k[z^k]F(z)^j\right)t^{-j-1}.
\]

Der Koeffizient von t⁻ʲ⁻¹ in M_y ist genau der Multiplikationsanteil der j-ten Komponente von A_Nᵀy, für 0 ≤ j ≤ N. Daher verschwinden die ersten N+1 Laurent-Koeffizienten von T_y genau dann, wenn A_Nᵀy=u. Das beweist (5).

Die einzigen möglichen Pole von T_y sind 0 und a, jeweils höchstens der Ordnung N+1. Nach Multiplikation mit dem Nenner in (6) entsteht ein Polynom. Aus dem Verschwinden bei ∞ zur Ordnung N+2 folgt Grad höchstens \(2N+2-(N+2)=N\). ∎

**Prüfung:** Die Identität und die Gradschranke wurden bei N=5 zusätzlich mit exakter rationaler Symbolrechnung geprüft, für Jets nichtpolynomialer analytischer Funktionen. Diese Gegenprüfung ersetzt nicht den obigen Beweis.

### Anwendung auf Kneser und die fehlende analytische Korrektur

**Status: bewiesen.** Setze \(F(z)=f_K(\rho z)/\rho\) und \(m(z)=f_K'(f_K(\rho z))\). Dann ist A_N genau die gewichtete Matrix J_N. Differenzieren der exakten Halbiterationsgleichung liefert in einer Umgebung von 0

\[
m(z)=\frac{e^{\rho z}}{f_K'(\rho z)}.
\tag{7}
\]

Damit benötigt J_N(f_K) nur die Taylor-Koeffizienten von f_K bis Index N+1. Formel (7) gilt hier für die **exakte** Funktion f_K; sie darf nicht ohne Defektkorrektur an beliebigen Fast-Lösungen eingesetzt werden.

Schreibe \(Y(t)=\sum_{k=0}^N y_k t^{-k}\). Für einen holomorphen Keim m ist

\[
\frac{m(t)Y(t)}t=M_y(t)+B_y(t),
\]

wobei B_y in 0 holomorph ist. Seine Taylor-Koeffizienten sind explizit

\[
[t^j]B_y(t)=\sum_{k=0}^N y_k m_{j+k+1},\qquad j\ge0.
\tag{8}
\]

Für ein polynomiales m ist B_y ein Polynom. Für ein allgemeines holomorphes m ist diese Schlussfolgerung nicht erlaubt. Die stets korrekte Identität lautet

\[
T_y(t)=R_y(t)+\frac{m(t)Y(t)}t-B_y(t)-U(t),
\tag{9}
\]

als lokale meromorphe Identität, mit rationaler Fortsetzung von T_y über (5).

**Zusätzliche bewiesene Folgerung:** Ist m kein rationaler Keim und y ≠ 0, kann B_y kein Polynom sein. Andernfalls wäre \(m(t)=t(M_y(t)+B_y(t))/Y(t)\) ein rationaler Keim, ein Widerspruch.

**Bewiesenes Gegenbeispiel zur unqualifizierten Verallgemeinerung:** Für N=0, y₀=1, F(z)=a+z, m(z)=eᶻ und u₀=2 ist

\[
T_y(t)=\frac1{t-a}-\frac1t=\frac a{t(t-a)},
\qquad B_y(t)=\frac{e^t-1}{t}.
\]

B_y ist kein Polynom. Kein Polynom P kann die nichtkorrigierte Funktion \(e^t/t+1/(t-a)-2/t+P(t)\) zu einer rationalen Funktion mit dem verlangten Abfall bei ∞ machen. Dieses Beispiel betrifft allgemeine Keime; es ist kein Gegenbeispiel zur Kneser-Konvergenzvermutung.

Der Bericht `w_grenzschicht/BERICHT.md`, §5, setzt für die dortige Formel mit der Summe über inverse Äste ausdrücklich ein Polynom f vom Grad d voraus. Das Problemblatt verwendet die kubische Version mit drei Urbildern. Die komprimierte Darstellung im Prompt sollte diese Voraussetzung nicht verlieren. Satz 2 gibt eine exakte Form für Kneser, ohne globale inverse Äste vorauszusetzen.

### Exakte lokale Form mit inversen Ästen

**Status: bewiesen unter den folgenden lokalen Fortsetzungshypothesen.** Sei Ω ein beschränktes Jordan-Gebiet, das 0 enthält; F und m seien in einer Umgebung seines Abschlusses holomorph. Für t ∈ Ω, t ≠ 0,a, t ∉ F(∂Ω), und einfache Urbilder von t in Ω folgt aus dem Residuensatz

\[
T_y(t)=\frac{m(t)Y(t)}t
 +\sum_{\substack{F(z)=t\\z\in\Omega}}\frac{Y(z)}{zF'(z)}
 -U(t)+I_F(t)-B_y(t),
\]

\[
I_F(t)=\frac1{2\pi i}\int_{\partial\Omega}
 \frac{Y(z)}{z(t-F(z))}\,dz.
\tag{10}
\]

Der Integrand hat in 0 Residuum R_y(t) und an einem einfachen Urbild z das Residuum \(-Y(z)/(zF'(z))\); Umstellen ergibt (10). Bei nicht einfachen Urbildern benutzt man stattdessen die entsprechenden Residuen. Dies formuliert die lokale nichtlokale Randwertgleichung exakt. Eine polynomiale Schranke für I_F−B_y auf der nötigen globalen Kurve ist damit noch nicht gewonnen.

## 3. Koeffizientenschranke (K2): Rechenversuch bleibt offen

**Status: kein neues Zertifikat.** Das vorhandene Verfahren `p32_rand_ableitungen.py` wurde untersucht. Seine alten Logs enthalten beim rechten Bogen bereits einen nichtendlichen Integralwert. Solche Werte dürfen nicht als Zertifikat gelten.

In dieser Bearbeitung wurden lokale Varianten mit expliziter Zurückweisung nichtendlicher Intervalle, skalierten Poincaré-Jets und einer schärferen Cauchy-Schranke getestet. Es entstand noch keine abgeschlossene Schranke für das gesamte Randintegral. Insbesondere wird weder C₆ ≤ 300 noch die vollständige Bulk-Schranke als neu bewiesen ausgegeben.

Die theoretische Implikation lautet weiterhin: Wenn \(\partial_\theta^5 f_K(\rho e^{i\theta})\) absolut stetig ist und
\[
\frac1{2\pi}\int_0^{2\pi}|\partial_\theta^6 f_K(\rho e^{i\theta})|\,d\theta\le C_6,
\]
dann liefern sechs partielle Integrationen \(|a_k\rho^k|\le C_6k^{-6}\). Diese Implikation stand bereits im Projekt und ist kein neues Ergebnis dieser Bearbeitung.

## 4. Reproduzierbarkeit und Literaturstatus

Die Datei `ergebnisse/lemma_checks.json` enthält die eigenen symbolischen und numerischen Gegenprüfungen. Das Skript `experimente/verify_lemmas.py` erzeugt sie erneut; `experimente/alpha_check.py` reproduziert die Intervallprüfung von α. Beide lesen `../gemeinsam/kneser_halbexp_koeff_zert_640.json` relativ zu diesem Forschungsordner und schreiben nur nach `ergebnisse/`. Der letzte, unvollständige K2-Code steht in `experimente/k2_boundary.py`; Versuchsgeschichte und offene Prüfungen stehen in `VERSUCHE.md`. Aufrufe und Abhängigkeiten sind in `README.md` dokumentiert. Der mathematische Beweis von Satz 1 benötigt keine numerische Rechnung.

Für die hier bewiesenen Sätze wurde kein externer Konvergenzsatz verwendet. Die Recherche zu Randwertproblemen mit Carleman-Verschiebung ergab beispielsweise [Guo–Zhang, Riemann–Hilbert problems with Carleman shift and conjugation for polynomial Beltrami equations](https://www.tandfonline.com/doi/abs/10.1080/17476933.2019.1689564). Geprüft wurde nur der Verlagsabstract, nicht ein auf unsere Geometrie anwendbarer Satz. Daraus wird ausdrücklich keine Übertragbarkeit auf (H₀) abgeleitet.

## Abschließende drei Punkte

1. **Neu hergeleitet und bewiesen:** Eine explizite, gleichmäßig auf \(\overline D\) konvergierende Folge reeller Polynome mit Grad ≤ 5N und exakt erfüllten Kompositionsgleichungen bis Ordnung N; außerdem die rationale adjungierte Darstellung für beliebige holomorphe Keime einschließlich ihrer analytischen Korrektur.
2. **Präzise verbleibende Lücke:** Die zusätzlichen Koeffizienten oberhalb von N sind in Satz 1 wesentlich für die Konstruktion. Ihre Beseitigung und die Kontrolle der ursprünglichen quadratischen Abschnitte J_N sind offen. Satz 2 liefert keine polynomiale Schranke für die Eckmoden. (H₀), seine kohärente Familienerweiterung und die Identifikation mit dem fortgesetzten Grad-N-Zweig bleiben unbewiesen.
3. **Ein aussichtsreicher nächster Schritt:** Für die exakte Gleichung (5)/(10) eine gewichtete Cauchy-Operatorformulierung auf einer festen Kontur aufbauen und den Term I_F−B_y dort mit polynomialer Genauigkeit kontrollieren. Das würde den Übergang vom kubischen Modell zur echten Kneser-Funktion präzisieren und die benötigte starke Eckabschätzung direkt adressieren.
