# Agent 10 — Neues Gebiet: der Erzeuger der Tetration (Fluss von exp) und zufällige Türme

## Kurzfassung (einfache Sprache)

**Teil A — Tetration als Strömung.** Eine Tetration mit gebrochenen Höhen ist dasselbe wie eine Strömung auf der Zahlengeraden:
Nach genau einer Zeiteinheit ist jeder Punkt x bei e^x angekommen. Die Strömung hat ein Geschwindigkeitsfeld J(x) = sexp′(slog x)
(der „Erzeuger“ oder die Julia-Funktion). Wir haben gefragt, wie dieses Feld in der komplexen Ebene aussieht. Ergebnisse:

1. **Es gibt kein überall glattes Feld.** Kein von null verschiedenes ganzes oder meromorphes J hat e^x als Zeit-1-Abbildung (bewiesen, für jede Basis).
2. **An jedem Fixpunkt ist das Feld starr:** Ist J dort holomorph, so ist es bis auf einen Faktor das reguläre Feld (bewiesen).
3. **Ein Feld kann nicht an zwei Fixpunkten gleichzeitig glatt sein**, ausser die Multiplikatoren sind „resonant“ (λ_p = λ_q^m). Bei Möbius-Abbildungen
   gilt die Resonanz immer (λ_pλ_q = 1) — das erklärt, warum dort beide Fixpunkte dieselbe Iteration liefern. Für b^x über η ist sie unmöglich,
   unter η gibt es genau eine Resonanzbasis b_k je k ≥ 2 (b₂ = 1.04480434792877928…). Dort weichen die beiden Iterationen trotzdem voneinander ab (numerisch).
4. **Folge für Kneser:** Jedes reelle Feld muss an L oder L̄ singulär sein. Für Knesers Feld haben wir die Singularität genau bestimmt:
   J_K − J_L ≈ C·(w−L)^{1+2πi/L}. Bei Basis e ist J_K am Fixpunkt **genau 5-mal stetig differenzierbar, aber nicht 6-mal** (Exponent 5.44695…, auf 4 Stellen gemessen).
5. **Das reguläre Feld lebt über C ohne die Bahn der Null** {0, 1, e, e^e, …}: Es lässt sich längs jedes Weges fortsetzen, der diese Punkte meidet (bewiesen);
   um 0 und um 1 verzweigt es tatsächlich (numerisch, 20 Stellen).

**Teil B — Türme mit zufälligen Basen.** Für Basen ≥ 1 divergiert der Zufallsturm fast sicher, sobald Basen über η möglich sind und Basen nahe 1
nicht zu häufig (bewiesen). Für Basen in (0,1) hat der Turm höchstens zwei Häufungspunkte, und Konvergenz ist ein 0-1-Ereignis (bewiesen).
Überraschend: Mit gleichverteilten Basen in (0, β) **konvergiert der Turm auch dann, wenn alle Basen unter Eulers Schranke e^{−e} liegen**
(β = 0.03 … 0.05, bei 0.03 erst nach sehr vielen Stufen; numerisch) — der Zufall zerstört Eulers Zweierschwingung. Die „Verschmelzung“ geschieht schlagartig mit einer Rate, die für kleines β
extrem schnell fällt. Vermutung: Konvergenz fast sicher für jedes β > 0.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| A1 | Kein ganzes J ≢ 0 mit J(b^z) = ln b · b^z · J(z) (jedes b ∉ {0,1}) | **bewiesen** | Wachstumsabschätzung M(r) ≤ \|c\| r M(√(log²r+π²)/\|c\|) | Aussage vermutlich bekannt (Folklore); dieser Beweis vermutlich neu |
| A1′ | Auch kein meromorphes J ≢ 0 | **bewiesen** | Pole sind rückwärts-invariant, Rückwärtsbahnen häufen sich auf der Julia-Menge | vermutlich neu in dieser Form |
| A2 | J holomorph an hyperbolischem Fixpunkt p ⇒ J = c·J_p (J_p = log λ·σ/σ′), J(p) = 0 | **bewiesen** | Koeffizientenvergleich in Schröder-Koordinaten | klassisch |
| A3 | Zwei hyperbolische Fixpunkte p, q (einer abstossend): ein gemeinsamer holomorpher Erzeuger auf zusammenhängendem D ∋ p, q ⇒ λ_p = λ_q^m für ein m ∈ Z | **bewiesen** | Periodenargument mit der ganzen Superfunktion | Karlin–McGregor 1968 zeigen Verwandtes (nur Möbius); Resonanzform und Beweis vermutlich neu |
| A3a | b > η: kein Erzeuger holomorph um L und L̄ zugleich; jeder reelle Erzeuger ist an L (bzw. L̄) singulär | **bewiesen** (Folgerung) | m = 1 und Im log λ ∉ πZ | vermutlich neu (als Satz über Erzeuger) |
| A3b | 1 < b < η: λ_aλ_A < 1 immer; Resonanzbasen b_k (λ_aλ_A^k = 1) existieren genau einmal je k ≥ 2: b₂ = 1.0448043479287792775, b₃ = 1.0015621748618726604, b₄ = 1.000038279179596386 | **bewiesen** (λ_aλ_A<1, Existenz) / **numerisch** (Werte, 20 St.) | Ungleichung ln(1+x) < x/√(1+x); findroot | vermutlich neu |
| A3c | An b₂, b₃ unterscheiden sich die regulären Iterationen an a und A trotzdem: \|ĉ₁\| = 1.4010369·10⁻⁴ bzw. 3.6650371·10⁻³, glatt durch die Resonanz | **numerisch** | dps 60/N 48 gegen dps 90/N 64: 10 St. gleich (rel. 10⁻⁵¹); Kontrolle √2: 3.626087773·10⁻²⁵ = Graphwert | neu |
| A3d | Möbius: λ_pλ_q = 1 ist genau die Resonanz m = −1 — passt zu it.satz.moebius_iteration | **bewiesen** | — | Einordnung neu |
| A4 | Knesers Erzeuger: J_K(w) = J_L(w)(1+θ′(s)), s = slog_K w; an L: J_K − J_L = 2πi a₁ e^{2πis} J_L (1+O(\|w−L\|^{4.45})) ~ C (w−L)^{α}, α = 1+2πi/L = 5.44695072+1.05794i | **bewiesen** (Formel) / **numerisch** (Exponent) | lokale Steigung 5.4466–5.4475 bei r = 10⁻³→10⁻⁴; dps 40/60 | vermutlich neu |
| A4a | J_K analytisch an L ⇔ θ konstant; also (da θ ≠ konst.) echter transzendenter Verzweigungspunkt | **bewiesen** | A2 + Periodizität von θ | vermutlich neu |
| A4b | Basis e: J_K ist am Fixpunkt C⁵, nicht C⁶; allgemein C^{⌊Re α⌋}, Tabelle b = 1.45: C⁴⁵ … b = 10⁶: C³ | **bewiesen** unter a₁(b) ≠ 0 (für e, 2 belegt) | s4_ausgabe.txt | neu |
| A4c | Umlauf um L multipliziert den singulären Koeffizienten mit e^{−4π²/L} (Betrag 1.2954·10⁻³) | **hergeleitet** + **numerisch** | gemessen 1.2975·10⁻³ (Abweichung = O(r)-Korrektur) | neu |
| A5 | σ und J_L = Lσ/σ′ lassen sich längs jedes Weges in C∖P fortsetzen, P = {0,1,e,e^e,…} | **bewiesen** | hyperbolische Metrik, iterierte Logarithmen | vermutlich neu (für exp als Satz über Erzeuger) |
| A5a | 0 und 1 sind echte Verzweigungspunkte von J_L; J_L(L+2πi) (längs der Senkrechten) = −14.322345969901284658 + 10.705167958195350901 i ≠ 0 | **numerisch** | zwei Einstellungen (dps 30/36, Tiefe 60/75, 4000/2500 Schritte), 20 St. gleich; Kontrollschleifen um L und 3: Sprung < 10⁻³⁰ | neu |
| A6 | Keine Resonanz λ_j = λ_k^m zwischen Fixpunkten von exp für \|j\|,\|k\| ≤ 40, m ∈ [−12,12]∖{0,1}: kleinster rel. Abstand 0.424 | **numerisch** | s4 | — |
| B1 | Basen in (0,1): Durchmesser D_N von f₁∘…∘f_N([0,1]) fällt monoton; Turm konvergiert ⇔ D_N → 0; sonst gerade/ungerade Teilfolgen gegen zwei Werte; bei u.i.v. Basen ist Konvergenz ein 0-1-Ereignis | **bewiesen** | geschachtelte Intervalle, Ergodizität | Teil 1 klassisch, 0-1 vermutlich bekannt |
| B2 | Robuste Zweierschwingung: zu jedem a₀ < e^{−e} gibt es ε > 0, so dass **jede** Basenfolge in [a₀−ε, a₀+ε] divergiert | **bewiesen** | Fallenpaar um den 2-Zyklus | vermutlich bekannt (Barrow-Umfeld) |
| B3 | Basen u.i.v. ≥ 1, P(ln a > 1/e) > 0, P(ln a ≤ ε) ≤ Cε^κ ⇒ Turm = ∞ fast sicher | **bewiesen** | Vorwärtskette, Martingal | vermutlich neu (Zufallsturm) |
| B4 | Basen ~ U(0,β): Turm konvergiert für β ≥ 0.03 (also auch ganz unter e^{−e} = 0.0660) in allen Läufen mit N ≤ 256 000; Verschmelzung schlagartig, Überlebensrate exponentiell in N: Rate ≈ 3.4·10⁻⁵ (β=0.035), 4·10⁻⁶ (0.03), 10⁻⁷ (0.025) | **numerisch** | 2 Saaten, M = 200–1000, Exponentialgesetz an 3 N-Werten | neu |
| B5 | Vermutung: für jedes β > 0 konvergiert der U(0,β)-Turm fast sicher (Zufall zerstört Eulers Periode-2-Verzweigung, Rate ≈ e^{−c/β}) | **Vermutung** | B4; Grund: Basen beliebig nahe 0 zerstören jedes Fallenpaar (B2 gilt nicht) | neu |
| B6 | Vermutung (Schwelle mit slog): Basen ≥ 1 mit P(ln a ≤ ε) = κ/slog(1/ε) nahe 0 und Masse über η: Konvergenz f.s. ⇔ κ > 1 | **Vermutung** (heuristisch hergeleitet) | Reduktion auf Kartenhaus-Kette in slog-Koordinaten | neu |

## Teil A — Der Erzeuger J (Julia-Funktion) von exp

Bezeichnungen: f = exp_b, c = ln b; ein **Erzeuger** ist eine holomorphe Lösung der Julia-Gleichung J(f(z)) = f′(z) J(z) (it.lem.julia).
Jede Tetration T liefert J_T = T′∘slog_T. Basis e: L = 0.3181315052…+1.3372357014…i, Multiplikator λ = L; σ Schröder-Funktion (σ∘exp = Lσ, σ′(L)=1),
Ψ = σ⁻¹ ganz, F_L(s) = Ψ(e^{Ls}) reguläre Superfunktion, h = F_L⁻¹, J_L = 1/h′ = Lσ/σ′.

**Warum das Gebiet zählt.** Die Frage „welche gebrochene Iteration?“ ist äquivalent zur Wahl von J, und die Julia-Gleichung ist *linear*.
Damit lassen sich Fragen der Fortsetzung als Fragen über ein einziges holomorphes Objekt stellen: Wo darf J singulär sein, wie stark, und welche Singularitäten
sind erzwungen? A1–A3 zeigen: Singularitäten sind unvermeidlich und an Fixpunktpaaren erzwungen; A4 misst Knesers erzwungene Singularität genau;
A5 verbindet den Erzeuger mit dem zentralen Objekt der transzendenten Dynamik, der Bahn des singulären Werts 0.

### A1 Kein ganzer / meromorpher Erzeuger

*Beweis (ganz).* Sei M(r) = max_{|w|=r}|J|. Der Kreis |w| = R ist {e^{cz} : cz = log R + iy, |y| ≤ π}, und |z| ≤ √(log²R+π²)/|c|.
Aus J(e^{cz}) = c e^{cz} J(z) folgt M(R) ≤ |c| R · M(√(log²R+π²)/|c|). Iteriert man r_{j+1} = √(log²r_j+π²)/|c|, bis r_N unter eine feste Schranke fällt
(N ≤ log* R + O(1)), ergibt sich M(R) ≤ |c|^N R (log R + C)^N M(C′) = O(R^{1+ε}). Also ist J ein Polynom vom Grad ≤ 1, J = αw + β;
α e^{cz} + β = c e^{cz}(αz + β) für alle z erzwingt α = β = 0 (lineare Unabhängigkeit von 1, e^{cz}, z e^{cz}). ∎

*Beweis (meromorph).* Ist p ein Pol und e^{cq} = p, so ist J(q) = J(e^{cq})/(c e^{cq}) = ∞: alle Urbilder eines Pols sind Pole.
Für p ≠ 0 häuft sich die Rückwärtsbahn von p an jedem Punkt der (nichtleeren) Julia-Menge (0 ist der einzige mögliche Ausnahmewert) — Widerspruch zur Diskretheit.
Ein Pol bei 0 erzwingt einen Pol bei f(0) = 1 ≠ 0. Also ist J ganz und nach dem ersten Teil null. ∎

### A2 Starrheit am Fixpunkt
J(p) = λJ(p) ⇒ J(p) = 0. In Schröder-Koordinaten u = σ(z) wird J zu K(u) mit K(λu) = λK(u); Koeffizientenvergleich k_n(λ^n − λ) = 0 ⇒ K = k₁u,
also J = c·J_p. (Die Julia-Gleichung allein lässt jedes c zu; Zeit-1 = f verlangt c ∈ {1 + 2πik/log λ}.) ∎

### A3 Zwei Fixpunkte schliessen sich aus (ausser bei Resonanz)
*Satz.* f ganz, p, q hyperbolische Fixpunkte, q abstossend (F_q ganz), Multiplikatoren λ_p, λ_q, ℓ = log λ. Ist J ≢ 0 holomorph auf zusammenhängendem D ∋ p, q
und erfüllt nahe p und q die Julia-Gleichung, dann gilt λ_p = λ_q^m für ein m ∈ Z∖{0}.

*Beweis.* Nach A2: J = c₁J_p nahe p, J = c₂J_q nahe q. Integriert man 1/J längs eines Weges in D (Nullstellen gemieden), wird die Fortsetzung von h_p
zu A h_q + B mit A = c₂/c₁. Umkehren: F_p(As + B) = F_q(s) auf einer offenen Menge, also überall (F_q ganz; F_p wird dadurch ganz fortgesetzt).
F_q hat die Periode 2πi/ℓ_q, und die Periodengruppe von F_p ist genau (2πi/ℓ_p)Z (bei Re(ℓ_p s) → −∞ ist Ψ_p injektiv). Also A = mℓ_q/ℓ_p.
Aus F_q(s+1) = f(F_q(s)) und F_p(w+1) = f(F_p(w)) folgt F_p(w + A) = F_p(w + 1), also A − 1 ∈ (2πi/ℓ_p)Z: mℓ_q − ℓ_p ∈ 2πiZ, d. h. λ_q^m = λ_p. ∎

*Folgerungen.*
- **b > η, Paar L, L̄:** ℓ̄ = mℓ + 2πin ⇒ Realteil (m−1)log|λ| = 0 ⇒ m = 1 ⇒ Im ℓ ∈ πZ, unmöglich für nichtreelles λ. Kein Erzeuger ist um L und L̄ zugleich holomorph.
  Ein reeller Erzeuger auf einem zusammenhängenden, zur reellen Achse symmetrischen Gebiet, das L enthält, ist durch Spiegelung auch um L̄ holomorph — also unmöglich.
  **Knesers Singularität an L ist erzwungen.**
- **1 < b < η, a < A:** λ_a = ln a, λ_A = ln A; mit a = (1+1/t)^t, A = (1+1/t)^{t+1} ist λ_aλ_A = t(t+1)ln²(1+1/t) < 1 (Ungleichung ln(1+x) < x/√(1+x)),
  also nur Resonanzen λ_aλ_A^k = 1 mit k ≥ 2; je k genau eine Lösung (s3: b₂, b₃, b₄ oben). Für alle anderen b ist die Verschiedenheit der beiden regulären Iterationen
  (kl.beob.zwei_abel, „Bummer“) damit **bewiesen**. An b₂, b₃ greift der Beweis nicht; die Rechnung (A3c) zeigt aber auch dort deutliche Abweichung.
  Nebenbei: |ĉ₁|·e^{2π²/|ln λ_a|} = 0.091 (√2), 0.085 (b₂), 0.078 (b₃) — Stütze für kl.verm.streifenbreite über 20 Grössenordnungen, mit langsam driftender Konstante.
- **Möbius:** λ_pλ_q = 1 immer — genau die Resonanz m = −1. Der Satz erklärt also, warum die Übereinstimmung bei Möbius (it.satz.moebius_iteration) möglich ist.

### A4 Knesers Erzeuger am Fixpunkt
Aus sexp_K = F_L∘(id + θ), θ(s) = Σ a_n e^{2πins}: slog_K = (id+θ)⁻¹∘h und **J_K(w) = J_L(w)·(1 + θ′(slog_K w))**.
Nahe L ist e^{2πis} ≈ e^{−2πia₀}σ(w)^{2πi/L}, also J_K = L(w−L)(1+O(w−L)) + 2πi a₁ e^{−2πia₀}L (w−L)^{1+2πi/L} + … (Zweig von a₀ passend zum Blatt).
Ist J_K an L analytisch, so ist J_K = cJ_L (A2), also 1+θ′ konstant, also θ konstant (periodisch) — falsch. Re α = 1 + 2π Im L/|L|² = 5.44695072…:
alle Ableitungen bis Ordnung 5 sind stetig, die 6. wächst wie |w−L|^{−0.553}. Numerik (s1): Verhältnis D/(erstes θ-Glied) − 1 = O(r^{4.45}), Steigungen 5.4466–5.4475.
Ein Blattwechsel (Umlauf um L) ersetzt s durch ≈ s + 2πi/L und multipliziert das Glied mit e^{−4π²/L}; gemessen 1.2975·10⁻³ gegen 1.2954·10⁻³.

### A5 Der reguläre Erzeuger über C∖P
J_L∘Ψ(u) = L u Ψ′(u) ist ganz, also genügt es, σ = Ψ⁻¹ fortzusetzen. Ψ hat keine kritischen Punkte (Ψ′(λu)λ = e^{Ψ(u)}Ψ′(u) würde Nullstellen gegen 0 häufen).
*Beweis der Fortsetzung:* Für einen Weg γ in Ω = C∖P lassen sich alle iterierten Logarithmen γ_n = log^n∘γ stetig wählen (log^k γ trifft nie 0, weil γ die Punkte exp^k(0) meidet).
exp: C∖exp⁻¹(P) → Ω ist eine Überlagerung, und C∖exp⁻¹(P) ⊊ Ω; also verkürzt jeder Logarithmus die hyperbolische Länge, auf Kompakta um einen festen Faktor c < 1.
Die Startpunkte log^n(γ(0)) → L, also γ_n → L gleichmässig, und σ(γ(t)) = L^nσ_loc(γ_n(t)). ∎
Numerik (s2): Umlauf um 0 ändert J_L(1/2) von 1.39814+0.02656i auf −6.56052+2.25811i, Umlauf um 1 ändert J_L(3/2) von 2.05008−0.30198i auf −8.27091+2.37225i;
Umläufe um L und um 3 (kein Punkt von P) ändern nichts (< 10⁻³⁰). Bei trivialer Monodromie um 0 müsste J_L(L+2πi) = e^{−L}J_L(L) = 0 sein; gemessen ≠ 0 (A5a).

## Teil B — Zufällige Türme

T_N = f_{a₁}∘…∘f_{a_N}(1), f_a(x) = a^x, Basen unabhängig gleichverteilt (u.i.v.).

- **B1** (bewiesen). Für a_n ∈ (0,1) liegen alle T_M (M ≥ N) in I_N = f₁∘…∘f_N([0,1]); I_{N+1} ⊂ I_N. F_N(1) ist ein Endpunkt, je nach Parität oben oder unten,
  also T_{2n} → v, T_{2n+1} → u mit [u,v] = ∩I_N. Konvergenz ⇔ D_∞ = 0. Da f₁ auf [0,1] Ableitungen zwischen a₁|ln a₁| > 0 und |ln a₁| hat, ist {D_∞ = 0} verschiebungsinvariant ⇒ 0-1.
- **B2** (bewiesen). Für a₀ < e^{−e} ist der 2-Zyklus p<q von f_{a₀} anziehend: |f′(p)f′(q)| < 1. Radien r_p, r_q mit |f′(p)|r_p < r_q, |f′(q)|r_q < r_p ergeben Intervalle
  I_p, I_q mit f_a(I_p) ⊂ I_q, f_a(I_q) ⊂ I_p für alle a nahe a₀. Die innere Bahn von 1 landet nach festen M Schritten in der Falle; T_N und T_{N+1} landen in verschiedenen Fallen,
  also |T_N − T_{N+1}| ≥ dist(I_p, I_q) > 0. ∎
- **B3** (bewiesen). Vorwärtskette X_{n+1} = exp(c_{n+1}X_n) hat dieselbe Verteilung wie T_n, und T_n ↑. Mit p = P(c ≥ c*) > 0, c* > 1/e: nach r guten Schritten ist X ≥ exp_{c*}^r(1);
  von dort fällt X_{k+1} ≥ X_k² nur mit Wahrscheinlichkeit ≤ C(2 log X/X)^κ aus, summierbar längs X^{2^k}. Daher h(x) = P_x(X_n → ∞) ≥ η₀ > 0 gleichmässig,
  Martingal-Konvergenz gibt X_n → ∞ fast sicher, also P(T_∞ ≤ M) = lim P(X_n ≤ M) = 0. ∎ (Mit Atom bei a = 1 konvergiert der Turm dagegen immer.)
- **B4/B5** (numerisch / Vermutung). z1–z3: Für U(0, β) konvergiert der Turm bei β = 0.035–0.05 schon nach N ≤ 64 000 in allen Läufen; bei β = 0.03 überleben
  92 % / 77 % / 40 % (N = 16k / 64k / 256k) — Exponentialgesetz mit Rate ≈ 4·10⁻⁶; bei 0.025 Rate ≈ 10⁻⁷; bei 0.015–0.02 keine Verschmelzung bis 256 000.
  Der Durchmesser springt von ≈ 1 direkt auf < 10⁻³⁰⁰ (Verschmelzung ist ein seltenes Einzelereignis). B2 greift hier nicht, weil Basen beliebig nahe 0 jedes Fallenpaar sprengen.
  Hinweis: der naive Ljapunow-Exponent χ = E log(T|ln a₁|) ist hier kein Kriterium (negativ auf beiden Seiten, wie bei einer Gabelverzweigung).
- **B6** (Vermutung). In slog-Koordinaten Y = slog X steigt Y je Schritt um 1 + o(1), solange c nicht kleiner als ≈ 1/X ist; bei langsam variierendem G(ε) = P(c ≤ ε)
  springt die Kette dann auf beschränkte Höhe (Kartenhaus-Kette mit Rückfallwahrscheinlichkeit q_Y ≈ G(1/exp^Y(1))). Positiv rekurrent ⇔ Σ_k Π_{j≤k}(1−q_j) < ∞;
  für G = κ/slog(1/ε) ist q_j ≈ κ/j, also Konvergenz ⇔ κ > 1. Lücke: Fehlerterme der slog-Näherung und Rückfallverteilung nicht streng kontrolliert.

## Vorschläge für Grapheinträge
- Kein ganzer oder meromorpher Erzeuger von b^z — satz — bewiesen
- Starrheit des Erzeugers an einem hyperbolischen Fixpunkt — lemma — bewiesen
- Zwei Fixpunkte, ein Erzeuger nur bei Resonanz λ_p = λ_q^m — satz — bewiesen
- Jeder reelle Erzeuger von exp_b (b > η) ist an L singulär — folgerung — bewiesen
- Resonanzbasen b_k unter η (λ_aλ_A^k = 1), λ_aλ_A < 1 — satz — bewiesen (Werte numerisch)
- Abweichung der regulären Iterationen an den Resonanzbasen b₂, b₃ — beobachtung — numerisch
- Knesers Erzeuger: J_K = J_L(1+θ′∘slog_K) und Verzweigung (w−L)^{1+2πi/L} — satz — bewiesen (Exponent numerisch bestätigt)
- Glattheitsgrad von J_K am Fixpunkt (Basis e: C⁵, nicht C⁶; Tabelle über b) — satz — bewiesen unter a₁ ≠ 0
- Regulärer Erzeuger setzt sich längs C∖{0,1,e,…} fort; 0 und 1 verzweigen — satz + beobachtung — bewiesen / numerisch
- Zufallsturm (0,1): zwei Häufungspunkte, 0-1-Gesetz — satz — bewiesen
- Robuste Zweierschwingung unter e^{−e} — satz — bewiesen
- Zufallsturm mit Basen ≥ 1 divergiert f.s. (Potenzschranke bei 1) — satz — bewiesen
- U(0,β)-Turm konvergiert auch unter e^{−e} — beobachtung — numerisch; f.s. für jedes β > 0 — vermutung
- slog-Schwelle κ > 1 für Zufallstürme mit Basen ≥ 1 — vermutung

## Offene Fragen
1. Monodromiegruppe von σ (bzw. J_L) über C∖P: welche Darstellung von π₁(C∖P)? Ist jeder Punkt exp^k(0) ein echter Verzweigungspunkt?
2. Wo genau liegen die Singularitäten von J_K auf seiner Riemannfläche (Bild der Rückwärtsbahn von L — dicht)? Ist J_K auf dem Hauptblatt nur an L, L̄ und den Schnitten singulär?
3. Resonanzbasen b_k: gibt es einen anderen Grund, warum die Iterationen an a und A dort verschieden sind? (Satz von Karlin–McGregor prüfen: deckt er den resonanten Fall?)
4. Fixpunktpaare L_j, L_k von exp: ist L_j = L_k^m je möglich? (Schanuel-artige Frage; numerisch ausgeschlossen für kleine Indizes.)
5. Beweis von B5 (Rate e^{−c/β}?) und von B6; komplexe Zufallsbasen (zufälliges Shell–Thron-Gebiet).

## Reproduktion
Python: `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`. Kein Motor nötig; θ-Koeffizienten werden nur gelesen aus
`forschung/identitaeten_2026_10_14/etappe2/daten/a_e_80_0.003_16384.json` (über `lib10.py`, Kopie von etappe3/lib3.py).
- `s1_erzeuger_bei_L.py` → `s1_ausgabe.txt` (A4; ~5 min)
- `s2_monodromie_regulaer.py 1` / `2` → `s2_ausgabe_1.txt`, `s2_ausgabe_2.txt` (A5a; je ~5–8 min)
- `s3_resonanz_unter_eta.py` → `s3_ausgabe.txt` (A3b/c)
- `s4_tabellen.py` → `s4_ausgabe.txt` (A4b, A6)
- `z1_zufallsturm.py`, `z2_schwelle.py`, `z3_durchmesser.py` → `z*_ausgabe.txt` (B4)
