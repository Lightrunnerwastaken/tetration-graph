# Gebiet 16 — Julia-Menge von eᶻ und wie Kneser-sexp die Ebene abbildet

**In einfacher Sprache.** Bei eᶻ ist die Julia-Menge die ganze Ebene; interessant ist darum die Menge I der *flüchtenden* Punkte
(Bahn geht gegen ∞), die aus „Haaren“ besteht. Ein Punkt z wird von Kneser-sexp F genau dann auf einen flüchtenden Punkt
abgebildet, wenn F(z+n) → ∞; diese Menge ist im ganzen Streifen dicht, hat aber Fläche null. Neu gefunden: Oben in der
Halbebene ist das Muster der Flucht **doppelt periodisch** — ausser der Periode 1 gibt es die schräge Periode
ω = 2πi/L − 4 ≈ 0.447 + 1.058i, und die Fluchtzeit wächst pro Schritt ω um genau 4. Ausserdem gezeigt (mit Literaturlücke):
**Kein Haar von eᶻ endet am Fixpunkt L oder häuft sich dort**; darum steigt kein zurückgezogenes Haar in einem senkrechten
Streifen bis Im z → ∞. Senkrechte Linien gehen unter F in Bögen von L̄ nach L über (keine Kreise), und der ganze Streifen
landet in einer Linse innerhalb |Im w| < 1.756.

Zeit: Start 13:22, Rechnungen fertig 13:43 (Uhr aus `date`).

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| 1 | E := F⁻¹(I(exp)) ist 1-invariant (z∈E ⇔ z+1∈E), und z∈E ⇔ F(z+n)→∞ | bewiesen | aus F(z+1)=exp F(z) | trivial, aber nirgends notiert |
| 2 | E ist dicht, hat Lebesgue-Mass 0; enthält ℝ∩Def(F) und alle F⁻¹(exp⁻ⁿ(ℝ)) | hergeleitet mit Lücke | Dichte: I dicht in J=ℂ; Mass 0: Rees/Lyubich (lit_fp_1); Lücke: Mass-0-Aussage nur sekundär geprüft | bekannt übertragen |
| 3 | Auf H gilt E∩H = h⁻¹(E_L) mit h(z)=z+θ(z) (Kneser), E_L = F_L⁻¹(I) ist **Λ-periodisch**, Λ = ℤ + (2πi/L)ℤ | bewiesen | F_L(u)=Ψ(e^{Lu}) erfüllt F_L(u+1)=exp F_L(u), F_L(u+2πi/L)=F_L(u) | vermutlich neu (als Aussage über sexp) |
| 4 | Zweiseitiger Eintrittsindex N(z)=min{n∈ℤ: Re F(z+n)>R}: N(z+1)=N(z)−1 exakt; N_L(u)+Re(Lu)/Re L ist Λ-periodisch; also N(z+ω)=N(z)+4 asymptotisch, ω=2πi/L−4 | bewiesen (exakt für F_L), numerisch für Kneser | Gitter 200×1201, R=50: Anteil N(z+ω)=N(z)+4 je Im-Band: 0.74 (0–0.5), 0.93 (0.5–1), 0.996 (1–1.5), 0.9999 (2–2.5), **1.0 ab Im z ≥ 3**; Gegenprobe mit Streifenaustritt T statt N: 0.91 / 0.9996 / 1.0 ab Im ≥ 2; zweite Einstellung R=20, Gitter 150×900: 0.80 / 0.95 / 0.9975 / 0.9999 / 1.0 ab Im ≥ 3 | vermutlich neu |
| 5 | Steigung: N(z) ≈ (Im L/Re L)·Im z − Re z + periodisch, Im L/Re L = 4.2034 | bewiesen (Formel), numerisch | Ausgleichsgerade über 1<Im z<7.4: 4.244 (1 % Abweichung durch periodischen Anteil) | vermutlich neu |
| 6 | F bildet den Streifen 0≤Re z<1 in die Linse Ω zwischen F(iℝ) und exp(F(iℝ)) ab; Ω ⊂ {|Im w|<1.756} ⊂ S₀ = {|Im w|<π} | numerisch | Taylor-Auswerter, Abweichung zu gp ≤ 2·10⁻¹¹; max |Im F| = 1.7558 | Zahl neu |
| 7 | Im F(z) > 0 auf dem oberen Grundstreifen, sogar Im F/Im z ≥ 0.178; daher F(z) ∉ ℝ für 0<Re z<2, Im z>0 | numerisch (+ bewiesen für den Schritt Streifen→Streifen+1: 0<Im F<π ⇒ Im eᶠ>0) | 400×800 Punkte bis Im z = 10⁻⁴ | neu |
| 8 | Die ersten nichtreellen Urbilder der reellen Achse (F(z) < 0) beginnen bei Re z ≈ 2.52 (Im F(z) = π zuerst bei z ≈ 1.516+0.635i) | numerisch | Gitterauflösung 0.0025 | neu |
| 9 | **Kein dynamischer Strahl (Haar) von exp hat L oder L̄ in seinem Abschluss**, auch keine zerbrochenen Urbild-Strahlen von ℝ | hergeleitet mit Lücke | Beweis unten; Lücke: Struktursätze von Schleicher–Zimmer (Strahlen disjunkt, Schwanz in S_{s₀}, g_{0^∞}=ℝ) im Fall „singulärer Wert flüchtet“ nicht am Original geprüft | vermutlich bekannt (Devaney 1993, Rempe 2006 schliessen diesen Fall aus) |
| 10 | Folgerung: Jedes F-Urbild eines Haars ist in jedem senkrechten Streifen a≤Re z≤b nach oben beschränkt; in u-Koordinaten bleibt Re(Lu) auf jedem Strahl-Urbild nach unten beschränkt | hergeleitet mit Lücke (wie 9) | F→L gleichmässig in senkrechten Streifen | neu |
| 11 | Bilder senkrechter Linien F(x₀+iy) sind Bögen von L̄ über F(x₀) nach L, symmetrisch zur reellen Achse, aber keine Kreise (beste Kreisanpassung weicht 2.9 % (x₀=0) bis 6.0 % (x₀=¾) des Radius ab); nahe L log-Spiralen (pro Einheit Höhe dreht das Argument um Re L = 0.318 rad, der Abstand zu L schrumpft um den Faktor e^{−Im L} = 0.263) | bewiesen (asymptotisch), Bild | F(z) = L + e^{L(z+θ∞)} + … | bekannt-artig |

## Beweise und Details

**Zu 3.** Ψ sei die Poincaré-Funktion am abstossenden Fixpunkt L (Ψ(0)=L, Ψ'(0)=1, Ψ(Lw)=exp Ψ(w)); sie ist ganz. Wegen e^L = L
ist log L = L, also F_L(u) := Ψ(e^{Lu}) ganz mit F_L(u+1) = exp F_L(u) und F_L(u+2πi/L) = F_L(u). Damit ist
E_L = {u : F_L(u) ∈ I} invariant unter beiden Verschiebungen. Kneser: F = F_L∘h auf H mit h(z) = z+θ(z), θ 1-periodisch,
θ(z) → θ∞ für Im z → ∞. Also E∩H = h⁻¹(E_L), und wegen h(z+ω) − h(z) − ω = θ(z+ω) − θ(z) → 0 ist das Muster oben
asymptotisch Λ-periodisch. Reduzierte Gitterbasis: 1 und ω = 2πi/L − 4 = 0.446951 + 1.057940i.

**Zu 4/5.** N_L(u) := min{n∈ℤ : Re F_L(u+n) > R} (das Minimum existiert, weil F_L(u+n) → L für n → −∞). Dann N_L(u+1) = N_L(u)−1
und N_L(u+P) = N_L(u), P = 2πi/L. Die lineare Funktion A(u) = −Re(Lu)/Re L erfüllt A(u+1) = A(u)−1 (da Re(L·1)/Re L = 1)
und A(u+P) = A(u) (da Re(L·2πi/L) = 0). Also ist N_L − A Λ-periodisch (unbeschränkt nur nahe der abgeschlossenen
Nullmenge der Bahnen, die nie Re > R erreichen). A(z) = −Re z + (Im L/Re L)·Im z. Deutung: Nahe L wächst der Abstand pro
Schritt um |L| = e^{Re L}, und |e^{Lu}| = e^{Re(Lu)}; die Zahl der Schritte bis zum Verlassen der Umgebung ist −Re(Lu)/Re L.
Insbesondere N(z+ω) = N(z) + 4 (weil ω = P − 4). Gleiches gilt für den Streifenaustritt T(z) = min{n : |Im F(z+n)| > π}.

**Zu 9 (Beweisskizze).** Strahlen g_s (Schleicher–Zimmer) sind Kurven aus flüchtenden Punkten mit Adresse s = s₀s₁…, Schwanz
g_s(t) ≈ t + 2πi s₀, exp(g_s) ⊂ g_{σs}, verschiedene Adressen ⇒ disjunkt; g_{0^∞} = ℝ (jeder reelle Punkt flieht, Adresse 0^∞).
Die Geraden Im z = (2k+1)π bestehen aus flüchtenden Punkten, deren Bahn nach zwei Schritten reell ist. (a) Ist s nicht
schliesslich 0, so trifft g_{σⁿs} diese Geraden nie (sonst g_{σⁿ⁺¹s} ∩ ℝ ≠ ∅, also σⁿ⁺¹s = 0^∞, Widerspruch); g_{σⁿs} ist
zusammenhängend, also ganz im offenen Streifen S_{sₙ}. Liegt L im Abschluss von g_s, dann L = expⁿ(L) im Abschluss von
expⁿ(g_s) ⊂ S_{sₙ}; da L im Innern von S₀ liegt, folgt sₙ = 0 für alle n, also s = 0^∞ — aber der Abschluss von ℝ enthält L
nicht. (b) Ist s schliesslich 0 (oder ein zerbrochenes Urbild von ℝ), dann expᵐ(g_s) ⊂ ℝ für ein m, und L = expᵐ(L) läge im
Abschluss von ℝ. Widerspruch. Gleich für L̄. Bemerkung: Das betrifft jeden periodischen Punkt, dessen Bahn ganz in S₀ liegt.
Die Aussage kann nicht auf beliebige zusammenhängende Teilmengen von I erweitert werden (I(exp) ist zusammenhängend und dicht, Rempe).

**Zu 10.** F(z) → L gleichmässig für Im z → ∞ in a ≤ Re z ≤ b (Kneser). Ein Haar-Urbild, das dort beliebig hoch steigt, hätte
ein Bild, das sich bei L häuft — ausgeschlossen durch 9.

**Numerik.** `kneser_feld.py`: Taylor-Reihen von F um 0.5+i·y_c (y_c = 0, 1.5, 3, 4.5, 6), Koeffizienten aus 96 gp-Werten auf Kreisen
(Radius 1.4 bzw. 2.0); Gegenprobe gegen direkte gp-Werte an 9 Punkten: Abweichung ≤ 1.5·10⁻¹¹. (Erster Versuch mit Radius 2.0 am
Zentrum 0.5 gab 2·10⁻⁵ wegen ungenauer gp-Werte nahe Re z = −1.5 — verworfen.) `e1_eintrittsindex.py`: Indizes N, T, Gittertest,
Bilder. `e2_reell.py`: Befunde 7 und 8. `e3_gegenprobe.py`: zweite Einstellung für Befund 4 und Kreistest (Befund 11). Alle Zahlen in `cache/e1_zahlen.json`. Rechnung in double precision; die Fluchtbilder
sind wie alle Escape-time-Bilder an den Rändern der Niveaumengen empfindlich (Chaos), die Gitterstatistik ist es nicht.

**Bilder.**
- `bilder/e1_zebene.png`: links N(z) mod 6 in der z-Ebene (Re z ∈ [−1,3], |Im z| ≤ 7.5), rechts T(z) mod 6. Die weissen Kreuze
  sind 3i + mω: das schräge Gitter ist sichtbar. Unten spiegelbildlich (Gitter mit L̄). Nahe der reellen Achse bricht das Gitter,
  dort ist h weit weg von einer Verschiebung (θ ist an den ganzen Zahlen unbeschränkt, siehe unten).
- `bilder/e1_webene.png`: dynamische Ebene von exp mit Fluchtzeit; farbig die Bilder der senkrechten Linien Re z = 0, ¼, ½, ¾, 1
  (Bögen von L̄ nach L, keine Kreise), gestrichelt Bilder waagrechter Strecken (Spiralstücke um L). Die Linse liegt weit innerhalb S₀.

**Nebenbefund (bewiesen, Kneser-bekannt).** θ ist an jeder ganzen Zahl unbeschränkt: Für z → 0 in H gilt F_L(h(z)−1) = F(z−1) → 0,
und F_L lässt 0 aus; also h(z) → ∞. Das Bild h((n,n+1)) der reellen Achse im u-Bild ist ein Bogen von ∞ nach ∞, ein Hochheben
des reellen Fluchtstücks (F(n),F(n+1)) — also selbst Teil von E_L.

## Vorschläge für Grapheinträge
- „Das Fluchtmuster von sexp ist oben doppelt periodisch (Gitter ℤ + 2πi/L·ℤ)“ — Satz — bewiesen (Gitter) / numerisch (Kneser-Nähe)
- „Fluchtzeit in Kneser-Koordinaten: N(z) = (Im L/Re L)·Im z − Re z + Λ-periodisch, N(z+ω)=N(z)+4“ — Satz — bewiesen für F_L, numerisch für sexp
- „Kneser bildet den Grundstreifen in eine Linse in |Im w| < 1.756 ab; Im F > 0 oben“ — Beobachtung — numerisch
- „Erste nichtreelle Urbilder der reellen Achse bei Re z ≈ 2.52“ — Beobachtung — numerisch
- „Kein Haar von eᶻ häuft sich bei L; Haar-Urbilder unter sexp steigen in senkrechten Streifen nicht ins Unendliche“ — Satz — hergeleitet mit Lücke (Schleicher–Zimmer)
- „F⁻¹(I(exp)) ist dicht, 1-invariant, vom Mass null“ — Satz — hergeleitet mit Lücke (Literatur)

## Offene Fragen
1. Hausdorff-Dimension von E in jeder offenen Menge = 2 (McMullen lokal übertragen) — Lücke: lokale Fassung von McMullens Satz.
2. Wie sieht E_L/Λ auf dem Torus ℂ/Λ aus — welche Homologieklassen tragen die Haar-Urbilder? (Nach 9: keine Klasse, die in Richtung L läuft.)
3. Gibt es eine geschlossene Beschreibung der Linse Ω (Kreise sind es nicht, Befund 11)?
4. Lage der Menge {T = ∞} (Devaneys unzerlegbares Kontinuum in Kneser-Koordinaten): ausser ℝ sichtbar nichts — ist F⁻¹(Λ₀) ∩ Streifen = ℝ-Stück?
5. Für andere Basen b > e^{1/e} dasselbe Gitter ℤ + 2πi/L_b·ℤ; wie ändert sich die Steigung Im L_b / Re L_b, und was passiert bei b → e^{1/e}?
