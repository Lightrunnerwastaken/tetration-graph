# Runde 18.10.2026 — 03 Basis-Ebene: Tetration als Funktion der Basis

## Kurz in einfacher Sprache

Unter der Schwelle η = e^{1/e} baut man die Tetration „regulär“ (am anziehenden Fixpunkt), über η nach Kneser. Das sind zwei völlig
verschiedene Bauweisen. Wir haben gemessen: **Kneser über η ist (auf 33 Stellen) einfach die Fortsetzung der regulären Tetration über η hinaus.**
Man nimmt die Taylor-Reihe der regulären Tetration in der Basis, rechnet sie über η hinaus aus — und trifft Knesers Werte.
Die gemeinsame Reihe bei η konvergiert bis zur Basis 1. Trotzdem ist η kein ganz gewöhnlicher Punkt: die reguläre Tetration hat bei η eine
winzige, extrem flache Singularität (Grösse etwa 10⁻³⁰ bis 10⁻⁵⁰), die man nur in sehr hohen Taylor-Koeffizienten sieht.
Die erste Ableitung nach der Basis bei η haben wir als geschlossene Summe (Nörlund-Hauptlösung) hergeleitet und auf 13 Stellen bestätigt.

Wichtiges Hilfsresultat: Der Motor rechnet für komplexe Basen im Shell–Thron-Inneren die reguläre Tetration (geprüft gegen eine eigene,
unabhängige Koenigs-Rechnung auf 15 Stellen); Kneser nur im dünnen „Horn“ rechts von η, das die Shell–Thron-Spitze bildet.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | Nahe η ist das Kneser-Gebiet ein dünnes Horn: Shell–Thron-Rand bei $\lambda=e^{it}$ gibt $b-\eta=\frac{\eta}{e}\bigl(\frac{t^2}{2}+\frac{it^3}{6}\bigr)+O(t^4)$, also $\lvert\mathrm{Im}\,b\rvert\approx\frac{\sqrt2\,e^{1/2}}{3\sqrt\eta}\,(\mathrm{Re}\,b-\eta)^{3/2}$; Spitze bei η, alle übrigen Basen um η sind im Shell–Thron-Inneren | bewiesen (Taylor von $\lambda e^{-\lambda}$) | Rechnung unten | Spitze bekannt (Kardioide in λ), die 3/2-Breite vermutlich nicht eingetragen |
| B2 | Der Motor (fatou_fork) liefert für komplexe Basen im Shell–Thron-Inneren die reguläre (Koenigs-)Tetration am anziehenden Fixpunkt | numerisch | b = 1.4±0.01i, 1.4+0.1i: Motor = eigene Koenigs-Rechnung (`regulaer.py`) auf 15 Stellen, auch bei z = 0.5+3i | Motorverhalten, für den Graphen neu |
| B3 | **Kneser = Fortsetzung der regulären Familie durch η.** Taylor-Reihe von $b\mapsto T^{\rm reg}_b(x)$ um $b_0=\eta-s$, ausgewertet bei $\eta+\delta$, trifft Kneser | numerisch | x=0.5: s=0.05 → Abweichung 5·10⁻³⁰ (δ=0.01), 3·10⁻²⁷ (δ=0.02). Mit Reihe bei η (H₂₀, s=0.03): $K-H$ = 1·10⁻³³ (δ=10⁻³ und 10⁻² nach Abzug des a₀-Artefakts), unter dem Restglied für δ=0.05/0.1/0.2. Gegenprobe x=−0.5: 1.8·10⁻³³. Kneser bei dps 40 und 60 gleich | vermutlich neu (Graph hat nur gleiche Steigung, kl.beob.stetig_eta) |
| B4 | Die Taylor-Koeffizienten $a_n$ von $T_b(1/2)$ bei η (Tabelle unten) sind von beiden Seiten gleich und wechseln das Vorzeichen; Konvergenzradius $\eta-1$, Singularität bei $b=1$, $a_n\sim C(-1)^n n^{-\gamma}(\eta-1)^{-n}$ mit $\gamma\approx2.4$ (erwartet 5/2 aus $T_b(1/2)\approx 1+\varepsilon-\varepsilon^{3/2}$, $\varepsilon=b-1$) | numerisch | a₀…a₁₂ aus s=0.02 und s=0.03 auf ≥13 Stellen gleich; a₀ = sexpeta(1/2) auf 10⁻³⁰ | neu |
| B5 | Die reguläre Familie ist bei η **nicht** analytisch: in der Multiplikator-Koordinate hat $g(\lambda)=T^{\rm reg}(1/2)$ um λ₀=0.8 neben dem Hauptteil (Radius ≈0.8, Singularität λ=0 ↔ b=1) einen winzigen Anteil mit Radius genau $1-\lambda_0$ = 0.200 (k=64…248, dps 150); Amplitude ~10⁻³⁶, um λ₀=0.9 unter 10⁻⁶⁵ — flache (wesentliche) Singularität bei λ=1, also b=η | numerisch | e2b: Radius 0.20±0.02 über k=64…248 (Rauschgrenze 10⁻¹⁵⁰); b-Koeffizienten um η−s fallen unter-geometrisch (≈e^{−c√k}) | neu |
| B6 | Stärke dieser Singularität: bei fester Ordnung k=64 wächst $\ln|c_k s^k|$ etwa linear in $-1/\sqrt s$ mit Steigung ≈ 10–15 | numerisch, grob | s = 0.02, 0.03, 0.05, 0.08 (e7) | — |
| B7 | Folgerung aus B3–B5: $b\mapsto T_b(x)$ (regulär unter, Kneser über η) ist bei η numerisch $C^\infty$ mit konvergenter gemeinsamer Taylor-Reihe $H$, aber nicht analytisch: $T^{\rm reg}=H+E^-$ mit flachem, bei η singulärem $E^-$, $T^{\rm K}=H$ (+ höchstens $10^{-33}$) | Vermutung | B3–B5 | beantwortet kl.frage.eta_analytisch numerisch: „glatt, nicht analytisch“, aber anders als vermutet sitzt der Defekt auf der regulären Seite |
| B8 | **Erste Basisableitung bei η geschlossen:** $c_1(x)=\partial_bT_b(x)|_\eta=T_\eta'(x)\,G(x)$, $G(x+1)-G(x)=\frac{e}{\eta}\frac{T_\eta(x)}{T_\eta'(x)}$, $G(0)=0$, $G$ = Nörlund-Hauptlösung | hergeleitet mit Lücke | Formel = Taylor-Wert: x=0.5: 0.438581800170784 vs 0.4385818001708944 (1·10⁻¹³, N=60), x=−0.5: 2.5·10⁻¹³; N=20: 2·10⁻¹¹ (Konvergenz in N) | neu |
| B9 | Ränder bei λ = e^{iπ/2}, e^{iπ/4}, e^{iπ/8} sind „starke“ Singularitäten (Schwanzamplitude ~10⁻³); zu λ=1 hin bricht die Amplitude ein (θ=π/64: ~10⁻³⁰) | numerisch | e6, je zwei Kreisradien (gleiche Koeffizienten) | — |
| B10 | Kneser-Grössen für grosse Basen (Tabelle unten): $T_b'(0)$ wächst langsam (2.33, 2.79, 3.14 bei 10², 10³, 10⁴), Wendepunkt $x_w(b)$ wandert nach links (−1.18, −1.25, −1.30) | numerisch (Werte) | x_w(2) = −0.0302445427536393 trifft den zertifizierten Wert | kein Gesetz gefunden |

## Details

### B1 Horn
Mit $b=\exp(\lambda e^{-\lambda})$ (fp.satz.multiplikator_rezept) und $\lambda=e^{it}$: $\lambda e^{-\lambda}=e^{-1}(1+u)e^{-u}$, $u=it-t^2/2-it^3/6+\dots$,
$(1+u)e^{-u}=1-u^2/2+u^3/3+O(u^4)=1+t^2/2+it^3/6+O(t^4)$. Also $\ln b-1/e=(t^2/2+it^3/6)/e$, Spitze nach links mit Öffnung $\propto(\mathrm{Re}\,b-\eta)^{3/2}$.
Bei δ = 0.01 ist das Horn nur ±6.5·10⁻⁴ breit (numerisch bestätigt: die reguläre Tetration wird bei Im b ≈ 7·10⁻⁴ wild, |λ| → 1).

### B3/B5 Methode
`regulaer.py`: $g(\lambda,x)=\sigma^{-1}(\lambda^x\sigma(1))$ über Rückwärts-Logarithmen ab $a+\lambda^{x+n}\sigma(1)$ mit Schröder-Reihe 3. Ordnung.
Warum g im punktierten Einheitskreis analytisch ist: für |λ|<1 konvergiert der Turm gegen a (Shell–Thron), σ ist dort analytisch in λ, Ψ=σ⁻¹ ist ganz.
Deshalb kann eine Singularität im Abstand 0.200 von 0.8 nur auf dem Einheitskreis nahe 1 liegen (der nächste Kreispunkt ist 1, Abstand 0.2).
Ob sie genau bei λ=1 oder bei nahen parabolischen Punkten $e^{2\pi ip/q}$ sitzt, trennt der Wurzeltest nicht (R = 0.200…0.206).
Das a₀-Artefakt: Umentwicklung der Reihe um η−s an der Stelle η summiert den Singulärschwanz mit (−1.5·10⁻³⁰ bei s=0.03, gegen sexpeta geprüft) — das ist selbst ein Zeichen der Singularität.

### Tabelle a_n (x = 1/2), $T_b(1/2)=\sum a_n(b-\eta)^n$
| n | a_n | n | a_n |
|---|---|---|---|
| 0 | 1.257153075054172622172 | 7 | 0.30199804837399430 |
| 1 | 0.438581800170894429326 | 8 | −0.48038294449992811 |
| 2 | −0.184336953718053020244 | 9 | 0.79867809308437451 |
| 3 | 0.129898167501740064395 | 10 | −1.3730210119323530 |
| 4 | −0.124702011169116832761 | 11 | 2.4235324371497700 |
| 5 | 0.147965640516363063694 | 12 | −4.370822084408 |
| 6 | −0.201950074503851512877 | 15 | 28.258895 |
x = −1/2: a₀…a₃ = 0.62207798171892641606, −0.22217528150401936021, 0.25915103947337897874, −0.36084872568475470149.

### B8 Herleitung (Lücke benannt)
Ableiten von $T_b(x+1)=b^{T_b(x)}$ nach b bei η: $c(x+1)=T(x+1)\,[c(x)/e+T(x)/\eta]$, mit $T'(x+1)=T(x+1)T'(x)/e$ und $c=T'G$ folgt
$G(x+1)-G(x)=r(x)=\frac{e}{\eta}\frac{T(x)}{T'(x)}$; $T_b(0)=1$ gibt $G(0)=0$. Bekannt ist auch $G=-\partial_b\,\mathrm{slog}_b\circ T_\eta$.
Die Lösungen unterscheiden sich um 1-periodische Funktionen; $r(x)\sim e x^2/(2\eta)$ wächst polynomial.
**Lücke:** dass gerade die Nörlund-Hauptlösung (asymptotische Entwicklung für x→∞ ohne periodischen Anteil) die Basisableitung ist; plausibel,
weil die Fatou-Koordinate der Störung eine asymptotische Entwicklung in $1/(y-e)$ und Logarithmen hat. Formel mit Euler–Maclaurin (Glieder bis r‴) in `e8_c1_noerlund.py`.

### Deutung (Vermutung)
Kneser für b>η ist die Douady–Lavaurs-Torkoordinate (Zylinder zwischen L und L̄); die Grenzwertsätze der parabolischen Implosion geben $T^{\rm K}_b\to T_\eta$
(untere η-Tetration). B3 sagt mehr: Kneser hängt so glatt von b ab, dass er die konvergente Reihe H fortsetzt; die reguläre Familie trägt dagegen
die Spur der dicht liegenden parabolischen Randpunkte (B9), die zu η hin flach wird. Ein Kandidat für die Grösse von $E^-$ ist $\exp(-2\pi N)$ mit der
Durchgangszeit $N=\pi\sqrt{2\eta/(e\,|b-\eta|)}$ (gibt Steigung ≈ 20 in B6 statt 10–15 bei effektivem Abstand s+ρ; nur grob).

### B10 Werte (dps 30)
| b | T'(0) | T''(0) | x_w |
|---|---|---|---|
| 2 | 0.889364954620976 | 0.0173530979307 | −0.0302445427536393 |
| 3 | 1.14909879145946 | 0.716664350009 | −0.603863714358098 |
| 5 | 1.40602362406508 | 1.62222718161 | −0.84692037022525 |
| 30 | 2.03120317345009 | 4.63385039301 | −1.11289628131011 |
| 100 | 2.33385048951909 | 6.47678151514 | −1.18188034423232 |
| 1000 | 2.78599024341557 | 9.67356140182 | −1.25479311682541 |
| 10⁴ | 3.1405112564132 | 12.5406025267 | −1.29677936033631 |
b = 10⁶ scheitert im Motor (invz0rfunc), b = 10 scheiterte an einer Dateisperre (nicht wiederholt).

### Fehlschläge / nicht gemacht
- Kneser-Werte unter b = η + 10⁻⁴ scheitern im Motor bei dps 30 (exp-Überlauf).
- Die Konstante in B6 ist nicht sauber bestimmt; Monodromie um η nicht gerechnet (Motor wechselt beim Kreis um η automatisch zwischen Kneser und regulär).
- Teil (2) des Auftrags (Knospengrössen, Misiurewicz-Punkte) nicht bearbeitet; θ-Welle a₀(b), κ_b nicht gerechnet.

## Reproduktion (Python: miniconda, `PYTHONIOENCODING=utf-8`)
- B2: `regulaer.py` gegen Motor (`motor.py`, eigener Speicher `daten/motor_cache.json`).
- B3: `e3_fortsetzung.py 50 128 0.05 0.8 0.5 1/100 1/50`; `e5_kneser_gegen_H.py e4_x0.5_s0.03_80_256_60.json 0.5 20 1/1000 1/100 1/50 1/20 1/10 1/5` (KDPS=40/60).
- B4: `e4_ableitungen_eta.py 80 256 0.03 0.8 0.5 60 30` und `... 0.02 ... 50 22`; x=−0.5 analog.
- B5: `e2b_schwanz.py 150 512 0.8 0.15 0.5`, `e2b_schwanz.py 100 384 0.9 0.075 0.5`; B6: `e7_schwanz_b.py 80 256 0.8 0.5 0.02 0.03 0.05 0.08`.
- B8: `e8_c1_noerlund.py 60 0.5 -0.5 1.5`; B9: `e6_randstaerke.py 60 256 0.8 0.15 0.5 0.5 0.25 0.125` (und 0.0625 0.03125 0.015625); B10: `e9_grosse_basis.py 30 …`.
- E1 (`e1_eta_grenze.py`): (T_b − T_η)/δ bei x = 0.5, −0.5, −1.5, 2 für δ = 10⁻², 10⁻³.

## Vorschläge für Grapheinträge
- kl.beob.kneser_fortsetzung_regulaer — „Kneser über η ist die Fortsetzung der regulären Tetration durch η (33 Stellen)“ — beobachtung — numerisch
- kl.beob.eta_taylor — „Taylor-Reihe von T_b(1/2) bei η: Koeffizienten, Radius η−1, Singularität bei b=1“ — beobachtung — numerisch
- kl.beob.regulaer_flach_singulaer — „Die reguläre Familie hat bei η eine flache Singularität“ — beobachtung — numerisch
- kl.verm.eta_glatt_nicht_analytisch — „T_b(x) ist bei η C^∞, nicht analytisch; Defekt auf der regulären Seite“ — vermutung — vermutung
- kl.satz.c1_noerlund — „Erste Basisableitung bei η als Nörlund-Hauptlösung“ — satz — skizze (Lücke: Wahl der Hauptlösung)
- fp.satz.horn_eta — „Das Kneser-Gebiet ist bei η ein Horn der Breite ∝ (b−η)^{3/2}“ — satz — bewiesen
- kn.beob.motor_komplex_regulaer — „Der Motor rechnet im Shell–Thron-Inneren regulär“ — beobachtung — numerisch

## Offene Fragen
1. Ist $T^{\rm K}=H$ exakt (Kneser analytisch in b auf einer Umgebung von η ohne Schlitz?) oder nur bis auf einen flachen Rest unter 10⁻³³?
2. Genaue Form von $E^-$: $\exp(-A/\sqrt{\eta-b})$? Konstante A; Zusammenhang mit der Durchgangszeit und den Écalle–Voronin-Invarianten von $e^u-1$.
3. Höhere Ableitungen $c_k$ ebenfalls als Hauptlösungen einer Hierarchie? (c₂ = 2a₂ = −0.36867390743610604049 bei x=1/2 als Prüfwert.)
4. Ist |λ|=1 eine natürliche Grenze der regulären Familie, mit Stärke → 0 zu λ=1 hin (B9)?
