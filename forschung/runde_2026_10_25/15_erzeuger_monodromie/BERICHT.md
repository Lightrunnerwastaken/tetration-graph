# Agent 15 — Monodromie des regulären Erzeugers J_L von exp

Zeit: Start 13:22, Ende ≈ 13:55 (Uhr aus `date`). Basis e. Bauen auf Runde 18, Bericht 10 (A5, A5a, offene Frage 1).

## Kurzfassung (einfache Sprache)

Der reguläre Erzeuger J_L (das Geschwindigkeitsfeld der Tetration, die am Fixpunkt L „regulär“ ist) ist mehrdeutig: Läuft man um 0, 1, e, e^e, … herum,
kommt man auf einem anderen „Blatt“ an. Wir haben **alle Blätter und alle Umläufe vollständig beschrieben**. Jedes Blatt hat eine Adresse, eine endliche
Folge ganzer Zahlen (m₁, m₂, m₃, …) — das sind die Zweige der wiederholten Logarithmen. Ein Umlauf um 0 erhöht m₁ um 1. Ein Umlauf um 1 erhöht m₂ um 1,
aber nur, wenn m₁ = 0 ist; sonst passiert gar nichts. Ein Umlauf um e erhöht m₃, nur wenn m₁ = m₂ = 0 ist, und so weiter — wie ein Zählwerk, bei dem
das k-te Rad sich nur dreht, wenn alle Räder davor auf null stehen. Daraus folgt: Die Gruppe aller Umläufe ist ein aus lauter „Lampenanzünder-Gruppen“
geschachteltes Kranzprodukt (Z ≀ Z ≀ Z ≀ …), sie ist **nicht frei** (es gibt echte Beziehungen zwischen den Umläufen), und jede Verzweigung ist
logarithmisch. Für Knesers Tetration heisst das: Ihre Korrekturfunktion θ hat an **jeder ganzen Zahl einen logarithmischen Verzweigungspunkt**
und wird dort unbeschränkt — genau so viel, dass Knesers Erzeuger auf der reellen Achse glatt bleibt.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| M1 | Die Riemannfläche von σ (und von J_L) ist die ganze u-Ebene C mit Projektion Ψ; Ψ: C∖Ψ⁻¹(P) → C∖P ist eine Überlagerung; die Blätter von J_L entsprechen eineindeutig denen von σ | **bewiesen** | Runde 18 A5 (Wegehebung) + Forster Satz 4.19; Ψ(cu) ≡ Ψ(u) ⇒ c = 1 | Überlagerungsaussage für Poincaré-Funktionen bekannt; für J_L vermutlich neu formuliert |
| M2 | Kodierung der Faser über L durch endliche Folgen (m₁, m₂, …) ∈ Z^(∞); Umläufe: a (um 0): m₁ ↦ m₁+1; b_k (um p_k = exp^k(0)): m_{k+1} ↦ m_{k+1}+1, falls m₁ = … = m_k = 0, sonst identisch | **bewiesen** + **numerisch** | Kranzrekursion a = Verschiebung (alle Abschnitte trivial), b_k = (b_{k−1} an der Stelle 0); 48/48 Proben (12 Blätter × a, b₁, b₂, b₃) in Gleitkomma, Gegenprobe mpmath dps 30: 12/12 (4 Blätter × a, b₁, b₂) | vermutlich neu |
| M3 | Auf dem Blatt m ist p_k genau dann Verzweigungspunkt, wenn m₁ = … = m_k = 0; dann logarithmisch; sonst sind σ, J_L dort regulär. 0 ist auf jedem Blatt logarithmisch. Auf dem Hauptblatt ist **jedes** exp^k(0) ein logarithmischer Verzweigungspunkt (offene Frage 1 aus Runde 18) | **bewiesen** | aus M2: Zyklen von b_k unendlich bzw. Fixpunkte | neu (Antwort auf R18-Frage 1) |
| M4 | Monodromiegruppe G = ⋃ G_k, G_k = ⟨a, b₁, …, b_k⟩ ≅ G_{k−1} ≀ Z (eingeschränktes Kranzprodukt), G₀ = Z, G₁ = Z ≀ Z (Lampenanzünder); G ist nicht frei, nicht endlich erzeugt, lokal auflösbar, nicht auflösbar; G_k hat Ableitungslänge k+1 | **bewiesen** (Ableitungslänge: Standardtatsache über Kranzprodukte) | Relationen [b₁, aⁿb₁a⁻ⁿ] = 1 (n ≠ 0), [b₂, b₁b₂b₁⁻¹] = 1, [b₂, ab₁a⁻¹] = 1 | neu |
| M4a | Die Relationen von M4 gelten geometrisch (Hebung des Kommutatorweges schliesst sich), [a,b₁], [b₁,b₂], [a,b₂] nicht | **numerisch** | 8 Probeblätter: Sprung < 4·10⁻¹⁴ bzw. bis 7.5 | — |
| M4b | G ist zugleich die iterierte Monodromiegruppe IMG(exp) (Basispunkt L, Rand C∖P) | **bewiesen** | Wirkung auf endliche Wörter und auf Z^(∞) hat denselben Kern | IMG-Sprache (Nekrashevych) bekannt; Ergebnis für exp vermutlich neu |
| M5 | Geschlossene Werte: nach b_k^j (Hauptblatt, Basispunkt L) ist σ = L^{k+1}·s_j und J_L = L^{k+1}·J^{(j)}, s_j, J^{(j)} = Werte längs ℓ_j bei L+2πij; s₁ = −0.95262397349458069216579816424 + 2.24193859082677437839204171742 i, J^{(1)} = −14.3223459699012846576275 + 10.7051679581953509008574 i | **bewiesen** (Formel) / **numerisch** (Werte, 30 St.) | J^{(1)} stimmt mit Runde 18 (Iterierte-Log-Methode) auf alle 20 dort angegebenen Stellen; zwei Wege ℓ₁ (senkrecht / links herum) gleich auf 2·10⁻³⁰ | neu |
| M6 | Kneser: J_K = J_L · Q(σ^{2πi/L}) mit Q(q) = 1/Φ′(s), q = e^{2πis}, Φ = (id+θ)⁻¹ — Knesers Erzeuger ist eine Funktion von σ allein (über die Potenz σ^{2πi/L}) | **bewiesen** (auf dem Definitionsgebiet von Φ) | aus J_K = J_L(1+θ′∘slog_K) (R18 A4) und Φ(s+1) = Φ(s)+1 | vermutlich neu |
| M7 | Knesers θ hat an jeder ganzen Zahl n einen logarithmischen Verzweigungspunkt; s + θ(s) → ∞ für s → n von oben | **bewiesen** | F_L lässt 0 aus, sexp_K(−1) = 0, a wirkt frei (M3); Periodizität | Singularität bei ganzen Zahlen wohl bekannt; Art (logarithmisch, unbeschränkt) und Begründung über die Monodromie vermutlich neu |
| M7a | θ(iy) wandert für y = 0.3 → 0.004 extrem langsam: Re(L·θ(iy)) = 1.098, 1.166, 1.264, 1.347, 1.382, 1.407 (monoton); \|F_L(−1+iy+θ(iy))\|/y → 1.09176 = sexp′(0) | **numerisch** | Reihen mit 2048 / 4096 Gliedern, 12 St. gleich (für y ≥ 0.004 ist die Gegenprobe schwach, da Glieder über 2048 < e⁻⁵¹) | — |
| M8 | Riemannfläche der regulären slog h = (log σ)/L ist die s-Ebene mit Projektion F_L; Überlagerung über C∖(P∪{L}); Faser Z^(∞) × Z; Monodromiegruppe enthält ⊕_{Z^(∞)} Z (Normalenhülle des Umlaufs um L) mit Quotient G | **bewiesen** (Überlagerung, Untergruppe, Quotient) / **hergeleitet mit Lücke** (ob die Erweiterung zerfällt, also ≅ Z ≀_{Z^(∞)} G, nicht geprüft) | — | vermutlich neu |
| M10 | Für jede Basis b > η gilt dieselbe Rekursion (a Verschiebung, b_k = (b_{k−1} an 0)), also dieselbe Gruppe G | **hergeleitet mit Lücke** | Rechteck- und Halbebenenargument übertragen sich wörtlich (Im log_b w ∈ (0, π/ln b) für Im w > 0, 0 < Im L_b < π/ln b); Lücke: Wegehebung (R18 A5) nur für b = e ausgeschrieben, keine Numerik für b ≠ e | vermutlich neu |
| M9 | Die Kodierung hängt von den Verbindungswegen ab: mit senkrechten Wegen ℓ_m (rechts an 0 vorbei) hätte a an der Stelle −1 den Abschnitt a⁻¹; erst die Wege links an 0 vorbei machen a zur reinen Verschiebung | **bewiesen** (Rechteckargument) | — | Hinweis |

## Beweise und Einzelheiten

**Bezeichnungen.** f = exp, L = 0.31813…+1.33723…i (e^L = L, Multiplikator λ = L, log λ = L), Ψ ganz mit Ψ(Lu) = e^{Ψ(u)}, Ψ(0) = L, Ψ′(0) = 1,
σ = Ψ⁻¹ (Keim bei L), J_L = Lσ/σ′, P = {p₀ = 0, p₁ = 1, p₂ = e, p₃ = e^e, …}. Basispunkt L.

### M1 — Riemannfläche
Ψ hat keine kritischen Punkte (R18), ist also lokal biholomorph. Runde 18 A5 zeigt: der Hauptkeim von σ lässt sich längs jedes Weges in C∖P fortsetzen.
Dasselbe gilt für jeden Keim: Ist Ψ(u₀) = z₀, so ist u₀/Lⁿ für grosses n in der Univalenzscheibe von Ψ um 0, also ist der Keim z ↦ Lⁿσ_loc(ℓ_n(z)) mit
ℓ_n = Ψ(·/Lⁿ)∘σ eine Verkettung von Log-Zweigen — genau die Form, für die A5 die Fortsetzung beweist. Also hat Ψ: Ψ⁻¹(C∖P) → C∖P die Wegehebungs-
eigenschaft und ist nach Forster (Riemannsche Flächen, Satz 4.19) eine Überlagerung. Die maximale Fortsetzung von σ ist also die u-Ebene selbst
(Punkte über P eingeschlossen, wo Ψ sie annimmt); die Faser über z ist Ψ⁻¹(z), abzählbar; die Monodromie ist transitiv (C∖Ψ⁻¹(P) zusammenhängend).
J_L auf dem Blatt u ist K(u) = L u Ψ′(u). Zwei Blätter u, u′ mit gleichem J-Keim hätten σ₂/σ₂′ = σ₁/σ₁′, also σ₂ = cσ₁, also Ψ(cu) = Ψ(u) auf einer
offenen Menge, also überall; Ableiten bei 0 gibt c = 1. Damit hat J_L dieselbe Monodromie wie σ. ∎

### M2 — Die Kranzrekursion
Erzeuger von π₁(C∖P, L): a = Lasso längs der Strecke [L, 0], b_k = Lasso längs [L, p_k] (Strecken von L zur reellen Achse, kreuzungsfrei).
Urbilder von L unter exp: L + 2πim. Verbindungswege ℓ_m: L → −1 + i·Im L → −1 + i(Im L + 2πm) → L + 2πim (links an 0 vorbei).
Abschnitt g|_m = ℓ_m · (Hebung von g ab L+2πim) · ℓ_{g(m)}⁻¹.
- **a:** Die Hebung ab L+2πim geht waagrecht nach links (Höhe Im L + 2πm), um 2π hinauf, zurück nach L+2πi(m+1); also a(m) = m+1.
  Weil ℓ_m genau die zusammengesetzten Hebungen von a sind (bis auf Homotopie, die Rechtecke enthalten keinen Punkt von P), ist jeder Abschnitt trivial.
- **b_k (k ≥ 1):** b_k umläuft 0 nicht, die Hebung ist geschlossen: b_k(m) = m. Für m ≠ 0 umläuft sie p_{k−1} + 2πim ∉ P, Abschnitt trivial.
  Für m = 0 ist sie das Lasso längs log([L, p_k]) um p_{k−1}; das Gebiet zwischen log([L,p_k]) und [L, p_{k−1}] liegt in der oberen Halbebene und berührt
  die reelle Achse nur in p_{k−1}, also b_k|₀ = b_{k−1} (mit b₀ = a).

Kodierung: code(u) = (m₁, code(u′)) mit Ψ(u/L) = L + 2πim₁ und u′ = Hebung von ℓ_{m₁}⁻¹ ab u/L. Für die Monodromie M_g gilt dann
code(M_g u) = g(m₁)·code(M_{g|m₁} u′) — die Wirkung auf den Adressen ist genau die Kranzrekursion. Nur u = 0 hat Adresse 0^∞ (Ψ ist nahe 0 injektiv);
da G transitiv auf der Faser wirkt und die endlichen Folgen eine G-Bahn bilden, sind alle Adressen endlich, und code ist eine G-äquivariante Bijektion
Faser → Z^(∞). Ausgerechnet: a(m₁, m₂, …) = (m₁+1, m₂, …); b_k ändert nur m_{k+1}, und nur wenn m₁ = … = m_k = 0. ∎

*Numerik (m1).* Hebung durch Newton-Verfolgung von Ψ(u(t)) = z(t) (unabhängig von der Iterierte-Log-Methode aus R18). Für 12 Blätter
((), (1), (−1), (2), (0,1), (0,−1), (3,−2), (0,0,1), (1,0,1), (0,1,−1), (−1,0,0,1), (0,0,0,2)) werden u aus der Adresse gebaut, zurückkodiert (12/12 gleich),
dann a, b₁, b₂, b₃ hochgehoben und neu kodiert: **48/48 wie vorhergesagt** (Gleitkomma, `m1_ausgabe_A.txt`). Gegenprobe mpmath dps 30, feinere Schritte
(`m1_ausgabe_B.txt`, `m1_ausgabe_B2.txt`; wegen des 5-Minuten-Limits nur die Blätter (), (1), (0,−1), (0,0,1) mit a, b₁, b₂; der erste Lauf wurde nach
2 Blättern vom Zeitlimit beendet, der zweite lief durch): **12/12 gleiche Adressen**, Sprünge bei trivialer Wirkung < 10⁻²⁹.

### M3 — Art der Verzweigungen
Über einer kleinen punktierten Scheibe D* um p_k zerfällt die Überlagerung in Komponenten; die Komponente des Blattes m ist endlich-zyklisch vom Grad =
Zykellänge von b_k durch m, oder universell bei unendlichem Zyklus. Nach M2 ist der Zyklus unendlich (m_{k+1} läuft durch Z), falls m₁ = … = m_k = 0,
sonst Länge 1. Grad 1: die Komponente ist eine punktierte Scheibe, deren Loch ein Punkt von C ist (∞ ist ausgeschlossen: sonst enthielte die Komponente eine punktierte Umgebung von ∞, auf der Ψ beschränkt wäre — Liouville),
also ist σ dort regulär mit Ψ(σ(p_k)) = p_k. Unendlich: logarithmische Singularität (asymptotischer Kanal, σ → ∞). Endliche Zyklen der Länge ≥ 2 kommen nie vor:
**es gibt keine algebraischen Verzweigungspunkte.** ∎

### M4 — Gruppenstruktur
b_j (j ≥ 1) = δ₀(b_{j−1}) (wirkt nur unter dem Buchstaben 0), a = Verschiebung. Konjugation mit aⁿ setzt δ₀ nach δ_n um. Also enthält G_k die Gruppe
⊕_n δ_n(G_{k−1}) und die Verschiebung, d. h. G_k = G_{k−1} ≀ Z; umgekehrt liegt jeder Erzeuger darin. G₁ = Z ≀ Z wirkt auf (m₁, m₂) als
(m₁, m₂) ↦ (m₁ + c, m₂ + f(m₁)), f endlich getragen. Jede endlich erzeugte Untergruppe liegt in einem G_k, das nur die ersten k+1 Buchstaben ändert,
b_{k+1} aber den (k+2)-ten: G ist nicht endlich erzeugt, lokal auflösbar, wegen unbeschränkter Ableitungslänge nicht auflösbar. Relationen:
aⁿb₁a⁻ⁿ wirkt nur bei m₁ = n, b₁ nur bei m₁ = 0 → sie vertauschen; in der freien Gruppe π₁(C∖P) ist [b₁, ab₁a⁻¹] ≠ 1. Also ist die Monodromie
**nicht treu**: π₁(C∖P) → G hat einen grossen Kern. ∎ IMG: Ein Element wirkt auf allen endlichen Wörtern trivial genau dann, wenn auf allen Wörtern v0^∞, also stimmt G mit IMG(exp) überein.

### M5 — Geschlossene Werte
Adresse 0^k j 0^∞ gehört zu u = L^{k+1}s_j, s_j = Endpunkt der Hebung von ℓ_j ab 0 (wegen ℓ₀ trivial). Aus Ψ(L s_j) = e^{L+2πij} = L und
Ψ′(Lv)L = e^{Ψ(v)}Ψ′(v) folgt Ψ′(L^i s_j) = Ψ′(s_j) für alle i ≥ 1, also J_L = L·L^{k+1}s_j·Ψ′(s_j) = L^{k+1}J^{(j)}.
Probe: |Δu| nach a, b₁, b₂, b₃ vom Hauptblatt = |L|^{1,2,3,4}·|s₁| = 3.348, 4.602, 6.326, 8.696 (gemessen genau so).

### M6/M7 — Verbindung zu θ
Mit h = slog_L = (log σ)/L, slog_K = Φ∘h, Φ = (id+θ)⁻¹, ist 1 + θ′(slog_K) = 1/Φ′(h); Φ′ ist 1-periodisch, also Φ′(s) = R(e^{2πis}) und
**J_K = J_L · Q(e^{2πih}) = J_L · Q(σ^{2πi/L}),** Q = 1/R (M6).
*M7, Beweis.* Für Im s > 0 gilt sexp_K(s) = F_L(s + θ(s)), F_L(w) = Ψ(e^{Lw}). Ψ lässt 0 aus (Ψ(u) = 0 ⇒ e^{Ψ(u/L)} = 0), also auch F_L. Für s → −1 von oben
geht sexp_K(s) → 0; bliebe s + θ(s) auf einer Folge beschränkt, hätte F_L am Häufungspunkt den Wert 0 — also |s + θ(s)| → ∞. θ lässt sich in eine punktierte
Umgebung von −1 fortsetzen (θ = h∘sexp_K − id, sexp_K ist dort analytisch mit einfacher Nullstelle, h längs jedes Weges in C∖P fortsetzbar). Ein Umlauf um −1
wird von sexp_K auf einen Umlauf um 0 abgebildet; dort wirkt a auf **jedem** Blatt mit unendlichem Zyklus (M2/M3). Also kehrt θ nie zurück:
logarithmischer Verzweigungspunkt. Wegen θ(s+1) = θ(s) gilt dasselbe an jeder ganzen Zahl. ∎
Bedeutung: J_K ist auf der reellen Achse analytisch (auch bei 0, 1, e, …), J_L auf dem Kneser-Blatt nicht; der Faktor (1 + θ′∘slog_K) hebt die Monodromie
von J_L an jedem p_k genau auf. **Die Singularitäten von θ an den ganzen Zahlen sind die Spur der Monodromie von J_L.**
Numerik (m3): die Divergenz ist extrem langsam (Grössenordnung slog(log(1/y))); Re(Lθ(iy)), also log|σ|, steigt monoton.

### M8 — Riemannfläche der regulären slog
F_L = Ψ∘exp_L hat die singulären Werte P ∪ {L} (L als asymptotischer Wert, Re(Ls) → −∞); F_L ist daher eine Überlagerung über C∖(P∪{L}) und die s-Ebene ist die
Riemannfläche von h. Faser = {(u, k)}: u ∈ Ψ⁻¹(z), k Zweig von log u. Der Umlauf δ um L ändert k nur auf dem Blatt, dessen u nahe 0 liegt; seine Konjugierten
gδg⁻¹ ändern k auf genau einem anderen Blatt; da G transitiv ist, entsteht ⊕_{Z^(∞)} Z, und der Quotient ist G. Lücke: Zerfall der Erweiterung nicht geprüft.

## Vorschläge für Grapheinträge
- Riemannfläche des regulären Erzeugers ist die u-Ebene (Überlagerung über C∖P) — satz — bewiesen
- Monodromie des regulären Erzeugers: a verschiebt m₁, b_k verschiebt m_{k+1} falls m₁…m_k = 0 — satz — bewiesen (48/48 numerisch bestätigt)
- Jeder Punkt exp^k(0) ist auf dem Hauptblatt logarithmischer Verzweigungspunkt; keine algebraischen Verzweigungen — folgerung — bewiesen
- Monodromiegruppe = IMG(exp) = iteriertes Kranzprodukt Z ≀ Z ≀ …, nicht frei, nicht endlich erzeugt — satz — bewiesen
- Werte nach einem Umlauf: J_L = L^{k+1}·J^{(1)}, J^{(1)} = −14.32234596990128466 + 10.70516795819535090 i — beobachtung — bewiesen (Formel) / numerisch (Wert)
- Knesers Erzeuger als Funktion von σ^{2πi/L} — satz — bewiesen
- Knesers θ: logarithmische Verzweigung und Unbeschränktheit an jeder ganzen Zahl — satz — bewiesen
- Riemannfläche der regulären slog: s-Ebene, Monodromie ⊕Z-Erweiterung von G — satz — bewiesen bis auf Zerfall

## Offene Fragen
1. Zerfällt die Erweiterung in M8, d. h. ist die Monodromie der regulären slog das Kranzprodukt Z ≀_{Z^(∞)} G?
2. Riemannfläche von slog_K selbst: Φ ist nur auf einem Halbstreifen-Gebiet erklärt. Wie sieht die maximale Fortsetzung von slog_K aus, und welche
   Blätter von h (Adressen) erreicht sie? Vermutung: genau die Blätter, deren Adresse an den Stellen gemäss der θ-Verzweigung von 0 abweicht.
3. Asymptotik von θ an der ganzen Zahl: gilt θ(n+iy) = slog_L(log(1/y)) + O(1) (passender Zweig)? Numerik reicht bis y = 0.004 nur für den Trend.
4. Basis-Abhängigkeit (siehe M10): Für b < η (die Bahn von 0 konvergiert gegen a) ist P nicht mehr abgeschlossen-diskret, es häuft sich bei a —
   dort ändert sich die Monodromie grundlegend. Nicht gerechnet.
5. Gibt es eine Wahl von θ (bzw. einer Abelfunktion h + p(h)), deren Erzeuger auf einer grösseren Teilfläche einwertig ist? M6 sagt: p muss die Monodromie
   von J_L an den p_k aufheben; die Knesersche Lösung tut das auf der reellen Achse.

## Reproduktion
Python `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`, kein Motor.
- `m1_monodromie.py A` → `m1_ausgabe_A.txt` (≈ 25 s; M2, M5); `m1_monodromie.py B` / `B2` → `m1_ausgabe_B.txt`, `m1_ausgabe_B2.txt` (mpmath, je ≤ 5 min, Teilmengen)
- `m2_relationen.py` → `m2_ausgabe.txt` (M4a, ≈ 10 s)
- `m3_theta_ganzzahl.py` → `m3_ausgabe.txt` (M7a; liest θ-Koeffizienten aus `identitaeten_2026_10_14/etappe2/daten/` über `runde_2026_10_18/10_neues_gebiet/lib10.py`, nur lesend)
