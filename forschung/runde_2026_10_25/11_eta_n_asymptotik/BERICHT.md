# 11 — Asymptotik von ηᵣ, eᵣ und Kenngrößen der Berührung

Runde 25.10.2026, Arbeitszeit 13:23–13:50 (aus `date`). Ein Agent, keine Unteragenten.

**Kurz in einfacher Sprache.** ηᵣ ist die Basis, bei der die Rechenart Hᵣ die Diagonale gerade berührt. eᵣ ist die Höhe, bei der das passiert. Ein neues Glied η₉ war in 30 Minuten nicht zu rechnen. Ich habe deshalb die vorhandenen Glieder r=3…8 neu ausgewertet und an jeder Berührstelle zwei Formzahlen gemessen: die Krümmung βᵣ und das „iterative Residuum“ ιᵣ. ιᵣ ist die Zahl, die man nicht wegtransformieren kann. Sie sagt, wie langsam die nächste Rechenart an der Grenzbasis gegen eᵣ läuft. ιᵣ fällt von 1/3 (bei der Potenz) schnell auf etwa −0.6 und bleibt dann fast stehen. Für den Abstand 2−ηᵣ passt die Form C/(r·ln r) mit C≈2.20 besser als ein reines Potenzgesetz. Dann wäre der bisher gemessene Exponent 1.46 nur eine Übergangszahl. Für eᵣ passt keine einfache Form; das ist ein ehrlicher Fehlschlag.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| B1 | Iteratives Residuum ιᵣ=1−a₃/a₂² der Berührung (aₖ=Fᵣ⁽ᵏ⁾(eᵣ)/k!, Fᵣ=Hᵣ(ηᵣ,·)): ι₃=1/3 exakt; ι₄=−0.3940091010835869, ι₅=−0.5523418565090349, ι₆=−0.6006192688985848, ι₇=−0.6057589653632122, ι₈=−0.5928654176566894 | numerisch | zwei Cauchy-Einstellungen je Rang (ρ=0.1/M=32/dps 25 und ρ=0.03/M=16/dps 30; r=8: ρ=0.03 und 0.05, dps 20/24) stimmen auf ≥17 Stellen überein. Gilt für die gewählte Fortsetzung; bei r=8 die experimentelle Zwei-Bogen-Auswahl mit etwa 13 belegten Stellen | vermutlich neu |
| B2 | ιᵣ ist nicht monoton: Minimum bei r=7, danach steigt es (ι₈>ι₇). Zuwächse −0.727, −0.158, −0.048, −0.005, +0.013 | numerisch | aus B1 | vermutlich neu |
| B3 | Verfeinerte kritische Rate: eᵣ−Hᵣ₊₁(ηᵣ,k)=1/(βᵣ(k+ιᵣ ln k+Cᵣ+o(1))) für k→∞ | bewiesen unter den Voraussetzungen von M3 (reell-analytische Falte mit βᵣ>0, Konvergenz von unten) | Fatou-Koordinate, Beweis unten; r=3 gibt die bekannte Form e−ᵏη≈2e/k−(2e/3)·ln k/k² | Methode klassisch (Écalle/Fatou); Anwendung auf die Ränge vermutlich neu |
| B4 | Krümmung βᵣ=Fᵣ″(eᵣ)/2: 0.1839397 (1/(2e)), 0.2934833259466268, 0.3769856540413645, 0.4429294090481398, 0.4961476451837973, 0.5368242373506136 | numerisch | zwei Einstellungen wie B1; β₇, β₈ bestätigen die Werte der Runde 12.10. (0.4961476452, 0.5368242374) | Kontrolle; β₄–β₆ als Liste neu |
| B5 | βᵣ wächst mit fallenden Zuwächsen (0.110, 0.084, 0.066, 0.053, 0.041). Grenzwert offen: Drei-Punkt-Grenzfits geben 1.72 (r=5,6,7) und 0.87 (r=6,7,8) | numerisch, Grenzwert offen | instabile Fits = keine Aussage | — |
| B6 | Vermutung: 2−ηᵣ ~ C/(r ln r), C≈2.20. vᵣ=(2−ηᵣ)·r·ln r = 1.8303, 2.0222, 2.1136, 2.1596, 2.1827, 2.1935; Aitken-Werte 2.1966, 2.2063, 2.2061, 2.2028 | Vermutung | Gegenprobe mit anderen Normierungen: r^1.5 (Aitken 2.96, 3.03, 3.50, 2.33, instabil), (r−1)ln(r−1) und r(ln r)² (driften gleichsinnig). Nur 6 Glieder | vermutlich neu (Forum und Runde 12.10. nennen nur Potenzfits) |
| B7 | Folgerung aus B6: der lokale Exponent strebt langsam gegen 1, nicht gegen 1.46. Die gemessenen 1.462, 1.471, 1.470, 1.466, 1.460 sind dann eine Übergangszone | hergeleitet mit Lücke | Lücke = B6 selbst | — |
| B8 | Vorhersagen für η₉: 1.888827 (Modell B6, Zuwachsquotient 0.464), 1.889137 (Potenzfit r=3…6 der Runde 12.10.), 1.88892–1.88932 (A(r+c)^−p durch drei Punkte). Schon 4 Stellen von η₉ trennen die Modelle | Vermutung / Prüfvorschlag | `analyse.py` | — |
| B9 | Freie Drei-Punkt-Grenzfits L−A r^−q: für η geben sie L=1.99978, 1.99766, 1.99693 und fallen leicht unter 2. Das passt zu einer Konvergenz, die langsamer ist als jedes Potenzgesetz, wie in B6 | numerisch | `analyse.py` | — |
| B10 | eᵣ folgt keinem einfachen Gesetz. Der freie Grenzfit liefert L=4.118, 4.178, 4.229, also über 4 und mit steigender Tendenz. Der lokale Exponent von 4−eᵣ steigt (1.185 → 1.386). A(r+c)^−p ist instabil (p=1.47 → 2.42), und keine der Normierungen r ln r, r^1.5, r ln²r ist stabil | negativ | Wenn eᵣ≤4 gilt, muss 4−eᵣ zuletzt schneller fallen als jede Anpassung an r≤8 | — |
| B11 | (4−eᵣ)/(2−ηᵣ) = 2.308, 2.499, 2.660, 2.777, 2.849, 2.877; die Zuwächse brechen ein (0.19, 0.16, 0.12, 0.07, 0.03). Das Verhältnis hat wahrscheinlich bald ein Maximum | numerisch | — | — |
| B12 | Effektiver Konvergenzradius der Taylorreihe von Fᵣ um eᵣ (Quotient aₖ/aₖ₊₁ bei k≤8): ≈2.9, 1.84, 1.47, 1.29, ≈1.15 für r=4…8. Er schrumpft | numerisch, nur effektiv | bei r=8 nur k≤5; nicht asymptotisch | Hinweis auf eine wandernde Singularität |
| B13 | η₉ war mit den Regeln dieser Runde nicht rechenbar | negativ (Aufwand) | H₈ oberhalb η₈ braucht die komplexe Zwei-Bogen-Octation (98–241 s je Basis), H₉ darauf 200–340 s. Nötig wären mindestens 3 Basen plus Cauchy-Jets, also deutlich über 5 min je Lauf und über 30 min insgesamt | — |
| B14 | Komplexe Abtastung der Octation um e₈ mit Radius 0.1 bricht ab („o8knhex: Newton bleibt nicht im Bogen“); 0.03 und 0.05 funktionieren | negativ (technisch) | `cache/jet_r8_d20_M32_rho0.1.json` | — |

## Beweis von B3

Setze z=x−eᵣ. Dann ist Fᵣ(eᵣ+z)=eᵣ+z+a₂z²+a₃z³+O(z⁴) mit a₂=βᵣ>0. Mit w=a₂z wird daraus g(w)=w+w²+c w³+O(w⁴), c=a₃/a₂², und ι:=1−c. Für Φ(w)=−1/w+ι·log(−w) gilt auf w<0:
−1/g(w)=−1/w+1+(c−1)w+O(w²) und log(−g(w))=log(−w)+w+O(w²). Also Φ(g(w))−Φ(w)=1+(c−1+ι)w+O(w²)=1+O(w²).
Für den Orbit wₖ↑0 (Konvergenz von unten, die Voraussetzung aus M3) ist wₖ=O(1/k). Die Fehler O(wₖ²) sind deshalb summierbar, und es folgt Φ(wₖ)=k+C+o(1). Daraus wird −1/wₖ=k+C−ι·log(−wₖ)+o(1)=k+ι ln k+C+o(1). Zurück übersetzt:
eᵣ−Hᵣ₊₁(ηᵣ,k)=−wₖ/βᵣ=1/(βᵣ(k+ιᵣ ln k+Cᵣ+o(1))). ∎

Prüfung bei r=3: β=1/(2e), ι=1/3. Das ergibt e−ᵏη=2e/k−(2e/3)ln k/k²+O(1/k²). Das ist die bekannte langsame Konvergenz an der Basis e^(1/e).
Bedeutung: ιᵣ≈−0.6 heißt, dass der Logarithmusterm bei hohen Rängen das Vorzeichen wechselt. Die nächste Rechenart erreicht eᵣ dort *schneller* als die reine 1/(βk)-Form. Das ist wichtig, wenn man ηᵣ aus endlichen Höhen schätzt (vgl. M3).

## Rechenweg und Reproduktion

- `runtime/`: unveränderte Kopie der eingefrorenen Engine aus `eta_hyperoperationen/rang_gamma/runtime` (Snapshot von fatou_fork.gp plus hyper.gp, even8.gp). Die Kombinationsdatei entsteht nur hier; es gibt keine geteilten Zwischenspeicher.
- `jets.py --rank r --dps D --M M --rho ρ`: initialisiert Hᵣ bei ηᵣ (`hinit`), tastet hval(r,·) auf dem Kreis |x−eᵣ|=ρ ab und berechnet aₖ, β und ι. Ausgabe in `cache/jet_r*_d*_M*_rho*.json`. Laufzeiten: r≤5 unter 5 s, r=6,7 13–42 s, r=8 75–108 s.
- `analyse.py`: alle Tabellen, Fits, Normierungen und η₉-Vorhersagen (Laufzeit < 2 s). Eingabe ist `eta_values.csv` (Kopie aus `forschung/eta_hyperoperationen`).
- Taylorkoeffizienten r=7 (ρ=0.1): 1, 0.49615, 0.39528, 0.29364, 0.22858, 0.17778, 0.13697, 0.10626. r=8: 1, 0.53682, 0.45903, 0.37837, 0.33146.

## Vorschläge für Grapheinträge

- `kn.satz.kritische_rate_log`: Verfeinerte kritische Rate mit ιᵣ ln k — Satz — bewiesen (Voraussetzungen wie kn.satz.kritische_rate)
- `kn.beob.iota_raenge`: Iteratives Residuum der Berührung, r=3…8, Minimum bei r=7 — Beobachtung — numerisch
- `kn.beob.beta_raenge`: βᵣ für r=3…8 mit unabhängiger Kontrolle von β₇, β₈ — Beobachtung — numerisch
- `kn.verm.eta_r_log`: 2−ηᵣ ~ C/(r ln r), C≈2.20, Prüfung über η₉≈1.88883 gegen 1.88914 — Vermutung
- `kn.gescheitert.e_potenzgesetz`: 4−eᵣ passt in kein Potenzgesetz und der freie Grenzwert driftet über 4 — gescheitert / negativ

## Offene Fragen

1. η₉ mit 4 Stellen würde B6 gegen das Potenzgesetz entscheiden. Dafür braucht es einen einzigen längeren Lauf (H₈ komplex, H₉ regulär, Sekante in b), am besten auf dem Server und nicht auf dem Laptop.
2. Hat ιᵣ einen Grenzwert ι_∞ (bei etwa −0.6)? Liegt nach dem Minimum bei r=7 ein Paritätseffekt vor (gerade Ränge komplex, ungerade regulär gebaut)?
3. Was ist die Singularität, die den Konvergenzradius um eᵣ schrumpfen lässt (B12)? Ist sie reell und liegt sie rechts von eᵣ, und läuft eᵣ+Rᵣ gegen 4?
4. Eine Herleitung für C/(r ln r): Nach der Heuristik ηᵣ₊₁−ηᵣ≈G/Sᵣ mit Sᵣ=∂_bHᵣ(b,eᵣ) bräuchte man Sᵣ~r² ln r. Sᵣ lässt sich für r≤7 billig über zwei nahe Basen messen; das wurde hier nicht mehr gemacht.
