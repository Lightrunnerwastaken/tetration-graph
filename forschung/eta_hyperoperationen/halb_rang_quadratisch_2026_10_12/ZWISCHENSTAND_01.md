# Zwischenstand: quadratischer halber Rang, 12.10.2026

Status: Forschung am fraktionalen Rang. Der alte H₀-Beweis bleibt pausiert.
Der gelieferte Forschungsplan wurde als Vorschlag gelesen, nicht als zusätzliche Arbeitsanweisung.

**Bewiesen im deklarierten endlichen Interpolationsmodell:** Ein nichtlinearer halber Rangschritt existiert bis zur zweiten Ordnung. Die quadratische Carleman-Matrix ist unipotent; ihre endliche Binomialpotenz liefert eine Ranghalbgruppe modulo kubischer Terme. Zwei halbe Schritte reproduzieren den vollen Schritt in diesem Sinn.

**Bewiesener skalarer Grenzwert:** Für v(x)=x(x−1) und Höhe x=3 ist der quadratische Koeffizient des halben Schritts

    Q_(1/2)(v,v)(3) = 341/12 − 35 log(2)
                    = 4.156515347068580837… .

Für Modellgrösse K handelt es sich genau um die ersten K Terme einer Newtonreihe. Der positive Rest ist für K≥4 höchstens

    8/(3√π) · (K−3)^(−9/2).

**Wichtige Begrenzung:** Der interpolierte quadratische Koeffizient bei x=0.5 driftet: K=16 ergibt +0.00700266, K=32 −0.00767058, K=64 −0.02808715. Das ist ein Warnsignal für die Fortsetzung zu einer Funktion und kein bewiesener Nichtexistenzsatz. Weitere Präzisionskontrollen laufen.

**Vergleich, keine Hyperoperation:** Die einfache Mittelung (b^x+T_b(x))/2 hat eine numerische Tangenz bei b≈1.54254852987747. Unser zuvor festgelegter Fit sagt bei Rang 3.5 etwa 1.55681727674051 voraus. Diese Differenz zeigt, dass der Fit eine Rangfortsetzung nicht auswählt. Der Mittelwert ist kein nachgewiesenes H_3.5.

Nächste Prüfung: Präzisionskontrolle der Drift, explizite mathematische Herleitung und dokumentierte Voraussetzungen für eine echte analytische Ranghalbgruppe.
