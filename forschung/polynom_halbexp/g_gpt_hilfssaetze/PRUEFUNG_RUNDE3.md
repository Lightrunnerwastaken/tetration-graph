# Prüfungen der dritten Runde

Zeitpunkt: 2026-10-11T12:09:09+02:00.

- Exakte SymPy-Prüfung: Rang-eins-Zuwachs mit voller Quelle, monische Biorthogonalität, Companion-Eigenbasis und Determinante.
- 512-Bit-Intervallgeometrie: beide kritischen Werte haben Betrag <0,476; die kubischen Eigenwerte in L sind getrennt.
- Monodromie: analytischer 3-Zykel nach der Deformation zu ∞; separate numerische Wurzelverfolgung mit 721 Punkten bestätigt ihn.
- Kubische Zuwächse: N=200,300,400,600,800,1000, feste γ*-Punkte und 24 Kreiswinkel, hochpräzise Approx-Solves.
- Ganze Eckkreise: dieselben sechs Grade, 900…2500 Bit, geprüfte Matrixarithmetik und vollständige endliche Taylor-Majoranten.
- Für jeden dieser sechs Grade ist √N sup|𝒟_N|<3,13 als echte Intervallaussage verifiziert.
- Zusätzlich direkte G-Majorante auf vollständigen Kreisen bei N=200 und N=1000; Status und Werte in den jeweiligen JSON-Dateien.
- Geprüfte Anfangsnorm bei N=199 auf der 200-Scheibe für das bedingte Summenargument.
- Kneser-Gegenprobe: N=200,300,400,600, gleiche festen γ*-Punkte; Approx-Solves, keine Propagation der Koeffizientenradien.
- Numerischer Präzisionsvergleich sowie aufeinanderfolgende Grade 200…210; Unterschiede stehen im Bericht.
- Verwendete W-Inversen/Singulärwerte wurden nur gelesen. Keine neue SVD.
- Syntaxprüfung aller Python-Dateien im Paket erfolgreich. Alle finalen JSON-Dateien lesbar.

Ein Versuch mit geprüfter ungepreconditionter LU bei N=200 und 900 Bit lieferte
`ZeroDivisionError: singular matrix in solve()` durch unzureichende Intervallschärfe.
Das ist kein Beleg einer tatsächlichen Singularität. Die anschließend erfolgreichen
preconditionierten Prüfungen desselben exakten Modells zertifizieren die Invertibilität.
Die bereits gespeicherten Ergebnisse wurden beim Wechsel vom doppelt berechneten
Preconditioner zur gemeinsamen geprüften Inversen erhalten; Zwischenläufe wurden gezielt beendet.
Eine direkte Float64-Auswertung der Randzeile in der Bucht war wegen Auslöschung unbrauchbar
(scheinbarer Zuwachs etwa 5,20·10^57 bei N=1000). Der verworfene Probe-Code ist ausdrücklich
als solcher unter `archiv/` erhalten; seine Ergebnisse werden nicht als Messwerte verwendet.

Keine All-n-Schlussfolgerung, kein globales Kurvenzertifikat für γ*, keine Kneser-Übertragung,
keine K2-Läufe. Unabhängige P-Prüfung bleibt offen.

Versionen: {"python-flint": "0.9.0", "sympy": "1.14.0", "numpy": "2.2.6", "matplotlib": "3.10.9"}.
