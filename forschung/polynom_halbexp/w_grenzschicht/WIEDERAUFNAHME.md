# Strang W — Wiederaufnahme (angehalten am 11.10.2026)

**Lesen zuerst:** `BERICHT.md` — Kurzfassungen ganz oben (neueste zuerst), dann §27 (W12-Stand), §22–§26 (Reduktion, Weg A),
§19–§20 (Kegel, Modell). Problemblatt: `OFFENES_PROBLEM.md` (deutsch), `OPEN_PROBLEM_en.md` (englisch).

**Stand in einem Satz:** (H₀) ist auf die globale Eckfrage reduziert. Die n-te-Wurzel-Stufe allein trägt nicht (Normierung = (E_L)); weiter geht es nur
auf der starken Ebene.

**Nächster Schritt:** Die Schweiss-RH-Aufgabe aufsetzen und mit O abstimmen (BRETT 04:35: Phase N·c₁ζ^{κ₁}, Amplitude |ζ|^{α}, Re α = −0.476 oder −0.642,
Skalen 1/N und N^{−0.894}). Dafür:
- (y, T) auf der geschweissten Fläche, mit u als Inhomogenität;
- äussere Parametrix aus H auf ℳ (Referenzmodell `w40_modell_ref.npz`, Loader `w40_lade.py`);
- lokale Parametrix in L (O).

**Wichtige Skripte:** w40_freier_rand_ecke.py (Modell; Argumente im Kopf); w41 (H*-Prüfung), w42/w43 (Nullstellen), w44–w47 (Testvektor,
Konstanten, Kreisbedingung, D6-Konstante), w50 (Skalenkreise), w51–w53 (biorthogonale Interpolanten).
Python: C:/Users/janis/miniconda3/python.exe; arb mit ctx.cap = N + 10 setzen.
