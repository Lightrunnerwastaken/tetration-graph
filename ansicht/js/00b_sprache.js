// ===================================================================== Sprache (DE/EN)
// Die Oberfläche ist deutsch geschrieben. Im englischen Modus werden (1) die Daten mit den vorhandenen `en`-Feldern überlagert,
// bevor die anderen Module sie lesen, und (2) alle Texte, die in die Seite gelangen (Textknoten, einige Attribute, Canvas-Text),
// über EN_WOERTER/EN_MUSTER ins Englische gesetzt. Umschalten = Wahl merken und neu laden.
const SPRACHE = (() => {
  try {
    const q = new URLSearchParams(location.search).get("lang");
    if (q === "en" || q === "de") return q;
  } catch (e) { /* ohne URL-Parameter */ }
  try {
    const s = localStorage.getItem("tetra-sprache");
    if (s === "en" || s === "de") return s;
  } catch (e) { /* localStorage nicht verfügbar */ }
  return "en";                     // Standard für alle (öffentliche Bibliothek); Deutsch über den Knopf oder ?lang=de, die Wahl wird gemerkt
})();
const EN = SPRACHE === "en";
document.documentElement.lang = SPRACHE;

// ------------------------------------------------------------------ Daten überlagern
function sprachDaten(G) {
  const ersetze = (ziel, quelle, felder) => { for (const f of felder) if (quelle && typeof quelle[f] === "string" && quelle[f].trim()) ziel[f] = quelle[f]; };
  for (const k of G.knoten) {
    const en = k.en;
    if (en) {
      ersetze(k, en, ["titel", "aussage", "bereich", "einfach", "naechste_frage", "gegenbeispielsuche", "bedingungen_text"]);
      if (en.beweis && k.beweis) {
        ersetze(k.beweis, en.beweis, ["text"]);
        for (const [i, l] of (en.beweis.luecken || []).entries()) if (typeof l === "string" && l.trim() && k.beweis.luecken) k.beweis.luecken[i] = l;
      }
      for (const [i, f] of (en.formen || []).entries()) if (k.formen && k.formen[i] && typeof k.formen[i] === "object") ersetze(k.formen[i], f, ["notiz"]);
      for (const [i, v] of (en.voraussetzungen || []).entries()) if (k.voraussetzungen && k.voraussetzungen[i]) ersetze(k.voraussetzungen[i], v, ["text", "notiz"]);
      for (const [i, v] of (en.versionen || []).entries()) if (k.versionen && k.versionen[i]) ersetze(k.versionen[i], v, ["aenderung", "grund"]);
      for (const [i, v] of (en.nachgeprueft || []).entries()) if (k.nachgeprueft && k.nachgeprueft[i]) ersetze(k.nachgeprueft[i], v, ["ergebnis"]);
      for (const [i, v] of (en.quellen || []).entries()) if (k.quellen && k.quellen[i] && typeof k.quellen[i] === "object") ersetze(k.quellen[i], v, ["stelle"]);
      if (en.pruefung) { k.pruefung = k.pruefung || {}; ersetze(k.pruefung, en.pruefung, Object.keys(en.pruefung)); }
      if (en.bewertung) { k.bewertung = k.bewertung || {}; ersetze(k.bewertung, en.bewertung, ["richtigkeit", "nutzen", "recherche"]); }
    }
    if (k.en_theorie && k.theorie) ersetze(k.theorie, k.en_theorie, ["worum", "genau", "warum", "beispiel", "einordnung"]);
  }
  for (const g of G.gebiete || []) ersetze(g, g.en, ["name", "kurz"]);
  for (const s of G.straenge || []) ersetze(s, s.en, ["name", "verlauf"]);
  for (const e of G.experimente || []) {
    ersetze(e, e.en, ["titel", "frage"]);
    const r = e.en && e.en.ergebnis;
    if (r && e.ergebnis) {
      ersetze(e.ergebnis, r, ["zusammenfassung"]);
      for (const [i, p] of (r.pruefungen || []).entries()) if (e.ergebnis.pruefungen && e.ergebnis.pruefungen[i]) ersetze(e.ergebnis.pruefungen[i], p, ["name"]);
    }
  }
  for (const q of Object.values(G.quellen || {})) ersetze(q, q.en, ["notiz", "titel"]);
  for (const o of Object.values(G.objekte || {})) ersetze(o, o.en, ["name"]);
  for (const x of [...((G.ideen || {}).haupt || []), ...((G.ideen || {}).unter || [])]) ersetze(x, x.en, ["name", "kurz"]);
  for (const kap of (G.lehrbuch || {}).kapitel || []) {
    if (!kap.en) continue;
    ersetze(kap, kap.en, ["titel", "untertitel"]);
    for (const [i, t] of (kap.en.teile || []).entries()) if (kap.teile && kap.teile[i]) ersetze(kap.teile[i], t, ["abschnitt", "text", "unterschrift"]);
  }
}
if (EN) sprachDaten(window.GRAPH);

// ------------------------------------------------------------------ Texte übersetzen
const WOERTER = new Map(Object.entries(EN_WOERTER));
const MUSTER = EN_MUSTER.map(([re, rep]) => [re instanceof RegExp ? re : new RegExp(re), rep]);
function uebersetzeKern(kern) {
  const treffer = WOERTER.get(kern);
  if (treffer !== undefined) return treffer;
  for (const [re, rep] of MUSTER) if (re.test(kern)) return kern.replace(re, rep);
  return kern;
}
function uebersetzeText(s) {
  if (typeof s !== "string" || !s) return s;
  const m = /^(\s*)([\s\S]*?)(\s*)$/.exec(s);
  const kern = m[2];
  if (!kern) return s;
  const treffer = WOERTER.get(kern);
  if (treffer !== undefined) return m[1] + treffer + m[3];
  for (const [re, rep] of MUSTER) if (re.test(kern)) return m[1] + kern.replace(re, rep) + m[3];
  // Aufzählungen "A · B · C" (Beziehungen, Zeilen der Tafel): Teile einzeln übersetzen
  if (kern.includes(" · ")) {
    const teile = kern.split(" · "), neu = teile.map((x) => uebersetzeKern(x));
    if (neu.some((x, i) => x !== teile[i])) return m[1] + neu.join(" · ") + m[3];
  }
  return s;
}
const t = (s) => (EN ? uebersetzeText(s) : s);

// Heuristik für den Test: sieht ein Text noch deutsch aus?
const NAMEN_MIT_UMLAUT = /Schr[öo]der|B[öo]ttcher|Erd[őo]s|M[üu]ller|Hölder|Weierstra[ßs]|Lévy|Écalle|Gödel|Dzyadyk|Kühn|Bürger|Gauß|Gauss|Möbius|Mobius|Cesàro|Pólya|Poincaré|Fréchet|Cauchy|Lambert|Thron|Löwner|Lowner|Kármán|Szegő|Szeg|Eötvös/g;
const DEUTSCHE_WOERTER = /(^|[^\p{L}_])(und|der|die|das|nicht|mit|für|von|zu|im|auf|alle|noch|keine|oder|wird|sind|eine|einer|zum|zur|wenn|dass|aber|auch|bei|nach|aus|über|unter|ohne|nur|wie|werden|haben|gibt|jede|jeder|jedes|Aussage|Aussagen|Beweis|Voraussetzung|Quelle|Quellen|Satz|Sätze|Vermutung)($|[^\p{L}_])/u;
function sieht_deutsch_aus(s) {
  if (typeof s !== "string") return false;
  const ohneNamen = s.replace(NAMEN_MIT_UMLAUT, " ");
  return /[äöüÄÖÜß]/.test(ohneNamen) || DEUTSCHE_WOERTER.test(ohneNamen.replace(/\$[^$]*\$/g, " "));
}
window.__sprache = { sprache: SPRACHE, uebersetzeText, sieht_deutsch_aus, canvas: [] };

if (EN) {
  const ATTR = ["title", "aria-label", "placeholder", "alt"];
  const bearbeite = (knoten) => {
    if (knoten.nodeType === 3) {
      if (knoten.parentNode && /^(SCRIPT|STYLE)$/.test(knoten.parentNode.nodeName)) return;
      const neu = uebersetzeText(knoten.nodeValue);
      if (neu !== knoten.nodeValue) knoten.nodeValue = neu;
    } else if (knoten.nodeType === 1) {
      if (/^(SCRIPT|STYLE)$/.test(knoten.nodeName)) return;
      for (const a of ATTR) if (knoten.hasAttribute(a)) { const alt = knoten.getAttribute(a), neu = uebersetzeText(alt); if (neu !== alt) knoten.setAttribute(a, neu); }
      for (const kind of knoten.childNodes) bearbeite(kind);
    }
  };
  const beobachter = new MutationObserver((liste) => {
    for (const m of liste) {
      if (m.type === "childList") m.addedNodes.forEach(bearbeite);
      else if (m.type === "characterData") bearbeite(m.target);
      else if (m.type === "attributes") bearbeite(m.target);
    }
  });
  const starte = () => {
    bearbeite(document.documentElement);
    document.title = uebersetzeText(document.title);
    beobachter.observe(document.documentElement, { childList: true, subtree: true, characterData: true, attributes: true, attributeFilter: ATTR });
  };
  // Canvas: Text wird nicht ins DOM gesetzt
  const ctx = CanvasRenderingContext2D.prototype;
  for (const fn of ["fillText", "strokeText", "measureText"]) {
    const orig = ctx[fn];
    ctx[fn] = function (text, ...rest) {
      const neu = uebersetzeText(String(text));
      if (sieht_deutsch_aus(neu)) window.__sprache.canvas.push(neu);
      return orig.call(this, neu, ...rest);
    };
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", starte); else starte();
}

// ------------------------------------------------------------------ Umschalter
function sprachKnopf() {
  const knopf = document.getElementById("sprache-knopf");
  if (!knopf) return;
  knopf.textContent = EN ? "DE" : "EN";
  knopf.title = EN ? "Auf Deutsch umschalten" : "Switch to English";
  knopf.setAttribute("aria-label", knopf.title);
  knopf.addEventListener("click", () => {
    try { localStorage.setItem("tetra-sprache", EN ? "de" : "en"); } catch (e) { /* nicht speicherbar */ }
    const u = new URL(location.href);
    u.searchParams.set("lang", EN ? "de" : "en");
    location.href = u.toString();
  });
}
if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", sprachKnopf); else sprachKnopf();
