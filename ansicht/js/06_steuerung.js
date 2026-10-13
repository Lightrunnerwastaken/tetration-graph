// ===================================================================== Steuerung
function pfadZeigen() {
  const p = $("#pfad");
  if (zustand.ansicht === "theorie") {
    const buchModus = zustand.theorieModus !== "graph";
    p.innerHTML = `<span class="modus-wahl" role="group" aria-label="Theorie-Ansicht">
      <button data-theoriemodus="buch" aria-pressed="${buchModus}">Lehrbuch</button>
      <button data-theoriemodus="graph" aria-pressed="${!buchModus}">Abhängigkeiten</button></span>
      <span class="klein">${buchModus ? "von den Grundlagen bis zur Forschung, Kapitel für Kapitel" : "Beweisstufen, links Grundlagen"}</span>`;
    return;
  }
  if (zustand.ansicht === "ideen") { const h = ideen.pfadHtml(); if (p.innerHTML !== h) p.innerHTML = h; return; }
  if (zustand.ansicht === "forschung") { p.innerHTML = `<span>Forschung: Lücken, Experimente, Korrekturen, nächste Fragen</span>`; return; }
  const o = karte.ort();
  const teile = [`<button data-karte="ganz">Karte</button>`];
  const k = zustand.gewaehlt && K.get(zustand.gewaehlt);
  const gid = k ? k.gebiet : o.gebiet, sid = k ? k._strang : o.strang;
  if (gid) teile.push(`<button data-gebiet="${gid}">${esc(GEB.get(gid).name)}</button>`);
  if (sid && STR.has(sid)) teile.push(`<button data-strang="${sid}">${esc(STR.get(sid).name)}</button>`);
  if (k) teile.push(`<span>${esc(kurz(reinText(k.titel), 48))}</span>`);
  const html = teile.join(" › ");
  if (p.innerHTML !== html) p.innerHTML = html;
}
setInterval(() => { if (zustand.ansicht === "karte") pfadZeigen(); }, 250);

function filterLeiste() {
  const f = $("#filter");
  if (zustand.ansicht === "forschung" || zustand.ansicht === "ideen" || (zustand.ansicht === "theorie" && zustand.theorieModus !== "graph")) { f.innerHTML = ""; return; }
  const aktiv = (GRUPPEN.length - zustand.filter.size) + (4 - zustand.ebenen.size);
  f.innerHTML = Object.entries(NIVEAU_NAME).map(([n, name]) => `<button data-niveau="${n}" class="niveau-knopf" aria-pressed="${zustand.niveaus.has(n)}" title="${esc(S.niveau[n])}">${name}</button>`).join("")
    + `<button data-filtermehr="1" aria-expanded="${!!zustand.filterOffen}">${zustand.filterOffen ? "Weniger" : "Weitere Filter"}${aktiv ? ` (${aktiv} aus)` : ""}</button>`
    + (zustand.filterOffen ? `<span class="filter-umbruch"></span>`
      + GRUPPEN.map((g) => `<button data-filter="${g.id}" aria-pressed="${zustand.filter.has(g.id)}"><span class="punkt" style="background:var(${g.v})"></span>${esc(g.name)}</button>`).join("")
      + ["A", "B", "C", "X"].map((e) => `<button data-ebene="${e}" aria-pressed="${zustand.ebenen.has(e)}" title="${esc(S.ebenen[e])}">Ebene ${e}</button>`).join("") : "");
}
function filterAnwenden() { theorie.filter(); karte.anstossen(); }

function setzeAnsicht(a) {
  zustand.ansicht = a;
  document.querySelectorAll(".ansichten button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.ansicht === a)));
  document.body.dataset.ansicht = a;
  document.body.dataset.modus = zustand.theorieModus || "buch";
  $("#tipp").hidden = true;
  if (a !== "ideen") ideen.verberge();
  if (a === "ideen") { karte.verberge(); theorie.verberge(); buch.verberge(); $("#forschung").hidden = true; ideen.zeige(); }
  else if (a === "karte") { theorie.verberge(); buch.verberge(); $("#forschung").hidden = true; karte.zeige(); }
  else if (a === "theorie") {
    karte.verberge(); $("#forschung").hidden = true;
    if (zustand.theorieModus === "graph") { buch.verberge(); theorie.zeige(); } else { theorie.verberge(); buch.zeige(); }
  }
  else { karte.verberge(); theorie.verberge(); buch.verberge(); $("#forschung").hidden = false; zeichneForschung(); }
  pfadZeigen();
  filterLeiste();
  hashSetzen();
}

function waehle(id, opt = {}) {
  zustand.gewaehlt = id;
  if (!id) {
    karte.fokus(null);
    theorie.hervorheben(null);
    tafelSchliessen();
    ideen.anstossen();
    pfadZeigen();
    hashSetzen();
    return;
  }
  if (zustand.ansicht === "ideen") { zeigeKnoten(id); ideen.zuKnoten(id); pfadZeigen(); hashSetzen(); return; }
  if (zustand.ansicht === "theorie" && zustand.theorieModus !== "graph") {
    zeigeKnoten(id); buch.zuKnoten(id); pfadZeigen(); hashSetzen(); return;
  }
  if (zustand.ansicht === "theorie" && !theorie.passt(id)) setzeAnsicht("karte");
  zeigeKnoten(id);
  if (zustand.ansicht === "karte") karte.fokus(id, opt.fliegen !== false);
  if (zustand.ansicht === "theorie") theorie.hervorheben(id);
  pfadZeigen();
  hashSetzen();
}

function hashSetzen() {
  let h = zustand.gewaehlt ? `#a=${zustand.gewaehlt}&v=${zustand.ansicht}` : `#v=${zustand.ansicht}`;
  if (zustand.ansicht === "theorie") h += zustand.theorieModus === "graph" ? "&m=graph" : (zustand.kapitel ? `&k=${zustand.kapitel}` : "");
  if (location.hash !== h) history.replaceState(null, "", h);
}

document.addEventListener("click", (e) => {
  const t = e.target.closest("[data-idee],[data-id],[data-exp],[data-gebiet],[data-strang],[data-ansicht],[data-filter],[data-ebene],[data-niveau],[data-filtermehr],[data-tafel],[data-karte],[data-theoriemodus],[data-lehrbuch]");
  if (!t || t.closest(".knoten") || t.closest(".inhaltskarte")) return;
  const d = t.dataset;
  if (d.idee) {
    if (zustand.ansicht !== "ideen") setzeAnsicht("ideen");
    if (zustand.gewaehlt) { zustand.gewaehlt = null; tafelSchliessen(); hashSetzen(); }
    ideen.oeffne(d.idee === "_alle" ? null : d.idee);
  } else if (d.id) { if (zustand.ansicht === "forschung") { zustand.gewaehlt = d.id; zeigeKnoten(d.id); hashSetzen(); } else waehle(d.id); }
  else if (d.exp) zeigeExperiment(d.exp);
  else if (d.gebiet) {
    if (zustand.ansicht !== "karte") setzeAnsicht("karte");
    waehleNichts(); zeigeGebiet(d.gebiet); karte.zuGebiet(d.gebiet);
  } else if (d.strang) {
    if (zustand.ansicht !== "karte") setzeAnsicht("karte");
    waehleNichts(); zeigeStrang(d.strang); karte.zuStrang(d.strang);
  } else if (d.karte === "ganz") { waehleNichts(); tafelSchliessen(); karte.ganz(); }
  else if (d.ansicht) setzeAnsicht(d.ansicht);
  else if (d.theoriemodus) { zustand.theorieModus = d.theoriemodus; setzeAnsicht("theorie"); }
  else if (d.lehrbuch) { zustand.theorieModus = "buch"; setzeAnsicht("theorie"); buch.zuKnoten(d.lehrbuch); }
  else if (d.tafel === "zu") tafelSchliessen();
  else if (d.filter) { const s = zustand.filter; s.has(d.filter) ? s.delete(d.filter) : s.add(d.filter); filterLeiste(); filterAnwenden(); }
  else if (d.ebene) { const s = zustand.ebenen; s.has(d.ebene) ? s.delete(d.ebene) : s.add(d.ebene); filterLeiste(); filterAnwenden(); }
  else if (d.filtermehr) { zustand.filterOffen = !zustand.filterOffen; filterLeiste(); }
  else if (d.niveau) { const s = zustand.niveaus; s.has(d.niveau) ? s.delete(d.niveau) : s.add(d.niveau); filterLeiste(); filterAnwenden(); }
});
function waehleNichts() { zustand.gewaehlt = null; karte.fokus(null, false); hashSetzen(); }

$("#tafel-knopf").onclick = () => {
  if (document.body.classList.contains("tafel-offen")) tafelSchliessen();
  else if (zustand.gewaehlt) zeigeKnoten(zustand.gewaehlt);
  else zeigeUebersicht();
};
$("#zoom-plus").onclick = () => (zustand.ansicht === "ideen" ? ideen.zoom(1.4) : zustand.ansicht === "karte" ? karte.zoom(1.6) : theorie.zoom(1.25));
$("#zoom-minus").onclick = () => (zustand.ansicht === "ideen" ? ideen.zoom(1 / 1.4) : zustand.ansicht === "karte" ? karte.zoom(1 / 1.6) : theorie.zoom(0.8));
$("#zoom-fit").onclick = () => (zustand.ansicht === "ideen" ? (waehle(null), ideen.oeffne(null)) : zustand.ansicht === "karte" ? (waehleNichts(), karte.ganz()) : setzeAnsicht("theorie"));

// ------------------------------------------------------------------ Suche
const suche = $("#suche"), treffer = $("#treffer");
let auswahl = 0, trefferListe = [];
function suchen(q) {
  q = q.trim().toLowerCase();
  if (!q) { treffer.hidden = true; return; }
  const w = q.split(/\s+/);
  trefferListe = G.knoten.map((k) => {
    const t = reinText(k.titel).toLowerCase(), body = (reinText(k.aussage) + " " + k.id).toLowerCase();
    let s = 0;
    for (const x of w) { if (t.includes(x)) s += 3; else if (body.includes(x)) s += 1; else return null; }
    return [s, k];
  }).filter(Boolean).sort((a, b) => b[0] - a[0]).slice(0, 14).map((x) => x[1]);
  auswahl = 0;
  treffer.innerHTML = trefferListe.length ? trefferListe.map((k, i) => `<li role="option" aria-selected="${i === 0}" data-wahl="${k.id}">${punkt(k._status_eff)}<span>${esc(reinText(k.titel))}</span><span class="id">${esc(GEB.get(k.gebiet).name)} · ${NIVEAU_NAME[k._niveau] || ""}</span></li>`).join("")
    : `<li class="klein">Keine Aussage gefunden. Andere Begriffe oder eine id wie „kn.satz…“ versuchen.</li>`;
  treffer.hidden = false;
}
suche.addEventListener("input", () => suchen(suche.value));
suche.addEventListener("keydown", (e) => {
  if (e.key === "ArrowDown" || e.key === "ArrowUp") {
    e.preventDefault();
    auswahl = Math.max(0, Math.min(trefferListe.length - 1, auswahl + (e.key === "ArrowDown" ? 1 : -1)));
    treffer.querySelectorAll("li").forEach((li, i) => li.setAttribute("aria-selected", String(i === auswahl)));
  } else if (e.key === "Enter" && trefferListe[auswahl]) { waehle(trefferListe[auswahl].id); treffer.hidden = true; suche.blur(); }
  else if (e.key === "Escape") { treffer.hidden = true; suche.blur(); }
});
treffer.addEventListener("click", (e) => { const li = e.target.closest("[data-wahl]"); if (li) { waehle(li.dataset.wahl); treffer.hidden = true; } });
document.addEventListener("pointerdown", (e) => { if (!e.target.closest(".suche")) treffer.hidden = true; });
document.addEventListener("keydown", (e) => {
  if (e.target.closest?.("input")) return;
  if (e.key === "/") { e.preventDefault(); suche.focus(); suche.select(); }
  else if (e.key === "Escape") {
    if (!$("#legende").hidden) { $("#legende").hidden = true; $("#legende-knopf").setAttribute("aria-expanded", "false"); }
    else if (zustand.ansicht === "ideen") ideen.hoch();
    else if (zustand.gewaehlt) waehle(null);
    else if (document.body.classList.contains("tafel-offen")) tafelSchliessen();
    else if (zustand.ansicht === "karte") karte.ganz();
  } else if (zustand.ansicht === "karte" && (e.key === "+" || e.key === "=")) karte.zoom(1.6);
  else if (zustand.ansicht === "karte" && e.key === "-") karte.zoom(1 / 1.6);
});

// ------------------------------------------------------------------ Legende
function legendeBauen() {
  const beispiel = (art, st) => {
    const kk = { art, _status_eff: st, ebene: "A", _tragweite: 1, titel: "" };
    const s = document.createElementNS(NS, "svg");
    s.setAttribute("viewBox", "-15 -11 30 22");
    const f = form(art, 8, el("g", {}, s));
    const stl = stilFuer(kk);
    f.setAttribute("style", `fill:${stl.fill};stroke:${stl.stroke};stroke-width:1.4;${stl.strich ? "stroke-dasharray:" + stl.strich : ""}`);
    return s.outerHTML;
  };
  $("#legende").innerHTML = `
    <h4>Zoomstufen der Karte</h4>
    <dl class="legende-stufen"><dt>Gebiete</dt><dd>Inseln mit Statusring; Bänder = Querverbindungen, Lichtpunkte fliessen darauf.</dd>
    <dt>Verläufe</dt><dd>Forschungsstränge als Flächen mit einer Zeile, was dort passiert; Licht fliesst von Voraussetzung zu Folgerung.</dd>
    <dt>Aussagen</dt><dd>Titel, wichtige zuerst; Überlappendes bleibt verborgen.</dd>
    <dt>Inhalt</dt><dd>Karten mit der gesetzten Aussage neben den Knoten.</dd></dl>
    <h4>Form = Art der Aussage</h4><div class="legende-gitter">
    ${[["satz", "Satz, Lemma, Korollar"], ["definition", "Definition"], ["konstruktion", "Konstruktion"], ["vermutung", "Vermutung"], ["annahme", "Annahme / Hypothese"],
      ["beobachtung", "Beobachtung (Messung)"], ["gegenbeispiel", "Gegenbeispiel"], ["korrektur", "Korrektur"], ["gescheitert", "Gescheiterter Ansatz"], ["frage", "Offene Frage"]]
      .map(([a, t]) => `${beispiel(a, a === "frage" ? "offen" : "bewiesen")}<span>${t}</span>`).join("")}</div>
    <h4>Farbe = Prüfstatus (abgeleitet entlang der Beweiskette)</h4><div class="legende-gitter">
    ${RINGFOLGE.map((st) => `${beispiel("satz", st)}<span>${esc(S.status[st].label)}</span>`).join("")}</div>
    <h4>Buchstabe im Knoten = Ebene</h4><div class="legende-gitter">${Object.entries(S.ebenen).map(([e, t]) => `<b style="text-align:center">${e}</b><span>${esc(t)}</span>`).join("")}</div>
    <h4>Grösse und Markierungen</h4>
    <p class="klein">Grösse: wie viele Aussagen (transitiv) darauf aufbauen. Roter pulsierender Ring: neu prüfen. Gestrichelter Rand: bedingt oder Skizze. Durchgestrichen: widerlegt.</p>
    <h4>Linien</h4><p class="klein">Blau: Beweis (Pfeil zur Folgerung) · Rot gestrichelt: widerlegt/korrigiert · Grün: entlastet oder stützt · Grau: Spezialfall, verwandt.</p>
    <p class="klein" style="margin-top:12px">Tasten: <b>/</b> Suche · <b>+ −</b> Zoom · <b>Esc</b> zurück · Doppelklick zoomt hinein.</p>`;
}
// Handy: Suche und Werkzeuge hinter je einem Knopf (Kopfzeile in einer Zeile)
function handyMenue(auf) {
  document.body.classList.toggle("menue-offen", auf);
  $("#menue-knopf").setAttribute("aria-expanded", String(auf));
}
function handySuche(auf) {
  document.body.classList.toggle("suche-offen", auf);
  $("#such-knopf").setAttribute("aria-expanded", String(auf));
  if (auf) { handyMenue(false); $("#suche").focus(); }
}
$("#menue-knopf").onclick = (e) => { e.stopPropagation(); handySuche(false); handyMenue(!document.body.classList.contains("menue-offen")); };
$("#such-knopf").onclick = (e) => { e.stopPropagation(); handySuche(!document.body.classList.contains("suche-offen")); };
document.addEventListener("click", (e) => {
  if (document.body.classList.contains("menue-offen") && (!e.target.closest(".werkzeuge") || e.target.closest("button"))) handyMenue(false);
  if (document.body.classList.contains("suche-offen") && !e.target.closest(".suche")) handySuche(false);
});
$("#treffer").addEventListener("click", () => handySuche(false));

$("#legende-knopf").onclick = () => {
  const l = $("#legende");
  l.hidden = !l.hidden;
  $("#legende-knopf").setAttribute("aria-expanded", String(!l.hidden));
};
$("#thema-knopf").onclick = () => {
  const dunkel = document.documentElement.dataset.theme
    ? document.documentElement.dataset.theme === "dark"
    : matchMedia("(prefers-color-scheme: dark)").matches;
  document.documentElement.dataset.theme = dunkel ? "light" : "dark";
  try { localStorage.setItem("tg-thema", document.documentElement.dataset.theme); } catch (e) { /* ohne Speicher */ }
};
try { const t = localStorage.getItem("tg-thema"); if (t) document.documentElement.dataset.theme = t; } catch (e) { /* ohne Speicher */ }

// ------------------------------------------------------------------ Start
legendeBauen();
{
  const h = new URLSearchParams(location.hash.slice(1));
  const v = ["karte", "theorie", "forschung"].includes(h.get("v")) ? h.get("v") : "ideen";
  zustand.theorieModus = h.get("m") === "graph" ? "graph" : "buch";
  if (h.get("k")) zustand.kapitel = h.get("k");
  setzeAnsicht(v);
  const a = h.get("a");
  if (a && K.has(a)) {
    if (v === "forschung") { zustand.gewaehlt = a; zeigeKnoten(a); }
    else waehle(a);
  }
}
