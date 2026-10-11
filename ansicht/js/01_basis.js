// ===================================================================== Basis
// Daten, Vokabular, Text/Formeln, SVG-Formen. Alle Module teilen diesen Gültigkeitsbereich.
const G = window.GRAPH;
const S = G.schema;
const NS = "http://www.w3.org/2000/svg";
const $ = (s, el = document) => el.querySelector(s);
const K = new Map(G.knoten.map((k) => [k.id, k]));
const GEB = new Map(G.gebiete.map((g) => [g.id, g]));
const EXP = new Map(G.experimente.map((e) => [e.id, e]));
const STR = new Map((G.straenge || []).map((s) => [s.id, s]));
const Q = G.quellen;
const BEWEGUNG = !matchMedia("(prefers-reduced-motion: reduce)").matches;

const STATUSVAR = {
  definition: "--s-definition", bewiesen: "--s-bewiesen", bewiesen_literatur: "--s-literatur", entlastet: "--s-bewiesen",
  bedingt: "--s-bedingt", skizze: "--s-skizze", numerisch: "--s-numerisch", vermutung: "--s-vermutung",
  annahme: "--s-annahme", offen: "--s-offen", gescheitert: "--s-gescheitert", widerlegt: "--s-widerlegt", ueberholt: "--s-ueberholt",
};
const farbe = (st) => `var(${STATUSVAR[st] || "--s-offen"})`;
const GRUPPEN = [
  { id: "bewiesen", name: "bewiesen", status: ["bewiesen", "bewiesen_literatur", "entlastet", "definition"], v: "--s-bewiesen" },
  { id: "bedingt", name: "bedingt / Skizze", status: ["bedingt", "skizze"], v: "--s-bedingt" },
  { id: "offen", name: "Vermutung / Annahme / offen", status: ["vermutung", "annahme", "offen"], v: "--s-vermutung" },
  { id: "messung", name: "Beobachtung", status: ["numerisch"], v: "--s-numerisch" },
  { id: "tot", name: "widerlegt / überholt / gescheitert", status: ["widerlegt", "ueberholt", "gescheitert"], v: "--s-widerlegt" },
];
const gruppeVon = (st) => (GRUPPEN.find((g) => g.status.includes(st)) || GRUPPEN[2]).id;
const RINGFOLGE = ["definition", "bewiesen_literatur", "bewiesen", "entlastet", "bedingt", "skizze", "numerisch",
  "vermutung", "annahme", "offen", "gescheitert", "ueberholt", "widerlegt"];
const KANTENKLASSE = {
  stuetzt_sich_auf: "beweis", folgt_unmittelbar_aus: "beweis", widerlegt: "rot", korrigiert: "rot",
  entlastet: "gruen", stuetzt: "evidenz", spezialfall_von: "struktur", verallgemeinert: "struktur",
  ersetzt: "struktur", uebertragbar_auf: "struktur", analog_zu: "leicht", motiviert: "leicht", verwandt: "leicht",
};
const THEORIE_ARTEN = new Set(["definition", "satz", "lemma", "korollar", "konstruktion", "vermutung", "annahme"]);

const zustand = {
  ansicht: "karte", gewaehlt: null,
  filter: new Set(GRUPPEN.map((g) => g.id)), ebenen: new Set(["A", "B", "C", "X"]),
  niveaus: new Set(["einstieg", "vertiefung", "forschung"]),
  theorieModus: "buch", kapitel: null,
};
const NIVEAU_NAME = { einstieg: "Einstieg", vertiefung: "Vertiefung", forschung: "Forschung" };
const sichtbar = (k) => zustand.filter.has(gruppeVon(k._status_eff)) && zustand.ebenen.has(k.ebene) && zustand.niveaus.has(k._niveau || "vertiefung");

// Nachbarschaft über alle Kantentypen
const NACHBARN = new Map(G.knoten.map((k) => [k.id, new Set()]));
for (const e of G.kanten) { NACHBARN.get(e.von)?.add(e.nach); NACHBARN.get(e.nach)?.add(e.von); }

// ------------------------------------------------------------------ Text + Mathe
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const ID_RE = /\b(?:gr|fg|it|fp|kl|kn|bw|ar|br|exp)\.[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*/g;
const LINK_RE = /\[\[([a-z]{1,3}\.[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*)(?:\|([^\]]+))?\]\]/g;
function verweisKnopf(id, label) {
  if (K.has(id)) return `<button class="verweis" data-id="${id}">${label ? mathe(label) : esc(reinText(K.get(id).titel))}</button>`;
  if (EXP.has(id)) return `<button class="verweis" data-exp="${id}">${label ? mathe(label) : esc(EXP.get(id).titel)}</button>`;
  return esc(label || id);
}
function textTeil(s) {
  // [[id|Text]] zuerst durch Platzhalter ersetzen, damit die id-Erkennung nicht in Attribute greift
  const platz = [];
  s = String(s ?? "").replace(LINK_RE, (m, id, label) => { platz.push(verweisKnopf(id, label)); return `\u0000${platz.length - 1}\u0000`; });
  let h = esc(s).replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>")
    .replace(/(^|[\s(„])\*(?!\s)([^*\n]+?)\*(?=[\s.,;:!?)“]|$)/g, "$1<em>$2</em>");
  h = h.replace(ID_RE, (m) => {
    if (K.has(m)) return `<button class="verweis" data-id="${m}">${esc(K.get(m).titel)}</button>`;
    if (EXP.has(m)) return `<button class="verweis" data-exp="${m}">${esc(EXP.get(m).titel)}</button>`;
    return m;
  });
  h = h.replace(/\u0000(\d+)\u0000/g, (m, n) => platz[+n]);
  return h.replace(/\n\s*\n/g, "<br><br>").replace(/\n/g, " ");
}
// Titel sind Klartext mit gelegentlichen Hochzahlen wie x^x, e^{1/e}, n^{−3/2}: diese hochstellen, $…$ wie üblich setzen
function titelHtml(t) {
  if (/\$/.test(t || "")) return mathe(t);
  return esc(t).replace(/([\w)])\^\{([^}]+)\}/g, "$1<sup>$2</sup>").replace(/([\w)])\^([\w−-]+)/g, "$1<sup>$2</sup>");
}
// Kleines Markdown für Lehrbuchtexte: Absätze, ### Zwischentitel, Listen, > Merke, $$Formel$$, ![[bild:name]] Unterschrift
function bildHtml(name, unterschrift) {
  const svg = (G.bilder || {})[name] || (EXP.has(name) ? EXP.get(name).svg : "");
  if (!svg) return `<p class="klein">Bild ${esc(name)} fehlt.</p>`;
  return `<figure class="lb-bild">${svg}${unterschrift ? `<figcaption>${mathe(unterschrift)}</figcaption>` : ""}</figure>`;
}
function md(s) {
  if (s == null) return "";
  const bloecke = [];
  let puffer = [];
  for (const zeile of String(s).replace(/\r/g, "").split("\n")) {
    const offen = (puffer.join("\n").match(/\$\$/g) || []).length % 2 === 1;   // $$…$$ über Leerzeilen zusammenhalten
    if (!zeile.trim() && !offen) { if (puffer.join("").trim()) bloecke.push(puffer.join("\n").trim()); puffer = []; continue; }
    if (/^#{2,4}\s/.test(zeile) && !offen) {                                       // Zwischentitel ist immer ein eigener Block
      if (puffer.join("").trim()) bloecke.push(puffer.join("\n").trim());
      bloecke.push(zeile.trim()); puffer = []; continue;
    }
    // Einleitungssatz direkt vor einer Liste und Text direkt nach einer Liste: getrennte Blöcke
    const istListe = (z) => /^([-*]|\d+\.) /.test(z);
    if (puffer.length && !offen && zeile.trim() && !/^\s/.test(zeile) && istListe(zeile) !== istListe(puffer[0])) {
      bloecke.push(puffer.join("\n").trim()); puffer = [];
    }
    puffer.push(zeile);
  }
  if (puffer.join("").trim()) bloecke.push(puffer.join("\n").trim());
  return bloecke.map((t) => {
    let m;
    if ((m = t.match(/^#{2,4}\s+(.*)$/))) return `<h4 class="lb-zwischen">${mathe(m[1])}</h4>`;
    if ((m = t.match(/^!\[\[bild:([A-Za-z0-9_.-]+)\]\]\s*([\s\S]*)$/))) return bildHtml(m[1], m[2]);
    if (/^[-*] /.test(t)) return `<ul class="lb-liste">${t.split(/\n(?=[-*] )/).map((li) => `<li>${mathe(li.replace(/^[-*] /, ""))}</li>`).join("")}</ul>`;
    if (/^\d+\. /.test(t)) return `<ol class="lb-liste">${t.split(/\n(?=\d+\. )/).map((li) => `<li>${mathe(li.replace(/^\d+\. /, ""))}</li>`).join("")}</ol>`;
    if (/^> /.test(t)) return `<div class="lb-merke">${mathe(t.replace(/^> ?/gm, ""))}</div>`;
    const zeilen = t.split("\n");
    if (zeilen.length >= 2 && zeilen.every((z) => /^\|.*\|\s*$/.test(z)) && /^\|[\s:|-]+\|\s*$/.test(zeilen[1])) {   // Tabelle
      const zellen = (z) => {                     // an | trennen, aber nicht in [[id|Text]] und nicht in $…$
        const s = z.trim().slice(1, -1), aus = [];
        let akt = "", tiefe = 0, inMathe = false;
        for (let i = 0; i < s.length; i++) {
          if (s.startsWith("[[", i)) { tiefe++; akt += "[["; i++; continue; }
          if (s.startsWith("]]", i) && tiefe) { tiefe--; akt += "]]"; i++; continue; }
          if (s[i] === "$") inMathe = !inMathe;
          if (s[i] === "|" && !tiefe && !inMathe) { aus.push(akt.trim()); akt = ""; } else akt += s[i];
        }
        aus.push(akt.trim());
        return aus;
      };
      const kopf = zellen(zeilen[0]).map((c) => `<th>${mathe(c)}</th>`).join("");
      const rumpf = zeilen.slice(2).map((z) => `<tr>${zellen(z).map((c) => `<td>${mathe(c)}</td>`).join("")}</tr>`).join("");
      return `<div class="lb-tabelle-rahmen"><table class="lb-tabelle"><thead><tr>${kopf}</tr></thead><tbody>${rumpf}</tbody></table></div>`;
    }
    if (/^\$\$[\s\S]*\$\$$/.test(t) && (t.match(/\$\$/g) || []).length === 2) return `<div class="lb-formel">${tex(t.slice(2, -2), true)}</div>`;
    return `<p>${mathe(t)}</p>`;
  }).join("");
}
function tex(s, gross) {
  try {
    const opt = { displayMode: gross, throwOnError: false, strict: false };
    const teile = gross ? teileAnQquad(s) : [];
    if (teile.length < 2) return katex.renderToString(s, opt);
    // Doppelformel „A,\qquad B“: nebeneinander, wenn Platz ist, sonst untereinander (formelnEinpassen entscheidet)
    return `<span class="tex-paar"><span class="tex-breit">${katex.renderToString(s, opt)}</span>` +
      `<span class="tex-stapel">${teile.map((t) => katex.renderToString(t, opt)).join("")}</span></span>`;
  }
  catch (e) { return esc(s); }
}
// an \qquad auf oberster Klammerebene trennen; gibt es keins, an „,\quad“ (Aufzählung mehrerer Formeln)
function teileAnQquad(s) {
  const teile = teileAnTrenner(s, "\\qquad", false);
  return teile.length > 1 ? teile : teileAnTrenner(s, "\\quad", true);
}
function teileAnTrenner(s, trenner, nachKomma) {
  const teile = [];
  let tiefe = 0, start = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === "\\") {
      const passt = tiefe === 0 && s.startsWith(trenner, i) && !/[a-zA-Z]/.test(s[i + trenner.length] || "") &&
        (!nachKomma || s.slice(start, i).trimEnd().endsWith(","));
      if (passt) { teile.push(s.slice(start, i)); start = i + trenner.length; i += trenner.length - 1; }
      else i++;
      continue;
    }
    if (c === "{") tiefe++;
    else if (c === "}") tiefe--;
  }
  teile.push(s.slice(start));
  return teile.map((t) => t.trim()).filter(Boolean);
}
function formelnEinpassen() {
  for (const p of document.querySelectorAll(".tex-paar")) {
    if (!p.offsetParent) continue;
    p.classList.remove("gestapelt");
    const d = p.querySelector(".tex-breit .katex-display");
    if (d && d.scrollWidth > d.clientWidth + 1) p.classList.add("gestapelt");
  }
}
{
  let geplant = false;
  const bald = () => { if (!geplant) { geplant = true; setTimeout(() => { geplant = false; formelnEinpassen(); }, 30); } };
  new MutationObserver(bald).observe(document.documentElement, { childList: true, subtree: true, attributes: true, attributeFilter: ["hidden"] });
  addEventListener("resize", bald);
}
function mathe(s) {
  if (s == null) return "";
  s = String(s);
  const re = /\$\$([\s\S]+?)\$\$|\$([^$]+?)\$/g;
  let out = "", i = 0, m;
  while ((m = re.exec(s))) {
    out += textTeil(s.slice(i, m.index));
    out += m[1] !== undefined ? tex(m[1], true) : tex(m[2], false);
    i = re.lastIndex;
  }
  return out + textTeil(s.slice(i));
}
const reinText = (s) => String(s ?? "").replace(/\$+/g, "").replace(/\\(mathrm|operatorname|mathbb|mathcal|mathbf)\{([^}]*)\}/g, "$2")
  .replace(/\\(eta|lambda|alpha|beta|theta|mu|pi|sigma|varepsilon|varphi|Phi|Gamma)\b/g, (m, w) => ({ eta: "η", lambda: "λ", alpha: "α", beta: "β", theta: "θ", mu: "µ", pi: "π", sigma: "σ", varepsilon: "ε", varphi: "φ", Phi: "Φ", Gamma: "Γ" }[w]))
  .replace(/\\([a-zA-Z]+)/g, "$1").replace(/[{}]/g, "").replace(/\s+/g, " ").trim();
const kurz = (s, n) => (s.length > n ? s.slice(0, n - 1).trimEnd() + "…" : s);
function umbruch(s, breite, max) {
  const w = s.split(/\s+/), z = [];
  let cur = "";
  for (const x of w) {
    if ((cur + " " + x).trim().length > breite && cur) { z.push(cur); cur = x; } else cur = (cur + " " + x).trim();
    if (z.length === max) break;
  }
  if (z.length < max && cur) z.push(cur);
  if (z.length === max && w.join(" ").length > z.join(" ").length) z[max - 1] = kurz(z[max - 1] + "…", breite);
  return z;
}
function hash(s) { let h = 2166136261; for (const c of s) h = Math.imul(h ^ c.charCodeAt(0), 16777619); return (h >>> 0) / 4294967296; }

// ------------------------------------------------------------------ SVG-Formen (Theorie, Legende)
function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs || {})) if (v != null) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}
function formPunkte(art, r) {
  const poly = (n, rr, dreh = 0) => Array.from({ length: n }, (_, i) => [rr * Math.cos(dreh + (2 * Math.PI * i) / n), rr * Math.sin(dreh + (2 * Math.PI * i) / n)]);
  switch (art) {
    case "konstruktion": return poly(6, r * 1.08);
    case "vermutung": return poly(4, r * 1.25, -Math.PI / 2);
    case "beobachtung": return poly(3, r * 1.25, -Math.PI / 2);
    case "korrektur": return poly(3, r * 1.25, Math.PI / 2);
    case "gescheitert": return poly(8, r * 1.02, Math.PI / 8);
    case "gegenbeispiel": return Array.from({ length: 8 }, (_, i) => {
      const a = -Math.PI / 2 + (Math.PI * i) / 4, rr = i % 2 ? r * 0.45 : r * 1.15;
      return [rr * Math.cos(a), rr * Math.sin(a)];
    });
    default: return null;
  }
}
function form(art, r, g) {
  if (art === "definition") return el("rect", { x: -r * 0.9, y: -r * 0.9, width: r * 1.8, height: r * 1.8, rx: r * 0.3, class: "form" }, g);
  if (art === "annahme") return el("rect", { x: -r * 1.3, y: -r * 0.75, width: r * 2.6, height: r * 1.5, rx: r * 0.75, class: "form" }, g);
  const p = formPunkte(art, r);
  if (p) return el("polygon", { points: p.map(([x, y]) => `${x.toFixed(1)},${y.toFixed(1)}`).join(" "), class: "form" }, g);
  return el("circle", { r, class: "form" }, g);
}
function stilFuer(k) {
  const st = k._status_eff, c = farbe(st);
  if (st === "offen") return { fill: "var(--flaeche)", stroke: c, schrift: c, strich: null };
  if (st === "bedingt") return { fill: `color-mix(in srgb, ${c} 30%, var(--flaeche))`, stroke: c, schrift: c, strich: "4 2.5" };
  if (st === "skizze") return { fill: `color-mix(in srgb, ${c} 22%, var(--flaeche))`, stroke: c, schrift: c, strich: "2 2.5" };
  if (st === "annahme") return { fill: `color-mix(in srgb, ${c} 28%, var(--flaeche))`, stroke: c, schrift: c, strich: null };
  if (st === "ueberholt") return { fill: "var(--flaeche-2)", stroke: c, schrift: c, strich: "1 2" };
  return { fill: c, stroke: c, schrift: "#fff", strich: null };
}
const radius = (k, faktor = 1) => Math.min(26, (8 + 3 * Math.sqrt(k._tragweite || 0)) * faktor);
