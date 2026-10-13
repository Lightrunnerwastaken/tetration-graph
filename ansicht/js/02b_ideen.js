// ===================================================================== Ideen-Landkarte
// Hauptideen → Unterideen → Aussagen, aufklappbar per Klick. Lage, Füllstand und Herkunft kommen fertig aus
// wissen/landkarte.py (G.ideen); hier wird nur gezeichnet und bedient.
const ideen = (() => {
  const D = G.ideen || { haupt: [], unter: [], kopplung_haupt: [], kopplung_unter: [], aussagen: {} };
  const leinwand = $("#ideen-flaeche");
  const ctx = leinwand.getContext("2d");
  const tipp = $("#tipp");
  const H = new Map(D.haupt.map((h) => [h.id, h]));
  const U = new Map(D.unter.map((u) => [u.id, u]));
  const A = D.aussagen;
  const PALETTE = ["#7fa8ff", "#3fd6b4", "#bb92ff", "#f2b84b", "#ff8f70", "#6fd0ff", "#e98fd0", "#a6c46a", "#c9d1e3"];
  const hFarbe = new Map(D.haupt.map((h, i) => [h.id, PALETTE[i % PALETTE.length]]));
  const GRUPPEN = [
    ["bewiesen", "#7fa8ff", "bewiesen"], ["bedingt", "#b9cdfa", "bedingt"], ["numerisch", "#3fd6b4", "numerisch"],
    ["vermutung", "#bb92ff", "Vermutungen"], ["offen", "#a3aec3", "offen"], ["verworfen", "#ff6a5c", "verworfen"],
  ];
  const HERKUNFT = [["literatur", "#93a6d4", "Literatur"], ["forum", "#e3ab3d", "Forum"], ["eigen", "#7fe3c4", "eigen"]];
  const STATUSFARBE = {
    definition: "#c9d1e3", bewiesen: "#7fa8ff", bewiesen_literatur: "#a6bbea", entlastet: "#7fa8ff", bedingt: "#6f9cff",
    skizze: "#9aa6bd", numerisch: "#3fd6b4", vermutung: "#bb92ff", annahme: "#f2b84b", offen: "#a3aec3",
    gescheitert: "#e3926f", widerlegt: "#ff6a5c", ueberholt: "#5d6880",
  };
  const rgba = (h, a) => { const n = parseInt(h.slice(1), 16); return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${Math.max(0, Math.min(1, a)).toFixed(3)})`; };
  const SERIF = `"KaTeX_Main", Cambria, Georgia, serif`;

  // Kopplungen nachschlagen
  const kopH = D.kopplung_haupt;
  const quer = new Map();          // unteridee -> Map(fremde Hauptidee -> Gewicht)
  const innen = [];                // [u, v, w] innerhalb einer Hauptidee
  for (const [a, b, w] of D.kopplung_unter) {
    const ea = U.get(a)?.eltern, eb = U.get(b)?.eltern;
    if (!ea || !eb) continue;
    if (ea === eb) { innen.push([a, b, w]); continue; }
    for (const [x, ey] of [[a, eb], [b, ea]]) {
      if (!quer.has(x)) quer.set(x, new Map());
      quer.get(x).set(ey, (quer.get(x).get(ey) || 0) + w);
    }
  }

  // ------------------------------------------------------------------ Kamera
  const cam = { x: 0, y: 0, k: 1 };
  let W = 0, Hh = 0, dpr = 1, flug = null, aktiv = false, laeuft = false, kFit = 1, erstesMal = true;
  const offen = { haupt: null, unter: null };
  let hover = null, chips = [];
  const zuSchirm = (x, y) => [(x - cam.x) * cam.k + W / 2, (y - cam.y) * cam.k + Hh / 2];
  const zuWelt = (sx, sy) => [(sx - W / 2) / cam.k + cam.x, (sy - Hh / 2) / cam.k + cam.y];

  function ziel() {                       // Rahmen der offenen Stufe in Weltkoordinaten
    if (offen.unter) { const u = U.get(offen.unter); return [u.x - u.r * 1.18, u.y - u.r * 1.18, u.x + u.r * 1.18, u.y + u.r * 1.18]; }
    if (offen.haupt) { const h = H.get(offen.haupt); return [h.x - h.r * 1.06, h.y - h.r * 1.06, h.x + h.r * 1.06, h.y + h.r * 1.06]; }
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const h of D.haupt) { x0 = Math.min(x0, h.x - h.r); y0 = Math.min(y0, h.y - h.r); x1 = Math.max(x1, h.x + h.r); y1 = Math.max(y1, h.y + h.r); }
    return [x0, y0, x1, y1];
  }
  function einpassK() {
    const [x0, y0, x1, y1] = ziel();
    return Math.max(1e-3, Math.min((W - 80) / Math.max(1, x1 - x0), (Hh - 110) / Math.max(1, y1 - y0)));
  }
  function groesse() {
    const r = leinwand.parentElement.getBoundingClientRect();
    if (!r.width || r.height < 40) return;          // Handy, Stufe „Text“: Karte zugeklappt
    const rel = kFit ? cam.k / kFit : 1;
    W = r.width; Hh = r.height;
    dpr = Math.min(2, window.devicePixelRatio || 1);
    leinwand.width = Math.round(W * dpr); leinwand.height = Math.round(Hh * dpr);
    kFit = einpassK();
    cam.k = kFit * rel;
    anstossen();
  }
  function fliege(x, y, k, ms = 760) {
    if (!BEWEGUNG) { Object.assign(cam, { x, y, k }); anstossen(); return; }
    flug = { von: { ...cam }, nach: { x, y, k }, t0: performance.now(), ms };
    anstossen();
  }
  function flugSchritt(now) {
    if (!flug) return;
    const t = Math.min(1, (now - flug.t0) / flug.ms);
    const e = t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
    const { von, nach } = flug;
    cam.k = Math.exp(Math.log(von.k) * (1 - e) + Math.log(nach.k) * e);
    cam.x = von.x + (nach.x - von.x) * e;
    cam.y = von.y + (nach.y - von.y) * e;
    if (t >= 1) flug = null;
  }
  function einpassen(sofort) {
    kFit = einpassK();
    const [x0, y0, x1, y1] = ziel();
    const x = (x0 + x1) / 2, y = (y0 + y1) / 2 + 12 / kFit;
    if (sofort) { Object.assign(cam, { x, y, k: kFit }); anstossen(); } else fliege(x, y, kFit);
  }

  // ------------------------------------------------------------------ Zeichnen
  function ringe(x, y, r, fuell, herkunft, alpha, dicke) {
    const n = Object.values(fuell).reduce((s, v) => s + v, 0);
    if (n) {
      let a0 = -Math.PI / 2;
      for (const [g, f] of GRUPPEN) {
        const v = fuell[g] || 0;
        if (!v) continue;
        const a1 = a0 + (2 * Math.PI * v) / n;
        ctx.beginPath(); ctx.arc(x, y, r + dicke / 2 + 1, a0 + 0.012, a1 - 0.012);
        ctx.strokeStyle = rgba(f, alpha); ctx.lineWidth = dicke;
        ctx.setLineDash(g === "offen" ? [3, 3] : []);
        ctx.stroke();
        a0 = a1;
      }
      ctx.setLineDash([]);
    }
    const m = Object.values(herkunft).reduce((s, v) => s + v, 0);
    if (m) {
      let a0 = -Math.PI / 2;
      for (const [g, f] of HERKUNFT) {
        const v = herkunft[g] || 0;
        if (!v) continue;
        const a1 = a0 + (2 * Math.PI * v) / m;
        ctx.beginPath(); ctx.arc(x, y, r - Math.max(2, dicke * 0.35), a0, a1);
        ctx.strokeStyle = rgba(f, alpha * 0.8); ctx.lineWidth = Math.max(1.5, dicke * 0.3);
        ctx.stroke();
        a0 = a1;
      }
    }
  }
  function zeilen(text, breite, max) {
    const w = String(text).split(/\s+/), z = [];
    let cur = "";
    for (const x of w) {
      if ((cur + " " + x).trim().length > breite && cur) { z.push(cur); cur = x; } else cur = (cur + " " + x).trim();
    }
    if (cur) z.push(cur);
    if (z.length > max) { z.length = max; z[max - 1] = kurz(z[max - 1] + "…", breite + 1); }
    return z;
  }
  function beschriftung(x, y, text, px, farbe, breite, max, gewicht = 600) {
    ctx.font = `${gewicht} ${px}px ${SERIF}`;
    ctx.fillStyle = farbe; ctx.textAlign = "center"; ctx.textBaseline = "middle";
    const z = zeilen(text, breite, max);
    z.forEach((t, i) => ctx.fillText(t, x, y + (i - (z.length - 1) / 2) * px * 1.18));
    return z.length;
  }
  function kreisIdee(o, farbe, alpha, gross, mitText) {
    const [sx, sy] = zuSchirm(o.x, o.y);
    const r = o.r * cam.k;
    if (!(r > 0.5)) return;
    if (sx + r < -40 || sx - r > W + 40 || sy + r < -40 || sy - r > Hh + 40) return;
    const gr = ctx.createRadialGradient(sx, sy, r * 0.1, sx, sy, r);
    gr.addColorStop(0, rgba(farbe, 0.16 * alpha)); gr.addColorStop(1, rgba(farbe, 0.05 * alpha));
    ctx.beginPath(); ctx.arc(sx, sy, r, 0, 2 * Math.PI); ctx.fillStyle = gr; ctx.fill();
    ctx.strokeStyle = rgba(farbe, (hover === o.id ? 0.85 : 0.4) * alpha); ctx.lineWidth = hover === o.id ? 2 : 1; ctx.stroke();
    ringe(sx, sy, r, o.fuell, o.herkunft, alpha, Math.max(3, Math.min(gross ? 10 : 7, r * 0.06)));
    if (!mitText || r < 22) return;
    const px = gross ? Math.max(13, Math.min(20, r * 0.13)) : Math.max(11, Math.min(15, r * 0.16));
    const breite = Math.max(8, Math.floor((r * 1.6) / (px * 0.52)));
    const zz = beschriftung(sx, sy - px * 0.4, o.name, px, rgba("#e8eef9", alpha), breite, 3);
    ctx.font = `400 ${Math.round(px * 0.72)}px ${SERIF}`;
    ctx.fillStyle = rgba("#a9b7cf", alpha);
    ctx.fillText(`${o.n} Aussage${o.n === 1 ? "" : "n"}`, sx, sy - px * 0.4 + ((zz + 1) / 2) * px * 1.18 + 4);
  }
  function linie(x0, y0, x1, y1, breite, farbe, alpha) {
    const [a, b] = zuSchirm(x0, y0), [c, d] = zuSchirm(x1, y1);
    ctx.beginPath(); ctx.moveTo(a, b); ctx.lineTo(c, d);
    ctx.strokeStyle = rgba(farbe, alpha); ctx.lineWidth = breite; ctx.stroke();
  }
  function randpunkt(h) {                  // Randchip-Position einer fremden Hauptidee
    const [sx, sy] = zuSchirm(h.x, h.y);
    const cx = W / 2, cy = Hh / 2, dx = sx - cx, dy = sy - cy;
    const m = 70, s = Math.min((W / 2 - m) / Math.abs(dx || 1e-9), (Hh / 2 - m) / Math.abs(dy || 1e-9));
    return s < 1 ? [cx + dx * s, cy + dy * s] : [sx, sy];
  }
  // Handy: die anderen Hauptideen als wischbare Leiste unten statt als Schildchen am Rand
  const leiste = $("#ideen-chips");
  let leisteSchluessel = "";
  function leisteSetzen(ids) {
    const s = ids.join(",");
    if (s === leisteSchluessel) return;
    leisteSchluessel = s;
    leiste.innerHTML = ids.map((id) => `<button data-idee="${id}" style="--f:${hFarbe.get(id)}">${esc(H.get(id).name)}</button>`).join("");
    leiste.hidden = !ids.length;
  }
  const schmal = () => W < 600;

  function zeichneChips() {
    chips = [];
    if (schmal()) { leisteSetzen(offen.haupt ? D.haupt.filter((h) => h.id !== offen.haupt).map((h) => h.id) : []); return; }
    leisteSetzen([]);
    if (!offen.haupt) return;
    ctx.font = `600 12px ${SERIF}`;
    for (const h of D.haupt) {
      if (h.id === offen.haupt) continue;
      const [sx, sy] = zuSchirm(h.x, h.y);
      if (sx > 60 && sx < W - 60 && sy > 60 && sy < Hh - 60) continue;
      const [px, py] = randpunkt(h);
      const t = kurz(h.name, 30), b = ctx.measureText(t).width + 22;
      const x = Math.max(8, Math.min(W - b - 8, px - b / 2));
      let y = Math.max(46, Math.min(Hh - 34, py - 12));
      const schritt = py > Hh / 2 ? -30 : 30;          // nicht auf andere Chips legen
      for (let n = 0; n < 12 && chips.some((c) => x < c.x + c.b + 6 && c.x < x + b + 6 && Math.abs(c.y - y) < 28); n++) y += schritt;
      const f = hFarbe.get(h.id);
      ctx.beginPath(); ctx.roundRect(x, y, b, 24, 12);
      ctx.fillStyle = hover === "chip:" + h.id ? "rgba(30,44,72,.95)" : "rgba(13,21,37,.88)"; ctx.fill();
      ctx.strokeStyle = rgba(f, 0.7); ctx.lineWidth = 1; ctx.stroke();
      ctx.fillStyle = "#dbe4f3"; ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.fillText(t, x + 11, y + 12.5);
      chips.push({ id: h.id, x, y, b, mx: x + b / 2, my: y + 12 });
    }
  }
  function zeichneAussagen(u, alpha) {
    const ids = u.mitglieder;
    const drin = new Set(ids);
    ctx.lineWidth = 1;
    for (const i of ids) {
      const k = K.get(i), a = A[i];
      for (const d of k._abhaengig || []) {
        if (!drin.has(d)) continue;
        linie(a.x, a.y, A[d].x, A[d].y, 1, "#9db8ff", 0.22 * alpha);
      }
    }
    for (const i of ids) {
      const k = K.get(i), a = A[i];
      const [sx, sy] = zuSchirm(a.x, a.y);
      const r = Math.max(a.haupt ? 6 : 3.5, a.r * cam.k * 0.55);
      const f = STATUSFARBE[k._status_eff] || "#a3aec3";
      const hohl = k._status_eff === "offen" || k._status_eff === "ueberholt";
      const al = (a.haupt ? 1 : 0.55) * alpha;
      ctx.beginPath(); ctx.arc(sx, sy, r, 0, 2 * Math.PI);
      if (hohl) { ctx.strokeStyle = rgba(f, al); ctx.lineWidth = 1.5; ctx.stroke(); }
      else { ctx.fillStyle = rgba(f, al); ctx.fill(); }
      if (zustand.gewaehlt === i || hover === i) {
        ctx.beginPath(); ctx.arc(sx, sy, r + 4, 0, 2 * Math.PI);
        ctx.strokeStyle = zustand.gewaehlt === i ? "#ffffff" : "rgba(255,255,255,.6)"; ctx.lineWidth = 1.5; ctx.stroke();
      }
      if (alpha > 0.5 && (a.haupt || hover === i || zustand.gewaehlt === i)) {
        ctx.font = `${a.haupt ? 600 : 400} ${a.haupt ? 13 : 12}px ${SERIF}`;
        ctx.fillStyle = a.haupt ? "#e8eef9" : "#a9b7cf"; ctx.textAlign = "left"; ctx.textBaseline = "middle";
        ctx.fillText(kurz(reinText(k.titel), 40), sx + r + 6, sy);
      }
    }
  }
  function zeichne() {
    if (Hh < 40) return;
    const c = ctx;
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    const bg = c.createRadialGradient(W / 2, Hh / 2, 0, W / 2, Hh / 2, Math.max(W, Hh) * 0.75);
    bg.addColorStop(0, "#0d1628"); bg.addColorStop(1, "#060a13");
    c.fillStyle = bg; c.fillRect(0, 0, W, Hh);
    const h0 = offen.haupt && H.get(offen.haupt), u0 = offen.unter && U.get(offen.unter);
    if (!h0) leisteSetzen([]);
    if (!h0) {
      for (const [a, b, w] of kopH) {
        const p = H.get(a), q = H.get(b);
        linie(p.x, p.y, q.x, q.y, Math.min(8, 0.8 + 1.1 * Math.sqrt(w)), "#8f9dbb", 0.18 + Math.min(0.35, w / 120));
      }
      for (const h of D.haupt) kreisIdee(h, hFarbe.get(h.id), 1, true, true);
    } else {
      for (const h of D.haupt) if (h.id !== h0.id) kreisIdee(h, hFarbe.get(h.id), 0.22, true, false);
      kreisIdee({ ...h0, fuell: {}, herkunft: {} }, hFarbe.get(h0.id), 0.55, true, false);
      // Kopplungen nach aussen: von der Unteridee zur fremden Hauptidee
      if (!u0) {
        for (const id of h0.kinder) {
          const u = U.get(id), m = quer.get(id);
          if (!m) continue;
          for (const [eh, w] of m) {
            const ziel = H.get(eh), [rx, ry] = randpunkt(ziel), [wx, wy] = zuWelt(rx, ry);
            linie(u.x, u.y, wx, wy, Math.min(5, 0.6 + 0.8 * Math.sqrt(w)), hFarbe.get(eh), 0.16 + Math.min(0.3, w / 25));
          }
        }
        for (const [a, b, w] of innen) {
          if (U.get(a).eltern !== h0.id) continue;
          const p = U.get(a), q = U.get(b);
          linie(p.x, p.y, q.x, q.y, Math.min(6, 0.8 + Math.sqrt(w)), "#c9d1e3", 0.15 + Math.min(0.35, w / 20));
        }
      }
      for (const id of h0.kinder) {
        const u = U.get(id);
        if (u0 && id !== u0.id) kreisIdee(u, hFarbe.get(h0.id), 0.3, false, true);
        else if (!u0) kreisIdee(u, hFarbe.get(h0.id), 1, false, true);
      }
      if (u0) {
        const [sx, sy] = zuSchirm(u0.x, u0.y), r = u0.r * cam.k;
        c.beginPath(); c.arc(sx, sy, r, 0, 2 * Math.PI);
        c.fillStyle = rgba(hFarbe.get(h0.id), 0.05); c.fill();
        c.strokeStyle = rgba(hFarbe.get(h0.id), 0.45); c.lineWidth = 1; c.stroke();
        ringe(sx, sy, r, u0.fuell, u0.herkunft, 1, Math.max(3, Math.min(9, r * 0.03)));
        c.font = `600 15px ${SERIF}`; c.fillStyle = "#e8eef9"; c.textAlign = "center"; c.textBaseline = "bottom";
        c.fillText(u0.name, sx, sy - r - 14);
        zeichneAussagen(u0, 1);
      }
      zeichneChips();
    }
  }

  // ------------------------------------------------------------------ Treffer
  function treffer(sx, sy) {
    for (const ch of chips) if (sx >= ch.x && sx <= ch.x + ch.b && sy >= ch.y && sy <= ch.y + 24) return { art: "haupt", id: ch.id, chip: true };
    const [wx, wy] = zuWelt(sx, sy);
    const drin = (o) => Math.hypot(wx - o.x, wy - o.y) <= o.r;
    if (offen.unter) {
      let best = null, bd = Infinity;
      for (const i of U.get(offen.unter).mitglieder) {
        const a = A[i], [ax, ay] = zuSchirm(a.x, a.y);
        const d = Math.hypot(sx - ax, sy - ay), r = Math.max(a.haupt ? 6 : 3.5, a.r * cam.k * 0.55) + 5;
        if (d <= r && d < bd) { bd = d; best = i; }
      }
      if (best) return { art: "aussage", id: best };
    }
    if (offen.haupt) {
      for (const id of H.get(offen.haupt).kinder) if (id !== offen.unter && drin(U.get(id))) return { art: "unter", id };
      if (offen.unter && drin(U.get(offen.unter))) return null;
      if (drin(H.get(offen.haupt))) return null;
    }
    for (const h of D.haupt) if (h.id !== offen.haupt && drin(h)) return { art: "haupt", id: h.id };
    return null;
  }
  function standText(o) {
    const t = GRUPPEN.filter(([g]) => o.fuell[g]).map(([g, , n]) => `${o.fuell[g]} ${n}`).join(", ");
    const q = HERKUNFT.filter(([g]) => o.herkunft[g]).map(([g, , n]) => `${n} ${o.herkunft[g]}`).join(" · ");
    return `${t}<br><span class="tipp-status">Quellen: ${q}</span>`;
  }
  function tippZeigen(sx, sy, html) {
    tipp.innerHTML = html; tipp.hidden = false;
    let x = sx + 14, y = sy + 12;
    if (x + 300 > W) x = sx - 310;
    if (y + 120 > Hh) y = sy - 120;
    tipp.style.left = x + "px"; tipp.style.top = y + "px";
  }

  // ------------------------------------------------------------------ Bedienung
  let zug = null;
  leinwand.addEventListener("pointerdown", (e) => {
    zug = { x: e.clientX, y: e.clientY, cx: cam.x, cy: cam.y, bewegt: false };
    leinwand.setPointerCapture(e.pointerId);
  });
  leinwand.addEventListener("pointermove", (e) => {
    const r = leinwand.getBoundingClientRect(), sx = e.clientX - r.left, sy = e.clientY - r.top;
    if (zug) {
      const dx = e.clientX - zug.x, dy = e.clientY - zug.y;
      if (Math.abs(dx) + Math.abs(dy) > 4) zug.bewegt = true;
      if (zug.bewegt) { flug = null; cam.x = zug.cx - dx / cam.k; cam.y = zug.cy - dy / cam.k; leinwand.classList.add("zieht"); tipp.hidden = true; anstossen(); return; }
    }
    const t = treffer(sx, sy);
    const neu = t ? (t.chip ? "chip:" + t.id : t.id) : null;
    if (neu !== hover) { hover = neu; anstossen(); }
    leinwand.style.cursor = t ? "pointer" : "grab";
    if (!t) { tipp.hidden = true; return; }
    if (t.art === "aussage") {
      const k = K.get(t.id);
      tippZeigen(sx, sy, `<b>${esc(kurz(reinText(k.titel), 80))}</b><br><span class="tipp-status">${esc(S.status[k._status_eff]?.label || k._status_eff)} · ${esc(S.arten[k.art] || k.art)}${A[t.id].haupt ? " · Hauptaussage" : ""}</span>`);
    } else {
      const o = t.art === "haupt" ? H.get(t.id) : U.get(t.id);
      tippZeigen(sx, sy, `<b>${esc(o.name)}</b><br>${esc(o.kurz)}<br>${standText(o)}`);
    }
  });
  function loslassen(e) {
    if (!zug) return;
    const war = zug; zug = null; leinwand.classList.remove("zieht");
    if (war.bewegt) return;
    const r = leinwand.getBoundingClientRect();
    const t = treffer(e.clientX - r.left, e.clientY - r.top);
    if (!t) return;
    tipp.hidden = true;
    if (t.art === "aussage") waehle(t.id);
    else oeffne(t.id);
  }
  leinwand.addEventListener("pointerup", loslassen);
  leinwand.addEventListener("pointercancel", () => { zug = null; leinwand.classList.remove("zieht"); });
  leinwand.addEventListener("pointerleave", () => { if (hover) { hover = null; anstossen(); } tipp.hidden = true; });
  leinwand.addEventListener("wheel", (e) => {
    e.preventDefault();
    const r = leinwand.getBoundingClientRect(), sx = e.clientX - r.left, sy = e.clientY - r.top;
    const [wx, wy] = zuWelt(sx, sy);
    const k = Math.max(kFit * 0.6, Math.min(kFit * 4, cam.k * Math.exp(-e.deltaY * 0.0015)));
    flug = null;
    cam.k = k; cam.x = wx - (sx - W / 2) / k; cam.y = wy - (sy - Hh / 2) / k;
    anstossen();
  }, { passive: false });
  new ResizeObserver(() => { if (aktiv) groesse(); }).observe(leinwand.parentElement);

  function oeffne(id) {
    if (!id) { offen.haupt = null; offen.unter = null; }
    else if (H.has(id)) { offen.haupt = id; offen.unter = null; }
    else if (U.has(id)) { offen.haupt = U.get(id).eltern; offen.unter = id; }
    hover = null;
    einpassen(false);
    if (typeof pfadZeigen === "function") pfadZeigen();
  }
  function hoch() {
    if (zustand.gewaehlt && offen.unter) { waehle(null); return true; }
    if (offen.unter) { oeffne(offen.haupt); return true; }
    if (offen.haupt) { oeffne(null); return true; }
    return false;
  }

  // ------------------------------------------------------------------ Schleife
  function anstossen() { if (!laeuft && aktiv) { laeuft = true; requestAnimationFrame(bild); } }
  function bild(now) {
    laeuft = false;
    if (!aktiv) return;
    flugSchritt(now);
    zeichne();
    if (flug) { laeuft = true; requestAnimationFrame(bild); }
  }

  return {
    zeige() {
      aktiv = true; leinwand.hidden = false;
      groesse();
      if (erstesMal) { erstesMal = false; einpassen(true); }
      document.fonts?.load(`600 20px "KaTeX_Main"`).then(anstossen);
      anstossen();
    },
    verberge() { aktiv = false; leinwand.hidden = true; tipp.hidden = true; leiste.hidden = true; leisteSchluessel = ""; },
    oeffne, hoch, anstossen,
    zuKnoten(id) { const a = A[id]; if (a && a.idee !== offen.unter) oeffne(a.idee); else anstossen(); },
    zoom(f) { const k = Math.max(kFit * 0.6, Math.min(kFit * 4, cam.k * f)); fliege(cam.x, cam.y, k, 380); },
    pfadHtml() {
      const t = [`<button data-idee="_alle">Ideen</button>`];
      if (offen.haupt) t.push(`<button data-idee="${offen.haupt}">${esc(H.get(offen.haupt).name)}</button>`);
      if (offen.unter) t.push(`<button data-idee="${offen.unter}">${esc(U.get(offen.unter).name)}</button>`);
      const k = zustand.gewaehlt && K.get(zustand.gewaehlt);
      if (k && A[k.id]?.idee === offen.unter) t.push(`<span>${esc(kurz(reinText(k.titel), 48))}</span>`);
      return t.join(" › ");
    },
    get aktiv() { return aktiv; },
  };
})();
