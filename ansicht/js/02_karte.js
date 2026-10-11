// ===================================================================== Karte
// Eine zoombare Karte auf Canvas. Stufen: Gebiete → Verläufe (Stränge) → Aussagen → Inhalt.
// Was sichtbar ist, hängt stetig vom Zoom z = k / k0 ab (k0 = ganze Karte passt).
const karte = (() => {
  const leinwand = $("#karte");
  const ctx = leinwand.getContext("2d");
  const kartenEbene = $("#karten");
  const tipp = $("#tipp");

  const KF = {
    bg1: "#0d1628", bg2: "#060a13", raster: "rgba(130,170,240,0.045)", rasterStark: "rgba(130,170,240,0.085)",
    text: "#e8eef9", text2: "#a9b7cf", text3: "#71819e",
    status: {
      definition: "#c9d1e3", bewiesen: "#7fa8ff", bewiesen_literatur: "#a6bbea", entlastet: "#7fa8ff", bedingt: "#6f9cff",
      skizze: "#9aa6bd", numerisch: "#3fd6b4", vermutung: "#bb92ff", annahme: "#f2b84b", offen: "#a3aec3",
      gescheitert: "#e3926f", widerlegt: "#ff6a5c", ueberholt: "#5d6880",
    },
    kante: { beweis: "#9db8ff", rot: "#ff6a5c", gruen: "#3fd6b4", evidenz: "#3fd6b4", struktur: "#8f9dbb", leicht: "#7d8aa6" },
  };
  const HOHL = new Set(["offen", "ueberholt"]);
  const STUFEN = [{ name: "Gebiete", z: 1 }, { name: "Verläufe", z: 2.5 }, { name: "Aussagen", z: 4.8 }, { name: "Inhalt", z: 8.6 }];
  const stufeVon = (z) => (z < 1.75 ? 0 : z < 3.6 ? 1 : z < 6.6 ? 2 : 3);

  const rgb = (h) => { const n = parseInt(h.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; };
  const rgba = (h, a) => { const [r, g, b] = rgb(h); return `rgba(${r},${g},${b},${Math.max(0, Math.min(1, a)).toFixed(3)})`; };
  const aufhellen = (h, f) => "#" + rgb(h).map((c) => Math.round(c + (255 - c) * f).toString(16).padStart(2, "0")).join("");
  const glatt = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
  const sprites = new Map();
  function glimmen(h) {
    if (sprites.has(h)) return sprites.get(h);
    const c = document.createElement("canvas");
    c.width = c.height = 128;
    const g = c.getContext("2d");
    const gr = g.createRadialGradient(64, 64, 0, 64, 64, 64);
    gr.addColorStop(0, rgba(h, 0.9)); gr.addColorStop(0.22, rgba(h, 0.38)); gr.addColorStop(0.55, rgba(h, 0.09)); gr.addColorStop(1, rgba(h, 0));
    g.fillStyle = gr;
    g.fillRect(0, 0, 128, 128);
    sprites.set(h, c);
    return c;
  }

  // ------------------------------------------------------------------ Layout
  const L = { knoten: new Map(), straenge: new Map(), gebiete: new Map(), kanten: [], baender: [], strangVon: new Map() };
  const rWelt = (id) => 7 + 2.6 * Math.sqrt(K.get(id)._tragweite || 0);

  function huelle(pkt) {
    if (pkt.length < 3) return pkt.slice();
    const p = pkt.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const kreuz = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const u = [], l = [];
    for (const q of p) { while (l.length >= 2 && kreuz(l[l.length - 2], l[l.length - 1], q) <= 0) l.pop(); l.push(q); }
    for (const q of p.slice().reverse()) { while (u.length >= 2 && kreuz(u[u.length - 2], u[u.length - 1], q) <= 0) u.pop(); u.push(q); }
    return l.slice(0, -1).concat(u.slice(0, -1));
  }

  function layout() {
    for (const s of G.straenge || []) L.straenge.set(s.id, { ...s, glieder: [] });
    for (const k of G.knoten) {
      let sid = k._strang;
      if (!sid || !L.straenge.has(sid)) {
        sid = k.gebiet + "._weitere";
        if (!L.straenge.has(sid)) L.straenge.set(sid, { id: sid, gebiet: k.gebiet, name: "Weitere", verlauf: "", glieder: [], pseudo: true });
      }
      L.straenge.get(sid).glieder.push(k.id);
      L.strangVon.set(k.id, sid);
    }
    // Gebiete: Startkreis, dann Abstossung nach Grösse und Anziehung nach Querverbindungen
    const geb = G.gebiete.map((g, i) => {
      const n = G.knoten.filter((k) => k.gebiet === g.id).length, a = -Math.PI / 2 + (2 * Math.PI * i) / G.gebiete.length;
      return { id: g.id, name: g.name, farbe: aufhellen(g.farbe, 0.38), n, R: 95 + 44 * Math.sqrt(n), x: Math.cos(a) * 1150, y: Math.sin(a) * 800 };
    });
    const gi = new Map(geb.map((g, i) => [g.id, i]));
    const W = geb.map(() => geb.map(() => 0));
    for (const e of G.kanten) {
      const a = gi.get(K.get(e.von).gebiet), b = gi.get(K.get(e.nach).gebiet);
      if (a !== b) { W[a][b]++; W[b][a]++; }
    }
    for (let it = 0; it < 600; it++) {
      for (let i = 0; i < geb.length; i++) for (let j = i + 1; j < geb.length; j++) {
        const A = geb[i], B = geb[j], dx = B.x - A.x, dy = B.y - A.y, d = Math.hypot(dx, dy) + 0.01;
        const min = A.R + B.R + 190;
        let f = 0;
        if (d < min) f = -(min - d) * 0.5;
        else if (W[i][j]) f = (d - min) * 0.006 * Math.sqrt(W[i][j]);
        A.x += (dx / d) * f * 0.5; A.y += (dy / d) * f * 0.5; B.x -= (dx / d) * f * 0.5; B.y -= (dy / d) * f * 0.5;
      }
      for (const g of geb) { g.x *= 0.996; g.y *= 0.996; }
    }
    // Stränge innerhalb der Gebiete
    for (const g of geb) {
      const st = [...L.straenge.values()].filter((s) => s.gebiet === g.id);
      st.forEach((s, i) => {
        s.r = 42 + 30 * Math.sqrt(s.glieder.length);
        const a = (2 * Math.PI * i) / st.length + hash(g.id) * 6;
        const rr = st.length > 1 ? g.R * 0.5 : 0;
        s.x = g.x + Math.cos(a) * rr; s.y = g.y + Math.sin(a) * rr;
      });
      for (let it = 0; it < 300; it++) {
        for (let i = 0; i < st.length; i++) for (let j = i + 1; j < st.length; j++) {
          const A = st[i], B = st[j], dx = B.x - A.x, dy = B.y - A.y, d = Math.hypot(dx, dy) + 0.01, min = A.r + B.r + 46;
          if (d < min) { const f = (min - d) * 0.5; A.x -= (dx / d) * f * 0.5; A.y -= (dy / d) * f * 0.5; B.x += (dx / d) * f * 0.5; B.y += (dy / d) * f * 0.5; }
        }
        for (const s of st) { s.x += (g.x - s.x) * 0.03; s.y += (g.y - s.y) * 0.03; }
      }
    }
    // Aussagen: gemeinsame Simulation mit Federn, Strangbindung und Kollision
    const ids = G.knoten.map((k) => k.id);
    const P = new Map();
    for (const s of L.straenge.values()) {
      s.glieder.forEach((id, j) => {
        const a = j * 2.39996 + hash(id) * 2, rr = 24 * Math.sqrt(j + 0.6);
        P.set(id, { x: s.x + Math.cos(a) * rr, y: s.y + Math.sin(a) * rr });
      });
    }
    const R = new Map(ids.map((id) => [id, rWelt(id)]));
    const E = G.kanten.filter((e) => P.has(e.von) && P.has(e.nach)).map((e) => {
      const sa = L.strangVon.get(e.von), sb = L.strangVon.get(e.nach);
      const gleichGebiet = K.get(e.von).gebiet === K.get(e.nach).gebiet;
      return [e.von, e.nach, sa === sb ? 74 : gleichGebiet ? 150 : 340, sa === sb ? 0.05 : gleichGebiet ? 0.018 : 0.004];
    });
    const SCHRITTE = 420;
    for (let s = 0; s < SCHRITTE; s++) {
      const temp = 34 * (1 - s / SCHRITTE) + 0.4;
      const D = new Map(ids.map((id) => [id, [0, 0]]));
      for (let i = 0; i < ids.length; i++) {
        const a = ids[i], pa = P.get(a);
        for (let j = i + 1; j < ids.length; j++) {
          const b = ids[j], pb = P.get(b);
          const dx = pa.x - pb.x, dy = pa.y - pb.y, d = Math.hypot(dx, dy) + 0.01;
          if (d > 230) continue;
          let f = 2400 / (d * d) * 8;
          const min = R.get(a) + R.get(b) + 26;
          if (d < min) f += (min - d) * 0.9;
          D.get(a)[0] += (dx / d) * f; D.get(a)[1] += (dy / d) * f;
          D.get(b)[0] -= (dx / d) * f; D.get(b)[1] -= (dy / d) * f;
        }
      }
      for (const [a, b, l0, c] of E) {
        const pa = P.get(a), pb = P.get(b), dx = pb.x - pa.x, dy = pb.y - pa.y, d = Math.hypot(dx, dy) + 0.01;
        const f = (d - l0) * c;
        D.get(a)[0] += (dx / d) * f; D.get(a)[1] += (dy / d) * f;
        D.get(b)[0] -= (dx / d) * f; D.get(b)[1] -= (dy / d) * f;
      }
      for (const id of ids) {
        const p = P.get(id), st = L.straenge.get(L.strangVon.get(id)), d = D.get(id);
        d[0] += (st.x - p.x) * 0.07; d[1] += (st.y - p.y) * 0.07;
        const m = Math.hypot(d[0], d[1]) + 1e-9, sch = Math.min(m, temp);
        p.x += (d[0] / m) * sch; p.y += (d[1] / m) * sch;
      }
    }
    ids.forEach((id) => L.knoten.set(id, { id, x: P.get(id).x, y: P.get(id).y, r: R.get(id), k: K.get(id) }));
    // Stränge: Hülle und Anker aus den tatsächlichen Positionen
    for (const s of L.straenge.values()) {
      const pkt = s.glieder.map((id) => [L.knoten.get(id).x, L.knoten.get(id).y]);
      s.x = pkt.reduce((u, p) => u + p[0], 0) / pkt.length;
      s.y = pkt.reduce((u, p) => u + p[1], 0) / pkt.length;
      s.huelle = huelle(pkt);
      s.rad = Math.max(30, ...pkt.map((p) => Math.hypot(p[0] - s.x, p[1] - s.y))) + 30;
      s.oben = Math.min(...pkt.map((p) => p[1])) - 34;
    }
    for (const g of geb) {
      const pkt = G.knoten.filter((k) => k.gebiet === g.id).map((k) => L.knoten.get(k.id));
      g.x = pkt.reduce((u, p) => u + p.x, 0) / pkt.length;
      g.y = pkt.reduce((u, p) => u + p.y, 0) / pkt.length;
      g.R = Math.max(...pkt.map((p) => Math.hypot(p.x - g.x, p.y - g.y) + p.r)) + 60;
      const z = {};
      pkt.forEach((p) => (z[p.k._status_eff] = (z[p.k._status_eff] || 0) + 1));
      g.ring = RINGFOLGE.filter((s) => z[s]).map((s) => [s, z[s] / pkt.length]);
      g.offen = pkt.filter((p) => ["vermutung", "offen", "annahme", "skizze"].includes(p.k._status_eff)).length;
      g.neu = pkt.filter((p) => p.k._neu_pruefen?.length).length;
      L.gebiete.set(g.id, g);
    }
    // Gebiete starr auseinanderschieben, bis sich ihre Ringe nicht berühren
    const verschiebe = (g, dx, dy) => {
      g.x += dx; g.y += dy;
      for (const k of G.knoten) if (k.gebiet === g.id) { const n = L.knoten.get(k.id); n.x += dx; n.y += dy; }
      for (const st of L.straenge.values()) if (st.gebiet === g.id) {
        st.x += dx; st.y += dy; st.oben += dy; st.huelle = st.huelle.map(([x, y]) => [x + dx, y + dy]);
      }
    };
    for (let it = 0; it < 400; it++) {
      let bewegt = false;
      for (let i = 0; i < geb.length; i++) for (let j = i + 1; j < geb.length; j++) {
        const A = geb[i], B = geb[j], dx = B.x - A.x, dy = B.y - A.y, d = Math.hypot(dx, dy) + 0.01, min = A.R + B.R + 70;
        if (d >= min) continue;
        const f = (min - d) / 2 + 1;
        verschiebe(A, (-dx / d) * f, (-dy / d) * f);
        verschiebe(B, (dx / d) * f, (dy / d) * f);
        bewegt = true;
      }
      // sanft zusammenhalten, damit nichts davontreibt
      const mx = geb.reduce((u, g) => u + g.x, 0) / geb.length, my = geb.reduce((u, g) => u + g.y, 0) / geb.length;
      for (const g of geb) verschiebe(g, (mx - g.x) * 0.01, (my - g.y) * 0.01);
      if (!bewegt && it > 20) break;
    }
    // Kanten als leicht gebogene Kurven; Beweise fliessen von der Voraussetzung zur Folgerung
    for (const e of G.kanten) {
      const kl = KANTENKLASSE[e.typ] || "leicht";
      const [va, vb] = kl === "beweis" ? [e.nach, e.von] : [e.von, e.nach];
      const a = L.knoten.get(va), b = L.knoten.get(vb);
      if (!a || !b) continue;
      const dx = b.x - a.x, dy = b.y - a.y, len = Math.hypot(dx, dy) + 0.01;
      const bieg = (hash(e.von + e.nach) - 0.5) * 0.28 * len;
      const c = [(a.x + b.x) / 2 - (dy / len) * bieg, (a.y + b.y) / 2 + (dx / len) * bieg];
      L.kanten.push({ e, kl, a, b, c, len, phase: hash(e.nach + e.von), quer: K.get(e.von).gebiet !== K.get(e.nach).gebiet });
    }
    for (let i = 0; i < geb.length; i++) for (let j = i + 1; j < geb.length; j++) {
      if (W[i][j]) L.baender.push({ a: geb[i], b: geb[j], w: W[i][j], phase: hash(geb[i].id + geb[j].id) });
    }
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const g of geb) { x0 = Math.min(x0, g.x - g.R); y0 = Math.min(y0, g.y - g.R); x1 = Math.max(x1, g.x + g.R); y1 = Math.max(y1, g.y + g.R); }
    L.box = [x0, y0, x1, y1];
  }
  const qPunkt = (kt, t) => {
    const u = 1 - t;
    return [u * u * kt.a.x + 2 * u * t * kt.c[0] + t * t * kt.b.x, u * u * kt.a.y + 2 * u * t * kt.c[1] + t * t * kt.b.y];
  };

  // ------------------------------------------------------------------ Kamera
  const cam = { x: 0, y: 0, k: 1 };
  let k0 = 1, W = 0, H = 0, dpr = 1, flug = null, aktiv = false, laeuft = false, erstesMal = false;
  const z = () => cam.k / k0;
  const zuSchirm = (x, y) => [(x - cam.x) * cam.k + W / 2, (y - cam.y) * cam.k + H / 2];
  const zuWelt = (sx, sy) => [(sx - W / 2) / cam.k + cam.x, (sy - H / 2) / cam.k + cam.y];
  function groesse() {
    const r = leinwand.parentElement.getBoundingClientRect();
    if (!r.width) return;
    W = r.width; H = r.height;
    dpr = Math.min(2, window.devicePixelRatio || 1);
    leinwand.width = Math.round(W * dpr); leinwand.height = Math.round(H * dpr);
    const [x0, y0, x1, y1] = L.box;
    const altZ = k0 ? z() : 1;
    k0 = Math.min((W - 60) / (x1 - x0), (H - 110) / (y1 - y0));
    cam.k = k0 * altZ;
    anstossen();
  }
  function fliege(x, y, k, ms = 820) {
    k = Math.max(k0 * 0.6, Math.min(k0 * 16, k));
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
  function rahmen(x0, y0, x1, y1, zMin, zMax) {
    groesse();                                  // Tafel kann gerade auf- oder zugegangen sein
    const k = Math.min((W - 160) / Math.max(1, x1 - x0), (H - 170) / Math.max(1, y1 - y0));
    fliege((x0 + x1) / 2, (y0 + y1) / 2, Math.max(k0 * zMin, Math.min(k0 * zMax, k)));
  }
  function ganz() { const [x0, y0, x1, y1] = L.box; fliege((x0 + x1) / 2, (y0 + y1) / 2 + 20, k0); }

  // ------------------------------------------------------------------ Fokus
  let fokus = null, zweiter = new Set(), hover = null;
  function setzeFokus(id) {
    fokus = id;
    zweiter = new Set();
    if (id) for (const n of NACHBARN.get(id)) for (const m of NACHBARN.get(n)) if (m !== id && !NACHBARN.get(id).has(m)) zweiter.add(m);
    anstossen();
  }
  function knotenAlpha(id) {
    if (!sichtbar(K.get(id))) return 0.07;
    if (!fokus || id === fokus || NACHBARN.get(fokus).has(id)) return 1;
    return zweiter.has(id) ? 0.4 : 0.09;
  }
  function zuKnoten(id) {
    const n = L.knoten.get(id);
    if (!n) return;
    const gr = [id, ...NACHBARN.get(id)].map((i) => L.knoten.get(i)).filter(Boolean);
    const x0 = Math.min(...gr.map((p) => p.x - p.r)), x1 = Math.max(...gr.map((p) => p.x + p.r));
    const y0 = Math.min(...gr.map((p) => p.y - p.r)), y1 = Math.max(...gr.map((p) => p.y + p.r + 30));
    // Knoten selbst zur Mitte: Rahmen symmetrisch um ihn legen
    const hx = Math.max(n.x - x0, x1 - n.x), hy = Math.max(n.y - y0, y1 - n.y);
    rahmen(n.x - hx, n.y - hy, n.x + hx, n.y + hy, 4.2, 9.5);
  }
  function zuStrang(sid) {
    const s = L.straenge.get(sid);
    const p = s.glieder.map((i) => L.knoten.get(i));
    rahmen(Math.min(...p.map((q) => q.x)) - 60, s.oben - 80, Math.max(...p.map((q) => q.x)) + 60, Math.max(...p.map((q) => q.y)) + 60, 2.6, 6.2);
  }
  function zuGebiet(gid) {
    const g = L.gebiete.get(gid);
    rahmen(g.x - g.R, g.y - g.R, g.x + g.R, g.y + g.R, 1.9, 3.4);
  }

  // ------------------------------------------------------------------ Zeichnen
  function formPfad(art, x, y, r) {
    ctx.beginPath();
    if (art === "definition") { const a = r * 0.9; ctx.roundRect(x - a, y - a, 2 * a, 2 * a, r * 0.3); return; }
    if (art === "annahme") { ctx.roundRect(x - r * 1.3, y - r * 0.75, r * 2.6, r * 1.5, r * 0.75); return; }
    const p = formPunkte(art, r);
    if (!p) { ctx.arc(x, y, r, 0, 2 * Math.PI); return; }
    p.forEach(([px, py], i) => (i ? ctx.lineTo(x + px, y + py) : ctx.moveTo(x + px, y + py)));
    ctx.closePath();
  }
  const welt = () => ctx.setTransform(dpr * cam.k, 0, 0, dpr * cam.k, dpr * (W / 2 - cam.x * cam.k), dpr * (H / 2 - cam.y * cam.k));
  const schirm = () => ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const SERIF = '"KaTeX_Main", Cambria, Georgia, serif';
  const SANS = '"Segoe UI Variable Text", "Segoe UI", system-ui, sans-serif';

  function raster() {
    let s = 24 * cam.k;
    while (s < 22) s *= 2;
    while (s > 70) s /= 2;
    const ox = ((W / 2 - cam.x * cam.k) % s + s) % s, oy = ((H / 2 - cam.y * cam.k) % s + s) % s;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let x = ox; x < W; x += s) { ctx.moveTo(Math.round(x) + 0.5, 0); ctx.lineTo(Math.round(x) + 0.5, H); }
    for (let y = oy; y < H; y += s) { ctx.moveTo(0, Math.round(y) + 0.5); ctx.lineTo(W, Math.round(y) + 0.5); }
    ctx.strokeStyle = KF.raster;
    ctx.stroke();
  }

  function zeichne(now) {
    const Z = z(), k = cam.k, t = now / 1000;
    schirm();
    const gr = ctx.createRadialGradient(W / 2, H * 0.45, 0, W / 2, H * 0.45, Math.max(W, H) * 0.8);
    gr.addColorStop(0, KF.bg1); gr.addColorStop(1, KF.bg2);
    ctx.fillStyle = gr;
    ctx.fillRect(0, 0, W, H);
    raster();
    welt();

    const aGeb = 1 - glatt(1.45, 2.4, Z);
    const aHalo = 0.3 + 0.7 * (1 - glatt(1.3, 4.5, Z));
    const aStrang = glatt(1.3, 2.1, Z) * (1 - 0.6 * glatt(5.5, 10, Z));
    const aKante = glatt(1.35, 2.2, Z);
    const aTeilchen = glatt(1.35, 2.2, Z);

    // Gebiets-Höfe
    for (const g of L.gebiete.values()) {
      const h = ctx.createRadialGradient(g.x, g.y, g.R * 0.1, g.x, g.y, g.R * 1.25);
      h.addColorStop(0, rgba(g.farbe, 0.16 * aHalo * (fokus ? 0.5 : 1)));
      h.addColorStop(1, rgba(g.farbe, 0));
      ctx.fillStyle = h;
      ctx.beginPath(); ctx.arc(g.x, g.y, g.R * 1.25, 0, 2 * Math.PI); ctx.fill();
      if (aGeb > 0.01) {
        let w0 = -Math.PI / 2;
        ctx.lineWidth = 7 / k;
        for (const [st, anteil] of g.ring) {
          const w1 = w0 + 2 * Math.PI * anteil;
          ctx.beginPath(); ctx.arc(g.x, g.y, g.R, w0 + 0.01, w1 - 0.01);
          ctx.strokeStyle = rgba(KF.status[st], 0.85 * aGeb);
          ctx.stroke();
          w0 = w1;
        }
      }
    }
    // Bänder zwischen Gebieten (grobe Struktur)
    if (aGeb > 0.01) {
      ctx.globalCompositeOperation = "lighter";
      for (const bd of L.baender) {
        const { a, b } = bd;
        const cx = (a.x + b.x) * 0.32, cy = (a.y + b.y) * 0.32;
        const lg = ctx.createLinearGradient(a.x, a.y, b.x, b.y);
        lg.addColorStop(0, rgba(a.farbe, 0.32 * aGeb)); lg.addColorStop(1, rgba(b.farbe, 0.32 * aGeb));
        ctx.strokeStyle = lg;
        ctx.lineWidth = (1.5 + 2.4 * Math.sqrt(bd.w)) / k;
        ctx.lineCap = "round";
        ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.quadraticCurveTo(cx, cy, b.x, b.y); ctx.stroke();
        if (BEWEGUNG) {
          const n = Math.min(5, Math.ceil(bd.w / 3));
          for (let i = 0; i < n; i++) {
            const tt = ((t * 0.07 + bd.phase + i / n) % 1), u = 1 - tt;
            const px = u * u * a.x + 2 * u * tt * cx + tt * tt * b.x, py = u * u * a.y + 2 * u * tt * cy + tt * tt * b.y;
            const s = 9 / k;
            ctx.globalAlpha = 0.8 * aGeb;
            ctx.drawImage(glimmen(aufhellen(a.farbe, 0.4)), px - s, py - s, 2 * s, 2 * s);
          }
          ctx.globalAlpha = 1;
        }
      }
      ctx.globalCompositeOperation = "source-over";
    }
    // Strang-Flächen
    if (aStrang > 0.01) {
      for (const s of L.straenge.values()) {
        if (s.pseudo) continue;
        const g = L.gebiete.get(s.gebiet);
        const betroffen = !fokus || s.glieder.includes(fokus) || s.glieder.some((i) => NACHBARN.get(fokus).has(i));
        const a = aStrang * (betroffen ? 1 : 0.35);
        ctx.lineJoin = "round"; ctx.lineCap = "round";
        ctx.beginPath();
        s.huelle.forEach(([x, y], i) => (i ? ctx.lineTo(x, y) : ctx.moveTo(x, y)));
        if (s.huelle.length > 2) ctx.closePath();
        if (s.huelle.length === 1) ctx.lineTo(s.huelle[0][0] + 0.01, s.huelle[0][1]);
        ctx.lineWidth = 64;
        ctx.strokeStyle = rgba(g.farbe, 0.075 * a);
        ctx.stroke();
        ctx.fillStyle = rgba(g.farbe, 0.075 * a);
        if (s.huelle.length > 2) ctx.fill();
      }
    }
    // Kanten
    ctx.lineCap = "round";
    for (const kt of L.kanten) {
      const { e, kl } = kt;
      let a = aKante * Math.min(knotenAlpha(e.von), knotenAlpha(e.nach));
      const anFokus = fokus && (e.von === fokus || e.nach === fokus);
      if (fokus && !anFokus) a *= 0.55;
      if (kl === "leicht" || kl === "struktur") a *= 0.55;
      const anHover = hover && (e.von === hover || e.nach === hover);
      if (kt.quer && !anFokus && !anHover) a *= 0.12 + 0.45 * glatt(3.6, 6.5, Z);   // Querverbindungen zeigen die Bänder
      if (a < 0.01) continue;
      const farbe = KF.kante[kl];
      ctx.beginPath(); ctx.moveTo(kt.a.x, kt.a.y); ctx.quadraticCurveTo(kt.c[0], kt.c[1], kt.b.x, kt.b.y);
      ctx.setLineDash(kl === "rot" ? [6 / k, 4 / k] : kl === "leicht" ? [2 / k, 4 / k] : []);
      if (anFokus || hover && (e.von === hover || e.nach === hover)) {
        ctx.globalCompositeOperation = "lighter";
        ctx.lineWidth = 6 / k; ctx.strokeStyle = rgba(farbe, 0.16 * a); ctx.stroke();
        ctx.globalCompositeOperation = "source-over";
        ctx.lineWidth = 1.9 / k; ctx.strokeStyle = rgba(farbe, 0.95 * a);
      } else {
        ctx.lineWidth = 1.15 / k; ctx.strokeStyle = rgba(farbe, 0.6 * a);
      }
      ctx.stroke();
      if (Z > 2.6 && (kl === "beweis" || kl === "rot" || kl === "gruen")) {          // Pfeilspitze
        const [px, py] = qPunkt(kt, 1 - (kt.b.r + 3) / kt.len), [qx, qy] = qPunkt(kt, 1 - (kt.b.r + 12) / kt.len);
        const dx = px - qx, dy = py - qy, d = Math.hypot(dx, dy) + 1e-9, s = 5.5 / k;
        ctx.setLineDash([]);
        ctx.beginPath();
        ctx.moveTo(px, py);
        ctx.lineTo(px - (dx / d) * s * 1.6 + (dy / d) * s, py - (dy / d) * s * 1.6 - (dx / d) * s);
        ctx.lineTo(px - (dx / d) * s * 1.6 - (dy / d) * s, py - (dy / d) * s * 1.6 + (dx / d) * s);
        ctx.closePath();
        ctx.fillStyle = rgba(farbe, 0.85 * a);
        ctx.fill();
      }
    }
    ctx.setLineDash([]);
    // Lichtfluss entlang der Beweise
    if (BEWEGUNG && aTeilchen > 0.01) {
      ctx.globalCompositeOperation = "lighter";
      for (const kt of L.kanten) {
        if (kt.kl === "leicht" || kt.kl === "struktur") continue;
        const anFokus = fokus && (kt.e.von === fokus || kt.e.nach === fokus);
        if (kt.quer && !anFokus && Z < 5) continue;
        const a = aTeilchen * Math.min(knotenAlpha(kt.e.von), knotenAlpha(kt.e.nach)) * (fokus && !anFokus ? 0.35 : 1);
        if (a < 0.02) continue;
        const n = anFokus ? 2 : 1;
        for (let i = 0; i < n; i++) {
          const tt = ((t * 55) / kt.len + kt.phase + i / n) % 1;
          const [px, py] = qPunkt(kt, tt);
          const s = (anFokus ? 7 : 5) / k;
          ctx.globalAlpha = a * (anFokus ? 1 : 0.75) * Math.sin(Math.PI * tt);
          ctx.drawImage(glimmen(KF.kante[kt.kl]), px - s, py - s, 2 * s, 2 * s);
        }
      }
      ctx.globalAlpha = 1;
      ctx.globalCompositeOperation = "source-over";
    }
    // Knoten
    const liste = [...L.knoten.values()].sort((p, q) => p.r - q.r);
    const stern = Z < 1.55;
    for (const n of liste) {
      const a = knotenAlpha(n.id);
      const st = n.k._status_eff, c = KF.status[st];
      const r = stern ? Math.max(n.r * 0.55, 1.8 / k) : n.r;
      const istFokus = n.id === fokus, istHover = n.id === hover;
      // Leuchten
      ctx.globalCompositeOperation = "lighter";
      const puls = istFokus && BEWEGUNG ? 0.15 * Math.sin(t * 2.4) : 0;
      const gs = r * (istFokus ? 7 : istHover ? 6 : stern ? 5.5 : 4.6);
      ctx.globalAlpha = a * (stern ? 0.6 : istFokus ? 0.95 + puls : istHover ? 0.85 : 0.55);
      ctx.drawImage(glimmen(c), n.x - gs, n.y - gs, 2 * gs, 2 * gs);
      ctx.globalAlpha = 1;
      ctx.globalCompositeOperation = "source-over";
      // Kern
      formPfad(stern ? "satz" : n.k.art, n.x, n.y, r);
      if (HOHL.has(st)) { ctx.fillStyle = rgba(KF.bg1, a); ctx.fill(); ctx.lineWidth = 1.6 / k; ctx.strokeStyle = rgba(c, a); ctx.stroke(); }
      else if (st === "bedingt" || st === "skizze" || st === "annahme") {
        ctx.fillStyle = rgba(c, 0.35 * a); ctx.fill();
        ctx.lineWidth = 1.6 / k; ctx.setLineDash(st === "annahme" ? [] : [3.5 / k, 2.5 / k]); ctx.strokeStyle = rgba(c, a); ctx.stroke(); ctx.setLineDash([]);
      } else { ctx.fillStyle = rgba(c, a); ctx.fill(); }
      if (!stern && st === "widerlegt") {
        ctx.beginPath(); const q = r * 0.5;
        ctx.moveTo(n.x - q, n.y - q); ctx.lineTo(n.x + q, n.y + q); ctx.moveTo(n.x + q, n.y - q); ctx.lineTo(n.x - q, n.y + q);
        ctx.lineWidth = 2 / k; ctx.strokeStyle = rgba("#ffffff", a); ctx.stroke();
      }
      if (!stern && n.k._neu_pruefen?.length) {
        const rr = r + (5 + (BEWEGUNG ? 2 * Math.sin(t * 3) : 0)) / k;
        ctx.beginPath(); ctx.arc(n.x, n.y, rr, 0, 2 * Math.PI);
        ctx.lineWidth = 1.5 / k; ctx.setLineDash([3 / k, 2.5 / k]); ctx.strokeStyle = rgba(KF.status.widerlegt, a); ctx.stroke(); ctx.setLineDash([]);
      }
      if (istFokus) {
        ctx.beginPath(); ctx.arc(n.x, n.y, r + 9 / k, 0, 2 * Math.PI);
        ctx.lineWidth = 1.4 / k; ctx.strokeStyle = rgba("#ffffff", 0.7); ctx.stroke();
      }
      if (!stern && Z > 3.2 && r * k >= 9 && st !== "widerlegt") {
        ctx.font = `600 ${Math.min(12, r * k * 0.95) / k}px ${SANS}`;
        ctx.textAlign = "center"; ctx.textBaseline = "middle";
        ctx.fillStyle = HOHL.has(st) || st === "bedingt" || st === "skizze" || st === "annahme" ? rgba(c, a) : rgba("#0b1220", a);
        ctx.fillText(n.k.art === "frage" ? "?" : n.k.ebene, n.x, n.y + 0.5 / k);
      }
    }
    schirm();
    beschriften(Z);
    karten(Z);
    stufeAnzeigen(Z);
  }

  // Bedienelemente über der Leinwand: dort keine Beschriftung und keine Karten
  function sperrflaechen() {
    const r0 = leinwand.getBoundingClientRect();
    return ["#pfad", "#filter", "#stufen", ".zoom"].map((q) => $(q)).filter((e) => e && e.offsetParent && e.style.visibility !== "hidden")
      .map((e) => { const r = e.getBoundingClientRect(); return [r.left - r0.left - 6, r.top - r0.top - 6, r.right - r0.left + 6, r.bottom - r0.top + 6]; });
  }
  let sperre = [], sperreZeit = 0;

  // ------------------------------------------------------------------ Beschriftung (Schirmraum, mit Kollisionsprüfung)
  function beschriften(Z) {
    if (performance.now() - sperreZeit > 500) { sperre = sperrflaechen(); sperreZeit = performance.now(); }
    const belegt = sperre.slice();
    if (Z >= 1.75) {                                       // Knoten selbst freihalten
      for (const n of L.knoten.values()) {
        if (knotenAlpha(n.id) < 0.3) continue;
        const [sx, sy] = zuSchirm(n.x, n.y), r = n.r * cam.k + 3;
        if (sx > -30 && sx < W + 30 && sy > -30 && sy < H + 30) belegt.push([sx - r, sy - r, sx + r, sy + r]);
      }
    }
    const frei = (x0, y0, x1, y1) => {
      for (const b of belegt) if (x0 < b[2] && x1 > b[0] && y0 < b[3] && y1 > b[1]) return false;
      belegt.push([x0, y0, x1, y1]);
      return true;
    };
    ctx.textBaseline = "alphabetic";
    // Gebiete
    const aGeb = 1 - glatt(1.45, 2.4, Z);
    if (aGeb > 0.01) {
      for (const g of L.gebiete.values()) {
        const [sx, sy] = zuSchirm(g.x, g.y);
        const fs = Math.max(15, Math.min(30, g.R * cam.k * 0.2));
        ctx.font = `600 ${fs}px ${SERIF}`;
        ctx.textAlign = "center";
        const zeilen = umbruch(g.name, 16, 2);
        const w = Math.max(...zeilen.map((zl) => ctx.measureText(zl).width));
        frei(sx - w / 2, sy - fs * 1.2, sx + w / 2, sy + fs * zeilen.length);
        ctx.shadowColor = "rgba(0,0,0,.8)"; ctx.shadowBlur = 8;
        zeilen.forEach((zl, i) => { ctx.fillStyle = rgba(KF.text, aGeb); ctx.fillText(zl, sx, sy + (i - (zeilen.length - 1) / 2) * fs * 1.1); });
        ctx.font = `13px ${SANS}`;
        ctx.fillStyle = rgba(KF.text2, aGeb);
        ctx.fillText(`${g.n} Aussagen · ${g.offen} offen${g.neu ? ` · ${g.neu} neu prüfen` : ""}`, sx, sy + (zeilen.length / 2) * fs * 1.1 + 12);
        ctx.shadowBlur = 0;
      }
    }
    // Aussagen-Kandidaten: Fokus, Hover und Nachbarn haben Vorrang vor allen anderen Texten
    const kandidaten = [];
    for (const n of L.knoten.values()) {
      if (sichtbarKarte.has(n.id)) continue;               // hat eine Inhaltskarte
      if (fokus && n.id !== fokus && !NACHBARN.get(fokus).has(n.id) && !zweiter.has(n.id)) continue;
      const wichtig = n.id === fokus || n.id === hover || (fokus && NACHBARN.get(fokus).has(n.id));
      const a = Math.max(glatt(3.1, 4.1, Z), (n.k._tragweite || 0) >= 4 ? glatt(2.2, 3.0, Z) : 0, wichtig ? glatt(1.6, 2.2, Z) : 0)
        * Math.max(knotenAlpha(n.id), wichtig ? 1 : 0);
      if (a < 0.03) continue;
      const prio = (n.id === fokus ? 1e6 : 0) + (n.id === hover ? 5e5 : 0) + (wichtig ? 1e5 : 0) + (n.k._tragweite || 0) * 10 + n.r;
      kandidaten.push([prio, n, a, wichtig]);
    }
    kandidaten.sort((p, q) => q[0] - p[0]);
    let anzahl = 0;
    const etikett = (n, a) => {
      const [sx, sy] = zuSchirm(n.x, n.y);
      if (sx < -100 || sx > W + 100 || sy < -40 || sy > H + 40) return;
      const fs = n.id === fokus ? 14 : 12.5;
      ctx.font = `${n.id === fokus ? "600 " : ""}${fs}px ${SERIF}`;
      ctx.textAlign = "center";
      const text = kurz(reinText(n.k.titel), n.id === fokus ? 60 : 36);
      const w = ctx.measureText(text).width;
      const y = sy + n.r * cam.k + 18;
      const rect = [sx - w / 2 - 3, y - fs, sx + w / 2 + 3, y + 4];
      if (!frei(...rect)) {
        if (n.id !== fokus) return;
        belegt.push(rect);                                 // Fokus wird immer beschriftet
      }
      ctx.shadowColor = "rgba(0,0,0,.9)"; ctx.shadowBlur = 5;
      ctx.fillStyle = rgba(KF.text, a);
      ctx.fillText(text, sx, y);
      ctx.shadowBlur = 0;
      anzahl++;
    };
    for (const [, n, a, wichtig] of kandidaten) if (wichtig) etikett(n, a);
    // Stränge: Name und eine Zeile Verlauf
    const aName = glatt(1.4, 2.0, Z) * (1 - glatt(5.6, 6.8, Z));
    const aVerlauf = glatt(1.6, 2.3, Z) * (1 - glatt(4.3, 5.6, Z));
    if (aName > 0.01) {
      const st = [...L.straenge.values()].filter((s) => !s.pseudo).sort((p, q) => q.glieder.length - p.glieder.length);
      for (const s of st) {
        const betroffen = !fokus || s.glieder.includes(fokus) || s.glieder.some((i) => NACHBARN.get(fokus).has(i));
        const a = betroffen ? 1 : 0.35;
        const g = L.gebiete.get(s.gebiet);
        const [sx, sy] = zuSchirm(s.x, s.oben);
        ctx.font = `600 15px ${SERIF}`;
        ctx.textAlign = "center";
        const w = ctx.measureText(s.name).width;
        const zeilen = aVerlauf > 0.02 ? umbruch(reinText(s.verlauf), 46, 3) : [];
        const hoehe = 18 + zeilen.length * 15;
        const breite = Math.max(w, zeilen.length ? 290 : 0);
        if (!frei(sx - breite / 2, sy - hoehe, sx + breite / 2, sy + 4)) continue;
        ctx.shadowColor = "rgba(0,0,0,.85)"; ctx.shadowBlur = 6;
        ctx.fillStyle = rgba(aufhellen(g.farbe, 0.35), aName * a);
        ctx.fillText(s.name, sx, sy - zeilen.length * 15);
        ctx.font = `12.5px ${SANS}`;
        zeilen.forEach((zl, i) => { ctx.fillStyle = rgba(KF.text2, aVerlauf * a); ctx.fillText(zl, sx, sy - (zeilen.length - 1 - i) * 15 + 2); });
        ctx.shadowBlur = 0;
      }
    }
    // übrige Aussagen, nur wo Platz ist
    for (const [, n, a, wichtig] of kandidaten) {
      if (anzahl > 90) break;
      if (!wichtig) etikett(n, a);
    }
  }

  // ------------------------------------------------------------------ Inhaltskarten (DOM, nur ganz nah)
  const kartenCache = new Map();
  const sichtbarKarte = new Set();
  function karteHtml(k) {
    return `<div class="kk-kopf"><span class="punkt" style="background:${KF.status[k._status_eff]}"></span>${esc(S.arten[k.art])} · ${esc(S.status[k._status_eff]?.label)} · Ebene ${k.ebene}</div>
      <div class="kk-titel">${mathe(k.titel)}</div><div class="kk-text">${mathe(k.aussage)}</div>`;
  }
  function karten(Z) {
    const a = glatt(6.4, 7.6, Z);
    sichtbarKarte.clear();
    if (a < 0.02) { kartenEbene.innerHTML = ""; kartenCache.clear(); return; }
    const kand = [...L.knoten.values()].filter((n) => knotenAlpha(n.id) > 0.3).map((n) => {
      const [sx, sy] = zuSchirm(n.x, n.y);
      return { n, sx, sy, d: Math.hypot(sx - W / 2, sy - H / 2) - (n.id === fokus ? 1e5 : 0) - (fokus && NACHBARN.get(fokus).has(n.id) ? 5e4 : 0) };
    }).filter((c) => c.sx > -50 && c.sx < W + 50 && c.sy > -50 && c.sy < H + 50).sort((p, q) => p.d - q.d).slice(0, 9);
    const belegt = sperre.slice();
    const behalten = new Set();
    for (const { n, sx, sy } of kand) {
      let el2 = kartenCache.get(n.id);
      if (!el2) {
        el2 = document.createElement("button");
        el2.className = "inhaltskarte";
        el2.dataset.id = n.id;
        el2.innerHTML = karteHtml(n.k);
        kartenEbene.appendChild(el2);
        kartenCache.set(n.id, el2);
      }
      const w = el2.offsetWidth || 270, h = el2.offsetHeight || 120, r = n.r * cam.k;
      const versuche = [[sx + r + 14, sy - h / 2], [sx - r - 14 - w, sy - h / 2], [sx - w / 2, sy + r + 12], [sx - w / 2, sy - r - 12 - h]];
      let ort = null;
      for (const [x, y] of versuche) {
        if (x < 4 || y < 40 || x + w > W - 4 || y + h > H - 4) continue;
        if (belegt.some((b) => x < b[2] && x + w > b[0] && y < b[3] && y + h > b[1])) continue;
        ort = [x, y];
        break;
      }
      if (!ort) { el2.style.display = "none"; continue; }
      belegt.push([ort[0], ort[1], ort[0] + w, ort[1] + h]);
      el2.style.display = "";
      el2.style.transform = `translate(${ort[0].toFixed(1)}px, ${ort[1].toFixed(1)}px)`;
      el2.style.opacity = a;
      el2.classList.toggle("fokus", n.id === fokus);
      behalten.add(n.id);
      sichtbarKarte.add(n.id);
    }
    for (const [id, el2] of kartenCache) if (!behalten.has(id)) { el2.remove(); kartenCache.delete(id); }
  }

  // ------------------------------------------------------------------ Stufenleiste
  const leiste = $("#stufen");
  leiste.innerHTML = STUFEN.map((s, i) => `<button data-stufe="${i}" aria-pressed="false"><span class="stufe-punkt"></span><span class="stufe-name">${s.name}</span></button>`).join("");
  let letzteStufe = -1;
  function stufeAnzeigen(Z) {
    const s = stufeVon(Z);
    if (s === letzteStufe) return;
    letzteStufe = s;
    leiste.querySelectorAll("button").forEach((b, i) => b.setAttribute("aria-pressed", String(i === s)));
    $("#filter").style.visibility = s === 0 ? "hidden" : "visible";
  }
  leiste.addEventListener("click", (e) => {
    const b = e.target.closest("[data-stufe]");
    if (!b) return;
    const i = +b.dataset.stufe;
    if (i === 0) { setzeFokus(null); ganz(); return; }
    const n = fokus && L.knoten.get(fokus);
    fliege(n ? n.x : cam.x, n ? n.y : cam.y, k0 * STUFEN[i].z);
  });

  // ------------------------------------------------------------------ Trefferprüfung und Eingabe
  function trefferKnoten(sx, sy) {
    const [wx, wy] = zuWelt(sx, sy);
    let best = null, bd = Infinity;
    for (const n of L.knoten.values()) {
      if (knotenAlpha(n.id) < 0.08) continue;
      const d = Math.hypot(n.x - wx, n.y - wy) * cam.k;
      const grenze = Math.max(9, n.r * cam.k + 5);
      if (d < grenze && d < bd) { bd = d; best = n.id; }
    }
    return best;
  }
  function trefferStrang(sx, sy) {
    const [wx, wy] = zuWelt(sx, sy);
    let best = null, bd = Infinity;
    for (const s of L.straenge.values()) {
      if (s.pseudo) continue;
      const d = Math.hypot(s.x - wx, s.y - wy);
      const dl = Math.hypot(s.x - wx, s.oben - wy);
      const m = Math.min(d / s.rad, dl * cam.k / 60);
      if (m < 1 && m < bd) { bd = m; best = s.id; }
    }
    return best;
  }
  function trefferGebiet(sx, sy) {
    const [wx, wy] = zuWelt(sx, sy);
    for (const g of L.gebiete.values()) if (Math.hypot(g.x - wx, g.y - wy) < g.R) return g.id;
    return null;
  }
  function klick(sx, sy) {
    const id = trefferKnoten(sx, sy);
    if (id && z() >= 1.55) { waehle(id); return; }
    const Z = z();
    if (Z < 1.75) {
      const g = trefferGebiet(sx, sy);
      if (g) { zeigeGebiet(g); zuGebiet(g); return; }
    } else if (Z < 6.6) {
      const s = trefferStrang(sx, sy);
      if (s) { zeigeStrang(s); zuStrang(s); return; }
    }
    if (fokus) waehle(null);
  }
  function zeigeTipp(sx, sy, html) {
    tipp.innerHTML = html;
    tipp.hidden = false;
    let x = sx + 14, y = sy + 12;
    if (x + 300 > W) x = sx - 310;
    tipp.style.left = x + "px"; tipp.style.top = y + "px";
  }
  const zeiger = new Map();
  let zug = null, pinch = null;
  leinwand.addEventListener("pointerdown", (e) => {
    leinwand.setPointerCapture(e.pointerId);
    zeiger.set(e.pointerId, [e.offsetX, e.offsetY]);
    if (zeiger.size === 2) {
      const [p, q] = [...zeiger.values()];
      pinch = { d: Math.hypot(p[0] - q[0], p[1] - q[1]), k: cam.k };
      zug = null;
    } else zug = { x: e.offsetX, y: e.offsetY, cx: cam.x, cy: cam.y, bewegt: false };
    flug = null;
  });
  leinwand.addEventListener("pointermove", (e) => {
    if (zeiger.has(e.pointerId)) zeiger.set(e.pointerId, [e.offsetX, e.offsetY]);
    if (pinch && zeiger.size === 2) {
      const [p, q] = [...zeiger.values()];
      cam.k = Math.max(k0 * 0.6, Math.min(k0 * 16, pinch.k * Math.hypot(p[0] - q[0], p[1] - q[1]) / pinch.d));
      anstossen();
      return;
    }
    if (zug) {
      const dx = e.offsetX - zug.x, dy = e.offsetY - zug.y;
      if (Math.abs(dx) + Math.abs(dy) > 3) zug.bewegt = true;
      if (zug.bewegt) { cam.x = zug.cx - dx / cam.k; cam.y = zug.cy - dy / cam.k; leinwand.classList.add("zieht"); tipp.hidden = true; anstossen(); }
      return;
    }
    const id = trefferKnoten(e.offsetX, e.offsetY);
    if (id !== hover) { hover = id; anstossen(); }
    if (id && z() >= 1.55) {
      const k = K.get(id);
      leinwand.style.cursor = "pointer";
      zeigeTipp(e.offsetX, e.offsetY, `<div>${esc(reinText(k.titel))}</div><div class="tipp-status">${esc(S.arten[k.art])} · Ebene ${k.ebene} · ${esc(S.status[k._status_eff]?.label)}${k._neu_pruefen?.length ? " · neu prüfen" : ""}</div>`);
    } else if (z() < 1.75 && trefferGebiet(e.offsetX, e.offsetY)) {
      const g = GEB.get(trefferGebiet(e.offsetX, e.offsetY));
      leinwand.style.cursor = "pointer";
      zeigeTipp(e.offsetX, e.offsetY, `<div><b>${esc(g.name)}</b></div><div class="tipp-status">${esc(reinText(g.kurz))}</div><div class="tipp-status">Klicken: hineinzoomen</div>`);
    } else if (z() >= 1.75 && z() < 6.6 && trefferStrang(e.offsetX, e.offsetY)) {
      const s = L.straenge.get(trefferStrang(e.offsetX, e.offsetY));
      leinwand.style.cursor = "pointer";
      zeigeTipp(e.offsetX, e.offsetY, `<div><b>${esc(s.name)}</b></div><div class="tipp-status">${esc(reinText(s.verlauf))}</div>`);
    } else { leinwand.style.cursor = ""; tipp.hidden = true; }
  });
  const loslassen = (e) => {
    zeiger.delete(e.pointerId);
    if (zeiger.size < 2) pinch = null;
    if (zug && !zug.bewegt) klick(zug.x, zug.y);
    zug = null;
    leinwand.classList.remove("zieht");
  };
  leinwand.addEventListener("pointerup", loslassen);
  leinwand.addEventListener("pointercancel", loslassen);
  leinwand.addEventListener("pointerleave", () => { if (hover) { hover = null; anstossen(); } tipp.hidden = true; });
  leinwand.addEventListener("wheel", (e) => {
    e.preventDefault();
    flug = null;
    const f = Math.exp(-e.deltaY * 0.0016);
    const [wx, wy] = zuWelt(e.offsetX, e.offsetY);
    cam.k = Math.max(k0 * 0.6, Math.min(k0 * 16, cam.k * f));
    cam.x = wx - (e.offsetX - W / 2) / cam.k; cam.y = wy - (e.offsetY - H / 2) / cam.k;
    anstossen();
  }, { passive: false });
  leinwand.addEventListener("dblclick", (e) => {
    const [wx, wy] = zuWelt(e.offsetX, e.offsetY);
    fliege(wx, wy, cam.k * 2, 500);
  });
  kartenEbene.addEventListener("click", (e) => { const b = e.target.closest(".inhaltskarte"); if (b) waehle(b.dataset.id); });
  new ResizeObserver(() => { if (aktiv) groesse(); }).observe(leinwand.parentElement);

  // ------------------------------------------------------------------ Schleife
  function anstossen() { if (!laeuft && aktiv) { laeuft = true; requestAnimationFrame(bild); } }
  function bild(now) {
    laeuft = false;
    if (!aktiv) return;
    flugSchritt(now);
    zeichne(now);
    if (flug || BEWEGUNG) { laeuft = true; requestAnimationFrame(bild); }
  }

  layout();
  return {
    zeige() {
      aktiv = true;
      leinwand.hidden = false; kartenEbene.hidden = false; leiste.hidden = false;
      groesse();
      if (!erstesMal) {
        erstesMal = true;
        const [x0, y0, x1, y1] = L.box;
        Object.assign(cam, { x: (x0 + x1) / 2, y: (y0 + y1) / 2 + 20 / k0, k: k0 });
      }
      document.fonts?.load(`600 20px "KaTeX_Main"`).then(anstossen);
      anstossen();
    },
    verberge() { aktiv = false; leinwand.hidden = true; kartenEbene.hidden = true; leiste.hidden = true; tipp.hidden = true; },
    fokus(id, fliegen = true) { setzeFokus(id); if (id && fliegen) zuKnoten(id); },
    zuGebiet, zuStrang, ganz, anstossen,
    zoom(f) { fliege(cam.x, cam.y, cam.k * f, 380); },
    ort() {                                   // wo die Kamera gerade hinschaut, für die Pfadleiste
      const Z = z();
      let g = null, s = null, bd = Infinity;
      if (Z >= 1.75) for (const q of L.gebiete.values()) { const d = Math.hypot(q.x - cam.x, q.y - cam.y) / q.R; if (d < bd) { bd = d; g = q.id; } }
      bd = Infinity;
      if (Z >= 3) for (const q of L.straenge.values()) { if (q.pseudo || q.gebiet !== g) continue; const d = Math.hypot(q.x - cam.x, q.y - cam.y); if (d < bd) { bd = d; s = q.id; } }
      return { gebiet: g, strang: s, stufe: stufeVon(Z) };
    },
    strangVon: (id) => L.strangVon.get(id),
    get aktiv() { return aktiv; },
  };
})();
