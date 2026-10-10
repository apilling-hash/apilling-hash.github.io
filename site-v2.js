// Motion + 3D for the portfolio (v2). Loaded from the page logic via dynamic import.
const THREE_URL = 'https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js';

const clamp = (v, a = 0, b = 1) => Math.min(b, Math.max(a, v));
const smooth = (a, b, x) => { const t = clamp((x - a) / (b - a)); return t * t * (3 - 2 * t); };
const easeInOut = t => t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
const easeOut = t => 1 - Math.pow(1 - t, 3);
const damp = (a, b, k, dt) => a + (b - a) * (1 - Math.exp(-k * dt));
const EASE = 'cubic-bezier(.2,.7,.2,1)';
const inView = (el, m = 100) => { const r = el.getBoundingClientRect(); return r.bottom > -m && r.top < innerHeight + m; };

const ATOMS = {
  Dusk: { Si: 0xd8d4ca, O: 0x8fb3ff, C: 0x5a6272, H: 0xeceae4, bond: 0x6b7383, body: '#0f1219' },
  Night: { Si: 0xcfc6b4, O: 0xe3a65b, C: 0x5d6470, H: 0xe9e6df, bond: 0x6e7580, body: '#10141b' },
  Paper: { Si: 0xc9ccd3, O: 0x2d55c8, C: 0x2a2d33, H: 0xfafaf7, bond: 0xa9adb5, body: '#f3f3ef' },
  Warm:  { Si: 0xd9ccb2, O: 0xb0503a, C: 0x35322e, H: 0xf3efe7, bond: 0xbab2a4, body: '#f5f2ec' }
};
const DEEP = {
  Dusk: { Si: 0xd8d4ca, O: 0x8fb3ff, C: 0x5a6272, H: 0xeceae4, bond: 0x5d6575 },
  Night: { Si: 0xcfc6b4, O: 0xe3a65b, C: 0x5d6470, H: 0xe9e6df, bond: 0x5f6672 },
  Paper: { Si: 0xc9ccd3, O: 0x93aaf2, C: 0x6c727c, H: 0xf3f3ef, bond: 0x565b64 },
  Warm:  { Si: 0xd9ccb2, O: 0xe08a72, C: 0x6f6a63, H: 0xf3efe7, bond: 0x5f5a53 }
};
const RAD = { Si: 0.52, O: 0.4, C: 0.37, H: 0.2 };
const BOND_R = 0.105;
const EL = ['Si', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'];
const INTRA = [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [2, 6], [3, 7], [3, 8], [3, 9]];
const NAMES = {
  3: 'Hexamethylcyclotrisiloxane', 4: 'Octamethylcyclotetrasiloxane', 5: 'Decamethylcyclopentasiloxane',
  6: 'Dodecamethylcyclohexasiloxane', 7: 'Tetradecamethylcycloheptasiloxane', 8: 'Hexadecamethylcyclooctasiloxane'
};
const CIDS = { 3: 10914, 4: 11169, 5: 10913, 6: 10911, 7: 7874, 8: 11170 };
const ringR = n => 1.64 / (2 * Math.sin(Math.PI / (2 * n)));

let T;
const opts = { palette: 'Dusk', showH: true };

const geomCache = {};
function siloxane(n) {
  if (geomCache[n]) return geomCache[n];
  const V = T.Vector3, R = ringR(n);
  const A = [0, 0.3, 0.55, 0.8, 1.0, 1.2][n - 3];
  const alt = n > 3 ? 0.16 * Math.min(1, (n - 3) / 2) : 0;
  const ring = (j, r) => {
    const phi = j * Math.PI / n;
    return new V(r * Math.cos(phi), A * 0.6 * Math.cos(2 * phi) + (j % 2 ? -alt : alt), r * Math.sin(phi));
  };
  const up = new V(0, 1, 0), tet = 54.75 * Math.PI / 180, hA = 70.5 * Math.PI / 180;
  const units = [];
  for (let k = 0; k < n; k++) {
    const si = ring(2 * k, R * 0.95), o = ring(2 * k + 1, R * 1.05);
    const phi = 2 * k * Math.PI / n, rad = new V(Math.cos(phi), 0, Math.sin(phi));
    const u = [si, o];
    const cs = [1, -1].map(sg => si.clone().add(rad.clone().multiplyScalar(Math.cos(tet) * 1.87)).add(up.clone().multiplyScalar(sg * Math.sin(tet) * 1.87)));
    u.push(cs[0], cs[1]);
    cs.forEach((c, ci) => {
      const d = c.clone().sub(si).normalize();
      const e1 = new V().crossVectors(d, Math.abs(d.y) < 0.9 ? up : new V(1, 0, 0)).normalize(), e2 = new V().crossVectors(d, e1).normalize();
      for (let m = 0; m < 3; m++) {
        const al = m * 2 * Math.PI / 3 + ci * 0.5 + k * 0.3;
        const dir = d.clone().multiplyScalar(Math.cos(hA)).add(e1.clone().multiplyScalar(Math.sin(hA) * Math.cos(al))).add(e2.clone().multiplyScalar(Math.sin(hA) * Math.sin(al)));
        u.push(c.clone().add(dir.multiplyScalar(1.09)));
      }
    });
    units.push(u);
  }
  return (geomCache[n] = units);
}

const mols = [];
class Mol {
  constructor(capUnits, lite, deep) {
    this.deep = deep;
    const seg = lite ? 18 : 32;
    const sph = new T.SphereGeometry(1, seg, Math.round(seg * 0.7));
    const cyl = new T.CylinderGeometry(1, 1, 1, lite ? 10 : 18, 1, true);
    const mat = () => new T.MeshPhysicalMaterial({ roughness: 0.46, metalness: 0, clearcoat: 0.35, clearcoatRoughness: 0.45 });
    this.group = new T.Group();
    this.m = {};
    const caps = { Si: capUnits, O: capUnits, C: 2 * capUnits, H: 6 * capUnits };
    for (const el in caps) {
      const mesh = new T.InstancedMesh(sph, mat(), caps[el]);
      mesh.instanceMatrix.setUsage(T.DynamicDrawUsage); mesh.frustumCulled = false;
      this.m[el] = mesh; this.group.add(mesh);
    }
    this.bond = new T.InstancedMesh(cyl, new T.MeshStandardMaterial({ roughness: 0.6 }), 11 * capUnits);
    this.bond.instanceMatrix.setUsage(T.DynamicDrawUsage); this.bond.frustumCulled = false;
    this.group.add(this.bond);
    this._m = new T.Matrix4(); this._q = new T.Quaternion(); this._s = new T.Vector3();
    this._v = new T.Vector3(); this._mid = new T.Vector3(); this._up = new T.Vector3(0, 1, 0);
    this.colors(); mols.push(this);
  }
  colors() {
    const c = (this.deep ? DEEP : ATOMS)[opts.palette] || ATOMS.Dusk;
    for (const el in this.m) this.m[el].material.color.setHex(c[el]);
    this.bond.material.color.setHex(c.bond);
  }
  // atoms: [{el,p,s}], bonds: [[p1,p2,s,isH]]
  draw(atoms, bonds) {
    const { _m, _q, _s, _v, _mid, _up } = this, showH = opts.showH;
    const c = { Si: 0, O: 0, C: 0, H: 0 };
    _q.identity();
    for (const a of atoms) {
      if (a.s <= 0.001 || (!showH && a.el === 'H')) continue;
      const r = RAD[a.el] * a.s; _s.set(r, r, r);
      _m.compose(a.p, _q, _s); this.m[a.el].setMatrixAt(c[a.el]++, _m);
    }
    for (const el in c) { this.m[el].count = c[el]; this.m[el].instanceMatrix.needsUpdate = true; }
    let b = 0;
    for (const [p1, p2, s, isH] of bonds) {
      if (s <= 0.001 || (isH && !showH)) continue;
      _v.subVectors(p2, p1); const len = _v.length(); if (len < 1e-4) continue;
      _q.setFromUnitVectors(_up, _v.divideScalar(len));
      _s.set(BOND_R * s, len, BOND_R * s);
      _mid.addVectors(p1, p2).multiplyScalar(0.5);
      _m.compose(_mid, _q, _s); this.bond.setMatrixAt(b++, _m);
    }
    this.bond.count = b; this.bond.instanceMatrix.needsUpdate = true;
  }
}

function makeStage(canvas, fov, lite) {
  const renderer = new T.WebGLRenderer({ canvas, antialias: !lite || devicePixelRatio < 2, alpha: true, powerPreference: 'high-performance' });
  renderer.setPixelRatio(Math.min(devicePixelRatio, lite ? 1.5 : 2));
  renderer.outputColorSpace = T.SRGBColorSpace;
  renderer.toneMapping = T.NeutralToneMapping;
  const scene = new T.Scene();
  const camera = new T.PerspectiveCamera(fov, 1, 0.1, 200);
  scene.add(new T.HemisphereLight(0xfff7ec, 0x8e95a0, 1.5));
  const key = new T.DirectionalLight(0xffffff, 2.1); key.position.set(5, 7, 9); scene.add(key);
  const fill = new T.DirectionalLight(0xdfe6ff, 0.7); fill.position.set(-7, -3, 4); scene.add(fill);
  const st = { renderer, scene, camera, w: 1, h: 1, aspect: 1 };
  st.resize = () => {
    const r = canvas.parentElement.getBoundingClientRect();
    st.w = Math.max(1, r.width); st.h = Math.max(1, r.height); st.aspect = st.w / st.h;
    renderer.setSize(st.w, st.h, false); camera.aspect = st.aspect; camera.updateProjectionMatrix();
  };
  st.resize();
  return st;
}

const flatten = units => { const a = []; units.forEach(u => u.forEach((p, i) => a.push({ el: EL[i], p: p.clone(), s: 1 }))); return a; };
function bondIdx(n) {
  const out = [];
  for (let k = 0; k < n; k++) {
    INTRA.forEach(([a, b]) => out.push([k * 10 + a, k * 10 + b, b >= 4]));
    out.push([k * 10 + 1, ((k + 1) % n) * 10, false]);
  }
  return out;
}

// Exploding / reassembling molecule (hero + contact)
function burstScene(ctx, cfg) {
  const canvas = document.getElementById(cfg.canvas), wrap = document.getElementById(cfg.wrap);
  if (!canvas || !wrap) return null;
  const st = makeStage(canvas, 30, ctx.lite);
  const n = cfg.n;
  const mol = new Mol(n, ctx.lite, cfg.deep); st.scene.add(mol.group);
  const base = flatten(siloxane(n)), bi = bondIdx(n);
  let seed = cfg.seed || 7; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const dirs = base.map(a => {
    const d = a.p.clone().normalize().multiplyScalar(0.8).add(new T.Vector3(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5));
    return d.normalize().multiplyScalar((2.4 + rnd() * 3.6 + (a.el === 'H' ? 1.2 : 0)) * (cfg.spread || 1));
  });
  const spin = base.map(() => (rnd() - 0.5) * 1.6);
  const atoms = base.map(a => ({ el: a.el, p: a.p.clone(), s: 1 }));
  const fitR = ringR(n) + 2.3;
  let rx = 0.5, ry = 0.6, e = cfg.startE || 0, time = 0;

  function layout() {
    const dist = 16, halfH = Math.tan(st.camera.fov * Math.PI / 360) * dist, halfW = halfH * st.aspect;
    st.camera.position.set(0, 0, dist);
    const l = cfg.layout(st.aspect, halfW, halfH);
    mol.group.scale.setScalar(l.size / fitR);
    mol.group.position.set(l.x, l.y, 0);
  }
  layout();

  function frame(dt, s) {
    time += dt;
    const target = ctx.still ? 0 : cfg.explode(wrap.getBoundingClientRect());
    e = ctx.still ? 0 : damp(e, target, 5, dt);
    const ee = easeInOut(clamp(e));
    for (let i = 0; i < atoms.length; i++) {
      const d = dirs[i], a = spin[i] * ee, ca = Math.cos(a), sa = Math.sin(a);
      atoms[i].p.set(base[i].p.x + (d.x * ca - d.z * sa) * ee, base[i].p.y + d.y * ee, base[i].p.z + (d.x * sa + d.z * ca) * ee);
      atoms[i].s = 1 - 0.12 * ee;
    }
    const bs = 1 - smooth(0.04, 0.32, ee);
    mol.draw(atoms, bi.map(([a, b, h]) => [atoms[a].p, atoms[b].p, bs, h]));
    const tx = 0.5 + s.my * 0.28, ty = (ctx.still ? 0.6 : time * 0.16) + s.mx * 0.55 + ee * 0.8;
    rx = damp(rx, tx, 3, dt); ry = ctx.still ? ty : damp(ry, ty, 3, dt);
    mol.group.rotation.set(rx, ry, 0.12);
    st.renderer.render(st.scene, st.camera);
  }
  return {
    frame(dt, s) { if (inView(wrap)) frame(dt, s); },
    resize() { st.resize(); layout(); if (ctx.still) frame(0, { mx: 0, my: 0 }); },
    renderOnce() { frame(0, { mx: 0, my: 0 }); },
    dispose() { st.renderer.dispose(); }
  };
}

// ---------- Research ----------
function setStep(i, animate) {
  const n = i + 3;
  document.querySelectorAll('[data-step]').forEach(b => {
    const on = +b.dataset.step === i;
    b.style.color = on ? 'var(--ink)' : 'var(--muted)';
    b.setAttribute('aria-current', on ? 'step' : 'false');
    const num = b.querySelector('[data-num]'); if (num) num.style.color = on ? 'var(--accent)' : 'var(--muted)';
  });
  document.querySelectorAll('[data-desc]').forEach(d => {
    const on = +d.dataset.desc === i;
    d.style.gridTemplateRows = on ? '1fr' : '0fr';
    d.style.opacity = on ? '1' : '0';
  });
  const tag = document.getElementById('rs-tag');
  if (tag) {
    tag.textContent = `D${n} · ${NAMES[n]}`;
    tag.href = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + CIDS[n];
    if (animate) tag.animate([{ opacity: 0, transform: 'translate3d(0,6px,0)' }, { opacity: 1, transform: 'none' }], { duration: 500, easing: EASE });
  }
}

function buildRing(nf, pool) {
  const n0 = Math.min(8, Math.floor(nf + 1e-6)), t = n0 >= 8 ? 0 : nf - n0;
  const A = siloxane(n0), B = t > 0 ? siloxane(n0 + 1) : null;
  const atoms = [], bonds = [], idx = [];
  let pi = 0; const P = () => pool[pi] || (pool[pi] = new T.Vector3()), nextP = () => { const p = P(); pi++; return p; };
  const tg = easeOut(clamp(t * 1.15));
  for (let k = 0; k < n0; k++) {
    idx[k] = [];
    for (let i = 0; i < 10; i++) {
      const p = nextP(); p.copy(A[k][i]); if (B) p.lerp(B[k][i], t);
      idx[k][i] = atoms.length; atoms.push({ el: EL[i], p, s: 1 });
    }
  }
  if (B) {
    const M = A[n0 - 1][1].clone().add(A[0][0]).multiplyScalar(0.5);
    idx[n0] = [];
    for (let i = 0; i < 10; i++) {
      const p = nextP(); p.copy(M).lerp(B[n0][i], tg);
      idx[n0][i] = atoms.length; atoms.push({ el: EL[i], p, s: tg });
    }
  }
  const units = B ? n0 + 1 : n0;
  for (let k = 0; k < units; k++) {
    const s = k === n0 ? tg : 1;
    INTRA.forEach(([a, b]) => bonds.push([atoms[idx[k][a]].p, atoms[idx[k][b]].p, s, b >= 4]));
  }
  for (let k = 0; k < n0 - 1; k++) bonds.push([atoms[idx[k][1]].p, atoms[idx[k + 1][0]].p, 1, false]);
  const last = atoms[idx[n0 - 1][1]].p, first = atoms[idx[0][0]].p;
  if (B) {
    bonds.push([last, first, 1 - smooth(0, 0.5, t), false]);
    bonds.push([last, atoms[idx[n0][0]].p, tg, false]);
    bonds.push([atoms[idx[n0][1]].p, first, tg, false]);
  } else bonds.push([last, first, 1, false]);
  return { atoms, bonds, R: B ? ringR(n0) + (ringR(n0 + 1) - ringR(n0)) * t : ringR(n0) };
}

function researchScene(ctx) {
  const canvas = document.getElementById('rs-canvas'), wrap = document.getElementById('rs-wrap');
  if (!canvas || !wrap) return null;
  const st = makeStage(canvas, 30, ctx.lite);
  const mol = new Mol(9, ctx.lite); st.scene.add(mol.group);
  const pool = [];
  let nf = 3, shown = -1, time = 0, yaw = 0, pitch = 0, vy = 0, vp = 0, dragging = false, lx = 0, ly = 0;
  setStep(0, false); shown = 0;

  const down = e => { dragging = true; lx = e.clientX; ly = e.clientY; vy = vp = 0; canvas.setPointerCapture?.(e.pointerId); canvas.style.cursor = 'grabbing'; };
  const move = e => {
    if (!dragging) return;
    const dx = e.clientX - lx, dy = e.clientY - ly; lx = e.clientX; ly = e.clientY;
    yaw += dx * 0.008; pitch = clamp(pitch + dy * 0.006, -0.9, 0.9); vy = dx * 0.48; vp = dy * 0.36;
  };
  const up = () => { dragging = false; canvas.style.cursor = 'grab'; };
  canvas.addEventListener('pointerdown', down); canvas.addEventListener('pointermove', move);
  canvas.addEventListener('pointerup', up); canvas.addEventListener('pointercancel', up);

  const steps = [...document.querySelectorAll('[data-step]')];
  const onStep = ev => {
    const i = +ev.currentTarget.dataset.step;
    const top = wrap.getBoundingClientRect().top + scrollY, span = wrap.offsetHeight - innerHeight;
    scrollTo({ top: top + span * (i + 0.2) / 6, behavior: 'smooth' });
  };
  steps.forEach(b => b.addEventListener('click', onStep));

  function frame(dt) {
    time += dt;
    const r = wrap.getBoundingClientRect(), span = Math.max(1, r.height - innerHeight);
    const p = clamp(-r.top / span) * 6, i = Math.min(5, Math.floor(p)), f = p - i;
    nf = damp(nf, 3 + i + (i < 5 ? easeInOut(smooth(0.6, 0.95, f)) : 0), 7, dt);
    const step = Math.round(nf) - 3;
    if (step !== shown) { shown = step; setStep(step, true); }
    const { atoms, bonds, R } = buildRing(nf, pool);
    mol.draw(atoms, bonds);
    st.camera.position.set(0, 0, 15 / Math.min(1, st.aspect * 1.05));
    const halfH = Math.tan(st.camera.fov * Math.PI / 360) * 15;
    mol.group.scale.setScalar(halfH * (0.62 + 0.035 * (nf - 3)) / (R + 2.3));
    if (!dragging) { yaw += vy * dt; pitch = clamp(pitch + vp * dt, -0.9, 0.9); vy = damp(vy, 0, 3, dt); vp = damp(vp, 0, 3, dt); pitch = damp(pitch, 0, 0.6, dt); }
    mol.group.rotation.set(0.62 + pitch, time * 0.12 + yaw, 0);
    st.renderer.render(st.scene, st.camera);
  }
  return {
    frame(dt) { if (inView(wrap)) frame(dt); },
    resize() { st.resize(); },
    dispose() {
      steps.forEach(b => b.removeEventListener('click', onStep));
      canvas.removeEventListener('pointerdown', down); canvas.removeEventListener('pointermove', move);
      canvas.removeEventListener('pointerup', up); canvas.removeEventListener('pointercancel', up);
      st.renderer.dispose();
    }
  };
}

function stillRing() {
  const img = document.getElementById('rs-still'); if (!img) return;
  const c = document.createElement('canvas'), holder = document.createElement('div');
  holder.style.cssText = 'position:fixed;left:-9999px;top:0;width:640px;height:640px';
  holder.appendChild(c); document.body.appendChild(holder);
  const st = makeStage(c, 30, false);
  const mol = new Mol(8, false); st.scene.add(mol.group);
  const { atoms, bonds, R } = buildRing(6, []);
  mol.draw(atoms, bonds);
  st.camera.position.set(0, 0, 15);
  mol.group.scale.setScalar(Math.tan(15 * Math.PI / 180) * 15 * 0.8 / (R + 2.3));
  mol.group.rotation.set(0.7, 0.4, 0);
  st.renderer.render(st.scene, st.camera);
  img.src = c.toDataURL('image/png');
  st.renderer.dispose(); holder.remove();
  mols.splice(mols.indexOf(mol), 1);
}

// ---------- Simple Industry: tilt-in, then a look at the module list ----------
function shotScene() {
  const wrap = document.getElementById('si-wrap'), fr = document.getElementById('si-frame'), sh = document.getElementById('si-shadow');
  const zoom = document.getElementById('si-zoom'), note = document.getElementById('si-note');
  if (!wrap || !fr) return null;
  let ps = 0, last = -1;
  function apply(p) {
    const a = easeInOut(clamp(p / 0.42));
    const z = easeInOut(smooth(0.5, 0.66, p)) * (1 - easeInOut(smooth(0.84, 0.97, p)));
    fr.style.transform = `perspective(1800px) translate3d(0, ${(1 - a) * 7}vh, 0) rotateX(${(1 - a) * 26}deg) rotateZ(${(1 - a) * -5}deg) scale(${0.5 + 0.5 * a})`;
    if (sh) {
      sh.style.transform = `perspective(1800px) translate3d(0, ${(1 - a) * 14 + 3}vh, 0) rotateX(${(1 - a) * 26}deg) rotateZ(${(1 - a) * -5}deg) scale(${(0.5 + 0.5 * a) * 0.94})`;
      sh.style.opacity = String(0.6 - 0.25 * a);
    }
    if (zoom) zoom.style.transform = `scale(${1 + 1.35 * z})`;
    if (note) { note.style.opacity = String(smooth(0.25, 0.9, z)); note.style.transform = `translate3d(${(1 - z) * 48}px,0,0)`; }
  }
  return {
    frame(dt) {
      if (!inView(wrap, 200)) return;
      const r = wrap.getBoundingClientRect(), span = Math.max(1, r.height - innerHeight);
      const p = clamp((-r.top + innerHeight * 0.25) / (span + innerHeight * 0.25));
      ps = damp(ps, p, 9, dt);
      if (Math.abs(ps - last) > 1e-4) { last = ps; apply(ps); }
    }
  };
}

// ---------- Reveals ----------
function reveals() {
  const els = [...document.querySelectorAll('[data-reveal]')];
  const lineSets = [...document.querySelectorAll('[data-lines]')];
  const rules = [...document.querySelectorAll('[data-rule]')];
  const all = [...els, ...lineSets, ...rules];
  function show(el) {
    if (el.dataset.shown) return; el.dataset.shown = '1';
    const d = +(el.dataset.delay || 0);
    if (el.hasAttribute('data-lines')) {
      [...el.querySelectorAll('[data-line] > span')].forEach((s, i) => {
        s.style.transition = `transform 1.15s ${EASE} ${d + i * 90}ms`; s.style.transform = 'none';
      });
    } else if (el.hasAttribute('data-rule')) {
      el.style.transition = `transform 1.4s ${EASE} ${d}ms`; el.style.transform = 'scaleX(1)';
    } else {
      el.style.transition = `opacity .9s ${EASE} ${d}ms, transform 1.1s ${EASE} ${d}ms`;
      el.style.opacity = '1'; el.style.transform = 'none';
    }
  }
  els.forEach(el => { el.style.opacity = '0'; el.style.transform = 'translate3d(0,22px,0)'; });
  rules.forEach(el => { el.style.transformOrigin = '0 50%'; el.style.transform = 'scaleX(0)'; });
  lineSets.forEach(el => el.querySelectorAll('[data-line] > span').forEach(s => { s.style.display = 'inline-block'; s.style.transform = 'translate3d(0,108%,0)'; }));
  const check = () => all.forEach(el => { if (!el.dataset.shown && el.getBoundingClientRect().top < innerHeight * 0.92) show(el); });
  requestAnimationFrame(() => requestAnimationFrame(check));
  addEventListener('scroll', check, { passive: true });
  return () => removeEventListener('scroll', check);
}

export async function init(o = {}) {
  Object.assign(opts, o);
  document.body.style.background = (ATOMS[opts.palette] || ATOMS.Dusk).body;
  const still = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const lite = innerWidth < 760 || matchMedia('(pointer: coarse)').matches;
  const ctx = { still, lite };
  const cleanups = [];
  if (!still) cleanups.push(reveals());
  const api = { dispose: () => cleanups.forEach(f => f && f()), setPalette() {}, setShowH() {} };

  try { T = await import(THREE_URL); } catch (e) { console.warn('three.js failed to load', e); return api; }

  const scenes = [];
  const add = f => { try { const s = f(); s && scenes.push(s); return s; } catch (e) { console.warn(e); } };
  const hero = add(() => burstScene(ctx, {
    canvas: 'hero-canvas', wrap: 'hero-wrap', n: 4,
    layout: (asp, hw, hh) => asp > 1.15 ? { size: Math.min(hh * 0.74, hw * 0.36), x: hw * 0.46, y: -hh * 0.02 } : { size: Math.min(hw * 0.72, hh * 0.32), x: 0, y: hh * 0.04 },
    explode: r => Math.sin(Math.PI * smooth(0.02, 0.98, clamp(-r.top / Math.max(1, r.height - innerHeight))))
  }));
  const contact = add(() => burstScene(ctx, {
    canvas: 'ct-canvas', wrap: 'contact', n: 6, deep: true, seed: 31, spread: 1.4, startE: 1,
    layout: (asp, hw, hh) => ({ size: Math.min(hw * 0.8, hh * 0.8), x: 0, y: 0 }),
    explode: r => 1 - smooth(0.05, 0.75, (innerHeight - r.top) / innerHeight)
  }));
  if (still) {
    hero && hero.renderOnce(); contact && contact.renderOnce();
    try { stillRing(); } catch (e) { console.warn(e); }
  } else {
    add(() => researchScene(ctx));
    add(() => shotScene());
  }

  const s = { mx: 0, my: 0, tx: 0, ty: 0 };
  const onMouse = e => { s.tx = (e.clientX / innerWidth) * 2 - 1; s.ty = (e.clientY / innerHeight) * 2 - 1; };
  addEventListener('pointermove', onMouse, { passive: true });
  let raf = 0, prev = performance.now();
  const loop = now => {
    const dt = Math.min(0.05, (now - prev) / 1000); prev = now;
    s.mx = damp(s.mx, s.tx, 4, dt); s.my = damp(s.my, s.ty, 4, dt);
    scenes.forEach(sc => sc.frame && sc.frame(dt, s));
    raf = requestAnimationFrame(loop);
  };
  if (!still) raf = requestAnimationFrame(loop);
  let rt = 0;
  const onResize = () => { clearTimeout(rt); rt = setTimeout(() => scenes.forEach(sc => sc.resize && sc.resize()), 120); };
  addEventListener('resize', onResize);
  const redrawStill = () => { if (still) { hero && hero.renderOnce(); contact && contact.renderOnce(); try { stillRing(); } catch (e) {} } };

  api.setPalette = name => {
    opts.palette = name; document.body.style.background = (ATOMS[name] || ATOMS.Dusk).body;
    mols.forEach(m => m.colors()); redrawStill();
  };
  api.setShowH = v => { opts.showH = v; redrawStill(); };
  api.dispose = () => {
    cancelAnimationFrame(raf);
    removeEventListener('pointermove', onMouse); removeEventListener('resize', onResize);
    scenes.forEach(sc => sc.dispose && sc.dispose());
    cleanups.forEach(f => f && f());
    mols.length = 0;
  };
  return api;
}
