// Curved-surface demo: a small civic hall whose architecture is driven by
// analytic curved panels rather than by adding curved objects to the office.
//
// CSG members provide the plinth, glass walls, mullions and ring beams. The
// barrel vault, rotunda dome and conical roof lantern are registered through
// the wasm curved-panel API and tessellated by the kernel adapter.

const P = (x, y, z) => ({ x, y, z });
const rect = (halfW, halfH) => ({ Rect: { half_w: halfW, half_h: halfH } });
const circle = (r) => ({ Circle: { radius: r } });
const hSection = (halfW, halfH, web, flange) => ({
  HSection: { half_w: halfW, half_h: halfH, web, flange },
});

const extrude = (profile, origin, axis, length) => ({
  Extrude: { profile, origin, axis, length },
});

const boxZ = (cx, cy, z0, sizeX, sizeY, h) =>
  extrude(rect(sizeY / 2, sizeX / 2), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);
const cylZ = (cx, cy, z0, r, h) => extrude(circle(r), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);
const prismXY = (cx, cy, z0, angle, length, thickness, height) => {
  const dx = Math.cos(angle);
  const dy = Math.sin(angle);
  return extrude(
    rect(thickness / 2, height / 2),
    P(cx - (dx * length) / 2, cy - (dy * length) / 2, z0 + height / 2),
    { x: dx, y: dy, z: 0 },
    length,
  );
};
function ringSegments(kind, label, cx, cy, z0, radius, segments, thickness, height, lengthScale = 0.94) {
  for (let i = 0; i < segments; i++) {
    const a = (Math.PI * 2 * (i + 0.5)) / segments;
    const length = 2 * radius * Math.sin(Math.PI / segments) * lengthScale;
    add(kind, `${label} ${i}`, prismXY(
      cx + Math.cos(a) * radius,
      cy + Math.sin(a) * radius,
      z0,
      a + Math.PI / 2,
      length,
      thickness,
      height,
    ));
  }
}

let nextId = 1;
const members = [];
function add(kind, label, node) {
  const id = nextId++;
  members.push({ id, kind, label, node });
  return id;
}

let nextCurvedId = 10001;
const curvedPanels = [];
function addCurved(kind, label, node) {
  const id = nextCurvedId++;
  curvedPanels.push({ id, kind, label, node });
  return id;
}

// ── Plinth and main hall ─────────────────────────────────────────────────────
add('plinth', 'stone plinth', boxZ(0.8, 0.0, 0.0, 28.4, 10.2, 0.35));
add('plinth', 'entry forecourt', boxZ(-13.4, 0.0, 0.02, 5.0, 6.6, 0.12));
add('glass', 'south glass wall', boxZ(0.0, -4.12, 0.35, 22.2, 0.12, 3.95));
add('glass', 'north glass wall', boxZ(0.0, 4.12, 0.35, 22.2, 0.12, 3.95));
add('glass', 'west glass end', boxZ(-10.9, 0.0, 0.35, 0.12, 8.0, 3.95));
add('roof-frame', 'south spring beam', boxZ(0.0, -4.16, 4.22, 22.6, 0.22, 0.28));
add('roof-frame', 'north spring beam', boxZ(0.0, 4.16, 4.22, 22.6, 0.22, 0.28));
add('roof-frame', 'west portal beam', boxZ(-10.9, 0.0, 4.2, 0.26, 8.35, 0.32));

for (let i = 0; i < 8; i++) {
  const x = -10.5 + i * 3.0;
  add('mullion', `south vault mullion ${i}`, boxZ(x, -4.24, 0.35, 0.16, 0.16, 4.05));
  add('mullion', `north vault mullion ${i}`, boxZ(x, 4.24, 0.35, 0.16, 0.16, 4.05));
}
for (const z of [1.7, 2.95]) {
  add('mullion', `south horizontal mullion z${z}`, boxZ(0.0, -4.26, z, 22.1, 0.08, 0.08));
  add('mullion', `north horizontal mullion z${z}`, boxZ(0.0, 4.26, z, 22.1, 0.08, 0.08));
}

// ── Rotunda and connection bay ───────────────────────────────────────────────
const ROT = { x: 12.9, y: 0.0, r: 4.3 };
add('glass', 'rotunda link glass bay', boxZ(10.7, 0.0, 0.35, 2.8, 7.3, 3.7));
add('plinth', 'rotunda circular plinth', cylZ(ROT.x, ROT.y, 0.0, 4.95, 0.42));
ringSegments('glass', 'rotunda glass panel', ROT.x, ROT.y, 0.38, 4.28, 20, 0.1, 3.28);
ringSegments('roof-frame', 'rotunda dome ring beam', ROT.x, ROT.y, 3.58, 4.25, 20, 0.32, 0.36, 1.0);
ringSegments('roof-frame', 'lantern curb segment', ROT.x, ROT.y, 6.78, 1.13, 12, 0.2, 0.24, 1.0);

for (let i = 0; i < 12; i++) {
  const a = (Math.PI * 2 * i) / 12;
  const x = ROT.x + Math.cos(a) * 4.55;
  const y = ROT.y + Math.sin(a) * 4.55;
  add('mullion', `rotunda perimeter column ${i}`, cylZ(x, y, 0.35, 0.075, 3.62));
}

for (let i = 0; i < 5; i++) {
  const y = -2.4 + i * 1.2;
  add(
    'mullion',
    `entry canopy fin ${i}`,
    extrude(hSection(0.055, 0.07, 0.006, 0.008), P(-14.6, y, 3.15), { x: 1, y: 0, z: 0 }, 5.5),
  );
}

// ── Analytic curved panels ──────────────────────────────────────────────────
addCurved('curved-roof', 'thick barrel vault with rectangular skylight trims', {
  cylinder: {
    axis_origin: P(-11.1, 0.0, 1.22),
    axis_dir: P(1.0, 0.0, 0.0),
    radius: 5.12,
    theta_min: -0.9,
    theta_max: 0.9,
    z_min: 0.0,
    z_max: 22.0,
    thickness: 0.18,
    holes: [
      { rectangle: { u_min: -0.16, u_max: 0.16, v_min: 2.1, v_max: 4.3, reverse: true } },
      { rectangle: { u_min: -0.16, u_max: 0.16, v_min: 6.1, v_max: 8.3, reverse: true } },
      { rectangle: { u_min: -0.16, u_max: 0.16, v_min: 10.1, v_max: 12.3, reverse: true } },
      { rectangle: { u_min: -0.16, u_max: 0.16, v_min: 14.1, v_max: 16.3, reverse: true } },
    ],
  },
});

addCurved('curved-dome', 'rotunda glass dome with UV oculi', {
  sphere: {
    center: P(ROT.x, ROT.y, 3.98),
    radius: 4.15,
    pole: P(0.0, 0.0, 1.0),
    theta_min: 0.08,
    theta_max: 6.2,
    phi_min: 0.12,
    phi_max: 1.48,
    holes: [
      { circle: { center: [0.9, 0.92], radius: 0.11, reverse: true } },
      { circle: { center: [2.1, 0.92], radius: 0.11, reverse: true } },
      { circle: { center: [3.3, 0.92], radius: 0.11, reverse: true } },
      { circle: { center: [4.5, 0.92], radius: 0.11, reverse: true } },
      { circle: { center: [5.7, 0.92], radius: 0.11, reverse: true } },
    ],
  },
});

addCurved('curved-cone', 'copper conical roof lantern with slot trims', {
  cone: {
    apex: P(ROT.x, ROT.y, 8.74),
    axis: P(0.0, 0.0, -1.0),
    half_angle: 0.48,
    theta_min: 0.08,
    theta_max: 6.2,
    height_min: 0.72,
    height_max: 1.95,
    holes: [
      { rectangle: { u_min: 0.7, u_max: 1.05, v_min: 1.05, v_max: 1.55, reverse: true } },
      { rectangle: { u_min: 2.05, u_max: 2.4, v_min: 1.05, v_max: 1.55, reverse: true } },
      { rectangle: { u_min: 3.4, u_max: 3.75, v_min: 1.05, v_max: 1.55, reverse: true } },
      { rectangle: { u_min: 4.75, u_max: 5.1, v_min: 1.05, v_max: 1.55, reverse: true } },
    ],
  },
});

export const KIND_OF = new Map(members.map((m) => [m.id, m.kind]));
export const CURVED_KIND_OF = new Map(curvedPanels.map((m) => [m.id, m.kind]));
export const LABEL_OF = new Map(members.map((m) => [m.id, m.label]));
export const MEMBERS = members;
export const CURVED_PANELS = curvedPanels;
export const DEMO = {
  id: 'curved',
  label: 'Curved hall',
  description: 'Curved-panel civic hall',
};
export const BOUNDS = {
  x: [-15.9, 18.0],
  y: [-5.25, 5.25],
  zMax: 8.95,
};
export const VIEW = {
  target: [2.5, 0.0, 3.8],
  heroCamera: [33.0, -25.0, 14.8],
  heroTarget: [1.5, 0.0, 3.9],
  sectionCamera: [30.0, -23.0, 13.6],
  sectionTarget: [1.6, 0.0, 3.4],
  sectionZ: 3.15,
};
