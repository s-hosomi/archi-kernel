// Curved-surface demo inspired by Hagia Sophia's massing: a central dome,
// east/west half domes, heavy buttresses and four minarets. This is not a
// historical reconstruction; it is a compact architectural demo that uses the
// current analytic curved-panel API without relying on fragile angled detail.

const P = (x, y, z) => ({ x, y, z });
const rect = (halfW, halfH) => ({ Rect: { half_w: halfW, half_h: halfH } });
const circle = (r) => ({ Circle: { radius: r } });

const extrude = (profile, origin, axis, length) => ({
  Extrude: { profile, origin, axis, length },
});
const openings = (base, list) =>
  list.length
    ? { OpeningSubtraction: { base, openings: list.map((shape, i) => [i + 1, { shape }]) } }
    : base;

// Vertical prism helpers (axis +Z: half_w spans Y, half_h spans X).
const boxZ = (cx, cy, z0, sizeX, sizeY, h) =>
  extrude(rect(sizeY / 2, sizeX / 2), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);
const cylZ = (cx, cy, z0, r, h) => extrude(circle(r), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);

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

const TAU = Math.PI * 2;
const FULL0 = 0.0;
const FULL1 = TAU;

const uvRectangle = (uMin, uMax, vMin, vMax) => ({
  rectangle: { u_min: uMin, u_max: uMax, v_min: vMin, v_max: vMax, reverse: true },
});
const uvPolygon = (points) => ({ polygon: { points, reverse: true } });
const pointedArchTrim = (u, halfWidth, baseV, shoulderV, apexV) =>
  uvPolygon([
    [u - halfWidth, baseV],
    [u + halfWidth, baseV],
    [u + halfWidth, shoulderV],
    [u, apexV],
    [u - halfWidth, shoulderV],
  ]);
const diamondTrim = (u, v, halfWidth, halfHeight) =>
  uvPolygon([
    [u, v - halfHeight],
    [u + halfWidth, v],
    [u, v + halfHeight],
    [u - halfWidth, v],
  ]);

const domeClerestoryTrims = Array.from({ length: 24 }, (_, i) =>
  pointedArchTrim((i + 0.5) * (TAU / 24), 0.066, 1.42, 1.25, 1.12),
);
const domeUpperJewelTrims = Array.from({ length: 12 }, (_, i) =>
  diamondTrim((i + 0.5) * (TAU / 12), 0.76, 0.055, 0.085),
);
const domeWindowTrims = [...domeClerestoryTrims, ...domeUpperJewelTrims];
const vaultSkylightTrims = [3.0, 7.0, 11.0, 15.0, 19.0, 23.0].map((z) =>
  uvRectangle(-0.18, 0.18, z - 0.72, z + 0.72),
);
const halfDomeTrims = (thetaCenter) =>
  [-0.72, -0.36, 0.0, 0.36, 0.72].map((offset) =>
    pointedArchTrim(thetaCenter + offset, 0.09, 1.06, 0.9, 0.74),
  );
const minaretCapTrims = Array.from({ length: 4 }, (_, i) => {
  const u = (i + 0.5) * (TAU / 4);
  return uvRectangle(u - 0.055, u + 0.055, 0.72, 1.28);
});

function wallWithWindowsY(y, label) {
  const body = boxZ(0.0, y, 0.45, 30.0, 0.65, 5.45);
  const holes = [-11.5, -7.7, -3.9, 3.9, 7.7, 11.5].map((x) =>
    boxZ(x, y, 2.0, 1.15, 0.95, 2.05),
  );
  add('stone', label, openings(body, holes));
  for (const x of [-11.5, -7.7, -3.9, 3.9, 7.7, 11.5]) {
    add('window', `${label} window ${x}`, boxZ(x, y + Math.sign(y) * 0.03, 2.05, 0.96, 0.05, 1.88));
  }
}

function endWallWithWindowsX(x, label) {
  const body = boxZ(x, 0.0, 0.45, 0.75, 17.4, 5.45);
  const holes = [-5.4, -2.0, 2.0, 5.4].map((y) =>
    boxZ(x, y, 2.0, 0.95, 1.12, 2.05),
  );
  add('stone', label, openings(body, holes));
  for (const y of [-5.4, -2.0, 2.0, 5.4]) {
    add('window', `${label} window ${y}`, boxZ(x + Math.sign(x) * 0.03, y, 2.05, 0.05, 0.94, 1.88));
  }
}

function minaret(x, y, label) {
  add('minaret', `${label} shaft`, cylZ(x, y, 0.4, 0.46, 14.2));
  add('minaret', `${label} lower balcony`, cylZ(x, y, 9.25, 0.72, 0.22));
  add('minaret', `${label} upper balcony`, cylZ(x, y, 12.0, 0.62, 0.2));
  addCurved('curved-cone', `${label} conical cap`, {
    cone: {
      apex: P(x, y, 16.7),
      axis: P(0.0, 0.0, -1.0),
      half_angle: 0.24,
      theta_min: FULL0,
      theta_max: FULL1,
      height_min: 0.28,
      height_max: 2.18,
      holes: minaretCapTrims,
    },
  });
}

// Base platform and enclosed nave.
add('plinth', 'stepped stone plinth lower', boxZ(0.0, 0.0, 0.0, 34.6, 20.2, 0.25));
add('plinth', 'stepped stone plinth upper', boxZ(0.0, 0.0, 0.25, 32.6, 18.2, 0.24));
add('stone', 'main prayer hall mass', boxZ(0.0, 0.0, 0.49, 26.0, 13.8, 4.8));
wallWithWindowsY(-8.2, 'south arcade wall');
wallWithWindowsY(8.2, 'north arcade wall');
endWallWithWindowsX(-15.2, 'west entry wall');
endWallWithWindowsX(15.2, 'east apse wall');

// Central square, drum and buttresses.
add('stone', 'central square roof base', boxZ(0.0, 0.0, 5.25, 14.4, 14.4, 1.25));
add('drum', 'central circular dome drum', cylZ(0.0, 0.0, 6.35, 7.1, 0.95));
for (const [x, y] of [
  [-6.65, -6.65],
  [6.65, -6.65],
  [-6.65, 6.65],
  [6.65, 6.65],
]) {
  add('buttress', `main pier ${x},${y}`, boxZ(x, y, 0.49, 2.45, 2.45, 7.15));
}
for (const x of [-11.7, 0.0, 11.7]) {
  add('buttress', `south exterior buttress ${x}`, boxZ(x, -9.15, 0.49, 1.35, 1.6, 6.1));
  add('buttress', `north exterior buttress ${x}`, boxZ(x, 9.15, 0.49, 1.35, 1.6, 6.1));
}
for (const y of [-5.7, 5.7]) {
  add('buttress', `west exterior buttress ${y}`, boxZ(-16.0, y, 0.49, 1.6, 1.35, 6.1));
  add('buttress', `east exterior buttress ${y}`, boxZ(16.0, y, 0.49, 1.6, 1.35, 6.1));
}

// Lower side aisle barrel vaults. They are opaque, simple and kept inside the
// masonry envelope so no curved surface protrudes beyond its support.
addCurved('curved-stone-vault', 'south side aisle barrel vault', {
  cylinder: {
    axis_origin: P(-13.0, -5.55, 3.4),
    axis_dir: P(1.0, 0.0, 0.0),
    radius: 1.75,
    theta_min: -0.78,
    theta_max: 0.78,
    z_min: 0.0,
    z_max: 26.0,
    thickness: 0.12,
    holes: vaultSkylightTrims,
  },
});
addCurved('curved-stone-vault', 'north side aisle barrel vault', {
  cylinder: {
    axis_origin: P(-13.0, 5.55, 3.4),
    axis_dir: P(1.0, 0.0, 0.0),
    radius: 1.75,
    theta_min: -0.78,
    theta_max: 0.78,
    z_min: 0.0,
    z_max: 26.0,
    thickness: 0.12,
    holes: vaultSkylightTrims,
  },
});

// Main dome and east/west half domes.
addCurved('curved-stone-dome', 'central dome with UV clerestory trims', {
  sphere: {
    center: P(0.0, 0.0, 7.15),
    radius: 6.8,
    pole: P(0.0, 0.0, 1.0),
    theta_min: FULL0,
    theta_max: FULL1,
    phi_min: 0.02,
    phi_max: 1.52,
    holes: domeWindowTrims,
    thickness: 0.18,
  },
});
addCurved('curved-stone-dome', 'east half dome', {
  sphere: {
    center: P(6.2, 0.0, 4.55),
    radius: 5.35,
    pole: P(1.0, 0.0, 0.0),
    theta_min: -1.26,
    theta_max: 1.26,
    phi_min: 0.02,
    phi_max: 1.2,
    holes: halfDomeTrims(0.0),
    thickness: 0.14,
  },
});
addCurved('curved-stone-dome', 'west half dome', {
  sphere: {
    center: P(-6.2, 0.0, 4.55),
    radius: 5.35,
    pole: P(-1.0, 0.0, 0.0),
    theta_min: Math.PI - 1.26,
    theta_max: Math.PI + 1.26,
    phi_min: 0.02,
    phi_max: 1.2,
    holes: halfDomeTrims(Math.PI),
    thickness: 0.14,
  },
});

// Minarets, kept outside the plinth corners.
minaret(-17.1, -9.8, 'south-west minaret');
minaret(17.1, -9.8, 'south-east minaret');
minaret(-17.1, 9.8, 'north-west minaret');
minaret(17.1, 9.8, 'north-east minaret');

export const KIND_OF = new Map(members.map((m) => [m.id, m.kind]));
export const CURVED_KIND_OF = new Map(curvedPanels.map((m) => [m.id, m.kind]));
export const LABEL_OF = new Map(members.map((m) => [m.id, m.label]));
export const MEMBERS = members;
export const CURVED_PANELS = curvedPanels;
export const DEMO = {
  id: 'curved',
  label: 'Ayasofya study',
  description: 'Hagia Sophia inspired curved-panel massing study',
};
export const BOUNDS = {
  x: [-18.2, 18.2],
  y: [-10.7, 10.7],
  zMax: 16.9,
};
export const VIEW = {
  target: [0.0, 0.0, 6.2],
  heroCamera: [46.0, -39.0, 23.0],
  heroTarget: [0.0, 0.0, 7.2],
  sectionCamera: [38.0, -32.0, 20.0],
  sectionTarget: [0.0, 0.0, 5.8],
  sectionZ: 6.25,
  shots: {
    front: { camera: [0.0, -52.0, 16.0], target: [0.0, 0.0, 6.8] },
    back: { camera: [0.0, 52.0, 16.0], target: [0.0, 0.0, 6.8] },
    east: { camera: [54.0, 0.0, 16.0], target: [0.0, 0.0, 6.8] },
    west: { camera: [-54.0, 0.0, 16.0], target: [0.0, 0.0, 6.8] },
    ne: { camera: [40.0, 38.0, 21.0], target: [0.0, 0.0, 7.0] },
    nw: { camera: [-40.0, 38.0, 21.0], target: [0.0, 0.0, 7.0] },
    se: { camera: [40.0, -38.0, 21.0], target: [0.0, 0.0, 7.0] },
    sw: { camera: [-40.0, -38.0, 21.0], target: [0.0, 0.0, 7.0] },
  },
};
