// Demo building: a three-storey RC office block, modelled entirely as
// archi-kernel CSG and evaluated by the kernel at load time.
//
// Architectural programme
//   - 3 × 2 structural grid (6.4 m × 5.6 m bays), storey heights 4.2 / 3.6 / 3.6
//   - top storey set back one bay → roof terrace with parapet and a steel pergola
//   - round-column colonnade along the front elevation (Y = 0)
//   - core (stair + elevator) with cast walls and door openings, voids stacked
//     through every slab
//   - glazed front: each storey's facade wall carries a window grid
//     (batched OpeningSubtraction), ground floor gets the entrance + canopy
//   - full deduction chain for quantity take-off: columns ▷ girders ▷ beams ▷ slabs,
//     walls deducted by their columns (公共建築数量積算基準の優先控除)
//
// Everything below is plain JSON in archi-kernel's serde representation.
// Units: metres, Z up.
//
// Profile-frame conventions (kernel's deterministic plane_basis):
//   axis +Z → (u, v) = (+Y, −X): half_w spans Y, half_h spans X
//   axis +X → (u, v) = (+Z, −Y): half_w spans Z, half_h spans Y
//   axis +Y →                    half_w spans Z, half_h spans X (verified)

const P = (x, y, z) => ({ x, y, z });
const rect = (halfW, halfH) => ({ Rect: { half_w: halfW, half_h: halfH } });
const circle = (r) => ({ Circle: { radius: r } });
const hSection = (halfW, halfH, web, flange) => ({
  HSection: { half_w: halfW, half_h: halfH, web, flange },
});

const extrude = (profile, origin, axis, length) => ({
  Extrude: { profile, origin, axis, length },
});
const clip = (base, clippers) =>
  clippers.length ? { Clip: { base, clippers, rule: 'Priority' } } : base;
const openings = (base, list) =>
  list.length
    ? { OpeningSubtraction: { base, openings: list.map((shape, i) => [i + 1, { shape }]) } }
    : base;

// Vertical prism helpers (axis +Z: half_w spans Y, half_h spans X).
const boxZ = (cx, cy, z0, sizeX, sizeY, h) =>
  extrude(rect(sizeY / 2, sizeX / 2), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);
const cylZ = (cx, cy, z0, r, h) => extrude(circle(r), P(cx, cy, z0), { x: 0, y: 0, z: 1 }, h);

// ── Grid & levels ────────────────────────────────────────────────────────────
const GX = [0, 6.4, 12.8, 19.2]; // column lines, X (3 bays)
const GY = [0, 5.6, 11.2]; // column lines, Y (2 bays)
const LEVELS_Z = [4.2, 7.8, 11.4]; // structural floor tops, storeys 1..3
const STOREY0 = [0, 4.2, 7.8]; // storey bottoms
const SETBACK_X = 12.8; // top storey exists only for x ≤ 12.8
const COL = 0.65;
const ROUND_R = 0.325;
const SLAB_T = 0.22;
const G = { w: 0.4, d: [0.75, 0.7, 0.7] }; // girders per storey
const B = { w: 0.3, d: 0.55 }; // secondary beams
const MARGIN = 0.35; // slab cantilever beyond the grid

// Core (stair + elevator) in the middle bay, rear side.
const CORE = { x0: 7.4, x1: 11.8, y0: 8.0, y1: 11.0 };
const STAIR = { x0: 7.7, x1: 10.2, y0: 8.6, y1: 10.4 }; // stair well void
const ELEV = { x0: 10.5, x1: 11.5, y0: 8.6, y1: 10.4 }; // elevator shaft void

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

// Storey exists at (x, s)? Top storey (s=2) is set back beyond SETBACK_X.
const inStorey = (s, x) => s < 2 || x <= SETBACK_X + 1e-9;

// ── Columns ──────────────────────────────────────────────────────────────────
// Front row (y = 0) is a round-column colonnade; the rest are square RC.
const colId = new Map(); // `${s}:${x}:${y}` -> id
for (let s = 0; s < 3; s++) {
  const z0 = STOREY0[s];
  const h = LEVELS_Z[s] - z0;
  for (const x of GX) {
    if (!inStorey(s, x)) continue;
    for (const y of GY) {
      const round = y === 0;
      const node = round ? cylZ(x, y, z0, ROUND_R, h) : boxZ(x, y, z0, COL, COL, h);
      const id = add(
        round ? 'column-round' : 'column',
        `${round ? 'round col' : 'col'} (${x},${y}) S${s + 1}`,
        node,
      );
      colId.set(`${s}:${x}:${y}`, id);
    }
  }
}
const colsOf = (s) =>
  [...colId.entries()].filter(([k]) => k.startsWith(`${s}:`)).map(([, id]) => id);

// ── Girders (top flush with the floor above, clipped by their columns) ──────
const girderIds = [[], [], []]; // per storey
for (let s = 0; s < 3; s++) {
  const zTop = LEVELS_Z[s];
  const d = G.d[s];
  // X-direction girders along every Y line.
  for (const y of GY) {
    for (let i = 0; i + 1 < GX.length; i++) {
      const [x0, x1] = [GX[i], GX[i + 1]];
      if (!inStorey(s, x1) || !inStorey(s, x0)) continue;
      const body = extrude(
        rect(d / 2, G.w / 2),
        P(x0, y, zTop - d / 2),
        { x: 1, y: 0, z: 0 },
        x1 - x0,
      );
      const clippers = [colId.get(`${s}:${x0}:${y}`), colId.get(`${s}:${x1}:${y}`)].filter(
        Boolean,
      );
      girderIds[s].push(add('girder', `girder X S${s + 1} y${y} bay${i}`, clip(body, clippers)));
    }
  }
  // Y-direction girders along every X line.
  for (const x of GX) {
    if (!inStorey(s, x)) continue;
    for (let j = 0; j + 1 < GY.length; j++) {
      const [y0, y1] = [GY[j], GY[j + 1]];
      const body = extrude(
        rect(d / 2, G.w / 2),
        P(x, y0, zTop - d / 2),
        { x: 0, y: 1, z: 0 },
        y1 - y0,
      );
      const clippers = [colId.get(`${s}:${x}:${y0}`), colId.get(`${s}:${x}:${y1}`)].filter(
        Boolean,
      );
      girderIds[s].push(add('girder', `girder Y S${s + 1} x${x} bay${j}`, clip(body, clippers)));
    }
  }
}

// ── Secondary beams (小梁, Y-direction at bay third-points, girders win) ─────
const beamIds = [[], [], []];
for (let s = 0; s < 3; s++) {
  const zTop = LEVELS_Z[s];
  for (let i = 0; i + 1 < GX.length; i++) {
    const [x0, x1] = [GX[i], GX[i + 1]];
    if (!inStorey(s, x1)) continue;
    for (const t of [1 / 3, 2 / 3]) {
      const bx = x0 + (x1 - x0) * t;
      for (let j = 0; j + 1 < GY.length; j++) {
        const [y0, y1] = [GY[j], GY[j + 1]];
        const body = extrude(
          rect(B.d / 2, B.w / 2),
          P(bx, y0, zTop - B.d / 2),
          { x: 0, y: 1, z: 0 },
          y1 - y0,
        );
        beamIds[s].push(
          add('beam', `beam S${s + 1} x${bx.toFixed(1)} bay${j}`, {
            __clipByGirders: { body, s },
          }),
        );
      }
    }
  }
}

// Resolve beam clips now that all girder ids exist. A beam is deducted by the
// X-girders of its storey (girders take priority over beams); the keep closure
// is idempotent, so listing girders the beam never touches changes nothing.
for (const m of members) {
  if (m.node.__clipByGirders) {
    const { body, s } = m.node.__clipByGirders;
    const xGirders = [];
    for (let k = 0; k < members.length; k++) {
      const g = members[k];
      if (g.kind === 'girder' && g.label.startsWith(`girder X S${s + 1}`)) xGirders.push(g.id);
    }
    m.node = clip(body, xGirders);
  }
}

// ── Slabs (deducted by columns, girders and beams; voids stacked) ───────────
function slab(s) {
  const zTop = LEVELS_Z[s];
  const xMax = s === 2 ? SETBACK_X : GX[GX.length - 1];
  const x0 = GX[0] - MARGIN;
  const x1 = xMax + MARGIN;
  const y0 = GY[0] - MARGIN;
  const y1 = GY[GY.length - 1] + MARGIN;
  const body = extrude(
    rect((y1 - y0) / 2, (x1 - x0) / 2),
    P((x0 + x1) / 2, (y0 + y1) / 2, zTop - SLAB_T),
    { x: 0, y: 0, z: 1 },
    SLAB_T,
  );
  const voids = [];
  const through = (vx0, vx1, vy0, vy1) =>
    boxZ((vx0 + vx1) / 2, (vy0 + vy1) / 2, zTop - SLAB_T - 0.1, vx1 - vx0, vy1 - vy0, SLAB_T + 0.2);
  // Stair + elevator voids on every floor (stacked).
  voids.push(through(STAIR.x0, STAIR.x1, STAIR.y0, STAIR.y1));
  voids.push(through(ELEV.x0, ELEV.x1, ELEV.y0, ELEV.y1));
  // Round MEP sleeves beside the core (two sizes), storeys 1..2 only.
  if (s < 2) {
    voids.push(cylZ(6.6, 9.6, zTop - SLAB_T - 0.1, 0.175, SLAB_T + 0.2));
    voids.push(cylZ(6.6, 8.7, zTop - SLAB_T - 0.1, 0.1, SLAB_T + 0.2));
  }
  const clippers = [...colsOf(s), ...girderIds[s], ...beamIds[s]];
  return clip(openings(body, voids), clippers);
}
add('slab', 'slab L2', slab(0));
add('slab', 'slab L3 (terrace)', slab(1));
add('slab', 'roof slab', slab(2));

// ── Facade wall (front, Y = 0): window grid + entrance, per storey ──────────
const WALL_T = 0.18;
for (let s = 0; s < 3; s++) {
  const z0 = STOREY0[s];
  const h = LEVELS_Z[s] - G.d[s] - z0; // stops under the girder
  const xEnd = s === 2 ? SETBACK_X : GX[GX.length - 1];
  const body = boxZ(xEnd / 2, 0, z0, xEnd, WALL_T, h);
  const voids = [];
  const win = (cx, sill, w, hh) => boxZ(cx, 0, z0 + sill, w, WALL_T + 0.2, hh);
  for (let i = 0; i + 1 < GX.length; i++) {
    const [x0, x1] = [GX[i], GX[i + 1]];
    if (x1 > xEnd + 1e-9) continue;
    const mid = (x0 + x1) / 2;
    if (s === 0 && i === 1) {
      // Entrance bay: double door + sidelights instead of windows.
      voids.push(win(mid, 0.0, 2.6, 3.0));
      voids.push(win(mid - 2.05, 0.45, 0.9, 2.55));
      voids.push(win(mid + 2.05, 0.45, 0.9, 2.55));
    } else {
      // Three windows per bay.
      for (const t of [0.25, 0.5, 0.75]) {
        voids.push(win(x0 + (x1 - x0) * t, s === 0 ? 1.0 : 0.85, 1.35, s === 0 ? 1.8 : 1.6));
      }
    }
  }
  const clippers = GX.filter((x) => x <= xEnd + 1e-9)
    .map((x) => colId.get(`${s}:${x}:0`))
    .filter(Boolean);
  add('wall', `facade wall S${s + 1}`, clip(openings(body, voids), clippers));
}

// ── Gable end walls (west x=0, east x=19.2 / setback) with window pairs ──────
for (let s = 0; s < 3; s++) {
  const z0 = STOREY0[s];
  const h = LEVELS_Z[s] - G.d[s] - z0;
  const ends = s === 2 ? [GX[0], SETBACK_X] : [GX[0], GX[GX.length - 1]];
  for (const [e, x] of ends.entries()) {
    const body = boxZ(x, GY[1], z0, WALL_T, GY[2] - GY[0], h);
    const voids = [];
    for (let j = 0; j + 1 < GY.length; j++) {
      const [y0, y1] = [GY[j], GY[j + 1]];
      for (const t of [0.32, 0.68]) {
        voids.push(boxZ(x, y0 + (y1 - y0) * t, z0 + (s === 0 ? 1.0 : 0.85), WALL_T + 0.2, 1.25, 1.6));
      }
    }
    const clippers = GY.map((y) => colId.get(`${s}:${x}:${y}`)).filter(Boolean);
    add('wall', `gable wall ${e === 0 ? 'W' : 'E'} S${s + 1}`, clip(openings(body, voids), clippers));
  }
}

// ── Core walls (stair + elevator enclosure, storeys 1..3) ────────────────────
for (let s = 0; s < 3; s++) {
  const z0 = STOREY0[s];
  const h = LEVELS_Z[s] - G.d[s] - z0;
  const t = 0.2;
  const cx = (CORE.x0 + CORE.x1) / 2;
  const cy = (CORE.y0 + CORE.y1) / 2;
  const door = (cx_, cy_, alongX) =>
    boxZ(cx_, cy_, z0, alongX ? 1.0 : t + 0.2, alongX ? t + 0.2 : 1.0, 2.1);
  // South wall (faces the floor plate) with the stair door and elevator door.
  add(
    'wall-core',
    `core wall S${s + 1} south`,
    openings(boxZ(cx, CORE.y0, z0, CORE.x1 - CORE.x0, t, h), [
      door(STAIR.x0 + 0.7, CORE.y0, true),
      door((ELEV.x0 + ELEV.x1) / 2, CORE.y0, true),
    ]),
  );
  // North wall, plain.
  add('wall-core', `core wall S${s + 1} north`, boxZ(cx, CORE.y1, z0, CORE.x1 - CORE.x0, t, h));
  // East / west walls.
  add('wall-core', `core wall S${s + 1} west`, boxZ(CORE.x0, cy, z0, t, CORE.y1 - CORE.y0, h));
  add('wall-core', `core wall S${s + 1} east`, boxZ(CORE.x1, cy, z0, t, CORE.y1 - CORE.y0, h));
  // Divider between stair and elevator.
  add(
    'wall-core',
    `core wall S${s + 1} divider`,
    boxZ((STAIR.x1 + ELEV.x0) / 2, cy, z0, t, CORE.y1 - CORE.y0, h),
  );
}

// ── Girder sleeves (MEP through the webs of two front girders) ──────────────
// Re-model: replace two storey-1 front X-girders with sleeved versions is
// complex; instead the long rear girder on line y=5.6 storey 1 gets sleeves.
{
  const s = 0;
  const zTop = LEVELS_Z[s];
  const d = G.d[s];
  const y = GY[1];
  const zC = zTop - d / 2;
  // Find the two rear X-girders of storey 1 on line y=5.6 and add sleeves.
  for (const m of members) {
    if (m.kind !== 'girder') continue;
    if (!m.label.startsWith('girder X S1 y5.6')) continue;
    const base = m.node.Clip ? m.node.Clip.base : m.node;
    const x0 = base.Extrude.origin.x;
    const sleeves = [0.3, 0.55].map((t) =>
      extrude(
        circle(0.1),
        P(x0 + (GX[1] - GX[0]) * t, y - 1.0, zC),
        { x: 0, y: 1, z: 0 },
        2.0,
      ),
    );
    const withSleeves = openings(base, sleeves);
    m.node = m.node.Clip
      ? { Clip: { base: withSleeves, clippers: m.node.Clip.clippers, rule: 'Priority' } }
      : withSleeves;
    m.label += ' (sleeved)';
  }
}

// ── Entrance canopy (cantilever slab over the entrance bay) ──────────────────
add('canopy', 'entrance canopy', boxZ((GX[1] + GX[2]) / 2, -1.1, 3.3, 7.2, 2.6, 0.18));

// ── Parapets ─────────────────────────────────────────────────────────────────
const PAR_T = 0.15;
const PAR_H = 0.65;
function parapet(label, cx, cy, lenX, lenY, zTop) {
  add('parapet', label, boxZ(cx, cy, zTop, lenX, lenY, PAR_H));
}
{
  // Roof parapet (over the set-back block, x 0..12.8).
  const z = LEVELS_Z[2];
  const [x0, x1] = [GX[0] - MARGIN, SETBACK_X + MARGIN];
  const [y0, y1] = [GY[0] - MARGIN, GY[2] + MARGIN];
  parapet('roof parapet S', (x0 + x1) / 2, y0 + PAR_T / 2, x1 - x0, PAR_T, z);
  parapet('roof parapet N', (x0 + x1) / 2, y1 - PAR_T / 2, x1 - x0, PAR_T, z);
  parapet('roof parapet W', x0 + PAR_T / 2, (y0 + y1) / 2, PAR_T, y1 - y0, z);
  parapet('roof parapet E', x1 - PAR_T / 2, (y0 + y1) / 2, PAR_T, y1 - y0, z);
}
{
  // Terrace parapet (level 2 top, around the bay beyond the setback).
  const z = LEVELS_Z[1];
  const [x0, x1] = [SETBACK_X + MARGIN, GX[3] + MARGIN];
  const [y0, y1] = [GY[0] - MARGIN, GY[2] + MARGIN];
  parapet('terrace parapet S', (x0 + x1) / 2, y0 + PAR_T / 2, x1 - x0, PAR_T, z);
  parapet('terrace parapet N', (x0 + x1) / 2, y1 - PAR_T / 2, x1 - x0, PAR_T, z);
  parapet('terrace parapet E', x1 - PAR_T / 2, (y0 + y1) / 2, PAR_T, y1 - y0, z);
}

// ── Steel pergola on the terrace ─────────────────────────────────────────────
{
  const z0 = LEVELS_Z[1];
  const PH = 2.5; // post height
  const posts = [
    [14.2, 1.6],
    [18.2, 1.6],
    [14.2, 9.6],
    [18.2, 9.6],
  ];
  const postIds = posts.map(([x, y], i) =>
    add('steel', `pergola post ${i}`, boxZ(x, y, z0, 0.15, 0.15, PH)),
  );
  // Main H-beams along Y on each post pair (clipped by their posts).
  for (const [i, x] of [14.2, 18.2].entries()) {
    const body = extrude(
      hSection(0.125, 0.1, 0.008, 0.012), // H-250×200 spanning Y
      P(x, 1.6, z0 + PH + 0.125),
      { x: 0, y: 1, z: 0 },
      8.0,
    );
    add('steel', `pergola beam ${i}`, clip(body, [postIds[i], postIds[i + 2]]));
  }
  // Purlins along X resting on top of the beams (no boolean needed).
  for (let k = 0; k < 7; k++) {
    const y = 1.8 + k * 1.3;
    add(
      'steel',
      `pergola purlin ${k}`,
      extrude(hSection(0.06, 0.05, 0.006, 0.008), P(13.7, y, z0 + PH + 0.31), { x: 1, y: 0, z: 0 }, 5.0),
    );
  }
}

// ── Analytic curved panel layer ──────────────────────────────────────────────
// These elements use the kernel's trimmed analytic surface tessellators rather
// than the CSG evaluator. They demonstrate the current curved-panel path:
// barrel-vault cylinder with punched skylights, a spherical lounge dome with
// UV oculi, and a conical entry canopy with slotted trims.
addCurved('curved-roof', 'thick barrel-vault roof with skylights', {
  cylinder: {
    axis_origin: P(0.0, 5.6, 8.95),
    axis_dir: P(1.0, 0.0, 0.0),
    radius: 4.2,
    theta_min: -0.94,
    theta_max: 0.94,
    z_min: -0.35,
    z_max: SETBACK_X + 0.35,
    thickness: 0.16,
    holes: [
      { rectangle: { u_min: -0.22, u_max: 0.22, v_min: 1.2, v_max: 2.75, reverse: true } },
      { rectangle: { u_min: -0.22, u_max: 0.22, v_min: 4.05, v_max: 5.6, reverse: true } },
      { rectangle: { u_min: -0.22, u_max: 0.22, v_min: 6.9, v_max: 8.45, reverse: true } },
      { rectangle: { u_min: -0.22, u_max: 0.22, v_min: 9.75, v_max: 11.3, reverse: true } },
    ],
  },
});

addCurved('curved-dome', 'spherical terrace lounge dome with round UV oculi', {
  sphere: {
    center: P(16.2, 5.6, 7.55),
    radius: 3.0,
    pole: P(0.0, 0.0, 1.0),
    theta_min: 0.25,
    theta_max: 6.03,
    phi_min: 0.18,
    phi_max: 1.23,
    holes: [
      { circle: { center: [1.45, 0.82], radius: 0.12, reverse: true } },
      { circle: { center: [3.15, 0.82], radius: 0.12, reverse: true } },
      { circle: { center: [4.85, 0.82], radius: 0.12, reverse: true } },
    ],
  },
});

addCurved('curved-cone', 'conical entry canopy with slotted trims', {
  cone: {
    apex: P(9.6, -1.25, 2.25),
    axis: P(0.0, 0.0, 1.0),
    half_angle: 0.62,
    theta_min: 0.15,
    theta_max: 4.55,
    height_min: 1.0,
    height_max: 2.35,
    holes: [
      { rectangle: { u_min: 0.85, u_max: 1.25, v_min: 1.2, v_max: 1.85, reverse: true } },
      { rectangle: { u_min: 2.05, u_max: 2.45, v_min: 1.2, v_max: 1.85, reverse: true } },
      { rectangle: { u_min: 3.25, u_max: 3.65, v_min: 1.2, v_max: 1.85, reverse: true } },
    ],
  },
});

// ── Exports ──────────────────────────────────────────────────────────────────
export const KIND_OF = new Map(members.map((m) => [m.id, m.kind]));
export const CURVED_KIND_OF = new Map(curvedPanels.map((m) => [m.id, m.kind]));
export const LABEL_OF = new Map(members.map((m) => [m.id, m.label]));
export const MEMBERS = members;
export const CURVED_PANELS = curvedPanels;
export const BOUNDS = {
  x: [GX[0] - MARGIN, GX[3] + MARGIN],
  y: [GY[0] - 3.2, GY[2] + MARGIN],
  zMax: 13.4,
};
