// Demo building: a two-storey RC frame exercising every kernel feature —
// rectangular & circular columns, clipped girders with circular sleeves,
// slabs with a stair opening and a round duct void, and a wall with windows.
//
// Everything is plain JSON in archi-kernel's serde representation.
// Units: metres, Z up.
//
// Profile-frame conventions (from the kernel's deterministic plane_basis):
//   axis +Z → (u, v) = (+Y, −X): half_w spans Y, half_h spans X
//   axis +X → (u, v) = (+Z, −Y): half_w spans Z, half_h spans Y

const P = (x, y, z) => ({ x, y, z });
const rectProfile = (halfW, halfH) => ({ Rect: { half_w: halfW, half_h: halfH } });
const circleProfile = (r) => ({ Circle: { radius: r } });

const extrude = (profile, origin, axis, length) => ({
  Extrude: { profile, origin, axis, length },
});

const clip = (base, clippers) => ({
  Clip: { base, clippers, rule: 'Priority' },
});

const openings = (base, list) => ({
  OpeningSubtraction: { base, openings: list.map((shape, i) => [i + 1, { shape }]) },
});

// ── Grid ─────────────────────────────────────────────────────────────────────
const GX = [0, 7, 14]; // column lines, X
const GY = [0, 6]; // column lines, Y
const STOREY = 3.5; // floor-to-floor
const SLAB_T = 0.25;
const COL = 0.6; // square column side
const G_W = 0.4; // girder width
const G_D = 0.7; // girder depth
const ROUND_R = 0.35; // round column radius
// The round column replaces the square one at (GX[2], GY[1]).
const ROUND_AT = { x: GX[2], y: GY[1] };

let nextId = 1;
const members = []; // { id, kind, node, label }

function add(kind, label, node) {
  members.push({ id: nextId++, kind, label, node });
  return nextId - 1;
}

// ── Columns (two storeys, continuous per storey) ────────────────────────────
const columnIds = []; // ids per storey for clipping girders
for (let s = 0; s < 2; s++) {
  const z0 = s * STOREY;
  for (const x of GX) {
    for (const y of GY) {
      const round = x === ROUND_AT.x && y === ROUND_AT.y;
      const profile = round ? circleProfile(ROUND_R) : rectProfile(COL / 2, COL / 2);
      const id = add(
        round ? 'column-round' : 'column',
        `${round ? 'round column' : 'column'} (${x},${y}) S${s + 1}`,
        extrude(profile, P(x, y, z0), { x: 0, y: 0, z: 1 }, STOREY),
      );
      columnIds.push({ id, x, y, storey: s });
    }
  }
}

const colsOf = (storey) => columnIds.filter((c) => c.storey === storey);
const colAt = (storey, x, y) =>
  columnIds.find((c) => c.storey === storey && c.x === x && c.y === y);

// ── Girders (top of each storey, clipped by their end columns) ──────────────
// Girder body: top flush with the floor level above (z = (s+1)·STOREY).
function girderX(s, y, x0, x1, withSleeves) {
  const zTop = (s + 1) * STOREY;
  const body = extrude(
    rectProfile(G_D / 2, G_W / 2), // axis +X: half_w spans Z (depth), half_h spans Y (width)
    P(x0, y, zTop - G_D / 2),
    { x: 1, y: 0, z: 0 },
    x1 - x0,
  );
  let node = body;
  if (withSleeves) {
    // Two ø200 sleeves through the web, horizontal (axis +Y), at mid-depth.
    const zc = zTop - G_D / 2;
    const sleeve = (sx) =>
      extrude(circleProfile(0.1), P(sx, y - 1.0, zc), { x: 0, y: 1, z: 0 }, 2.0);
    node = openings(body, [sleeve(x0 + 1.6), sleeve(x0 + 3.1)]);
  }
  const clippers = [colAt(s, x0, y), colAt(s, x1, y)].filter(Boolean).map((c) => c.id);
  return clip(node, clippers);
}

function girderY(s, x, y0, y1) {
  const zTop = (s + 1) * STOREY;
  // axis +Y: by plane_basis the frame spans (Z, X) — symmetric enough to state:
  // half_w spans X? To stay convention-safe we build Y girders as +Y extrusions
  // with a square-ish asymmetry that we verified visually; width along X.
  const body = extrude(
    rectProfile(G_D / 2, G_W / 2),
    P(x, y0, zTop - G_D / 2),
    { x: 0, y: 1, z: 0 },
    y1 - y0,
  );
  const clippers = [colAt(s, x, y0), colAt(s, x, y1)].filter(Boolean).map((c) => c.id);
  return clip(body, clippers);
}

const girderIds = [];
for (let s = 0; s < 2; s++) {
  for (const y of GY) {
    for (let i = 0; i + 1 < GX.length; i++) {
      const withSleeves = s === 0 && y === GY[0] && i === 0; // showcase sleeves once
      girderIds.push(
        add('girder', `girder X S${s + 1} y=${y} bay${i}`, girderX(s, y, GX[i], GX[i + 1], withSleeves)),
      );
    }
  }
  for (const x of GX) {
    girderIds.push(add('girder', `girder Y S${s + 1} x=${x}`, girderY(s, x, GY[0], GY[1])));
  }
}

// ── Slabs (clipped by columns and girders; stair + duct voids on level 1) ───
function slab(s) {
  const zTop = (s + 1) * STOREY;
  const margin = 0.4; // cantilevered edge beyond the column grid
  const x0 = GX[0] - margin;
  const x1 = GX[GX.length - 1] + margin;
  const y0 = GY[0] - margin;
  const y1 = GY[1] + margin;
  // Vertical extrusion (axis +Z): half_w spans Y, half_h spans X.
  const body = extrude(
    rectProfile((y1 - y0) / 2, (x1 - x0) / 2),
    P((x0 + x1) / 2, (y0 + y1) / 2, zTop - SLAB_T),
    { x: 0, y: 0, z: 1 },
    SLAB_T,
  );
  let node = body;
  if (s === 0) {
    const stair = extrude(
      rectProfile(1.4 / 2, 3.0 / 2), // 3.0 (X) × 1.4 (Y) stair well
      P(10.5, 3.0, zTop - SLAB_T - 0.1),
      { x: 0, y: 0, z: 1 },
      SLAB_T + 0.2,
    );
    const duct = extrude(
      circleProfile(0.3),
      P(2.2, 4.4, zTop - SLAB_T - 0.1),
      { x: 0, y: 0, z: 1 },
      SLAB_T + 0.2,
    );
    node = openings(body, [stair, duct]);
  }
  const clippers = [
    ...colsOf(s).map((c) => c.id),
    // Girders of this storey deduct the slab (slab has lowest priority).
    ...girderIds.slice(s * 9, s * 9 + 9),
  ];
  return clip(node, clippers);
}
add('slab', 'slab L2', slab(0));
add('slab', 'roof slab', slab(1));

// ── Wall (storey 1, along Y=0 between the first two columns) ────────────────
{
  const t = 0.18;
  const h = STOREY - G_D; // stops under the girder
  // Vertical extrusion: half_w spans Y (=t/2), half_h spans X (=length/2).
  const body = extrude(
    rectProfile(t / 2, 7 / 2),
    P(3.5, 0, 0),
    { x: 0, y: 0, z: 1 },
    h,
  );
  // Two windows: vertical prisms through the wall thickness (over-thick in Y).
  const win = (cx) =>
    extrude(rectProfile(0.3, 1.6 / 2), P(cx, 0, 0.9), { x: 0, y: 0, z: 1 }, 1.2);
  const node = clip(openings(body, [win(2.0), win(5.0)]), [
    colAt(0, GX[0], 0).id,
    colAt(0, GX[1], 0).id,
  ]);
  add('wall', 'wall S1', node);
}

export const KIND_OF = new Map(members.map((m) => [m.id, m.kind]));
export const LABEL_OF = new Map(members.map((m) => [m.id, m.label]));
export const MEMBERS = members;
export const LEVELS = { storey: STOREY, slabT: SLAB_T };
