// archi-kernel viewer — Three.js + wasm.
//
// Visual language: a light "architectural drawing" theme. Warm paper
// background, matte concrete materials with subtle per-kind tones, dark
// feature edges, and a vermilion section plane whose caps and outlines come
// from the kernel's own closed-form section() — not from screen-space tricks.
// Z is up, as in the kernel.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import init, { KernelModel } from './pkg/archi_kernel_wasm.js';
import { MEMBERS, KIND_OF, LEVELS } from './building.js';

THREE.Object3D.DEFAULT_UP.set(0, 0, 1);

const params = new URLSearchParams(location.search);
const SHOT = params.get('shot'); // 'hero' | 'section' | null

// ── Palette ──────────────────────────────────────────────────────────────────
const PAPER = 0xedebe6;
const INK = 0x23262b;
const VERMILION = 0xe04f2e;
const TONES = {
  column: 0xb6b1a7,
  'column-round': 0xa9b3b8,
  girder: 0xa9aeb6,
  slab: 0xc7c3ba,
  wall: 0xd2cec6,
};

// ── Renderer / scene ─────────────────────────────────────────────────────────
const canvas = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.shadowMap.enabled = true;
renderer.shadowMap.type = THREE.PCFSoftShadowMap;
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.05;
renderer.localClippingEnabled = true;

const scene = new THREE.Scene();
scene.background = new THREE.Color(PAPER);
scene.fog = new THREE.Fog(PAPER, 60, 140);

const camera = new THREE.PerspectiveCamera(34, 2, 0.1, 500);
camera.up.set(0, 0, 1);

const controls = new OrbitControls(camera, canvas);
controls.target.set(7, 3, 2.6);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.maxPolarAngle = Math.PI * 0.55;

const pmrem = new THREE.PMREMGenerator(renderer);
scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
scene.environmentIntensity = 0.55;

// Key light: low warm sun from the south-west, long soft shadows.
const sun = new THREE.DirectionalLight(0xfff2e2, 2.6);
sun.position.set(-14, -22, 24);
sun.castShadow = true;
sun.shadow.mapSize.set(4096, 4096);
sun.shadow.camera.left = -16;
sun.shadow.camera.right = 16;
sun.shadow.camera.top = 16;
sun.shadow.camera.bottom = -16;
sun.shadow.camera.far = 80;
sun.shadow.bias = -2e-4;
sun.shadow.normalBias = 0.02;
sun.target.position.set(7, 3, 0);
scene.add(sun, sun.target);
scene.add(new THREE.HemisphereLight(0xdfe8f2, 0xb8b0a2, 0.5));

// Ground: shadow catcher + radially fading grid, drawing-board style.
{
  const ground = new THREE.Mesh(
    new THREE.CircleGeometry(70, 64),
    new THREE.ShadowMaterial({ color: 0x4a4338, opacity: 0.22 }),
  );
  ground.receiveShadow = true;
  scene.add(ground);

  const grid = new THREE.GridHelper(80, 80, 0xb9b3a8, 0xd6d1c7);
  grid.rotation.x = Math.PI / 2; // GridHelper is XZ; rotate into XY (Z-up)
  grid.position.z = -0.002;
  grid.material.transparent = true;
  grid.material.opacity = 0.45;
  scene.add(grid);
}

// ── Kernel ───────────────────────────────────────────────────────────────────
const hud = document.getElementById('hud-stats');
const memberGroup = new THREE.Group();
scene.add(memberGroup);
const sectionGroup = new THREE.Group();
scene.add(sectionGroup);

const clipPlane = new THREE.Plane(new THREE.Vector3(0, 0, -1), 0); // z ≤ h keeps below
let sectionOn = false;
let sectionZ = 1.5;

const edgeMaterial = new THREE.LineBasicMaterial({
  color: INK,
  transparent: true,
  opacity: 0.34,
});
const lineMaterials = []; // LineMaterials needing resolution updates
const meshMaterials = [];

function memberMaterial(kind) {
  const m = new THREE.MeshStandardMaterial({
    color: TONES[kind] ?? 0xbdb9b0,
    roughness: 0.88,
    metalness: 0.0,
  });
  // Clip the shadow casters along with the geometry, or the removed upper
  // storeys would keep casting onto the ground in section mode.
  m.clipShadows = true;
  meshMaterials.push(m);
  return m;
}

let model;
let triangleTotal = 0;

async function build() {
  const t0 = performance.now();
  await init();
  model = new KernelModel();
  for (const m of MEMBERS) model.insert(BigInt(m.id), JSON.stringify(m.node));

  const statuses = JSON.parse(model.evaluate_all());
  const failed = statuses.filter((s) => !s.ok);
  let concrete = 0;
  for (const s of statuses) concrete += s.volume ?? 0;

  for (const m of MEMBERS) {
    let data;
    try {
      data = model.mesh(BigInt(m.id), 0.0015);
    } catch (e) {
      console.warn(`member ${m.id} (${m.label}) failed:`, e);
      continue;
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
    geo.setIndex(new THREE.BufferAttribute(data.indices, 1));
    geo.computeVertexNormals();
    triangleTotal += data.indices.length / 3;

    const mesh = new THREE.Mesh(geo, memberMaterial(KIND_OF.get(m.id)));
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    mesh.userData.id = m.id;
    memberGroup.add(mesh);

    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(geo, 24),
      edgeMaterial,
    );
    mesh.add(edges);
  }

  const dt = performance.now() - t0;
  hud.innerHTML =
    `${MEMBERS.length} members · ${triangleTotal.toLocaleString()} tris · ` +
    `evaluated in ${dt.toFixed(0)} ms<br>` +
    `concrete ${concrete.toFixed(2)} m³` +
    (failed.length ? ` · <span style="color:#c0392b">${failed.length} failed</span>` : '');
  document.getElementById('hud').classList.add('ready');
  window.__model = model; // debug hook

  applyShotMode();
}

// ── Section plane (kernel-computed caps + outlines) ─────────────────────────
function setSection(on, z = sectionZ) {
  sectionOn = on;
  sectionZ = z;
  clipPlane.constant = z;
  for (const m of meshMaterials) m.clippingPlanes = on ? [clipPlane] : null;
  edgeMaterial.clippingPlanes = on ? [clipPlane] : null;
  rebuildSectionGraphics();
}

function rebuildSectionGraphics() {
  sectionGroup.clear();
  lineMaterials.length = 0;
  if (!sectionOn || !model) return;

  const { members: out, errors } = JSON.parse(model.section_all(0, 0, sectionZ, 0, 0, 1));
  if (errors.length) console.warn('section errors:', errors);
  const capMaterial = new THREE.MeshBasicMaterial({
    color: VERMILION,
    side: THREE.DoubleSide,
    polygonOffset: true,
    polygonOffsetFactor: -1,
    polygonOffsetUnits: -2,
  });

  for (const member of out) {
    for (const profile of member.profiles) {
      // Filled cap: THREE.Shape in the (world XY) section plane.
      const toV2 = (pts) => pts.map((p) => new THREE.Vector2(p[0], p[1]));
      const shape = new THREE.Shape(toV2(profile.outer.points));
      for (const hole of profile.holes) {
        shape.holes.push(new THREE.Path(toV2(hole.points)));
      }
      const cap = new THREE.Mesh(new THREE.ShapeGeometry(shape, 24), capMaterial);
      cap.position.z = sectionZ;
      sectionGroup.add(cap);

      // Crisp outline on top of the cap (fat lines).
      for (const loop of [profile.outer, ...profile.holes]) {
        const flat = [];
        for (const p of loop.points) flat.push(p[0], p[1], p[2] + 0.004);
        const first = loop.points[0];
        flat.push(first[0], first[1], first[2] + 0.004);
        const lg = new LineGeometry();
        lg.setPositions(flat);
        const lm = new LineMaterial({
          color: 0x7e2613,
          linewidth: 2.2,
          worldUnits: false,
        });
        lineMaterials.push(lm);
        sectionGroup.add(new Line2(lg, lm));
      }
    }
  }
  syncLineResolution();
}

// ── UI ───────────────────────────────────────────────────────────────────────
const toggle = document.getElementById('section-toggle');
const slider = document.getElementById('section-z');
toggle.addEventListener('change', () => setSection(toggle.checked));
slider.addEventListener('input', () => {
  document.getElementById('section-z-label').textContent =
    `z = ${Number(slider.value).toFixed(2)} m`;
  if (sectionOn) setSection(true, Number(slider.value));
  else sectionZ = Number(slider.value);
});

// ── Shots (deterministic states for screenshots) ─────────────────────────────
function applyShotMode() {
  if (SHOT === 'section') {
    camera.position.set(14.5, -11.0, 10.5);
    controls.target.set(6.2, 3.0, 0.9);
    toggle.checked = true;
    slider.value = '3.3';
    setSection(true, 3.3);
  } else {
    camera.position.set(19.5, -14.0, 9.0);
    controls.target.set(7.0, 3.1, 2.9);
  }
  controls.update();
  if (SHOT) {
    controls.enableDamping = false;
    document.getElementById('hud-hint').style.display = 'none';
  }
  // Signal readiness for headless capture after a few settled frames.
  let frames = 0;
  const tick = () => {
    if (++frames > 8) {
      window.__READY = true;
      return;
    }
    requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
}

// ── Resize / loop ────────────────────────────────────────────────────────────
function syncLineResolution() {
  for (const lm of lineMaterials) {
    lm.resolution.set(canvas.clientWidth * renderer.getPixelRatio(),
      canvas.clientHeight * renderer.getPixelRatio());
  }
}

function resize() {
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (canvas.width !== Math.floor(w * renderer.getPixelRatio()) ||
      canvas.height !== Math.floor(h * renderer.getPixelRatio())) {
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    syncLineResolution();
  }
}

renderer.setAnimationLoop(() => {
  resize();
  controls.update();
  renderer.render(scene, camera);
});

build().catch((e) => {
  hud.textContent = `failed to start: ${e}`;
  console.error(e);
});
