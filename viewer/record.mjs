// Records the section-plane sweep as PNG frames for the README GIF.
// Usage: node record.mjs [outdir]   (server must be running on :8741)
//
// Drives the real UI: sets the section slider and dispatches `input`, so every
// frame's vermilion caps are honest kernel `section_all()` output.

import { chromium } from 'playwright';
import { mkdirSync } from 'node:fs';

const outdir = process.argv[2] ?? '/tmp/ak-frames';
mkdirSync(outdir, { recursive: true });

const browser = await chromium.launch();
const page = await browser.newPage({
  viewport: { width: 1100, height: 680 },
  deviceScaleFactor: 1,
});
page.on('console', (m) => {
  const t = m.text();
  if (t.includes('section errors') || t.includes('failed')) console.log('[console]', t.slice(0, 200));
});
await page.goto('http://localhost:8741/index.html?shot=section');
await page.waitForFunction('window.__READY === true', { timeout: 60000 });

// Sweep: pause at the top (full floor plate), glide down through the wall and
// windows, pause at the bottom, then the frame list is ping-ponged by the
// caller. Eased steps so the motion reads as deliberate, not linear.
const Z_TOP = 11.15;
const Z_BOT = 0.55;
const STEPS = 56;
const ease = (t) => 0.5 - 0.5 * Math.cos(Math.PI * t); // smooth in/out

const zs = [];
for (let i = 0; i < 6; i++) zs.push(Z_TOP); // hold on the floor plate
for (let i = 0; i <= STEPS; i++) zs.push(Z_TOP + (Z_BOT - Z_TOP) * ease(i / STEPS));
for (let i = 0; i < 6; i++) zs.push(Z_BOT); // hold at the base

const setZ = async (zv) => {
  await page.evaluate((v) => {
    const slider = document.getElementById('section-z');
    slider.value = String(v);
    slider.dispatchEvent(new Event('input'));
  }, zv);
  // Two settled frames so the rebuilt caps and clipping are on screen.
  await page.evaluate(
    () => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))),
  );
  return page.evaluate(() => (window.__lastSectionErrors ?? []).length);
};

let n = 0;
for (const z of zs) {
  const errs = await setZ(z);
  if (errs) console.log(`frame ${n} z=${z.toFixed(4)}: ${errs} section errors`);
  await page.screenshot({ path: `${outdir}/f${String(n).padStart(3, '0')}.png` });
  n++;
}
console.log(`captured ${n} frames in ${outdir}`);
await browser.close();
