import { chromium } from 'playwright';
const [,, url, file, w = '1600', h = '1000'] = process.argv;
const b = await chromium.launch({ args: ['--use-angle=metal'] });
const p = await b.newPage({ viewport: { width: +w, height: +h }, deviceScaleFactor: 2 });
await p.goto(url);
await p.waitForFunction('window.__READY === true', { timeout: 60000 });
await p.waitForTimeout(600);
await p.screenshot({ path: file });
await b.close();
