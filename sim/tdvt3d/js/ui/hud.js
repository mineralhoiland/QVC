// Tab bar, layer badges, sliders, camera presets, auto-orbit, FPS, play, sim presets, sparkline charts.

const TABS = [
  { id: 'hopf', title: 'Hopf fibration', layers: ['L2'] },
  { id: 'whitehead', title: 'Whitehead core', layers: ['L2'] },
  { id: 'lattice', title: 'BCC lattice', layers: ['L2'] },
  { id: 'horizon', title: 'Acoustic horizon', layers: ['L2'] },
  { id: 'torsion', title: 'Teleparallel torsion', layers: ['L2', 'L3'] },
  { id: 'spectra', title: 'Solitons & spectra', layers: ['L1', 'L2'] },
  { id: 'hypersphere', title: 'S³ flow', layers: ['L2', 'L3'] },
];

export function createHUD({ onTab, onPreset, onOrbit, onPlay }) {
  const root = document.getElementById('hud');
  root.innerHTML = `
    <header class="hud-top">
      <div class="brand">
        <div class="title">TDVT 3D</div>
        <div class="sub">Topological Dynamic Vacuum — Three.js / WebGL</div>
      </div>
      <button type="button" id="playBtn" class="play-btn" aria-pressed="false">Pause</button>
      <nav class="tabs" id="tabs"></nav>
      <div class="layers" id="layers"></div>
    </header>
    <aside class="hud-left">
      <div class="panel" id="sliders"></div>
      <div class="panel cam">
        <button data-preset="iso">iso</button>
        <button data-preset="top">top</button>
        <button data-preset="core">core</button>
        <button id="orbitBtn" class="toggle">auto-orbit</button>
      </div>
      <div class="panel sim-presets" id="simPresets"></div>
      <div class="panel stats">
        <div><span class="k">FPS</span> <span id="fps">—</span></div>
        <div><span class="k">particles</span> <span id="npart">—</span></div>
        <div><span class="k">pixelRatio</span> <span id="pr">—</span></div>
      </div>
      <canvas id="charts" class="hud-charts" width="640" height="120"></canvas>
    </aside>
    <aside class="hud-right">
      <div class="panel" id="numbers"></div>
      <div class="panel" id="dots"></div>
    </aside>
    <div class="banner" id="banner"></div>
    <div class="selftest" id="selftest" hidden></div>
  `;

  let playing = true;
  const playBtn = root.querySelector('#playBtn');
  playBtn.addEventListener('click', () => {
    setPlaying(!playing);
    onPlay?.(playing);
  });

  const tabsEl = root.querySelector('#tabs');
  for (const t of TABS) {
    const b = document.createElement('button');
    b.textContent = t.title;
    b.dataset.tab = t.id;
    b.addEventListener('click', () => onTab(t.id));
    tabsEl.appendChild(b);
  }
  root.querySelectorAll('[data-preset]').forEach((b) => {
    b.addEventListener('click', () => onPreset(b.dataset.preset));
  });
  const orbitBtn = root.querySelector('#orbitBtn');
  orbitBtn.addEventListener('click', () => {
    const on = orbitBtn.classList.toggle('on');
    onOrbit(on);
  });

  function setPlaying(v) {
    playing = !!v;
    playBtn.textContent = playing ? 'Pause' : 'Play';
    playBtn.classList.toggle('paused', !playing);
    playBtn.setAttribute('aria-pressed', playing ? 'false' : 'true');
  }

  function setTab(id) {
    tabsEl.querySelectorAll('button').forEach((b) => b.classList.toggle('on', b.dataset.tab === id));
    const t = TABS.find((x) => x.id === id);
    const layers = root.querySelector('#layers');
    layers.innerHTML = '';
    for (const L of (t?.layers || [])) {
      const s = document.createElement('span');
      s.className = `badge ${L}`;
      s.textContent = L;
      s.title = L === 'L1' ? 'laboratory lock' : L === 'L2' ? 'analogue geometry' : 'spacetime ontology (schematic)';
      layers.appendChild(s);
    }
  }

  function setSimPresets(presets = []) {
    const el = root.querySelector('#simPresets');
    el.innerHTML = '';
    for (const p of presets) {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = p.label;
      b.dataset.id = p.id;
      b.addEventListener('click', () => p.onClick?.());
      el.appendChild(b);
    }
  }

  function setSliders(defs) {
    const el = root.querySelector('#sliders');
    el.innerHTML = '';
    for (const d of defs) {
      const wrap = document.createElement('label');
      wrap.className = 'slider';
      if (d.type === 'checkbox') {
        wrap.classList.toggle('parked', !!d.parked);
        wrap.innerHTML = `<span>${d.label}</span><input type="checkbox" ${d.value ? 'checked' : ''} ${d.parked || d.disabled ? 'disabled' : ''}>`;
        wrap.querySelector('input').addEventListener('change', (e) => d.oninput(e.target.checked));
      } else if (d.type === 'select') {
        const opts = (d.options || []).map((o) =>
          `<option value="${o.value}" ${String(o.value) === String(d.value) ? 'selected' : ''}>${o.label}</option>`).join('');
        wrap.innerHTML = `<span>${d.label}</span><select>${opts}</select><em></em>`;
        wrap.querySelector('select').addEventListener('change', (e) => d.oninput(e.target.value));
      } else {
        wrap.innerHTML = `<span>${d.label}</span><input type="${d.type || 'range'}" min="${d.min}" max="${d.max}" step="${d.step}" value="${d.value}"><em id="sv-${d.id}">${d.fmt ? d.fmt(d.value) : d.value}</em>`;
        const inp = wrap.querySelector('input');
        inp.addEventListener('input', () => {
          const v = parseFloat(inp.value);
          wrap.querySelector('em').textContent = d.fmt ? d.fmt(v) : v;
          d.oninput(v);
        });
      }
      el.appendChild(wrap);
    }
  }

  function setNumbers(rows) {
    const el = root.querySelector('#numbers');
    el.innerHTML = rows.map((r) =>
      `<div class="nv ${r.kind || 'live'}"><span class="k">${r.label}</span><span class="v">${r.value}</span></div>`
    ).join('');
  }

  function setDots(items) {
    const el = root.querySelector('#dots');
    if (!items?.length) { el.innerHTML = ''; return; }
    el.innerHTML = items.map((d) =>
      `<div class="dot ${d.pass ? 'ok' : 'fail'}"><i></i><span>${d.label}</span></div>`
    ).join('');
  }

  function setBanner(text, kind = 'note') {
    const el = root.querySelector('#banner');
    el.className = `banner ${kind} ${text ? '' : 'empty'}`;
    el.textContent = text || '';
  }

  function setStats({ fps, particles, pixelRatio }) {
    root.querySelector('#fps').textContent = fps != null ? fps.toFixed(0) : '—';
    root.querySelector('#npart').textContent = particles != null ? particles.toLocaleString() : '—';
    root.querySelector('#pr').textContent = pixelRatio != null ? pixelRatio.toFixed(2) : '—';
  }

  function setSelftest(html, visible) {
    const el = root.querySelector('#selftest');
    el.hidden = !visible;
    if (html != null) el.innerHTML = html;
  }

  function setCharts(series = []) {
    const canvas = root.querySelector('#charts');
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = 'rgba(6, 10, 18, 0.92)';
    ctx.fillRect(0, 0, w, h);
    if (!series.length) return;
    const n = series.length;
    const pad = 6;
    const rowH = (h - pad * (n + 1)) / n;
    series.forEach((s, i) => {
      const y0 = pad + i * (rowH + pad);
      const data = s.data || s.values || [];
      const color = s.color || '#7ee0ff';
      ctx.fillStyle = 'rgba(255,255,255,0.05)';
      ctx.fillRect(pad, y0, w - pad * 2, rowH);
      if (data.length > 1) {
        let min = Infinity, max = -Infinity;
        for (let k = 0; k < data.length; k++) {
          const v = data[k];
          if (v < min) min = v;
          if (v > max) max = v;
        }
        if (!(max > min)) { min -= 1; max += 1; }
        const span = max - min;
        const left = pad + 2;
        const innerW = w - pad * 2 - 4;
        const innerH = rowH - 14;
        ctx.beginPath();
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.4;
        for (let k = 0; k < data.length; k++) {
          const x = left + (k / (data.length - 1)) * innerW;
          const y = y0 + 12 + (1 - (data[k] - min) / span) * innerH;
          if (k === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        }
        ctx.stroke();
      }
      const last = data.length ? data[data.length - 1] : (s.value ?? 0);
      const unit = s.unit ? ` ${s.unit}` : '';
      ctx.fillStyle = color;
      ctx.font = '10px ui-monospace, SFMono-Regular, Menlo, monospace';
      ctx.textBaseline = 'top';
      ctx.fillText(`${s.label}  ${fmtShort(last)}${unit}`, pad + 4, y0 + 1);
    });
  }

  setPlaying(true);

  return {
    setTab, setSliders, setNumbers, setDots, setBanner, setStats, setSelftest,
    setSimPresets, setCharts, setPlaying, TABS,
  };
}

export function fmtSelftestHTML(result, progressMsg) {
  if (!result) {
    return `<h2>Self-test</h2><p class="busy">${progressMsg || 'Running live compute…'}</p>`;
  }
  const rows = result.rows.map((r) => {
    const g = typeof r.got === 'number' ? r.got.toPrecision(5) : String(r.got);
    const e = typeof r.expected === 'number' ? r.expected.toPrecision(5) : String(r.expected);
    return `<tr class="${r.pass ? 'ok' : 'fail'}"><td>${r.pass ? 'PASS' : 'FAIL'}</td><td>${r.label}</td><td>${g}</td><td>${e}</td><td>${r.extra || ''}</td></tr>`;
  }).join('');
  return `
    <h2>Self-test ${result.allPass ? 'PASS' : 'FAIL'}</h2>
    <p>${result.ms.toFixed(0)} ms live JS. Whitehead 32³ Q=${result.grids[32].Q.toFixed(4)} E=${result.grids[32].E.toFixed(4)}; 48³ Q=${result.grids[48].Q.toFixed(4)}. Richardson ${result.richardson.toFixed(4)}.</p>
    <table><thead><tr><th></th><th>check</th><th>got</th><th>expected</th><th></th></tr></thead><tbody>${rows}</tbody></table>
  `;
}

function fmtShort(v) {
  const n = Number(v);
  if (!Number.isFinite(n)) return '—';
  const a = Math.abs(n);
  if (a === 0) return '0';
  if (a >= 1000 || a < 0.01) return n.toExponential(2);
  return n.toFixed(3);
}
