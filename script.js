// ── PROTEIN ANALYSIS ──
async function runProtein() {
  const sequence = document.getElementById('protein-input').value.trim();
  if (!sequence) return;

  hide('protein-results');
  hide('protein-error');

  const res  = await fetch('/analyze_protein', {
    method:  'POST',
    headers: { 'Content-Type': 'application/json' },
    body:    JSON.stringify({ sequence })
  });
  const data = await res.json();

  if (data.error) { showError('protein-error', data.error); return; }

  setText('p-length', data.length + ' aa');
  setText('p-mw',     data.mol_weight.toLocaleString() + ' Da');
  setText('p-pi',     data.isoelectric);
  setText('p-charge', (data.charge_at_7 > 0 ? '+' : '') + data.charge_at_7);
  setText('p-pos',    data.pos_residues);
  setText('p-neg',    data.neg_residues);

  // Hydrophobicity bar (scale −4.5 to +4.5 → 0% to 100%)
  const hydro = data.hydrophobicity;
  const pct   = ((hydro + 4.5) / 9) * 100;
  setText('p-hydro-val',   hydro);
  setText('p-hydro-label', data.hydro_label);
  document.getElementById('p-hydro-bar').style.width = pct.toFixed(1) + '%';

  // Composition tags
  const compBox = document.getElementById('p-composition');
  compBox.innerHTML = '';
  const sorted = Object.entries(data.composition).sort((a, b) => b[1] - a[1]);
  sorted.forEach(([aa, count]) => {
    const tag = document.createElement('div');
    tag.className = 'aa-tag';
    tag.innerHTML = `<span>${aa}</span> × ${count}`;
    compBox.appendChild(tag);
  });

  show('protein-results');
}

function clearProtein() {
  document.getElementById('protein-input').value = '';
  hide('protein-results');
  hide('protein-error');
}

// ── HELPERS ──
function setText(id, val) { document.getElementById(id).textContent = val; }
function show(id)  { document.getElementById(id).classList.remove('hidden'); }
function hide(id)  { document.getElementById(id).classList.add('hidden'); }
function showError(id, msg) {
  const el = document.getElementById(id);
  el.textContent = '⚠ ' + msg;
  el.classList.remove('hidden');
}
