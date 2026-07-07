async function runProtein() {
  const seq = document.getElementById('protein-input').value.trim();
  if (!seq) return;

  hide('protein-results');
  hide('protein-error');

  const res = await fetch('/analyze_protein', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ sequence: seq })
  });
  const data = await res.json();

  if (data.error) {
    showError('protein-error', data.error);
    return;
  }

  setText('p-length', data.length + ' aa');
  setText('p-mw', data.mol_weight.toLocaleString() + ' Da');
  setText('p-pi', data.isoelectric);
  setText('p-charge', (data.charge_at_7 > 0 ? '+' : '') + data.charge_at_7);
  setText('p-pos', data.pos_residues);
  setText('p-neg', data.neg_residues);

  const h = data.hydrophobicity;
  const pct = (h + 4.5) / 9 * 100;
  setText('p-hydro-val', h);
  setText('p-hydro-label', data.hydro_label);
  document.getElementById('p-hydro-bar').style.width = pct.toFixed(1) + '%';

  const box = document.getElementById('p-composition');
  box.innerHTML = '';
  Object.entries(data.composition)
    .sort((a, b) => b[1] - a[1])
    .forEach(pair => {
      const tag = document.createElement('div');
      tag.className = 'aa-tag';
      tag.innerHTML = `<span>${pair[0]}</span> × ${pair[1]}`;
      box.appendChild(tag);
    });

  show('protein-results');
}

function clearProtein() {
  document.getElementById('protein-input').value = '';
  hide('protein-results');
  hide('protein-error');
}

function setText(id, val) {
  document.getElementById(id).textContent = val;
}

function show(id) {
  document.getElementById(id).classList.remove('hidden');
}

function hide(id) {
  document.getElementById(id).classList.add('hidden');
}

function showError(id, msg) {
  const el = document.getElementById(id);
  el.textContent = '⚠ ' + msg;
  el.classList.remove('hidden');
}
