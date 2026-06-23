from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

# ─────────────────────────────────────────
#  PROTEIN ANALYZER LOGIC
# ─────────────────────────────────────────

# Residue weights (g/mol) — derived from atomic composition using IUPAC 2021
# atomic weights. Each value = free amino acid MW - H2O (18.015 g/mol).
# Free amino acid weights are computed from molecular formulas, not looked up,
# so these are as accurate as the underlying atomic weights allow.
AA_WEIGHTS = {
    'A':  71.079, 'R': 156.189, 'N': 114.104, 'D': 115.088, 'C': 103.139,
    'E': 129.115, 'Q': 128.131, 'G':  57.052, 'H': 137.142, 'I': 113.160,
    'L': 113.160, 'K': 128.175, 'M': 131.193, 'F': 147.177, 'P':  97.117,
    'S':  87.078, 'T': 101.105, 'W': 186.214, 'Y': 163.176, 'V':  99.133,
}

# pKa values — Lehninger Principles of Biochemistry
# Arg pKa set to 12.0 to match standard bioinformatics tool conventions
# (e.g. Peptide2.0, ExPASy), rather than Lehninger's 12.5.
PKA_NTERM   = 8.0    # N-terminus alpha-amino group
PKA_CTERM   = 3.1    # C-terminus alpha-carboxyl group
AA_PKA_SIDE = {
    'D': 3.9,   # Aspartate side chain
    'E': 4.1,   # Glutamate side chain
    'H': 6.0,   # Histidine side chain
    'C': 8.3,   # Cysteine side chain
    'Y': 10.1,  # Tyrosine side chain
    'K': 10.5,  # Lysine side chain
    'R': 12.0,  # Arginine side chain (convention: 12.0, not Lehninger's 12.5)
}

# Kyte-Doolittle hydrophobicity scale (1982)
AA_HYDROPHOBICITY = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

def calculate_charge_at_ph(seq, ph):
    """
    Net charge at a given pH using Henderson-Hasselbalch.
    Basic groups contribute positive charge below their pKa.
    Acidic groups contribute negative charge above their pKa.
    """
    charge = 0.0
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))   # N-terminus
    charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))   # C-terminus
    for aa, pka in AA_PKA_SIDE.items():
        count = seq.count(aa)
        if count == 0:
            continue
        if aa in ('D', 'E', 'C', 'Y'):
            charge -= count / (1.0 + 10 ** (pka - ph))
        else:
            charge += count / (1.0 + 10 ** (ph - pka))
    return charge

def calculate_pi(seq):
    """
    Find pI by binary search over pH 0–14.
    pI is the pH where net charge = 0.
    """
    ph_low, ph_high = 0.0, 14.0
    for _ in range(1000):
        ph_mid = (ph_low + ph_high) / 2
        charge = calculate_charge_at_ph(seq, ph_mid)
        if abs(charge) < 0.0001:
            break
        if charge > 0:
            ph_low = ph_mid
        else:
            ph_high = ph_mid
    return round(ph_mid, 2)

def analyze_protein(sequence):
    seq = sequence.upper().strip()

    invalid = [ch for ch in seq if ch not in AA_WEIGHTS]
    if invalid:
        return {"error": f"Invalid characters: {set(invalid)}. Use standard single-letter amino acid codes."}

    n = len(seq)

    # Molecular weight: sum residue weights + one water (18.015) for N/C termini
    mw = sum(AA_WEIGHTS[aa] for aa in seq) + 18.015

    # Isoelectric point
    pi = calculate_pi(seq)

    # Net charge at pH 7
    charge_at_7 = round(calculate_charge_at_ph(seq, 7.0), 2)

    # Hydrophobicity
    hydro = sum(AA_HYDROPHOBICITY[aa] for aa in seq) / n
    if hydro > 1.5:
        hydro_label = "Highly hydrophobic — likely membrane-associated"
    elif hydro > 0:
        hydro_label = "Moderately hydrophobic"
    elif hydro >= -1:
        hydro_label = "Slightly hydrophilic"
    else:
        hydro_label = "Highly hydrophilic — likely soluble/cytoplasmic"

    # Amino acid composition
    composition = {}
    for aa in seq:
        composition[aa] = composition.get(aa, 0) + 1

    return {
        "length":         n,
        "mol_weight":     round(mw, 2),
        "isoelectric":    pi,
        "charge_at_7":    charge_at_7,
        "pos_residues":   seq.count('K') + seq.count('R') + seq.count('H'),
        "neg_residues":   seq.count('D') + seq.count('E'),
        "hydrophobicity": round(hydro, 3),
        "hydro_label":    hydro_label,
        "composition":    composition,
    }


# ─────────────────────────────────────────
#  FLASK ROUTES
# ─────────────────────────────────────────

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze_protein', methods=['POST'])
def protein_route():
    data = request.get_json()
    result = analyze_protein(data.get('sequence', ''))
    return jsonify(result)


if __name__ == '__main__':
    app.run(debug=True)
