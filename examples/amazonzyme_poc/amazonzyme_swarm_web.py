import os
import json
import urllib.request
import numpy as np
import scipy.spatial.distance
import psycopg
from flask import Flask, render_template_string, jsonify, request

app = Flask(__name__)
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

# Read the base HTML template from swarm_web.py and customize title, controls, and top-scenarios bar
with open("/Users/jonatas/code/pg_bio/examples/swarm_discovery/swarm_web.py", "r") as f:
    src = f.read()

HTML_TEMPLATE = src.split('HTML_TEMPLATE = """')[1].split('"""\n\nimport os')[0]
HTML_TEMPLATE = HTML_TEMPLATE.replace(
    "<title>Swarm Discovery Heatmap</title>",
    "<title>AmazonZyme × pg_bio — Globular Enzyme Discovery Swarm</title>"
).replace(
    "<h2>Protein Correlation Swarm</h2>",
    "<h2>🌿 AmazonZyme × pg_bio Globular Swarm</h2>"
).replace(
    '<span class="subtitle">Pending (gray) vs Explored (blue). Click cell to view 3D-aligned models &amp; full metadata.</span>',
    '<span class="subtitle">All 8 AmazonZyme Organism Scenarios · 3D Globular Gate (Rg &le; 1.45 · pLDDT &ge; 81.5)</span>'
).replace(
    '<option value="40" selected>Latest 40</option>\n                <option value="best-40">Top Best 40 Matches</option>',
    '<option value="best-40" selected>Top Best 40 Matches (All Scenarios)</option>\n                <option value="40">Latest 40</option>'
).replace(
    'visualMap: {\n                    min: 0, max: 0.2,',
    'visualMap: {\n                    min: 0, max: 0.72,'
).replace(
    '<div id="heatmap"></div>',
    '<div id="top-scenarios-bar" style="display:flex; gap:8px; overflow-x:auto; padding:7px 20px; background:#0d1117; border-bottom:1px solid #263238; font-size:0.78em; align-items:center;"></div>\n    <div id="heatmap"></div>'
).replace(
    'var labels = data.nodes.map(n => n.id);',
    '''var labels = data.nodes.map(n => n.id);
                var tsBar = document.getElementById('top-scenarios-bar');
                if (tsBar && data.top_scenarios) {
                    tsBar.innerHTML = '<b style="color:#ffaa00; white-space:nowrap; margin-right:4px;">🏆 Top Match per Scenario:</b>' +
                        data.top_scenarios.map(s => {
                            return `<span onclick="triggerPair('${s.bait_id}','${s.orphan_id}')" style="cursor:pointer; white-space:nowrap; background:#1e293b; border:1px solid #334155; padding:3px 8px; border-radius:4px; color:#e2e8f0;" title="${s.research_line}: ${s.bait_org} -> ${s.orphan_org}">
                                <b style="color:#38bdf8">${s.short_line}</b>: <span style="color:cyan">${s.bait_id}</span>→<span style="color:magenta">${s.orphan_id}</span> <span style="color:#4ade80">(${s.dist})</span>
                            </span>`;
                        }).join('');
                }
                window.triggerPair = function(bId, oId) {
                    let n1 = data.nodes.find(n => n.id === bId);
                    let n2 = data.nodes.find(n => n.id === oId);
                    let i1 = labels.indexOf(bId);
                    let i2 = labels.indexOf(oId);
                    if (n1 && n2 && i1 >= 0 && i2 >= 0) {
                        openModal(bId, oId, data.matrix[i2][i1], n1, n2);
                    }
                };'''
)


def compute_3d_shape_metrics(ca_coords):
    """Returns (rg_angstrom, compactness_ratio, elongation_ratio, fold_class)."""
    n = len(ca_coords)
    if n < 20:
        return 0.0, 9.99, 9.99, "TOO_SHORT"
    centroid = np.mean(ca_coords, axis=0)
    rg = float(np.sqrt(np.mean(np.sum((ca_coords - centroid) ** 2, axis=1))))
    ideal_rg = 2.25 * (n ** 0.38)
    comp_ratio = rg / ideal_rg

    centered = ca_coords - centroid
    cov = np.dot(centered.T, centered) / n
    eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
    elongation = float(np.sqrt(eigvals[0] / max(eigvals[1], 1e-4)))

    if comp_ratio <= 1.45 and elongation <= 3.0:
        fold_class = "GLOBULAR_CATALYTIC_CORE"
    else:
        fold_class = "EXTENDED_COIL_OR_ROD"
    return round(rg, 2), round(comp_ratio, 3), round(elongation, 2), fold_class


def globular_kabsch_superpose_pdb(pdb1_str, pdb2_str):
    """
    Proper right-handed Kabsch 3D superposition (det(R) = +1, no reflection/mirroring):
      - For row-vector coordinates P, Q (N x 3) where Q_aligned = (Q - cQ) @ R + cP,
        H = (Q - cQ).T @ (P - cP) = U @ S @ Vt  ==>  R = U @ diag(1, 1, det(U @ Vt)) @ Vt.
      - Tests both full-length global alignment (0, 0, L) AND multi-secondary-structure
        domain windows (w = min(140, L)) so homologous proteins never flip upside-down.
    """
    def extract_ca_and_plddt(pdb_text):
        coords = []
        plddts = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    b = float(line[60:66])
                    coords.append([x, y, z])
                    plddts.append(b * 100.0 if b <= 1.0 else b)
                except ValueError:
                    pass
        return np.array(coords), (float(np.mean(plddts)) if plddts else 0.0)

    c1, mean_plddt1 = extract_ca_and_plddt(pdb1_str)
    c2, mean_plddt2 = extract_ca_and_plddt(pdb2_str)
    n1, n2 = len(c1), len(c2)
    L = min(n1, n2)
    if L < 10:
        return pdb2_str, 99.9, 0.0, 0, mean_plddt1, mean_plddt2

    def fit_kabsch(P, Q):
        cP = np.mean(P, axis=0)
        cQ = np.mean(Q, axis=0)
        P0 = P - cP
        Q0 = Q - cQ
        H = Q0.T @ P0
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(U @ Vt))
        if d == 0:
            d = 1.0
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        return R, cP, cQ

    d0 = max(1.24 * ((L - 15.0) ** (1.0 / 3.0)) - 1.8, 0.5) if L > 15 else 0.5

    # Always include the full-length global seed (0, 0, L) first, plus domain-sized windows
    w = min(140, L)
    seeds = [(0, 0, L)]
    if n1 > w or n2 > w:
        step1 = max(20, (n1 - w) // 4) if n1 > w else 1
        step2 = max(20, (n2 - w) // 4) if n2 > w else 1
        for i1 in range(0, n1 - w + 1, step1):
            for i2 in range(max(0, i1 - 60), min(n2 - w + 1, i1 + 61), max(15, step2 // 2)):
                seeds.append((i1, i2, w))

    best_tm = -1.0
    best_transform = (np.eye(3), np.zeros(3), np.zeros(3))
    best_stats = (9.99, 0.0, L)

    for i1, i2, win in seeds:
        R, cP, cQ = fit_kabsch(c1[i1:i1 + win], c2[i2:i2 + win])
        for cutoff in (8.0, 5.0, 3.5):
            c2_al = (c2 - cQ) @ R + cP
            D = scipy.spatial.distance.cdist(c1, c2_al)
            idx2 = np.argmin(D, axis=1)
            min_d = D[np.arange(n1), idx2]
            mask = min_d < cutoff
            if int(np.sum(mask)) >= max(20, int(0.20 * L)):
                R, cP, cQ = fit_kabsch(c1[mask], c2[idx2[mask]])

        c2_al = (c2 - cQ) @ R + cP
        D = scipy.spatial.distance.cdist(c1, c2_al)
        min_d = np.min(D, axis=1)
        tm = float(np.sum(1.0 / (1.0 + (min_d / d0) ** 2)) / L)
        core_mask = min_d < 4.0
        if int(np.sum(core_mask)) >= 15:
            n_core = int(np.sum(core_mask))
            core_rmsd = float(np.sqrt(np.mean(min_d[core_mask] ** 2)))
        else:
            top_k = np.sort(min_d)[:max(15, int(0.5 * L))]
            n_core = int(len(top_k))
            core_rmsd = float(np.sqrt(np.mean(top_k ** 2)))

        if tm > best_tm:
            best_tm = tm
            best_transform = (R, cP, cQ)
            best_stats = (core_rmsd, tm, n_core)

    best_R, best_cP, best_cQ = best_transform
    best_rmsd, best_tm, best_n_core = best_stats

    out_lines = []
    for line in pdb2_str.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                pt_new = (xyz - best_cQ) @ best_R + best_cP
                new_line = f"{line[:30]}{pt_new[0]:8.3f}{pt_new[1]:8.3f}{pt_new[2]:8.3f}{line[54:]}"
                out_lines.append(new_line)
            except ValueError:
                out_lines.append(line)
        else:
            out_lines.append(line)

    return "\n".join(out_lines), best_rmsd, best_tm, best_n_core, mean_plddt1, mean_plddt2


@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)


def fetch_pdb_string(uniprot_id):
    if uniprot_id == "A0A8E4UR43" and os.path.exists("/Users/jonatas/code/pg_bio/examples/amazonzyme_poc/8V16.pdb"):
        with open("/Users/jonatas/code/pg_bio/examples/amazonzyme_poc/8V16.pdb", "r") as f:
            return f.read(), 200
    path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{uniprot_id}-F1-model_v4.pdb"
    if os.path.exists(path):
        with open(path, 'r') as f:
            content = f.read()
            if "ATOM" in content:
                return content, 200

    for ver in (6, 4):
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{ver}.pdb"
            pdb_data = urllib.request.urlopen(url, timeout=10).read().decode('utf-8')
            if "ATOM" in pdb_data:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, 'w') as f:
                    f.write(pdb_data)
                return pdb_data, 200
        except Exception:
            pass

    try:
        with psycopg.connect(DB_URI) as conn:
            with conn.cursor() as cur:
                cur.execute("SELECT sequence FROM proteins WHERE uniprot_id = %s", (uniprot_id,))
                res = cur.fetchone()
                if res and res[0]:
                    seq = res[0]
                    cur.execute("SELECT bio_fold_sequence(%s)", (seq[:380],))
                    pdb_data = cur.fetchone()[0]
                    if pdb_data and "ATOM" in pdb_data:
                        os.makedirs(os.path.dirname(path), exist_ok=True)
                        with open(path, 'w') as f:
                            f.write(pdb_data)
                        return pdb_data, 200
    except Exception as e:
        return f"Error folding {uniprot_id}: {e}", 500
    return "Structure not found", 404


@app.route('/pdb/<uniprot_id>')
def get_pdb(uniprot_id):
    return fetch_pdb_string(uniprot_id)


@app.route('/api/compare/<id1>/<id2>')
def compare_proteins(id1, id2):
    pdb1, s1 = fetch_pdb_string(id1)
    pdb2, s2 = fetch_pdb_string(id2)
    if s1 != 200 or s2 != 200:
        return jsonify({'error': f'Could not load 3D structure ({id1}: {s1}, {id2}: {s2})'}), 500

    pdb2_aligned, core_rmsd, tm_score, n_core, plddt1, plddt2 = globular_kabsch_superpose_pdb(pdb1, pdb2)

    def extract_ca(pdb_text):
        return np.array([
            [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            for l in pdb_text.splitlines()
            if l.startswith("ATOM  ") and l[12:16].strip() == "CA"
        ])

    rg1, comp1, elong1, fclass1 = compute_3d_shape_metrics(extract_ca(pdb1))
    rg2, comp2, elong2, fclass2 = compute_3d_shape_metrics(extract_ca(pdb2))

    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            def build_meta(uid, plddt, rg, comp, fclass):
                cur.execute("""
                    SELECT name, sequence, molecular_weight(sequence), sequence_entropy(sequence),
                           isoelectric_point(sequence), gravy_score(sequence)
                    FROM proteins WHERE uniprot_id = %s
                """, (uid,))
                row = cur.fetchone()
                name, seq, mw, ent, pi, gravy = row if row else ("Unknown", "", 0.0, 0.0, 7.0, 0.0)
                L = max(len(seq), 1)
                acidic = (seq.count('D') + seq.count('E')) / L * 100.0
                basic = (seq.count('K') + seq.count('R') + seq.count('H')) / L * 100.0
                hydro = sum(seq.count(c) for c in 'AILMFWV') / L * 100.0

                cur.execute("""
                    SELECT research_line, reactor_profile, wetlab_solubility, orphan_organism
                    FROM amazonzyme_swarm_discoveries WHERE orphan_id = %s LIMIT 1
                """, (uid,))
                az = cur.fetchone()
                domains = [f"3D Fold: {fclass} (Rg={rg}Å, Ratio={comp})"]
                organism = az[3] if az and az[3] else ""
                if az:
                    domains.append(f"Line: {az[0]}")
                    domains.append(f"pI: {pi:.2f} ({az[1]})")
                    domains.append(f"GRAVY: {gravy:.3f}")
                else:
                    domains.append(f"pI: {pi:.2f} | GRAVY: {gravy:.3f}")

                return {
                    'id': uid,
                    'name': name,
                    'organism': organism,
                    'length': len(seq),
                    'mw_kda': f"{mw/1000.0:.1f}",
                    'entropy': f"{ent:.2f}",
                    'acidic_pct': f"{acidic:.1f}",
                    'basic_pct': f"{basic:.1f}",
                    'hydro_pct': f"{hydro:.1f}",
                    'mean_plddt': f"{plddt:.1f}",
                    'domains': domains,
                    'seq': seq
                }

            m1 = build_meta(id1, plddt1, rg1, comp1, fclass1)
            m2 = build_meta(id2, plddt2, rg2, comp2, fclass2)

            cur.execute(
                "SELECT sequence_alignment_score(%s, %s), embedding_cosine_distance(get_esm_embedding(%s), get_esm_embedding(%s))",
                (m1['seq'], m2['seq'], m1['seq'], m2['seq'])
            )
            sw_score, esm_kmer_dist = cur.fetchone()
            m1.pop('seq', None)
            m2.pop('seq', None)

    return jsonify({
        'pdb1': pdb1,
        'pdb2_aligned': pdb2_aligned,
        'meta1': m1,
        'meta2': m2,
        'alignment': {
            'core_rmsd': f"{core_rmsd:.2f}",
            'tm_score': f"{tm_score:.3f}",
            'aligned_residues': n_core,
            'sw_score': sw_score,
            'esm_kmer_dist': f"{esm_kmer_dist:.4f}"
        }
    })


@app.route('/api/data')
def get_data():
    limit_param = request.args.get('limit', 'best-40')
    n1 = request.args.get('n1')
    n2 = request.args.get('n2')

    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            # 1. Load metadata for all nodes in amazonzyme_swarm_discoveries
            cur.execute("""
                SELECT orphan_id, status, research_line, orphan_organism, isoelectric_pi,
                       gravy_hydropathy, mean_plddt, rg_angstrom, compactness_ratio
                FROM amazonzyme_swarm_discoveries
            """)
            all_meta = {
                r[0]: {
                    'status': r[1],
                    'raw_line': r[2],
                    'family': f"{r[2]} | pI={r[4]:.2f} | GRAVY={r[5]:.2f} | pLDDT={r[6]:.1f} | Rg={r[7]:.1f}A ({r[8]:.2f})",
                    'organism': r[3]
                }
                for r in cur.fetchall()
            }

            # 2. Fetch the #1 best match per research_line for the Top Scenarios bar
            cur.execute("""
                WITH best_per_line AS (
                    SELECT
                        research_line, bait_id, bait_organism, orphan_id, orphan_organism,
                        vector_distance, hybrid_score,
                        ROW_NUMBER() OVER (PARTITION BY research_line ORDER BY vector_distance ASC, hybrid_score ASC) AS rn
                    FROM amazonzyme_swarm_discoveries
                    WHERE bait_id != orphan_id
                )
                SELECT research_line, bait_id, bait_organism, orphan_id, orphan_organism, vector_distance
                FROM best_per_line
                WHERE rn = 1
                ORDER BY vector_distance ASC;
            """)
            top_scenarios = []
            for r_line, b_id, b_org, o_id, o_org, v_dist in cur.fetchall():
                short_line = r_line.split(" (")[0]
                top_scenarios.append({
                    'research_line': r_line,
                    'short_line': short_line,
                    'bait_id': b_id,
                    'bait_org': b_org,
                    'orphan_id': o_id,
                    'orphan_org': o_org,
                    'dist': f"{v_dist:.4f}"
                })

            # 3. Build nodes_info according to limit_param, ALWAYS pairing (bait_id, orphan_id)
            nodes_info = {}
            if limit_param.lower() == 'best-40':
                # Stratified selection across all 8 research lines so every scenario's top matches + baits are included
                cur.execute("""
                    WITH ranked_edges AS (
                        SELECT
                            id, research_line, bait_id, bait_organism, orphan_id, orphan_organism,
                            status, vector_distance, hybrid_score,
                            ROW_NUMBER() OVER (
                                PARTITION BY research_line
                                ORDER BY vector_distance ASC, hybrid_score ASC, id ASC
                            ) AS scenario_rank
                        FROM amazonzyme_swarm_discoveries
                        WHERE bait_id != orphan_id
                    )
                    SELECT bait_id, orphan_id, research_line, bait_organism, orphan_organism, status
                    FROM ranked_edges
                    ORDER BY scenario_rank ASC, vector_distance ASC
                    LIMIT 80;
                """)
                edge_rows = cur.fetchall()
                for b_id, o_id, r_line, b_org, o_org, st in edge_rows:
                    if b_id not in nodes_info and len(nodes_info) < 40:
                        nodes_info[b_id] = all_meta.get(b_id, {
                            'status': 'Explored',
                            'raw_line': r_line,
                            'family': r_line,
                            'organism': b_org
                        })
                    if o_id not in nodes_info and len(nodes_info) < 40:
                        nodes_info[o_id] = all_meta.get(o_id, {
                            'status': st,
                            'raw_line': r_line,
                            'family': r_line,
                            'organism': o_org
                        })
                    if len(nodes_info) >= 40:
                        break
            else:
                if limit_param.lower() == 'all':
                    max_nodes = 800
                else:
                    try:
                        max_nodes = int(limit_param)
                    except Exception:
                        max_nodes = 40

                cur.execute("""
                    SELECT bait_id, orphan_id, research_line, bait_organism, orphan_organism, status
                    FROM amazonzyme_swarm_discoveries
                    WHERE bait_id != orphan_id
                    ORDER BY id DESC
                    LIMIT %s;
                """, (max_nodes * 2,))
                for b_id, o_id, r_line, b_org, o_org, st in cur.fetchall():
                    if b_id not in nodes_info and len(nodes_info) < max_nodes:
                        nodes_info[b_id] = all_meta.get(b_id, {
                            'status': 'Explored',
                            'raw_line': r_line,
                            'family': r_line,
                            'organism': b_org
                        })
                    if o_id not in nodes_info and len(nodes_info) < max_nodes:
                        nodes_info[o_id] = all_meta.get(o_id, {
                            'status': st,
                            'raw_line': r_line,
                            'family': r_line,
                            'organism': o_org
                        })
                    if len(nodes_info) >= max_nodes:
                        break

            # Ensure URL-requested n1 and n2 are always included in the matrix
            for extra_id in (n1, n2):
                if extra_id and extra_id != 'undefined' and extra_id not in nodes_info and extra_id in all_meta:
                    nodes_info[extra_id] = all_meta[extra_id]

            node_ids = list(nodes_info.keys())
            if not node_ids:
                return jsonify({'nodes': [], 'matrix': [], 'top_scenarios': []})

            cur.execute("""
                SELECT uniprot_id, name, length(sequence), embedding
                FROM proteins
                WHERE uniprot_id = ANY(%s)
            """, (node_ids,))
            prot_rows = cur.fetchall()

    emb_dict = {}
    nodes = []
    for r in prot_rows:
        uid, name, length, emb_str = r
        emb_dict[uid] = np.array(eval(emb_str))
        info = nodes_info[uid]
        nodes.append({
            'id': uid,
            'status': info['status'],
            'raw_line': info.get('raw_line', ''),
            'family': info['family'],
            'organism': info['organism'],
            'name': name,
            'length': length
        })

    # Sort by research_line then id so each AmazonZyme organism scenario forms a coherent diagonal block
    nodes.sort(key=lambda x: (x.get('raw_line', ''), x['id']))
    emb_matrix = np.array([emb_dict[n['id']] for n in nodes])
    dist_matrix = scipy.spatial.distance.cdist(emb_matrix, emb_matrix, metric='cosine')
    np.fill_diagonal(dist_matrix, 0.0)

    return jsonify({
        'nodes': nodes,
        'matrix': dist_matrix.tolist(),
        'top_scenarios': top_scenarios
    })


if __name__ == '__main__':
    app.run(port=8085)
