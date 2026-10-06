from flask import Flask, render_template_string, jsonify
import psycopg
import numpy as np
from scipy.spatial.distance import cdist

app = Flask(__name__)
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Swarm Discovery Heatmap</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.0/dist/echarts.min.js"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        body { background-color: #121212; color: #fff; font-family: sans-serif; text-align: center; margin: 0; overflow: hidden; }
        #heatmap { width: 100vw; height: 100vh; }
        .info { position: absolute; top: 10px; left: 10px; z-index: 10; text-align: left; background: rgba(0,0,0,0.7); padding: 15px; border-radius: 8px; pointer-events: none; }
        
        #modal { display: none; position: fixed; top: 10%; left: 10%; width: 80%; height: 80%; background: #1e1e1e; border: 2px solid #4facfe; box-shadow: 0 0 20px #000; z-index: 100; border-radius: 10px; overflow: hidden; }
        .modal-header { display: flex; justify-content: space-between; padding: 10px 20px; background: #000; font-weight: bold; font-size: 1.2em; }
        .close-btn { cursor: pointer; color: #ff5555; }
        .viewers-container { display: flex; width: 100%; height: calc(100% - 45px); }
        .viewer-box { width: 50%; height: 100%; position: relative; }
        .viewer-title { position: absolute; top: 10px; left: 10px; z-index: 10; background: rgba(0,0,0,0.8); padding: 8px 12px; border-radius: 5px; }
        .viewer-title a { text-decoration: none; }
        .viewer-title a:hover { text-decoration: underline; }
        .controls-bar { margin-top: 5px; font-size: 0.8em; }
        .btn { background: #333; color: #fff; border: 1px solid #555; padding: 3px 8px; cursor: pointer; border-radius: 3px; }
        .btn:hover { background: #555; }
        
        /* Custom tooltip for axis labels */
        #axis-tooltip { display: none; position: absolute; background: rgba(0,0,0,0.9); border: 1px solid #4facfe; padding: 10px; border-radius: 5px; z-index: 50; pointer-events: none; font-size: 12px; text-align: left; max-width: 300px; }
    </style>
</head>
<body>
    <div class="info">
        <h2>Protein Correlation Swarm</h2>
        <p>Pending (gray) vs Explored (blue).</p>
        <p style="color: #4facfe;"><b>Click a correlation cell to view synced 3D models!</b></p>
        <p style="color: #ffaa00; font-size: 0.9em;">Hover over axis labels for protein details.</p>
    </div>
    <div id="heatmap"></div>
    <div id="axis-tooltip"></div>

    <div id="modal">
        <div class="modal-header">
            <div>
                <span id="modal-title">Comparing Proteins</span>
                <div class="controls-bar">
                    <button class="btn" onclick="applyDefaultColoring()">Color by Identity</button>
                    <button class="btn" onclick="applyConfidenceColoring()">Color by AI Confidence</button>
                    <span style="border-left: 1px solid #555; margin: 0 5px;"></span>
                    <button class="btn" onclick="toggleMergeView()">🧬 Merge View (Superposition)</button>
                    <button id="pocket-btn" class="btn" onclick="togglePocketSelection()">🔍 Spatial Pocket Selection (pg_bio)</button>
                </div>
            </div>
            <span class="close-btn" onclick="closeModal()">✖ Close</span>
        </div>
        <div class="viewers-container">
            <div class="viewer-box">
                <div id="v1-title" class="viewer-title"></div>
                <div id="viewer1" style="width: 100%; height: 100%; position: relative;"></div>
            </div>
            <div class="viewer-box">
                <div id="v2-title" class="viewer-title"></div>
                <div id="viewer2" style="width: 100%; height: 100%; position: relative;"></div>
            </div>
        </div>
    </div>

    <script>
        var chart = echarts.init(document.getElementById('heatmap'));
        var syncId = null;
        var master = null;
        var globalV1 = null;
        var globalV2 = null;

        function closeModal() {
            document.getElementById('modal').style.display = 'none';
            if (syncId) cancelAnimationFrame(syncId);
        }
        
        function applyDefaultColoring() {
            if(!globalV1 || !globalV2) return;
            globalV1.setStyle({}, {cartoon: {color: "cyan"}});
            globalV2.setStyle({}, {cartoon: {color: "magenta"}});
            globalV1.render();
            globalV2.render();
        }
        
        function getpLDDTColor(atom) {
            // ESMFold stores pLDDT as 0.0-1.0. AlphaFold uses 0.0-100.0.
            let b = atom.b <= 1.0 ? atom.b * 100 : atom.b;
            if (b > 90) return '#0053d6'; // Dark blue (Very high)
            if (b > 70) return '#65cbff'; // Light blue (Confident)
            if (b > 50) return '#ffdb13'; // Yellow (Low)
            return '#ff7d45';             // Orange/Red (Very low)
        }
        
        function applyConfidenceColoring() {
            if(!globalV1 || !globalV2) return;
            globalV1.setStyle({}, {cartoon: {colorfunc: getpLDDTColor}});
            globalV2.setStyle({}, {cartoon: {colorfunc: getpLDDTColor}});
            globalV1.render();
            globalV2.render();
        }

        fetch('/api/data').then(r => r.json()).then(data => {
            var labels = data.nodes.map(n => n.id);
            var heatmapData = [];
            
            for(let i=0; i<data.nodes.length; i++) {
                for(let j=0; j<data.nodes.length; j++) {
                    heatmapData.push([j, i, data.matrix[i][j]]);
                }
            }
            
            var option = {
                tooltip: { 
                    position: 'top', 
                    formatter: function(p) { 
                        let n1 = data.nodes[p.value[0]]; // Y axis
                        let n2 = data.nodes[p.value[1]]; // X axis
                        return `<div style="text-align: left;">
                                <b style="color: #4facfe; font-size: 1.1em;">Vector Distance: ${p.value[2].toFixed(4)}</b><hr style="border-color:#333;">
                                <b style="color: cyan;">Node X: ${n2.id} (${n2.status})</b><br>
                                Family: ${n2.family}<br>Organism: ${n2.organism}<br>Protein: ${n2.name}<br>Length: ${n2.length} aa<br><br>
                                <b style="color: magenta;">Node Y: ${n1.id} (${n1.status})</b><br>
                                Family: ${n1.family}<br>Organism: ${n1.organism}<br>Protein: ${n1.name}<br>Length: ${n1.length} aa
                                </div>`;
                    } 
                },
                grid: { height: '80%', top: '10%' },
                xAxis: { 
                    type: 'category', data: labels, triggerEvent: true,
                    axisLabel: { 
                        formatter: function(val) {
                            var node = data.nodes.find(n => n.id === val);
                            return node.status === 'Pending' ? '{pending|' + val + '}' : '{explored|' + val + '}';
                        },
                        rich: { pending: { color: '#888' }, explored: { color: '#4facfe', fontWeight: 'bold' } },
                        rotate: 90, fontSize: 10
                    }
                },
                yAxis: { 
                    type: 'category', data: labels, triggerEvent: true,
                    axisLabel: {
                        formatter: function(val) {
                            var node = data.nodes.find(n => n.id === val);
                            return node.status === 'Pending' ? '{pending|' + val + '}' : '{explored|' + val + '}';
                        },
                        rich: { pending: { color: '#888' }, explored: { color: '#4facfe', fontWeight: 'bold' } },
                        fontSize: 10
                    }
                },
                visualMap: {
                    min: 0, max: 0.2, calculable: true, orient: 'vertical', right: '2%', top: 'center',
                    precision: 3,
                    textStyle: { color: '#ffffff' },
                    inRange: { color: ['#00fa9a', '#4facfe', '#444466'] }
                },
                series: [{
                    name: 'Distance', type: 'heatmap', data: heatmapData,
                    emphasis: { itemStyle: { shadowBlur: 10, shadowColor: 'rgba(0, 0, 0, 0.5)' } }
                }]
            };
            chart.setOption(option);
            
            // Custom Tooltip for Axis Labels
            var axisTooltip = document.getElementById('axis-tooltip');
            
            // Track the active node from mouseover to be used by mousemove
            var hoveredNode = null;

            chart.on('mouseover', function(params) {
                if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                    hoveredNode = data.nodes.find(n => n.id === params.value);
                    if (hoveredNode) {
                        axisTooltip.innerHTML = `<b style="color:#4facfe;">${hoveredNode.id} (${hoveredNode.status})</b><br>
                                                 <b>Family:</b> ${hoveredNode.family}<br>
                                                 <b>Organism:</b> ${hoveredNode.organism}<br>
                                                 <b>Name:</b> ${hoveredNode.name}<br>
                                                 <b>Length:</b> ${hoveredNode.length} aa`;
                        axisTooltip.style.display = 'block';
                    }
                }
            });

            // Make tooltip follow the mouse and stay within bounds
            document.addEventListener('mousemove', function(e) {
                if (hoveredNode && axisTooltip.style.display === 'block') {
                    let tooltipRect = axisTooltip.getBoundingClientRect();
                    let x = e.clientX + 15;
                    let y = e.clientY + 15;

                    // Prevent clipping on the right edge
                    if (x + tooltipRect.width > window.innerWidth) {
                        x = e.clientX - tooltipRect.width - 10;
                    }
                    // Prevent clipping on the bottom edge
                    if (y + tooltipRect.height > window.innerHeight) {
                        y = e.clientY - tooltipRect.height - 10;
                    }

                    axisTooltip.style.left = x + 'px';
                    axisTooltip.style.top = y + 'px';
                }
            });

            chart.on('mouseout', function(params) {
                if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                    axisTooltip.style.display = 'none';
                    hoveredNode = null;
                }
            });
            
            chart.on('click', function(params) {
                if (params.componentType !== 'series') return;
                let id1 = labels[params.value[0]];
                let id2 = labels[params.value[1]];
                
                let node1 = data.nodes.find(n => n.id === id1);
                let node2 = data.nodes.find(n => n.id === id2);
                
                document.getElementById('modal-title').innerText = `Comparing Distance: ${params.value[2].toFixed(4)}`;
                
                document.getElementById('v1-title').innerHTML = `<b style="font-size:1.1em;">${node1.organism}</b><br>
                    <a href="https://www.uniprot.org/uniprotkb/${node1.id}/entry" target="_blank" style="color:cyan; font-weight:bold;">${node1.id} 🔗</a> 
                    <a href="https://alphafold.ebi.ac.uk/entry/${node1.id}" target="_blank" style="font-size:0.7em; color:#aaa;">(AF DB)</a><br>
                    <span style="font-size:0.85em">${node1.name}</span><br>
                    <span style="font-size:0.8em; color:#aaa">${node1.length} aa</span>`;
                    
                document.getElementById('v2-title').innerHTML = `<b style="font-size:1.1em;">${node2.organism}</b><br>
                    <a href="https://www.uniprot.org/uniprotkb/${node2.id}/entry" target="_blank" style="color:magenta; font-weight:bold;">${node2.id} 🔗</a> 
                    <a href="https://alphafold.ebi.ac.uk/entry/${node2.id}" target="_blank" style="font-size:0.7em; color:#aaa;">(AF DB)</a><br>
                    <span style="font-size:0.85em">${node2.name}</span><br>
                    <span style="font-size:0.8em; color:#aaa">${node2.length} aa</span>`;
                
                document.getElementById('modal').style.display = 'block';
                
                document.getElementById('viewer1').innerHTML = '';
                document.getElementById('viewer2').innerHTML = '';
                
                if (syncId) cancelAnimationFrame(syncId);
                
                let v1 = $3Dmol.createViewer("viewer1", {backgroundColor: "0x1e1e1e"});
                let v2 = $3Dmol.createViewer("viewer2", {backgroundColor: "0x1e1e1e"});
                
                globalV1 = v1;
                globalV2 = v2;

                
                Promise.all([
                    fetch('/pdb/' + id1),
                    fetch('/pdb/' + id2)
                ]).then(async ([res1, res2]) => {
                    let pdb1 = await res1.text();
                    let pdb2 = await res2.text();
                    
                    if(!res1.ok) {
                        document.getElementById('viewer1').innerHTML = `<div style="padding:20px; color:#ff5555;"><b>Failed to fold ${id1}</b><br><br>${pdb1}</div>`;
                    } else {
                        v1.addModel(pdb1, "pdb");
                        v1.setStyle({}, {cartoon: {color: "cyan"}});
                        v1.zoomTo();
                        v1.render();
                    }
                    
                    if(!res2.ok) {
                        document.getElementById('viewer2').innerHTML = `<div style="padding:20px; color:#ff5555;"><b>Failed to fold ${id2}</b><br><br>${pdb2}</div>`;
                    } else {
                        v2.addModel(pdb2, "pdb");
                        v2.setStyle({}, {cartoon: {color: "magenta"}});
                        v2.zoomTo();
                        v2.render();
                    }
                    
                    document.getElementById('viewer1').onmouseenter = () => master = 1;
                    document.getElementById('viewer2').onmouseenter = () => master = 2;
                    
                    function sync() {
                        if (master === 1 && res1.ok && res2.ok) {
                            v2.setView(v1.getView());
                        } else if (master === 2 && res1.ok && res2.ok) {
                            v1.setView(v2.getView());
                        }
                        syncId = requestAnimationFrame(sync);
                    }
                    sync();
                });
            });
            // Cross-link over the Heatmap Axis Labels!
            chart.on('click', function(params) {
                if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                    window.open('https://www.uniprot.org/uniprotkb/' + params.value + '/entry', '_blank');
                }
            });
        });

        // Advanced UI Toggles
        var isMerged = false;
        var pocketSelectionMode = false;
        
        function toggleMergeView() {
            if(!globalV1 || !globalV2 || !window.currentPdb1 || !window.currentPdb2) return;
            isMerged = !isMerged;
            
            if (isMerged) {
                // Hide viewer 2, expand viewer 1
                document.querySelector('#viewer2').parentNode.style.display = 'none';
                document.querySelector('#viewer1').parentNode.style.width = '100%';
                
                globalV1.clear();
                // Load both into viewer 1 for Superposition
                globalV1.addModel(window.currentPdb1, "pdb");
                globalV1.addModel(window.currentPdb2, "pdb");
                
                globalV1.setStyle({model: 0}, {cartoon: {color: "cyan"}});
                globalV1.setStyle({model: 1}, {cartoon: {color: "magenta"}});
                globalV1.zoomTo();
                globalV1.render();
            } else {
                // Restore split view
                document.querySelector('#viewer2').parentNode.style.display = 'block';
                document.querySelector('#viewer1').parentNode.style.width = '50%';
                
                globalV1.clear();
                globalV1.addModel(window.currentPdb1, "pdb");
                globalV1.setStyle({}, {cartoon: {color: "cyan"}});
                globalV1.zoomTo();
                globalV1.render();
                
                globalV2.zoomTo();
                globalV2.render();
            }
        }
        
        function togglePocketSelection() {
            pocketSelectionMode = !pocketSelectionMode;
            let btn = document.getElementById('pocket-btn');
            
            if (pocketSelectionMode) {
                btn.style.background = '#4facfe';
                btn.style.color = '#000';
                btn.innerText = '🔴 Pocket Selection ON (Click Atoms)';
                
                globalV1.setClickable({}, true, function(atom, viewer, event, container) {
                    viewer.addSphere({ center: {x:atom.x, y:atom.y, z:atom.z}, radius: 1.5, color: 'yellow' });
                    viewer.addLabel(atom.resn + atom.resi, {position: {x:atom.x, y:atom.y, z:atom.z}, backgroundColor: 'black', fontColor:'white'});
                    viewer.render();
                    console.log("Selected Atom for Z-Order search:", atom.x, atom.y, atom.z);
                    alert(`Selected ${atom.resn} ${atom.resi} at (${atom.x.toFixed(2)}, ${atom.y.toFixed(2)}, ${atom.z.toFixed(2)})\n\n(Backend pg_bio Z-Order search integration coming next!)`);
                });
            } else {
                btn.style.background = '#333';
                btn.style.color = '#fff';
                btn.innerText = '🔍 Spatial Pocket Selection (pg_bio)';
                globalV1.removeAllShapes();
                globalV1.removeAllLabels();
                globalV1.setClickable({}, false);
                globalV1.render();
            }
        }
    </script>
</body>
</html>
"""

import os

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/pdb/<uniprot_id>')
def get_pdb(uniprot_id):
    path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{uniprot_id}-F1-model_v4.pdb"
    failed_path = f"{path}.failed"
    
    if os.path.exists(path):
        with open(path, 'r') as f:
            return f.read()
    if os.path.exists(failed_path):
        with open(failed_path, 'r') as f:
            return f.read(), 413
            
    import urllib.request
    try:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        return urllib.request.urlopen(url).read().decode('utf-8')
    except:
        # Fallback to local ESMFold if AlphaFold doesn't have it
        try:
            conn = psycopg.connect(DB_URI)
            cur = conn.cursor()
            cur.execute("SELECT sequence FROM proteins WHERE uniprot_id = %s", (uniprot_id,))
            res = cur.fetchone()
            if res and res[0]:
                seq = res[0]
                
                def attempt_fold(sequence):
                    cur.execute("SELECT bio_fold_sequence(%s)", (sequence,))
                    return cur.fetchone()[0]

                try:
                    pdb_data = attempt_fold(seq)
                    if pdb_data and "ATOM" in pdb_data:
                        os.makedirs(os.path.dirname(path), exist_ok=True)
                        with open(path, 'w') as f:
                            f.write(pdb_data)
                        return pdb_data
                    return "ESMFold returned empty structure.", 500
                except psycopg.errors.InternalError as e:
                    err_msg = str(e)
                    if "413" in err_msg:
                        # CHUNKING LOGIC FOR LARGE PROTEINS
                        chunk_size = 350
                        chunks = [seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size)]
                        merged_lines = []
                        
                        try:
                            for idx, chunk in enumerate(chunks):
                                conn.rollback() # reset transaction after 413 error
                                chunk_pdb = attempt_fold(chunk)
                                
                                for line in chunk_pdb.split("\n"):
                                    if line.startswith("ATOM  "):
                                        x = float(line[30:38])
                                        new_x = x + (60.0 * idx) # Shift 60 angstroms per chunk
                                        res_seq = int(line[22:26])
                                        new_res_seq = res_seq + (chunk_size * idx)
                                        atom_serial = int(line[6:11])
                                        new_atom_serial = (atom_serial + (4000 * idx)) % 100000
                                        
                                        new_line = line[:6] + f"{new_atom_serial:5d}" + line[11:22] + f"{new_res_seq:4d}" + line[26:30] + f"{new_x:8.3f}" + line[38:]
                                        merged_lines.append(new_line)
                            
                            merged_lines.append("TER")
                            final_pdb = "\n".join(merged_lines)
                            
                            os.makedirs(os.path.dirname(path), exist_ok=True)
                            with open(path, 'w') as f:
                                f.write(final_pdb)
                            return final_pdb
                            
                        except Exception as inner_e:
                            msg = f"Protein too large and chunking failed: {inner_e}"
                            with open(failed_path, 'w') as f: f.write(msg)
                            return msg, 413

                    msg = f"PostgreSQL ESMFold Error: {err_msg}"
                    with open(failed_path, 'w') as f: f.write(msg)
                    return msg, 500
            return "Sequence not found in local DB.", 404
        except Exception as e:
            return f"Error connecting to DB: {e}", 500

@app.route('/api/data')
def get_data():
    conn = psycopg.connect(DB_URI)
    cur = conn.cursor()
    
    # Get top 40 nodes from the graph, grouped by orphan_id
    cur.execute("""
        SELECT orphan_id, MAX(status), MAX(family_name), MAX(orphan_organism)
        FROM orphan_discoveries 
        GROUP BY orphan_id
        ORDER BY MAX(id) DESC 
        LIMIT 40
    """)
    rows = cur.fetchall()
    
    nodes_info = {r[0]: {'status': r[1], 'family': r[2], 'organism': r[3]} for r in rows}
    
    if not nodes_info:
        return jsonify({'nodes': [], 'matrix': []})
        
    # Get embeddings and extra protein info
    node_ids = list(nodes_info.keys())
    cur.execute("""
        SELECT uniprot_id, name, length(sequence), embedding 
        FROM proteins 
        WHERE uniprot_id = ANY(%s)
    """, (node_ids,))
    
    emb_dict = {}
    nodes = []
    
    for r in cur.fetchall():
        uid, name, length, emb_str = r
        emb_dict[uid] = np.array(eval(emb_str))
        
        info = nodes_info[uid]
        nodes.append({
            'id': uid,
            'status': info['status'],
            'family': info['family'],
            'organism': info['organism'],
            'name': name,
            'length': length
        })
    
    # Not all nodes might have embeddings instantly, so sort them predictably
    nodes.sort(key=lambda x: x['id'])
    
    matrix = []
    for i in range(len(nodes)):
        row = []
        emb1 = emb_dict[nodes[i]['id']]
        for j in range(len(nodes)):
            emb2 = emb_dict[nodes[j]['id']]
            # Cosine distance
            dist = scipy.spatial.distance.cosine(emb1, emb2)
            row.append(dist)
        matrix.append(row)
        
    return jsonify({'nodes': nodes, 'matrix': matrix})

import scipy.spatial.distance

if __name__ == '__main__':
    app.run(port=8080)
