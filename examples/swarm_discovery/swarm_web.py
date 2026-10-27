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
        body { background-color: #121212; color: #fff; font-family: sans-serif; margin: 0; overflow: hidden; display: flex; flex-direction: column; height: 100vh; }
        .topbar { display: flex; align-items: center; justify-content: space-between; padding: 10px 20px; background: #000; border-bottom: 1px solid #333; z-index: 10; }
        .topbar h2 { margin: 0; font-size: 1.2em; color: #fff; display: inline-block; }
        .topbar .subtitle { color: #ffaa00; font-size: 0.9em; margin-left: 15px; }
        .topbar .controls { display: flex; align-items: center; gap: 15px; }
        #heatmap { flex: 1; width: 100%; position: relative; }
        
        #modal { display: none; position: fixed; top: 6%; left: 6%; width: 88%; height: 88%; background: #1e1e1e; border: 2px solid #4facfe; box-shadow: 0 0 25px #000; z-index: 100; border-radius: 10px; overflow: hidden; }
        .modal-header { display: flex; justify-content: space-between; align-items: flex-start; padding: 10px 20px; background: #000; font-weight: bold; font-size: 1.1em; border-bottom: 1px solid #333; }
        .close-btn { cursor: pointer; color: #ff5555; white-space: nowrap; }
        #alignment-metrics-bar { display: flex; gap: 15px; flex-wrap: wrap; align-items: center; padding: 6px 20px; background: #141820; border-bottom: 1px solid #2b3548; font-size: 0.82em; color: #cbd5e1; }
        .metric-pill { background: #1e293b; border: 1px solid #334155; padding: 3px 8px; border-radius: 4px; }
        .metric-pill b { color: #4facfe; }
        .viewers-container { display: flex; width: 100%; height: calc(100% - 95px); }
        .viewer-box { width: 50%; height: 100%; position: relative; border-right: 1px solid #2a2a2a; }
        .viewer-title { position: absolute; top: 10px; left: 10px; z-index: 10; background: rgba(0,0,0,0.85); border: 1px solid #333; padding: 10px 14px; border-radius: 6px; max-width: 85%; font-size: 0.88em; line-height: 1.4; }
        .viewer-title a { text-decoration: none; }
        .viewer-title a:hover { text-decoration: underline; }
        .domain-badge { display: inline-block; background: #263238; border: 1px solid #4facfe; color: #80d8ff; font-size: 0.78em; padding: 1px 6px; border-radius: 3px; margin: 2px 3px 0 0; }
        .controls-bar { margin-top: 6px; font-size: 0.8em; display: flex; flex-wrap: wrap; gap: 5px; align-items: center; }
        .btn { background: #333; color: #fff; border: 1px solid #555; padding: 4px 9px; cursor: pointer; border-radius: 4px; }
        .btn:hover { background: #555; }
        
        /* Custom tooltip for axis labels */
        #axis-tooltip { display: none; position: absolute; background: rgba(0,0,0,0.9); border: 1px solid #4facfe; padding: 10px; border-radius: 5px; z-index: 50; pointer-events: none; font-size: 12px; text-align: left; max-width: 300px; }
    </style>
</head>
<body>
    <div class="topbar">
        <div>
            <h2>Protein Correlation Swarm</h2>
            <span class="subtitle">Pending (gray) vs Explored (blue). Click cell to view 3D-aligned models &amp; full metadata.</span>
        </div>
        <div class="controls">
            <label style="color: #fff; font-size: 0.9em;">Grid View Limit: </label>
            <select id="grid-limit" onchange="loadData()" style="background: #333; color: #fff; border: 1px solid #555; padding: 4px; border-radius: 4px;">
                <option value="40" selected>Latest 40</option>
                <option value="best-40">Top Best 40 Matches</option>
                <option value="100">Latest 100</option>
                <option value="250">Latest 250</option>
                <option value="all">Latest 800 (Max)</option>
            </select>
        </div>
    </div>
    <div id="heatmap"></div>
    <div id="axis-tooltip"></div>

    <div id="modal">
        <div class="modal-header">
            <div>
                <span id="modal-title">Comparing Proteins</span>
                <div class="controls-bar">
                    <button class="btn" onclick="applyDefaultColoring()">🎨 Color by Identity</button>
                    <button class="btn" onclick="applyConfidenceColoring()">🎯 Color by pLDDT Confidence</button>
                    <button class="btn" onclick="applySecondaryStructureColoring()">🌀 Color by Secondary Structure</button>
                    <button class="btn" onclick="applyAcidicSurfaceColoring()">🧪 Highlight Acidic Shell (D/E)</button>
                    <span style="border-left: 1px solid #555; margin: 0 5px;"></span>
                    <button class="btn" onclick="toggleMergeView()">🧬 Merge View (3D Superposition)</button>
                    <button id="pocket-btn" class="btn" onclick="togglePocketSelection()">🔍 Spatial Pocket Selection (pg_bio)</button>
                    <span style="border-left: 1px solid #555; margin: 0 5px;"></span>
                    <button class="btn" onclick="copyShareLink()">🔗 Share Link</button>
                </div>
            </div>
            <span class="close-btn" onclick="closeModal()">✖ Close (Esc)</span>
        </div>
        <div id="alignment-metrics-bar">Loading 3D structural superposition &amp; pg_bio biochemical metadata...</div>
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
            let newUrl = new URL(window.location);
            newUrl.searchParams.delete('node1');
            newUrl.searchParams.delete('node2');
            window.history.pushState({}, '', newUrl);
            if (window.chart) {
                window.chart.dispatchAction({ type: 'downplay' });
            }
        }
        
        function copyShareLink() {
            navigator.clipboard.writeText(window.location.href).then(() => {
                alert('Correlation link copied to clipboard! You can share this URL.');
            }).catch(err => {
                alert('Failed to copy: ' + window.location.href);
            });
        }
        
        document.addEventListener('keydown', function(e) {
            if (e.key === 'Escape' && document.getElementById('modal').style.display === 'block') {
                closeModal();
            }
        });
        
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

        function applySecondaryStructureColoring() {
            if(!globalV1 || !globalV2) return;
            [globalV1, globalV2].forEach(v => {
                v.getModels().forEach(m => m.calcSS());
                v.setStyle({}, {cartoon: {color: 'white', opacity: 0.35}});
                v.setStyle({ss: 'h'}, {cartoon: {color: '#ff6ec7'}});
                v.setStyle({ss: 's'}, {cartoon: {color: '#00f2fe'}});
                v.render();
            });
        }

        function applyAcidicSurfaceColoring() {
            if(!globalV1 || !globalV2) return;
            [globalV1, globalV2].forEach(v => {
                v.setStyle({}, {cartoon: {color: 'white', opacity: 0.35}});
                v.setStyle({resn: ['ASP', 'GLU']}, {cartoon: {color: '#ff4444'}, stick: {colorscheme: 'redCarbon'}});
                v.render();
            });
        }

        function loadData() {
            var limit = document.getElementById('grid-limit') ? document.getElementById('grid-limit').value : '40';
            let urlParams = new URLSearchParams(window.location.search);
            let p1 = urlParams.get('node1') === 'undefined' ? '' : (urlParams.get('node1') || '');
            let p2 = urlParams.get('node2') === 'undefined' ? '' : (urlParams.get('node2') || '');
            chart.showLoading({text: 'Fetching Matrix...', color: '#4facfe', textColor: '#fff', maskColor: 'rgba(0,0,0,0.8)'});
            fetch(`/api/data?limit=${limit}&n1=${p1}&n2=${p2}`).then(r => r.json()).then(data => {
                chart.hideLoading();
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
                        backgroundColor: 'rgba(0,0,0,0.9)',
                        borderColor: '#4facfe',
                        borderWidth: 1,
                        textStyle: { color: '#fff' },
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
                if (params.componentType !== 'series' || !params.value) return;
                
                let id1 = labels[params.value[0]];
                let id2 = labels[params.value[1]];
                if (!id1 || !id2) return;
                
                let node1 = data.nodes.find(n => n.id === id1);
                let node2 = data.nodes.find(n => n.id === id2);
                if (!node1 || !node2) return;
                
                let distance = params.value[2];
                if (distance === '-') {
                    distance = data.matrix[params.value[0]][params.value[1]];
                }
                
                if (distance === '-' || distance === undefined || distance === null) return;
                
                openModal(id1, id2, distance, node1, node2);
            });
            
            // Auto-open modal from URL params
            let initialParams = new URLSearchParams(window.location.search);
            let init_p1 = initialParams.get('node1');
            let init_p2 = initialParams.get('node2');
            if (init_p1 && init_p2 && init_p1 !== 'undefined' && init_p2 !== 'undefined') {
                let n1 = data.nodes.find(n => n.id === init_p1);
                let n2 = data.nodes.find(n => n.id === init_p2);
                if (n1 && n2) {
                    let idx1 = labels.indexOf(init_p1);
                    let idx2 = labels.indexOf(init_p2);
                    if (idx1 >= 0 && idx2 >= 0) {
                        let distance = data.matrix[idx2][idx1];
                        
                        setTimeout(() => {
                            chart.dispatchAction({
                                type: 'highlight',
                                seriesIndex: 0,
                                dataIndex: idx2 * labels.length + idx1
                            });
                        }, 500);
                        
                        openModal(init_p1, init_p2, distance, n1, n2);
                    }
                }
            }
            
            function formatMetadataCard(node, meta, color) {
                let domainsHtml = (meta.domains && meta.domains.length > 0)
                    ? '<div style="margin-top:4px;">' + meta.domains.map(d => `<span class="domain-badge">${d}</span>`).join('') + '</div>'
                    : '<div style="margin-top:4px;"><span class="domain-badge" style="border-color:#555;color:#aaa;">Dark Proteome (No Pfam/InterPro)</span></div>';
                return `<b style="font-size:1.05em; color:#fff;">${meta.organism || node.organism}</b>
                    <span style="font-size:0.75em; background:#222; border:1px solid #444; padding:1px 5px; border-radius:3px; margin-left:6px; color:#ffaa00;">${node.family || 'Orphan'}</span><br>
                    <a href="https://www.uniprot.org/uniprotkb/${node.id}/entry" target="_blank" style="color:${color}; font-weight:bold;">${node.id} 🔗</a> 
                    <a href="https://alphafold.ebi.ac.uk/entry/${node.id}" target="_blank" style="font-size:0.75em; color:#aaa;">(AlphaFold DB)</a>
                    <span style="font-size:0.75em; color:#94a3b8; margin-left:6px;">pLDDT: <b style="color:#38bdf8">${meta.mean_plddt}</b></span><br>
                    <span style="font-size:0.8em; color:#cbd5e1;">${node.name}</span><br>
                    <span style="font-size:0.78em; color:#94a3b8;">
                        <b>${meta.length} aa</b> | <b>${meta.mw_kda} kDa</b> | Entropy: <b>${meta.entropy}</b> |
                        Acidic(D/E): <b style="color:#f87171">${meta.acidic_pct}%</b> |
                        Basic: <b style="color:#60a5fa">${meta.basic_pct}%</b> |
                        Hydrophobic: <b>${meta.hydro_pct}%</b>
                    </span>
                    ${domainsHtml}`;
            }

            function openModal(id1, id2, distance, node1, node2) {
                // Update URL
                let newUrl = new URL(window.location);
                newUrl.searchParams.set('node1', id1);
                newUrl.searchParams.set('node2', id2);
                window.history.pushState({}, '', newUrl);

                document.getElementById('modal-title').innerHTML = `Comparing: <span style="color:cyan">${id1}</span> vs <span style="color:magenta">${id2}</span> <span style="color:#aaa; font-size:0.8em; margin-left:10px;">Vector Dist: ${distance.toFixed(4)}</span>`;
                document.getElementById('alignment-metrics-bar').innerHTML = `⏳ Running 3D Kabsch Superposition &amp; querying <code>pg_bio</code> + UniProt metadata...`;

                document.getElementById('v1-title').innerHTML = `<b style="font-size:1.1em;">${node1.organism}</b><br><span style="color:cyan">${node1.id}</span> (${node1.length} aa)`;
                document.getElementById('v2-title').innerHTML = `<b style="font-size:1.1em;">${node2.organism}</b><br><span style="color:magenta">${node2.id}</span> (${node2.length} aa)`;
                
                document.getElementById('modal').style.display = 'block';
                
                document.getElementById('viewer1').innerHTML = '';
                document.getElementById('viewer2').innerHTML = '';
                
                if (syncId) cancelAnimationFrame(syncId);
                
                let v1 = $3Dmol.createViewer("viewer1", {backgroundColor: "0x1e1e1e"});
                let v2 = $3Dmol.createViewer("viewer2", {backgroundColor: "0x1e1e1e"});
                
                globalV1 = v1;
                globalV2 = v2;

                fetch(`/api/compare/${id1}/${id2}`).then(r => r.json()).then(comp => {
                    if (comp.error) {
                        document.getElementById('alignment-metrics-bar').innerHTML = `<span style="color:#ff5555">Error: ${comp.error}</span>`;
                        return;
                    }
                    let al = comp.alignment;
                    document.getElementById('alignment-metrics-bar').innerHTML = `
                        <span class="metric-pill">3D Core RMSD: <b>${al.core_rmsd} Å</b> (${al.aligned_residues} CA atoms)</span>
                        <span class="metric-pill">TM-Score: <b>${al.tm_score}</b></span>
                        <span class="metric-pill">Smith-Waterman (<code>pg_bio</code>): <b>${al.sw_score}</b></span>
                        <span class="metric-pill">Tri-peptide Cosine Dist: <b>${al.esm_kmer_dist}</b></span>
                        <span class="metric-pill">3D Orientation: <b style="color:#4ade80">✔ Kabsch Superimposed</b></span>
                    `;
                    document.getElementById('v1-title').innerHTML = formatMetadataCard(node1, comp.meta1, 'cyan');
                    document.getElementById('v2-title').innerHTML = formatMetadataCard(node2, comp.meta2, 'magenta');

                    let pdb1 = comp.pdb1;
                    let pdb2 = comp.pdb2_aligned;

                    v1.addModel(pdb1, "pdb");
                    v1.setStyle({}, {cartoon: {color: "cyan"}});
                    v1.zoomTo();
                    v1.render();
                    window.currentPdb1 = pdb1;

                    v2.addModel(pdb2, "pdb");
                    v2.setStyle({}, {cartoon: {color: "magenta"}});
                    v2.setView(v1.getView());
                    v2.render();
                    window.currentPdb2 = pdb2;

                    document.getElementById('viewer1').onmouseenter = () => master = 1;
                    document.getElementById('viewer2').onmouseenter = () => master = 2;
                    
                    function sync() {
                        if (master === 1) {
                            v2.setView(v1.getView());
                        } else if (master === 2) {
                            v1.setView(v2.getView());
                        }
                        syncId = requestAnimationFrame(sync);
                    }
                    sync();
                });
            }
            // Cross-link over the Heatmap Axis Labels!
            chart.on('click', function(params) {
                if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                    window.open('https://www.uniprot.org/uniprotkb/' + params.value + '/entry', '_blank');
                }
            });
        }); // End of fetch

        } // End of loadData

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
        
        // Initial load
        loadData();
    </script>
</body>
</html>
"""

import os
import json
import urllib.request
import numpy as np

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

def fetch_pdb_string(uniprot_id):
    path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{uniprot_id}-F1-model_v4.pdb"
    if os.path.exists(path):
        with open(path, 'r') as f:
            content = f.read()
            if "ATOM" in content:
                return content, 200

    # Try AlphaFold v6 then v4
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
                    return pdb_data, 200
                return "ESMFold returned empty structure.", 500
            except psycopg.errors.InternalError as e:
                err_msg = str(e)
                if "413" in err_msg:
                    chunk_size = 350
                    chunks = [seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size)]
                    merged_lines = []
                    try:
                        for idx, chunk in enumerate(chunks):
                            conn.rollback()
                            chunk_pdb = attempt_fold(chunk)
                            for line in chunk_pdb.split("\n"):
                                if line.startswith("ATOM  "):
                                    x = float(line[30:38])
                                    new_x = x + (60.0 * idx)
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
                        return final_pdb, 200
                    except Exception as inner_e:
                        return f"Protein too large and chunking failed: {inner_e}", 413
                return f"PostgreSQL ESMFold Error: {err_msg}", 500
        return "Sequence not found in local DB.", 404
    except Exception as e:
        return f"Error connecting to DB: {e}", 500

@app.route('/pdb/<uniprot_id>')
def get_pdb(uniprot_id):
    pdb_str, status = fetch_pdb_string(uniprot_id)
    return pdb_str, status

def kabsch_superpose_pdb(pdb_ref, pdb_mov):
    def get_ca(pdb_str):
        coords, bfacs = [], []
        for l in pdb_str.splitlines():
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                b = float(l[60:66])
                bfacs.append(b * 100.0 if b <= 1.0 else b)
        return np.array(coords), (np.mean(bfacs) if bfacs else 0.0)

    c1, plddt1 = get_ca(pdb_ref)
    c2, plddt2 = get_ca(pdb_mov)
    if len(c1) < 5 or len(c2) < 5:
        return pdb_mov, 0.0, 0.0, 0, plddt1, plddt2

    def fit_kabsch(P, Q):
        cP, cQ = P.mean(axis=0), Q.mean(axis=0)
        P0, Q0 = P - cP, Q - cQ
        H = Q0.T @ P0
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(U @ Vt))
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        return R, cP, cQ

    n1, n2 = len(c1), len(c2)
    L = min(n1, n2)
    d0 = max(1.24 * ((L - 15) ** (1.0 / 3.0)) - 1.8, 0.5) if L > 15 else 0.5
    best_tm = -1.0
    best_transform = (np.eye(3), np.zeros(3), np.zeros(3))
    best_stats = (9.99, 0.0, L)

    import scipy.spatial.distance
    w = min(60, L)
    seeds = [(0, 0, L)]
    for i in range(0, n1 - w + 1, max(15, w // 3)):
        for j in range(max(0, i - 30), min(n2 - w + 1, i + 31), 5):
            seeds.append((i, j, w))

    for i1, i2, win in seeds:
        R, cP, cQ = fit_kabsch(c1[i1:i1+win], c2[i2:i2+win])
        for cutoff in (7.0, 4.5, 3.0):
            c2_al = (c2 - cQ) @ R + cP
            D = scipy.spatial.distance.cdist(c1, c2_al)
            idx2 = np.argmin(D, axis=1)
            min_d = D[np.arange(n1), idx2]
            mask = min_d < cutoff
            if mask.sum() >= max(15, int(0.15 * L)):
                R, cP, cQ = fit_kabsch(c1[mask], c2[idx2[mask]])
        c2_al = (c2 - cQ) @ R + cP
        D = scipy.spatial.distance.cdist(c1, c2_al)
        min_d = np.min(D, axis=1)
        tm = float(np.sum(1.0 / (1.0 + (min_d / d0) ** 2)) / L)
        core_mask = min_d < 3.5
        if core_mask.sum() >= 10:
            n_core = int(core_mask.sum())
            core_rmsd = float(np.sqrt(np.mean(min_d[core_mask] ** 2)))
        else:
            top_k = np.sort(min_d)[:max(10, int(0.4 * L))]
            n_core = len(top_k)
            core_rmsd = float(np.sqrt(np.mean(top_k ** 2)))
        if tm > best_tm:
            best_tm = tm
            best_transform = (R, cP, cQ)
            best_stats = (core_rmsd, tm, n_core)

    R, cP, cQ = best_transform
    out_lines = []
    for line in pdb_mov.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            new_xyz = (xyz - cQ) @ R + cP
            line = f"{line[:30]}{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}{line[54:]}"
        out_lines.append(line)
    return "\n".join(out_lines), best_stats[0], best_stats[1], best_stats[2], plddt1, plddt2

@app.route('/api/compare/<id1>/<id2>')
def compare_proteins(id1, id2):
    pdb1, s1 = fetch_pdb_string(id1)
    pdb2, s2 = fetch_pdb_string(id2)
    if s1 != 200 or s2 != 200:
        return jsonify({'error': f'Could not load 3D structure ({id1}: {s1}, {id2}: {s2})'}), 500

    pdb2_aligned, core_rmsd, tm_score, n_core, plddt1, plddt2 = kabsch_superpose_pdb(pdb1, pdb2)

    conn = psycopg.connect(DB_URI)
    cur = conn.cursor()

    def build_meta(uid, plddt):
        cur.execute("SELECT name, sequence, molecular_weight(sequence), sequence_entropy(sequence) FROM proteins WHERE uniprot_id = %s", (uid,))
        row = cur.fetchone()
        name, seq, mw, ent = row if row else ("Unknown", "", 0.0, 0.0)
        L = max(len(seq), 1)
        acidic = (seq.count('D') + seq.count('E')) / L * 100.0
        basic = (seq.count('K') + seq.count('R') + seq.count('H')) / L * 100.0
        hydro = sum(seq.count(c) for c in 'AILMFWV') / L * 100.0
        domains = []
        organism = ""
        try:
            u = json.loads(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=4).read().decode())
            organism = u.get("organism", {}).get("scientificName", "")
            for x in u.get("uniProtKBCrossReferences", []):
                if x.get("database") in ("Pfam", "InterPro"):
                    val = x.get("properties", [{}])[0].get("value", "")
                    domains.append(f"{x['id']}: {val}" if val else x['id'])
            for feat in u.get("features", []):
                if feat.get("type") == "Transmembrane":
                    s = feat.get("location", {}).get("start", {}).get("value", "?")
                    e = feat.get("location", {}).get("end", {}).get("value", "?")
                    domains.append(f"TM Helix ({s}-{e})")
        except Exception:
            pass
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
            'domains': domains[:6],
            'seq': seq
        }

    m1 = build_meta(id1, plddt1)
    m2 = build_meta(id2, plddt2)

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
    from flask import request
    import scipy.spatial.distance
    conn = psycopg.connect(DB_URI)
    cur = conn.cursor()
    
    limit_param = request.args.get('limit', 'best-40')
    if limit_param.lower() == 'best-40':
        cur.execute("""
            WITH ranked_edges AS (
                SELECT
                    id, bait_id, orphan_id, vector_distance, family_name,
                    bait_organism, orphan_organism, status,
                    ROW_NUMBER() OVER (PARTITION BY family_name ORDER BY vector_distance ASC, id DESC) AS fam_rank
                FROM orphan_discoveries
                WHERE bait_id != orphan_id
            )
            SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism, vector_distance
            FROM ranked_edges
            ORDER BY fam_rank ASC, vector_distance ASC, id DESC
            LIMIT 80
        """)
        edge_rows = cur.fetchall()
        nodes_info = {}
        for b_id, o_id, st, fam, b_org, o_org, v_dist in edge_rows:
            if b_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[b_id] = {'status': 'Explored', 'family': fam, 'organism': b_org}
            if o_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[o_id] = {'status': st, 'family': fam, 'organism': o_org}
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
            SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism
            FROM orphan_discoveries
            WHERE bait_id != orphan_id
            ORDER BY id DESC
            LIMIT %s
        """, (max_nodes * 2,))
        edge_rows = cur.fetchall()
        nodes_info = {}
        for b_id, o_id, st, fam, b_org, o_org in edge_rows:
            if b_id not in nodes_info and len(nodes_info) < max_nodes:
                nodes_info[b_id] = {'status': 'Explored', 'family': fam, 'organism': b_org}
            if o_id not in nodes_info and len(nodes_info) < max_nodes:
                nodes_info[o_id] = {'status': st, 'family': fam, 'organism': o_org}
            if len(nodes_info) >= max_nodes:
                break
    
    n1 = request.args.get('n1')
    n2 = request.args.get('n2')
    extra_nodes = []
    if n1 and n1 != 'undefined' and n1 not in nodes_info: extra_nodes.append(n1)
    if n2 and n2 != 'undefined' and n2 not in nodes_info: extra_nodes.append(n2)
    
    if extra_nodes:
        cur.execute("""
            SELECT orphan_id, MAX(status), MAX(family_name), MAX(orphan_organism)
            FROM orphan_discoveries 
            WHERE orphan_id = ANY(%s)
            GROUP BY orphan_id
        """, (extra_nodes,))
        for r in cur.fetchall():
            nodes_info[r[0]] = {'status': r[1], 'family': r[2], 'organism': r[3]}
        missing_still = [x for x in extra_nodes if x not in nodes_info]
        if missing_still:
            cur.execute("""
                SELECT bait_id, 'Explored', MAX(family_name), MAX(bait_organism)
                FROM orphan_discoveries
                WHERE bait_id = ANY(%s)
                GROUP BY bait_id
            """, (missing_still,))
            for r in cur.fetchall():
                nodes_info[r[0]] = {'status': r[1], 'family': r[2], 'organism': r[3]}
    
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
    
    # Sort them predictably
    nodes.sort(key=lambda x: x['id'])
    
    if not nodes:
        return jsonify({'nodes': [], 'matrix': []})

    # Vectorized cosine distance matrix
    emb_list = [emb_dict[n['id']] for n in nodes]
    emb_matrix = np.array(emb_list)
    dist_matrix = scipy.spatial.distance.cdist(emb_matrix, emb_matrix, metric='cosine')
    
    # Fill diagonal with 0.0 to fix floating point epsilon (like 1.1e-16)
    np.fill_diagonal(dist_matrix, 0.0)

    matrix = dist_matrix.tolist()
        
    return jsonify({'nodes': nodes, 'matrix': matrix})

if __name__ == '__main__':
    app.run(port=8080)
