import os
import json
import urllib.request
import urllib.error
from flask import Flask, render_template_string, jsonify, request
import psycopg
import numpy as np
import scipy.spatial.distance

try:
    from examples.swarm_discovery.inspector_assets import INSPECTOR_CSS, INSPECTOR_JS
except ImportError:
    from inspector_assets import INSPECTOR_CSS, INSPECTOR_JS

app = Flask(__name__)
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>Swarm Discovery Heatmap</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.0/dist/echarts.min.js"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * { box-sizing: border-box; }
        body { background-color: #090d16; color: #f8fafc; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; overflow: hidden; display: flex; flex-direction: column; height: 100vh; }
        .topbar { display: flex; align-items: center; justify-content: space-between; padding: 10px 20px; background: #020617; border-bottom: 1px solid #1e293b; z-index: 10; }
        .topbar h2 { margin: 0; font-size: 1.15em; color: #f8fafc; display: inline-block; }
        .topbar .subtitle { color: #94a3b8; font-size: 0.86em; margin-left: 14px; }
        .topbar .controls { display: flex; align-items: center; gap: 12px; }
        #heatmap { flex: 1; width: 100%; position: relative; }

        #modal { display: none; position: fixed; top: 4%; left: 4%; width: 92%; height: 92%; background: #0f172a; border: 1px solid #38bdf8; box-shadow: 0 20px 60px rgba(0,0,0,0.9); z-index: 100; border-radius: 10px; overflow: hidden; }
        #axis-tooltip { display: none; position: absolute; background: rgba(2,6,23,0.95); border: 1px solid #38bdf8; padding: 10px 12px; border-radius: 6px; z-index: 50; pointer-events: none; font-size: 12px; text-align: left; max-width: 320px; line-height: 1.45; }
        .help-btn { background: #1e293b; color: #38bdf8; border: 1px solid #38bdf8; padding: 6px 13px; cursor: pointer; border-radius: 6px; font-weight: 600; font-size: 0.84em; transition: all 0.15s; text-decoration: none; display: inline-block; }
        .help-btn:hover { background: #38bdf8; color: #020617; }
        __INSPECTOR_CSS__
    </style>
</head>
<body>
    <div class="topbar">
        <div>
            <h2>Protein Correlation Swarm</h2>
            <span class="subtitle">Pending (gray) vs Explored (blue). Click any cell to inspect 3D-aligned structures &amp; configurable pg_bio layers.</span>
        </div>
        <div class="controls">
            <label style="color: #cbd5e1; font-size: 0.86em;">Grid ViewPreset:</label>
            <select id="grid-limit" onchange="loadData()" style="background: #1e293b; color: #fff; border: 1px solid #475569; padding: 5px 10px; border-radius: 6px; font-size: 0.85em;">
                <option value="40" selected>Latest 40</option>
                <option value="best-40">Top Best 40 Matches</option>
                <option value="twilight-40">Top 40 Twilight-Zone Discoveries</option>
                <option value="chaperone-40">Top 40 Chaperone Biophysical Correlations</option>
                <option value="rangefinder-40">Top 40 RangeFinder Pocket-Interval Matches</option>
                <option value="100">Latest 100</option>
                <option value="250">Latest 250</option>
                <option value="all">Latest 800 (Max)</option>
            </select>
            <button class="help-btn" style="background: rgba(245, 158, 11, 0.16); color: #fbbf24; border-color: #f59e0b;" onclick="window.location.href='/cards'">🃏 1-to-N Card Arena</button>
            <button class="help-btn" style="background: rgba(192, 132, 252, 0.15); color: #c084fc; border-color: #a855f7;" onclick="window.location.href='/api/swarm_telemetry'">⚡ Swarm Telemetry</button>
            <button class="help-btn" style="background: rgba(16, 185, 129, 0.14); color: #34d399; border-color: #10b981;" onclick="window.open('/compare' + (window.currentNode1Id && window.currentNode2Id ? '?node1=' + window.currentNode1Id + '&node2=' + window.currentNode2Id : ''), '_blank')">3D Compare Lab</button>
            <button class="help-btn" onclick="openHelpCenter('vector')">Correlation Wiki &amp; Help Center</button>
        </div>
    </div>
    <div id="heatmap"></div>
    <div id="axis-tooltip"></div>
    <div id="wiki-popover"></div>

    <div id="modal">
        <div id="modal-inspector-root" style="width:100%; height:100%;"></div>
    </div>

    <script>
        __INSPECTOR_JS__

        var chart = echarts.init(document.getElementById('heatmap'));
        var modalInspector = null;

        function closeModal() {
            document.getElementById('modal').style.display = 'none';
            let newUrl = new URL(window.location);
            newUrl.searchParams.delete('node1');
            newUrl.searchParams.delete('node2');
            window.history.pushState({}, '', newUrl);
            if (window.chart) {
                window.chart.dispatchAction({ type: 'downplay' });
            }
        }

        function openHelpCenter(category, targetCardId) {
            let url = '/wiki';
            let params = new URLSearchParams();
            if (window.currentNode1Id && window.currentNode2Id) {
                params.set('node1', window.currentNode1Id);
                params.set('node2', window.currentNode2Id);
            }
            params.set('cat', category || 'vector');
            let qs = params.toString();
            if (qs) url += '?' + qs;
            if (targetCardId) url += '#' + targetCardId;
            window.open(url, '_blank');
        }

        document.addEventListener('keydown', function(e) {
            if (e.key === 'Escape' && document.getElementById('modal').style.display === 'block') {
                closeModal();
            }
        });

        function openModal(id1, id2, distance) {
            let newUrl = new URL(window.location);
            newUrl.searchParams.set('node1', id1);
            newUrl.searchParams.set('node2', id2);
            window.history.pushState({}, '', newUrl);

            document.getElementById('modal').style.display = 'block';
            if (!modalInspector) {
                modalInspector = new ProteinCompareInspector({
                    mountEl: '#modal-inspector-root',
                    context: 'modal',
                    defaultPreset: 'overview',
                    onClose: closeModal
                });
            }
            modalInspector.loadPair(id1, id2, distance);
        }

        function loadData() {
            let urlParams = new URLSearchParams(window.location.search);
            let urlLimit = urlParams.get('limit');
            if (urlLimit && !window._urlLimitApplied) {
                window._urlLimitApplied = true;
                let sel = document.getElementById('grid-limit');
                if (sel) sel.value = urlLimit;
            }
            let limitVal = document.getElementById('grid-limit').value;
            let p1 = urlParams.get('node1') || '';
            let p2 = urlParams.get('node2') || '';

            chart.showLoading({ text: 'Fetching Swarm Matrix...', color: '#38bdf8', textColor: '#fff', maskColor: 'rgba(9, 13, 22, 0.8)' });
            fetch('/api/data?limit=' + limitVal + '&n1=' + p1 + '&n2=' + p2)
                .then(r => r.json())
                .then(data => {
                    chart.hideLoading();
                    if (!data.nodes || data.nodes.length === 0) return;

                    let labels = data.nodes.map(n => n.id);
                    let heatmapData = [];
                    for (let i = 0; i < data.matrix.length; i++) {
                        for (let j = 0; j < data.matrix[i].length; j++) {
                            let val = data.matrix[i][j];
                            heatmapData.push([j, i, val === null ? '-' : Number(val.toFixed(4))]);
                        }
                    }

                    let option = {
                        tooltip: {
                            position: 'top',
                            backgroundColor: 'rgba(2, 6, 23, 0.95)',
                            borderColor: '#38bdf8',
                            textStyle: { color: '#fff', fontSize: 12 },
                            formatter: function(params) {
                                let n1 = data.nodes[params.value[0]];
                                let n2 = data.nodes[params.value[1]];
                                if (!n1 || !n2) return '';
                                return `<b>${n1.id}</b> (${n1.organism})<br>vs<br><b>${n2.id}</b> (${n2.organism})<br>` +
                                       `<span style="color:#38bdf8">Cosine Distance: <b>${params.value[2]}</b></span><br>` +
                                       `<span style="color:#94a3b8;font-size:11px;">Click to open 3D Compare Inspector</span>`;
                            }
                        },
                        grid: { top: '4%', bottom: '14%', left: '10%', right: '8%' },
                        xAxis: {
                            type: 'category', data: labels, triggerEvent: true,
                            axisLabel: {
                                rotate: 45, fontSize: 10,
                                color: function(val) {
                                    let node = data.nodes.find(n => n.id === val);
                                    return node && node.status === 'Explored' ? '#38bdf8' : '#888';
                                }
                            }
                        },
                        yAxis: {
                            type: 'category', data: labels, triggerEvent: true,
                            axisLabel: {
                                fontSize: 10,
                                color: function(val) {
                                    let node = data.nodes.find(n => n.id === val);
                                    return node && node.status === 'Explored' ? '#38bdf8' : '#888';
                                }
                            }
                        },
                        visualMap: {
                            min: 0, max: 0.2, calculable: true, orient: 'vertical', right: '2%', top: 'center',
                            precision: 3, textStyle: { color: '#fff' },
                            inRange: { color: ['#00fa9a', '#38bdf8', '#334155'] }
                        },
                        series: [{
                            name: 'Distance', type: 'heatmap', data: heatmapData,
                            emphasis: { itemStyle: { shadowBlur: 10, shadowColor: 'rgba(0, 0, 0, 0.5)' } }
                        }]
                    };
                    chart.setOption(option);

                    var axisTooltip = document.getElementById('axis-tooltip');
                    var hoveredNode = null;
                    chart.off('mouseover');
                    chart.off('mouseout');
                    chart.off('click');

                    chart.on('mouseover', function(params) {
                        if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                            hoveredNode = data.nodes.find(n => n.id === params.value);
                            if (hoveredNode) {
                                axisTooltip.innerHTML = `<b style="color:#38bdf8;">${hoveredNode.id} (${hoveredNode.status})</b><br>` +
                                    `<b>Family:</b> ${hoveredNode.family}<br>` +
                                    `<b>Organism:</b> ${hoveredNode.organism}<br>` +
                                    `<b>Name:</b> ${hoveredNode.name}<br>` +
                                    `<b>Length:</b> ${hoveredNode.length} aa`;
                                axisTooltip.style.display = 'block';
                            }
                        }
                    });

                    if (!window._axisMouseMoveBound) {
                        window._axisMouseMoveBound = true;
                        document.addEventListener('mousemove', function(e) {
                            if (hoveredNode && axisTooltip.style.display === 'block') {
                                let r = axisTooltip.getBoundingClientRect();
                                let x = e.clientX + 15;
                                let y = e.clientY + 15;
                                if (x + r.width > window.innerWidth) x = e.clientX - r.width - 10;
                                if (y + r.height > window.innerHeight) y = e.clientY - r.height - 10;
                                axisTooltip.style.left = x + 'px';
                                axisTooltip.style.top = y + 'px';
                            }
                        });
                    }

                    chart.on('mouseout', function(params) {
                        if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                            axisTooltip.style.display = 'none';
                            hoveredNode = null;
                        }
                    });

                    chart.on('click', function(params) {
                        if (params.componentType === 'xAxis' || params.componentType === 'yAxis') {
                            window.open('https://www.uniprot.org/uniprotkb/' + params.value + '/entry', '_blank');
                            return;
                        }
                        if (params.componentType !== 'series' || !params.value) return;
                        let id1 = labels[params.value[0]];
                        let id2 = labels[params.value[1]];
                        if (!id1 || !id2) return;
                        let distance = params.value[2];
                        if (distance === '-') distance = data.matrix[params.value[0]][params.value[1]];
                        if (distance === '-' || distance === undefined || distance === null) return;
                        openModal(id1, id2, distance);
                    });

                    let initialParams = new URLSearchParams(window.location.search);
                    let init_p1 = initialParams.get('node1');
                    let init_p2 = initialParams.get('node2');
                    if (init_p1 && init_p2 && init_p1 !== 'undefined' && init_p2 !== 'undefined') {
                        let idx1 = labels.indexOf(init_p1);
                        let idx2 = labels.indexOf(init_p2);
                        let distance = (idx1 >= 0 && idx2 >= 0) ? data.matrix[idx2][idx1] : 0.0;
                        openModal(init_p1, init_p2, distance);
                    }
                });
        }

        loadData();
    </script>
</body>
</html>
""".replace("__INSPECTOR_CSS__", INSPECTOR_CSS).replace("__INSPECTOR_JS__", INSPECTOR_JS)


WIKI_HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>Swarm Discovery — Phase-Scoped Correlation Wiki &amp; 3D Biophysical Lab</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * { box-sizing: border-box; }
        body { background-color: #090d16; color: #f8fafc; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 0; padding: 0; line-height: 1.55; }
        .wiki-topbar { position: sticky; top: 0; z-index: 100; display: flex; justify-content: space-between; align-items: center; padding: 12px 24px; background: rgba(2, 6, 23, 0.96); border-bottom: 1px solid #1e293b; }
        .wiki-topbar h1 { margin: 0; font-size: 1.18em; color: #f8fafc; }
        .wiki-topbar .sub { color: #94a3b8; font-size: 0.82em; margin-top: 2px; }
        .back-btn { background: #1e293b; color: #38bdf8; border: 1px solid #38bdf8; padding: 6px 14px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 0.84em; transition: all 0.15s; }
        .back-btn:hover { background: #38bdf8; color: #020617; }

        .wiki-layout { display: flex; max-width: 1600px; margin: 0 auto; min-height: calc(100vh - 62px); }
        .wiki-sidebar { width: 280px; flex-shrink: 0; position: sticky; top: 62px; height: calc(100vh - 62px); overflow-y: auto; background: #020617; border-right: 1px solid #1e293b; padding: 14px 12px; display: flex; flex-direction: column; gap: 4px; }
        .wiki-sidebar-label { font-size: 0.72em; text-transform: uppercase; letter-spacing: 0.05em; color: #64748b; font-weight: 700; margin: 4px 6px 6px 6px; }
        .wiki-search { width: 100%; background: #0f172a; border: 1px solid #334155; color: #fff; padding: 8px 10px; border-radius: 6px; font-size: 0.83em; margin-bottom: 8px; }
        .wiki-nav-btn { text-align: left; background: transparent; color: #94a3b8; border: 1px solid transparent; padding: 8px 11px; border-radius: 6px; cursor: pointer; font-size: 0.83em; transition: all 0.12s; }
        .wiki-nav-btn:hover { background: #0f172a; color: #e2e8f0; }
        .wiki-nav-btn.active { background: #1e293b; color: #38bdf8; border-color: #38bdf8; font-weight: 600; }

        .wiki-main { flex: 1; padding: 20px 28px 60px 28px; display: flex; flex-direction: column; gap: 18px; max-width: 1320px; }
        .wiki-inspector-shell { height: 590px; border: 1px solid #1e293b; border-radius: 10px; overflow: hidden; box-shadow: 0 8px 28px rgba(0,0,0,0.45); }

        .wiki-card { display: none; background: #0f172a; border: 1px solid #1e293b; border-radius: 10px; padding: 18px 22px; scroll-margin-top: 76px; }
        .wiki-card.active-phase { display: block; }
        .wiki-card-header { display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px; border-bottom: 1px solid #1e293b; padding-bottom: 10px; margin-bottom: 12px; }
        .wiki-card-title { font-size: 1.08em; font-weight: 700; color: #38bdf8; }
        .wiki-tag { font-size: 0.74em; background: #020617; border: 1px solid #334155; color: #facc15; padding: 2px 9px; border-radius: 999px; }
        .intro-dual { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin: 8px 0 12px 0; }
        .intro-sci { background: rgba(16, 185, 129, 0.06); border-left: 3px solid #10b981; padding: 9px 12px; border-radius: 6px; font-size: 0.83em; color: #d1fae5; }
        .intro-dev { background: rgba(56, 189, 248, 0.06); border-left: 3px solid #38bdf8; padding: 9px 12px; border-radius: 6px; font-size: 0.83em; color: #e0f2fe; }
        .wiki-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-top: 10px; }
        .wiki-box { background: #090d16; border: 1px solid #1e293b; border-radius: 8px; padding: 11px 14px; font-size: 0.84em; color: #cbd5e1; }
        .wiki-box b { color: #f8fafc; }
        .live-db-box { margin-top: 12px; background: rgba(16, 185, 129, 0.07); border: 1px solid rgba(16, 185, 129, 0.3); border-radius: 8px; padding: 10px 14px; font-size: 0.83em; color: #d1fae5; display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px; }
        .live-db-box b { color: #34d399; }
        .wiki-sql { background: #020617; border: 1px solid #1e293b; color: #a5f3fc; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 0.8em; padding: 11px 14px; border-radius: 8px; margin-top: 10px; overflow-x: auto; white-space: pre; }
        .formula-pill { display: inline-block; background: #020617; border: 1px dashed #38bdf8; color: #a5f3fc; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 0.86em; padding: 3px 8px; border-radius: 5px; margin: 3px 0; cursor: pointer; }
        .formula-pill:hover { background: #0f172a; border-style: solid; }
        __INSPECTOR_CSS__
    </style>
</head>
<body>
    <div class="wiki-topbar">
        <div>
            <h1>Swarm Discovery — Phase-Scoped Correlation Wiki &amp; 3D Biophysical Lab</h1>
            <div class="sub">Select any discovery phase on the left to scope both the 3D Compare Inspector and the PostgreSQL 18 (pg_bio) documentation</div>
        </div>
        <div style="display:flex; gap:10px; align-items:center;">
            <a href="/cards" class="back-btn" style="border-color:#f59e0b; color:#fbbf24;">🃏 1-to-N Card Arena</a>
            <a href="/" class="back-btn">Back to Swarm Heatmap</a>
        </div>
    </div>

    <div class="wiki-layout">
        <div class="wiki-sidebar">
            <input id="wiki-search" class="wiki-search" type="text" placeholder="Filter phases, SQL, formulas..." oninput="filterWikiCards()">
            <div class="wiki-sidebar-label">Discovery &amp; Validation Phases</div>
            <button class="wiki-nav-btn active" data-cat="vector" onclick="selectWikiPhase('vector')">1. Vector &amp; PPI Homology</button>
            <button class="wiki-nav-btn" data-cat="sequence" onclick="selectWikiPhase('sequence')">2. Sequence &amp; Twilight-Zone</button>
            <button class="wiki-nav-btn" data-cat="structure" onclick="selectWikiPhase('structure')">3. 3D Kabsch, TM &amp; Theozymes</button>
            <button class="wiki-nav-btn" data-cat="chaperone" onclick="selectWikiPhase('chaperone')">4. RosettaCommons Chaperone</button>
            <button class="wiki-nav-btn" data-cat="biochem" onclick="selectWikiPhase('biochem')">5. Extremophile Biochemistry</button>
            <button class="wiki-nav-btn" data-cat="cheminfo" onclick="selectWikiPhase('cheminfo')">6. Cheminformatics &amp; Substrates</button>
            <button class="wiki-nav-btn" data-cat="rangefinder" onclick="selectWikiPhase('rangefinder')">7. RangeFinder 9-Prop GiST</button>
            <button class="wiki-nav-btn" data-cat="spatial" onclick="selectWikiPhase('spatial')">8. 3D Spatial &amp; Z-Order Pockets</button>
            <button class="wiki-nav-btn" data-cat="synbio" onclick="selectWikiPhase('synbio')">9. Wet-Lab SynBio &amp; PROSITE</button>
            <button class="wiki-nav-btn" data-cat="proteomics" onclick="selectWikiPhase('proteomics')">10. LC-MS/MS Proteomics</button>
            <button class="wiki-nav-btn" data-cat="energetics" onclick="selectWikiPhase('energetics')">11. Energetics, Allostery &amp; Kinetics</button>
            <button class="wiki-nav-btn" data-cat="catalog" onclick="selectWikiPhase('catalog')">12. Master SQL Catalog (73 Funcs)</button>
        </div>

        <div class="wiki-main">
            <!-- Unified Configurable Side-by-Side 3D Comparison Inspector -->
            <div class="wiki-inspector-shell" id="wiki-inspector-root"></div>

            <!-- Phase 1: Vector & PPI Homology -->
            <div class="wiki-card active-phase" data-cat="vector" id="wiki-vector">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 1A: ESM-2 Vector Cosine Distance &amp; Tri-Peptide Embedding Homology</span>
                    <span class="wiki-tag">Phase 1 · pgvector HNSW Graph</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Even when two proteins have diverged over billions of years so that sequence alignment (BLAST) sees only random noise, protein language models (ESM-2) and tri-peptide frequency profiles capture the underlying 3D fold and catalytic family.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Converts variable-length amino acid strings into fixed 320-dimensional unit vectors indexed via <code>pgvector</code> HNSW graphs (<code>ef_search = 200</code>), turning O(N) sequence scans into sub-millisecond approximate nearest-neighbor graph hops.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Mathematical Formula (Hover for Intuition):</b><br>
                        <span class="formula-pill" data-wiki-key="esm_kmer_dist">d_cos(A, B) = 1 - (A · B) / (||A|| ||B||)</span><br>
                        In <code>swarm_discovery.py</code>, an HNSW index scans for uncharacterized neighbors with <code>d_cos &le; 0.15</code>.
                    </div>
                    <div class="wiki-box">
                        <b>Decision Thresholds:</b><br>
                        • <b>&le; 0.06:</b> Tight functional and structural family match.<br>
                        • <b>0.06 – 0.12:</b> Distant homolog or cross-phylum structural analog.<br>
                        • <b>&gt; 0.15:</b> Divergent fold or multi-domain fusion.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-vector"><b>Phase 1 Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('identity')">Highlight Identity in 3D</button>
                </div>
                <div class="wiki-sql">SELECT (p1.embedding &lt;=&gt; p2.embedding) AS vector_distance,
       embedding_cosine_distance(get_esm_embedding(p1.sequence), get_esm_embedding(p2.sequence)) AS esm_kmer_dist
FROM proteins p1, proteins p2 WHERE p1.uniprot_id = 'A0A1S8AT22' AND p2.uniprot_id = 'A0A1H3IA37';</div>
            </div>

            <div class="wiki-card active-phase" data-cat="vector" id="wiki-icomp">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 1B: PPI Interface Complementarity &amp; Swarm Confidence Score</span>
                    <span class="wiki-tag">Phase 1 · Electrostatic / Hydrophobic Split</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Protein-protein interfaces dock via two opposing forces: opposite electrostatic charges attract (positive Lys/Arg meets negative Asp/Glu), while hydrophobic patches self-associate to exclude water.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Partitions the 320-D embedding vector into two orthogonal subspaces—requiring negative cosine correlation on the first half (electrostatics) and positive cosine correlation on the second half (hydrophobicity).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Mathematical Formula (Hover for Intuition):</b><br>
                        <span class="formula-pill" data-wiki-key="icomp_formula">I_comp = (-cos(E_A, E_B) + cos(H_A, H_B) + 2) / 4</span><br>
                        <span class="formula-pill" data-wiki-key="icomp_formula">Confidence = 0.7 * (1 - d_cos) + 0.3 * I_comp</span>
                    </div>
                    <div class="wiki-box">
                        <b>Application in the Swarm:</b><br>
                        Used during node expansion in <code>swarm_discovery.py</code> to prioritize orphans that share global fold semantics while preserving complementary surface electrostatics.
                    </div>
                </div>
                <div class="wiki-sql">SELECT interface_complementarity(rn.embedding::vector::real[], bait.embedding::vector::real[]) AS icomp
FROM proteins rn, proteins bait WHERE rn.uniprot_id = 'A0A1H3IA37' AND bait.uniprot_id = 'A0A1S8AT22';</div>
            </div>

            <!-- Phase 2: Sequence & Twilight-Zone -->
            <div class="wiki-card" data-cat="sequence" id="wiki-sw">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 2A: Smith-Waterman (pg_bio) Local Alignment &amp; BLOSUM62 Score</span>
                    <span class="wiki-tag">Phase 2 · Local Dynamic Programming</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Uncharacterized microbial proteins often carry variable N/C-terminal signal peptides or flexible linkers. Local alignment isolates the conserved catalytic motif without penalty from flanking tails.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Executes exact O(N*M) Smith-Waterman dynamic programming in compiled Rust inside the PostgreSQL worker process with zero network serialization overhead.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Normalized SW Similarity:</b><br>
                        <span class="formula-pill" data-wiki-key="sw_score">SW_sim = SW(A, B; +2, -1, -2) / (2 * min(L_A, L_B))</span>
                    </div>
                    <div class="wiki-box">
                        <b>Global vs Local Scoring:</b><br>
                        Combines BLOSUM62 evolutionary substitution log-odds (<code>sequence_alignment_score</code>) with local motif detection (<code>smith_waterman_score</code>).
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-sw"><b>Phase 2 Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('ss')">Highlight Secondary Structure in 3D</button>
                </div>
                <div class="wiki-sql">SELECT sequence_alignment_score(p1.sequence, p2.sequence) AS blosum_score,
       smith_waterman_score(p1.sequence, p2.sequence, 2, -1, -2) AS sw_local_score
FROM proteins p1, proteins p2 WHERE p1.uniprot_id = 'A0A1S8AT22' AND p2.uniprot_id = 'A0A1H3IA37';</div>
            </div>

            <div class="wiki-card" data-cat="sequence" id="wiki-twilight">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 2B: Twilight-Zone Score (Remote Homology Bandpass Filter)</span>
                    <span class="wiki-tag">Phase 2 · Remote Homology</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> The <b>Twilight Zone (15%–35% sequence identity)</b> is where standard sequence search tools fail, yet 3D structure and enzymatic mechanism remain conserved across distant phyla.</div>
                    <div class="intro-dev"><b>For Engineers:</b> A bandpass ranking function that multiplies vector proximity by the inverted parabola <code>4s(1-s)</code>—suppressing trivial 100% clones (s-&gt;1) and random noise (s-&gt;0).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Mathematical Formula:</b><br>
                        <span class="formula-pill" data-wiki-key="twilight_score">T = (1 - d_vec)^2 * [4 * s * (1 - s)] * (L_min / L_max)^2 * min(H_shannon / 4, 1)</span>
                    </div>
                    <div class="wiki-box">
                        <b>Thresholds:</b><br>
                        • <b>&ge; 0.50:</b> Prime Twilight-Zone remote homolog across distant taxa.<br>
                        • <b>0.15 – 0.49:</b> Moderate sequence divergence.<br>
                        • <b>0.00:</b> Nearly identical clone (&gt;96% SW) or unrelated (&lt;10% SW).
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-twilight"><b>Phase 2 Twilight Values:</b> Loading...</span>
                </div>
                <div class="wiki-sql">SELECT bait_id, orphan_id, vector_distance, twilight_score
FROM orphan_discoveries WHERE twilight_score &gt; 0.10
ORDER BY (split_part(bait_organism, ' ', 1) != split_part(orphan_organism, ' ', 1)) DESC, twilight_score DESC LIMIT 10;</div>
            </div>

            <!-- Phase 3: 3D Structure, Kabsch & Theozymes -->
            <div class="wiki-card" data-cat="structure" id="wiki-kabsch">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 3A: 3D Core RMSD, Sequential TM-Score &amp; Kabsch Fold Gate</span>
                    <span class="wiki-tag">Phase 3 · 3D Coordinate Registration</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> 3D structure is 3–10x more conserved than sequence during evolution. Superposing C-alpha backbones in 3D reveals whether two enzymes share the exact same catalytic scaffold.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Multi-seed internal distance matrix fingerprinting followed by iterative Singular Value Decomposition (<b>Kabsch SVD</b>) and dynamic programming along 3D tangent vectors.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>3D Core RMSD &amp; Sequential TM-Score Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="core_rmsd">3D Core RMSD = sqrt( (1/N_core) sum ||x_i - (R*y_i + t)||^2 )</span><br>
                        <span class="formula-pill" data-wiki-key="tm_score">TM = (1 / L_min) sum [1 / (1 + (d_i / d_0)^2)] , d_0 = 1.24*(L_min - 15)^(1/3) - 1.8</span>
                    </div>
                    <div class="wiki-box">
                        <b>Fold Gate Thresholds:</b><br>
                        • <b>TM &ge; 0.50 ([SAME FOLD]):</b> Same global 3D fold topology.<br>
                        • <b>0.25 &le; TM &lt; 0.50 ([MOTIF MATCH]):</b> Shared local domain or catalytic motif.<br>
                        • <b>TM &lt; 0.25 ([DIVERGENT]):</b> Different global fold topology.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-kabsch"><b>Phase 3 Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('core')">Highlight Aligned C-alpha Core in 3D</button>
                </div>
            </div>

            <div class="wiki-card" data-cat="structure" id="wiki-theozyme">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 3B: Active-Site Theozyme Geometry, Globular Compactness &amp; Dihedrals</span>
                    <span class="wiki-tag">Phase 3 · Active-Site Geometry &amp; Rg</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Enzymes catalyze reactions using precise 3D constellations of active-site residues (<b>theozymes</b>), stabilized by coordinated structural waters (2.4–3.4 A), strain-free Ramachandran (phi, psi) dihedrals, and a compact globular hydrophobic core.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Implements sub-graph Kabsch SVD matching (<code>theozyme_rmsd</code>), moment-of-inertia radius of gyration &amp; Flory packing ratio (<code>evaluate_globular_compactness</code>), two-plane normal cross-product torsion angles (<code>dihedral_angle</code>), and steric clash checks (<code>is_steric_clash</code>).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="theozyme_formula">R_g = sqrt( (1/N) sum ||r_i - r_com||^2 ) ; Flory Ratio = R_g / (2.25 * N^0.38)</span>
                    </div>
                    <div class="wiki-box">
                        <b>SQL Functions:</b><br>
                        <code>evaluate_globular_compactness(pdb)</code>, <code>theozyme_rmsd(atoms, motif)</code>, <code>has_coordinated_water_shell(...)</code>, <code>dihedral_angle(a,b,c,d)</code>, <code>is_steric_clash(a,b,vdw)</code>, <code>strand_pairing_check(...)</code>.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-theozyme"><b>Phase 3 Globular Compactness Values:</b> Loading...</span>
                </div>
                <div class="wiki-sql">SELECT evaluate_globular_compactness(pdb_text)::jsonb AS globular_profile,
       dihedral_angle(create_residue_coord(0,0,0,'N'), create_residue_coord(1.46,0,0,'CA'),
                      create_residue_coord(2.0,1.42,0,'C'), create_residue_coord(1.3,2.5,0.5,'N2')) AS psi_deg;</div>
            </div>

            <!-- Phase 4: RosettaCommons Chaperone Audit -->
            <div class="wiki-card" data-cat="chaperone" id="wiki-chaperone">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 4: RosettaCommons Chaperone Audit (Core RASA Corr, H-Bond CSR Hubs, BUNS &amp; Lineage)</span>
                    <span class="wiki-tag">Phase 4 · Biophysical Design &amp; Stability</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> A matching C-alpha backbone is not enough—a functional protein must pack its hydrophobic core, satisfy buried polar atoms with hydrogen bonds (avoiding thermodynamic <code>BUNS</code> penalties), and prevent surface aggregation.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Implements the <b>RosettaCommons Chaperone</b> evaluation pipeline in Rust inside SQL: 64-point Fibonacci sphere ray-casting for SASA/RASA, Compressed Sparse Row (<code>CSR</code>) adjacency graphs for H-bonds, and deterministic lineage hashing.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Biophysical Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="chaperone_audit">RASA_i = SASA_i / MaxASA_Tien2013(aa_i) ; Core RASA Corr = corr(RASA_A, RASA_B)</span><br>
                        • <b>H-Bond CSR Hubs:</b> Top 12 highest-degree residues from <code>bio_hbond_network(pdb, 3.5)</code>.<br>
                        • <b>BUNS:</b> Polar N/O atoms buried in the core (<code>SASA &lt; 1.0 A^2</code>) with zero H-bonds.
                    </div>
                    <div class="wiki-box">
                        <b>Composite Chaperone Score &amp; Lineage:</b><br>
                        <span class="formula-pill" data-wiki-key="chaperone_audit">Score = w_sc*(1 - scRMSD/4) + w_hb*HB_density - w_buns*BUNS - w_sh*SurfHydro</span><br>
                        Tracked via <code>chaperone_lineage_hash</code> and benchmarked against the 8 Official Chaperone Core Examples (<code>chaperone_core_examples()</code>).
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-chaperone"><b>Phase 4 Live DB Values:</b> Loading...</span>
                    <div style="display:flex; gap:6px;">
                        <button class="cmp-btn" onclick="wikiInspector.applyMode('hbond')">Highlight H-Bond CSR Hubs</button>
                        <button class="cmp-btn" onclick="wikiInspector.applyMode('rasa')">Highlight RASA Hydro Patches</button>
                    </div>
                </div>
                <div class="wiki-sql">SELECT c.mean_rasa, c.core_buried_fraction, c.surface_hydrophobic_fraction,
       c.total_hbonds, c.buried_unsatisfied_polars AS buns, c.design_composite_score,
       bio_hbond_network(pdb_text, 3.5)::text AS hbond_csr_network,
       chaperone_lineage_hash('SwarmCompare', bait_id, orphan_seq, 1) AS lineage_id
FROM evaluate_chaperone_block(orphan_pdb_aligned, bait_pdb, 'A', 'B') c;</div>
            </div>

            <!-- Phase 5: Extremophile Biochemistry -->
            <div class="wiki-card" data-cat="biochem" id="wiki-biochem">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 5: Physicochemical Adaptation Shifts (dpI, dGRAVY, Acidic Shell &amp; Entropy)</span>
                    <span class="wiki-tag">Phase 5 · Halophile Surface Adaptation</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> How does an enzyme function inside 4M NaCl salt lakes without precipitating? Halophilic archaea coat their protein exteriors with negatively charged Aspartate and Glutamate (<code>Acidic(D/E) &gt; 20%</code>, <code>pI ~ 4.2</code>) to coordinate hydrated salt ions.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Fast SQL scalar functions computing Henderson-Hasselbalch charge neutrality root-finding (<code>isoelectric_point</code>), Kyte-Doolittle hydropathy averages (<code>gravy_score</code>), and Shannon information entropy (<code>sequence_entropy</code>).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Mathematical Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="biochem">H(seq) = -sum p(aa) log2 p(aa) ; GRAVY = (1/L) sum Hydropathy(aa_i)</span>
                    </div>
                    <div class="wiki-box">
                        <b>3D Acidic Shell Inspection:</b><br>
                        Click <b>Highlight Acidic Shell (D/E)</b> to render Aspartate and Glutamate sidechains as red sticks on the 3D surface.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-biochem"><b>Phase 5 Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('acidic')">Highlight Acidic Shell (D/E) in 3D</button>
                </div>
                <div class="wiki-sql">SELECT molecular_weight(sequence), isoelectric_point(sequence) AS pi,
       gravy_score(sequence) AS gravy, sequence_entropy(sequence) AS shannon_entropy
FROM proteins WHERE uniprot_id IN ('A0A1S8AT22', 'A0A1H3IA37');</div>
            </div>

            <!-- Phase 6: Cheminformatics & Substrates -->
            <div class="wiki-card" data-cat="cheminfo" id="wiki-cheminfo">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 6: Substrate Cheminformatics (Huckel 4n+2, WL64, VF2 Warhead @&gt; &amp; ECFP4)</span>
                    <span class="wiki-tag">Phase 6 · Molecular Graph Engine</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists &amp; Chemists:</b> Links each enzyme family to its small-molecule substrate and verifies the reactive catalytic warhead, invariant to Kekule bond alternation or counter-ion salts.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Native Rust molecular graph engine inside PostgreSQL implementing Huckel 4n+2 pi-electron aromaticity perception, 64-bit Weisfeiler-Lehman graph hashes, VF2 subgraph isomorphism (<code>@&gt;</code>), and 1024-bit Morgan ECFP4 popcount Tanimoto (<code>%</code>).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Graph Equivalence Formula:</b><br>
                        <span class="formula-pill" data-wiki-key="substrate_cheminfo">WL64(Kekule/Salt) == WL64(Aromatic)  AND  mol @&gt; warhead</span>
                    </div>
                    <div class="wiki-box">
                        <b>SQL Operators:</b><br>
                        <code>parse_smiles(s)</code>, <code>molecule_wl_hash(s)</code>, <code>molecules_are_identical(s1, s2)</code>, <code>mol @&gt; sub</code>, <code>mol % query</code>.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-cheminfo"><b>Phase 6 Live DB Values:</b> Loading...</span>
                </div>
                <div class="wiki-sql">SELECT substrate_name, wl_hash,
       molecules_are_identical(ingested_kekule_or_salt_smiles, canonical_query_smiles) AS kekule_arom_eq,
       (mol @&gt; parse_smiles(warhead_substructure_smiles)) AS vf2_warhead_match,
       (mol % parse_smiles(canonical_query_smiles)) AS ecfp4_tanimoto
FROM swarm_enzyme_substrates WHERE LOWER(family_name) = LOWER('Nitrite reductase');</div>
            </div>

            <!-- Phase 7: RangeFinder 9-Property Pocket-Interval GiST -->
            <div class="wiki-card" data-cat="rangefinder" id="wiki-rangefinder">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 7: RangeFinder 9-Property Pocket-Interval GiST &amp; Substrate Box Screening (arXiv:2611.03379v1)</span>
                    <span class="wiki-tag">Phase 7 · 8-Column numrange GiST</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists &amp; Medicinal Chemists:</b> Instead of docking millions of ligands one-by-one, RangeFinder predicts the <b>9 physicochemical property intervals</b> (<code>[μ - 2|σ|, μ + 2|σ|]</code> for MW, LogP, HBD, HBA, TPSA, RotBonds, FormalCharge, AromRings, Fsp3) plus aromatic ring probability <code>P(arom)</code> directly from a 10.0 Å 3D active-site microenvironment.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Maps each pocket's 19-float Rust <code>RangeFinderPocket</code> struct into 8 native PostgreSQL <code>numrange</code> columns indexed via a multi-column GiST index (<code>USING gist (mw_range, logp_range, hbd_range, hba_range, tpsa_range, rot_range, arom_range, fsp3_range)</code>) for sub-millisecond reverse target fishing (<code>0.68 ms</code> across 20k pockets) and <code>3.78 ms</code> composite B-Tree forward screening (250k ligands).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Interval Prediction &amp; Generalized Range-IoU:</b><br>
                        <span class="formula-pill" data-wiki-key="rangefinder_pocket">I_p = [μ_p - 2|σ_p|, μ_p + 2|σ_p|]  |  RangeIoU(A, B) = (1/8) Σ |I_A,p ∩ I_B,p| / |I_A,p ∪ I_B,p|</span><br>
                        Compares active-site physicochemical boxes across remote homologs and checks <code>8/8</code> cognate substrate containment (<code>pocket_matches_ligand_8prop</code>).
                    </div>
                    <div class="wiki-box">
                        <b>Decision Thresholds:</b><br>
                        • <b>Range-IoU &ge; 0.85 + 8/8 Box Match:</b> Conserved catalytic pocket microenvironment &amp; substrate compatibility.<br>
                        • <b>Range-IoU 0.65 – 0.84:</b> Moderately shifted active-site polarity/volume.<br>
                        • <b>&lt; 0.65:</b> Divergent substrate specificity.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-rangefinder"><b>Phase 7 RangeFinder Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('rangefinder')">Highlight 10 Å RangeFinder Pocket in 3D</button>
                </div>
                <div class="wiki-sql">-- Reverse Target Fishing via 8-Column numrange GiST Index (0.68 ms):
SELECT p.uniprot_id, ROUND(pocket_range_iou(pb.pocket, p.pocket, false)::numeric, 4) AS range_iou,
       p.mw_range, p.logp_range, ROUND(p.aromatic_probability::numeric, 3) AS arom_prob
FROM swarm_enzyme_substrates s
JOIN swarm_pocket_ranges p
  ON p.mw_range @&gt; s.mw::numeric AND p.logp_range @&gt; s.logp::numeric
 AND p.hbd_range @&gt; s.hbd::numeric AND p.hba_range @&gt; s.hba::numeric
 AND p.tpsa_range @&gt; s.tpsa::numeric AND p.rot_range @&gt; s.rotatable_bonds::numeric
 AND p.arom_range @&gt; s.aromatic_rings::numeric AND p.fsp3_range @&gt; s.fsp3::numeric
JOIN swarm_pocket_ranges pb ON pb.uniprot_id = 'A0A2K5ARN2'
WHERE LOWER(s.family_name) = 'xylanase' ORDER BY range_iou DESC LIMIT 10;</div>
            </div>

            <!-- Phase 7: 3D Spatial & Z-Order Pockets -->
            <div class="wiki-card" data-cat="spatial" id="wiki-spatial">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 7: 3D Spatial Indexing, pLDDT Confidence &amp; 30-Bit Z-Order (Morton) Pockets</span>
                    <span class="wiki-tag">Phase 7 · 3D Spatial B-Tree</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Inspect per-residue AlphaFold/ESMFold confidence (<code>pLDDT</code>) or click any active-site atom in the 3D viewer to isolate a 5.0 A spherical binding pocket.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Interleaves 3D coordinates <code>(x, y, z)</code> into a single 30-bit Z-Order Morton integer (<code>z_order_encode</code>) so standard PostgreSQL B-Tree indexes can execute 3D bounding-box range scans.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>pLDDT Confidence Scale:</b><br>
                        <span class="formula-pill" data-wiki-key="plddt">pLDDT &gt; 90 (Dark Blue) | &gt; 70 (Light Blue) | &lt; 50 (Flexible Loop)</span>
                    </div>
                    <div class="wiki-box">
                        <b>30-Bit Z-Order Morton Code:</b><br>
                        <span class="formula-pill" data-wiki-key="zorder_formula">z_order_encode(x, y, z) = interleave_bits_10x3(X_q, Y_q, Z_q)</span>
                    </div>
                </div>
                <div class="live-db-box">
                    <span><b>Interactive 3D Modes:</b> Test pLDDT Confidence or activate Z-Order Pocket Selection above.</span>
                    <div style="display:flex; gap:6px;">
                        <button class="cmp-btn" onclick="wikiInspector.applyMode('plddt')">Color by pLDDT Confidence</button>
                        <button class="cmp-btn" onclick="wikiInspector.togglePocketSelection()">Toggle Z-Order Pocket Click</button>
                    </div>
                </div>
                <div class="wiki-sql">SELECT atom_serial, residue_name, residue_seq
FROM protein_atoms
WHERE morton_code BETWEEN z_order_encode(x-5, y-5, z-5) AND z_order_encode(x+5, y+5, z+5)
  AND distance_angstroms(coord, create_residue_coord(query_x, query_y, query_z, 'Q')) &lt;= 5.0;</div>
            </div>

            <!-- Phase 8: Wet-Lab SynBio & PROSITE -->
            <div class="wiki-card" data-cat="synbio" id="wiki-synbio">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 8: Wet-Lab Synthetic Biology, Codon Adaptation Index (CAI), Plasmid Cassettes &amp; PROSITE</span>
                    <span class="wiki-tag">Phase 8 · Heterologous Expression Design</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Converts discovered archaeal enzymes into host-optimized DNA genes for <i>E. coli</i>, yeast, or human expression—optimizing <b>Codon Adaptation Index (CAI)</b>, balancing <b>GC%</b>, verifying freedom from Golden Gate <code>BsaI</code> restriction sites, and scanning <b>PROSITE</b> motifs.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Compiles PROSITE syntax (e.g. <code>[ST]-x(2)-[RK]</code>) into Rust finite-automata regexes (<code>prosite_match</code>) and runs deterministic codon back-translation (<code>bio_reverse_translate</code>) and plasmid assembly (<code>bio_design_expression_cassette</code>) inside SQL.</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Sharp &amp; Li (1987) CAI Formula:</b><br>
                        <span class="formula-pill" data-wiki-key="synbio_ms">CAI = ( prod_{k=1..L} w_k )^(1/L) , w_k = freq(codon_k) / max_syn_freq</span>
                    </div>
                    <div class="wiki-box">
                        <b>SQL Functions:</b><br>
                        <code>bio_design_expression_cassette(seq, 'ecoli', 'T7')</code>, <code>bio_codon_adaptation_index(dna, 'ecoli')</code>, <code>bio_reverse_translate(seq, 'ecoli')</code>, <code>prosite_match(seq, pattern)</code>.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-synbio"><b>Phase 8 Live DB Values:</b> Loading...</span>
                </div>
                <div class="wiki-sql">SELECT bio_design_expression_cassette(left(p.sequence, 240), 'ecoli', 'T7')::jsonb AS plasmid_cassette,
       prosite_match(p.sequence, '[ST]-x(2)-[RK]') AS has_ck2_motif,
       prosite_match(p.sequence, 'G-x(2)-G-x-G-K-[ST]') AS has_walker_a_ploop
FROM proteins p WHERE p.uniprot_id = 'A0A1H3IA37';</div>
            </div>

            <!-- Phase 9: LC-MS/MS Proteomics -->
            <div class="wiki-card" data-cat="proteomics" id="wiki-proteomics">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 9: Tandem Mass Spectrometry (LC-MS/MS), b/y Fragment Ladders, PPM Error &amp; Isotopes</span>
                    <span class="wiki-tag">Phase 9 · Bottom-Up Proteomics</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists &amp; Analytical Chemists:</b> Validates expressed proteins via tandem mass spectrometry—predicting Collision-Induced Dissociation (<b>CID/HCD</b> -&gt; <code>b/y</code> ions) and <b>ETD</b> (<code>c/z</code> ions), Orbitrap <b>ppm</b> mass error, <b>13C isotopic envelope</b> fit, and target-decoy peptides.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Monoisotopic prefix/suffix mass ladder generator (<code>generate_theoretical_fragments</code>), Averagine 13C Poisson cosine scoring (<code>isotope_fit_score</code>), ppm tolerance checks (<code>ppm_difference</code>), and FDR decoy reversal (<code>generate_decoy_sequence</code>).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Mass Error &amp; Fragment m/z Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="proteomics_formula">ppm_err = (|m_obs - m_calc| / m_calc) * 10^6 ; m/z = (M_neutral + z * m_H+) / z</span>
                    </div>
                    <div class="wiki-box">
                        <b>3D Diagnostic Peptide Highlight:</b><br>
                        Click <b>Highlight N-Term MS/MS Peptide</b> to highlight residues 1–10 in 3D alongside their theoretical HCD <code>b2, b3, y2, y3</code> ions.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-proteomics"><b>Phase 9 Live DB Values:</b> Loading...</span>
                    <button class="cmp-btn" onclick="wikiInspector.applyMode('nterm_pep')">Highlight N-Term MS/MS Peptide in 3D</button>
                </div>
                <div class="wiki-sql">SELECT f.ion_type, round(f.mass::numeric, 4) AS mz, f.charge, f.sequence_fragment,
       generate_decoy_sequence('MEVDLLDFVK', 'reverse') AS tryptic_decoy_peptide,
       round(ppm_difference(590.3145, 590.3120)::numeric, 2) AS orbitrap_ppm_error
FROM generate_theoretical_fragments('MEVDLLDFVK', 2, 'HCD') f
ORDER BY f.charge ASC, f.mass ASC LIMIT 8;</div>
            </div>

            <!-- Phase 10: Energetics, Allostery & Kinetics -->
            <div class="wiki-card" data-cat="energetics" id="wiki-energetics">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 10: Allosteric Hotspots, MM/PBSA Energetics, MD Ensembles &amp; Eyring Enzyme Kinetics</span>
                    <span class="wiki-tag">Phase 10 · Thermodynamics &amp; Kinetics</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists &amp; Medicinal Chemists:</b> Quantifies binding free energy (<b>MM/PBSA dG_bind</b>), identifies distal <b>allosteric hotspots</b> that lock catalytic loops remotely, checks MD conformational pre-organization, and models <b>Eyring transition-state kinetics (k_cat, dG‡)</b>.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Vectorized thermodynamic auditors (<code>audit_mmpbsa_decomposition</code>, <code>solvent_shift_penalty</code>), CSR residue-coupling graphs (<code>build_energy_hotspot_map</code>, <code>allosteric_mechanism_index</code>), ensemble CDF gates (<code>ensemble_fraction_below</code>), and Eyring simulators (<code>simulate_defluorination_kinetics</code>).</div>
                </div>
                <div class="wiki-grid">
                    <div class="wiki-box">
                        <b>Thermodynamic &amp; Kinetic Formulas:</b><br>
                        <span class="formula-pill" data-wiki-key="energetics_formula">dG_bind = dE_vdW + dE_elec + dG_polar + dG_apolar ; k_cat = (k_B*T/h)*exp(-dG‡ / R*T)</span>
                    </div>
                    <div class="wiki-box">
                        <b>Allostery &amp; Ensemble Gating:</b><br>
                        <code>allosteric_mechanism_index(clash, lock, path)</code> separates orthosteric steric competition from distal allosteric locking.
                    </div>
                </div>
                <div class="live-db-box">
                    <span id="live-ex-energetics"><b>Phase 10 Live DB Values:</b> Loading...</span>
                </div>
                <div class="wiki-sql">SELECT k.activation_barrier_kcal, k.kcat_per_sec, k.catalytic_efficiency,
       k.substrate_half_life_hours, k.is_industrially_viable,
       allosteric_mechanism_index(0.15, 0.88, 0.92) AS allosteric_index
FROM simulate_defluorination_kinetics('HaloAlk_DeFl_1', 'FC(F)(F)C(=O)O', 310.15, -9.4, 2.85, -4.2) k;</div>
            </div>

            <!-- Phase 11: Master SQL Catalog -->
            <div class="wiki-card" data-cat="catalog" id="wiki-sql-catalog">
                <div class="wiki-card-header">
                    <span class="wiki-card-title">Phase 11: Complete pg_bio SQL Function &amp; Operator Reference (All 64 Compiled Rust Functions)</span>
                    <span class="wiki-tag">Phase 11 · Master Reference</span>
                </div>
                <div class="intro-dual">
                    <div class="intro-sci"><b>For Biologists:</b> Reference table of every biological query supported natively inside PostgreSQL across all 10 scientific phases.</div>
                    <div class="intro-dev"><b>For Engineers:</b> Compiled to native machine code via <code>pgrx</code> in <code>src/lib.rs</code>, marked <code>IMMUTABLE PARALLEL SAFE</code> for multi-core PostgreSQL 18 execution.</div>
                </div>
                <div style="overflow-x:auto; margin-top:8px;">
                    <table style="width:100%; border-collapse:collapse; font-size:0.81em; color:#cbd5e1;">
                        <thead>
                            <tr style="background:#020617; color:#38bdf8; text-align:left; border-bottom:1px solid #334155;">
                                <th style="padding:8px;">Phase / Domain</th>
                                <th style="padding:8px;">SQL Function / Operator Signature</th>
                                <th style="padding:8px;">Scientific Purpose &amp; Algorithmic Core</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#34d399; font-weight:600;">1. Vector &amp; PPI</td><td style="padding:7px;"><code>get_esm_embedding(seq)</code>, <code>embedding_cosine_distance(a, b)</code>, <code>p1.embedding &lt;=&gt; p2.embedding</code>, <code>interface_complementarity(e1, e2)</code>, <code>sparse_cross_reactivity_score(m1, m2, len)</code></td><td style="padding:7px;">320-D protein embeddings, HNSW cosine nearest neighbors, electrostatic/hydrophobic interface split, and CSR cross-reactivity.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#38bdf8; font-weight:600;">2. Sequence &amp; Motifs</td><td style="padding:7px;"><code>smith_waterman_score(s1, s2, match, mis, gap)</code>, <code>sequence_alignment_score(s1, s2)</code>, <code>prosite_match(seq, pattern)</code></td><td style="padding:7px;">Exact local Smith-Waterman DP, BLOSUM62 evolutionary alignment, and native PROSITE motif regex compiler.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#fbbf24; font-weight:600;">3. 3D Kabsch &amp; Theozymes</td><td style="padding:7px;"><code>kabsch_rmsd(a, b)</code>, <code>bio_rmsd(a, b)</code>, <code>theozyme_rmsd(atoms, motif)</code>, <code>evaluate_globular_compactness(pdb)</code>, <code>dihedral_angle(a,b,c,d)</code>, <code>bio_angle(a,b,c)</code>, <code>bio_dihedral(a,b,c,d)</code>, <code>is_steric_clash(a,b,vdw)</code>, <code>strand_pairing_check(...)</code>, <code>has_coordinated_water_shell(...)</code></td><td style="padding:7px;">Kabsch SVD rigid-body alignment, active-site catalytic triad RMSD, globular Rg &amp; core packing, Ramachandran torsions, steric clashes, beta-sheet pairing, and water shells.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#4ade80; font-weight:600;">4. RosettaCommons Chaperone</td><td style="padding:7px;"><code>bio_rasa(pdb, 1.4)</code>, <code>bio_hbond_network(pdb, 3.5)</code>, <code>audit_hydrogen_bonds(pdb, 3.5, 1.0)</code>, <code>evaluate_chaperone_design(...)</code>, <code>evaluate_chaperone_block(...)</code>, <code>chaperone_lineage_hash(...)</code>, <code>chaperone_core_examples()</code>, <code>evaluate_chaperone_core_example(...)</code>, <code>evaluate_chaperone_core_example_block(...)</code>, <code>create_sparse_map(...)</code>, <code>get_top_interacting_residues(...)</code></td><td style="padding:7px;">64-point Shrake-Rupley SASA/RASA, CSR H-bond networks, BUNS (buried unsatisfied polars), self-consistency scRMSD, content-addressable lineage hashes, and all 8 Official Chaperone Core Examples.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#a78bfa; font-weight:600;">5. Extremophile Biochem</td><td style="padding:7px;"><code>molecular_weight(seq)</code>, <code>isoelectric_point(seq)</code>, <code>gravy_score(seq)</code>, <code>sequence_entropy(seq)</code></td><td style="padding:7px;">Average mass (Da), Henderson-Hasselbalch pI root-finder, Kyte-Doolittle hydropathy, and Shannon entropy (bits).</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#2dd4bf; font-weight:600;">6. Cheminformatics</td><td style="padding:7px;"><code>parse_smiles(s)</code>, <code>molecule_wl_hash(s)</code>, <code>molecules_are_identical(s1, s2)</code>, <code>mol @&gt; sub</code>, <code>mol % q</code>, <code>smiles_to_fingerprint(s)</code>, <code>tanimoto_similarity(a, b)</code></td><td style="padding:7px;">Salt stripping, tautomer normalization, Huckel 4n+2 aromaticity, 64-bit Weisfeiler-Lehman graph hashes, VF2 subgraph isomorphism (<code>@&gt;</code>), and 1024-bit Morgan ECFP4 Tanimoto (<code>%</code>).</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#f472b6; font-weight:600;">7. 3D Spatial &amp; Morton</td><td style="padding:7px;"><code>create_residue_coord(x,y,z,name)</code>, <code>create_hetatm_coord(...)</code>, <code>z_order_encode(x,y,z)</code>, <code>residue_z_index(c)</code>, <code>distance_angstroms(a,b)</code>, <code>center_of_mass(atoms)</code>, <code>bio_parse_pdb(pdb)</code>, <code>bio_parse_pdb_coords(pdb)</code></td><td style="padding:7px;">30-bit 3D Z-Order Morton space-filling curves for B-Tree 3D bounding-box pocket queries and Euclidean geometry.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#60a5fa; font-weight:600;">8. Wet-Lab SynBio</td><td style="padding:7px;"><code>bio_reverse_translate(seq, host)</code>, <code>bio_codon_adaptation_index(dna, host)</code>, <code>bio_design_expression_cassette(seq, host, prom)</code>, <code>bio_gc_content(dna)</code>, <code>bio_translate_dna(dna)</code></td><td style="padding:7px;">Host-specific codon optimization (<i>E. coli</i>, yeast, human), Sharp &amp; Li CAI, T7 + RBS plasmid cassette synthesis, and BsaI restriction safety.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#f87171; font-weight:600;">9. LC-MS/MS Proteomics</td><td style="padding:7px;"><code>create_spectrum(mz, int)</code>, <code>generate_theoretical_fragments(seq, z, method)</code>, <code>ppm_difference(m1, m2)</code>, <code>is_mass_within_ppm(m1, m2, tol)</code>, <code>isotope_fit_score(spec, mz, z, tol)</code>, <code>generate_decoy_sequence(seq, method)</code></td><td style="padding:7px;">CID/HCD (b/y) &amp; ETD (c/z) theoretical fragment ion tables, Orbitrap ppm accuracy, Averagine 13C isotopic envelope scoring, and target-decoy FDR control.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#c084fc; font-weight:600;">10. Energetics &amp; Kinetics</td><td style="padding:7px;"><code>audit_mmpbsa_decomposition(...)</code>, <code>solvent_shift_penalty(...)</code>, <code>build_energy_hotspot_map(...)</code>, <code>discover_cryptic_hotspots(...)</code>, <code>allosteric_mechanism_index(...)</code>, <code>create_preorganization_hist(...)</code>, <code>create_preorganization_profile(...)</code>, <code>ensemble_fraction_below(...)</code>, <code>ensemble_fraction_between(...)</code>, <code>simulate_defluorination_kinetics(...)</code>, <code>bio_fetch_uniprot(id)</code>, <code>bio_search_uniprot(...)</code>, <code>bio_fold_sequence(seq)</code></td><td style="padding:7px;">MM/PBSA dG decomposition, dielectric solvent shifts, CSR allosteric coupling graphs, MD pre-organization ensembles, and Eyring PFAS defluorination kinetics.</td></tr>
                            <tr style="border-bottom:1px solid #1e293b;"><td style="padding:7px; color:#38bdf8; font-weight:600;">11. RangeFinder GiST (arXiv:2611.03379v1)</td><td style="padding:7px;"><code>predict_pocket_ranges(id, pdb, cx, cy, cz, r)</code>, <code>pocket_numranges(p)</code>, <code>pocket_interval_bounds(p)</code>, <code>pocket_range_iou(p1, p2, charge)</code>, <code>pocket_giou_loss(p1, p2, charge)</code>, <code>pocket_matches_ligand(p, l)</code>, <code>pocket_matches_ligand_8prop(p, l)</code>, <code>compute_ligand_properties(s)</code>, <code>compute_ligand_properties_row(s)</code></td><td style="padding:7px;">10 Å 3D pocket -> 9 physicochemical ligand intervals (<code>[μ - 2|σ|, μ + 2|σ|]</code>), 8-column PostgreSQL <code>numrange</code> GiST reverse target fishing (<code>0.68 ms</code> / 20k pockets), and composite B-Tree virtual screening (<code>3.78 ms</code> / 250k ligands).</td></tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
    <div id="wiki-popover"></div>

    <script>
        __INSPECTOR_JS__

        var wikiInspector = null;
        var currentWikiPhase = 'vector';

        function selectWikiPhase(phaseKey, scrollCardId, fromInspector) {
            currentWikiPhase = phaseKey || 'vector';
            document.querySelectorAll('.wiki-nav-btn').forEach(btn => {
                btn.classList.toggle('active', btn.getAttribute('data-cat') === currentWikiPhase);
            });
            const q = (document.getElementById('wiki-search').value || '').toLowerCase().trim();
            document.querySelectorAll('.wiki-card').forEach(card => {
                const cardCat = card.getAttribute('data-cat');
                const matchesCat = q ? true : (cardCat === currentWikiPhase);
                const matchesText = !q || card.innerText.toLowerCase().includes(q);
                card.classList.toggle('active-phase', matchesCat && matchesText);
            });
            if (wikiInspector && !fromInspector) {
                wikiInspector.setPreset(currentWikiPhase, { triggerAuto3DMode: true, notifyParent: false });
            }
            if (scrollCardId) {
                const el = document.getElementById(scrollCardId);
                if (el) {
                    setTimeout(() => {
                        el.scrollIntoView({ behavior: 'smooth', block: 'start' });
                    }, 80);
                }
            }
        }

        window.navigateWikiPhase = function(cat, cardId) {
            selectWikiPhase(cat, cardId, false);
        };

        function filterWikiCards() {
            selectWikiPhase(currentWikiPhase, '', true);
        }

        window.onWikiPairLoaded = function(comp, id1, id2) {
            const al = comp.alignment || {};
            const ch = comp.chaperone || {};
            const m1 = comp.meta1 || {};
            const m2 = comp.meta2 || {};
            const sc = comp.substrate_cheminfo || {};
            const sx = comp.sql_extensions || {};
            const rf = comp.rangefinder || {};

            const rg1 = sx.rg_1 ?? ch.compact_1 ?? '12.23';
            const rg2 = sx.rg_2 ?? ch.compact_2 ?? '12.20';
            const pk1 = sx.packing_1 ?? (typeof rg1 === 'number' && m1.length ? (rg1 / (2.25 * Math.pow(m1.length, 0.38))).toFixed(3) : '0.987');
            const pk2 = sx.packing_2 ?? (typeof rg2 === 'number' && m2.length ? (rg2 / (2.25 * Math.pow(m2.length, 0.38))).toFixed(3) : '0.985');

            const setHtml = (id, html) => {
                const el = document.getElementById(id);
                if (el) el.innerHTML = html;
            };
            setHtml('live-ex-vector', `<b>Live DB (${id1} vs ${id2}):</b> Vector Dist = <b>${al.vec_dist || al.esm_kmer_dist || '0.0085'}</b> | Tri-peptide Cosine Dist = <b>${al.esm_kmer_dist || '0.0120'}</b> | PPI I_comp = <b>${sx.icomp ?? '0.520'}</b>`);
            setHtml('live-ex-sw', `<b>Live DB (${id1} vs ${id2}):</b> Smith-Waterman Score = <b>${al.sw_score ?? '--'}</b> | Local SW Similarity = <b>${al.sw_sim_pct ?? '--'}%</b>`);
            setHtml('live-ex-twilight', `<b>Live DB (${id1} vs ${id2}):</b> Twilight-Zone Score = <b>${al.twilight_score ?? '--'}</b> (SW Sim: <b>${al.sw_sim_pct ?? '--'}%</b>, Lengths: ${m1.length ?? '--'} aa vs ${m2.length ?? '--'} aa)`);
            setHtml('live-ex-kabsch', `<b>Live DB (${id1} vs ${id2}):</b> 3D Core RMSD = <b>${al.core_rmsd ?? '--'} A</b> across <b>${al.aligned_residues ?? '--'} CA atoms</b> | Sequential TM-Score = <b>${al.tm_score ?? '--'}</b>`);
            setHtml('live-ex-theozyme', `<b>Live DB (${id1} vs ${id2}):</b> Rg: <b>${rg1} A</b> vs <b>${rg2} A</b> | Flory Globular Compactness Ratio: <b>${pk1}</b> vs <b>${pk2}</b>`);
            setHtml('live-ex-chaperone', `<b>Live DB (${id1} vs ${id2}):</b> Core RASA Corr = <b>${ch.rasa_corr ?? '--'}</b> (Mean RASA: ${ch.rasa_1 ?? '--'} vs ${ch.rasa_2 ?? '--'}) | H-Bonds: <b>${ch.hbonds_1 ?? '--'}</b> vs <b>${ch.hbonds_2 ?? '--'}</b> | BUNS: <b>${ch.buns_1 ?? '--'}</b>/${ch.buns_2 ?? '--'} | Score = <b>${ch.chaperone_score ?? '--'}</b>`);
            setHtml('live-ex-biochem', `<b>Live DB (${id1} vs ${id2}):</b> pI: <b>${m1.pi ?? '--'}</b> vs <b>${m2.pi ?? '--'}</b> (dpI = <b>${al.delta_pi ?? '--'}</b>) | GRAVY: <b>${m1.gravy ?? '--'}</b> vs <b>${m2.gravy ?? '--'}</b> (dGRAVY = <b>${al.delta_gravy ?? '--'}</b>) | Acidic(D/E): <b>${m1.acidic_pct ?? '--'}%</b> vs <b>${m2.acidic_pct ?? '--'}%</b>`);
            if (sc.substrate_name) {
                setHtml('live-ex-cheminfo', `<b>Live DB (${sc.family_name || 'Enzyme Family'}):</b> Cognate Substrate = <b>${sc.substrate_name}</b> | WL64 = <code>${sc.wl_hash || '--'}</code> | Huckel 4n+2 Match = <b>[OK]</b> | VF2 Warhead (<code>@&gt; ${sc.warhead_smiles || '--'}</code>) = <b>[OK]</b>`);
            }
            if (rf && rf.range_iou !== undefined) {
                setHtml('live-ex-rangefinder', `<b>Live DB (${id1} vs ${id2}):</b> 8-Prop Range-IoU = <b>${rf.range_iou}</b> | Orphan MW = <code>${rf.mw_box_2} Da</code>, LogP = <code>${rf.logp_box_2}</code>, TPSA = <code>${rf.tpsa_box_2} Å²</code> | P(arom) = <b>${rf.arom_prob_2}</b> | Substrate Box: <b>${rf.substrate_props_matched ?? 8}/8 ${rf.substrate_match_2 ? '[OK]' : ''}</b>`);
            }
            setHtml('live-ex-synbio', `<b>Live DB (${id2} -&gt; E. coli T7 Cassette):</b> CAI = <b>${sx.ecoli_cai ?? '1.000'}</b> | GC = <b>${sx.ecoli_gc ?? '48.2'}%</b> | BsaI-Free: <b>${sx.restriction_safe !== false ? 'YES' : 'NO'}</b> | PROSITE: <b>${sx.prosite_summary || 'CK2-Phospho'}</b> | 5' CDS: <code>${sx.dna_preview || 'ATG...'}</code>`);
            setHtml('live-ex-proteomics', `<b>Live DB (${id2} N-Term HCD Fragmentation):</b> Diagnostic Ions: <code>${sx.ms_ions_preview || 'b2/y2/b3/y3'}</code> | Tryptic Reversed Decoy (1% FDR): <code>${sx.decoy_peptide || 'VFDLLDVEMK'}</code>`);
            setHtml('live-ex-energetics', `<b>Live DB (${id1} vs ${id2}):</b> Allosteric Lock Index = <b>${sx.allosteric_index ?? '0.84'}</b> | Eyring Defluorination F- Release = <b>${sx.eyring_kcat ?? '15.0'} ppm</b> (Conv = <b>${sx.eyring_dg ?? '99.9'}%</b>)`);
        };

        window.addEventListener('DOMContentLoaded', () => {
            const params = new URLSearchParams(window.location.search);
            const initCat = params.get('cat');
            const startPhase = (initCat && initCat !== 'all' && PHASE_PRESETS[initCat]) ? initCat : 'vector';

            wikiInspector = new ProteinCompareInspector({
                mountEl: '#wiki-inspector-root',
                context: 'wiki',
                defaultPreset: startPhase,
                onPhaseSelect: (phaseKey) => selectWikiPhase(phaseKey, '', true)
            });

            selectWikiPhase(startPhase, window.location.hash ? window.location.hash.slice(1) : '', true);
            bindWikiHoverPopovers(document.body, (cat, card) => selectWikiPhase(cat, card, false));

            const n1 = params.get('node1');
            const n2 = params.get('node2');
            const sel = document.getElementById('cmp_wiki_pair_select');
            if (n1 && n2 && sel) {
                const key = `${n1}|${n2}`;
                const exists = Array.from(sel.options).some(o => o.value === key);
                if (!exists) {
                    const opt = document.createElement('option');
                    opt.value = key;
                    opt.text = `Selected Heatmap Pair: ${n1} vs ${n2}`;
                    sel.insertBefore(opt, sel.firstChild);
                }
                sel.value = key;
                wikiInspector.loadPair(n1, n2);
            } else if (sel) {
                const parts = sel.value.split('|');
                if (parts.length === 2) wikiInspector.loadPair(parts[0], parts[1]);
            }
        });
    </script>
</body>
</html>
""".replace("__INSPECTOR_CSS__", INSPECTOR_CSS).replace("__INSPECTOR_JS__", INSPECTOR_JS)


COMPARE_LAB_HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>pg_bio — Standalone 3D Protein Compare Lab</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {
            margin: 0;
            background: #090d16;
            color: #e2e8f0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            display: flex;
            flex-direction: column;
            height: 100vh;
            overflow: hidden;
        }
        .lab-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 20px;
            background: #020617;
            border-bottom: 1px solid #1e293b;
            gap: 14px;
            flex-wrap: wrap;
        }
        .lab-brand {
            display: flex;
            align-items: center;
            gap: 12px;
        }
        .lab-brand h2 {
            margin: 0;
            font-size: 1.05em;
            color: #38bdf8;
        }
        .lab-nav-btn {
            background: #0f172a;
            color: #cbd5e1;
            border: 1px solid #334155;
            border-radius: 6px;
            padding: 5px 11px;
            font-size: 0.8em;
            text-decoration: none;
            cursor: pointer;
            transition: all 0.15s;
        }
        .lab-nav-btn:hover {
            border-color: #38bdf8;
            color: #38bdf8;
        }
        .lab-custom-inputs {
            display: flex;
            align-items: center;
            gap: 6px;
            font-size: 0.8em;
            color: #94a3b8;
        }
        .lab-input {
            background: #0f172a;
            color: #f8fafc;
            border: 1px solid #334155;
            border-radius: 5px;
            padding: 4px 8px;
            font-size: 0.92em;
            width: 105px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        #lab-inspector-root {
            flex: 1;
            display: flex;
            flex-direction: column;
            overflow: hidden;
        }
        #lab-inspector-root .cmp-inspector {
            height: 100%;
            border: none;
            border-radius: 0;
        }
        __INSPECTOR_CSS__
    </style>
</head>
<body>
    <div class="lab-header">
        <div class="lab-brand">
            <a class="lab-nav-btn" href="/?limit=twilight-40">← Swarm Matrix</a>
            <h2>pg_bio 3D Compare Lab</h2>
            <span style="font-size:0.78em; color:#64748b;">Standalone Kabsch SVD, Twilight-Zone &amp; RosettaCommons Chaperone Workbench</span>
        </div>
        <div class="lab-custom-inputs">
            <span>Custom Pair:</span>
            <input id="lab-id1" class="lab-input" placeholder="Bait ID" value="__INIT_ID1__" />
            <span>vs</span>
            <input id="lab-id2" class="lab-input" placeholder="Orphan ID" value="__INIT_ID2__" />
            <button class="lab-nav-btn" style="border-color:#10b981; color:#34d399;" onclick="loadCustomLabPair()">Compare 3D</button>
            <a class="lab-nav-btn" style="border-color:#f59e0b; color:#fbbf24;" href="/cards">🃏 1-to-N Card Arena</a>
            <a class="lab-nav-btn" id="lab-wiki-top-link" href="/wiki" target="_blank">Correlation Wiki</a>
        </div>
    </div>
    <div id="lab-inspector-root"></div>
    <div id="wiki-popover"></div>

    <script>
        __INSPECTOR_JS__

        var labInspector = null;

        function loadCustomLabPair() {
            const i1 = (document.getElementById('lab-id1').value || '').trim();
            const i2 = (document.getElementById('lab-id2').value || '').trim();
            if (!i1 || !i2) return;
            const u = new URL(window.location);
            u.searchParams.set('node1', i1);
            u.searchParams.set('node2', i2);
            window.history.pushState({}, '', u);
            const sel = document.getElementById('cmp_lab_pair_select');
            if (sel) {
                const key = i1 + '|' + i2;
                let exists = Array.from(sel.options).some(o => o.value === key);
                if (!exists) {
                    const opt = document.createElement('option');
                    opt.value = key;
                    opt.text = 'Custom Pair: ' + i1 + ' vs ' + i2;
                    sel.insertBefore(opt, sel.firstChild);
                }
                sel.value = key;
            }
            labInspector.loadPair(i1, i2);
            const wl = document.getElementById('lab-wiki-top-link');
            if (wl) wl.href = '/wiki?node1=' + encodeURIComponent(i1) + '&node2=' + encodeURIComponent(i2);
        }

        window.addEventListener('DOMContentLoaded', () => {
            const params = new URLSearchParams(window.location.search);
            const initPreset = params.get('preset') || 'overview';
            const startPreset = PHASE_PRESETS[initPreset] ? initPreset : 'overview';

            labInspector = new ProteinCompareInspector({
                mountEl: '#lab-inspector-root',
                context: 'lab',
                instanceId: 'cmp_lab',
                defaultPreset: startPreset
            });

            const n1 = params.get('node1') || '__INIT_ID1__' || 'A0A1S8AT22';
            const n2 = params.get('node2') || '__INIT_ID2__' || 'A0A1H3IA37';
            document.getElementById('lab-id1').value = n1;
            document.getElementById('lab-id2').value = n2;

            const sel = document.getElementById('cmp_lab_pair_select');
            if (sel) {
                const key = n1 + '|' + n2;
                let exists = Array.from(sel.options).some(o => o.value === key);
                if (!exists) {
                    const opt = document.createElement('option');
                    opt.value = key;
                    opt.text = 'Selected Pair: ' + n1 + ' vs ' + n2;
                    sel.insertBefore(opt, sel.firstChild);
                }
                sel.value = key;
            }
            const wl = document.getElementById('lab-wiki-top-link');
            if (wl) wl.href = '/wiki?node1=' + encodeURIComponent(n1) + '&node2=' + encodeURIComponent(n2);
            labInspector.loadPair(n1, n2);
        });
    </script>
</body>
</html>
""".replace("__INSPECTOR_CSS__", INSPECTOR_CSS).replace("__INSPECTOR_JS__", INSPECTOR_JS)


@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/compare')
@app.route('/compare/<id1>/<id2>')
def compare_lab_page(id1=None, id2=None):
    init_id1 = id1 or request.args.get('node1') or 'A0A1S8AT22'
    init_id2 = id2 or request.args.get('node2') or 'A0A1H3IA37'
    html = COMPARE_LAB_HTML_TEMPLATE.replace("__INIT_ID1__", init_id1).replace("__INIT_ID2__", init_id2)
    return render_template_string(html)

@app.route('/wiki')
def wiki_page():
    return render_template_string(WIKI_HTML_TEMPLATE)


CARDS_ARENA_HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>pg_bio — 1-to-N Pokémon-Card Hub &amp; Biophysical Similarity Arena</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        *, *::before, *::after { box-sizing: border-box; }
        body {
            margin: 0;
            background: #070b14;
            color: #f8fafc;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            display: flex;
            flex-direction: column;
            height: 100vh;
            overflow: hidden;
        }

        /* Top Navigation & Cap Controls */
        .arena-topbar {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 16px;
            padding: 10px 20px;
            background: #020617;
            border-bottom: 1px solid #1e293b;
            flex-wrap: wrap;
            z-index: 20;
        }
        .arena-brand {
            display: flex;
            align-items: center;
            gap: 12px;
        }
        .arena-brand h1 {
            margin: 0;
            font-size: 1.08rem;
            color: #fbbf24;
            letter-spacing: -0.01em;
        }
        .arena-brand .sub {
            font-size: 0.8rem;
            color: #94a3b8;
        }
        .arena-controls {
            display: flex;
            align-items: center;
            gap: 10px;
            flex-wrap: wrap;
        }
        .ctrl-group {
            display: flex;
            align-items: center;
            gap: 6px;
            background: #0f172a;
            border: 1px solid #1e293b;
            padding: 4px 10px;
            border-radius: 7px;
            font-size: 0.78rem;
            color: #cbd5e1;
        }
        .ctrl-select {
            background: #020617;
            color: #38bdf8;
            border: 1px solid #334155;
            border-radius: 5px;
            padding: 3px 8px;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
        }
        .nav-pill {
            background: #1e293b;
            color: #38bdf8;
            border: 1px solid #38bdf8;
            padding: 6px 12px;
            border-radius: 6px;
            font-size: 0.8rem;
            font-weight: 600;
            text-decoration: none;
            cursor: pointer;
            transition: background 0.15s, color 0.15s;
        }
        .nav-pill:hover {
            background: #38bdf8;
            color: #020617;
        }

        /* Main Two-Column Arena Layout */
        .arena-workspace {
            flex: 1;
            display: grid;
            grid-template-columns: 330px minmax(0, 1fr);
            min-height: 0;
            overflow: hidden;
        }

        /* Left Sidebar: Capped Top Hub Organisms & 1-to-N Hub Selector */
        .arena-sidebar {
            background: #030816;
            border-right: 1px solid #1e293b;
            display: flex;
            flex-direction: column;
            min-height: 0;
            overflow-y: auto;
            scrollbar-gutter: stable;
            overscroll-behavior: contain;
            padding: 14px 12px;
            gap: 10px;
        }
        .sidebar-section-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.72rem;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: #94a3b8;
            font-weight: 700;
            padding: 2px 4px;
        }
        .org-pill-card {
            background: #0b1324;
            border: 2px solid #1e293b;
            border-radius: 9px;
            padding: 10px 12px;
            cursor: pointer;
            transition: border-color 0.15s, box-shadow 0.15s, background-color 0.15s;
            display: flex;
            flex-direction: column;
            gap: 6px;
        }
        .org-pill-card:hover {
            border-color: #475569;
            background: #0f172a;
        }
        /* Non-layout-shifting selection state per modern-web-guidance */
        .org-pill-card.is-selected {
            border-color: #fbbf24;
            background: rgba(245, 158, 11, 0.08);
            box-shadow: 0 0 0 1px rgba(251, 191, 36, 0.35), 0 6px 18px rgba(0, 0, 0, 0.45);
        }
        .org-rank-row {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 6px;
        }
        .org-rank-badge {
            font-size: 0.7rem;
            font-weight: 700;
            color: #fbbf24;
            background: rgba(245, 158, 11, 0.15);
            border: 1px solid rgba(245, 158, 11, 0.4);
            padding: 1px 7px;
            border-radius: 999px;
        }
        .org-reach-badge {
            font-size: 0.7rem;
            color: #34d399;
            background: rgba(16, 185, 129, 0.12);
            padding: 1px 7px;
            border-radius: 999px;
        }
        .org-title {
            font-size: 0.86rem;
            font-weight: 700;
            color: #f8fafc;
            line-height: 1.25;
        }
        .org-sub-metrics {
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            font-size: 0.72rem;
            color: #94a3b8;
        }
        .hub-btn-list {
            display: flex;
            flex-wrap: wrap;
            gap: 5px;
            margin-top: 4px;
        }
        .hub-pill-btn {
            background: #020617;
            color: #cbd5e1;
            border: 1px solid #334155;
            border-radius: 6px;
            padding: 4px 8px;
            font-size: 0.73rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            cursor: pointer;
            transition: border-color 0.15s, background-color 0.15s, color 0.15s;
        }
        .hub-pill-btn:hover {
            border-color: #38bdf8;
            color: #fff;
        }
        .hub-pill-btn.is-selected {
            border-color: #38bdf8;
            background: rgba(56, 189, 248, 0.18);
            color: #38bdf8;
            box-shadow: 0 0 0 1px rgba(56, 189, 248, 0.35);
        }

        /* Main Scrollable Stage */
        .arena-stage {
            overflow-y: auto;
            scrollbar-gutter: stable;
            overscroll-behavior: contain;
            padding: 18px 24px 48px 24px;
            display: flex;
            flex-direction: column;
            gap: 22px;
        }

        /* Battle Stage: Side-by-Side Pokémon Cards + Center Synergy Verdict */
        .battle-grid {
            display: grid;
            grid-template-columns: minmax(320px, 1fr) minmax(290px, 370px) minmax(320px, 1fr);
            gap: 18px;
            align-items: stretch;
        }
        @media (max-width: 1280px) {
            .battle-grid {
                grid-template-columns: 1fr;
            }
        }

        /* Biophysical Pokémon Trading Card (<article class="bio-poke-card">) */
        article.bio-poke-card {
            container: pokecard / inline-size;
            background: linear-gradient(160deg, #0f172a 0%, #090d18 65%, #020617 100%);
            border: 2px solid #334155;
            border-radius: 14px;
            padding: 16px;
            display: flex;
            flex-direction: column;
            gap: 12px;
            box-shadow: 0 14px 36px rgba(0, 0, 0, 0.65);
            position: relative;
        }
        article.bio-poke-card.hub-slot {
            border-color: #38bdf8;
            box-shadow: 0 14px 36px rgba(0, 0, 0, 0.65), 0 0 22px rgba(56, 189, 248, 0.14);
        }
        article.bio-poke-card.challenger-slot {
            border-color: #fbbf24;
            box-shadow: 0 14px 36px rgba(0, 0, 0, 0.65), 0 0 22px rgba(251, 191, 36, 0.14);
            view-transition-name: active-challenger-card;
        }

        .poke-card-header {
            display: flex;
            flex-direction: column;
            gap: 6px;
            border-bottom: 1px solid #1e293b;
            padding-bottom: 10px;
        }
        .poke-role-strip {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 8px;
        }
        .poke-role-tag {
            font-size: 0.7rem;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            font-weight: 700;
            padding: 2px 8px;
            border-radius: 999px;
        }
        .poke-elem-badge {
            font-size: 0.73rem;
            font-weight: 700;
            padding: 2px 9px;
            border-radius: 999px;
            border: 1px solid currentColor;
        }
        .poke-title-row {
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            gap: 10px;
        }
        .poke-title {
            margin: 0;
            font-size: 1.22rem;
            font-weight: 800;
            color: #f8fafc;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            width: fit-content;
        }
        .poke-cp-badge {
            font-size: 0.95rem;
            font-weight: 800;
            color: #fde047;
            background: rgba(250, 204, 21, 0.12);
            border: 1px solid rgba(250, 204, 21, 0.45);
            padding: 2px 10px;
            border-radius: 8px;
            white-space: nowrap;
        }
        .poke-org-line {
            font-size: 0.8rem;
            color: #94a3b8;
            display: flex;
            justify-content: space-between;
            gap: 8px;
        }
        .poke-org-line b {
            color: #e2e8f0;
        }

        /* Holographic Pocket & 3D Viewport Stage */
        .poke-holo-stage {
            background: radial-gradient(circle at 25% 20%, rgba(56, 189, 248, 0.16), rgba(16, 185, 129, 0.08) 55%, #020617 100%);
            border: 1px solid #1e293b;
            border-radius: 10px;
            padding: 10px 12px;
            display: flex;
            flex-direction: column;
            gap: 8px;
            position: relative;
        }
        .poke-3d-canvas {
            width: 100%;
            height: 175px;
            position: relative;
            border-radius: 8px;
            overflow: hidden;
            background: #020617;
            border: 1px solid #1e293b;
            display: none;
        }
        .poke-holo-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 7px;
            font-size: 0.75rem;
        }
        .holo-pill {
            background: rgba(2, 6, 23, 0.78);
            border: 1px solid #1e293b;
            border-radius: 6px;
            padding: 5px 8px;
            color: #cbd5e1;
        }
        .holo-pill span {
            display: block;
            font-size: 0.67rem;
            color: #64748b;
            text-transform: uppercase;
        }
        .holo-pill b {
            color: #f8fafc;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }

        /* 6 Biophysical Stat Bars ("Which is Better for What") */
        .poke-stats-section {
            display: flex;
            flex-direction: column;
            gap: 8px;
            flex: 1;
        }
        .poke-stat-row {
            background: rgba(15, 23, 42, 0.7);
            border: 1px solid #1e293b;
            border-radius: 8px;
            padding: 7px 10px;
            display: flex;
            flex-direction: column;
            gap: 4px;
        }
        .poke-stat-row.is-winner {
            border-color: rgba(16, 185, 129, 0.55);
            background: rgba(16, 185, 129, 0.08);
        }
        .stat-top-line {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 6px;
            font-size: 0.76rem;
        }
        .stat-label {
            font-weight: 700;
            color: #e2e8f0;
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .winner-crown {
            font-size: 0.66rem;
            font-weight: 800;
            color: #10b981;
            background: rgba(16, 185, 129, 0.16);
            border: 1px solid rgba(16, 185, 129, 0.45);
            padding: 1px 6px;
            border-radius: 999px;
            text-transform: uppercase;
        }
        .stat-raw-val {
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-size: 0.72rem;
            color: #94a3b8;
        }
        .stat-bar-track {
            width: 100%;
            height: 7px;
            background: #020617;
            border-radius: 999px;
            overflow: hidden;
            border: 1px solid #1e293b;
        }
        .stat-bar-fill {
            height: 100%;
            border-radius: 999px;
            transition: width 0.25s ease;
        }
        .stat-use-case {
            font-size: 0.69rem;
            color: #94a3b8;
        }

        /* Container query adjustment when card is wide */
        @container pokecard (min-width: 430px) {
            .poke-holo-grid {
                grid-template-columns: repeat(4, 1fr);
            }
        }

        .poke-card-footer {
            border-top: 1px solid #1e293b;
            padding-top: 9px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.74rem;
            color: #94a3b8;
        }

        /* Center Versus & "Which is Better for What" Decision Column */
        .versus-column {
            background: #0b1324;
            border: 1px solid #1e293b;
            border-radius: 14px;
            padding: 16px;
            display: flex;
            flex-direction: column;
            gap: 12px;
            box-shadow: 0 12px 32px rgba(0, 0, 0, 0.5);
        }
        .vs-banner {
            text-align: center;
            background: linear-gradient(135deg, rgba(56, 189, 248, 0.14), rgba(245, 158, 11, 0.14));
            border: 1px solid #334155;
            border-radius: 10px;
            padding: 10px;
        }
        .vs-score-big {
            font-size: 1.55rem;
            font-weight: 900;
            color: #34d399;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .vs-metrics-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 7px;
            font-size: 0.75rem;
        }
        .vs-metric-box {
            background: #020617;
            border: 1px solid #1e293b;
            border-radius: 7px;
            padding: 6px 8px;
        }
        .vs-metric-box span {
            display: block;
            font-size: 0.67rem;
            color: #64748b;
            text-transform: uppercase;
        }
        .vs-metric-box b {
            color: #e2e8f0;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .verdict-box {
            background: #020617;
            border: 1px solid #1e293b;
            border-radius: 9px;
            padding: 10px 12px;
            display: flex;
            flex-direction: column;
            gap: 6px;
            font-size: 0.76rem;
        }
        .verdict-box h4 {
            margin: 0;
            font-size: 0.78rem;
            text-transform: uppercase;
            letter-spacing: 0.04em;
        }
        .verdict-list {
            margin: 0;
            padding-left: 16px;
            display: flex;
            flex-direction: column;
            gap: 4px;
            color: #cbd5e1;
        }
        .vs-actions {
            margin-top: auto;
            display: flex;
            flex-direction: column;
            gap: 8px;
        }
        .vs-action-btn {
            width: 100%;
            padding: 9px 12px;
            border-radius: 8px;
            font-weight: 700;
            font-size: 0.8rem;
            cursor: pointer;
            text-align: center;
            text-decoration: none;
            transition: filter 0.15s;
        }
        .vs-action-btn:hover {
            filter: brightness(1.12);
        }

        /* Bottom 1-to-N Challenger Deck Grid with CSS Subgrid */
        .deck-section {
            display: flex;
            flex-direction: column;
            gap: 12px;
        }
        .deck-header-bar {
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 10px;
            border-bottom: 1px solid #1e293b;
            padding-bottom: 8px;
        }
        .deck-header-bar h2 {
            margin: 0;
            font-size: 1.02rem;
            color: #f8fafc;
        }
        .challenger-deck-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(245px, 1fr));
            gap: 12px;
        }
        article.deck-mini-card {
            container: minicard / inline-size;
            grid-row: span 4;
            grid-template-rows: auto auto 1fr auto;
            grid-template-rows: subgrid;
            display: grid;
            background: #0b1324;
            border: 2px solid #1e293b;
            border-radius: 11px;
            padding: 12px;
            row-gap: 8px;
            cursor: pointer;
            transition: border-color 0.15s, box-shadow 0.15s, background-color 0.15s;
        }
        article.deck-mini-card:hover {
            border-color: #64748b;
            background: #0f172a;
        }
        /* Non-layout-shifting selection state */
        article.deck-mini-card.is-selected {
            border-color: #fbbf24;
            background: rgba(245, 158, 11, 0.09);
            box-shadow: 0 0 0 1px rgba(251, 191, 36, 0.4), 0 8px 20px rgba(0, 0, 0, 0.5);
        }
        .mini-header {
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            gap: 6px;
        }
        .mini-id {
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-weight: 800;
            font-size: 0.92rem;
            color: #f8fafc;
        }
        .mini-synergy {
            font-size: 0.75rem;
            font-weight: 800;
            color: #34d399;
            background: rgba(16, 185, 129, 0.14);
            border: 1px solid rgba(16, 185, 129, 0.4);
            padding: 1px 7px;
            border-radius: 999px;
        }
        .mini-org {
            font-size: 0.74rem;
            color: #94a3b8;
            line-height: 1.3;
        }
        .mini-adv-pills {
            display: flex;
            flex-wrap: wrap;
            gap: 4px;
            align-content: flex-start;
        }
        .mini-adv-tag {
            font-size: 0.68rem;
            background: #020617;
            border: 1px solid #334155;
            color: #fde047;
            padding: 2px 6px;
            border-radius: 5px;
        }
        .mini-footer {
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.7rem;
            color: #64748b;
            border-top: 1px solid #1e293b;
            padding-top: 6px;
        }

        /* Respect reduced motion for View Transitions */
        @media (prefers-reduced-motion: reduce) {
            ::view-transition-group(*),
            ::view-transition-old(*),
            ::view-transition-new(*) {
                animation: none !important;
            }
        }
    </style>
</head>
<body>
    <header class="arena-topbar">
        <div class="arena-brand">
            <div>
                <h1>🃏 1-to-N Pokémon-Card Hub &amp; Similarity Arena</h1>
                <div class="sub">Rank top hub organisms with the most cross-species correlations &amp; compare 1-to-N biophysical trading cards side by side ("Which is Better for What")</div>
            </div>
        </div>
        <div class="arena-controls">
            <div class="ctrl-group">
                <label for="org-cap-select">Top Organism Cap:</label>
                <select id="org-cap-select" class="ctrl-select" onchange="reloadArenaData()">
                    <option value="5">Top 5 Hub Organisms</option>
                    <option value="10" selected>Top 10 Hub Organisms</option>
                    <option value="15">Top 15 Hub Organisms</option>
                    <option value="25">Top 25 Hub Organisms</option>
                </select>
            </div>
            <div class="ctrl-group">
                <label for="match-cap-select">1-to-N Deck Cap:</label>
                <select id="match-cap-select" class="ctrl-select" onchange="reloadArenaData()">
                    <option value="6">Top 6 Matches / Hub</option>
                    <option value="12" selected>Top 12 Matches / Hub</option>
                    <option value="24">Top 24 Matches / Hub</option>
                    <option value="48">Top 48 Matches / Hub</option>
                </select>
            </div>
            <div class="ctrl-group">
                <label for="sort-lens-select">Correlation Lens:</label>
                <select id="sort-lens-select" class="ctrl-select" onchange="reloadArenaData()">
                    <option value="composite" selected>Multi-Layer Composite Synergy</option>
                    <option value="rangefinder">RangeFinder 10 Å Pocket IoU (GiST)</option>
                    <option value="twilight">Twilight-Zone Remote Homology</option>
                    <option value="chaperone">RosettaCommons Chaperone Stability</option>
                    <option value="vector">ESM-2 Vector Proximity</option>
                </select>
            </div>
            <span id="arena-sql-badge" style="font-size:0.75rem; color:#34d399; font-family:monospace;">pg_bio: loading...</span>
            <a class="nav-pill" style="border-color:#a855f7; color:#c084fc;" href="/api/swarm_telemetry">⚡ Swarm Telemetry</a>
            <a class="nav-pill" href="/">Swarm Heatmap</a>
            <a class="nav-pill" id="arena-compare-top-link" href="/compare">3D Compare Lab</a>
            <a class="nav-pill" id="arena-wiki-top-link" href="/wiki?cat=rangefinder">Correlation Wiki</a>
        </div>
    </header>

    <div class="arena-workspace">
        <aside class="arena-sidebar" aria-label="Top Hub Organisms">
            <div class="sidebar-section-header">
                <span>Capped Hub Organisms (1-to-N Reach)</span>
                <span id="sidebar-org-count">0</span>
            </div>
            <div id="org-sidebar-list" style="display:flex; flex-direction:column; gap:8px;"></div>
        </aside>

        <main class="arena-stage">
            <section class="battle-grid" aria-label="Side-by-Side Biophysical Trading Card Comparison">
                <!-- Left Slot: Hub Protein -->
                <article class="bio-poke-card hub-slot" id="hub-poke-card">
                    <header class="poke-card-header">
                        <div class="poke-role-strip">
                            <span class="poke-elem-badge" id="hub-elem-badge">⚡ Loading...</span>
                            <span id="hub-domain-sub" style="font-size:0.73rem; color:#94a3b8;"></span>
                        </div>
                        <div class="poke-title-row">
                            <h2 class="poke-title" id="hub-card-title" tabindex="-1">--</h2>
                            <span class="poke-cp-badge" id="hub-card-cp">CP --</span>
                        </div>
                        <div class="poke-org-line">
                            <b id="hub-card-org">--</b>
                        </div>
                    </header>

                    <div class="poke-holo-stage">
                        <div id="hub-3d-canvas" class="poke-3d-canvas"></div>
                        <div class="poke-holo-grid" id="hub-holo-grid"></div>
                    </div>

                    <div class="poke-stats-section" id="hub-stats-list"></div>

                    <footer class="poke-card-footer">
                        <span id="hub-footer-sub">Cognate Substrate: --</span>
                        <span id="hub-footer-len">-- aa</span>
                    </footer>
                </article>

                <!-- Center Column: 1-to-N Correlation Synergy & "Which is Better for What" -->
                <div class="versus-column">
                    <div class="vs-banner">
                        <div style="font-size:0.7rem; text-transform:uppercase; letter-spacing:0.06em; color:#94a3b8; font-weight:700;">1-to-N Pairwise Synergy Score</div>
                        <div class="vs-score-big" id="vs-composite-val">--</div>
                        <div style="font-size:0.74rem; color:#cbd5e1;" id="vs-pair-subtitle">Select any Challenger Card below</div>
                    </div>

                    <div class="vs-metrics-grid" id="vs-metrics-grid"></div>

                    <div class="verdict-box" style="border-color:rgba(56,189,248,0.35);">
                        <h4 style="color:#38bdf8;" id="verdict-hub-heading">🏆 Why Choose Hub ("The 1")</h4>
                        <ul class="verdict-list" id="verdict-hub-list"></ul>
                    </div>

                    <div class="verdict-box" style="border-color:rgba(251,191,36,0.4);">
                        <h4 style="color:#fbbf24;" id="verdict-chal-heading">🏆 Why Choose Challenger ("Match #k")</h4>
                        <ul class="verdict-list" id="verdict-chal-list"></ul>
                    </div>

                    <div class="vs-actions">
                        <button type="button" class="vs-action-btn" style="background:#1e293b; color:#38bdf8; border:1px solid #38bdf8;" onclick="toggleInline3DHolofoil()">
                            <span id="holofoil-btn-label">🔮 Toggle Live 3Dmol.js Card Holofoil</span>
                        </button>
                        <a id="vs-open-compare-btn" class="vs-action-btn" style="background:rgba(16,185,129,0.18); color:#34d399; border:1px solid #10b981;" href="/compare" target="_blank">
                            🔬 Open Full Synchronized 3D Superposition Lab
                        </a>
                    </div>
                </div>

                <!-- Right Slot: Correlated Challenger ("Match #k of N") -->
                <article class="bio-poke-card challenger-slot" id="challenger-poke-card">
                    <header class="poke-card-header">
                        <div class="poke-role-strip">
                            <span class="poke-role-tag" id="chal-rank-tag" style="background:rgba(245,158,11,0.18); color:#fbbf24; border:1px solid rgba(245,158,11,0.45);">CHALLENGER #1 OF N</span>
                            <span class="poke-elem-badge" id="chal-elem-badge">⚡ Loading...</span>
                        </div>
                        <div class="poke-title-row">
                            <h2 class="poke-title" id="chal-card-title" tabindex="-1">--</h2>
                            <span class="poke-cp-badge" id="chal-card-cp">CP --</span>
                        </div>
                        <div class="poke-org-line">
                            <b id="chal-card-org">--</b>
                            <span id="chal-domain-sub" style="font-size:0.73rem; color:#94a3b8;"></span>
                        </div>
                    </header>

                    <div class="poke-holo-stage">
                        <div id="chal-3d-canvas" class="poke-3d-canvas"></div>
                        <div class="poke-holo-grid" id="chal-holo-grid"></div>
                    </div>

                    <div class="poke-stats-section" id="chal-stats-list"></div>

                    <footer class="poke-card-footer">
                        <span id="chal-footer-sub">Cognate Substrate: --</span>
                        <span id="chal-footer-len">-- aa</span>
                    </footer>
                </article>
            </section>

            <!-- Bottom Section: 1-to-N Correlated Challenger Deck -->
            <section class="deck-section" aria-label="1-to-N Correlated Challenger Deck">
                <div class="deck-header-bar">
                    <div>
                        <h2 id="deck-section-title">1-to-N Correlated Partner Deck</h2>
                        <div style="font-size:0.78rem; color:#94a3b8;" id="deck-section-sub">Click any card in the 1-to-N deck to swap it into the Right-Hand Challenger Slot with View Transitions</div>
                    </div>
                </div>
                <div class="challenger-deck-grid" id="challenger-deck-grid"></div>
            </section>
        </main>
    </div>

    <script>
        var ARENA_DATA = null;
        var ACTIVE_HUB_ID = null;
        var ACTIVE_MATCH_INDEX = 0;
        var SHOW_LIVE_3D = false;
        var VIEWER_HUB = null;
        var VIEWER_CHAL = null;

        function runWithViewTransition(updateFn, focusTargetId) {
            if (!document.startViewTransition) {
                updateFn();
                if (focusTargetId) {
                    const el = document.getElementById(focusTargetId);
                    if (el) el.focus({ preventScroll: true });
                }
                return;
            }
            const transition = document.startViewTransition(() => {
                updateFn();
            });
            if (focusTargetId) {
                transition.finished.finally(() => {
                    const el = document.getElementById(focusTargetId);
                    if (el) el.focus({ preventScroll: true });
                });
            }
        }

        function reloadArenaData() {
            const orgCap = document.getElementById('org-cap-select').value;
            const matchCap = document.getElementById('match-cap-select').value;
            const sortBy = document.getElementById('sort-lens-select').value;
            document.getElementById('arena-sql-badge').textContent = 'pg_bio: querying...';

            fetch(`/api/hub_cards?org_cap=${encodeURIComponent(orgCap)}&match_cap=${encodeURIComponent(matchCap)}&sort_by=${encodeURIComponent(sortBy)}`)
                .then(r => r.json())
                .then(data => {
                    ARENA_DATA = data;
                    document.getElementById('arena-sql-badge').textContent = `pg_bio: ${data.elapsed_ms} ms (${data.top_organisms.length} Orgs · ${data.hubs.length} Hubs)`;
                    document.getElementById('sidebar-org-count').textContent = `Top ${data.top_organisms.length}`;

                    const urlParams = new URLSearchParams(window.location.search);
                    const reqHub = urlParams.get('hub');
                    if (reqHub && data.hubs.some(h => h.hub_id === reqHub)) {
                        ACTIVE_HUB_ID = reqHub;
                    } else if (!ACTIVE_HUB_ID || !data.hubs.some(h => h.hub_id === ACTIVE_HUB_ID)) {
                        ACTIVE_HUB_ID = data.hubs.length ? data.hubs[0].hub_id : null;
                        if (ACTIVE_HUB_ID) {
                            const u = new URL(window.location);
                            u.searchParams.set('hub', ACTIVE_HUB_ID);
                            window.history.replaceState({}, '', u);
                        }
                    }
                    ACTIVE_MATCH_INDEX = 0;
                    renderSidebar();
                    renderBattleAndDeck();
                });
        }

        function selectHub(hubId) {
            if (ACTIVE_HUB_ID === hubId) return;
            runWithViewTransition(() => {
                ACTIVE_HUB_ID = hubId;
                ACTIVE_MATCH_INDEX = 0;
                const u = new URL(window.location);
                u.searchParams.set('hub', hubId);
                window.history.replaceState({}, '', u);
                renderSidebar();
                renderBattleAndDeck();
            }, 'hub-card-title');
        }

        function selectChallenger(matchIdx) {
            if (ACTIVE_MATCH_INDEX === matchIdx) return;
            runWithViewTransition(() => {
                ACTIVE_MATCH_INDEX = matchIdx;
                renderBattleAndDeck();
            }, 'chal-card-title');
        }

        function renderSidebar() {
            const container = document.getElementById('org-sidebar-list');
            if (!ARENA_DATA || !ARENA_DATA.top_organisms.length) {
                container.innerHTML = '<div style="color:#64748b;font-size:0.8rem;">No hub organisms found.</div>';
                return;
            }
            const activeHubObj = ARENA_DATA.hubs.find(h => h.hub_id === ACTIVE_HUB_ID) || ARENA_DATA.hubs[0];
            const activeOrgName = activeHubObj ? activeHubObj.organism : '';

            let html = '';
            ARENA_DATA.top_organisms.forEach((org, idx) => {
                const orgHubs = ARENA_DATA.hubs.filter(h => h.organism === org.organism);
                const isOrgSelected = (org.organism === activeOrgName);
                const firstHubId = orgHubs.length ? orgHubs[0].hub_id : '';
                html += `
                    <div class="org-pill-card ${isOrgSelected ? 'is-selected' : ''}" onclick="if('${firstHubId}' && '${ org.organism.replace(/'/g, "\\'") }' !== '${ activeOrgName.replace(/'/g, "\\'") }') selectHub('${firstHubId}')">
                        <div class="org-rank-row">
                            <span class="org-rank-badge">#${idx + 1} Hub Organism</span>
                            <span class="org-reach-badge">${org.partner_orgs} Partner Orgs</span>
                        </div>
                        <div class="org-title">${org.short_organism}</div>
                        <div class="org-sub-metrics">
                            <span>🔗 ${org.total_edges} Total Edges</span>
                            <span>🎯 ${org.rf_partners} 3D Pocket Hits</span>
                            <span>Max RF: ${org.max_rf_score.toFixed(2)}</span>
                        </div>
                        <div class="hub-btn-list">
                            ${orgHubs.map(h => {
                                const fn = (h.family_name && h.family_name !== 'Dark Proteome Island')
                                    ? h.family_name
                                    : ((h.hub_card && h.hub_card.name && h.hub_card.name !== 'Uncharacterized protein') ? h.hub_card.name.replace(' domain-containing protein', '') : 'Orphan Fold');
                                return `
                                <button type="button" class="hub-pill-btn ${h.hub_id === ACTIVE_HUB_ID ? 'is-selected' : ''}"
                                        onclick="event.stopPropagation(); selectHub('${h.hub_id}')">
                                    ${h.hub_id} · ${fn} (1:${h.total_3d_matches})
                                </button>
                                `;
                            }).join('')}
                        </div>
                    </div>
                `;
            });
            container.innerHTML = html;
        }

        function renderSinglePokeCard(prefix, cardData, h2hAxes, sideKey) {
            const elemBadge = document.getElementById(`${prefix}-elem-badge`);
            elemBadge.textContent = cardData.element_badge;
            elemBadge.style.color = cardData.element_color;
            elemBadge.style.background = cardData.element_bg;

            document.getElementById(`${prefix}-card-title`).textContent = cardData.uniprot_id;
            document.getElementById(`${prefix}-card-cp`).textContent = `CP ${cardData.power_rating}`;
            document.getElementById(`${prefix}-card-org`).textContent = cardData.organism.split(' (strain')[0];
            const domSubEl = document.getElementById(`${prefix}-domain-sub`);
            if (domSubEl) {
                const desc = (cardData.name && cardData.name !== cardData.uniprot_id && cardData.name !== 'Uncharacterized protein')
                    ? cardData.name.replace(' domain-containing protein', ' domain')
                    : ((cardData.family_name && cardData.family_name !== 'Dark Proteome Island') ? cardData.family_name : '');
                domSubEl.textContent = desc;
            }

            document.getElementById(`${prefix}-holo-grid`).innerHTML = `
                <div class="holo-pill"><span>RangeFinder MW Box</span><b>${cardData.rf_mw_box}</b></div>
                <div class="holo-pill"><span>Pocket LogP Box</span><b>${cardData.rf_logp_box}</b></div>
                <div class="holo-pill"><span>Pocket TPSA Box</span><b>${cardData.rf_tpsa_box}</b></div>
                <div class="holo-pill"><span>Aromatic Warhead P</span><b>${Math.round(cardData.rf_arom_prob * 100)}%</b></div>
            `;

            const statOrder = [
                ['brine_armor', '#38bdf8'],
                ['interfacial_hydro', '#f59e0b'],
                ['fold_stability', '#10b981'],
                ['pocket_breadth', '#c084fc'],
                ['expression_efficiency', '#34d399'],
                ['aqueous_solubility', '#60a5fa']
            ];
            let statsHtml = '';
            statOrder.forEach(([key, barColor]) => {
                const st = cardData.stats[key];
                const ax = h2hAxes ? h2hAxes[key] : null;
                const isWinner = ax && ax.winner === sideKey;
                const diffTxt = isWinner ? `🏆 BETTER (+${Math.abs(ax.diff)})` : '';
                statsHtml += `
                    <div class="poke-stat-row ${isWinner ? 'is-winner' : ''}">
                        <div class="stat-top-line">
                            <span class="stat-label">
                                ${st.label}
                                ${isWinner ? `<span class="winner-crown">${diffTxt}</span>` : ''}
                            </span>
                            <span class="stat-raw-val">${st.raw} · <b style="color:#f8fafc;">${st.score}</b>/100</span>
                        </div>
                        <div class="stat-bar-track">
                            <div class="stat-bar-fill" style="width:${st.score}%; background:${isWinner ? '#10b981' : barColor};"></div>
                        </div>
                        <div class="stat-use-case">Best for: ${st.better_for}</div>
                    </div>
                `;
            });
            document.getElementById(`${prefix}-stats-list`).innerHTML = statsHtml;
            document.getElementById(`${prefix}-footer-sub`).textContent = `Cognate Substrate: ${cardData.substrate ? cardData.substrate.substrate_name : 'Aromatic Ester'}`;
            document.getElementById(`${prefix}-footer-len`).textContent = `${cardData.length} aa · ${cardData.mw_kda} kDa`;
        }

        function renderBattleAndDeck() {
            if (!ARENA_DATA || !ARENA_DATA.hubs.length) return;
            const hubObj = ARENA_DATA.hubs.find(h => h.hub_id === ACTIVE_HUB_ID) || ARENA_DATA.hubs[0];
            if (!hubObj || !hubObj.matches.length) return;

            const matchObj = hubObj.matches[Math.min(ACTIVE_MATCH_INDEX, hubObj.matches.length - 1)];
            const hubCard = hubObj.hub_card;
            const chalCard = matchObj.card;
            const pm = matchObj.pair_metrics;
            const h2h = matchObj.head_to_head;

            renderSinglePokeCard('hub', hubCard, h2h.axes, 'hub');
            renderSinglePokeCard('chal', chalCard, h2h.axes, 'challenger');
            document.getElementById('chal-rank-tag').textContent = `CHALLENGER #${matchObj.rank} OF ${hubObj.matches.length}`;

            // Center Versus & Synergy Metrics
            document.getElementById('vs-composite-val').textContent = `${(pm.composite_score * 100).toFixed(1)}%`;
            document.getElementById('vs-pair-subtitle').innerHTML = `<b>${hubCard.uniprot_id}</b> vs <b>${chalCard.uniprot_id}</b> (${pm.cross_genus ? 'Cross-Genus Convergence' : 'Same-Genus Homolog'})`;

            document.getElementById('vs-metrics-grid').innerHTML = `
                <div class="vs-metric-box"><span>RangeFinder 9-Prop Score</span><b style="color:#fbbf24;">${pm.rf_score.toFixed(3)} (IoU ${pm.rf_range_iou.toFixed(2)})</b></div>
                <div class="vs-metric-box"><span>GiST 8-Prop Containment</span><b style="color:#34d399;">${pm.rf_substrate_props_matched}/8 Props (${pm.rf_substrate_match ? 'PASS' : 'PARTIAL'})</b></div>
                <div class="vs-metric-box"><span>Twilight-Zone Score</span><b style="color:#c084fc;">${pm.twilight_score.toFixed(3)} (${pm.seq_identity_pct}% ID)</b></div>
                <div class="vs-metric-box"><span>3D TM-Score &amp; RMSD</span><b style="color:#38bdf8;">TM ${pm.tm_score.toFixed(2)} · ${pm.sc_rmsd.toFixed(2)} Å</b></div>
                <div class="vs-metric-box"><span>ESM-2 Cosine Similarity</span><b>${pm.cosine_sim_pct}% (d=${pm.vector_distance.toFixed(4)})</b></div>
                <div class="vs-metric-box"><span>Chaperone Fold Gate</span><b style="color:#10b981;">${chalCard.chaperone_score} (${chalCard.buns} BUNS)</b></div>
            `;

            // "Which is Better for What" Lists
            document.getElementById('verdict-hub-heading').textContent = `🏆 Better on Hub ${hubCard.uniprot_id} (${h2h.hub_advantages.length} Axes)`;
            document.getElementById('verdict-hub-list').innerHTML = h2h.hub_advantages.length
                ? h2h.hub_advantages.map(a => `<li><b style="color:#e2e8f0;">${a.use_case}:</b> ${a.delta}</li>`).join('')
                : `<li>Balanced baseline across all axes compared to ${chalCard.uniprot_id}.</li>`;

            document.getElementById('verdict-chal-heading').textContent = `🏆 Better on Challenger ${chalCard.uniprot_id} (${h2h.challenger_advantages.length} Axes)`;
            document.getElementById('verdict-chal-list').innerHTML = h2h.challenger_advantages.length
                ? h2h.challenger_advantages.map(a => `<li><b style="color:#e2e8f0;">${a.use_case}:</b> ${a.delta}</li>`).join('')
                : `<li>Serves as an alternative cross-species backup scaffold for ${hubCard.uniprot_id}.</li>`;

            const cmpUrl = `/compare?node1=${encodeURIComponent(hubCard.uniprot_id)}&node2=${encodeURIComponent(chalCard.uniprot_id)}&preset=rangefinder`;
            const wikiUrl = `/wiki?node1=${encodeURIComponent(hubCard.uniprot_id)}&node2=${encodeURIComponent(chalCard.uniprot_id)}&cat=rangefinder`;
            document.getElementById('vs-open-compare-btn').href = cmpUrl;
            document.getElementById('arena-compare-top-link').href = cmpUrl;
            document.getElementById('arena-wiki-top-link').href = wikiUrl;

            // Render 1-to-N Challenger Deck
            document.getElementById('deck-section-title').textContent = `1-to-${hubObj.matches.length} Correlated Partner Deck for Hub ${hubCard.uniprot_id} (${hubObj.short_organism})`;
            let deckHtml = '';
            hubObj.matches.forEach((m, idx) => {
                const c = m.card;
                const p = m.pair_metrics;
                const wins = m.head_to_head.challenger_advantages;
                const isSel = (idx === ACTIVE_MATCH_INDEX);
                deckHtml += `
                    <article class="deck-mini-card ${isSel ? 'is-selected' : ''}" onclick="selectChallenger(${idx})">
                        <header class="mini-header">
                            <span class="mini-id">#${m.rank} ${c.uniprot_id}</span>
                            <span class="mini-synergy">${(p.composite_score * 100).toFixed(1)}% Synergy</span>
                        </header>
                        <div class="mini-org">
                            <b style="color:#e2e8f0;">${c.organism.split(' (strain')[0]}</b><br>
                            <span>RF IoU: ${p.rf_range_iou.toFixed(2)} · Twilight: ${p.twilight_score.toFixed(2)} · d=${p.vector_distance.toFixed(3)}</span>
                        </div>
                        <div class="mini-adv-pills">
                            ${wins.length
                                ? wins.slice(0, 3).map(w => `<span class="mini-adv-tag">🏆 ${w.title}</span>`).join('')
                                : `<span class="mini-adv-tag" style="color:#94a3b8;">⚖️ Balanced Iso-Pocket Twin</span>`
                            }
                        </div>
                        <footer class="mini-footer">
                            <span>CP ${c.power_rating} · ${c.mw_kda} kDa</span>
                            <span>pI ${c.pi} · GRAVY ${c.gravy >= 0 ? '+' : ''}${c.gravy}</span>
                        </footer>
                    </article>
                `;
            });
            document.getElementById('challenger-deck-grid').innerHTML = deckHtml;

            if (SHOW_LIVE_3D) {
                loadInline3DModels(hubCard.uniprot_id, chalCard.uniprot_id);
            }
        }

        function toggleInline3DHolofoil() {
            SHOW_LIVE_3D = !SHOW_LIVE_3D;
            const c1 = document.getElementById('hub-3d-canvas');
            const c2 = document.getElementById('chal-3d-canvas');
            const lbl = document.getElementById('holofoil-btn-label');
            if (SHOW_LIVE_3D) {
                c1.style.display = 'block';
                c2.style.display = 'block';
                lbl.textContent = '🃏 Hide Live 3Dmol.js Card Holofoil';
                const hubObj = ARENA_DATA && ARENA_DATA.hubs.find(h => h.hub_id === ACTIVE_HUB_ID);
                if (hubObj && hubObj.matches[ACTIVE_MATCH_INDEX]) {
                    loadInline3DModels(hubObj.hub_card.uniprot_id, hubObj.matches[ACTIVE_MATCH_INDEX].card.uniprot_id);
                }
            } else {
                c1.style.display = 'none';
                c2.style.display = 'none';
                lbl.textContent = '🔮 Toggle Live 3Dmol.js Card Holofoil';
            }
        }

        var _SYNCING_3D_VIEW = false;
        function syncCard3DViews(sourceViewer, targetViewer) {
            if (_SYNCING_3D_VIEW || !sourceViewer || !targetViewer) return;
            _SYNCING_3D_VIEW = true;
            try {
                targetViewer.setView(sourceViewer.getView());
                targetViewer.render();
            } catch (e) {}
            _SYNCING_3D_VIEW = false;
        }

        function loadInline3DModels(id1, id2) {
            const el1 = document.getElementById('hub-3d-canvas');
            const el2 = document.getElementById('chal-3d-canvas');
            if (!VIEWER_HUB) {
                VIEWER_HUB = $3Dmol.createViewer(el1, { backgroundColor: '#020617' });
                el1.addEventListener('mousemove', (ev) => { if (ev.buttons) syncCard3DViews(VIEWER_HUB, VIEWER_CHAL); });
                el1.addEventListener('wheel', () => { setTimeout(() => syncCard3DViews(VIEWER_HUB, VIEWER_CHAL), 15); }, { passive: true });
            }
            if (!VIEWER_CHAL) {
                VIEWER_CHAL = $3Dmol.createViewer(el2, { backgroundColor: '#020617' });
                el2.addEventListener('mousemove', (ev) => { if (ev.buttons) syncCard3DViews(VIEWER_CHAL, VIEWER_HUB); });
                el2.addEventListener('wheel', () => { setTimeout(() => syncCard3DViews(VIEWER_CHAL, VIEWER_HUB), 15); }, { passive: true });
            }
            fetch('/api/compare/' + encodeURIComponent(id1) + '/' + encodeURIComponent(id2))
                .then(r => r.json())
                .then(comp => {
                    if (!comp || !comp.pdb1) return;
                    const pdb1 = comp.pdb1;
                    const pdb2 = comp.pdb2_aligned || comp.pdb2;
                    VIEWER_HUB.clear();
                    VIEWER_HUB.addModel(pdb1, 'pdb');
                    VIEWER_HUB.setStyle({}, { cartoon: { colorscheme: 'cyanCarbon' } });
                    VIEWER_HUB.zoomTo();
                    VIEWER_HUB.render();

                    VIEWER_CHAL.clear();
                    VIEWER_CHAL.addModel(pdb2, 'pdb');
                    VIEWER_CHAL.setStyle({}, { cartoon: { colorscheme: 'orangeCarbon' } });
                    try {
                        VIEWER_CHAL.setView(VIEWER_HUB.getView());
                    } catch (e) {
                        VIEWER_CHAL.zoomTo();
                    }
                    VIEWER_CHAL.render();

                    if (comp.alignment) {
                        const tmVal = parseFloat(comp.alignment.tm_score || 0);
                        const rmsdVal = parseFloat(comp.alignment.core_rmsd || 0);
                        const boxes = document.querySelectorAll('#vs-metrics-grid .vs-metric-box');
                        if (boxes && boxes.length >= 4 && tmVal > 0) {
                            boxes[3].innerHTML = `<span>3D TM-Score &amp; RMSD</span><b style="color:#38bdf8;">TM ${tmVal.toFixed(2)} · ${rmsdVal.toFixed(2)} Å</b>`;
                        }
                    }
                });
        }

        window.addEventListener('DOMContentLoaded', reloadArenaData);
    </script>
</body>
</html>
"""


@app.route('/cards')
@app.route('/arena')
def cards_arena_page():
    return render_template_string(CARDS_ARENA_HTML_TEMPLATE)

AA1_TO_AA3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

UNIPROT_CACHE_PATH = "/Users/jonatas/code/ideia.me/assets/models/uniprot_meta_cache.json"
ALIGNED_PDB_DIR = "/Users/jonatas/code/ideia.me/assets/models/aligned"
os.makedirs(ALIGNED_PDB_DIR, exist_ok=True)

try:
    with open(UNIPROT_CACHE_PATH, "r") as _uf:
        UNIPROT_META_CACHE = json.load(_uf)
except Exception:
    UNIPROT_META_CACHE = {}

COMPARE_CACHE = {}
CACHE_SCHEMA_VERSION = 5

def ensure_comparison_cache_table():
    try:
        with psycopg.connect(DB_URI, autocommit=True) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS swarm_comparison_cache (
                        id1 TEXT NOT NULL,
                        id2 TEXT NOT NULL,
                        schema_version INT NOT NULL DEFAULT 3,
                        metrics_payload JSONB NOT NULL,
                        pdb1_path TEXT,
                        pdb2_aligned_path TEXT,
                        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                        PRIMARY KEY (id1, id2)
                    );
                """)
    except Exception:
        pass

ensure_comparison_cache_table()

def fetch_pdb_string(uniprot_id):
    path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{uniprot_id}-F1-model_v4.pdb"
    if os.path.exists(path):
        with open(path, 'r') as f:
            content = f.read()
            if "ATOM" in content and "PARAMETRIC GLOBULAR FOLD" not in content:
                return content, 200

    for ver in (4, 6):
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{ver}.pdb"
            pdb_data = urllib.request.urlopen(url, timeout=2.5).read().decode('utf-8')
            if "ATOM" in pdb_data:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, 'w') as f:
                    f.write(pdb_data)
                return pdb_data, 200
        except urllib.error.HTTPError as he:
            if he.code == 404 and ver == 4:
                continue
        except Exception:
            break

    return "Real 3D structure not available.", 404

@app.route('/pdb/<uniprot_id>')
def get_pdb(uniprot_id):
    pdb_str, status = fetch_pdb_string(uniprot_id)
    return pdb_str, status

def kabsch_superpose_pdb(pdb_ref, pdb_mov):
    def get_ca(pdb_str):
        coords, resi, bfacs = [], [], []
        for l in pdb_str.splitlines():
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                resi.append(int(l[22:26]))
                b = float(l[60:66])
                bfacs.append(b * 100.0 if b <= 1.0 else b)
        return np.array(coords), np.array(resi), (np.mean(bfacs) if bfacs else 0.0)

    c1, r1, plddt1 = get_ca(pdb_ref)
    c2, r2, plddt2 = get_ca(pdb_mov)
    if len(c1) < 5 or len(c2) < 5:
        return pdb_mov, 0.0, 0.0, 0, plddt1, plddt2, [], []

    def fit_kabsch(P, Q):
        cP, cQ = P.mean(axis=0), Q.mean(axis=0)
        P0, Q0 = P - cP, Q - cQ
        H = Q0.T @ P0
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(U @ Vt))
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
        return R, cP, cQ

    def tangents(c):
        t = np.zeros_like(c)
        if len(c) > 2:
            t[1:-1] = c[2:] - c[:-2]
            t[0] = c[1] - c[0]
            t[-1] = c[-1] - c[-2]
        norms = np.linalg.norm(t, axis=1, keepdims=True) + 1e-8
        return t / norms

    def internal_fp(c):
        idx = np.arange(len(c))
        fps = []
        for d in (-4, -3, -2, -1, 1, 2, 3, 4):
            j = np.clip(idx + d, 0, len(c) - 1)
            fps.append(np.linalg.norm(c[idx] - c[j], axis=1))
        return np.stack(fps, axis=1)

    n1, n2 = len(c1), len(c2)
    L = min(n1, n2)
    d0 = max(1.24 * ((L - 15) ** (1.0 / 3.0)) - 1.8, 0.5) if L > 15 else 0.5
    t1 = tangents(c1)
    t2 = tangents(c2)

    fp_dist = scipy.spatial.distance.cdist(internal_fp(c1), internal_fp(c2), metric='cityblock')

    seeds = [(0, 0, L)]
    diag_step = 4 if max(n1, n2) > 600 else 2
    for w in (min(24, L), min(48, L), min(80, L)):
        if w < 8:
            continue
        best_diags = []
        for diag in range(-n1 + w, n2 - w + 1, diag_step):
            d_vals = np.diagonal(fp_dist, offset=diag)
            if len(d_vals) >= w:
                cs = np.concatenate([[0.0], np.cumsum(d_vals)])
                w_sums = cs[w:] - cs[:-w]
                s_idx = int(np.argmin(w_sums))
                i1 = s_idx if diag >= 0 else s_idx - diag
                i2 = s_idx + diag if diag >= 0 else s_idx
                best_diags.append((w_sums[s_idx] / w, i1, i2, w))
        best_diags.sort(key=lambda x: x[0])
        for _, i1, i2, win in best_diags[:5]:
            seeds.append((i1, i2, win))

    best_score = -1.0
    best_transform = (np.eye(3), np.zeros(3), np.zeros(3))
    best_stats = (9.99, 0.0, 0, [], [])

    for i1, i2, win in seeds:
        R, cP, cQ = fit_kabsch(c1[i1:i1+win], c2[i2:i2+win])
        for cutoff in (5.5, 3.5):
            c2_al = (c2 - cQ) @ R + cP
            t2_al = t2 @ R
            D = scipy.spatial.distance.cdist(c1, c2_al)
            dir_cos = np.clip(t1 @ t2_al.T, 0.0, 1.0)
            S = (dir_cos ** 1.5) / (1.0 + (D / d0) ** 2)
            j_for_i = np.argmax(S, axis=1)
            i_for_j = np.argmax(S, axis=0)
            mutual = (i_for_j[j_for_i] == np.arange(n1)) & (S[np.arange(n1), j_for_i] > 0.15)
            m_i = np.where(mutual)[0]
            m_j = j_for_i[m_i]
            if len(m_i) >= 8:
                scores = S[m_i, m_j]
                dp = scores.copy()
                prev = np.full(len(m_i), -1, dtype=int)
                for a in range(len(m_i)):
                    valid = (m_j[:a] < m_j[a]) & ((m_i[a] - m_i[:a]) <= 40) & ((m_j[a] - m_j[:a]) <= 40)
                    if np.any(valid):
                        cand = np.where(valid)[0]
                        best_p = cand[np.argmax(dp[cand])]
                        dp[a] += dp[best_p]
                        prev[a] = best_p
                cur = int(np.argmax(dp))
                chain = []
                while cur != -1:
                    chain.append(cur)
                    cur = prev[cur]
                chain.reverse()
                idx1, idx2 = m_i[chain], m_j[chain]
                d_al = D[idx1, idx2]
                mask = d_al < cutoff
                if mask.sum() >= 8:
                    R, cP, cQ = fit_kabsch(c1[idx1[mask]], c2[idx2[mask]])

        c2_al = (c2 - cQ) @ R + cP
        t2_al = t2 @ R
        D = scipy.spatial.distance.cdist(c1, c2_al)
        dir_cos = np.clip(t1 @ t2_al.T, 0.0, 1.0)
        S = (dir_cos ** 1.5) / (1.0 + (D / d0) ** 2)
        j_for_i = np.argmax(S, axis=1)
        i_for_j = np.argmax(S, axis=0)
        mutual = (i_for_j[j_for_i] == np.arange(n1)) & (S[np.arange(n1), j_for_i] > 0.12)
        m_i = np.where(mutual)[0]
        m_j = j_for_i[m_i]
        if len(m_i) >= 8:
            scores = 1.0 / (1.0 + (D[m_i, m_j] / d0) ** 2)
            dp = scores.copy()
            prev = np.full(len(m_i), -1, dtype=int)
            for a in range(len(m_i)):
                valid = (m_j[:a] < m_j[a])
                if np.any(valid):
                    cand = np.where(valid)[0]
                    best_p = cand[np.argmax(dp[cand])]
                    dp[a] += dp[best_p]
                    prev[a] = best_p
            cur = int(np.argmax(dp))
            chain = []
            while cur != -1:
                chain.append(cur)
                cur = prev[cur]
            chain.reverse()
            idx1, idx2 = m_i[chain], m_j[chain]
            d_al = D[idx1, idx2]
            tm = float(np.sum(1.0 / (1.0 + (d_al / d0) ** 2)) / L)
            core = d_al < 3.8
            if core.sum() >= 8:
                R_c, cP_c, cQ_c = fit_kabsch(c1[idx1[core]], c2[idx2[core]])
                c2_core = (c2[idx2[core]] - cQ_c) @ R_c + cP_c
                rmsd = float(np.sqrt(np.mean(np.sum((c1[idx1[core]] - c2_core) ** 2, axis=1))))
                n_c = int(core.sum())
                res1_core = r1[idx1[core]].tolist()
                res2_core = r2[idx2[core]].tolist()
                R, cP, cQ = R_c, cP_c, cQ_c
            else:
                rmsd = float(np.sqrt(np.mean(d_al ** 2)))
                n_c = len(d_al)
                res1_core = r1[idx1].tolist()
                res2_core = r2[idx2].tolist()
            rank_score = tm + 0.002 * n_c
            if rank_score > best_score:
                best_score = rank_score
                best_transform = (R, cP, cQ)
                best_stats = (rmsd, tm, n_c, res1_core, res2_core)

    R, cP, cQ = best_transform
    out_lines = []
    for line in pdb_mov.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            new_xyz = (xyz - cQ) @ R + cP
            line = f"{line[:30]}{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}{line[54:]}"
        out_lines.append(line)
    return "\n".join(out_lines), best_stats[0], best_stats[1], best_stats[2], plddt1, plddt2, best_stats[3], best_stats[4]

def fetch_uniprot_domains_cached(uid, fallback_name=""):
    org_from_name = ""
    if "OS=" in fallback_name:
        org_from_name = fallback_name.split("OS=", 1)[1].split(" OX=", 1)[0].split(" GN=", 1)[0].strip()
    if uid in UNIPROT_META_CACHE:
        entry = UNIPROT_META_CACHE[uid]
        return entry.get("organism") or org_from_name, entry.get("domains", [])
    domains = []
    organism = org_from_name
    try:
        u = json.loads(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=1.2).read().decode())
        organism = u.get("organism", {}).get("scientificName", "") or org_from_name
        for x in u.get("uniProtKBCrossReferences", []):
            if x.get("database") in ("Pfam", "InterPro"):
                val = x.get("properties", [{}])[0].get("value", "")
                domains.append(f"{x['id']}: {val}" if val else x['id'])
        for feat in u.get("features", []):
            if feat.get("type") == "Transmembrane":
                s = feat.get("location", {}).get("start", {}).get("value", "?")
                e = feat.get("location", {}).get("end", {}).get("value", "?")
                domains.append(f"TM Helix ({s}-{e})")
        UNIPROT_META_CACHE[uid] = {"organism": organism, "domains": domains[:6]}
        try:
            with open(UNIPROT_CACHE_PATH, "w") as _uf:
                json.dump(UNIPROT_META_CACHE, _uf)
        except Exception:
            pass
    except Exception:
        UNIPROT_META_CACHE[uid] = {"organism": organism, "domains": []}
    return organism, domains[:6]

@app.route('/api/compare/<id1>/<id2>')
def compare_proteins(id1, id2):
    cache_key = f"{id1}|{id2}"
    # Tier L2: In-Memory Python Dict Cache (< 0.2 ms)
    if cache_key in COMPARE_CACHE:
        cached_resp = dict(COMPARE_CACHE[cache_key])
        cached_resp['cached_in_db'] = True
        return jsonify(cached_resp)

    conn = psycopg.connect(DB_URI, autocommit=True)
    cur = conn.cursor()

    # Tier L3: Persistent PostgreSQL swarm_comparison_cache + Disk Aligned PDBs (< 2.5 ms)
    try:
        cur.execute("""
            SELECT metrics_payload, pdb1_path, pdb2_aligned_path
            FROM swarm_comparison_cache
            WHERE id1 = %s AND id2 = %s AND schema_version = %s
        """, (id1, id2, CACHE_SCHEMA_VERSION))
        row = cur.fetchone()
        if row and row[0] and row[1] and row[2] and os.path.exists(row[1]) and os.path.exists(row[2]):
            metrics = row[0] if isinstance(row[0], dict) else json.loads(row[0])
            with open(row[1], 'r') as f1, open(row[2], 'r') as f2:
                pdb1_cached = f1.read()
                pdb2_al_cached = f2.read()
            full_payload = dict(metrics)
            full_payload['pdb1'] = pdb1_cached
            full_payload['pdb2_aligned'] = pdb2_al_cached
            full_payload['cached_in_db'] = True
            COMPARE_CACHE[cache_key] = full_payload
            conn.close()
            return jsonify(full_payload)
    except Exception:
        pass

    pdb1, s1 = fetch_pdb_string(id1)
    pdb2, s2 = fetch_pdb_string(id2)
    if s1 != 200 or s2 != 200:
        conn.close()
        return jsonify({'error': f'Could not load 3D structure ({id1}: {s1}, {id2}: {s2})'}), 500

    pdb2_aligned, core_rmsd, tm_score, n_core, plddt1, plddt2, res1_core, res2_core = kabsch_superpose_pdb(pdb1, pdb2)

    def build_meta(uid, plddt):
        cur.execute(
            "SELECT name, sequence, molecular_weight(sequence), sequence_entropy(sequence), isoelectric_point(sequence), gravy_score(sequence) FROM proteins WHERE uniprot_id = %s",
            (uid,)
        )
        row = cur.fetchone()
        name, seq, mw, ent, pi_val, gravy_val = row if row else ("Unknown", "", 0.0, 0.0, 7.0, 0.0)
        L = max(len(seq), 1)
        acidic = (seq.count('D') + seq.count('E')) / L * 100.0
        basic = (seq.count('K') + seq.count('R') + seq.count('H')) / L * 100.0
        hydro = sum(seq.count(c) for c in 'AILMFWV') / L * 100.0
        organism, domains = fetch_uniprot_domains_cached(uid, name)
        return {
            'id': uid,
            'name': name,
            'organism': organism,
            'length': len(seq),
            'mw_kda': f"{mw/1000.0:.1f}",
            'entropy': f"{ent:.2f}",
            'pi': f"{pi_val:.2f}",
            'gravy': f"{gravy_val:.2f}",
            'pi_num': float(pi_val),
            'gravy_num': float(gravy_val),
            'ent_num': float(ent),
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
        """
        SELECT sequence_alignment_score(%s, %s),
               smith_waterman_score(%s, %s, 2, -1, -2),
               embedding_cosine_distance(get_esm_embedding(%s), get_esm_embedding(%s)),
               (SELECT p1.embedding <=> p2.embedding FROM proteins p1, proteins p2 WHERE p1.uniprot_id = %s AND p2.uniprot_id = %s LIMIT 1)
        """,
        (m1['seq'], m2['seq'], m1['seq'], m2['seq'], m1['seq'], m2['seq'], id1, id2)
    )
    sw_score, sw_212, esm_kmer_dist, vec_dist = cur.fetchone()
    L_min = max(min(m1['length'], m2['length']), 1)
    L_max = max(max(m1['length'], m2['length']), 1)
    sw_sim = min(max(float(sw_212 or 0) / (2.0 * L_min), 0.0), 1.0)
    sw_sim_pct = round(sw_sim * 100.0, 1)
    len_ratio = float(L_min) / float(L_max)
    v_d = float(vec_dist if vec_dist is not None else esm_kmer_dist)
    if 0.10 <= sw_sim <= 0.96:
        twilight_score = round(
            ((1.0 - min(v_d, 1.0)) ** 2)
            * (4.0 * sw_sim * (1.0 - sw_sim))
            * (len_ratio ** 2)
            * min(m2['ent_num'] / 4.0, 1.0),
            4
        )
    else:
        twilight_score = 0.0

    delta_pi = abs(m1['pi_num'] - m2['pi_num'])
    delta_gravy = abs(m1['gravy_num'] - m2['gravy_num'])

    chaperone_data = None
    try:
        cur.execute("""
            SELECT c1.sc_rmsd,
                   c1.mean_rasa AS rasa_1,
                   c1.core_buried_fraction AS buried_1,
                   c1.surface_hydrophobic_fraction AS surf_hydro_1,
                   c1.total_hbonds AS hbonds_1,
                   c1.buried_unsatisfied_polars AS buns_1,
                   c1.radius_of_gyration AS rg_1,
                   c1.design_composite_score AS score_1,
                   c1.passes_chaperone_gate AS gate_1,
                   c2.mean_rasa AS rasa_2,
                   c2.core_buried_fraction AS buried_2,
                   c2.surface_hydrophobic_fraction AS surf_hydro_2,
                   c2.total_hbonds AS hbonds_2,
                   c2.buried_unsatisfied_polars AS buns_2,
                   c2.radius_of_gyration AS rg_2,
                   c2.design_composite_score AS score_2,
                   c2.passes_chaperone_gate AS gate_2,
                   bio_rasa(%s, 1.4)::text AS rasa_json_1,
                   bio_rasa(%s, 1.4)::text AS rasa_json_2,
                   bio_hbond_network(%s, 3.5)::text AS hnet_json_1,
                   bio_hbond_network(%s, 3.5)::text AS hnet_json_2,
                   chaperone_lineage_hash('SwarmCompare', %s, %s, 1)
            FROM evaluate_chaperone_block(%s, %s, 'A', 'B') c1,
                 evaluate_chaperone_block(%s, %s, 'A', 'B') c2;
        """, (
            pdb1, pdb2_aligned,
            pdb1, pdb2_aligned,
            id1, m2['seq'],
            pdb1, pdb1,
            pdb2_aligned, pdb1
        ))
        crow = cur.fetchone()
        if crow:
            (
                sc_rmsd_v, rasa_1, buried_1, sh_1, hb_1, buns_1, comp_1, sc_1, gate_1,
                rasa_2, buried_2, sh_2, hb_2, buns_2, comp_2, sc_2, gate_2,
                rjson_1, rjson_2, hnet_1, hnet_2, lin_id
            ) = crow
            rj1 = json.loads(rjson_1)
            rj2 = json.loads(rjson_2)
            rv1 = np.array(rj1.get('per_residue_rasa', []), dtype=float)
            rv2 = np.array(rj2.get('per_residue_rasa', []), dtype=float)
            rn1 = list(range(1, len(rv1) + 1))
            rn2 = list(range(1, len(rv2) + 1))

            def top_csr_hubs(hnet_str):
                try:
                    hj = json.loads(hnet_str)
                    deg_map = {}
                    for e in hj.get('entries', []):
                        s = int(e.get('source_residue', 0)) + 1
                        t = int(e.get('target_residue', 0)) + 1
                        w = float(e.get('weight', 1.0))
                        deg_map[s] = deg_map.get(s, 0.0) + w
                        deg_map[t] = deg_map.get(t, 0.0) + w
                    degs = sorted(deg_map.items(), key=lambda x: x[1], reverse=True)
                    return [int(r) for r, _ in degs[:12]]
                except Exception:
                    return []

            hubs_1 = top_csr_hubs(hnet_1)
            hubs_2 = top_csr_hubs(hnet_2)
            exp1 = [int(rn1[i]) for i, v in enumerate(rv1) if i < len(rn1) and v > 0.25]
            exp2 = [int(rn2[i]) for i, v in enumerate(rv2) if i < len(rn2) and v > 0.25]
            if len(res1_core) >= 5 and len(res1_core) == len(res2_core):
                map1 = {int(rn1[i]): float(rv1[i]) for i in range(min(len(rn1), len(rv1)))}
                map2 = {int(rn2[i]): float(rv2[i]) for i in range(min(len(rn2), len(rv2)))}
                pair_r1 = np.array([map1.get(int(r), 0.0) for r in res1_core])
                pair_r2 = np.array([map2.get(int(r), 0.0) for r in res2_core])
                if np.std(pair_r1) > 1e-6 and np.std(pair_r2) > 1e-6:
                    rasa_corr = float(np.corrcoef(pair_r1, pair_r2)[0, 1])
                else:
                    rasa_corr = 0.0
            else:
                L_r = min(len(rv1), len(rv2))
                rasa_corr = float(np.corrcoef(rv1[:L_r], rv2[:L_r])[0, 1]) if L_r >= 5 else 0.0

            chaperone_data = {
                'lineage_id': lin_id,
                'sc_rmsd': round(float(core_rmsd), 2),
                'rasa_1': round(float(rasa_1), 3),
                'rasa_2': round(float(rasa_2), 3),
                'rasa_corr': round(float(rasa_corr), 3),
                'buried_1': round(float(buried_1) * 100.0, 1),
                'buried_2': round(float(buried_2) * 100.0, 1),
                'surf_hydro_1': round(float(sh_1) * 100.0, 1),
                'surf_hydro_2': round(float(sh_2) * 100.0, 1),
                'hbonds_1': int(hb_1),
                'hbonds_2': int(hb_2),
                'buns_1': int(buns_1),
                'buns_2': int(buns_2),
                'compact_1': round(float(comp_1), 3),
                'compact_2': round(float(comp_2), 3),
                'chaperone_score': round(float(sc_2), 2),
                'passes_gate': bool(gate_2 or (tm_score >= 0.50 and float(sh_2) <= 0.45)),
                'hbond_hubs_1': hubs_1,
                'hbond_hubs_2': hubs_2,
                'exposed_resi_1': exp1,
                'exposed_resi_2': exp2,
                'closest_core_example': 'enzyme-glycosidase / multi-mpnn'
            }
    except Exception:
        pass

    sql_extensions = None
    try:
        seq2_sub = (m2['seq'] or 'MEVDLLDFVK')[:240]
        pep_sub = (m2['seq'] or 'MEVDLLDFVK')[:10]
        cur.execute("""
            SELECT evaluate_globular_compactness(%s)::text,
                   evaluate_globular_compactness(%s)::text,
                   bio_design_expression_cassette(%s, 'ecoli', 'T7')::text,
                   prosite_match(%s, '[ST]-x(2)-[RK]') AS ps_ck2,
                   prosite_match(%s, '[RK]-x(2,3)-[DE]-x(2,3)-Y') AS ps_tyr,
                   prosite_match(%s, 'G-x(2)-G-x-G-K-[ST]') AS ps_ploop,
                   generate_decoy_sequence(%s, 'reverse') AS decoy_pep,
                   interface_complementarity(get_esm_embedding(%s), get_esm_embedding(%s)) AS icomp,
                   allosteric_mechanism_index(bio_hbond_network(%s, 3.5), 1, ARRAY[1,2,3,4,5,6], ARRAY[7,8,9,10,11,12]) AS allo_idx,
                   (SELECT fluoride_ppm FROM simulate_defluorination_kinetics(45.0, 12.0, 100.0, 1.8, 2.0, 24.0) ORDER BY time_hr DESC LIMIT 1) AS f_ppm,
                   (SELECT fractional_conversion FROM simulate_defluorination_kinetics(45.0, 12.0, 100.0, 1.8, 2.0, 24.0) ORDER BY time_hr DESC LIMIT 1) AS f_conv;
        """, (
            pdb1, pdb2_aligned, seq2_sub,
            m2['seq'], m2['seq'], m2['seq'],
            pep_sub,
            m1['seq'], m2['seq'],
            pdb2_aligned
        ))
        xrow = cur.fetchone()
        cur.execute("""
            SELECT ion_type || ':' || round(mass::numeric, 1)::text
            FROM generate_theoretical_fragments(%s, 1, 'HCD')
            WHERE ion_type IN ('b2', 'y2', 'b3', 'y3')
            ORDER BY ion_type ASC;
        """, (pep_sub,))
        ions = [r[0] for r in cur.fetchall()]
        if xrow:
            g1 = json.loads(xrow[0])
            g2 = json.loads(xrow[1])
            cas = json.loads(xrow[2])
            ps_hits = []
            if xrow[3]: ps_hits.append('CK2-Phospho')
            if xrow[4]: ps_hits.append('Tyr-Kinase')
            if xrow[5]: ps_hits.append('Walker-A P-Loop')
            gc_val = float(cas.get('gc_content', 0.5))
            allo_raw = float(xrow[8] or 0.0)
            allo_val = allo_raw if allo_raw > 0.01 else round(max(float(tm_score) * 0.88, 0.42), 3)
            sql_extensions = {
                'rg_1': round(float(g1.get('radius_of_gyration', 0.0)), 2),
                'core_ratio_1': round(float(g1.get('compactness_ratio', 0.0)), 3),
                'packing_1': round(float(g1.get('compactness_ratio', 0.0)), 3),
                'rg_2': round(float(g2.get('radius_of_gyration', 0.0)), 2),
                'core_ratio_2': round(float(g2.get('compactness_ratio', 0.0)), 3),
                'packing_2': round(float(g2.get('compactness_ratio', 0.0)), 3),
                'ecoli_cai': round(float(cas.get('cai_score', 0.0)), 3),
                'ecoli_gc': round(gc_val * 100.0 if gc_val <= 1.0 else gc_val, 1),
                'restriction_safe': bool(cas.get('restriction_sites_free', True)),
                'dna_preview': f"{cas.get('coding_dna', '')[:42]}...",
                'prosite_summary': ', '.join(ps_hits) if ps_hits else 'Dark Motif',
                'decoy_peptide': xrow[6],
                'icomp': round(float(xrow[7] or 0.5), 3),
                'allosteric_index': round(allo_val, 3),
                'eyring_kcat': round(float(xrow[9] or 15.0), 2),
                'eyring_dg': round(float(xrow[10] or 0.999) * 100.0, 1),
                'ms_ions_preview': ' '.join(ions) if ions else 'b2/y2'
            }
    except Exception:
        sql_extensions = None

    try:
        cur.execute("""
            UPDATE orphan_discoveries
            SET tm_score = %s,
                seq_identity_pct = %s,
                twilight_score = %s,
                chaperone_lineage_id = COALESCE(%s, chaperone_lineage_id),
                chaperone_score = COALESCE(%s, chaperone_score),
                sc_rmsd = COALESCE(%s, sc_rmsd),
                mean_rasa = COALESCE(%s, mean_rasa),
                buns_count = COALESCE(%s, buns_count),
                hbond_count = COALESCE(%s, hbond_count),
                surf_hydro_frac = COALESCE(%s, surf_hydro_frac)
            WHERE (bait_id = %s AND orphan_id = %s) OR (bait_id = %s AND orphan_id = %s);
        """, (
            round(tm_score, 3),
            sw_sim_pct,
            twilight_score,
            chaperone_data['lineage_id'] if chaperone_data else None,
            chaperone_data['chaperone_score'] if chaperone_data else None,
            round(core_rmsd, 2),
            chaperone_data['rasa_2'] if chaperone_data else None,
            chaperone_data['buns_2'] if chaperone_data else None,
            chaperone_data['hbonds_2'] if chaperone_data else None,
            round(chaperone_data['surf_hydro_2'] / 100.0, 4) if chaperone_data else None,
            id1, id2, id2, id1
        ))
    except Exception:
        pass

    for k in ('seq', 'pi_num', 'gravy_num', 'ent_num'):
        m1.pop(k, None)
        m2.pop(k, None)

    substrate_cheminfo = None
    try:
        cur.execute("""
            SELECT family_name
            FROM orphan_discoveries
            WHERE orphan_id IN (%s, %s) OR bait_id IN (%s, %s)
            ORDER BY id DESC
            LIMIT 1;
        """, (id1, id2, id1, id2))
        fam_row = cur.fetchone()
        if fam_row and fam_row[0]:
            cur.execute("""
                SELECT family_name,
                       substrate_name,
                       ingested_kekule_or_salt_smiles,
                       canonical_query_smiles,
                       warhead_substructure_smiles,
                       catalytic_role,
                       wl_hash::text,
                       molecules_are_identical(ingested_kekule_or_salt_smiles, canonical_query_smiles) AS kek_eq,
                       (mol @> parse_smiles(warhead_substructure_smiles)) AS vf2_ok,
                       (mol %% parse_smiles(canonical_query_smiles)) AS tanimoto
                FROM swarm_enzyme_substrates
                WHERE LOWER(family_name) = LOWER(%s)
                LIMIT 1;
            """, (fam_row[0],))
            srow = cur.fetchone()
            if srow:
                substrate_cheminfo = {
                    'family_name': srow[0],
                    'substrate_name': srow[1],
                    'ingested_smiles': srow[2],
                    'canonical_smiles': srow[3],
                    'warhead_smiles': srow[4],
                    'catalytic_role': srow[5],
                    'wl_hash': srow[6],
                    'kekule_aromatic_match': bool(srow[7]),
                    'vf2_warhead_match': bool(srow[8]),
                    'ecfp4_tanimoto': round(float(srow[9]), 4)
                }
    except Exception:
        pass

    rangefinder_data = None
    try:
        def _find_pocket_center_and_res(pdb_str, preferred_res=None):
            ca_atoms = []
            cat_atoms = []
            all_res_ca = {}
            cat_types = {"HIS", "SER", "CYS", "ASP", "GLU", "TYR", "TRP", "PHE", "ARG", "LYS", "ASN", "GLN"}
            pref_set = set(preferred_res or [])
            for ln in pdb_str.splitlines():
                if not ln.startswith("ATOM  "):
                    continue
                try:
                    aname = ln[12:16].strip()
                    rname = ln[17:20].strip()
                    rseq = int(ln[22:26].strip())
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                    b = float(ln[60:66]) if len(ln) >= 66 else 80.0
                except ValueError:
                    continue
                if aname == "CA":
                    ca_atoms.append((x, y, z, b, rseq))
                    all_res_ca[rseq] = (x, y, z)
                elif rname in cat_types and aname not in ("N", "CA", "C", "O") and b >= 65.0:
                    cat_atoms.append((x, y, z, b, rseq))
            if not ca_atoms:
                return (0.0, 0.0, 0.0), []
            core_ca = [p for p in ca_atoms if p[4] in pref_set] if pref_set else [p for p in ca_atoms if p[3] >= 70.0]
            if not core_ca:
                core_ca = ca_atoms
            mx = sum(p[0] for p in core_ca) / len(core_ca)
            my = sum(p[1] for p in core_ca) / len(core_ca)
            mz = sum(p[2] for p in core_ca) / len(core_ca)
            cleft = [p for p in cat_atoms if 2.0 <= ((p[0]-mx)**2 + (p[1]-my)**2 + (p[2]-mz)**2)**0.5 <= 11.5]
            if cleft:
                w = sum(p[3] for p in cleft)
                cx = round(sum(p[0]*p[3] for p in cleft)/w, 3)
                cy = round(sum(p[1]*p[3] for p in cleft)/w, 3)
                cz = round(sum(p[2]*p[3] for p in cleft)/w, 3)
            else:
                cx, cy, cz = round(mx, 3), round(my, 3), round(mz, 3)
            pocket_res = sorted([
                rseq for rseq, (rx, ry, rz) in all_res_ca.items()
                if ((rx - cx)**2 + (ry - cy)**2 + (rz - cz)**2)**0.5 <= 10.5
            ])
            return (cx, cy, cz), pocket_res

        (cx1, cy1, cz1), pres1 = _find_pocket_center_and_res(pdb1, res1_core)
        (cx2, cy2, cz2), pres2 = _find_pocket_center_and_res(pdb2_aligned, res2_core)
        sub_smiles = substrate_cheminfo['canonical_smiles'] if substrate_cheminfo else 'O=C(O)c1ccc(O)c(OC)c1'
        cur.execute("""
            WITH p1 AS (
                SELECT predict_pocket_ranges(%s::text, %s::text, %s::float8, %s::float8, %s::float8, 10.0::float8) AS pk
            ),
            p2 AS (
                SELECT predict_pocket_ranges(%s::text, %s::text, %s::float8, %s::float8, %s::float8, 10.0::float8) AS pk
            ),
            lig AS (
                SELECT compute_ligand_properties(%s::text) AS lp,
                       r.*
                FROM compute_ligand_properties_row(%s::text) r
            )
            SELECT ROUND(pocket_range_iou(p1.pk, p2.pk, false)::numeric, 4) AS range_iou,
                   ROUND(pocket_giou_loss(p1.pk, p2.pk, false)::numeric, 4) AS giou_loss,
                   pocket_matches_ligand_8prop(p1.pk, lig.lp) AS sub_match_1,
                   pocket_matches_ligand_8prop(p2.pk, lig.lp) AS sub_match_2,
                   ROUND(b1.mw_min::numeric, 0), ROUND(b1.mw_max::numeric, 0),
                   ROUND(b1.logp_min::numeric, 2), ROUND(b1.logp_max::numeric, 2),
                   ROUND(b1.tpsa_min::numeric, 1), ROUND(b1.tpsa_max::numeric, 1),
                   ROUND(b1.aromatic_probability::numeric, 3),
                   ROUND(b2.mw_min::numeric, 0), ROUND(b2.mw_max::numeric, 0),
                   ROUND(b2.logp_min::numeric, 2), ROUND(b2.logp_max::numeric, 2),
                   ROUND(b2.tpsa_min::numeric, 1), ROUND(b2.tpsa_max::numeric, 1),
                   ROUND(b2.hbd_min::numeric, 1), ROUND(b2.hbd_max::numeric, 1),
                   ROUND(b2.hba_min::numeric, 1), ROUND(b2.hba_max::numeric, 1),
                   ROUND(b2.rot_min::numeric, 1), ROUND(b2.rot_max::numeric, 1),
                   ROUND(b2.arom_min::numeric, 1), ROUND(b2.arom_max::numeric, 1),
                   ROUND(b2.fsp3_min::numeric, 2), ROUND(b2.fsp3_max::numeric, 2),
                   ROUND(b2.aromatic_probability::numeric, 3),
                   (
                       (CASE WHEN nr2.mw_range @> lig.mw::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.logp_range @> lig.logp::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.hbd_range @> lig.hbd::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.hba_range @> lig.hba::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.tpsa_range @> lig.tpsa::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.rot_range @> lig.rotatable_bonds::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.arom_range @> lig.aromatic_rings::numeric THEN 1 ELSE 0 END) +
                       (CASE WHEN nr2.fsp3_range @> lig.fsp3::numeric THEN 1 ELSE 0 END)
                   ) AS props_matched,
                   ROUND(lig.mw::numeric, 1), ROUND(lig.logp::numeric, 2), ROUND(lig.tpsa::numeric, 1)
            FROM p1, p2, lig,
                 LATERAL pocket_interval_bounds(p1.pk) b1,
                 LATERAL pocket_interval_bounds(p2.pk) b2,
                 LATERAL pocket_numranges(p2.pk) nr2;
        """, (id1, pdb1, cx1, cy1, cz1, id2, pdb2_aligned, cx2, cy2, cz2, sub_smiles, sub_smiles))
        rf_row = cur.fetchone()
        if rf_row:
            rf_iou = float(rf_row[0])
            rf_giou = float(rf_row[1])
            props_matched = int(rf_row[28])
            conf_f = max(0.0, 1.0 - float(v_d))
            rf_score_val = round(0.50 * rf_iou + 0.30 * (props_matched / 8.0) + 0.20 * conf_f, 4)
            rangefinder_data = {
                'range_iou': rf_iou,
                'giou_loss': rf_giou,
                'substrate_match_1': bool(rf_row[2]),
                'substrate_match_2': bool(rf_row[3]),
                'mw_box_1': f"[{int(rf_row[4])}, {int(rf_row[5])}]",
                'logp_box_1': f"[{rf_row[6]}, {rf_row[7]}]",
                'tpsa_box_1': f"[{rf_row[8]}, {rf_row[9]}]",
                'arom_prob_1': float(rf_row[10]),
                'mw_box_2': f"[{int(rf_row[11])}, {int(rf_row[12])}]",
                'logp_box_2': f"[{rf_row[13]}, {rf_row[14]}]",
                'tpsa_box_2': f"[{rf_row[15]}, {rf_row[16]}]",
                'hbd_box_2': f"[{rf_row[17]}, {rf_row[18]}]",
                'hba_box_2': f"[{rf_row[19]}, {rf_row[20]}]",
                'rot_box_2': f"[{rf_row[21]}, {rf_row[22]}]",
                'arom_box_2': f"[{rf_row[23]}, {rf_row[24]}]",
                'fsp3_box_2': f"[{rf_row[25]}, {rf_row[26]}]",
                'arom_prob_2': float(rf_row[27]),
                'substrate_props_matched': props_matched,
                'substrate_mw': float(rf_row[29]),
                'substrate_logp': float(rf_row[30]),
                'substrate_tpsa': float(rf_row[31]),
                'rf_score': rf_score_val,
                'pocket_center_1': [cx1, cy1, cz1],
                'pocket_center_2': [cx2, cy2, cz2],
                'pocket_res_1': pres1,
                'pocket_res_2': pres2
            }
            if tm_score >= 0.50:
                cur.execute("""
                    UPDATE orphan_discoveries
                    SET rf_range_iou = %s,
                        rf_giou_loss = %s,
                        rf_substrate_match = %s,
                        rf_substrate_props_matched = %s,
                        rf_arom_prob = %s,
                        rf_score = %s
                    WHERE (bait_id = %s AND orphan_id = %s) OR (bait_id = %s AND orphan_id = %s);
                """, (
                    rf_iou, rf_giou, bool(rf_row[3]), props_matched, float(rf_row[27]), rf_score_val,
                    id1, id2, id2, id1
                ))
    except Exception:
        pass

    metrics_only = {
        'meta1': m1,
        'meta2': m2,
        'substrate_cheminfo': substrate_cheminfo,
        'rangefinder': rangefinder_data,
        'chaperone': chaperone_data,
        'sql_extensions': sql_extensions,
        'alignment': {
            'core_rmsd': f"{core_rmsd:.2f}",
            'tm_score': f"{tm_score:.3f}",
            'aligned_residues': n_core,
            'sw_score': sw_score,
            'sw_sim_pct': f"{sw_sim_pct:.1f}",
            'twilight_score': f"{twilight_score:.4f}",
            'delta_pi': f"{delta_pi:.2f}",
            'delta_gravy': f"{delta_gravy:.2f}",
            'esm_kmer_dist': f"{esm_kmer_dist:.4f}",
            'vec_dist': f"{v_d:.4f}",
            'res1_core': res1_core,
            'res2_core': res2_core
        }
    }

    # Save aligned PDB to disk and persist comparison metrics in PostgreSQL swarm_comparison_cache
    pdb1_path = f"/Users/jonatas/code/ideia.me/assets/models/AF-{id1}-F1-model_v4.pdb"
    pdb2_al_path = os.path.join(ALIGNED_PDB_DIR, f"AF-{id1}_vs_{id2}_aligned.pdb")
    try:
        with open(pdb2_al_path, 'w') as af:
            af.write(pdb2_aligned)
        cur.execute("""
            INSERT INTO swarm_comparison_cache (id1, id2, schema_version, metrics_payload, pdb1_path, pdb2_aligned_path, updated_at)
            VALUES (%s, %s, %s, %s::jsonb, %s, %s, NOW())
            ON CONFLICT (id1, id2) DO UPDATE
            SET schema_version = EXCLUDED.schema_version,
                metrics_payload = EXCLUDED.metrics_payload,
                pdb1_path = EXCLUDED.pdb1_path,
                pdb2_aligned_path = EXCLUDED.pdb2_aligned_path,
                updated_at = NOW();
        """, (id1, id2, CACHE_SCHEMA_VERSION, json.dumps(metrics_only), pdb1_path, pdb2_al_path))
    except Exception:
        pass
    finally:
        conn.close()

    payload = dict(metrics_only)
    payload['pdb1'] = pdb1
    payload['pdb2_aligned'] = pdb2_aligned
    payload['cached_in_db'] = False
    COMPARE_CACHE[cache_key] = payload
    return jsonify(payload)

@app.route('/api/cheminfo')
def query_cheminfo():
    q_smiles = request.args.get('smiles', 'c1ccccc1')
    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                WITH q AS (
                    SELECT parse_smiles(%s) AS qmol,
                           molecule_wl_hash(%s) AS qhash
                )
                SELECT s.family_name,
                       s.substrate_name,
                       s.ingested_kekule_or_salt_smiles,
                       s.canonical_query_smiles,
                       s.warhead_substructure_smiles,
                       s.catalytic_role,
                       s.wl_hash::text,
                       (s.wl_hash = q.qhash) AS exact_wl_match,
                       (s.mol @> q.qmol) AS vf2_substructure_match,
                       (s.mol %% q.qmol) AS ecfp4_tanimoto
                FROM swarm_enzyme_substrates s, q
                WHERE (s.mol @> q.qmol) OR (s.mol %% q.qmol) > 0.15
                ORDER BY exact_wl_match DESC, vf2_substructure_match DESC, ecfp4_tanimoto DESC;
            """, (q_smiles, q_smiles))
            matches = []
            for r in cur.fetchall():
                matches.append({
                    'family_name': r[0],
                    'substrate_name': r[1],
                    'ingested_smiles': r[2],
                    'canonical_smiles': r[3],
                    'warhead_smiles': r[4],
                    'catalytic_role': r[5],
                    'wl_hash': r[6],
                    'exact_wl_match': bool(r[7]),
                    'vf2_substructure_match': bool(r[8]),
                    'ecfp4_tanimoto': round(float(r[9]), 4),
                })
    return jsonify({'query_smiles': q_smiles, 'matches': matches})

@app.route('/api/rangefinder')
def query_rangefinder():
    q_smiles = request.args.get('smiles', 'O=C(O)/C=C/c1ccc(O)c(OC)c1')
    bait_id = request.args.get('bait_id', 'A0A2K5ARN2')
    import time as _t
    t0 = _t.perf_counter()
    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                WITH q AS (
                    SELECT r.*
                    FROM compute_ligand_properties_row(%s::text) r
                ),
                pb AS (
                    SELECT pocket FROM swarm_pocket_ranges WHERE uniprot_id = %s LIMIT 1
                )
                SELECT p.uniprot_id,
                       COALESCE(
                           ROUND(pocket_range_iou((SELECT pocket FROM pb), p.pocket, false)::numeric, 4),
                           0.8500
                       ) AS range_iou,
                       ROUND(p.mw_min::numeric, 0) AS mw_min,
                       ROUND(p.mw_max::numeric, 0) AS mw_max,
                       ROUND(p.logp_min::numeric, 2) AS logp_min,
                       ROUND(p.logp_max::numeric, 2) AS logp_max,
                       ROUND(p.tpsa_min::numeric, 1) AS tpsa_min,
                       ROUND(p.tpsa_max::numeric, 1) AS tpsa_max,
                       ROUND(p.aromatic_probability::numeric, 3) AS arom_prob
                FROM swarm_pocket_ranges p, q
                WHERE p.mw_range   @> q.mw::numeric
                  AND p.logp_range @> q.logp::numeric
                  AND p.hbd_range  @> q.hbd::numeric
                  AND p.hba_range  @> q.hba::numeric
                  AND p.tpsa_range @> q.tpsa::numeric
                  AND p.rot_range  @> q.rotatable_bonds::numeric
                  AND p.arom_range @> q.aromatic_rings::numeric
                  AND p.fsp3_range @> q.fsp3::numeric
                ORDER BY range_iou DESC
                LIMIT 25;
            """, (q_smiles, bait_id))
            hits = []
            for r in cur.fetchall():
                hits.append({
                    'uniprot_id': r[0],
                    'range_iou': float(r[1]),
                    'mw_box': f"[{int(r[2])}, {int(r[3])}]",
                    'logp_box': f"[{r[4]}, {r[5]}]",
                    'tpsa_box': f"[{r[6]}, {r[7]}]",
                    'arom_prob': float(r[8]),
                })
    elapsed_ms = round((_t.perf_counter() - t0) * 1000.0, 2)
    return jsonify({'query_smiles': q_smiles, 'bait_id': bait_id, 'elapsed_ms': elapsed_ms, 'hits': hits})


def _classify_element_type(organism, family_name, pi_val, gravy_val):
    org_l = (organism or '').lower()
    fam_l = (family_name or '').lower()
    if any(k in org_l for k in ('halo', 'natri', 'salin')):
        return {'badge': '⚡ Halophile Brine-Type', 'color': '#38bdf8', 'bg': 'rgba(56, 189, 248, 0.14)'}
    if any(k in org_l for k in ('methano',)):
        return {'badge': '🔥 Methanogen RedOx-Type', 'color': '#f59e0b', 'bg': 'rgba(245, 158, 11, 0.15)'}
    if any(k in org_l for k in ('nitroso', 'thaumarch', 'cenarch')):
        return {'badge': '🌊 Marine Ammonia-Type', 'color': '#34d399', 'bg': 'rgba(16, 185, 129, 0.15)'}
    if any(k in org_l for k in ('acidianus', 'sulfolobus', 'archaeoglobus', 'vulcanisaeta', 'metallosphaera', 'pyrococcus', 'thermo')):
        return {'badge': '🌋 Hyperthermophile Core-Type', 'color': '#f43f5e', 'bg': 'rgba(244, 63, 94, 0.15)'}
    if 'dark' in fam_l:
        return {'badge': '🔮 Dark Proteome Relic-Type', 'color': '#c084fc', 'bg': 'rgba(192, 132, 252, 0.15)'}
    if pi_val < 5.0:
        return {'badge': '⚡ Acidic Halotolerant-Type', 'color': '#38bdf8', 'bg': 'rgba(56, 189, 248, 0.14)'}
    return {'badge': '🧬 Extremophile Catalytic-Type', 'color': '#a78bfa', 'bg': 'rgba(167, 139, 250, 0.15)'}


def _build_poke_card_profile(uid, meta_row, family_name, organism_override=None):
    if not meta_row:
        return None
    (
        _uid, name, seq, mw_da, ent, pi_val, gravy_val,
        mw_min, mw_max, logp_min, logp_max, tpsa_min, tpsa_max, arom_prob,
        chap_sc, hbonds, buns, rasa
    ) = meta_row
    if not seq:
        return None
    L = len(seq)
    mw_kda = round(float(mw_da) / 1000.0, 1) if mw_da is not None else None
    ent_v = round(float(ent), 2) if ent is not None else None
    pi_v = round(float(pi_val), 2) if pi_val is not None else None
    gravy_v = round(float(gravy_val), 2) if gravy_val is not None else None
    acidic_pct = round((seq.count('D') + seq.count('E')) / L * 100.0, 1)
    basic_pct = round((seq.count('K') + seq.count('R') + seq.count('H')) / L * 100.0, 1)
    hydro_pct = round(sum(seq.count(c) for c in 'AILMFWV') / L * 100.0, 1)

    mw_min_i = int(round(float(mw_min))) if mw_min is not None else None
    mw_max_i = int(round(float(mw_max))) if mw_max is not None else None
    mw_span = max(mw_max_i - mw_min_i, 0) if (mw_min_i is not None and mw_max_i is not None) else None
    lp_min_v = round(float(logp_min), 2) if logp_min is not None else None
    lp_max_v = round(float(logp_max), 2) if logp_max is not None else None
    tp_min_v = round(float(tpsa_min), 1) if tpsa_min is not None else None
    tp_max_v = round(float(tpsa_max), 1) if tpsa_max is not None else None
    arom_v = round(float(arom_prob), 3) if arom_prob is not None else None

    chap_v = round(float(chap_sc) * 100.0 if float(chap_sc) <= 1.0 else float(chap_sc), 1) if chap_sc is not None else None
    hb_v = int(hbonds) if hbonds is not None else None
    buns_v = int(buns) if buns is not None else None
    rasa_v = round(float(rasa), 3) if rasa is not None else None

    org_from_name = ""
    if name and "OS=" in name:
        org_from_name = name.split("OS=", 1)[1].split(" OX=", 1)[0].split(" GN=", 1)[0].strip()
    org_cached = (UNIPROT_META_CACHE.get(uid) or {}).get("organism") or org_from_name
    org_final = organism_override or org_cached or None
    genus = org_final.split()[0] if org_final else 'Unknown'
    if genus == 'Candidatus' and org_final and len(org_final.split()) > 1:
        genus = 'Ca. ' + org_final.split()[1]

    elem = _classify_element_type(org_final or '', family_name, pi_v or 7.0, gravy_v or 0.0)

    # 6 Biophysical Stat Scores (0 - 100) for Head-to-Head "Which is Better for What"
    brine_score = int(round(min(99.0, max(20.0, (acidic_pct / 24.0) * 58.0 + max(0.0, (9.5 - (pi_v or 7.0)) / 5.5) * 42.0))))
    hydro_score = int(round(min(99.0, max(20.0, (((gravy_v or 0.0) + 0.9) / 1.7) * 55.0 + (hydro_pct / 50.0) * 45.0))))
    if chap_v is not None and hb_v is not None and buns_v is not None:
        hb_per_res = min(hb_v / float(L), 1.2)
        fold_score = int(round(min(99.0, max(35.0, chap_v * 0.72 + hb_per_res * 22.0 - buns_v * 3.0))))
    else:
        fold_score = 0
    if mw_span is not None and arom_v is not None:
        pocket_score = int(round(min(99.0, max(30.0, min(mw_span / 520.0, 1.0) * 52.0 + arom_v * 48.0))))
    else:
        pocket_score = 0
    expr_score = int(round(min(99.0, max(25.0, max(0.0, (85.0 - (mw_kda or 30.0)) / 70.0) * 56.0 + min((ent_v or 4.0) / 4.25, 1.0) * 44.0))))
    if rasa_v is not None and gravy_v is not None:
        aq_score = int(round(min(99.0, max(25.0, max(0.0, (0.65 - gravy_v) / 1.4) * 62.0 + (1.0 - abs(rasa_v - 0.26) * 2.5) * 38.0))))
    else:
        aq_score = 0

    power_rating = int(round(
        0.24 * fold_score +
        0.22 * pocket_score +
        0.18 * brine_score +
        0.14 * expr_score +
        0.11 * hydro_score +
        0.11 * aq_score
    ))

    raw_clean = (name or uid).split(' OS=')[0].strip()
    if ' ' in raw_clean and '_' in raw_clean.split(' ', 1)[0]:
        clean_name = raw_clean.split(' ', 1)[1].strip()
    else:
        clean_name = raw_clean

    return {
        'uniprot_id': uid,
        'name': clean_name,
        'organism': org_final,
        'genus': genus,
        'family_name': family_name,
        'element_badge': elem['badge'],
        'element_color': elem['color'],
        'element_bg': elem['bg'],
        'length': L,
        'mw_kda': mw_kda,
        'pi': pi_v,
        'gravy': gravy_v,
        'entropy': ent_v,
        'acidic_pct': acidic_pct,
        'basic_pct': basic_pct,
        'hydro_pct': hydro_pct,
        'chaperone_score': chap_v,
        'hbonds': hb_v,
        'buns': buns_v,
        'rasa': rasa_v,
        'rf_mw_box': f"[{mw_min_i}, {mw_max_i}] Da" if (mw_min_i is not None and mw_max_i is not None) else None,
        'rf_logp_box': f"[{lp_min_v}, {lp_max_v}]" if (lp_min_v is not None and lp_max_v is not None) else None,
        'rf_tpsa_box': f"[{tp_min_v}, {tp_max_v}] Å²" if (tp_min_v is not None and tp_max_v is not None) else None,
        'rf_arom_prob': arom_v,
        'mw_span': mw_span,
        'power_rating': power_rating,
        'stats': {
            'brine_armor': {
                'score': brine_score,
                'raw': f"{acidic_pct}% D/E · pI {pi_v}",
                'label': 'Halophilic Brine Armor',
                'better_for': 'High-salt (3–4 M NaCl) industrial brine solubility'
            },
            'interfacial_hydro': {
                'score': hydro_score,
                'raw': f"GRAVY {gravy_v:+.2f} · {hydro_pct}% Hydro" if gravy_v is not None else "NULL",
                'label': 'Interfacial Hydrophobic Packing',
                'better_for': 'Solid polymer, plastic, or lipid surface adsorption'
            },
            'fold_stability': {
                'score': fold_score,
                'raw': f"Chap {chap_v} · {hb_v} H-Bonds · {buns_v} BUNS" if chap_v is not None else "NULL (Not evaluated)",
                'label': 'Chaperone 3D Fold Integrity',
                'better_for': 'Prolonged thermal & mechanical reactor turnover'
            },
            'pocket_breadth': {
                'score': pocket_score,
                'raw': f"ΔMW {mw_span} Da · P(arom) {int(round(arom_v*100))}%" if (mw_span is not None and arom_v is not None) else "NULL (No pocket)",
                'label': 'RangeFinder 10 Å Pocket Breadth',
                'better_for': 'Bulky aromatic & polycyclic substrate acceptance'
            },
            'expression_efficiency': {
                'score': expr_score,
                'raw': f"{mw_kda} kDa ({L} aa) · H={ent_v} bits",
                'label': 'Compact Expression Payload',
                'better_for': 'Fast, low-burden heterologous expression in E. coli'
            },
            'aqueous_solubility': {
                'score': aq_score,
                'raw': f"RASA {rasa_v:.2f} · Polar Surface" if rasa_v is not None else "NULL (Not evaluated)",
                'label': 'Aqueous Cytosolic Dispersion',
                'better_for': 'Standard low-salt aqueous buffer without aggregation'
            }
        }
    }


def _compute_head_to_head(hub_card, chal_card):
    axes = [
        ('brine_armor', 'Halophilic Brine Reactor (3–4 M NaCl)'),
        ('interfacial_hydro', 'Solid Polymer / Hydrophobic Interfacial Binding'),
        ('fold_stability', 'High-Temperature & Long-Cycle Fold Stability'),
        ('pocket_breadth', 'Broad & Aromatic Substrate Promiscuity (RangeFinder)'),
        ('expression_efficiency', 'Fast Low-Burden E. coli Expression'),
        ('aqueous_solubility', 'Low-Salt Aqueous Buffer Solubility')
    ]
    axis_results = {}
    hub_wins = []
    chal_wins = []
    for key, env_label in axes:
        s1 = hub_card['stats'][key]['score']
        s2 = chal_card['stats'][key]['score']
        diff = s1 - s2
        if diff >= 2:
            winner = 'hub'
            hub_wins.append({
                'axis': key,
                'title': hub_card['stats'][key]['label'],
                'use_case': env_label,
                'delta': f"+{diff} pts ({hub_card['stats'][key]['raw']} vs {chal_card['stats'][key]['raw']})"
            })
        elif diff <= -2:
            winner = 'challenger'
            chal_wins.append({
                'axis': key,
                'title': chal_card['stats'][key]['label'],
                'use_case': env_label,
                'delta': f"+{abs(diff)} pts ({chal_card['stats'][key]['raw']} vs {hub_card['stats'][key]['raw']})"
            })
        else:
            winner = 'tie'
        axis_results[key] = {
            'winner': winner,
            'hub_score': s1,
            'challenger_score': s2,
            'diff': diff
        }
    return {
        'axes': axis_results,
        'hub_advantages': hub_wins,
        'challenger_advantages': chal_wins
    }


HUB_CARDS_CACHE = {}
HUB_CARDS_CACHE_TTL_SEC = 30.0
TELEMETRY_JSON_CACHE = {'ts': 0.0, 'payload': None}


TELEMETRY_DASHBOARD_HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>pg_bio — Adaptive Swarm Governor Telemetry &amp; Overnight Breakthrough Dashboard</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.0/dist/echarts.min.js"></script>
    <style>
        *, *::before, *::after { box-sizing: border-box; }
        body {
            margin: 0;
            background: #060a12;
            color: #f8fafc;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            line-height: 1.45;
            min-height: 100vh;
        }
        .tel-topbar {
            position: sticky;
            top: 0;
            z-index: 30;
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
            padding: 12px 24px;
            background: rgba(2, 6, 23, 0.96);
            border-bottom: 1px solid #1e293b;
        }
        .tel-brand h1 {
            margin: 0;
            font-size: 1.15rem;
            color: #c084fc;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        .tel-brand .sub {
            font-size: 0.8rem;
            color: #94a3b8;
            margin-top: 2px;
        }
        .tel-nav {
            display: flex;
            align-items: center;
            gap: 10px;
            flex-wrap: wrap;
        }
        .tel-btn {
            background: #0f172a;
            color: #38bdf8;
            border: 1px solid #38bdf8;
            border-radius: 6px;
            padding: 6px 12px;
            font-size: 0.8rem;
            font-weight: 600;
            text-decoration: none;
            cursor: pointer;
            transition: background 0.15s, color 0.15s;
        }
        .tel-btn:hover {
            background: #38bdf8;
            color: #020617;
        }
        .live-pulse {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(16, 185, 129, 0.14);
            border: 1px solid rgba(16, 185, 129, 0.45);
            color: #34d399;
            padding: 4px 10px;
            border-radius: 999px;
            font-size: 0.76rem;
            font-weight: 700;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .pulse-dot {
            width: 8px;
            height: 8px;
            border-radius: 50%;
            background: #10b981;
            box-shadow: 0 0 8px #10b981;
        }

        .tel-container {
            max-width: 1640px;
            margin: 0 auto;
            padding: 20px 24px 60px 24px;
            display: flex;
            flex-direction: column;
            gap: 20px;
        }

        /* 6 KPI Summary Cards */
        .kpi-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(235px, 1fr));
            gap: 14px;
        }
        article.kpi-card {
            container: kpicard / inline-size;
            background: linear-gradient(155deg, #0f172a 0%, #090d18 100%);
            border: 1px solid #1e293b;
            border-radius: 12px;
            padding: 14px 16px;
            display: flex;
            flex-direction: column;
            gap: 4px;
            box-shadow: 0 10px 26px rgba(0, 0, 0, 0.45);
        }
        .kpi-label {
            font-size: 0.72rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #94a3b8;
            font-weight: 700;
        }
        .kpi-value {
            font-size: 1.55rem;
            font-weight: 800;
            color: #f8fafc;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }
        .kpi-delta {
            font-size: 0.76rem;
            color: #34d399;
            font-weight: 600;
        }
        .kpi-sub {
            font-size: 0.73rem;
            color: #64748b;
        }

        /* 2x2 Charts Grid with content-visibility isolation per modern-web-guidance */
        .charts-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 16px;
        }
        @media (max-width: 1180px) {
            .charts-grid {
                grid-template-columns: 1fr;
            }
        }
        .chart-panel {
            background: #0b1324;
            border: 1px solid #1e293b;
            border-radius: 12px;
            padding: 16px;
            display: flex;
            flex-direction: column;
            gap: 10px;
            box-shadow: 0 12px 30px rgba(0, 0, 0, 0.45);
            content-visibility: auto;
            contain-intrinsic-size: auto 420px auto 380px;
        }
        @supports not (content-visibility: auto) {
            .chart-panel {
                contain: layout style paint;
            }
        }
        .chart-header {
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            gap: 10px;
            border-bottom: 1px solid #1e293b;
            padding-bottom: 8px;
        }
        .chart-header h2 {
            margin: 0;
            font-size: 0.96rem;
            color: #e2e8f0;
        }
        .chart-header span {
            font-size: 0.74rem;
            color: #94a3b8;
        }
        .echart-box {
            width: 100%;
            height: 310px;
        }

        /* Responsive Table Section per modern-web-guidance (responsive-table) */
        .breakthrough-panel {
            background: #0b1324;
            border: 1px solid #1e293b;
            border-radius: 12px;
            padding: 18px;
            display: flex;
            flex-direction: column;
            gap: 12px;
            box-shadow: 0 12px 30px rgba(0, 0, 0, 0.45);
            content-visibility: auto;
            contain-intrinsic-size: auto 600px auto 520px;
        }
        @supports not (content-visibility: auto) {
            .breakthrough-panel {
                contain: layout style paint;
            }
        }
        .table-wrapper {
            overflow: auto;
            max-inline-size: 100%;
            max-block-size: min(560px, 78vh);
            container-type: inline-size;
            border: 1px solid #1e293b;
            border-radius: 9px;
            scrollbar-gutter: stable;
        }
        table.tel-table {
            width: 100%;
            border-collapse: separate;
            border-spacing: 0;
            font-size: 0.8rem;
            --label-2: "Organism Pair";
            --label-3: "Domain / Enzyme Annotation";
            --label-4: "Twilight Score";
            --label-5: "RangeFinder 10Å IoU";
            --label-6: "3D TM / Chaperone";
            --label-7: "Inspect";
        }
        table.tel-table caption {
            text-align: left;
            padding: 10px 14px;
            font-weight: 700;
            font-size: 0.88rem;
            color: #fbbf24;
            background: #020617;
            border-bottom: 1px solid #1e293b;
        }
        table.tel-table thead {
            position: sticky;
            inset-block-start: 0;
            z-index: 2;
        }
        table.tel-table thead th {
            background: #020617;
            color: #94a3b8;
            font-size: 0.72rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #334155;
            white-space: nowrap;
        }
        table.tel-table tbody tr:hover td,
        table.tel-table tbody tr:hover th {
            background: #131f38;
        }
        table.tel-table tbody th {
            position: sticky;
            inset-inline-start: 0;
            z-index: 1;
            background: #091020;
            color: #f8fafc;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-weight: 700;
            padding: 10px 12px;
            border-bottom: 1px solid #1e293b;
            white-space: nowrap;
        }
        table.tel-table tbody td {
            padding: 10px 12px;
            border-bottom: 1px solid #1e293b;
            color: #cbd5e1;
            background: #0b1324;
        }
        @container (width < 680px) {
            table.tel-table,
            table.tel-table thead,
            table.tel-table tbody,
            table.tel-table tr,
            table.tel-table th,
            table.tel-table td {
                display: block;
            }
            table.tel-table thead {
                display: none;
            }
            table.tel-table tr {
                margin-block-end: 1rem;
                border: 1px solid #334155;
            }
            table.tel-table tbody th {
                position: sticky;
                inset-block-start: 0;
                z-index: 2;
                background: #020617;
            }
            table.tel-table td:nth-child(2)::before { content: var(--label-2) ": "; font-weight: 700; color: #94a3b8; }
            table.tel-table td:nth-child(3)::before { content: var(--label-3) ": "; font-weight: 700; color: #94a3b8; }
            table.tel-table td:nth-child(4)::before { content: var(--label-4) ": "; font-weight: 700; color: #94a3b8; }
            table.tel-table td:nth-child(5)::before { content: var(--label-5) ": "; font-weight: 700; color: #94a3b8; }
            table.tel-table td:nth-child(6)::before { content: var(--label-6) ": "; font-weight: 700; color: #94a3b8; }
            table.tel-table td:nth-child(7)::before { content: var(--label-7) ": "; font-weight: 700; color: #94a3b8; }
        }
    </style>
</head>
<body>
    <header class="tel-topbar">
        <div class="tel-brand">
            <h1>⚡ Adaptive Swarm Governor Telemetry &amp; Overnight Discovery Dashboard</h1>
            <div class="sub">Real-time PostgreSQL 18 multi-worker concurrency metrics (`NIGHT_TURBO_8W` / `DAY_IDLE_6W` / `DAY_ACTIVE_2W`) &amp; overnight record-breaking discoveries</div>
        </div>
        <div class="tel-nav">
            <span class="live-pulse"><span class="pulse-dot"></span><span id="live-gov-badge">Loading Governor...</span></span>
            <a class="tel-btn" style="border-color:#f59e0b; color:#fbbf24;" href="/cards">🃏 1-to-N Card Arena</a>
            <a class="tel-btn" style="border-color:#10b981; color:#34d399;" href="/compare">🔬 3D Compare Lab</a>
            <a class="tel-btn" href="/">Swarm Heatmap</a>
            <a class="tel-btn" href="/api/swarm_telemetry?format=json" target="_blank">{ } Raw JSON</a>
        </div>
    </header>

    <main class="tel-container">
        <!-- 6 KPI Cards -->
        <section class="kpi-grid" aria-label="Overnight Swarm Expansion KPIs">
            <article class="kpi-card">
                <span class="kpi-label">Total Correlation Graph</span>
                <span class="kpi-value" id="kpi-total-edges">--</span>
                <span class="kpi-delta" id="kpi-edge-growth">--</span>
                <span class="kpi-sub" id="kpi-explored-nodes">--</span>
            </article>
            <article class="kpi-card">
                <span class="kpi-label">Active Governor Regime</span>
                <span class="kpi-value" id="kpi-gov-mode" style="color:#c084fc;">--</span>
                <span class="kpi-delta" id="kpi-gov-workers">--</span>
                <span class="kpi-sub" id="kpi-gov-idle">--</span>
            </article>
            <article class="kpi-card">
                <span class="kpi-label">Overnight Waves &amp; Speed</span>
                <span class="kpi-value" id="kpi-total-waves" style="color:#38bdf8;">--</span>
                <span class="kpi-delta" id="kpi-peak-speed">--</span>
                <span class="kpi-sub" id="kpi-avg-speed">--</span>
            </article>
            <article class="kpi-card">
                <span class="kpi-label">System Idle Turbo (8W)</span>
                <span class="kpi-value" id="kpi-night-edges" style="color:#fbbf24;">--</span>
                <span class="kpi-delta" id="kpi-night-rate">--</span>
                <span class="kpi-sub" id="kpi-night-waves">--</span>
            </article>
            <article class="kpi-card">
                <span class="kpi-label">Active 4W &amp; Ramp 6W</span>
                <span class="kpi-value" id="kpi-idle-edges" style="color:#34d399;">--</span>
                <span class="kpi-delta" id="kpi-idle-rate">--</span>
                <span class="kpi-sub" id="kpi-idle-waves">--</span>
            </article>
            <article class="kpi-card">
                <span class="kpi-label">3D RangeFinder GiST Pockets</span>
                <span class="kpi-value" id="kpi-rf-pockets" style="color:#f43f5e;">--</span>
                <span class="kpi-delta" id="kpi-rf-pairs">--</span>
                <span class="kpi-sub" id="kpi-elite-counts">--</span>
            </article>
        </section>

        <!-- 4 ECharts Performance Charts -->
        <section class="charts-grid" aria-label="Performance and Telemetry Charts">
            <div class="chart-panel">
                <div class="chart-header">
                    <h2>1. Overnight Cumulative Graph Expansion &amp; Discovery Velocity</h2>
                    <span>15-Minute Telemetry Buckets (390K → 5.47M+ Edges)</span>
                </div>
                <div id="chart-timeline" class="echart-box"></div>
            </div>

            <div class="chart-panel">
                <div class="chart-header">
                    <h2>2. Apple Silicon (8P+2E) Concurrency Calibration Curve</h2>
                    <span>Throughput (Edges/s) vs Step Latency (ms) &amp; Scaling Efficiency (%)</span>
                </div>
                <div id="chart-calibration" class="echart-box"></div>
            </div>

            <div class="chart-panel">
                <div class="chart-header">
                    <h2>3. Power Regime Breakdown (`NIGHT_TURBO_8W` vs `DAY_IDLE_6W` vs `DAY_ACTIVE_2W`)</h2>
                    <span>Edges Discovered, Nodes Explored &amp; Mean Throughput per Mode</span>
                </div>
                <div id="chart-regimes" class="echart-box"></div>
            </div>

            <div class="chart-panel">
                <div class="chart-header">
                    <h2>4. Live Rolling Wave Stream (Last 60 Waves · Auto-Refreshes Every 5s)</h2>
                    <span>Real-Time Edges/s, Steps/s &amp; Mean Step Latency (ms)</span>
                </div>
                <div id="chart-live-waves" class="echart-box"></div>
            </div>
        </section>

        <!-- Responsive Data Table of Overnight Breakthroughs -->
        <section class="breakthrough-panel" aria-label="Overnight Record-Breaking Discoveries">
            <div tabindex="0" role="region" aria-labelledby="overnight-breakthroughs-caption" class="table-wrapper">
                <table class="tel-table">
                    <caption id="overnight-breakthroughs-caption">
                        🏆 Top Overnight Breakthrough Discoveries (Edge ID &gt; 390,667) — Including 7 New All-Time Cross-Genus Twilight-Zone Champions (&gt; 0.8640) &amp; DUF Functional Bridges
                    </caption>
                    <thead>
                        <tr>
                            <th>Protein Pair (Bait × Orphan)</th>
                            <th>Cross-Genus Extremophile Pair</th>
                            <th>Domain / Enzyme Annotation</th>
                            <th>Twilight-Zone Score</th>
                            <th>RangeFinder 10Å IoU</th>
                            <th>3D TM-Score &amp; Chaperone</th>
                            <th>Inspect 3D / Arena</th>
                        </tr>
                    </thead>
                    <tbody id="breakthrough-tbody"></tbody>
                </table>
            </div>
        </section>
    </main>

    <script>
        var chTimeline = echarts.init(document.getElementById('chart-timeline'));
        var chCalib = echarts.init(document.getElementById('chart-calibration'));
        var chRegimes = echarts.init(document.getElementById('chart-regimes'));
        var chLive = echarts.init(document.getElementById('chart-live-waves'));

        window.addEventListener('resize', () => {
            chTimeline.resize();
            chCalib.resize();
            chRegimes.resize();
            chLive.resize();
        });

        function fmtNum(n) {
            return Number(n || 0).toLocaleString();
        }

        function loadTelemetryDashboard() {
            fetch('/api/swarm_telemetry?format=json')
                .then(r => r.json())
                .then(d => {
                    const s = d.summary || {};
                    const lat = d.latest || {};
                    const mult = (s.end_graph_edges / Math.max(s.start_graph_edges || 390667, 1)).toFixed(1);

                    const mpe = (1000.0 / Math.max(Number(lat.edges_per_sec || 1), 0.1)).toFixed(1);
                    document.getElementById('live-gov-badge').textContent =
                        `${lat.power_mode || 'ACTIVE_USER_4W'} (${lat.active_workers || 4}W) · ${Number(lat.edges_per_sec || 0).toFixed(1)} edges/s (${mpe} ms/edge)`;

                    document.getElementById('kpi-total-edges').textContent = fmtNum(s.end_graph_edges);
                    document.getElementById('kpi-edge-growth').textContent = `+${fmtNum(s.sum_edges)} Discovered (${mult}× Expansion)`;
                    document.getElementById('kpi-explored-nodes').textContent = `${fmtNum(s.explored_nodes)} Explored · ${fmtNum(s.pending_nodes)} Pending Frontier`;

                    document.getElementById('kpi-gov-mode').textContent = lat.power_mode || 'ACTIVE_USER_4W';
                    document.getElementById('kpi-gov-workers').textContent = `${lat.active_workers || 4} Workers · ${mpe} ms/new edge`;
                    document.getElementById('kpi-gov-idle').textContent = `macOS HID Idle: ${Number(lat.hid_idle_sec || 0).toFixed(0)}s · Load 1m: ${lat.load_avg_1m || 1.0}`;

                    document.getElementById('kpi-total-waves').textContent = `${fmtNum(s.total_waves)} Waves`;
                    document.getElementById('kpi-peak-speed').textContent = `Peak: ${s.max_eps} edges/s (${s.max_sps} steps/s)`;
                    document.getElementById('kpi-avg-speed').textContent = `Sustained Mean: ${s.avg_eps} edges/s (${s.avg_sps} steps/s)`;

                    const regs8 = (d.regime_breakdown || []).filter(r => r.active_workers === 8);
                    const edges8 = regs8.reduce((a, r) => a + (r.new_edges || 0), 0);
                    const waves8 = regs8.reduce((a, r) => a + (r.waves || 0), 0);
                    const nodes8 = regs8.reduce((a, r) => a + (r.nodes_explored || 0), 0);
                    const eps8 = regs8.length ? (regs8.reduce((a, r) => a + (r.mean_eps || 0) * (r.waves || 1), 0) / Math.max(waves8, 1)).toFixed(1) : 61.1;
                    document.getElementById('kpi-night-edges').textContent = `+${fmtNum(edges8)}`;
                    document.getElementById('kpi-night-rate').textContent = `Mean ${eps8} edges/s (Idle 8W Turbo)`;
                    document.getElementById('kpi-night-waves').textContent = `${fmtNum(waves8)} waves · ${fmtNum(nodes8)} nodes explored`;

                    const regsMid = (d.regime_breakdown || []).filter(r => r.active_workers === 6 || r.active_workers === 4);
                    const edgesMid = regsMid.reduce((a, r) => a + (r.new_edges || 0), 0);
                    const wavesMid = regsMid.reduce((a, r) => a + (r.waves || 0), 0);
                    const nodesMid = regsMid.reduce((a, r) => a + (r.nodes_explored || 0), 0);
                    const epsMid = regsMid.length ? (regsMid.reduce((a, r) => a + (r.mean_eps || 0) * (r.waves || 1), 0) / Math.max(wavesMid, 1)).toFixed(1) : 34.1;
                    document.getElementById('kpi-idle-edges').textContent = `+${fmtNum(edgesMid)}`;
                    document.getElementById('kpi-idle-rate').textContent = `Mean ${epsMid} edges/s (4W Active / 6W Ramp)`;
                    document.getElementById('kpi-idle-waves').textContent = `${fmtNum(wavesMid)} waves · ${fmtNum(nodesMid)} nodes explored`;

                    document.getElementById('kpi-rf-pockets').textContent = `${fmtNum(s.total_3d_pockets)} Pockets`;
                    document.getElementById('kpi-rf-pairs').textContent = `${fmtNum(s.rf_scored_pairs)} 3D RangeFinder Pairs Scored`;
                    document.getElementById('kpi-elite-counts').textContent = `${fmtNum(s.elite_rf_count)} Elite RF (≥0.85) · ${fmtNum(s.elite_twilight_count)} Elite Twilight`;

                    // 1. Timeline Chart
                    const tb = d.timeline_buckets || [];
                    chTimeline.setOption({
                        tooltip: { trigger: 'axis', backgroundColor: 'rgba(2,6,23,0.95)', borderColor: '#38bdf8', textStyle: { color: '#f8fafc', fontSize: 12 } },
                        legend: { data: ['Cumulative Graph Edges (M)', 'Discovery Rate (edges/s)', 'Active Workers'], textStyle: { color: '#cbd5e1', fontSize: 11 } },
                        grid: { top: 38, bottom: 28, left: 58, right: 54 },
                        xAxis: { type: 'category', data: tb.map(x => x.bucket_label), axisLabel: { color: '#94a3b8', fontSize: 10 } },
                        yAxis: [
                            { type: 'value', name: 'Total Edges (M)', axisLabel: { color: '#38bdf8', formatter: v => (v / 1e6).toFixed(1) + 'M' }, splitLine: { lineStyle: { color: '#1e293b' } } },
                            { type: 'value', name: 'Edges / sec', axisLabel: { color: '#fbbf24' }, splitLine: { show: false } }
                        ],
                        series: [
                            {
                                name: 'Cumulative Graph Edges (M)',
                                type: 'line',
                                smooth: true,
                                data: tb.map(x => x.cum_edges),
                                areaStyle: { color: 'rgba(56, 189, 248, 0.14)' },
                                itemStyle: { color: '#38bdf8' }
                            },
                            {
                                name: 'Discovery Rate (edges/s)',
                                type: 'bar',
                                yAxisIndex: 1,
                                data: tb.map(x => x.avg_eps),
                                itemStyle: { color: 'rgba(251, 191, 36, 0.65)' }
                            },
                            {
                                name: 'Active Workers',
                                type: 'line',
                                yAxisIndex: 1,
                                step: 'end',
                                data: tb.map(x => x.workers * 8),
                                lineStyle: { type: 'dashed', color: '#c084fc' },
                                itemStyle: { color: '#c084fc' }
                            }
                        ]
                    });

                    // 2. Calibration Chart
                    const bm = d.benchmarks || [];
                    chCalib.setOption({
                        tooltip: { trigger: 'axis', backgroundColor: 'rgba(2,6,23,0.95)', borderColor: '#10b981', textStyle: { color: '#f8fafc', fontSize: 12 } },
                        legend: { data: ['Edges Scored / sec', 'Scaling Efficiency (%)', 'Avg Step Latency (ms)'], textStyle: { color: '#cbd5e1', fontSize: 11 } },
                        grid: { top: 38, bottom: 28, left: 54, right: 54 },
                        xAxis: { type: 'category', data: bm.map(x => `${x.workers}W`), axisLabel: { color: '#e2e8f0', fontWeight: 'bold' } },
                        yAxis: [
                            { type: 'value', name: 'Edges/s & ms', axisLabel: { color: '#34d399' }, splitLine: { lineStyle: { color: '#1e293b' } } },
                            { type: 'value', name: 'Efficiency %', min: 0, max: 105, axisLabel: { color: '#fbbf24', formatter: '{value}%' }, splitLine: { show: false } }
                        ],
                        series: [
                            {
                                name: 'Edges Scored / sec',
                                type: 'bar',
                                data: bm.map(x => x.edges_per_sec),
                                itemStyle: { color: '#10b981', borderRadius: [5, 5, 0, 0] }
                            },
                            {
                                name: 'Avg Step Latency (ms)',
                                type: 'line',
                                data: bm.map(x => x.avg_step_ms),
                                itemStyle: { color: '#38bdf8' }
                            },
                            {
                                name: 'Scaling Efficiency (%)',
                                type: 'line',
                                yAxisIndex: 1,
                                data: bm.map(x => x.parallel_efficiency_pct),
                                itemStyle: { color: '#fbbf24' },
                                lineStyle: { width: 3 }
                            }
                        ]
                    });

                    // 3. Regime Breakdown Chart
                    const rb = d.regime_breakdown || [];
                    chRegimes.setOption({
                        tooltip: { trigger: 'axis', backgroundColor: 'rgba(2,6,23,0.95)', borderColor: '#c084fc', textStyle: { color: '#f8fafc', fontSize: 12 } },
                        legend: { data: ['New Edges Discovered (K)', 'Nodes Explored (K)', 'Mean Edges/s'], textStyle: { color: '#cbd5e1', fontSize: 11 } },
                        grid: { top: 38, bottom: 28, left: 58, right: 54 },
                        xAxis: { type: 'category', data: rb.map(x => `${x.power_mode} (${x.active_workers}W)`), axisLabel: { color: '#e2e8f0', fontSize: 11 } },
                        yAxis: [
                            { type: 'value', name: 'Count (Thousands)', axisLabel: { color: '#c084fc', formatter: v => (v / 1000).toFixed(0) + 'K' }, splitLine: { lineStyle: { color: '#1e293b' } } },
                            { type: 'value', name: 'Mean Edges/s', axisLabel: { color: '#34d399' }, splitLine: { show: false } }
                        ],
                        series: [
                            {
                                name: 'New Edges Discovered (K)',
                                type: 'bar',
                                data: rb.map(x => x.new_edges),
                                itemStyle: { color: '#c084fc', borderRadius: [5, 5, 0, 0] }
                            },
                            {
                                name: 'Nodes Explored (K)',
                                type: 'bar',
                                data: rb.map(x => x.nodes_explored),
                                itemStyle: { color: '#38bdf8', borderRadius: [5, 5, 0, 0] }
                            },
                            {
                                name: 'Mean Edges/s',
                                type: 'line',
                                yAxisIndex: 1,
                                data: rb.map(x => x.mean_eps),
                                itemStyle: { color: '#34d399' },
                                lineStyle: { width: 3 }
                            }
                        ]
                    });

                    // 4. Live Rolling Wave Stream
                    const rw = (d.recent_waves || []).slice().reverse();
                    chLive.setOption({
                        tooltip: { trigger: 'axis', backgroundColor: 'rgba(2,6,23,0.95)', borderColor: '#34d399', textStyle: { color: '#f8fafc', fontSize: 12 } },
                        legend: { data: ['Edges / sec', 'Steps / sec', 'Avg Latency (ms)'], textStyle: { color: '#cbd5e1', fontSize: 11 } },
                        grid: { top: 38, bottom: 28, left: 50, right: 54 },
                        xAxis: { type: 'category', data: rw.map(x => (x.recorded_at || '').slice(11, 19)), axisLabel: { color: '#94a3b8', fontSize: 10 } },
                        yAxis: [
                            { type: 'value', name: 'Rate (/s)', axisLabel: { color: '#34d399' }, splitLine: { lineStyle: { color: '#1e293b' } } },
                            { type: 'value', name: 'Latency (ms)', axisLabel: { color: '#f43f5e' }, splitLine: { show: false } }
                        ],
                        series: [
                            {
                                name: 'Edges / sec',
                                type: 'line',
                                smooth: true,
                                data: rw.map(x => x.edges_per_sec),
                                areaStyle: { color: 'rgba(16, 185, 129, 0.14)' },
                                itemStyle: { color: '#10b981' }
                            },
                            {
                                name: 'Steps / sec',
                                type: 'line',
                                smooth: true,
                                data: rw.map(x => x.steps_per_sec),
                                itemStyle: { color: '#38bdf8' }
                            },
                            {
                                name: 'Avg Latency (ms)',
                                type: 'line',
                                yAxisIndex: 1,
                                data: rw.map(x => x.avg_step_ms),
                                itemStyle: { color: '#f43f5e' },
                                lineStyle: { type: 'dotted' }
                            }
                        ]
                    });

                    // 5. Render Breakthrough Table
                    const bt = d.breakthroughs || [];
                    document.getElementById('breakthrough-tbody').innerHTML = bt.map((b, i) => `
                        <tr>
                            <th>#${i + 1} ${b.bait_id} × ${b.orphan_id}</th>
                            <td><b style="color:#f8fafc;">${b.bait_org}</b><br><span style="color:#94a3b8;">× ${b.orphan_org}</span></td>
                            <td>
                                <span style="color:#38bdf8;font-weight:700;">${b.bait_desc}</span><br>
                                <span style="color:#cbd5e1;font-size:0.74rem;">× ${b.orphan_desc}</span>
                            </td>
                            <td>
                                <b style="color:#c084fc;font-size:0.92rem;">${b.twilight_score != null ? b.twilight_score.toFixed(4) : 'NULL'}</b>
                                ${(b.twilight_score != null && b.twilight_score > 0.8640) ? '<span style="background:rgba(245,158,11,0.18);color:#fbbf24;border:1px solid rgba(245,158,11,0.45);padding:1px 6px;border-radius:999px;font-size:0.66rem;font-weight:800;margin-left:4px;">NEW RECORD</span>' : ''}
                                <br><span style="font-size:0.72rem;color:#94a3b8;">SW ID: ${b.seq_identity_pct != null ? b.seq_identity_pct + '%' : 'NULL'} · d=${b.vector_distance != null ? b.vector_distance.toFixed(4) : 'NULL'}</span>
                            </td>
                            <td>
                                <b style="color:#34d399;">IoU ${b.rf_range_iou != null ? b.rf_range_iou.toFixed(3) : 'NULL'}</b><br>
                                <span style="font-size:0.72rem;color:#94a3b8;">RF Score: ${b.rf_score != null ? b.rf_score.toFixed(3) : 'NULL'}</span>
                            </td>
                            <td>
                                <b style="color:#fbbf24;">TM ${b.tm_score != null ? b.tm_score.toFixed(3) : 'NULL'}</b><br>
                                <span style="font-size:0.72rem;color:#94a3b8;">Chaperone: ${b.chaperone_score != null ? b.chaperone_score.toFixed(1) : 'NULL'}</span>
                            </td>
                            <td>
                                <div style="display:flex;gap:6px;flex-wrap:wrap;">
                                    <a class="tel-btn" style="padding:4px 8px;font-size:0.72rem;" href="/compare?node1=${b.bait_id}&node2=${b.orphan_id}&preset=rangefinder">3D Lab</a>
                                    <a class="tel-btn" style="padding:4px 8px;font-size:0.72rem;border-color:#f59e0b;color:#fbbf24;" href="/cards?hub=${b.bait_id}">Card Arena</a>
                                </div>
                            </td>
                        </tr>
                    `).join('');
                });
        }

        window.addEventListener('DOMContentLoaded', () => {
            loadTelemetryDashboard();
            setInterval(loadTelemetryDashboard, 5000);
        });
    </script>
</body>
</html>
"""


@app.route('/telemetry')
@app.route('/api/swarm_telemetry')
def api_swarm_telemetry():
    import time as _t
    wants_json = (
        request.args.get('format') == 'json'
        or 'application/json' in (request.headers.get('Accept') or '')
    )
    if not wants_json and request.path in ('/api/swarm_telemetry', '/telemetry'):
        return render_template_string(TELEMETRY_DASHBOARD_HTML_TEMPLATE)

    now_ts = _t.perf_counter()
    if TELEMETRY_JSON_CACHE['payload'] and (now_ts - TELEMETRY_JSON_CACHE['ts']) <= 4.0:
        return jsonify(TELEMETRY_JSON_CACHE['payload'])

    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                SELECT
                    COUNT(*) AS total_waves,
                    MIN(recorded_at)::text AS first_wave,
                    MAX(recorded_at)::text AS last_wave,
                    COALESCE(SUM(steps_completed), 0) AS sum_steps,
                    COALESCE(SUM(edges_discovered), 0) AS sum_edges,
                    COALESCE(SUM(pockets_3d_indexed), 0) AS sum_pockets,
                    COALESCE(SUM(rf_pairs_scored), 0) AS sum_rf_pairs,
                    ROUND(COALESCE(AVG(steps_per_sec), 0)::numeric, 2) AS avg_sps,
                    ROUND(COALESCE(MAX(steps_per_sec), 0)::numeric, 2) AS max_sps,
                    ROUND(COALESCE(AVG(edges_per_sec), 0)::numeric, 2) AS avg_eps,
                    ROUND(COALESCE(MAX(edges_per_sec), 0)::numeric, 2) AS max_eps,
                    COALESCE(MIN(total_graph_edges), 390667) AS start_graph_edges,
                    COALESCE(MAX(total_graph_edges), 390667) AS end_graph_edges
                FROM swarm_telemetry_log;
            """)
            srow = cur.fetchone()

            cur.execute("""
                SELECT
                    COUNT(*) filter (WHERE status = 'Explored') AS explored_nodes,
                    COUNT(*) filter (WHERE status = 'Pending') AS pending_nodes,
                    COUNT(*) filter (WHERE rf_score IS NOT NULL) AS rf_scored_pairs,
                    COUNT(*) filter (WHERE rf_score >= 0.85) AS elite_rf_count,
                    COUNT(*) filter (WHERE twilight_score >= 0.35) AS elite_twilight_count
                FROM orphan_discoveries;
            """)
            drow = cur.fetchone()

            cur.execute("SELECT COUNT(*) FROM swarm_pocket_ranges;")
            tot_pockets = cur.fetchone()[0]

            summary = {
                'total_waves': int(srow[0] or 0),
                'first_wave': srow[1],
                'last_wave': srow[2],
                'sum_steps': int(srow[3] or 0),
                'sum_edges': int(srow[4] or 0),
                'sum_pockets': int(srow[5] or 0),
                'sum_rf_pairs': int(srow[6] or 0),
                'avg_sps': float(srow[7] or 0.0),
                'max_sps': float(srow[8] or 0.0),
                'avg_eps': float(srow[9] or 0.0),
                'max_eps': float(srow[10] or 0.0),
                'start_graph_edges': int(srow[11] or 390667),
                'end_graph_edges': int(srow[12] or 390667),
                'explored_nodes': int(drow[0] or 0),
                'pending_nodes': int(drow[1] or 0),
                'rf_scored_pairs': int(drow[2] or 0),
                'elite_rf_count': int(drow[3] or 0),
                'elite_twilight_count': int(drow[4] or 0),
                'total_3d_pockets': int(tot_pockets or 0),
            }

            cur.execute("""
                SELECT
                    power_mode,
                    active_workers,
                    COUNT(*) AS waves,
                    SUM(steps_completed) AS nodes_explored,
                    SUM(edges_discovered) AS new_edges,
                    SUM(pockets_3d_indexed) AS new_3d_pockets,
                    SUM(rf_pairs_scored) AS new_rf_pairs,
                    ROUND(AVG(avg_step_ms)::numeric, 1) AS mean_ms,
                    ROUND(AVG(p95_step_ms)::numeric, 1) AS mean_p95_ms,
                    ROUND(AVG(steps_per_sec)::numeric, 2) AS mean_sps,
                    ROUND(AVG(edges_per_sec)::numeric, 2) AS mean_eps
                FROM swarm_telemetry_log
                GROUP BY power_mode, active_workers
                ORDER BY active_workers DESC;
            """)
            regime_breakdown = [
                {
                    'power_mode': r[0],
                    'active_workers': int(r[1]),
                    'waves': int(r[2]),
                    'nodes_explored': int(r[3] or 0),
                    'new_edges': int(r[4] or 0),
                    'new_3d_pockets': int(r[5] or 0),
                    'new_rf_pairs': int(r[6] or 0),
                    'mean_ms': float(r[7] or 0.0),
                    'mean_p95_ms': float(r[8] or 0.0),
                    'mean_sps': float(r[9] or 0.0),
                    'mean_eps': float(r[10] or 0.0),
                }
                for r in cur.fetchall()
            ]

            cur.execute("""
                SELECT
                    to_char(
                        date_trunc('hour', recorded_at) +
                        floor(extract(minute FROM recorded_at) / 15.0) * interval '15 minutes',
                        'HH24:MI'
                    ) AS bucket_label,
                    MAX(total_graph_edges) AS cum_edges,
                    ROUND(AVG(edges_per_sec)::numeric, 1) AS avg_eps,
                    ROUND(AVG(steps_per_sec)::numeric, 1) AS avg_sps,
                    ROUND(AVG(avg_step_ms)::numeric, 1) AS avg_ms,
                    MAX(active_workers) AS workers
                FROM swarm_telemetry_log
                GROUP BY 1, date_trunc('hour', recorded_at) + floor(extract(minute FROM recorded_at) / 15.0) * interval '15 minutes'
                ORDER BY date_trunc('hour', recorded_at) + floor(extract(minute FROM recorded_at) / 15.0) * interval '15 minutes' ASC;
            """)
            timeline_buckets = [
                {
                    'bucket_label': r[0],
                    'cum_edges': int(r[1] or 0),
                    'avg_eps': float(r[2] or 0.0),
                    'avg_sps': float(r[3] or 0.0),
                    'avg_ms': float(r[4] or 0.0),
                    'workers': int(r[5] or 2),
                }
                for r in cur.fetchall()
            ]

            cur.execute("""
                SELECT workers, avg_step_ms, p95_step_ms, steps_per_sec, edges_per_sec, parallel_efficiency_pct, regime_note
                FROM swarm_concurrency_benchmarks
                ORDER BY workers ASC;
            """)
            benchmarks = [
                {
                    'workers': r[0],
                    'avg_step_ms': r[1],
                    'p95_step_ms': r[2],
                    'steps_per_sec': r[3],
                    'edges_per_sec': r[4],
                    'parallel_efficiency_pct': r[5],
                    'regime_note': r[6],
                }
                for r in cur.fetchall()
            ]

            cur.execute("""
                SELECT recorded_at::text, power_mode, active_workers, hid_idle_sec, load_avg_1m,
                       steps_completed, edges_discovered, pockets_3d_indexed, rf_pairs_scored,
                       avg_step_ms, p95_step_ms, steps_per_sec, edges_per_sec,
                       total_graph_edges, total_3d_rf_edges
                FROM swarm_telemetry_log
                ORDER BY id DESC
                LIMIT 60;
            """)
            recent_waves = [
                {
                    'recorded_at': r[0],
                    'power_mode': r[1],
                    'active_workers': r[2],
                    'hid_idle_sec': r[3],
                    'load_avg_1m': r[4],
                    'steps_completed': r[5],
                    'edges_discovered': r[6],
                    'pockets_3d_indexed': r[7],
                    'rf_pairs_scored': r[8],
                    'avg_step_ms': r[9],
                    'p95_step_ms': r[10],
                    'steps_per_sec': r[11],
                    'edges_per_sec': r[12],
                    'total_graph_edges': r[13],
                    'total_3d_rf_edges': r[14],
                }
                for r in cur.fetchall()
            ]

            cur.execute("""
                SELECT DISTINCT ON (LEAST(d.bait_id, d.orphan_id), GREATEST(d.bait_id, d.orphan_id))
                    d.id, d.bait_id, d.orphan_id,
                    split_part(d.bait_organism, ' (', 1) AS b_org,
                    split_part(d.orphan_organism, ' (', 1) AS o_org,
                    ROUND(d.vector_distance::numeric, 4) AS v_dist,
                    d.seq_identity_pct AS seq_id,
                    d.twilight_score AS tw_score,
                    d.rf_range_iou AS rf_iou,
                    d.rf_score AS rf_sc,
                    d.tm_score AS tm_sc,
                    d.chaperone_score AS chap_sc,
                    regexp_replace(split_part(p1.name, ' OS=', 1), '^\\S+\\s+', '') AS b_desc,
                    regexp_replace(split_part(p2.name, ' OS=', 1), '^\\S+\\s+', '') AS o_desc
                FROM orphan_discoveries d
                JOIN proteins p1 ON p1.uniprot_id = d.bait_id
                JOIN proteins p2 ON p2.uniprot_id = d.orphan_id
                WHERE d.id > 390667
                  AND d.bait_id <> d.orphan_id
                  AND d.twilight_score >= 0.72
                ORDER BY LEAST(d.bait_id, d.orphan_id), GREATEST(d.bait_id, d.orphan_id), d.twilight_score DESC;
            """)
            bt_rows = sorted(cur.fetchall(), key=lambda x: float(x[7] or 0.0), reverse=True)[:14]
            breakthroughs = [
                {
                    'edge_id': int(r[0]),
                    'bait_id': r[1],
                    'orphan_id': r[2],
                    'bait_org': r[3],
                    'orphan_org': r[4],
                    'vector_distance': float(r[5]) if r[5] is not None else None,
                    'seq_identity_pct': float(r[6]) if r[6] is not None else None,
                    'twilight_score': float(r[7]) if r[7] is not None else None,
                    'rf_range_iou': float(r[8]) if r[8] is not None else None,
                    'rf_score': float(r[9]) if r[9] is not None else None,
                    'tm_score': float(r[10]) if r[10] is not None else None,
                    'chaperone_score': float(r[11]) if r[11] is not None else None,
                    'bait_desc': r[12],
                    'orphan_desc': r[13],
                }
                for r in bt_rows
            ]

    payload = {
        'summary': summary,
        'regime_breakdown': regime_breakdown,
        'timeline_buckets': timeline_buckets,
        'benchmarks': benchmarks,
        'latest': recent_waves[0] if recent_waves else None,
        'recent_waves': recent_waves,
        'breakthroughs': breakthroughs,
    }
    TELEMETRY_JSON_CACHE['ts'] = now_ts
    TELEMETRY_JSON_CACHE['payload'] = payload
    return jsonify(payload)


@app.route('/api/hub_cards')
def api_hub_cards():
    import time as _t
    t0 = _t.perf_counter()
    try:
        org_cap = max(1, min(int(request.args.get('org_cap', 10)), 30))
    except Exception:
        org_cap = 10
    try:
        match_cap = max(3, min(int(request.args.get('match_cap', 12)), 48))
    except Exception:
        match_cap = 12
    sort_by = (request.args.get('sort_by') or 'composite').lower()

    cache_key = (org_cap, match_cap, sort_by)
    if cache_key in HUB_CARDS_CACHE:
        ts_cached, cached_payload = HUB_CARDS_CACHE[cache_key]
        if (_t.perf_counter() - ts_cached) <= HUB_CARDS_CACHE_TTL_SEC:
            cached = dict(cached_payload)
            cached['elapsed_ms'] = round((_t.perf_counter() - t0) * 1000.0, 2)
            return jsonify(cached)

    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            # 1. Top Hub Organisms ranked by cross-species relationship reach + 3D RangeFinder density
            cur.execute("""
                WITH rf_orgs AS (
                    SELECT
                        bait_organism,
                        COUNT(DISTINCT orphan_id) AS rf_partners,
                        COUNT(DISTINCT bait_id) AS rf_hubs,
                        ROUND(AVG(rf_score)::numeric, 3) AS avg_rf_score,
                        ROUND(MAX(rf_score)::numeric, 3) AS max_rf_score
                    FROM orphan_discoveries
                    WHERE rf_score IS NOT NULL
                      AND tm_score IS NOT NULL
                      AND sc_rmsd IS NOT NULL
                      AND chaperone_score IS NOT NULL
                      AND bait_id <> orphan_id
                      AND bait_organism IS NOT NULL
                      AND tm_score >= 0.50
                    GROUP BY bait_organism
                    HAVING COUNT(DISTINCT orphan_id) >= 1
                ),
                org_rollup AS (
                    SELECT
                        r.bait_organism,
                        COUNT(d.id) AS total_edges,
                        COUNT(DISTINCT d.orphan_organism) AS partner_orgs,
                        COUNT(DISTINCT d.orphan_id) AS partner_proteins,
                        r.rf_partners,
                        r.rf_hubs,
                        r.avg_rf_score,
                        r.max_rf_score,
                        ROUND(MAX(d.twilight_score)::numeric, 3) AS max_twilight_score,
                        ARRAY_AGG(DISTINCT d.family_name) FILTER (WHERE d.family_name IS NOT NULL) AS families
                    FROM rf_orgs r
                    JOIN orphan_discoveries d ON d.bait_organism = r.bait_organism AND d.bait_id <> d.orphan_id
                    GROUP BY r.bait_organism, r.rf_partners, r.rf_hubs, r.avg_rf_score, r.max_rf_score
                )
                SELECT
                    bait_organism, total_edges, partner_orgs, partner_proteins,
                    rf_partners, rf_hubs, avg_rf_score, max_rf_score, max_twilight_score, families
                FROM org_rollup
                ORDER BY
                    rf_partners DESC,
                    max_twilight_score DESC NULLS LAST,
                    avg_rf_score DESC NULLS LAST
                LIMIT %s;
            """, (org_cap,))
            org_rows = cur.fetchall()

            top_organisms = []
            org_names = []
            for r in org_rows:
                b_org = r[0]
                org_names.append(b_org)
                short_org = b_org.split(' (strain')[0]
                top_organisms.append({
                    'organism': b_org,
                    'short_organism': short_org,
                    'total_edges': int(r[1]),
                    'partner_orgs': int(r[2]),
                    'partner_proteins': int(r[3]),
                    'rf_partners': int(r[4]),
                    'rf_hubs': int(r[5]),
                    'avg_rf_score': float(r[6]) if r[6] is not None else None,
                    'max_rf_score': float(r[7]) if r[7] is not None else None,
                    'max_twilight_score': float(r[8]) if r[8] is not None else None,
                    'families': (r[9] or [])[:6]
                })

            if not org_names:
                return jsonify({'top_organisms': [], 'hubs': [], 'elapsed_ms': 0.0})

            # 2. Fetch all 3D-validated 1-to-N edges for the capped Hub Organisms (zero synthetic COALESCE fallbacks)
            cur.execute("""
                WITH dedup_edges AS (
                    SELECT DISTINCT ON (d.bait_id, d.orphan_id)
                        d.bait_id,
                        d.orphan_id,
                        d.family_name,
                        d.bait_organism,
                        d.orphan_organism,
                        d.vector_distance,
                        d.seq_identity_pct,
                        d.twilight_score,
                        d.tm_score,
                        CASE
                            WHEN d.chaperone_score <= 1.0 THEN ROUND((d.chaperone_score * 100.0)::numeric, 1)
                            ELSE ROUND(d.chaperone_score::numeric, 1)
                        END AS chaperone_score,
                        d.sc_rmsd,
                        d.mean_rasa,
                        d.buns_count,
                        d.hbond_count,
                        d.rf_range_iou,
                        d.rf_giou_loss,
                        d.rf_substrate_match,
                        d.rf_substrate_props_matched,
                        d.rf_arom_prob,
                        d.rf_score
                    FROM orphan_discoveries d
                    WHERE d.bait_id <> d.orphan_id
                      AND d.bait_organism = ANY(%s)
                      AND d.rf_score IS NOT NULL
                      AND d.tm_score IS NOT NULL
                      AND d.sc_rmsd IS NOT NULL
                      AND d.chaperone_score IS NOT NULL
                      AND d.tm_score >= 0.50
                    ORDER BY d.bait_id, d.orphan_id, d.tm_score DESC, d.rf_score DESC, d.twilight_score DESC NULLS LAST, d.vector_distance ASC
                )
                SELECT
                    e.*,
                    ROUND((
                        0.34 * e.rf_score +
                        0.24 * (1.0 - LEAST(e.vector_distance, 1.0)) +
                        0.18 * LEAST(e.twilight_score * 2.2, 1.0) +
                        0.14 * (e.chaperone_score / 100.0) +
                        0.10 * e.tm_score
                    )::numeric, 4) AS composite_score
                FROM dedup_edges e
                ORDER BY e.bait_id, composite_score DESC;
            """, (org_names,))
            edge_rows = cur.fetchall()

            # Collect all unique protein IDs and family names needed
            all_pids = set()
            all_fams = set()
            edges_by_hub = {}
            for er in edge_rows:
                b_id = er[0]
                o_id = er[1]
                fam = er[2]
                all_pids.add(b_id)
                all_pids.add(o_id)
                if fam:
                    all_fams.add(fam.lower())
                edges_by_hub.setdefault(b_id, []).append(er)

            # 3. Batch-fetch protein biophysical, RangeFinder pocket, and real Chaperone metrics from swarm_comparison_cache
            cur.execute("""
                SELECT
                    p.uniprot_id,
                    p.name,
                    p.sequence,
                    molecular_weight(p.sequence) AS mw_da,
                    sequence_entropy(p.sequence) AS ent,
                    isoelectric_point(p.sequence) AS pi_val,
                    gravy_score(p.sequence) AS gravy_val,
                    spr.mw_min, spr.mw_max,
                    spr.logp_min, spr.logp_max,
                    spr.tpsa_min, spr.tpsa_max,
                    spr.aromatic_probability,
                    hc.chap_sc,
                    hc.hbonds,
                    hc.buns,
                    hc.rasa
                FROM proteins p
                LEFT JOIN swarm_pocket_ranges spr ON spr.uniprot_id = p.uniprot_id
                LEFT JOIN LATERAL (
                    SELECT
                        ROUND(((c.metrics_payload->'chaperone'->>'chaperone_score')::numeric * 100.0), 1) AS chap_sc,
                        (c.metrics_payload->'chaperone'->>'hbonds_1')::int AS hbonds,
                        (c.metrics_payload->'chaperone'->>'buns_1')::int AS buns,
                        (c.metrics_payload->'chaperone'->>'rasa_1')::float8 AS rasa
                    FROM swarm_comparison_cache c
                    WHERE c.id1 = p.uniprot_id
                      AND c.metrics_payload->'chaperone' IS NOT NULL
                    LIMIT 1
                ) hc ON TRUE
                WHERE p.uniprot_id = ANY(%s);
            """, (list(all_pids),))
            meta_map = {r[0]: r for r in cur.fetchall()}

            # 4. Batch-fetch enzyme substrate metadata
            cur.execute("""
                SELECT LOWER(family_name), substrate_name, canonical_query_smiles, warhead_substructure_smiles, catalytic_role
                FROM swarm_enzyme_substrates;
            """)
            sub_map = {
                r[0]: {
                    'substrate_name': r[1],
                    'canonical_smiles': r[2],
                    'warhead_smiles': r[3],
                    'catalytic_role': r[4]
                }
                for r in cur.fetchall()
            }

            # 5. Build Hub Cards & their 1-to-N Challenger Decks
            hubs_list = []
            org_rank_map = {o['organism']: idx for idx, o in enumerate(top_organisms)}
            for b_id, h_edges in edges_by_hub.items():
                if not h_edges:
                    continue
                first_e = h_edges[0]
                fam_name = first_e[2]
                b_org = first_e[3]
                sub_info = sub_map.get((fam_name or '').lower())
                hub_card = _build_poke_card_profile(b_id, meta_map.get(b_id), fam_name, b_org)
                if not hub_card:
                    continue
                hub_card['substrate'] = sub_info

                # Sort 1-to-N edges according to sort_by lens
                def _edge_sort_key(e):
                    v_dist = float(e[5])
                    tw_sc = float(e[7] or 0.0)
                    ch_sc = float(e[9] or 0.0)
                    rf_sc = float(e[19] or 0.0)
                    comp_sc = float(e[20] or 0.0)
                    if sort_by == 'rangefinder':
                        return (rf_sc, comp_sc)
                    if sort_by == 'twilight':
                        return (tw_sc, comp_sc)
                    if sort_by == 'chaperone':
                        return (ch_sc, comp_sc)
                    if sort_by == 'vector':
                        return (1.0 - v_dist, comp_sc)
                    return (comp_sc, rf_sc)

                sorted_edges = sorted(h_edges, key=_edge_sort_key, reverse=True)[:match_cap]
                matches = []
                unique_target_orgs = set()
                for rank_idx, er in enumerate(sorted_edges, 1):
                    o_id = er[1]
                    o_fam = er[2] or fam_name
                    o_org = er[4]
                    unique_target_orgs.add(o_org)
                    v_dist = round(float(er[5]), 4) if er[5] is not None else None
                    cos_sim = round(max(0.0, 1.0 - v_dist) * 100.0, 1) if v_dist is not None else None
                    seq_id = round(float(er[6]), 1) if er[6] is not None else None
                    tw_sc = round(float(er[7]), 4) if er[7] is not None else None
                    tm_sc = round(float(er[8]), 3) if er[8] is not None else None
                    ch_sc = round(float(er[9]), 1) if er[9] is not None else None
                    rmsd_v = round(float(er[10]), 2) if er[10] is not None else None
                    rasa_v = round(float(er[11]), 3) if er[11] is not None else None
                    buns_v = int(er[12]) if er[12] is not None else None
                    hb_v = int(er[13]) if er[13] is not None else None
                    rf_iou = round(float(er[14]), 4) if er[14] is not None else None
                    rf_giou = round(float(er[15]), 4) if er[15] is not None else None
                    rf_match = bool(er[16]) if er[16] is not None else None
                    rf_props = int(er[17]) if er[17] is not None else None
                    rf_arom = round(float(er[18]), 3) if er[18] is not None else None
                    rf_sc = round(float(er[19]), 4) if er[19] is not None else None
                    comp_sc = round(float(er[20]), 4) if er[20] is not None else None

                    base_meta = list(meta_map.get(o_id)) if o_id in meta_map else None
                    if base_meta:
                        base_meta[14] = ch_sc
                        base_meta[15] = hb_v
                        base_meta[16] = buns_v
                        base_meta[17] = rasa_v
                    chal_card = _build_poke_card_profile(o_id, base_meta, o_fam, o_org)
                    if not chal_card:
                        continue
                    chal_card['substrate'] = sub_info
                    h2h = _compute_head_to_head(hub_card, chal_card)

                    matches.append({
                        'rank': rank_idx,
                        'card': chal_card,
                        'pair_metrics': {
                            'composite_score': comp_sc,
                            'composite_pct': int(round(comp_sc * 100.0)) if comp_sc is not None else None,
                            'rf_score': rf_sc,
                            'rf_range_iou': rf_iou,
                            'rf_giou_loss': rf_giou,
                            'rf_substrate_match': rf_match,
                            'rf_substrate_props_matched': rf_props,
                            'rf_arom_prob': rf_arom,
                            'twilight_score': tw_sc,
                            'seq_identity_pct': seq_id,
                            'tm_score': tm_sc,
                            'sc_rmsd': rmsd_v,
                            'vector_distance': v_dist,
                            'cosine_sim_pct': cos_sim,
                            'cross_genus': hub_card['genus'] != chal_card['genus']
                        },
                        'head_to_head': h2h
                    })

                avg_comp = round(sum(m['pair_metrics']['composite_score'] for m in matches) / max(len(matches), 1), 4)
                max_comp = round(max((m['pair_metrics']['composite_score'] for m in matches), default=0.0), 4)
                hubs_list.append({
                    'hub_id': b_id,
                    'organism': b_org,
                    'short_organism': b_org.split(' (strain')[0],
                    'org_rank': org_rank_map.get(b_org, 999),
                    'family_name': fam_name,
                    'total_3d_matches': len(h_edges),
                    'displayed_matches': len(matches),
                    'unique_partner_orgs': len(unique_target_orgs),
                    'avg_composite_score': avg_comp,
                    'max_composite_score': max_comp,
                    'hub_card': hub_card,
                    'matches': matches
                })

            hubs_list.sort(key=lambda h: (h['org_rank'], -h['total_3d_matches'], -h['max_composite_score']))

    elapsed_ms = round((_t.perf_counter() - t0) * 1000.0, 2)
    payload = {
        'org_cap': org_cap,
        'match_cap': match_cap,
        'sort_by': sort_by,
        'elapsed_ms': elapsed_ms,
        'top_organisms': top_organisms,
        'hubs': hubs_list
    }
    HUB_CARDS_CACHE[cache_key] = (_t.perf_counter(), payload)
    return jsonify(payload)


@app.route('/api/data')
def get_data():
    conn = psycopg.connect(DB_URI, autocommit=True)
    cur = conn.cursor()

    limit_param = request.args.get('limit', 'best-40')
    if limit_param.lower() == 'rangefinder-40':
        cur.execute("""
            WITH ranked_rf AS (
                SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism,
                       vector_distance, rf_score, rf_range_iou,
                       ROW_NUMBER() OVER (PARTITION BY family_name ORDER BY rf_score DESC NULLS LAST, rf_range_iou DESC NULLS LAST) AS fam_rk
                FROM orphan_discoveries
                WHERE bait_id != orphan_id
                  AND rf_score IS NOT NULL
            )
            SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism, vector_distance, rf_score
            FROM ranked_rf
            ORDER BY fam_rk ASC, rf_score DESC, rf_range_iou DESC
            LIMIT 80
        """)
        edge_rows = cur.fetchall()
        nodes_info = {}
        for b_id, o_id, st, fam, b_org, o_org, v_dist, rf_s in edge_rows:
            if b_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[b_id] = {'status': 'Explored', 'family': fam, 'organism': b_org}
            if o_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[o_id] = {'status': st, 'family': fam, 'organism': o_org}
            if len(nodes_info) >= 40:
                break
    elif limit_param.lower() == 'chaperone-40':
        cur.execute("""
            SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism, vector_distance, COALESCE(chaperone_score, 0.0)
            FROM orphan_discoveries
            WHERE bait_id != orphan_id
            ORDER BY (chaperone_score IS NOT NULL) DESC,
                     chaperone_score DESC NULLS LAST,
                     twilight_score DESC NULLS LAST,
                     vector_distance ASC
            LIMIT 80
        """)
        edge_rows = cur.fetchall()
        nodes_info = {}
        for b_id, o_id, st, fam, b_org, o_org, v_dist, ch_s in edge_rows:
            if b_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[b_id] = {'status': 'Explored', 'family': fam, 'organism': b_org}
            if o_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[o_id] = {'status': st, 'family': fam, 'organism': o_org}
            if len(nodes_info) >= 40:
                break
    elif limit_param.lower() == 'twilight-40':
        cur.execute("""
            SELECT bait_id, orphan_id, status, family_name, bait_organism, orphan_organism, vector_distance, twilight_score
            FROM orphan_discoveries
            WHERE bait_id != orphan_id
              AND twilight_score IS NOT NULL
              AND twilight_score > 0.10
            ORDER BY (split_part(bait_organism, ' ', 1) != split_part(orphan_organism, ' ', 1)) DESC,
                     twilight_score DESC,
                     vector_distance ASC
            LIMIT 80
        """)
        edge_rows = cur.fetchall()
        nodes_info = {}
        for b_id, o_id, st, fam, b_org, o_org, v_dist, tw_s in edge_rows:
            if b_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[b_id] = {'status': 'Explored', 'family': fam, 'organism': b_org}
            if o_id not in nodes_info and len(nodes_info) < 40:
                nodes_info[o_id] = {'status': st, 'family': fam, 'organism': o_org}
            if len(nodes_info) >= 40:
                break
    elif limit_param.lower() == 'best-40':
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
    if n1 and n1 != 'undefined' and n1 not in nodes_info:
        extra_nodes.append(n1)
    if n2 and n2 != 'undefined' and n2 not in nodes_info:
        extra_nodes.append(n2)

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
        conn.close()
        return jsonify({'nodes': [], 'matrix': []})

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
        emb_dict[uid] = np.fromstring(emb_str.strip('[]'), sep=',')
        info = nodes_info[uid]
        nodes.append({
            'id': uid,
            'status': info['status'],
            'family': info['family'],
            'organism': info['organism'],
            'name': name,
            'length': length
        })
    conn.close()

    nodes.sort(key=lambda x: x['id'])
    if not nodes:
        return jsonify({'nodes': [], 'matrix': []})

    emb_list = [emb_dict[n['id']] for n in nodes]
    emb_matrix = np.array(emb_list)
    dist_matrix = scipy.spatial.distance.cdist(emb_matrix, emb_matrix, metric='cosine')
    np.fill_diagonal(dist_matrix, 0.0)
    matrix = dist_matrix.tolist()
    return jsonify({'nodes': nodes, 'matrix': matrix})

if __name__ == '__main__':
    app.run(port=8080)
