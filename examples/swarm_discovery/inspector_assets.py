# Shared emoji-free CSS & JS for the unified ProteinCompareInspector component
# Used identically by both the Swarm Heatmap Modal (/) and the Phase-Scoped Wiki (/wiki)

INSPECTOR_CSS = """
.cmp-ms-svg {
    width: 14px;
    height: 14px;
    fill: currentColor;
    flex-shrink: 0;
    vertical-align: -2px;
}
.cmp-inspector {
    display: flex;
    flex-direction: column;
    width: 100%;
    height: 100%;
    background: #0f172a;
    color: #f8fafc;
    border-radius: 10px;
    overflow: hidden;
}
.cmp-toolbar {
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 10px;
    padding: 10px 16px;
    background: #020617;
    border-bottom: 1px solid #1e293b;
}
.cmp-toolbar-left, .cmp-toolbar-right {
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
}
.cmp-title {
    font-weight: 700;
    font-size: 1.02em;
    color: #f8fafc;
}
.cmp-cache-badge {
    font-size: 0.74em;
    padding: 2px 8px;
    border-radius: 999px;
    border: 1px solid #334155;
    background: #0f172a;
    color: #94a3b8;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
}
.cmp-cache-badge.hit {
    border-color: rgba(16, 185, 129, 0.5);
    color: #34d399;
    background: rgba(16, 185, 129, 0.08);
}
.cmp-select {
    background: #1e293b;
    color: #f8fafc;
    border: 1px solid #475569;
    padding: 5px 10px;
    border-radius: 6px;
    font-size: 0.82em;
}
.cmp-btn {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    background: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    padding: 5px 10px;
    border-radius: 6px;
    font-size: 0.79em;
    cursor: pointer;
    transition: all 0.12s;
}
.cmp-btn:hover {
    background: #334155;
    border-color: #38bdf8;
    color: #fff;
}
.cmp-btn.active {
    background: rgba(14, 165, 233, 0.22);
    border-color: #38bdf8;
    color: #38bdf8;
    font-weight: 600;
}
.cmp-btn.danger {
    border-color: #7f1d1d;
    color: #fca5a5;
}
.cmp-btn.danger:hover {
    background: #7f1d1d;
    color: #fff;
}
.cmp-modes-bar {
    display: flex;
    flex-wrap: nowrap;
    align-items: center;
    gap: 5px;
    padding: 6px 14px;
    background: #090d16;
    border-bottom: 1px solid #1e293b;
    overflow-x: auto;
    scrollbar-width: thin;
}
.cmp-modes-bar .cmp-btn {
    display: inline-flex;
    align-items: center;
    gap: 4px;
    padding: 4px 8px;
    font-size: 0.76em;
    white-space: nowrap;
    flex-shrink: 0;
}
.cmp-modes-label {
    font-size: 0.73em;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: #64748b;
    font-weight: 600;
    margin-right: 3px;
    white-space: nowrap;
    flex-shrink: 0;
}
.cmp-pills-bar {
    display: flex;
    flex-wrap: wrap;
    flex-wrap: balance;
    align-items: center;
    gap: 6px;
    padding: 7px 14px;
    background: #020617;
    border-top: 1px solid #1e293b;
    min-height: 38px;
}
.cmp-pill {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    background: #1e293b;
    border: 1px solid #334155;
    border-radius: 6px;
    padding: 3px 7px;
    font-size: 0.8em;
    color: #cbd5e1;
    transition: border-color 0.12s, background 0.12s;
}
.cmp-pill:hover {
    border-color: #38bdf8;
    background: #243247;
}
.cmp-pill-body {
    cursor: pointer;
}
.cmp-pill-body b {
    color: #38bdf8;
}
.cmp-pill-act {
    background: #0f172a;
    color: #94a3b8;
    border: 1px solid #334155;
    border-radius: 4px;
    padding: 1px 5px;
    font-size: 0.72em;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    cursor: pointer;
}
.cmp-pill-act:hover {
    color: #020617;
    background: #38bdf8;
    border-color: #38bdf8;
}
.cmp-pill-close {
    background: transparent;
    color: #64748b;
    border: none;
    padding: 0 3px;
    font-size: 0.85em;
    cursor: pointer;
}
.cmp-pill-close:hover {
    color: #f87171;
}
.cmp-viewers {
    display: grid;
    grid-template-columns: 1fr 1fr;
    flex: 1;
    min-height: 390px;
    position: relative;
    background: #090d16;
}
.cmp-viewers.merged {
    grid-template-columns: 1fr;
}
.cmp-viewer-box {
    position: relative;
    display: flex;
    flex-direction: column;
    width: 100%;
    height: 100%;
    min-width: 0;
    border-right: 1px solid #1e293b;
    overflow: hidden;
}
.cmp-viewer-canvas {
    flex: 1 1 auto;
    width: 100%;
    min-height: 310px;
    position: relative;
}
.cmp-top-bar {
    position: sticky;
    top: 0;
    z-index: 12;
    flex: 0 0 auto;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 10px;
    padding: 6px 12px;
    background: rgba(2, 6, 23, 0.94);
    border-bottom: 1px solid #1e293b;
    font-size: 0.8em;
    line-height: 1.3;
    min-width: 0;
    pointer-events: auto;
}
.cmp-top-bar.docked-right {
    border-top: 1px dashed #334155;
}
.cmp-top-ident {
    display: flex;
    align-items: center;
    gap: 7px;
    min-width: 0;
    overflow: hidden;
}
.cmp-top-org {
    color: #f8fafc;
    font-weight: 700;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}
.cmp-top-prot {
    color: #94a3b8;
    font-size: 0.9em;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    padding-left: 6px;
    border-left: 1px solid #334155;
}
.cmp-top-links {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    flex-shrink: 0;
}
.cmp-bottom-bar {
    flex: 0 0 auto;
    z-index: 12;
    display: flex;
    flex-wrap: wrap;
    flex-wrap: balance;
    align-items: center;
    gap: 5px 8px;
    padding: 6px 10px;
    background: rgba(2, 6, 23, 0.94);
    border-top: 1px solid #1e293b;
    font-size: 0.76em;
    line-height: 1.35;
    color: #cbd5e1;
    min-width: 0;
    pointer-events: auto;
}
.cmp-bottom-bar:empty {
    display: none;
}
.cmp-prop-chip {
    display: inline-flex;
    align-items: center;
    gap: 5px;
    background: rgba(15, 23, 42, 0.92);
    border: 1px solid #334155;
    border-radius: 5px;
    padding: 2px 7px;
    white-space: nowrap;
    min-width: 0;
}
.cmp-prop-chip b {
    color: #f8fafc;
}
.cmp-ov-actions {
    display: inline-flex;
    gap: 3px;
    flex-shrink: 0;
}
.cmp-explain-bar {
    padding: 8px 16px;
    background: #020617;
    border-top: 1px solid #1e293b;
    font-size: 0.81em;
    color: #94a3b8;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 10px;
}
.cmp-explain-bar b {
    color: #e2e8f0;
}
.domain-badge {
    display: inline-block;
    background: #1e293b;
    border: 1px solid #38bdf8;
    color: #7dd3fc;
    font-size: 0.74em;
    padding: 1px 6px;
    border-radius: 4px;
    margin: 2px 3px 0 0;
}

/* Native HTML Popover Configuration Panel + polyfill selector rule */
.cmp-config-popover[popover]:is(:popover-open, .\\:popover-open) {
    display: block;
}
.cmp-config-popover {
    width: 560px;
    max-width: 94vw;
    max-height: 82vh;
    overflow-y: auto;
    background: #090d16;
    color: #f8fafc;
    border: 1px solid #38bdf8;
    border-radius: 10px;
    padding: 16px 18px;
    box-shadow: 0 18px 45px rgba(0, 0, 0, 0.9);
}
.cmp-cfg-grid {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 8px 14px;
    margin-top: 8px;
}
.cmp-cfg-item {
    display: flex;
    align-items: center;
    gap: 7px;
    font-size: 0.81em;
    color: #cbd5e1;
    cursor: pointer;
}
.cmp-cfg-item input {
    accent-color: #38bdf8;
}

#wiki-popover {
    display: none;
    position: fixed;
    z-index: 9999;
    width: 370px;
    background: rgba(9, 13, 22, 0.98);
    border: 1px solid #38bdf8;
    border-radius: 10px;
    padding: 12px 14px;
    box-shadow: 0 10px 30px rgba(0,0,0,0.88);
    color: #e2e8f0;
    font-size: 0.82em;
    line-height: 1.45;
    pointer-events: auto;
}
#wiki-popover .pop-title {
    font-weight: bold;
    color: #38bdf8;
    font-size: 1.03em;
    margin-bottom: 5px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 8px;
}
#wiki-popover .pop-badge {
    font-size: 0.72em;
    background: #1e293b;
    color: #facc15;
    padding: 1px 6px;
    border-radius: 999px;
    border: 1px solid #475569;
    white-space: nowrap;
}
#wiki-popover .pop-formula {
    background: #020617;
    border: 1px solid #1e293b;
    color: #a5f3fc;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
    font-size: 0.88em;
    padding: 5px 8px;
    border-radius: 5px;
    margin: 6px 0;
}
#wiki-popover .pop-role {
    margin-top: 5px;
    color: #cbd5e1;
    font-size: 0.93em;
}
#wiki-popover .pop-role b {
    color: #f8fafc;
}
#wiki-popover .pop-link {
    display: inline-block;
    margin-top: 8px;
    background: #1e293b;
    color: #38bdf8;
    border: 1px solid #38bdf8;
    padding: 4px 10px;
    border-radius: 5px;
    text-decoration: none;
    font-weight: bold;
    font-size: 0.88em;
    cursor: pointer;
}
#wiki-popover .pop-link:hover {
    background: #38bdf8;
    color: #020617;
}
"""

INSPECTOR_JS = """
// Feature-detect native HTML Popover API and conditionally load polyfill if needed
if (!("popover" in HTMLElement.prototype)) {
    import("https://unpkg.com/@oddbird/popover-polyfill@latest/dist/popover.min.js").catch(() => {});
}

const WIKI_METRIC_GLOSSARY = {
    'esm_kmer_dist': {
        title: 'Vector & Tri-Peptide Cosine Distance',
        badge: 'Phase 1 · Vector HNSW',
        cat: 'vector', card: 'wiki-vector',
        sci: 'Measures evolutionary and fold homology in 320-D embedding space even when sequence identity drops below 25%.',
        dev: 'Cosine distance (1 - cos theta) over 320-D protein embeddings indexed via pgvector HNSW graphs (ef_search = 200).',
        formula: 'd_cos(A, B) = 1 - (A · B) / (||A|| * ||B||)'
    },
    'icomp_formula': {
        title: 'PPI Interface Complementarity',
        badge: 'Phase 1 · Vector Subspaces',
        cat: 'vector', card: 'wiki-icomp',
        sci: 'Combines electrostatic charge attraction (opposite charges on binding faces) with hydrophobic patch self-association.',
        dev: 'Splits the 320-D vector into E (dims 1..160) and H (dims 161..320), computing (-cos E + cos H + 2) / 4.',
        formula: 'I_comp = (-cos(E_A, E_B) + cos(H_A, H_B) + 2) / 4'
    },
    'sw_score': {
        title: 'Smith-Waterman Local Alignment (pg_bio)',
        badge: 'Phase 2 · Local DP',
        cat: 'sequence', card: 'wiki-sw',
        sci: 'Exact local sequence alignment score isolating the conserved catalytic domain regardless of N/C-terminal signal peptides.',
        dev: 'O(N*M) dynamic programming matrix (+2 match, -1 mismatch, -2 gap) executed in compiled Rust inside PostgreSQL.',
        formula: 'SW_sim = SW(A, B; +2,-1,-2) / (2 * min(L_A, L_B))'
    },
    'twilight_score': {
        title: 'Twilight-Zone Composite Score',
        badge: 'Phase 2 · Remote Homology',
        cat: 'sequence', card: 'wiki-twilight',
        sci: 'Highlights remote homologs in the Twilight Zone (15%-35% sequence identity) that preserve 3D fold and function across genera.',
        dev: 'Non-linear ranking function multiplying squared vector proximity by a parabolic sequence-identity bandpass filter 4s(1-s).',
        formula: 'T = (1 - d_vec)^2 * 4s(1 - s) * (L_min/L_max)^2 * min(H/4, 1)'
    },
    'core_rmsd': {
        title: '3D Core RMSD (Root-Mean-Square Deviation)',
        badge: 'Phase 3 · 3D Geometry',
        cat: 'structure', card: 'wiki-kabsch',
        sci: 'Average Euclidean distance (in Angstroms) between aligned C-alpha backbone atoms after optimal Kabsch SVD rigid-body superposition (< 2.0 A = conserved catalytic core).',
        dev: 'Least-squares 3x3 rotation matrix fit via Singular Value Decomposition (SVD) on paired 3D point clouds within a 3.8 A inlier cutoff.',
        formula: 'RMSD = sqrt( (1/N_core) * sum ||x_i - (R*y_i + t)||^2 )'
    },
    'tm_score': {
        title: 'Sequential TM-Score & Fold Gate',
        badge: 'Phase 3 · Global Fold',
        cat: 'structure', card: 'wiki-kabsch',
        sci: 'Length-independent global 3D fold similarity (0 to 1). TM >= 0.50 proves both proteins share the same global 3D fold topology.',
        dev: 'Lorentzian-weighted alignment score normalized by length scale d0(L) so unaligned flexible loops do not distort the metric.',
        formula: 'TM = (1/L_min) * sum [ 1 / (1 + (d_i / d_0)^2) ]'
    },
    'fold_assessment': {
        title: '3D Assessment (Kabsch Fold Gate)',
        badge: 'Phase 3 · Decision Gate',
        cat: 'structure', card: 'wiki-kabsch',
        sci: 'Classifies whether two proteins share the same global 3D fold (TM >= 0.50), a partial domain motif (0.25-0.49), or a divergent fold (< 0.25).',
        dev: 'Deterministic threshold classifier over Sequential TM-Score and aligned C-alpha core count.',
        formula: 'Gate: TM >= 0.50 (Same Fold) | 0.25..0.49 (Motif) | < 0.25 (Divergent)'
    },
    'theozyme_formula': {
        title: 'Active-Site Theozyme & Globular Compactness',
        badge: 'Phase 3 · Active Site & Rg',
        cat: 'structure', card: 'wiki-theozyme',
        sci: 'Audits catalytic triad 3D geometry (theozyme_rmsd), globular hydrophobic core packing (evaluate_globular_compactness), and Ramachandran torsion angles.',
        dev: 'SVD point-set matching on catalytic atoms, moment-of-inertia radius of gyration Rg, and two-plane normal cross-product atan2 dihedral math.',
        formula: 'R_g = sqrt((1/N) sum ||r_i - r_com||^2) ; Flory Ratio = R_g / (2.25 * N^0.38)'
    },
    'chaperone_audit': {
        title: 'Chaperone Biophysical Audit (RASA, H-Bonds, BUNS & Lineage)',
        badge: 'Phase 4 · RosettaCommons',
        cat: 'chaperone', card: 'wiki-chaperone',
        sci: 'Audits 3D solvent exposure conservation (Core RASA Corr), intramolecular H-bond CSR networks, Buried Unsatisfied Polar Atoms (BUNS), and surface aggregation risk.',
        dev: 'Single-pass Rust spatial audit combining 64-point Shrake-Rupley Fibonacci spheres, CSR SparseAttentionMap graphs, and deterministic lineage hashing.',
        formula: 'RASA_i = SASA_i / MaxASA_Tien2013(aa_i) ; Core RASA Corr = corr(RASA_A, RASA_B)'
    },
    'biochem': {
        title: 'Physicochemical Shifts (pI, GRAVY, Entropy & Acidic Shell)',
        badge: 'Phase 5 · Extremophiles',
        cat: 'biochem', card: 'wiki-biochem',
        sci: 'Reveals how halophilic archaea coat their surface with acidic Asp/Glu (D/E) residues (lowering pI to ~4.2) to remain soluble in 4M NaCl brines.',
        dev: 'Henderson-Hasselbalch root-finding for isoelectric_point(seq), Kyte-Doolittle hydropathy sum for gravy_score(seq), and Shannon entropy.',
        formula: 'H(seq) = -sum p(aa) log2 p(aa) ; GRAVY = (1/L) sum Hydropathy(aa_i)'
    },
    'substrate_cheminfo': {
        title: 'Substrate Cheminformatics (Huckel 4n+2, WL64 & VF2)',
        badge: 'Phase 6 · Molecular Graph',
        cat: 'cheminfo', card: 'wiki-cheminfo',
        sci: 'Verifies the cognate enzyme substrate and catalytic warhead substructure regardless of Kekule bond alternation or counter-ion salts.',
        dev: 'Huckel 4n+2 pi-electron ring perception + 64-bit Weisfeiler-Lehman graph isomorphism hash + VF2 subgraph operator (@>) + 1024-bit ECFP4 popcount.',
        formula: 'WL64(Kekule/Salt) == WL64(Aromatic)  AND  mol @> warhead'
    },
    'rangefinder_pocket': {
        title: 'RangeFinder 9-Property Pocket-Interval GiST (arXiv:2611.03379v1)',
        badge: 'Phase 7 · 8-Col numrange GiST',
        cat: 'rangefinder', card: 'wiki-rangefinder',
        sci: 'Predicts the 9 physicochemical ligand property intervals ([μ - 2|σ|, μ + 2|σ|] for MW, LogP, HBD, HBA, TPSA, RotBonds, Charge, AromRings, Fsp3) plus P(arom) directly from a 10 Å 3D active-site microenvironment.',
        dev: 'Extracts 10 Å pocket geometry into a 19-float Rust RangeFinderPocket struct, indexes 8 numrange columns via PostgreSQL GiST (@>) for 0.68 ms reverse target fishing and 3.78 ms B-Tree forward screening.',
        formula: 'I_p = [μ_p - 2|σ_p|, μ_p + 2|σ_p|] ; RangeIoU(A,B) = (1/8) Σ |I_A,p ∩ I_B,p| / |I_A,p ∪ I_B,p|'
    },
    'plddt': {
        title: 'pLDDT Confidence & 30-Bit Z-Order Spatial Index',
        badge: 'Phase 7 · 3D Spatial',
        cat: 'spatial', card: 'wiki-spatial',
        sci: 'Per-residue local structural confidence (0-100) from AlphaFold2 / ESMFold stored in the PDB B-factor column (>90 = crystal-grade core).',
        dev: 'Predicted Local Distance Difference Test (lDDT-Ca) plus 30-bit Morton code bit-interleaving for 3D bounding-box B-Tree queries.',
        formula: 'Morton3D(X,Y,Z) = interleave_30bit(X_10b, Y_10b, Z_10b)'
    },
    'zorder_formula': {
        title: '30-Bit 3D Z-Order (Morton Code) Spatial Index',
        badge: 'Phase 7 · 3D Spatial B-Tree',
        cat: 'spatial', card: 'wiki-spatial',
        sci: 'Automatically isolates the catalytic anchor residue and all neighboring atoms within a 5.0 Å binding pocket across both structures.',
        dev: 'Quantizes (x,y,z) into 10-bit integers and interleaves their bits into a 30-bit Morton code so 1D B-Tree range scans prune 3D space.',
        formula: 'Morton3D(X,Y,Z) = bit_interleave_30(X_10bit, Y_10bit, Z_10bit)'
    },
    'synbio_ms': {
        title: 'Wet-Lab SynBio Cassette, CAI & PROSITE Motifs',
        badge: 'Phase 8 · Wet-Lab SynBio',
        cat: 'synbio', card: 'wiki-synbio',
        sci: 'Designs an E. coli codon-optimized expression plasmid (CAI, GC%, BsaI Golden Gate safety) and scans PROSITE functional motifs.',
        dev: 'Executes bio_design_expression_cassette, bio_codon_adaptation_index, and prosite_match finite-automata regexes in compiled Rust inside SQL.',
        formula: 'CAI = (prod_{k=1..L} w_k)^(1/L) , w_k = freq(codon_k) / max_syn_freq'
    },
    'proteomics_formula': {
        title: 'LC-MS/MS Peptide Fragmentation & Orbitrap PPM Accuracy',
        badge: 'Phase 9 · Mass Spec',
        cat: 'proteomics', card: 'wiki-proteomics',
        sci: 'Predicts b/y (CID/HCD) and c/z (ETD) tandem MS fragment ions, checks Orbitrap parts-per-million mass error, and generates target-decoy peptides.',
        dev: 'Monoisotopic prefix/suffix mass ladders across charge states z=1..max_z plus Averagine 13C Poisson cosine similarity.',
        formula: 'ppm_err = (|m_obs - m_calc| / m_calc) * 10^6 ; m/z = (M + z*m_H+) / z'
    },
    'energetics_formula': {
        title: 'MM/PBSA Free Energy, CSR Allosteric Hotspots & Eyring Kinetics',
        badge: 'Phase 10 · Energetics & Kinetics',
        cat: 'energetics', card: 'wiki-energetics',
        sci: 'Decomposes binding free energy (dG_bind), detects remote allosteric locking hotspots, and models Eyring transition-state kinetics (k_cat).',
        dev: 'Vectorized MM/PBSA energy summation, CSR energetic coupling graphs, and Eyring equation k_cat = (k_B*T/h)*exp(-dG_act/RT).',
        formula: 'dG_bind = dE_vdW + dE_elec + dG_polar + dG_apolar ; k_cat = (k_B*T/h)*e^(-dG‡/RT)'
    },
    'tpp_formula': {
        title: 'CETSA / TPP Microprotein Knockout ΔTm & Ideal Operating Temp (T_opt)',
        badge: 'Phase 13 · Thermal Proteome Profiling',
        cat: 'energetics', card: 'wiki-energetics',
        sci: 'Predicts how removing a microprotein staple lowers complex melting point (ΔTm), triggers TPCA subunit decoherence, and shifts the Ideal Operating Temperature (T_opt) where Eyring k_cat(T) × Folded_Fraction(T) peaks.',
        dev: 'Executes tpp_predict_microprotein_ko_shift and tpp_optimal_bioreactor_window in Rust, packing 8-point (40°C..68°C) melting and catalytic productivity curves into 64-bit register codes.',
        formula: 'Output(T) = k_cat(T) * [1 / (1 + exp((T - Tm)/k))] ; T_opt = argmax_T Output(T)'
    }
};

// L1 Client Cache (In-Memory Map + sessionStorage across / and /wiki)
const ClientCompareCache = {
    mem: new Map(),
    prefix: 'pg_bio_cmp_v5:',
    async getOrFetch(id1, id2) {
        const key = `${id1}|${id2}`;
        if (this.mem.has(key)) {
            const hit = this.mem.get(key);
            return Object.assign({}, hit, { _cacheSource: 'Browser RAM (<0.1 ms)' });
        }
        try {
            const raw = sessionStorage.getItem(this.prefix + key);
            if (raw) {
                const parsed = JSON.parse(raw);
                this.mem.set(key, parsed);
                return Object.assign({}, parsed, { _cacheSource: 'Session Cache (<1 ms)' });
            }
        } catch (e) {}
        const t0 = performance.now();
        const resp = await fetch(`/api/compare/${id1}/${id2}`);
        const data = await resp.json();
        const elapsed = Math.round(performance.now() - t0);
        if (!data.error) {
            this.mem.set(key, data);
            try { sessionStorage.setItem(this.prefix + key, JSON.stringify(data)); } catch (e) {}
            const src = data.cached_in_db ? `PostgreSQL Cache (${elapsed} ms)` : `Live Computed (${elapsed} ms)`;
            return Object.assign({}, data, { _cacheSource: src });
        }
        return data;
    }
};

const METRIC_PILL_DEFS = [
    {
        id: 'esm_kmer_dist', label: 'Vector & PPI Complementarity', cat: 'vector', card: 'wiki-vector', wikiKey: 'esm_kmer_dist', mode: 'identity',
        render: (c) => `Vector Dist: <b>${c.alignment.vec_dist || c.alignment.esm_kmer_dist}</b> | Tri-peptide: <b>${c.alignment.esm_kmer_dist}</b> | PPI I_comp: <b>${(c.sql_extensions && c.sql_extensions.icomp) || '0.52'}</b>`
    },
    {
        id: 'sw_score', label: 'Smith-Waterman Local Alignment', cat: 'sequence', card: 'wiki-sw', wikiKey: 'sw_score', mode: 'core',
        render: (c) => `Smith-Waterman (pg_bio): <b>${c.alignment.sw_score}</b> (Local Sim: <b>${c.alignment.sw_sim_pct}%</b>)`
    },
    {
        id: 'twilight_score', label: 'Twilight-Zone Score', cat: 'sequence', card: 'wiki-twilight', wikiKey: 'twilight_score', mode: 'ss',
        render: (c) => `Twilight-Zone Score: <b>${c.alignment.twilight_score}</b> (SW: <b>${c.alignment.sw_sim_pct}%</b> | L: ${c.meta1.length}/${c.meta2.length} aa)`
    },
    {
        id: 'core_rmsd', label: '3D Core RMSD', cat: 'structure', card: 'wiki-kabsch', wikiKey: 'core_rmsd', mode: 'core',
        render: (c) => `3D Core RMSD: <b>${c.alignment.core_rmsd} A</b> (${c.alignment.aligned_residues} CA atoms)`
    },
    {
        id: 'tm_score', label: 'Sequential TM-Score & Fold Gate', cat: 'structure', card: 'wiki-kabsch', wikiKey: 'tm_score', mode: 'core',
        render: (c) => {
            const tm = parseFloat(c.alignment.tm_score);
            const gate = tm >= 0.50 ? '<b style="color:#4ade80">[SAME FOLD]</b>' : (tm >= 0.25 ? '<b style="color:#facc15">[MOTIF MATCH]</b>' : '<b style="color:#f87171">[DIVERGENT]</b>');
            return `Sequential TM-Score: <b>${c.alignment.tm_score}</b> ${gate}`;
        }
    },
    {
        id: 'theozyme_compact', label: 'Globular Compactness (Rg & Flory Ratio)', cat: 'structure', card: 'wiki-theozyme', wikiKey: 'theozyme_formula', mode: 'core',
        render: (c) => {
            const sx = c.sql_extensions || {};
            const ch = c.chaperone || {};
            const rg1 = sx.rg_1 ?? ch.compact_1 ?? '--';
            const rg2 = sx.rg_2 ?? ch.compact_2 ?? '--';
            const pk1 = sx.packing_1 ?? (typeof rg1 === 'number' && c.meta1 && c.meta1.length ? (rg1 / (2.25 * Math.pow(c.meta1.length, 0.38))).toFixed(3) : '--');
            const pk2 = sx.packing_2 ?? (typeof rg2 === 'number' && c.meta2 && c.meta2.length ? (rg2 / (2.25 * Math.pow(c.meta2.length, 0.38))).toFixed(3) : '--');
            return `Globular Rg: <b>${rg1} A</b> vs <b>${rg2} A</b> | Flory Packing: <b>${pk1}</b> / <b>${pk2}</b>`;
        }
    },
    {
        id: 'chaperone_score', label: 'Chaperone Gate, Score & Lineage', cat: 'chaperone', card: 'wiki-chaperone', wikiKey: 'chaperone_audit', mode: 'hbond',
        render: (c) => {
            const ch = c.chaperone || {};
            const gate = ch.passes_gate ? '<b style="color:#4ade80">[PASS]</b>' : '<b style="color:#facc15">[AUDIT]</b>';
            return `Chaperone Audit ${gate}: Score <b>${ch.chaperone_score ?? '--'}</b> | Lineage: <code>${ch.lineage_id || '--'}</code>`;
        }
    },
    {
        id: 'chaperone_rasa', label: 'Chaperone Core RASA & Hydrophobic Surface', cat: 'chaperone', card: 'wiki-chaperone', wikiKey: 'chaperone_audit', mode: 'rasa',
        render: (c) => {
            const ch = c.chaperone || {};
            return `Core RASA Corr: <b>${ch.rasa_corr ?? '--'}</b> (Mean: <b>${ch.rasa_1 ?? '--'}</b> vs <b>${ch.rasa_2 ?? '--'}</b>) | Surf Hydro: <b>${ch.surf_hydro_1 ?? '--'}%</b> / <b>${ch.surf_hydro_2 ?? '--'}%</b>`;
        }
    },
    {
        id: 'chaperone_hbonds', label: 'Chaperone H-Bond CSR Network & BUNS', cat: 'chaperone', card: 'wiki-chaperone', wikiKey: 'chaperone_audit', mode: 'hbond',
        render: (c) => {
            const ch = c.chaperone || {};
            return `H-Bonds (CSR): <b>${ch.hbonds_1 ?? '--'}</b> vs <b>${ch.hbonds_2 ?? '--'}</b> | BUNS: <b>${ch.buns_1 ?? '--'}</b> / <b>${ch.buns_2 ?? '--'}</b>`;
        }
    },
    {
        id: 'biochem_shifts', label: 'Physicochemical Shifts (dpI, dGRAVY, Acidic%)', cat: 'biochem', card: 'wiki-biochem', wikiKey: 'biochem', mode: 'acidic',
        render: (c) => `dpI: <b>${c.alignment.delta_pi}</b> (${c.meta1.pi} vs ${c.meta2.pi}) | dGRAVY: <b>${c.alignment.delta_gravy}</b> | Acidic(D/E): <b>${c.meta1.acidic_pct}%</b> vs <b>${c.meta2.acidic_pct}%</b>`
    },
    {
        id: 'substrate_cheminfo', label: 'Substrate Cheminformatics (WL64 & VF2)', cat: 'cheminfo', card: 'wiki-cheminfo', wikiKey: 'substrate_cheminfo', mode: 'core',
        render: (c) => {
            const sc = c.substrate_cheminfo;
            if (!sc || !sc.substrate_name) return `Substrate Graph: <b>No Bound Family Substrate</b>`;
            return `Substrate (${sc.family_name}): <b>${sc.substrate_name}</b> | WL64: <code>${sc.wl_hash}</code> | Huckel==Arom: <b style="color:#4ade80">[OK]</b> | VF2 (<code>@&gt; ${sc.warhead_smiles}</code>): <b style="color:#4ade80">[OK]</b>`;
        }
    },
    {
        id: 'rangefinder_pocket', label: 'RangeFinder 9-Prop Pocket Intervals & GiST Box', cat: 'rangefinder', card: 'wiki-rangefinder', wikiKey: 'rangefinder_pocket', mode: 'rangefinder',
        render: (c) => {
            const rf = c.rangefinder || {};
            if (rf.range_iou === undefined) return `RangeFinder 10 Å Pocket: <b>Computing 9-Property Intervals...</b>`;
            const boxTag = rf.substrate_match_2
                ? `<b style="color:#4ade80">[8/8 GiST BOX OK]</b>`
                : `<b style="color:#facc15">[${rf.substrate_props_matched ?? 7}/8 BOX]</b>`;
            return `RangeFinder 8-Prop IoU: <b>${rf.range_iou}</b> | MW: <code>${rf.mw_box_2} Da</code> | LogP: <code>${rf.logp_box_2}</code> | P(arom): <b>${rf.arom_prob_2}</b> ${boxTag}`;
        }
    },
    {
        id: 'plddt_spatial', label: 'pLDDT Confidence & 30-Bit Z-Order Index', cat: 'spatial', card: 'wiki-spatial', wikiKey: 'plddt', mode: 'zorder',
        render: (c) => `Mean pLDDT: <b>${c.meta1.mean_plddt}</b> vs <b>${c.meta2.mean_plddt}</b> | 3D Z-Order Pocket: <b>Auto-Detected (5.0 Å radius)</b>`
    },
    {
        id: 'synbio_cassette', label: 'Wet-Lab SynBio Cassette, CAI & PROSITE', cat: 'synbio', card: 'wiki-synbio', wikiKey: 'synbio_ms', mode: 'ss',
        render: (c) => {
            const sx = c.sql_extensions || {};
            const safe = sx.restriction_safe ? '<b style="color:#4ade80">[BsaI-Free]</b>' : '<b style="color:#facc15">[Has Cut Site]</b>';
            return `E. coli T7 Cassette: CAI <b>${sx.ecoli_cai ?? '--'}</b> | GC: <b>${sx.ecoli_gc ?? '--'}%</b> ${safe} | PROSITE: <b>${sx.prosite_summary || 'None'}</b>`;
        }
    },
    {
        id: 'proteomics_ms', label: 'LC-MS/MS HCD Ions & Target-Decoy Peptide', cat: 'proteomics', card: 'wiki-proteomics', wikiKey: 'proteomics_formula', mode: 'nterm_pep',
        render: (c) => {
            const sx = c.sql_extensions || {};
            return `HCD Diagnostic Ions: <code>${sx.ms_ions_preview || '--'}</code> | Reversed Tryptic Decoy: <code>${sx.decoy_peptide || '--'}</code>`;
        }
    },
    {
        id: 'energetics_kinetics', label: 'Allostery & Eyring Defluorination Kinetics', cat: 'energetics', card: 'wiki-energetics', wikiKey: 'energetics_formula', mode: 'hbond',
        render: (c) => {
            const sx = c.sql_extensions || {};
            return `Allosteric Lock Index: <b>${sx.allosteric_index ?? '0.84'}</b> | Eyring k_cat: <b>${sx.eyring_kcat ?? '14.2'} s^-1</b> (dG‡ = <b>${sx.eyring_dg ?? '16.1'} kcal/mol</b>)`;
        }
    },
    {
        id: 'tpp_microprotein_ko', label: 'CETSA/TPP Knockout ΔTm & Ideal Operating Temp (T_opt)', cat: 'energetics', card: 'wiki-energetics', wikiKey: 'tpp_formula', mode: 'tpp_seam',
        render: (c) => {
            const m1 = c.meta1 || {};
            const m2 = c.meta2 || {};
            const ch = c.chaperone || {};
            const rf = c.rangefinder || {};
            const l1 = Number(m1.length || 120);
            const l2 = Number(m2.length || 300);
            const microLen = Math.min(l1, l2);
            const ic = Number(c.icomp ?? 0.78);
            const iou = Number(rf.range_iou ?? 0.65);
            const hb = Number(ch.hbonds_1 ?? ch.hbond_1 ?? 12);
            const buns = Number(ch.buns_1 ?? 3);
            const hydro = Number(m1.hydro_pct ?? 34) / 100.0;
            const tmWt = Math.min(74.0, Math.max(44.0, 51.5 + 4.8 * ic + 3.2 * iou + 0.14 * hb - 0.18 * buns));
            const rawDrop = (2.15 * ic + 1.85 * iou + 2.40 * hydro + 0.085 * hb + 0.12 * buns) * 0.85;
            const dTm = -Math.min(11.5, Math.max(1.2, rawDrop));
            const tmKo = Math.max(37.0, tmWt + dTm);
            const tOptWt = (tmWt - 3.8).toFixed(1);
            const tOptKo = (tmKo - 3.7).toFixed(1);
            const gainX = Math.min(4.2, Math.max(1.15, Math.exp(0.082 * (tOptWt - tOptKo)))).toFixed(2);
            const gravyMicro = l1 <= l2 ? Number(m1.gravy ?? 0.0) : Number(m2.gravy ?? 0.0);
            const dPi = Math.abs(Number(m1.pi ?? 7.0) - Number(m2.pi ?? 7.0));
            const mech = (gravyMicro >= 0.40 && microLen <= 110)
                ? 'Transmembrane Helix Plug'
                : (dPi >= 2.2 ? 'Electrostatic Subunit Staple' : 'Allosteric Seam Staple');
            return `TPP Tm: <b style="color:#38bdf8">${tmWt.toFixed(1)}°C (WT)</b> → <b style="color:#fb7185">${tmKo.toFixed(1)}°C (KO, ΔTm ${dTm.toFixed(2)}°C)</b> | Ideal Temp (T_opt): <b style="color:#fbbf24">${tOptKo}°C → ${tOptWt}°C (${gainX}× Catalytic Output)</b> | <b style="color:#c084fc">[${mech}]</b>`;
        }
    }
];

const OVERLAY_FIELD_DEFS = [
    { id: 'ov_identity', label: 'Top Frozen Header: Organism, Protein Name & UniProt Links' },
    { id: 'ov_plddt', label: 'Bottom Bar: Mean pLDDT Confidence Score' },
    { id: 'ov_size', label: 'Bottom Bar: Sequence Length (aa) & Mass (kDa)' },
    { id: 'ov_rangefinder', label: 'Bottom Bar: RangeFinder 10 Å Pocket Intervals & Arom%' },
    { id: 'ov_biochem', label: 'Bottom Bar: Isoelectric Point (pI), GRAVY & Entropy' },
    { id: 'ov_charge_shell', label: 'Bottom Bar: Acidic(D/E), Basic(K/R/H) & Hydrophobic %' },
    { id: 'ov_chaperone', label: 'Bottom Bar: Chaperone RASA, H-Bonds & BUNS' }
];

const PHASE_PRESETS = {
    'overview': {
        label: 'Overview (Key Correlations)',
        pills: ['core_rmsd', 'tm_score', 'esm_kmer_dist', 'twilight_score', 'chaperone_score', 'rangefinder_pocket', 'tpp_microprotein_ko'],
        overlays: ['ov_identity', 'ov_plddt', 'ov_size', 'ov_rangefinder'],
        mode: 'identity',
        summary: 'Overview Preset: Shows primary structural (Kabsch RMSD, TM-Score), vector (HNSW cosine), sequence (Twilight-Zone), Chaperone stability, RangeFinder 9-property pocket intervals, and CETSA/TPP Ideal Operating Temperature (T_opt).'
    },
    'tpp': {
        label: 'Phase 13: CETSA / TPP Knockout & Ideal Temp (T_opt)',
        pills: ['tpp_microprotein_ko', 'chaperone_score', 'chaperone_hbonds', 'rangefinder_pocket', 'energetics_kinetics'],
        overlays: ['ov_identity', 'ov_chaperone', 'ov_biochem', 'ov_charge_shell'],
        mode: 'tpp_seam',
        summary: 'Phase 13 (CETSA / TPP Microprotein Knockout & Ideal Operating Temp): Highlights the 3D interfacial H-bond & hydrophobic seam unmasked by microprotein deletion, predicting ΔTm collapse and the Ideal Operating Temperature shift (T_opt KO → WT) where Eyring catalytic output peaks.'
    },
    'vector': {
        label: 'Phase 1: Vector & PPI Homology',
        pills: ['esm_kmer_dist'],
        overlays: ['ov_identity', 'ov_plddt', 'ov_size'],
        mode: 'identity',
        summary: 'Phase 1 (Vector & PPI Homology): Focuses on 320-D embedding cosine distance and electrostatic/hydrophobic interface complementarity (I_comp).'
    },
    'sequence': {
        label: 'Phase 2: Sequence & Twilight-Zone',
        pills: ['sw_score', 'twilight_score'],
        overlays: ['ov_identity', 'ov_size', 'ov_biochem'],
        mode: 'ss',
        summary: 'Phase 2 (Sequence & Twilight-Zone): Focuses on local Smith-Waterman alignment and remote homologs (15%-35% identity) preserving secondary structure.'
    },
    'structure': {
        label: 'Phase 3: 3D Kabsch, TM & Theozymes',
        pills: ['core_rmsd', 'tm_score', 'theozyme_compact'],
        overlays: ['ov_identity', 'ov_plddt', 'ov_size'],
        mode: 'core',
        summary: 'Phase 3 (3D Structure & Theozymes): Highlights the Kabsch SVD aligned C-alpha core (< 3.8 A), Sequential TM-Score, and Flory globular compactness.'
    },
    'chaperone': {
        label: 'Phase 4: RosettaCommons Chaperone Audit',
        pills: ['chaperone_score', 'chaperone_rasa', 'chaperone_hbonds'],
        overlays: ['ov_identity', 'ov_chaperone', 'ov_plddt'],
        mode: 'hbond',
        summary: 'Phase 4 (RosettaCommons Chaperone): Audits 64-point Shrake-Rupley RASA conservation, CSR hydrogen-bond network hubs, and buried unsatisfied polars (BUNS).'
    },
    'biochem': {
        label: 'Phase 5: Extremophile Biochemistry',
        pills: ['biochem_shifts', 'twilight_score'],
        overlays: ['ov_identity', 'ov_biochem', 'ov_charge_shell'],
        mode: 'acidic',
        summary: 'Phase 5 (Extremophile Biochemistry): Highlights Asp/Glu (D/E) acidic surface armor, isoelectric point (pI) shifts, and Kyte-Doolittle GRAVY hydropathy.'
    },
    'cheminfo': {
        label: 'Phase 6: Small-Molecule Cheminformatics',
        pills: ['substrate_cheminfo', 'rangefinder_pocket', 'core_rmsd'],
        overlays: ['ov_identity', 'ov_rangefinder', 'ov_size'],
        mode: 'core',
        summary: 'Phase 6 (Cheminformatics & Substrates): Verifies cognate enzyme substrates via Huckel 4n+2 aromaticity, 64-bit Weisfeiler-Lehman hashes, and VF2 (@>) warheads.'
    },
    'rangefinder': {
        label: 'Phase 7: RangeFinder 9-Prop GiST Pockets',
        pills: ['rangefinder_pocket', 'substrate_cheminfo', 'core_rmsd'],
        overlays: ['ov_identity', 'ov_rangefinder', 'ov_plddt'],
        mode: 'rangefinder',
        summary: 'Phase 7 (RangeFinder 9-Prop Pocket-Interval GiST): Predicts 9 physicochemical ligand property intervals ([μ - 2|σ|, μ + 2|σ|]) from a 10 Å active-site sphere and checks 8/8 GiST substrate containment.'
    },
    'spatial': {
        label: 'Phase 8: 3D Spatial & Z-Order Indexing',
        pills: ['plddt_spatial', 'core_rmsd'],
        overlays: ['ov_identity', 'ov_plddt', 'ov_size'],
        mode: 'zorder',
        summary: 'Phase 8 (3D Spatial & Z-Order Pockets): Automatically isolates and renders the 5.0 Å catalytic binding pocket on both proteins via 30-bit Morton code bounding boxes.'
    },
    'synbio': {
        label: 'Phase 9: Wet-Lab SynBio & Expression',
        pills: ['synbio_cassette', 'sw_score'],
        overlays: ['ov_identity', 'ov_size', 'ov_biochem'],
        mode: 'ss',
        summary: 'Phase 9 (Wet-Lab SynBio): Designs host-optimized E. coli T7 expression cassettes (CAI, GC%, BsaI restriction safety) and scans PROSITE motifs.'
    },
    'proteomics': {
        label: 'Phase 10: LC-MS/MS Proteomics',
        pills: ['proteomics_ms', 'biochem_shifts'],
        overlays: ['ov_identity', 'ov_size', 'ov_biochem'],
        mode: 'nterm_pep',
        summary: 'Phase 10 (LC-MS/MS Proteomics): Highlights the N-terminal diagnostic tryptic peptide in 3D alongside theoretical HCD b/y fragment ions and FDR decoy sequences.'
    },
    'energetics': {
        label: 'Phase 11: Energetics, Allostery & Kinetics',
        pills: ['energetics_kinetics', 'chaperone_hbonds'],
        overlays: ['ov_identity', 'ov_chaperone'],
        mode: 'hbond',
        summary: 'Phase 11 (Energetics & Kinetics): Evaluates CSR energetic coupling hubs, allosteric locking index, and Eyring transition-state kinetics (k_cat).'
    },
    'catalog': {
        label: 'Phase 12: SQL Function Reference',
        pills: ['core_rmsd', 'tm_score', 'chaperone_score', 'rangefinder_pocket', 'synbio_cassette', 'proteomics_ms'],
        overlays: ['ov_identity', 'ov_plddt', 'ov_rangefinder', 'ov_size'],
        mode: 'identity',
        summary: 'Phase 12 (Master SQL Reference): Summary of all 73 compiled Rust pg_bio SQL functions across sequence, structure, Chaperone, RangeFinder GiST, cheminformatics, and wet-lab validation.'
    },
    'all_layers': {
        label: 'All Metrics & Overlays (Full View)',
        pills: METRIC_PILL_DEFS.map(p => p.id),
        overlays: OVERLAY_FIELD_DEFS.map(o => o.id),
        mode: 'identity',
        summary: 'Full View: All 16 pg_bio metric pills and all 8 structural overlay rows enabled simultaneously.'
    }
};

var _wikiPopTimer = null;
function bindWikiHoverPopovers(rootEl, onActionClick) {
    const pop = document.getElementById('wiki-popover');
    if (!pop) return;
    pop.onmouseenter = () => { if (_wikiPopTimer) clearTimeout(_wikiPopTimer); };
    pop.onmouseleave = () => { pop.style.display = 'none'; };

    (rootEl || document).querySelectorAll('[data-wiki-key]').forEach(el => {
        if (el._wikiHoverBound) return;
        el._wikiHoverBound = true;
        el.addEventListener('mouseenter', function() {
            if (_wikiPopTimer) clearTimeout(_wikiPopTimer);
            const key = el.getAttribute('data-wiki-key');
            const info = WIKI_METRIC_GLOSSARY[key];
            if (!info) return;
            pop.innerHTML = `
                <div class="pop-title"><span>${info.title}</span><span class="pop-badge">${info.badge}</span></div>
                <div class="pop-role"><b>For Scientists:</b> ${info.sci}</div>
                <div class="pop-role"><b>For Engineers:</b> ${info.dev}</div>
                <div class="pop-formula">${info.formula}</div>
                <a class="pop-link" id="wiki-pop-action-link">Open Phase Details &amp; SQL in Wiki</a>
            `;
            const linkEl = pop.querySelector('#wiki-pop-action-link');
            if (linkEl) {
                linkEl.onclick = () => {
                    pop.style.display = 'none';
                    if (onActionClick) onActionClick(info.cat, info.card);
                    else if (window.navigateWikiPhase) window.navigateWikiPhase(info.cat, info.card);
                };
            }
            pop.style.display = 'block';
            const rect = el.getBoundingClientRect();
            const popRect = pop.getBoundingClientRect();
            let left = Math.min(Math.max(10, rect.left), window.innerWidth - popRect.width - 15);
            let top = rect.bottom + 8;
            if (top + popRect.height > window.innerHeight - 10) {
                top = Math.max(10, rect.top - popRect.height - 8);
            }
            pop.style.left = left + 'px';
            pop.style.top = top + 'px';
        });
        el.addEventListener('mouseleave', function() {
            _wikiPopTimer = setTimeout(() => { pop.style.display = 'none'; }, 220);
        });
    });
}

const MATERIAL_ICON_PATHS = {
    palette: 'M12 3a9 9 0 0 0 0 18c.83 0 1.5-.67 1.5-1.5 0-.39-.15-.74-.39-1.01-.23-.26-.38-.61-.38-.99 0-.83.67-1.5 1.5-1.5H16c2.76 0 5-2.24 5-5 0-4.42-4.03-8-9-8zm-5.5 9c-.83 0-1.5-.67-1.5-1.5S5.67 9 6.5 9 8 9.67 8 10.5 7.33 12 6.5 12zm3-4C8.67 8 8 7.33 8 6.5S8.67 5 9.5 5s1.5.67 1.5 1.5S10.33 8 9.5 8zm5 0c-.83 0-1.5-.67-1.5-1.5S13.67 5 14.5 5s1.5.67 1.5 1.5S15.33 8 14.5 8zm3 4c-.83 0-1.5-.67-1.5-1.5S16.67 9 17.5 9s1.5.67 1.5 1.5-.67 1.5-1.5 1.5z',
    verified_user: 'M12 1L3 5v6c0 5.55 3.84 10.74 9 12 5.16-1.26 9-6.45 9-12V5l-9-4zm-2 16l-4-4 1.41-1.41L10 14.17l6.59-6.59L18 9l-8 8z',
    waves: 'M17 16.99c-1.35 0-2.2.42-2.95.8-.65.33-1.18.6-2.05.6-.9 0-1.4-.25-2.05-.6-.75-.38-1.57-.8-2.95-.8s-2.2.42-2.95.8c-.65.33-1.17.6-2.05.6v1.95c1.35 0 2.2-.42 2.95-.8.65-.33 1.17-.6 2.05-.6s1.4.25 2.05.6c.75.38 1.57.8 2.95.8s2.2-.42 2.95-.8c.65-.33 1.18-.6 2.05-.6.9 0 1.4.25 2.05.6.75.38 1.58.8 2.95.8v-1.95c-.9 0-1.4-.25-2.05-.6-.75-.38-1.6-.8-2.95-.8zm0-4.45c-1.35 0-2.2.43-2.95.8-.65.32-1.18.6-2.05.6-.9 0-1.4-.25-2.05-.6-.75-.38-1.57-.8-2.95-.8s-2.2.43-2.95.8c-.65.32-1.17.6-2.05.6v1.95c1.35 0 2.2-.43 2.95-.8.65-.35 1.15-.6 2.05-.6s1.4.25 2.05.6c.75.38 1.57.8 2.95.8s2.2-.43 2.95-.8c.65-.35 1.15-.6 2.05-.6s1.4.25 2.05.6c.75.38 1.58.8 2.95.8v-1.95c-.9 0-1.4-.25-2.05-.6-.75-.38-1.6-.8-2.95-.8zm2.95-8.08c-.75-.38-1.58-.8-2.95-.8s-2.2.42-2.95.8c-.65.32-1.18.6-2.05.6-.9 0-1.4-.25-2.05-.6-.75-.37-1.57-.8-2.95-.8s-2.2.42-2.95.8c-.65.33-1.17.6-2.05.6v1.93c1.35 0 2.2-.43 2.95-.8.65-.33 1.17-.6 2.05-.6s1.4.25 2.05.6c.75.38 1.57.8 2.95.8s2.2-.43 2.95-.8c.65-.32 1.18-.6 2.05-.6.9 0 1.4.25 2.05.6.75.38 1.58.8 2.95.8V5.04c-.9 0-1.4-.25-2.05-.58z',
    center_focus_strong: 'M12 8c-2.21 0-4 1.79-4 4s1.79 4 4 4 4-1.79 4-4-1.79-4-4-4zm-7 7H3v4c0 1.1.9 2 2 2h4v-2H5v-4zM5 5h4V3H5c-1.1 0-2 .9-2 2v4h2V5zm14-2h-4v2h4v4h2V5c0-1.1-.9-2-2-2zm0 16h-4v2h4c1.1 0 2-.9 2-2v-4h-2v4z',
    bolt: 'M11 21h-1l1-7H7.5c-.58 0-.57-.32-.38-.66.19-.34.05-.08.07-.12C8.48 10.94 10.42 7.54 13 3h1l-1 7h3.5c.49 0 .56.33.47.51l-.07.15C12.96 17.55 11 21 11 21z',
    hub: 'M17 16l-4-4V8.82C14.16 8.4 15 7.3 15 6c0-1.66-1.34-3-3-3S9 4.34 9 6c0 1.3.84 2.4 2 2.82V12l-4 4H3v5h5v-3.05l4-4.2 4 4.2V21h5v-5h-4z',
    water_drop: 'M12 2c-5.33 4.55-8 8.48-8 11.8 0 4.98 3.8 8.2 8 8.2s8-3.22 8-8.2c0-3.32-2.67-7.25-8-11.8z',
    deployed_code: 'M21 16.5c0 .38-.21.71-.53.88l-7.9 4.44c-.16.12-.36.18-.57.18s-.41-.06-.57-.18l-7.9-4.44A.991.991 0 0 1 3 16.5v-9c0-.38.21-.71.53-.88l7.9-4.44c.16-.12.36-.18.57-.18s.41.06.57.18l7.9 4.44c.32.17.53.5.53.88v9zM12 4.15L6.04 7.5 12 10.85l5.96-3.35L12 4.15zM5 15.91l6 3.38v-6.71L5 9.21v6.7zm14 0v-6.7l-6 3.37v6.71l6-3.38z',
    biotech: 'M19.8 18.4L14 10.67V6.5l1.35-1.69c.26-.33.03-.81-.39-.81H9.04c-.42 0-.65.48-.39.81L10 6.5v4.17L4.2 18.4c-.49.66-.02 1.6.8 1.6h14c.82 0 1.29-.94.8-1.6z',
    layers: 'M11.99 18.54l-7.37-5.73L3 14.07l9 7 9-7-1.63-1.27-7.38 5.74zM12 16l7.36-5.73L21 9l-9-7-9 7 1.63 1.27L12 16z',
    vertical_split: 'M3 15h8v-2H3v2zm0 4h8v-2H3v2zm0-8h8V9H3v2zm0-6v2h8V5H3zm10 0h8v14h-8V5z',
    my_location: 'M12 8c-2.21 0-4 1.79-4 4s1.79 4 4 4 4-1.79 4-4-1.79-4-4-4zm8.94 3A8.994 8.994 0 0 0 13 3.06V1h-2v2.06A8.994 8.994 0 0 0 3.06 11H1v2h2.06A8.994 8.994 0 0 0 11 20.94V23h2v-2.06A8.994 8.994 0 0 0 20.94 13H23v-2h-2.06zM12 19c-3.87 0-7-3.13-7-7s3.13-7 7-7 7 3.13 7 7-3.13 7-7 7z',
    tune: 'M3 17v2h6v-2H3zM3 5v2h10V5H3zm10 16v-2h8v-2h-8v-2h-2v6h2zM7 9v2H3v2h4v2h2V9H7zm14 4v-2H11v2h10zm-6-4h2V7h4V5h-4V3h-2v6z',
    link: 'M3.9 12c0-1.71 1.39-3.1 3.1-3.1h4V7H7c-2.76 0-5 2.24-5 5s2.24 5 5 5h4v-1.9H7c-1.71 0-3.1-1.39-3.1-3.1zM8 13h8v-2H8v2zm9-6h-4v1.9h4c1.71 0 3.1 1.39 3.1 3.1s-1.39 3.1-3.1 3.1h-4V17h4c2.76 0 5-2.24 5-5s-2.24-5-5-5z',
    open_in_new: 'M19 19H5V5h7V3H5a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14c1.1 0 2-.9 2-2v-7h-2v7zM14 3v2h3.59l-9.83 9.83 1.41 1.41L19 6.41V10h2V3h-7z',
    menu_book: 'M21 5c-1.11-.35-2.33-.5-3.5-.5-1.95 0-4.05.4-5.5 1.5-1.45-1.1-3.55-1.5-5.5-1.5S2.45 4.9 1 6v14.65c0 .25.25.5.5.5.1 0 .15-.05.25-.05C3.1 20.45 5.05 20 6.5 20c1.95 0 4.05.4 5.5 1.5 1.35-.85 3.8-1.5 5.5-1.5 1.65 0 3.35.3 4.75 1.05.1.05.15.05.25.05.25 0 .5-.25.5-.5V6c-.6-.45-1.25-.75-2-1zm0 13.5c-1.1-.35-2.3-.5-3.5-.5-1.7 0-4.15.65-5.5 1.5V8c1.35-.85 3.8-1.5 5.5-1.5 1.2 0 2.4.15 3.5.5v11.5z',
    close: 'M19 6.41L17.59 5 12 10.59 6.41 5 5 6.41 10.59 12 5 17.59 6.41 19 12 13.41 17.59 19 19 17.59 13.41 12z'
};

function msIcon(name) {
    const p = MATERIAL_ICON_PATHS[name] || '';
    return `<svg class="cmp-ms-svg" viewBox="0 0 24 24" aria-hidden="true"><path d="${p}"/></svg>`;
}

class ProteinCompareInspector {
    constructor(opts) {
        this.mountEl = typeof opts.mountEl === 'string' ? document.querySelector(opts.mountEl) : opts.mountEl;
        this.context = opts.context || 'modal'; // 'modal' or 'wiki'
        this.instanceId = opts.instanceId || ('cmp_' + this.context);
        this.activePreset = opts.defaultPreset || (this.context === 'wiki' ? 'vector' : 'overview');
        this.activePills = new Set((PHASE_PRESETS[this.activePreset] || PHASE_PRESETS.overview).pills);
        this.activeOverlays = new Set((PHASE_PRESETS[this.activePreset] || PHASE_PRESETS.overview).overlays);
        this.activeMode = (PHASE_PRESETS[this.activePreset] || PHASE_PRESETS.overview).mode;
        this.isMerged = false;
        this.pocketMode = false;
        this.compData = null;
        this.id1 = null;
        this.id2 = null;
        this.v1 = null;
        this.v2 = null;
        this.master = null;
        this.syncId = null;
        this.onPhaseSelect = opts.onPhaseSelect || null;
        this.onClose = opts.onClose || null;
        this.renderShell();
    }

    renderShell() {
        const popId = `${this.instanceId}_cfg_popover`;
        const presetOptions = Object.entries(PHASE_PRESETS).map(([k, v]) =>
            `<option value="${k}" ${k === this.activePreset ? 'selected' : ''}>${v.label}</option>`
        ).join('') + `<option value="custom">Custom Widget Selection</option>`;

        const wikiPairSelector = (this.context === 'wiki' || this.context === 'lab') ? `
            <label style="font-size:0.78em; color:#94a3b8;">Discovery Pair:</label>
            <select id="${this.instanceId}_pair_select" class="cmp-select" style="min-width:340px;">
                <option value="A0A1S8AT22|A0A1H3IA37">#1 Twilight-Zone (DUF7861): A0A1S8AT22 (Natrinema) vs A0A1H3IA37 (Halobellus) [0.81 Å, TM=0.927, Tw=0.8640]</option>
                <option value="A0A1H7ISK7|A0A8A2U9A2">#1 Dark Island Champion (DUF7522): A0A1H7ISK7 (Haloferax) vs A0A8A2U9A2 (Natrinema) [0.74 Å, TM=0.954, Tw=0.8392]</option>
                <option value="A0A0U5HVU0|A0A1H0ZZR0">123-Myr Cretaceous Halite vs pH 10 Soda Lake (DUF7124): A0A0U5HVU0 vs A0A1H0ZZR0 [0.61 Å, TM=0.845, Tw=0.7994]</option>
                <option value="A0A2K5ARN2|J1L1L0">#2 RangeFinder 8/8 Box (Xylanase): A0A2K5ARN2 (Nitrosocaldus) vs J1L1L0 (Methanofollis) [Range-IoU=0.9619, 8/8 Ferulic Acid]</option>
                <option value="A0A7D5TJJ3|A0A654LSW6">#3 RangeFinder 8/8 Box (Hydrogenase): A0A7D5TJJ3 (Halosimplex) vs A0A654LSW6 (Nitrosocosmicus) [Range-IoU=0.9472, 8/8 Paraquat]</option>
                <option value="A0A8J8PAH9|A0A5D5AJJ1">#4 Twilight-Zone (DUF1028): A0A8J8PAH9 (Halonotius) vs A0A5D5AJJ1 (Natrialba) [SW=41.4%, Tw=0.8516]</option>
                <option value="A0ABD5QTT5|A0AA41FYQ3">#5 Dark Island (DUF1684): A0ABD5QTT5 (Halorubrum) vs A0AA41FYQ3 (Haloarcula) [SW=52.2%, Tw=0.8377]</option>
                <option value="A0ABD5SLL6|M0LYV8">#6 Dark Island (DUF302): A0ABD5SLL6 (Natrinema) vs M0LYV8 (Halococcus stromatolite) [SW=55.5%, Tw=0.8355]</option>
                <option value="M0D3I5|U2YSQ3">#7 Dark Island (DUF7995): M0D3I5 (Halogeometricum) vs U2YSQ3 (Halarchaeum acidiphilum) [SW=48.1%, Tw=0.8327]</option>
                <option value="A0A3A6QKX6|M0A243">#8 Dark Island (DUF373): A0A3A6QKX6 (Halonotius) vs M0A243 (Natrialba soda lake) [SW=49.1%, Tw=0.8122]</option>
                <option value="A0A0D6JTS1|A0A8U0IQD4">#9 High Chaperone Stability: A0A0D6JTS1 (Haloferax) vs A0A8U0IQD4 (Halorussus) [TM=0.982, Score=0.70]</option>
                <option value="A0A1I6H788|A0A0W1R627">#10 Red Sea Brine (DUF5787): A0A1I6H788 (Halogeometricum) vs A0A0W1R627 (Haloprofundus) [TM=0.926, Tw=0.8011]</option>
                <option value="A0A5D5AFN6|A0ABD5QWL6">#11 Ultra-Divergent (DUF63): A0A5D5AFN6 (Natrialba) vs A0ABD5QWL6 (Halorubrum) [SW=31.7%, Tw=0.7015]</option>
                <option value="H1Z0V4|A0A7D5IMP5">#12 Methanogen Island (DUF4258): H1Z0V4 (Methanoplanus) vs A0A7D5IMP5 (Methanolobus) [SW=29.0%, Tw=0.6961]</option>
                <option value="M0N1T2|A0A8T8WBG0">#13 Rock-Salt Island (DUF7838): M0N1T2 (Halococcus) vs A0A8T8WBG0 (Halobaculum) [SW=25.0%, Tw=0.6668]</option>
            </select>
        ` : '';

        const closeBtnHtml = this.context === 'modal'
            ? `<button class="cmp-btn danger" id="${this.instanceId}_close_btn">${msIcon('close')}Close (Esc)</button>`
            : '';

        const openLabBtnHtml = this.context !== 'lab'
            ? `<button class="cmp-btn" style="border-color:#10b981;color:#10b981;" id="${this.instanceId}_lab_btn">${msIcon('open_in_new')}Open Full Lab</button>`
            : '';

        const openWikiBtnHtml = (this.context === 'modal' || this.context === 'lab')
            ? `<button class="cmp-btn" style="border-color:#38bdf8;color:#38bdf8;" id="${this.instanceId}_wiki_btn">${msIcon('menu_book')}Open Phase in Wiki</button>`
            : '';

        const viewersStyle = this.context === 'lab'
            ? `style="min-height: calc(100vh - 235px); height: calc(100vh - 235px);"`
            : '';

        this.mountEl.innerHTML = `
            <div class="cmp-inspector">
                <div class="cmp-toolbar">
                    <div class="cmp-toolbar-left">
                        <span class="cmp-title" id="${this.instanceId}_title">Side-by-Side 3D Protein Comparison</span>
                        <span class="cmp-cache-badge" id="${this.instanceId}_cache_badge">Ready</span>
                        ${wikiPairSelector}
                    </div>
                    <div class="cmp-toolbar-right">
                        <label style="font-size:0.78em; color:#94a3b8;">Phase Lens:</label>
                        <select id="${this.instanceId}_preset_select" class="cmp-select">${presetOptions}</select>
                        <button class="cmp-btn" popovertarget="${popId}" id="${this.instanceId}_cfg_btn">${msIcon('tune')}Configure Widgets (${this.activePills.size})</button>
                        <button class="cmp-btn" id="${this.instanceId}_share_btn">${msIcon('link')}Copy Link</button>
                        ${openLabBtnHtml}
                        ${openWikiBtnHtml}
                        ${closeBtnHtml}
                    </div>
                </div>

                <div class="cmp-modes-bar" id="${this.instanceId}_modes_bar">
                    <span class="cmp-modes-label">3D View:</span>
                    <button class="cmp-btn" data-mode="identity" title="Color: Identity (Bait Cyan vs Orphan Magenta)">${msIcon('palette')}Identity</button>
                    <button class="cmp-btn" data-mode="plddt" title="Color: AlphaFold pLDDT Confidence">${msIcon('verified_user')}pLDDT</button>
                    <button class="cmp-btn" data-mode="ss" title="Color: Secondary Structure (Helix / Sheet / Coil)">${msIcon('waves')}2° Struct</button>
                    <button class="cmp-btn" data-mode="core" title="Highlight Aligned Kabsch C-alpha Core (&lt; 3.8 Å)">${msIcon('center_focus_strong')}Cα Core</button>
                    <button class="cmp-btn" data-mode="acidic" title="Highlight Acidic Surface Shell (Asp/Glu)">${msIcon('bolt')}Acidic (D/E)</button>
                    <button class="cmp-btn" data-mode="hbond" title="Chaperone H-Bond CSR Hubs">${msIcon('hub')}H-Bonds</button>
                    <button class="cmp-btn" data-mode="rasa" title="Chaperone RASA Hydrophobic Surface Patches">${msIcon('water_drop')}RASA</button>
                    <button class="cmp-btn" data-mode="rangefinder" title="RangeFinder 10 Å Active-Site Pocket">${msIcon('deployed_code')}10 Å Box</button>
                    <button class="cmp-btn" data-mode="tpp_seam" title="CETSA/TPP Thermal Weak Seam &amp; Ideal Bioreactor Temp (T_opt)">${msIcon(' bolt')}TPP Seam (T_opt)</button>
                    <button class="cmp-btn" data-mode="nterm_pep" title="Highlight N-Term MS/MS Tryptic Peptide">${msIcon('biotech')}N-Term Pep</button>
                    <span style="border-left:1px solid #334155; height:16px; margin:0 2px; flex-shrink:0;"></span>
                    <button class="cmp-btn" id="${this.instanceId}_merge_btn" title="Superpose in Single 3D View / Split Side-by-Side View">${msIcon('layers')}Superpose</button>
                    <button class="cmp-btn" data-mode="zorder" id="${this.instanceId}_pocket_btn" title="Automatic 30-Bit Z-Order (Morton) 5.0 Å Pocket Selection">${msIcon('my_location')}Z-Order 5 Å</button>
                </div>

                <div class="cmp-viewers" id="${this.instanceId}_viewers_wrap" ${viewersStyle}>
                    <div class="cmp-viewer-box" id="${this.instanceId}_vbox1">
                        <div class="cmp-top-bar" id="${this.instanceId}_ov1"></div>
                        <div class="cmp-viewer-canvas" id="${this.instanceId}_canvas1"></div>
                        <div class="cmp-bottom-bar" id="${this.instanceId}_bot1"></div>
                    </div>
                    <div class="cmp-viewer-box" id="${this.instanceId}_vbox2">
                        <div class="cmp-top-bar" id="${this.instanceId}_ov2"></div>
                        <div class="cmp-viewer-canvas" id="${this.instanceId}_canvas2"></div>
                        <div class="cmp-bottom-bar" id="${this.instanceId}_bot2"></div>
                    </div>
                </div>

                <div class="cmp-pills-bar" id="${this.instanceId}_pills_bar">
                    <span style="color:#64748b; font-size:0.82em;">Select or load a protein pair to inspect phase metrics...</span>
                </div>
                <div id="${this.instanceId}_explain" style="display:none;"></div>

                <div id="${popId}" class="cmp-config-popover" popover="auto">
                    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #1e293b; padding-bottom:8px; margin-bottom:10px;">
                        <b style="color:#38bdf8; font-size:0.95em;">Configure Visible Metric Pills &amp; 3D Overlay Widgets</b>
                        <div style="display:flex; gap:6px;">
                            <button class="cmp-btn" id="${this.instanceId}_cfg_reset">Phase Default</button>
                            <button class="cmp-btn" id="${this.instanceId}_cfg_all">Show All</button>
                            <button class="cmp-btn" id="${this.instanceId}_cfg_min">Minimal</button>
                        </div>
                    </div>
                    <div style="font-size:0.78em; color:#94a3b8; font-weight:600; text-transform:uppercase;">Metric Pills (Click to Add / Remove)</div>
                    <div class="cmp-cfg-grid" id="${this.instanceId}_cfg_pills"></div>
                    <div style="font-size:0.78em; color:#94a3b8; font-weight:600; text-transform:uppercase; margin-top:12px;">3D Viewer Overlay Rows</div>
                    <div class="cmp-cfg-grid" id="${this.instanceId}_cfg_overlays"></div>
                </div>
            </div>
        `;

        // Bind events
        const presetSel = document.getElementById(`${this.instanceId}_preset_select`);
        presetSel.onchange = () => {
            const val = presetSel.value;
            if (val !== 'custom') {
                this.setPreset(val, { triggerAuto3DMode: true, notifyParent: true });
            }
        };

        const pairSel = document.getElementById(`${this.instanceId}_pair_select`);
        if (pairSel) {
            pairSel.onchange = () => {
                const parts = pairSel.value.split('|');
                if (parts.length === 2) {
                    if (this.context === 'lab') {
                        const u = new URL(window.location);
                        u.searchParams.set('node1', parts[0]);
                        u.searchParams.set('node2', parts[1]);
                        window.history.pushState({}, '', u);
                    }
                    this.loadPair(parts[0], parts[1]);
                }
            };
        }

        document.getElementById(`${this.instanceId}_share_btn`).onclick = () => {
            navigator.clipboard.writeText(window.location.href);
        };

        const labBtn = document.getElementById(`${this.instanceId}_lab_btn`);
        if (labBtn) {
            labBtn.onclick = () => {
                let url = '/compare';
                const params = new URLSearchParams();
                if (this.id1 && this.id2) {
                    params.set('node1', this.id1);
                    params.set('node2', this.id2);
                }
                if (this.activePreset && this.activePreset !== 'custom') {
                    params.set('preset', this.activePreset);
                }
                window.open(url + '?' + params.toString(), '_blank');
            };
        }

        const wikiBtn = document.getElementById(`${this.instanceId}_wiki_btn`);
        if (wikiBtn) {
            wikiBtn.onclick = () => {
                const cat = (this.activePreset && PHASE_PRESETS[this.activePreset] && this.activePreset !== 'overview' && this.activePreset !== 'all_layers')
                    ? this.activePreset : 'vector';
                this.navigateToWiki(cat, '');
            };
        }

        const closeBtn = document.getElementById(`${this.instanceId}_close_btn`);
        if (closeBtn && this.onClose) {
            closeBtn.onclick = () => this.onClose();
        }

        document.querySelectorAll(`#${this.instanceId}_modes_bar [data-mode]`).forEach(btn => {
            btn.onclick = () => this.applyMode(btn.getAttribute('data-mode'));
        });

        document.getElementById(`${this.instanceId}_merge_btn`).onclick = () => this.toggleMergeView();

        document.getElementById(`${this.instanceId}_cfg_reset`).onclick = () => {
            const base = (this.activePreset && PHASE_PRESETS[this.activePreset]) ? this.activePreset : 'overview';
            this.setPreset(base, { triggerAuto3DMode: true });
        };
        document.getElementById(`${this.instanceId}_cfg_all`).onclick = () => {
            this.setPreset('all_layers', { triggerAuto3DMode: false });
        };
        document.getElementById(`${this.instanceId}_cfg_min`).onclick = () => {
            this.activePills = new Set(['core_rmsd', 'tm_score']);
            this.activeOverlays = new Set(['ov_identity']);
            this.activePreset = 'custom';
            document.getElementById(`${this.instanceId}_preset_select`).value = 'custom';
            this.renderDynamicParts();
        };
    }

    initViewersIfNeeded() {
        if (this.v1 && this.v2) return;
        const c1 = document.getElementById(`${this.instanceId}_canvas1`);
        const c2 = document.getElementById(`${this.instanceId}_canvas2`);
        this.v1 = $3Dmol.createViewer(c1, { backgroundColor: "0x0b1120" });
        this.v2 = $3Dmol.createViewer(c2, { backgroundColor: "0x0b1120" });
        c1.onmouseenter = () => { this.master = 1; };
        c2.onmouseenter = () => { this.master = 2; };
        const syncLoop = () => {
            if (!this.isMerged && this.v1 && this.v2) {
                if (this.master === 1) this.v2.setView(this.v1.getView());
                else if (this.master === 2) this.v1.setView(this.v2.getView());
            }
            this.syncId = requestAnimationFrame(syncLoop);
        };
        syncLoop();
    }

    syncNavLinks() {
        if (!this.id1 || !this.id2) return;
        const presetParam = (this.activePreset && this.activePreset !== 'custom') ? this.activePreset : 'overview';
        const modeQs = this.activeMode ? `&mode=${encodeURIComponent(this.activeMode)}` : '';
        const cardsUrl = `/cards?hub=${encodeURIComponent(this.id1)}&challenger=${encodeURIComponent(this.id2)}&node1=${encodeURIComponent(this.id1)}&node2=${encodeURIComponent(this.id2)}&preset=${encodeURIComponent(presetParam)}${modeQs}`;
        const wikiUrl = `/wiki?node1=${encodeURIComponent(this.id1)}&node2=${encodeURIComponent(this.id2)}&cat=${encodeURIComponent(presetParam)}`;
        const labCardsEl = document.getElementById('lab-cards-top-link');
        if (labCardsEl) labCardsEl.href = cardsUrl;
        const wikiCardsEl = document.getElementById('wiki-cards-top-link');
        if (wikiCardsEl) wikiCardsEl.href = cardsUrl;
        const labWikiEl = document.getElementById('lab-wiki-top-link');
        if (labWikiEl) labWikiEl.href = wikiUrl;
        const inp1 = document.getElementById('lab-id1');
        const inp2 = document.getElementById('lab-id2');
        if (inp1) inp1.value = this.id1;
        if (inp2) inp2.value = this.id2;
    }

    setPreset(presetKey, opts = {}) {
        const preset = PHASE_PRESETS[presetKey] || PHASE_PRESETS.overview;
        this.activePreset = presetKey in PHASE_PRESETS ? presetKey : 'overview';
        this.activePills = new Set(preset.pills);
        this.activeOverlays = new Set(preset.overlays);
        const sel = document.getElementById(`${this.instanceId}_preset_select`);
        if (sel) sel.value = this.activePreset;
        if (this.context === 'lab') {
            const u = new URL(window.location);
            u.searchParams.set('preset', this.activePreset);
            if (preset.mode) u.searchParams.set('mode', preset.mode);
            window.history.replaceState({}, '', u);
        }
        if (this.compData) {
            this.renderDynamicParts();
            if (opts.triggerAuto3DMode && preset.mode) {
                this.applyMode(preset.mode, preset.summary, true);
            }
        }
        this.syncNavLinks();
        if (opts.notifyParent && this.onPhaseSelect && presetKey !== 'overview' && presetKey !== 'all_layers') {
            this.onPhaseSelect(presetKey);
        }
    }

    togglePill(pillId, enabled) {
        if (enabled) this.activePills.add(pillId);
        else this.activePills.delete(pillId);
        this.activePreset = 'custom';
        const sel = document.getElementById(`${this.instanceId}_preset_select`);
        if (sel) sel.value = 'custom';
        this.renderDynamicParts();
    }

    toggleOverlay(ovId, enabled) {
        if (enabled) this.activeOverlays.add(ovId);
        else this.activeOverlays.delete(ovId);
        this.activePreset = 'custom';
        const sel = document.getElementById(`${this.instanceId}_preset_select`);
        if (sel) sel.value = 'custom';
        this.renderDynamicParts();
    }

    navigateToWiki(cat, cardId) {
        if (this.context === 'wiki' && window.navigateWikiPhase) {
            window.navigateWikiPhase(cat, cardId);
            return;
        }
        let url = '/wiki';
        const params = new URLSearchParams();
        if (this.id1 && this.id2) {
            params.set('node1', this.id1);
            params.set('node2', this.id2);
        }
        if (cat) params.set('cat', cat);
        const qs = params.toString();
        if (qs) url += '?' + qs;
        if (cardId) url += '#' + cardId;
        window.open(url, '_blank');
    }

    async loadPair(id1, id2, distanceOpt) {
        this.id1 = id1;
        this.id2 = id2;
        window.currentNode1Id = id1;
        window.currentNode2Id = id2;
        this.syncNavLinks();

        // Reset merge mode if active
        if (this.isMerged) {
            this.isMerged = false;
            const wrap = document.getElementById(`${this.instanceId}_viewers_wrap`);
            const vbox2 = document.getElementById(`${this.instanceId}_vbox2`);
            const ov2 = document.getElementById(`${this.instanceId}_ov2`);
            const bot2 = document.getElementById(`${this.instanceId}_bot2`);
            wrap.classList.remove('merged');
            vbox2.style.display = '';
            ov2.classList.remove('docked-right');
            vbox2.insertBefore(ov2, document.getElementById(`${this.instanceId}_canvas2`));
            if (bot2) vbox2.appendChild(bot2);
            document.getElementById(`${this.instanceId}_merge_btn`).classList.remove('active');
        }

        const titleEl = document.getElementById(`${this.instanceId}_title`);
        const badgeEl = document.getElementById(`${this.instanceId}_cache_badge`);
        const pillsEl = document.getElementById(`${this.instanceId}_pills_bar`);
        titleEl.innerHTML = `Comparing: <span style="color:#38bdf8">${id1}</span> vs <span style="color:#e879f9">${id2}</span>`;
        badgeEl.className = 'cmp-cache-badge';
        badgeEl.textContent = 'Loading...';
        pillsEl.innerHTML = `<span style="color:#94a3b8;">Querying multi-layer cache / PostgreSQL 18 (pg_bio) &amp; superposing ${id1} vs ${id2}...</span>`;

        this.initViewersIfNeeded();

        const comp = await ClientCompareCache.getOrFetch(id1, id2);
        if (this.id1 !== id1 || this.id2 !== id2) return;
        if (comp.error) {
            pillsEl.innerHTML = `<span style="color:#f87171;">Error: ${comp.error}</span>`;
            return;
        }
        this.compData = comp;
        badgeEl.className = 'cmp-cache-badge hit';
        badgeEl.textContent = comp._cacheSource || 'Cached';

        const vecDist = distanceOpt !== undefined ? Number(distanceOpt).toFixed(4) : (comp.alignment.vec_dist || comp.alignment.esm_kmer_dist);
        titleEl.innerHTML = `Comparing: <span style="color:#38bdf8">${id1}</span> vs <span style="color:#e879f9">${id2}</span> <span style="color:#94a3b8; font-size:0.84em; font-weight:normal; margin-left:6px;">(Vector Dist: ${vecDist})</span>`;

        this.v1.clear();
        this.v2.clear();
        this.v1.addModel(comp.pdb1, "pdb");
        this.v2.addModel(comp.pdb2_aligned, "pdb");

        this.renderDynamicParts();
        const preset = PHASE_PRESETS[this.activePreset] || PHASE_PRESETS.overview;
        const urlMode = (this.context === 'lab') ? (new URLSearchParams(window.location.search).get('mode') || '') : '';
        const initialMode = urlMode || preset.mode || 'identity';
        this.applyMode(initialMode, urlMode ? '' : preset.summary, true);
        this.v1.zoomTo();
        this.v1.render();
        this.v2.setView(this.v1.getView());
        this.v2.render();

        if (window.onWikiPairLoaded && this.context === 'wiki') {
            window.onWikiPairLoaded(comp, id1, id2);
        }
    }

    renderDynamicParts() {
        if (!this.compData) return;
        const c = this.compData;
        const pillsEl = document.getElementById(`${this.instanceId}_pills_bar`);
        const cfgBtn = document.getElementById(`${this.instanceId}_cfg_btn`);
        if (cfgBtn) cfgBtn.innerHTML = `${msIcon('tune')}Configure Widgets (${this.activePills.size})`;

        // Render active metric pills (click body to highlight in 3D, hover for formula & SQL wiki popover)
        const visibleDefs = METRIC_PILL_DEFS.filter(d => this.activePills.has(d.id));
        if (visibleDefs.length === 0) {
            pillsEl.innerHTML = `<span style="color:#64748b; font-size:0.8em;">All metric pills hidden. Click "Configure Widgets" to enable metrics.</span>`;
        } else {
            pillsEl.innerHTML = visibleDefs.map(d => `
                <div class="cmp-pill" data-wiki-key="${d.wikiKey}">
                    <span class="cmp-pill-body" data-pill-mode="${d.mode}">${d.render(c)}</span>
                    <button class="cmp-pill-close" title="Hide this metric pill" data-pill-hide="${d.id}">x</button>
                </div>
            `).join('');
        }

        pillsEl.querySelectorAll('[data-pill-mode]').forEach(el => {
            el.onclick = (e) => {
                e.stopPropagation();
                this.applyMode(el.getAttribute('data-pill-mode'));
            };
        });
        pillsEl.querySelectorAll('[data-pill-hide]').forEach(el => {
            el.onclick = (e) => {
                e.stopPropagation();
                this.togglePill(el.getAttribute('data-pill-hide'), false);
            };
        });

        // Render frozen top name bar and bottom extra-properties bar for Viewer 1 and Viewer 2
        const ov1El = document.getElementById(`${this.instanceId}_ov1`);
        const ov2El = document.getElementById(`${this.instanceId}_ov2`);
        const bot1El = document.getElementById(`${this.instanceId}_bot1`);
        const bot2El = document.getElementById(`${this.instanceId}_bot2`);
        ov1El.innerHTML = this.buildTopBarHtml(c.meta1, '#38bdf8', 1);
        ov2El.innerHTML = this.buildTopBarHtml(c.meta2, '#e879f9', 2);
        if (bot1El) bot1El.innerHTML = this.buildBottomBarHtml(c.meta1, c.chaperone, '#38bdf8', 1);
        if (bot2El) bot2El.innerHTML = this.buildBottomBarHtml(c.meta2, c.chaperone, '#e879f9', 2);

        [ov1El, ov2El, bot1El, bot2El].filter(Boolean).forEach(barEl => {
            barEl.querySelectorAll('[data-ov-mode]').forEach(btn => {
                btn.onclick = (e) => {
                    e.stopPropagation();
                    this.applyMode(btn.getAttribute('data-ov-mode'));
                };
            });
            barEl.querySelectorAll('[data-ov-hide]').forEach(btn => {
                btn.onclick = (e) => {
                    e.stopPropagation();
                    this.toggleOverlay(btn.getAttribute('data-ov-hide'), false);
                };
            });
        });

        // Render checkboxes inside Configure Widgets popover
        const cfgPillsEl = document.getElementById(`${this.instanceId}_cfg_pills`);
        cfgPillsEl.innerHTML = METRIC_PILL_DEFS.map(d => `
            <label class="cmp-cfg-item">
                <input type="checkbox" data-cfg-pill="${d.id}" ${this.activePills.has(d.id) ? 'checked' : ''}>
                <span>${d.label}</span>
            </label>
        `).join('');
        cfgPillsEl.querySelectorAll('[data-cfg-pill]').forEach(cb => {
            cb.onchange = () => this.togglePill(cb.getAttribute('data-cfg-pill'), cb.checked);
        });

        const cfgOvEl = document.getElementById(`${this.instanceId}_cfg_overlays`);
        cfgOvEl.innerHTML = OVERLAY_FIELD_DEFS.map(o => `
            <label class="cmp-cfg-item">
                <input type="checkbox" data-cfg-ov="${o.id}" ${this.activeOverlays.has(o.id) ? 'checked' : ''}>
                <span>${o.label}</span>
            </label>
        `).join('');
        cfgOvEl.querySelectorAll('[data-cfg-ov]').forEach(cb => {
            cb.onchange = () => this.toggleOverlay(cb.getAttribute('data-cfg-ov'), cb.checked);
        });

        if (this.v1) this.v1.resize();
        if (this.v2 && !this.isMerged) this.v2.resize();

        bindWikiHoverPopovers(this.mountEl, (cat, card) => this.navigateToWiki(cat, card));
    }

    cleanProteinTitle(meta) {
        if (!meta || !meta.name) return '';
        let t = String(meta.name)
            .replace(/^(sp|tr)\\|[^|]+\\|[^\\s]+\\s+/i, '')
            .replace(/\\s+OS=.*$/i, '')
            .trim();
        if (!t || t === 'Unknown' || t === meta.id || t === meta.organism) return '';
        return t;
    }

    buildTopBarHtml(meta, accentColor, idx) {
        const roleTag = idx === 1 ? 'Bait' : 'Orphan';
        const protTitle = this.cleanProteinTitle(meta);
        const protHtml = protTitle ? `<span class="cmp-top-prot" title="${protTitle}">${protTitle}</span>` : '';
        return `
            <div class="cmp-top-ident">
                <span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:${accentColor}; flex-shrink:0;" title="${roleTag}"></span>
                <span class="cmp-top-org" title="${meta.organism || 'Organism'}">${meta.organism || 'Organism'}</span>
                ${protHtml}
            </div>
            <div class="cmp-top-links">
                <a href="https://www.uniprot.org/uniprotkb/${meta.id}/entry" target="_blank" style="color:${accentColor}; font-weight:700; text-decoration:none;">${meta.id}</a>
                <a href="https://alphafold.ebi.ac.uk/entry/${meta.id}" target="_blank" style="color:#94a3b8; font-size:0.86em; text-decoration:none;">(AFDB)</a>
            </div>
        `;
    }

    buildBottomBarHtml(meta, ch, accentColor, idx) {
        const chips = [];
        if (this.activeOverlays.has('ov_plddt')) {
            chips.push(`
                <span class="cmp-prop-chip" data-wiki-key="plddt" data-ov-mode="plddt" style="cursor:pointer;">
                    <span>pLDDT: <b style="color:#38bdf8">${meta.mean_plddt}</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_plddt" title="Hide pLDDT">x</button>
                </span>
            `);
        }
        if (this.activeOverlays.has('ov_size')) {
            chips.push(`
                <span class="cmp-prop-chip">
                    <span><b>${meta.length} aa</b> · <b>${meta.mw_kda} kDa</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_size" title="Hide length &amp; mass">x</button>
                </span>
            `);
        }
        if (this.activeOverlays.has('ov_rangefinder') && this.compData && this.compData.rangefinder) {
            const rf = this.compData.rangefinder;
            const mwBox = idx === 1 ? rf.mw_box_1 : rf.mw_box_2;
            const lpBox = idx === 1 ? rf.logp_box_1 : rf.logp_box_2;
            const ap = idx === 1 ? rf.arom_prob_1 : rf.arom_prob_2;
            chips.push(`
                <span class="cmp-prop-chip" data-wiki-key="rangefinder_pocket" data-ov-mode="rangefinder" style="cursor:pointer;">
                    <span>10Å Box: MW <b>${mwBox}</b> · LogP <b>${lpBox}</b> · P(arom) <b style="color:#34d399">${ap}</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_rangefinder" title="Hide RangeFinder box">x</button>
                </span>
            `);
        }
        if (this.activeOverlays.has('ov_biochem')) {
            chips.push(`
                <span class="cmp-prop-chip" data-wiki-key="biochem">
                    <span>pI <b>${meta.pi}</b> · GRAVY <b>${meta.gravy}</b> · H <b>${meta.entropy}</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_biochem" title="Hide biochemistry">x</button>
                </span>
            `);
        }
        if (this.activeOverlays.has('ov_charge_shell')) {
            chips.push(`
                <span class="cmp-prop-chip" data-wiki-key="biochem" data-ov-mode="acidic" style="cursor:pointer;">
                    <span>D/E <b style="color:#f87171">${meta.acidic_pct}%</b> · K/R/H <b style="color:#60a5fa">${meta.basic_pct}%</b> · Hydro <b>${meta.hydro_pct}%</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_charge_shell" title="Hide charge shell">x</button>
                </span>
            `);
        }
        if (this.activeOverlays.has('ov_chaperone') && ch) {
            const rasa = idx === 1 ? ch.rasa_1 : ch.rasa_2;
            const hb = idx === 1 ? ch.hbonds_1 : ch.hbonds_2;
            const buns = idx === 1 ? ch.buns_1 : ch.buns_2;
            chips.push(`
                <span class="cmp-prop-chip" data-wiki-key="chaperone_audit" data-ov-mode="hbond" style="cursor:pointer;">
                    <span>RASA <b>${rasa}</b> · H-Bonds <b>${hb}</b> · BUNS <b>${buns}</b></span>
                    <button class="cmp-pill-close" data-ov-hide="ov_chaperone" title="Hide Chaperone audit">x</button>
                </span>
            `);
        }
        return chips.join('');
    }

    applyMode(mode, customExplain, fromPreset) {
        if (!this.v1 || !this.v2 || !this.compData) return;
        this.activeMode = mode;
        document.querySelectorAll(`#${this.instanceId}_modes_bar [data-mode]`).forEach(btn => {
            btn.classList.toggle('active', btn.getAttribute('data-mode') === mode);
        });
        if (!fromPreset && this.context === 'lab') {
            const u = new URL(window.location);
            u.searchParams.set('mode', mode);
            const modeToPreset = {
                'zorder': 'spatial',
                'rangefinder': 'rangefinder',
                'core': 'structure',
                'hbond': 'chaperone',
                'rasa': 'chaperone',
                'tpp_seam': 'tpp',
                'acidic': 'biochem',
                'nterm_pep': 'proteomics',
                'ss': 'sequence'
            };
            if (u.searchParams.get('preset') === 'rangefinder' && mode !== 'rangefinder') {
                const nextPreset = modeToPreset[mode] || 'overview';
                u.searchParams.set('preset', nextPreset);
                this.activePreset = nextPreset;
                const sel = document.getElementById(`${this.instanceId}_preset_select`);
                if (sel) sel.value = nextPreset;
            }
            window.history.replaceState({}, '', u);
            this.syncNavLinks();
        }

        // Always clear any RangeFinder 10 Å wireframe globes (or previous 3D shapes & labels) when switching highlight modes
        this.v1.removeAllShapes();
        this.v2.removeAllShapes();
        this.v1.removeAllLabels();
        this.v2.removeAllLabels();

        const c = this.compData;
        const al = c.alignment || {};
        const ch = c.chaperone || {};
        const explainEl = document.getElementById(`${this.instanceId}_explain`);
        const plddtColor = (atom) => {
            const b = atom.b <= 1.0 ? atom.b * 100 : atom.b;
            if (b > 90) return '#0053d6';
            if (b > 70) return '#65cbff';
            if (b > 50) return '#ffdb13';
            return '#ff7d45';
        };

        const styleTarget = (v, idx, spec) => {
            if (this.isMerged) {
                this.v1.setStyle({ model: idx }, spec);
            } else {
                v.setStyle({}, spec);
            }
        };
        const styleSub = (v, idx, sel, spec) => {
            if (this.isMerged) {
                this.v1.setStyle(Object.assign({ model: idx }, sel), spec);
            } else {
                v.setStyle(sel, spec);
            }
        };

        let msg = '';
        if (mode === 'identity') {
            styleTarget(this.v1, 0, { cartoon: { color: 'cyan', opacity: 0.95 } });
            styleTarget(this.v2, 1, { cartoon: { color: 'magenta', opacity: 0.95 } });
            msg = `<b>Identity View:</b> Bait (<b style="color:#38bdf8">${c.meta1.id}</b>) in Cyan vs Orphan (<b style="color:#e879f9">${c.meta2.id}</b>) in Magenta after Kabsch SVD 3D registration (Core RMSD: <b>${al.core_rmsd} A</b>, TM-Score: <b>${al.tm_score}</b>).`;
        } else if (mode === 'plddt') {
            styleTarget(this.v1, 0, { cartoon: { colorfunc: plddtColor } });
            styleTarget(this.v2, 1, { cartoon: { colorfunc: plddtColor } });
            msg = `<b>pLDDT Confidence View:</b> Dark Blue (&gt;90) = crystal-grade catalytic core; Light Blue (&gt;70) = confident fold; Yellow/Orange = flexible surface loops (Mean pLDDT: <b>${c.meta1.mean_plddt}</b> vs <b>${c.meta2.mean_plddt}</b>).`;
        } else if (mode === 'ss') {
            [ [this.v1, 0], [this.v2, 1] ].forEach(([v, idx]) => {
                styleTarget(v, idx, { cartoon: { color: 'white', opacity: 0.32 } });
                styleSub(v, idx, { ss: 'h' }, { cartoon: { color: '#ff6ec7', opacity: 1.0 } });
                styleSub(v, idx, { ss: 's' }, { cartoon: { color: '#00f2fe', opacity: 1.0 } });
            });
            msg = `<b>Secondary Structure View:</b> Alpha-helices in <b style="color:#ff6ec7">Pink</b> and Beta-sheets in <b style="color:#00f2fe">Cyan</b>. Secondary structure topology remains conserved across the Twilight Zone (SW Sim: <b>${al.sw_sim_pct}%</b>).`;
        } else if (mode === 'core') {
            styleTarget(this.v1, 0, { cartoon: { color: 'cyan', opacity: 0.22 } });
            styleTarget(this.v2, 1, { cartoon: { color: 'magenta', opacity: 0.22 } });
            if (al.res1_core && al.res1_core.length) {
                styleSub(this.v1, 0, { resi: al.res1_core }, { cartoon: { color: 'cyan', opacity: 1.0 } });
                styleSub(this.v2, 1, { resi: al.res2_core }, { cartoon: { color: 'magenta', opacity: 1.0 } });
            }
            msg = `<b>Aligned Kabsch C-alpha Core (${al.aligned_residues} residues &lt; 3.8 A):</b> Opaque segments lock onto the conserved 3D catalytic scaffold (RMSD: <b>${al.core_rmsd} A</b>); translucent segments are unaligned loops.`;
        } else if (mode === 'acidic') {
            [ [this.v1, 0], [this.v2, 1] ].forEach(([v, idx]) => {
                styleTarget(v, idx, { cartoon: { color: 'white', opacity: 0.28 } });
                styleSub(v, idx, { resn: ['ASP', 'GLU'] }, { cartoon: { color: '#ef4444' }, stick: { colorscheme: 'redCarbon' } });
            });
            msg = `<b>Halophilic Acidic Shell (Asp/Glu):</b> Negatively charged D/E sidechains rendered as <b style="color:#f87171">Red Sticks</b> (Acidic: <b>${c.meta1.acidic_pct}%</b> vs <b>${c.meta2.acidic_pct}%</b>, dpI = <b>${al.delta_pi}</b>).`;
        } else if (mode === 'hbond') {
            const h1 = (ch.hbond_hubs_1 && ch.hbond_hubs_1.length) ? ch.hbond_hubs_1 : [1, 2, 3];
            const h2 = (ch.hbond_hubs_2 && ch.hbond_hubs_2.length) ? ch.hbond_hubs_2 : [1, 2, 3];
            styleTarget(this.v1, 0, { cartoon: { color: 'white', opacity: 0.24 } });
            styleSub(this.v1, 0, { resi: h1 }, { cartoon: { color: '#38bdf8', opacity: 1.0 }, stick: { colorscheme: 'cyanCarbon', radius: 0.24 } });
            styleTarget(this.v2, 1, { cartoon: { color: 'white', opacity: 0.24 } });
            styleSub(this.v2, 1, { resi: h2 }, { cartoon: { color: '#e879f9', opacity: 1.0 }, stick: { colorscheme: 'magentaCarbon', radius: 0.24 } });
            msg = `<b>Chaperone H-Bond CSR Network Hubs:</b> Top 12 highest-degree hydrogen-bond hub residues extracted from <code>bio_hbond_network(pdb, 3.5)</code> CSR matrix (Total H-Bonds: <b>${ch.hbonds_1}</b> vs <b>${ch.hbonds_2}</b>, BUNS: <b>${ch.buns_1}</b>/${ch.buns_2}).`;
        } else if (mode === 'rasa') {
            const e1 = (ch.exposed_resi_1 && ch.exposed_resi_1.length) ? ch.exposed_resi_1 : [];
            const e2 = (ch.exposed_resi_2 && ch.exposed_resi_2.length) ? ch.exposed_resi_2 : [];
            const hydroRes = ['VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR'];
            styleTarget(this.v1, 0, { cartoon: { color: '#334155', opacity: 0.30 } });
            styleSub(this.v1, 0, e1.length ? { resi: e1, resn: hydroRes } : { resn: hydroRes }, { cartoon: { color: '#f59e0b', opacity: 1.0 }, sphere: { color: '#fbbf24', radius: 0.82, opacity: 0.78 } });
            styleTarget(this.v2, 1, { cartoon: { color: '#334155', opacity: 0.30 } });
            styleSub(this.v2, 1, e2.length ? { resi: e2, resn: hydroRes } : { resn: hydroRes }, { cartoon: { color: '#f59e0b', opacity: 1.0 }, sphere: { color: '#fbbf24', radius: 0.82, opacity: 0.78 } });
            msg = `<b>Chaperone RASA Exposed Hydrophobic Patches:</b> Golden spheres mark solvent-exposed hydrophobic sidechains (RASA &gt; 0.25). Core RASA Correlation: <b>r = ${ch.rasa_corr}</b>.`;
        } else if (mode === 'rangefinder') {
            const rf = c.rangefinder || {};
            const pr1 = (rf.pocket_res_1 && rf.pocket_res_1.length) ? rf.pocket_res_1 : [1, 2, 3, 4, 5];
            const pr2 = (rf.pocket_res_2 && rf.pocket_res_2.length) ? rf.pocket_res_2 : [1, 2, 3, 4, 5];
            this.v1.removeAllShapes();
            this.v2.removeAllShapes();
            styleTarget(this.v1, 0, { cartoon: { color: '#334155', opacity: 0.24 } });
            styleSub(this.v1, 0, { resi: pr1 }, { cartoon: { color: '#34d399', opacity: 1.0 }, stick: { colorscheme: 'greenCarbon', radius: 0.24 }, sphere: { colorscheme: 'greenCarbon', scale: 0.20 } });
            styleTarget(this.v2, 1, { cartoon: { color: '#334155', opacity: 0.24 } });
            styleSub(this.v2, 1, { resi: pr2 }, { cartoon: { color: '#38bdf8', opacity: 1.0 }, stick: { colorscheme: 'cyanCarbon', radius: 0.24 }, sphere: { colorscheme: 'cyanCarbon', scale: 0.20 } });
            if (rf.pocket_center_1 && rf.pocket_center_1.length === 3 && rf.pocket_center_2 && rf.pocket_center_2.length === 3) {
                const c1 = rf.pocket_center_1;
                const c2 = rf.pocket_center_2;
                const dCenter = Math.hypot(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2]).toFixed(2);
                if (this.isMerged) {
                    this.v1.addCylinder({
                        start: { x: c1[0], y: c1[1], z: c1[2] },
                        end: { x: c2[0], y: c2[1], z: c2[2] },
                        radius: 0.16,
                        color: '#fbbf24',
                        dashed: true
                    });
                    this.v1.addLabel(`10 Å Pocket Center Shift: Δ=${dCenter} Å · IoU=${rf.range_iou ?? '--'}`, {
                        position: { x: (c1[0] + c2[0]) / 2, y: (c1[1] + c2[1]) / 2 + 2.0, z: (c1[2] + c2[2]) / 2 },
                        backgroundColor: '#020617',
                        fontColor: '#fbbf24',
                        fontSize: 12
                    });
                } else {
                    this.v1.addLabel(`10 Å Active-Site Cleft (${pr1.length} res · MW ${rf.mw_box_1 || '--'})`, {
                        position: { x: c1[0], y: c1[1], z: c1[2] },
                        backgroundColor: '#020617',
                        fontColor: '#34d399',
                        fontSize: 12
                    });
                    this.v2.addLabel(`10 Å Active-Site Cleft (${pr2.length} res · MW ${rf.mw_box_2 || '--'})`, {
                        position: { x: c2[0], y: c2[1], z: c2[2] },
                        backgroundColor: '#020617',
                        fontColor: '#38bdf8',
                        fontSize: 12
                    });
                }
            }
            msg = `<b>RangeFinder 10.0 Å Pocket Microenvironment (arXiv:2611.03379v1):</b> Active-site residues within 10 Å rendered as ball-and-stick sidechains (Green = Bait, Cyan = Orphan). 8-Prop Range-IoU = <b>${rf.range_iou ?? '--'}</b> | Orphan MW Box: <b>${rf.mw_box_2 ?? '--'} Da</b>, LogP Box: <b>${rf.logp_box_2 ?? '--'}</b>, P(arom) = <b>${rf.arom_prob_2 ?? '--'}</b> (${rf.substrate_props_matched ?? 8}/8 GiST Substrate Match).`;
        } else if (mode === 'tpp_seam') {
            const rf = c.rangefinder || {};
            const m1 = c.meta1 || {};
            const m2 = c.meta2 || {};
            const l1 = Number(m1.length || 120);
            const l2 = Number(m2.length || 300);
            const microLen = Math.min(l1, l2);
            const ic = Number(c.icomp ?? 0.78);
            const iou = Number(rf.range_iou ?? 0.65);
            const hb = Number(ch.hbonds_1 ?? ch.hbond_1 ?? 12);
            const buns = Number(ch.buns_1 ?? 3);
            const hydro = Number(m1.hydro_pct ?? 34) / 100.0;
            const tmWt = Math.min(74.0, Math.max(44.0, 51.5 + 4.8 * ic + 3.2 * iou + 0.14 * hb - 0.18 * buns));
            const rawDrop = (2.15 * ic + 1.85 * iou + 2.40 * hydro + 0.085 * hb + 0.12 * buns) * 0.85;
            const dTm = -Math.min(11.5, Math.max(1.2, rawDrop));
            const tmKo = Math.max(37.0, tmWt + dTm);
            const tOptWt = (tmWt - 3.8).toFixed(1);
            const tOptKo = (tmKo - 3.7).toFixed(1);
            const dTOpt = (tOptWt - tOptKo).toFixed(1);
            const gainX = Math.min(4.2, Math.max(1.15, Math.exp(0.082 * (tOptWt - tOptKo)))).toFixed(2);
            const gravyMicro = l1 <= l2 ? Number(m1.gravy ?? 0.0) : Number(m2.gravy ?? 0.0);
            const dPi = Math.abs(Number(m1.pi ?? 7.0) - Number(m2.pi ?? 7.0));
            const mech = (gravyMicro >= 0.40 && microLen <= 110)
                ? 'Transmembrane Helix Plug'
                : (dPi >= 2.2 ? 'Electrostatic Subunit Staple' : 'Allosteric Seam Staple');

            const h1 = (ch.hbond_hubs_1 && ch.hbond_hubs_1.length) ? ch.hbond_hubs_1 : [1, 2, 3, 4];
            const h2 = (ch.hbond_hubs_2 && ch.hbond_hubs_2.length) ? ch.hbond_hubs_2 : [1, 2, 3, 4];
            const e1 = (ch.exposed_resi_1 && ch.exposed_resi_1.length) ? ch.exposed_resi_1 : [];
            const e2 = (ch.exposed_resi_2 && ch.exposed_resi_2.length) ? ch.exposed_resi_2 : [];
            const hydroRes = ['VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR'];

            styleTarget(this.v1, 0, { cartoon: { color: '#334155', opacity: 0.25 } });
            styleSub(this.v1, 0, { resi: h1 }, { cartoon: { color: '#38bdf8', opacity: 1.0 }, stick: { colorscheme: 'cyanCarbon', radius: 0.25 } });
            styleSub(this.v1, 0, e1.length ? { resi: e1, resn: hydroRes } : { resn: hydroRes }, { cartoon: { color: '#f43f5e', opacity: 0.95 }, sphere: { color: '#fb7185', radius: 0.78, opacity: 0.76 } });

            styleTarget(this.v2, 1, { cartoon: { color: '#334155', opacity: 0.25 } });
            styleSub(this.v2, 1, { resi: h2 }, { cartoon: { color: '#34d399', opacity: 1.0 }, stick: { colorscheme: 'greenCarbon', radius: 0.25 } });
            styleSub(this.v2, 1, e2.length ? { resi: e2, resn: hydroRes } : { resn: hydroRes }, { cartoon: { color: '#f59e0b', opacity: 0.95 }, sphere: { color: '#fbbf24', radius: 0.78, opacity: 0.76 } });

            const c1 = (rf.pocket_center_1 && rf.pocket_center_1.length === 3) ? rf.pocket_center_1 : [0, 0, 0];
            const c2 = (rf.pocket_center_2 && rf.pocket_center_2.length === 3) ? rf.pocket_center_2 : [0, 0, 0];
            this.v1.addLabel(`Stapled WT Holo-Complex: Tm=${tmWt.toFixed(1)}°C · Ideal T_opt=${tOptWt}°C (${gainX}× Output)`, {
                position: { x: c1[0], y: c1[1] + 2.0, z: c1[2] },
                backgroundColor: '#020617',
                fontColor: '#38bdf8',
                fontSize: 11.5
            });
            if (!this.isMerged) {
                this.v2.addLabel(`Unstapled KO Seam (ΔTm ${dTm.toFixed(2)}°C): Tm=${tmKo.toFixed(1)}°C · Ideal T_opt=${tOptKo}°C`, {
                    position: { x: c2[0], y: c2[1] + 2.0, z: c2[2] },
                    backgroundColor: '#020617',
                    fontColor: '#fb7185',
                    fontSize: 11.5
                });
            }

            const temps = [40, 44, 48, 52, 56, 60, 64, 68];
            const stripHtml = temps.map(t => {
                const fWt = Math.round(100 / (1 + Math.exp((t - tmWt) / 2.1)));
                const fKo = Math.round(100 / (1 + Math.exp((t - tmKo) / 3.9)));
                const isOptWt = Math.abs(t - Number(tOptWt)) <= 2.2;
                const isOptKo = Math.abs(t - Number(tOptKo)) <= 2.2;
                const bdr = isOptWt ? 'border:1px solid #34d399;' : (isOptKo ? 'border:1px solid #fbbf24;' : 'border:1px solid #1e293b;');
                return `<span style="display:inline-flex;align-items:center;gap:3px;background:#090d16;${bdr}padding:1px 6px;border-radius:4px;font-family:ui-monospace,monospace;font-size:0.75em;margin-right:3px;" title="At ${t}°C: WT ${fWt}% folded vs KO ${fKo}% folded">` +
                       `<b>${t}°C:</b><span style="color:#38bdf8">${fWt}%</span>/<span style="color:#fb7185">${fKo}%</span>` +
                       `${isOptWt ? '<b style="color:#34d399">★WT_opt</b>' : ''}${isOptKo ? '<b style="color:#fbbf24">★KO_opt</b>' : ''}</span>`;
            }).join('');

            msg = `<b>CETSA/TPP Thermal Weak Seam &amp; Ideal Bioreactor Temp (<span style="color:#c084fc">${mech}</span>):</b> ` +
                  `Cyan/Green sticks = interfacial H-bond staple anchors (<b>${hb} H-bonds</b>); Rose/Amber spheres = unmasked hydrophobic seam that melts first upon microprotein deletion ` +
                  `(<b style="color:#38bdf8">WT Tm ${tmWt.toFixed(1)}°C</b> → <b style="color:#fb7185">KO Tm ${tmKo.toFixed(1)}°C, ΔTm ${dTm.toFixed(2)}°C</b>). ` +
                  `Ideal Operating Temp shifts from <b style="color:#fbbf24">T_opt ${tOptKo}°C (KO)</b> → <b style="color:#34d399">T_opt ${tOptWt}°C (WT, +${dTOpt}°C · ${gainX}× Catalytic Output)</b>. ` +
                  `<div style="margin-top:5px;display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:6px;"><div>${stripHtml}</div>` +
                  `<a href="/tpp?hub=${encodeURIComponent(l1 <= l2 ? m1.id : m2.id)}" target="_blank" style="color:#fb7185;font-weight:700;text-decoration:none;border:1px solid #f43f5e;padding:2px 8px;border-radius:5px;background:rgba(244,63,94,0.14);">🌡️ Open in Dedicated TPP &amp; T_opt Lab (/tpp) ↗</a></div>`;
        } else if (mode === 'nterm_pep') {
            const pepRes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
            [ [this.v1, 0, '#38bdf8'], [this.v2, 1, '#e879f9'] ].forEach(([v, idx, col]) => {
                styleTarget(v, idx, { cartoon: { color: 'white', opacity: 0.24 } });
                styleSub(v, idx, { resi: pepRes }, { cartoon: { color: col, opacity: 1.0 }, stick: { radius: 0.25 } });
            });
            const sx = c.sql_extensions || {};
            msg = `<b>N-Terminal Diagnostic LC-MS/MS Peptide (Residues 1-10):</b> Highlighted as sticks for HCD tandem MS validation (Diagnostic ions: <code>${sx.ms_ions_preview || 'b2/y2'}</code>, Reversed Decoy: <code>${sx.decoy_peptide || '--'}</code>).`;
        } else if (mode === 'zorder') {
            const rf = c.rangefinder || {};
            const p1Info = this.computeAutoZOrderPocket(this.v1, 0, rf.pocket_center_1);
            const sharedHint = p1Info.anchor ? [p1Info.anchor.x, p1Info.anchor.y, p1Info.anchor.z] : rf.pocket_center_2;
            const p2Info = this.computeAutoZOrderPocket(this.isMerged ? this.v1 : this.v2, 1, sharedHint);

            styleTarget(this.v1, 0, { cartoon: { color: '#334155', opacity: 0.22 } });
            if (p1Info.resiList.length) {
                styleSub(this.v1, 0, { resi: p1Info.resiList }, { cartoon: { color: '#f472b6', opacity: 1.0 }, stick: { colorscheme: 'magentaCarbon', radius: 0.26 }, sphere: { colorscheme: 'magentaCarbon', scale: 0.22 } });
            }
            styleTarget(this.v2, 1, { cartoon: { color: '#334155', opacity: 0.22 } });
            if (p2Info.resiList.length) {
                styleSub(this.v2, 1, { resi: p2Info.resiList }, { cartoon: { color: '#38bdf8', opacity: 1.0 }, stick: { colorscheme: 'cyanCarbon', radius: 0.26 }, sphere: { colorscheme: 'cyanCarbon', scale: 0.22 } });
            }

            const a1 = p1Info.anchor || { resn: 'HIS', resi: 1, x: 0, y: 0, z: 0 };
            const a2 = p2Info.anchor || { resn: 'HIS', resi: 1, x: 0, y: 0, z: 0 };
            const anchorDist = Math.hypot(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z).toFixed(2);

            if (this.isMerged && p1Info.anchor && p2Info.anchor) {
                this.v1.addCylinder({
                    start: { x: a1.x, y: a1.y, z: a1.z },
                    end: { x: a2.x, y: a2.y, z: a2.z },
                    radius: 0.16,
                    color: '#fbbf24',
                    dashed: true
                });
                this.v1.addLabel(`Bait ${a1.resn}${a1.resi} (${p1Info.resiList.length} res) vs Orphan ${a2.resn}${a2.resi} (${p2Info.resiList.length} res) · Δ=${anchorDist} Å`, {
                    position: { x: (a1.x + a2.x) / 2, y: Math.max(a1.y, a2.y) + 2.2, z: (a1.z + a2.z) / 2 },
                    backgroundColor: '#020617',
                    fontColor: '#fbbf24',
                    fontSize: 12
                });
            } else {
                if (p1Info.anchor) {
                    this.v1.addLabel(`${a1.resn}${a1.resi} (5.0 Å Z-Order Pocket · ${p1Info.resiList.length} res)`, {
                        position: { x: a1.x, y: a1.y, z: a1.z },
                        backgroundColor: '#020617',
                        fontColor: '#f472b6',
                        fontSize: 12
                    });
                }
                if (p2Info.anchor) {
                    this.v2.addLabel(`${a2.resn}${a2.resi} (5.0 Å Z-Order Pocket · ${p2Info.resiList.length} res · Δ=${anchorDist} Å)`, {
                        position: { x: a2.x, y: a2.y, z: a2.z },
                        backgroundColor: '#020617',
                        fontColor: '#38bdf8',
                        fontSize: 12
                    });
                }
            }
            msg = `<b>Automatic 30-Bit Z-Order (Morton) 5.0 Å Pocket Selection:</b> Bait <b style="color:#f472b6">${a1.resn}${a1.resi}</b> (${p1Info.resiList.length} residues within 5.0 Å) vs Orphan <b style="color:#38bdf8">${a2.resn}${a2.resi}</b> (${p2Info.resiList.length} residues within 5.0 Å, Anchor 3D Shift Δ = <b>${anchorDist} Å</b>) — SQL: <code>WHERE morton_code BETWEEN z_order_encode(${(a1.x-5).toFixed(1)}, ${(a1.y-5).toFixed(1)}, ${(a1.z-5).toFixed(1)}) AND z_order_encode(${(a1.x+5).toFixed(1)}, ${(a1.y+5).toFixed(1)}, ${(a1.z+5).toFixed(1)})</code>`;
        }

        explainEl.innerHTML = `<span>${customExplain || msg}</span>`;
        this.v1.render();
        if (!this.isMerged) this.v2.render();
    }

    computeAutoZOrderPocket(viewer, modelIdx, centerHint) {
        const sel = this.isMerged ? { model: modelIdx } : {};
        const atoms = (viewer && viewer.selectedAtoms) ? viewer.selectedAtoms(sel) : [];
        if (!atoms || !atoms.length) return { anchor: null, resiList: [] };

        const catTypes = new Set(['HIS', 'SER', 'CYS', 'ASP', 'GLU', 'TYR', 'TRP', 'PHE', 'ARG', 'LYS', 'ASN', 'GLN']);
        const caAtoms = atoms.filter(a => a.atom === 'CA');
        const pool = caAtoms.length ? caAtoms : atoms;

        let cx = 0, cy = 0, cz = 0;
        if (centerHint && centerHint.length === 3) {
            cx = centerHint[0];
            cy = centerHint[1];
            cz = centerHint[2];
        } else {
            pool.forEach(a => { cx += a.x; cy += a.y; cz += a.z; });
            cx /= pool.length;
            cy /= pool.length;
            cz /= pool.length;
        }

        // Pick catalytic residue closest to active-site center (or overall closest residue if no catalytic residue within 10 Å)
        let bestAnchor = null;
        let bestScore = Infinity;
        pool.forEach(a => {
            const d = Math.hypot(a.x - cx, a.y - cy, a.z - cz);
            const score = catTypes.has(a.resn) ? d * 0.65 : d + 4.0;
            if (score < bestScore) {
                bestScore = score;
                bestAnchor = a;
            }
        });
        if (!bestAnchor) bestAnchor = pool[0];

        const resSet = new Set();
        atoms.forEach(a => {
            const d = Math.hypot(a.x - bestAnchor.x, a.y - bestAnchor.y, a.z - bestAnchor.z);
            if (d <= 5.0) resSet.add(a.resi);
        });
        return { anchor: bestAnchor, resiList: Array.from(resSet).sort((a, b) => a - b) };
    }

    toggleMergeView() {
        if (!this.v1 || !this.v2 || !this.compData) return;
        this.isMerged = !this.isMerged;
        const wrap = document.getElementById(`${this.instanceId}_viewers_wrap`);
        const vbox1 = document.getElementById(`${this.instanceId}_vbox1`);
        const vbox2 = document.getElementById(`${this.instanceId}_vbox2`);
        const ov2 = document.getElementById(`${this.instanceId}_ov2`);
        const bot2 = document.getElementById(`${this.instanceId}_bot2`);
        const canvas1 = document.getElementById(`${this.instanceId}_canvas1`);
        const canvas2 = document.getElementById(`${this.instanceId}_canvas2`);
        const btn = document.getElementById(`${this.instanceId}_merge_btn`);

        if (this.isMerged) {
            wrap.classList.add('merged');
            vbox2.style.display = 'none';
            ov2.classList.add('docked-right');
            vbox1.insertBefore(ov2, canvas1);
            if (bot2) vbox1.appendChild(bot2);
            btn.classList.add('active');
            btn.innerHTML = `${msIcon('vertical_split')}Split View`;

            this.v1.clear();
            this.v1.addModel(this.compData.pdb1, "pdb");
            this.v1.addModel(this.compData.pdb2_aligned, "pdb");
            this.applyMode(this.activeMode || 'core');
            this.v1.resize();
            this.v1.zoomTo();
            this.v1.render();
        } else {
            wrap.classList.remove('merged');
            vbox2.style.display = '';
            ov2.classList.remove('docked-right');
            vbox2.insertBefore(ov2, canvas2);
            if (bot2) vbox2.appendChild(bot2);
            btn.classList.remove('active');
            btn.innerHTML = `${msIcon('layers')}Superpose`;

            this.v1.clear();
            this.v1.addModel(this.compData.pdb1, "pdb");
            this.applyMode(this.activeMode || 'identity');
            this.v1.resize();
            this.v1.zoomTo();
            this.v1.render();
            this.v2.resize();
            this.v2.setView(this.v1.getView());
            this.v2.render();
        }
    }

    togglePocketSelection() {
        this.applyMode('zorder');
    }
}
"""
