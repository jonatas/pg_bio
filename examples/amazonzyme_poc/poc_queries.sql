
-- ============================================================================
-- AMAZONZYME x pg_bio: IN-SILICO BIOPROSPECTING & ENZYME ENGINEERING POC
-- Targeted directly at the published work of AmazonZyme's founders (UFAM):
--   1. Dra. Pamella Suely Santa Rosa Pimentel (CEO)
--   2. Prof. Dr. Edson Junior do Carmo (Co-founder, LAB-BIOSIN / UFAM)
--   3. MSc/PhD Cand. Samanta Gabriela Souza da Silva (Co-founder, PPGBIOTEC / UFAM)
-- ============================================================================

-- MODULE 1: 3D Z-ORDER ACTIVE-SITE INDEXING OF ADE1 / LIP4 (CLONE FOS4, PDB: 8V16)
-- Isolated from Amazonian Dark Earth (Terra Preta de Índio, Presidente Figueiredo, AM)
-- Reference: Batalha Jr., Carmo, Astolfi-Filho et al. (REMUNOM, 2026) & Vinces et al. (Biochemistry, 2024)
DROP TABLE IF EXISTS amazonzyme_ade1_atoms;
CREATE TABLE amazonzyme_ade1_atoms AS
SELECT 
    atom_id,
    residue,
    x, y, z,
    create_residue_coord(x, y, z, residue) AS coord,
    z_order_encode(x, y, z) AS z_index
FROM parse_pdb_file('/Users/jonatas/code/pg_bio/examples/amazonzyme_poc/8V16_chainA.pdb');

CREATE INDEX idx_ade1_z_order ON amazonzyme_ade1_atoms (z_index);

-- 1A. Verify 8V16 Chain A ingestion and 3D Center of Mass in Postgres:
SELECT 
    count(*) AS chain_a_atoms_8v16,
    center_of_mass(array_agg(coord))::text AS center_of_mass_3d
FROM amazonzyme_ade1_atoms;

-- 1B. Automatically detect the 3D Catalytic Triad (Ser94 - Asp217 - His245) in 8V16 via SQL:
-- Demonstrates how pg_bio detects active-site 3D geometry in <5ms even when BLASTP sequence
-- identity is only 33% (as reported for Clone FosLip4 / Ade1 vs PDB 3OM8).
SELECT 
    s.residue || ' (Atom ' || s.atom_id || ')' AS nucleophile_ser94,
    h.residue || ' (Atom ' || h.atom_id || ')' AS general_base_his245,
    d.residue || ' (Atom ' || d.atom_id || ')' AS acid_asp217,
    round(distance_angstroms(s.coord, h.coord)::numeric, 2) AS ser_to_his_angstroms,
    round(distance_angstroms(h.coord, d.coord)::numeric, 2) AS his_to_asp_angstroms,
    is_steric_clash(h.coord, d.coord, 2.5) AS has_steric_clash
FROM amazonzyme_ade1_atoms h
JOIN amazonzyme_ade1_atoms d ON h.residue = 'HIS' AND d.residue = 'ASP' AND distance_angstroms(h.coord, d.coord) < 3.1
JOIN amazonzyme_ade1_atoms s ON s.residue = 'SER' AND distance_angstroms(s.coord, h.coord) < 4.5
ORDER BY his_to_asp_angstroms ASC, ser_to_his_angstroms ASC
LIMIT 3;

-- MODULE 2: LIVE UNIPROT ENRICHMENT OF FOUNDERS' BENCHMARK ENZYMES INSIDE SQL
-- Matches:
--   - Ade1 / Lip4 (PDB 8V16, 29.3 kDa) from Terra Preta de Índio
--   - Penicillium citrinum LMI01 35 kDa Xylanase/Cellulase band (Pamella's PhD Thesis, UFAM 2019)
--   - Trichoderma reesei QM9414 Cellulase benchmark (UFAM1 & UFAM2 cocktails)
--   - Aspergillus awamori Glucoamylase (PDB 1GLM, Edson's MSc Thesis & Samanta's 2026 Ethanol paper)
--   - Bacillus subtilis Alpha-Amylase (PDB 1BAG, Samanta & Edson's 2026 Ethanol paper)
--   - Ktedonobacter racemifer Carboxylesterase B (Clone FosLip3 ORF 00017, Terra Preta de Índio 2026)
SELECT 
    target_role,
    id AS uniprot_id,
    left(organism, 36) AS organism,
    len AS aa_len,
    kda AS calc_mw_kda,
    entropy AS shannon_entropy,
    has_active_motif
FROM (
    SELECT '1. Ade1 / Lip4 (TPI Clone Fos4, PDB 8V16)' AS target_role, id, organism,
           length(sequence) AS len,
           round(molecular_weight(sequence)::numeric/1000.0, 1) AS kda,
           round(sequence_entropy(sequence)::numeric, 3) AS entropy,
           prosite_match(sequence, 'G-x-S-x-G') AS has_active_motif
    FROM bio_search_uniprot('xref:pdb-8V16')
    UNION ALL
    SELECT '2. P. citrinum 35kDa Band (Cocktail UFAM2)', id, organism,
           length(sequence),
           round(molecular_weight(sequence)::numeric/1000.0, 1),
           round(sequence_entropy(sequence)::numeric, 3),
           prosite_match(sequence, '[LIVMF]-x(2)-[DNE]')
    FROM (SELECT * FROM bio_search_uniprot('accession:B1B533') LIMIT 1) p
    UNION ALL
    SELECT '3. T. reesei QM9414 Endoglucanase (UFAM1/2)', id, organism,
           length(sequence),
           round(molecular_weight(sequence)::numeric/1000.0, 1),
           round(sequence_entropy(sequence)::numeric, 3),
           prosite_match(sequence, '[LIVMF]-x(2)-[DNE]')
    FROM (SELECT * FROM bio_search_uniprot('accession:P07982') LIMIT 1) t
    UNION ALL
    SELECT '4. A. awamori Glucoamylase (PDB 1GLM, SSF)', id, organism,
           length(sequence),
           round(molecular_weight(sequence)::numeric/1000.0, 1),
           round(sequence_entropy(sequence)::numeric, 3),
           prosite_match(sequence, 'W-E-E')
    FROM bio_search_uniprot('xref:pdb-1GLM')
    UNION ALL
    SELECT '5. B. subtilis Alpha-Amylase (PDB 1BAG)', id, organism,
           length(sequence),
           round(molecular_weight(sequence)::numeric/1000.0, 1),
           round(sequence_entropy(sequence)::numeric, 3),
           prosite_match(sequence, '[LIVMF]-x-D-x-V-x-N-H')
    FROM bio_search_uniprot('xref:pdb-1BAG')
    UNION ALL
    SELECT '6. K. racemifer Esterase B (TPI FosLip3)', id, organism,
           length(sequence),
           round(molecular_weight(sequence)::numeric/1000.0, 1),
           round(sequence_entropy(sequence)::numeric, 3),
           prosite_match(sequence, 'G-x-S-x-G')
    FROM (SELECT * FROM bio_search_uniprot('accession:D6TJC3') LIMIT 1) k
) founders_bench;

-- MODULE 3: OVERCOMING THE "HARPIA + BLAST" BOTTLENECK IN DARK PROTEOME BIOPROSPECTING
-- Using Hybrid Structural Embedding + Smith-Waterman Re-ranking (<~>) across 2,347,310 proteins
-- to recover extremophile Cellulases, Xylanases (FosXil), and Acid-Stable Alpha-Amylases (SSF/RSDE)
SELECT 
    'Cellulase (pH 4.2 / 60C LMI01-like)' AS enzyme_class,
    uniprot_id AS discovered_orphan_id,
    left(name, 75) AS orphan_annotation,
    round(distance::numeric, 4) AS esm_vector_dist,
    round(hybrid_score::numeric, 4) AS hybrid_bio_score
FROM pg_bio_search_homologs('A0A2U9ICY9', 0.35, 0.3, 2.0, 2, 50)
UNION ALL
SELECT 
    'Acid Alpha-Amylase (Single-Tank SSF)',
    uniprot_id,
    left(name, 75),
    round(distance::numeric, 4),
    round(hybrid_score::numeric, 4)
FROM pg_bio_search_homologs('A0A2U9IBL3', 0.35, 0.3, 2.0, 2, 50)
UNION ALL
SELECT 
    'Xylanase (Rescuing FosXil Dropout)',
    uniprot_id,
    left(name, 75),
    round(distance::numeric, 4),
    round(hybrid_score::numeric, 4)
FROM pg_bio_search_homologs('A0A2U1S825', 0.35, 0.3, 2.0, 2, 50);
