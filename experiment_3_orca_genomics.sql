-- Experiment 3: Genomic Architecture Validation (ORCA Deep Learning Model)
-- This script demonstrates how pg_bio can ingest and query Hi-C contact maps 
-- predicted by DNA/RNA deep learning models (like the Zhou Lab ORCA model).

CREATE EXTENSION IF NOT EXISTS pg_bio;

CREATE TABLE genomic_predictions (
    plasmid_id serial PRIMARY KEY,
    sequence_name text,
    -- We reuse our SparseAttentionMap type to act as a Compressed Sparse Row (CSR) 
    -- for the ORCA Hi-C DNA Contact Map. It perfectly models genomic looping!
    orca_hic_map SparseAttentionMap
);

-- Imagine we ran the teaflon_orca_input.fasta through the local PyTorch ORCA model.
-- The model predicts a 3D structural contact map for our 1,352 base pair plasmid.
-- We compress the probabilities and only store interactions with a high looping score (> 0.7).
INSERT INTO genomic_predictions (sequence_name, orca_hic_map) VALUES
(
    'TeaFlon Synthetic Plasmid',
    '{"sequence_length": 1352, "entries": [
        {"source_residue": 25, "target_residue": 1300, "weight": 0.92},  -- Promoter looping to Terminator!
        {"source_residue": 25, "target_residue": 500,  "weight": 0.15},  -- Weak interaction
        {"source_residue": 800, "target_residue": 810, "weight": 0.99}   -- Local tight folding (knot)
    ]}'
);

-- ====================================================================
-- THE BIOINFORMATICS QUERY
-- ====================================================================

-- Experiment:
-- We need to ensure our T7 Promoter (located around base pair 25) is structurally 
-- accessible and optimally folded to initiate transcription of the TeaFlon enzyme.
-- We ask Postgres to instantly scan the ORCA 3D genomic prediction and find the 
-- top DNA regions physically interacting with the promoter in 3D space.

SELECT 
    sequence_name as synthetic_construct,
    get_top_interacting_residues(orca_hic_map, 25, 2) as closest_physical_dna_contacts
FROM 
    genomic_predictions;

-- Expected Output:
-- {1300, 500}
-- Interpretation: Base pair 1300 (the T7 Terminator) is physically touching the T7 Promoter (Base pair 25).
-- This indicates ORCA predicts a healthy DNA loop structure, allowing RNA Polymerase to efficiently 
-- cycle transcription of our enzyme!
