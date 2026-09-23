-- Experiment 2: Sparse Attention Matrices & Structural Alignment
-- This script demonstrates querying AI Attention Weights natively in PostgreSQL.

CREATE EXTENSION IF NOT EXISTS pg_bio;

CREATE TABLE protein_attention_maps (
    protein_id serial PRIMARY KEY,
    name text,
    attention_data SparseAttentionMap
);

-- Imagine we extract the structural attention weights from ESMFold's neural network.
-- We only save interactions with an attention weight > 0.5 to compress the matrix.
INSERT INTO protein_attention_maps (name, attention_data) VALUES
(
    'Fluoroacetate Dehalogenase (Active Site Interactions)',
    create_sparse_map(
        298,
        ARRAY[105, 105, 150],
        ARRAY[150, 20, 88],
        ARRAY[0.98, 0.65, 0.88]::real[]
    )
);

-- Experiment:
-- A biologist wants to know which amino acids are physically forcing 
-- Residue 105 (Aspartate, part of our active site) into its unique shape.
-- By querying the compressed AI attention map, Postgres instantly returns 
-- the top 2 residues with the strongest evolutionary/spatial interaction.

SELECT 
    name as protein,
    get_top_interacting_residues(attention_data, 105, 2) as strongest_interactions
FROM 
    protein_attention_maps;
