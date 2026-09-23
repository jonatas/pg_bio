-- Experiment 1: Biological Homology & Spatial Search inside PostgreSQL
-- This script demonstrates how pg_bio allows scientists to perform complex 
-- bioinformatics workflows natively in the database.

-- Enable the extension
CREATE EXTENSION IF NOT EXISTS pg_bio;

-- ====================================================================
-- PART A: VECTOR EMBEDDING SEARCH (Finding functional homologues)
-- ====================================================================
CREATE TABLE protein_database (
    id serial PRIMARY KEY,
    name text,
    sequence text,
    embedding real[]
);

-- Insert the TeaFlon domains (using our mocked ESM embedding generator)
INSERT INTO protein_database (name, sequence, embedding) VALUES
('Fluoroacetate Dehalogenase', 'MFEGFERRLVD', get_esm_embedding('MFEGFERRLVD')),
('Class II Hydrophobin', 'MQFFAVALFA', get_esm_embedding('MQFFAVALFA')),
('Amelogenin Tag', 'WPSTDKTKREEVD', get_esm_embedding('WPSTDKTKREEVD'));

-- Run the experiment: 
-- Which protein in our database is most functionally similar to a mutated Dehalogenase?
WITH synthetic_query AS (
    SELECT get_esm_embedding('MFEGFERRLV') as emb -- Slight mutation
)
SELECT 
    p.name, 
    embedding_cosine_distance(p.embedding, q.emb) as distance
FROM 
    protein_database p, synthetic_query q
ORDER BY 
    distance ASC;

-- ====================================================================
-- PART B: 3D SPATIAL INDEXING (Finding atoms in a binding pocket)
-- ====================================================================
CREATE TABLE active_site_atoms (
    atom_id serial PRIMARY KEY,
    coord ResidueCoord,
    z_index bigint
);

-- Insert 3D coordinates for the Dehalogenase active site
INSERT INTO active_site_atoms (coord, z_index) VALUES
( '{"x": 10.5, "y": 12.0, "z": 8.1, "name": "TRP150"}', residue_z_index('{"x": 10.5, "y": 12.0, "z": 8.1, "name": "TRP"}') ),
( '{"x": 11.0, "y": 12.5, "z": 8.5, "name": "ASP105"}', residue_z_index('{"x": 11.0, "y": 12.5, "z": 8.5, "name": "ASP"}') ),
( '{"x": 50.0, "y": -10.0, "z": 100.0, "name": "GLY10"}', residue_z_index('{"x": 50.0, "y": -10.0, "z": 100.0, "name": "GLY"}') ); -- Far away atom

-- Create a B-Tree index on the Z-Order curve!
-- This instantly maps 3D spatial space into a highly optimized 1D binary tree.
CREATE INDEX idx_spatial_z_order ON active_site_atoms (z_index);

-- Experiment: Find all atoms near a specific spatial coordinate (e.g., a docked Teflon molecule)
-- By querying the 1D Z-index, Postgres can jump directly to the correct spatial sector without Math.sqrt()!
SELECT 
    (coord).name, 
    distance_angstroms(coord, '{"x": 10.0, "y": 11.5, "z": 8.0, "name": "PTFE"}') as true_distance
FROM 
    active_site_atoms
WHERE 
    z_index BETWEEN residue_z_index('{"x": 9.0, "y": 10.0, "z": 7.0, "name": ""}') 
                AND residue_z_index('{"x": 12.0, "y": 13.0, "z": 9.0, "name": ""}');
