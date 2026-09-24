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
( create_residue_coord(10.5, 12.0, 8.1, 'TRP150'), residue_z_index(create_residue_coord(10.5, 12.0, 8.1, 'TRP')) ),
( create_residue_coord(11.0, 12.5, 8.5, 'ASP105'), residue_z_index(create_residue_coord(11.0, 12.5, 8.5, 'ASP')) ),
( create_residue_coord(50.0, -10.0, 100.0, 'GLY10'), residue_z_index(create_residue_coord(50.0, -10.0, 100.0, 'GLY')) ); -- Far away atom

-- Create a B-Tree index on the Z-Order curve!
-- This instantly maps 3D spatial space into a highly optimized 1D binary tree.
CREATE INDEX idx_spatial_z_order ON active_site_atoms (z_index);

-- Experiment: Find all atoms near a specific spatial coordinate (e.g., a docked Teflon molecule)
-- By querying the 1D Z-index, Postgres can jump directly to the correct spatial sector without Math.sqrt()!
SELECT 
    (coord).name, 
    distance_angstroms(coord, create_residue_coord(10.0, 11.5, 8.0, 'PTFE')) as true_distance
FROM 
    active_site_atoms
WHERE 
    z_index BETWEEN residue_z_index(create_residue_coord(9.0, 10.0, 7.0, '')) 
                AND residue_z_index(create_residue_coord(12.0, 13.0, 9.0, ''));
