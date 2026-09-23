---
name: pg_bio_querying
description: Teaches the agent how to construct SQL queries using the custom pg_bio PostgreSQL extension for Z-Order spatial searches and AI vector embeddings.
---

# Querying the pg_bio Extension

When asked to perform bioinformatics searches, structural queries, or spatial alignment inside the database, use the `pg_bio` extension capabilities.

## 1. Z-Order Spatial Indexing (Bounding Boxes)
To find atoms within a radius of a specific target, do NOT use brute-force `sqrt()` math. Use the Z-Order curve B-Tree index.

**SQL Pattern:**
```sql
WITH bounds AS (
    SELECT 
        residue_z_index('{"x": MIN_X, "y": MIN_Y, "z": MIN_Z, "name": ""}') as min_z,
        residue_z_index('{"x": MAX_X, "y": MAX_Y, "z": MAX_Z, "name": ""}') as max_z
)
SELECT * FROM your_table, bounds 
WHERE z_index BETWEEN min_z AND max_z
AND distance_angstroms(coord, '{"x": TGT_X, "y": TGT_Y, "z": TGT_Z, "name": ""}') <= RADIUS;
```

## 2. ESM Vector Embeddings
To find structurally or functionally similar proteins, use the `embedding_cosine_distance` operator. The distance represents `1 - similarity`. Lower is closer.

**SQL Pattern:**
```sql
SELECT name, embedding_cosine_distance(embedding, get_esm_embedding('QUERY_SEQUENCE')) as distance
FROM protein_database
ORDER BY distance ASC LIMIT 5;
```

## 3. Sparse Attention Maps
To find which amino acids have the highest AI attention weights against a target residue, use the `get_top_interacting_residues` function against a `SparseAttentionMap` type.

**SQL Pattern:**
```sql
SELECT get_top_interacting_residues(attention_data, TARGET_RESIDUE_ID, LIMIT_INT) 
FROM protein_attention_maps;
```

## Testing
Always run queries using standard `psql` connected to the `pgrx` instance (e.g., via `cargo pgrx run pg18` or `psycopg` scripts).
