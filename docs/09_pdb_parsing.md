# Tutorial 9: In-Database 3D Structure Parsing

Molecular visualization and structural modeling pipelines usually read PDB or mmCIF files in Python/C++, serialize them into arrays, and upload them slowly into databases.

With `pg_bio` and the `pdbtbx` parser, you can ingest raw PDB coordinate files natively in SQL, map them to Z-Order Morton curves for spatial indexing, and dock them—all in one query.

## 1. Reading PDB Files natively

```sql
SELECT atom_id, residue, x, y, z 
FROM parse_pdb_file('/path/to/protein.pdb')
LIMIT 5;
```

## 2. Instant Z-Order Spatial Indexing

You can immediately transform these raw files into highly optimized 3D structures.

```sql
CREATE TABLE my_protein_atoms (
    atom_id int,
    residue text,
    coord ResidueCoord,
    z_index bigint
);

INSERT INTO my_protein_atoms (atom_id, residue, coord, z_index)
SELECT 
    atom_id, 
    residue,
    ROW(x, y, z, residue)::ResidueCoord,
    residue_z_index(ROW(x, y, z, residue)::ResidueCoord)
FROM parse_pdb_file('/path/to/protein.pdb');
```

This bypasses massive network overhead, letting PostgreSQL act directly as a bioinformatics computing engine!
