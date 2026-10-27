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

## 3. Parsing PDB directly from Strings (APIs)
If you are fetching PDB structures dynamically from an API (like ESMFold) or storing them in a TEXT column, you can use `bio_parse_pdb` to parse them entirely in memory without hitting the disk:

```sql
SELECT atom_name, res_name, res_id, x, y, z 
FROM bio_parse_pdb('ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 49.33           N  ...');
```

## 4. Structural Alignment & RMSD
Measuring the structural similarity between two protein folds is fundamental in structural biology. You can compute the Root Mean Square Deviation (RMSD) of two coordinate arrays natively:

```sql
WITH conf1 AS (
    SELECT ARRAY[create_residue_coord(27.34, 24.43, 2.61, 'N')] as atoms
),
conf2 AS (
    SELECT ARRAY[create_residue_coord(27.50, 24.40, 2.60, 'N')] as atoms
)
SELECT bio_rmsd(conf1.atoms, conf2.atoms) as rmsd_angstroms
FROM conf1, conf2;
```
