# Tutorial 6: Building SMILES Expressions & Pattern Matching

SMILES (Simplified Molecular-Input Line-Entry System) is a line notation for describing chemical structures. This guide will help you understand the rules of building SMILES expressions and how `pg_bio` natively matches and compares these structures.

## 1. Core SMILES Rules

### Atoms and Bonds
- **Atoms** are represented by their standard chemical symbols (e.g., `C` for Carbon, `O` for Oxygen, `N` for Nitrogen).
- **Aromatic Atoms** are represented by lowercase letters (e.g., `c`, `n`, `o`).
- **Bonds** between atoms are implicit if they are single bonds. 
  - Double bonds are represented by `=` (e.g., `C=O` is a carbonyl group).
  - Triple bonds are represented by `#` (e.g., `C#N` is a nitrile group).

### Branches
Branches from a main molecular chain are enclosed in parentheses `()`.
- Example: **Isobutane** `CC(C)C`
  - A central Carbon `C` is attached to a methyl group `(C)`, while the main chain continues.

### Rings
Rings are built by breaking one bond in the ring and labeling the two adjacent atoms with a matching number.
- Example: **Cyclohexane** `C1CCCCC1`
  - The `1` after the first and last `C` indicates they are bonded together, forming a 6-member ring.
- Example: **Benzene** `c1ccccc1` (or `C1=CC=CC=C1`)
  - A 6-member aromatic ring.

## 2. Structural Matching (`smiles_contains`)

The `smiles_contains` function allows you to query if a molecular graph contains a specific substructure.

### How it Works
`pg_bio` parses both the target compound and the query substructure into Native Rust graphs (via `purr`). It performs a subgraph isomorphism check. This means it doesn't just look for substring matches in text; it mathematically ensures the exact atoms and bonds exist in the same topological arrangement.

```sql
-- Query: Does Aspirin contain a Benzene ring?
SELECT smiles_contains(
    parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O'), -- Target
    'C1=CC=CC=C1'                            -- Substructure
) AS has_benzene;
-- Result: true
```

### Advanced Matching Example
If you look for a phenol group (`OC1=CC=CC=C1`), it will correctly identify whether the `O` is directly attached to the aromatic ring.

```sql
-- Query: Does Aspirin contain a Phenol group?
SELECT smiles_contains(
    parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O'), 
    'OC1=CC=CC=C1'
) AS has_phenol;
-- Result: true (Aspirin is technically an ester of phenol!)
```

## 3. Morgan Fingerprints & Topological Matching

Sometimes you don't want an *exact* subgraph match. You want a statistical **Similarity Score**. This is where Morgan Fingerprints (ECFPs) come into play.

### How the Expression builds the Fingerprint:
1. **Initial Invariants**: `pg_bio` assigns an integer "hash" to each atom based on its atomic number, charge, and degree (number of neighbors).
2. **Radius Expansion**: It iteratively expands outwards. A radius of `2` means it hashes the paths linking an atom to its neighbors, and its neighbors' neighbors.
3. **Folding**: The millions of unique path hashes are folded (modulo operations) into a fixed 1024-bit vector `Vec<f32>`.

### Tanimoto Similarity Operator (`%`)
The `%` operator compares two 1024-dimensional continuous arrays.
The mathematical formula for Tanimoto (Jaccard) similarity is:
`Intersection(A, B) / Union(A, B)`

If two molecules share the exact same chemical fragments, their similarity approaches `1.0`.

```sql
-- Compare Toluene (Benzene with a Methyl group) to Benzene
SELECT (smiles_to_fingerprint('CC1=CC=CC=C1') % smiles_to_fingerprint('C1=CC=CC=C1')) AS similarity;
-- Result: ~0.1143
```

Because `pg_bio` natively manages these vectors inside Postgres, it completely avoids the network I/O penalty of Python-based dataframe serializations, executing chemical comparisons 5-10x faster natively in bulk.
