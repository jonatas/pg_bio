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

## 2. Aromaticity Perception (Hückel's $4n+2$ Rule) & VF2 Subgraph Matching (`@>` / `smiles_contains`)

A classic trap in cheminformatics is inserting a molecule in **Kekulé form** (`C1=CC=CC=C1`) and querying it in **aromatic form** (`c1ccccc1`), or storing a molecule as a **salt** (`[Na+].O=C([O-])C1=CC=CC=C1`).

### How `pg_bio` Normalizes Molecular Graphs
Immediately after lexing a SMILES string, `pg_bio` runs:
1. **Salt Stripping & Counter-Ion Neutralization:** Splits on `.` connected components, keeps the largest heavy-atom component, and neutralizes counter-ion salt charges (`[O-].[Na+]` $\rightarrow$ `OH`, `[NH3+].[Cl-]` $\rightarrow$ `NH2`).
2. **1,3-Prototropic Tautomer Canonicalization:** Normalizes aliphatic keto-enol (`C=CO` $\rightarrow$ `CC=O`) and imidic-amide (`CC(=N)O` $\rightarrow$ `CC(=O)N`) tautomers.
3. **Hückel’s $4n+2$ Cycle Perception:** Evaluates $\pi$-electrons around every candidate cycle (`1` for `=C-`/`=N-`, `2` for lone-pair heteroatoms `[nH]`/`O`/`S` and `[C-]`, `0` for exocyclic `C(=O)` in Pyridone/Caffeine, `-1` abort for saturated $sp^3$ atoms) and promotes cycles with $4n+2$ ($6, 10, 14, \dots$) $\pi$-electrons to `BOND_AROMATIC`.
4. **VF2 Subgraph Isomorphism (`@>` and `smiles_contains`):** Matches query vertices `Q[i]` to target vertices `T[j]` by `(atomic_num, is_aromatic, degree)` and verifies bond-type consistency across all incident edges.

```sql
-- Query Kekulé Aspirin with aromatic Benzene (`c1ccccc1`) using the `@>` operator:
SELECT parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') @> parse_smiles('c1ccccc1') AS has_benzene,
       smiles_contains(parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O'), 'c1ccccc1C(=O)O') AS has_salicylate;
-- Result: true, true
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
