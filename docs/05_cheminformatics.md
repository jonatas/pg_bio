# Tutorial 5: Cheminformatics & Small Molecules

While proteins are the biological machines of a cell, **small molecules** (like pharmaceutical drugs, metabolites, and ligands) are what interact with them. 

The universal language used to represent small chemical structures in a database is called **SMILES** (Simplified Molecular-Input Line-Entry System). 

In this tutorial, we will use our new `pg_bio` native `Molecule` type to parse, store, and query these chemical structures directly inside PostgreSQL, laying the foundation for in-database docking algorithms.

## 1. The Normalized `Molecule` Data Type
Instead of storing raw chemistry strings in `TEXT` columns where Kekulé (`C1=CC=CC=C1`) and aromatic (`c1ccccc1`) forms fail to match, `pg_bio` provides the native `Molecule` PostgreSQL type backed by a 4-stage normalization pipeline:
1. **Connected-Component Salt Stripping:** Splits disconnected `.` fragments (e.g., `CC(=O)[O-].[Na+]`) and retains the largest heavy-atom parent molecule while neutralizing counter-ion charges (`CC(=O)O`).
2. **Prototropic Tautomer Normalization:** Canonicalizes aliphatic 1,3-keto-enol (`C=CO` $\rightarrow$ `CC=O`) and 1,3-imidic-amide (`CC(=N)O` $\rightarrow$ `CC(=O)N`) tautomers.
3. **Hückel’s $4n+2$ Aromaticity Perception:** Evaluates every simple and fused cycle ($3..=10$) using pure integer $\pi$-electron arithmetic (`(pi_sum >= 6) && ((pi_sum - 2) & 3) == 0`), promoting conjugated $4n+2$ rings to `BOND_AROMATIC`.
4. **Weisfeiler-Lehman (WL) 64-Bit Hash & 1024-Bit ECFP4 Bitset:** Computes vertex-permutation-invariant `wl_hash` (`i64`) and `ecfp_bits` (`16 x u64`) via `fmix64` and sorted neighbor hashes.

```sql
-- Kekulé Benzene and Aromatic Benzene compile to the exact same WL hash and graph!
SELECT molecule_wl_hash('C1=CC=CC=C1') AS kekule_hash,
       molecule_wl_hash('c1ccccc1')    AS aromatic_hash,
       molecules_are_identical('C1=CC=CC=C1', 'c1ccccc1') AS identical;
```
**Output:**
```text
    kekule_hash     |   aromatic_hash    | identical 
--------------------+--------------------+-----------
 466190659499846863 | 466190659499846863 | t
```

If we pass an invalid chemical graph (like an unclosed parenthesis `C(C`), the parser flags `is_valid = false`:
```sql
SELECT (parse_smiles('C(C')).is_valid;
-- Result: f
```

## 2. Exact VF2 Substructure Searching (`@>` and `smiles_contains`)
Finding specific chemical scaffolds (like an aromatic benzene or salicylate core) across candidate drugs and enzyme substrates requires exact **VF2 Subgraph Isomorphism** after Hückel $4n+2$ normalization.

You can use either `smiles_contains(target_mol, query_smiles)` or the native **`@>`** containment operator (`Molecule @> Molecule`):

```sql
-- Query Kekulé Aspirin using an aromatic Benzene query (`c1ccccc1`)
SELECT parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') @> parse_smiles('c1ccccc1') AS contains_benzene,
       smiles_contains(parse_smiles('CC(=O)Oc1ccccc1C(=O)O'), 'C1=CC=CC=C1') AS reverse_kekule_match;
```
**Output:**
```text
 contains_benzene | reverse_kekule_match 
------------------+----------------------
 t                | t
```

## 3. Hardware-Popcount ECFP4 Fingerprints & Tanimoto Similarity (`%`)
`pg_bio` supports the **`%` Tanimoto Similarity Operator** on both `Molecule % Molecule` (using hardware 64-bit `POPCNT` across the 16 `u64` words of `ecfp_bits`) and `real[] % real[]` (`smiles_to_fingerprint`):

```sql
-- 1. Direct hardware-popcount Tanimoto on Molecule (Kekulé vs. Aromatic = 1.0)
SELECT parse_smiles('C1=CC=CC=C1') % parse_smiles('c1ccccc1') AS benzene_identity,
       ROUND((parse_smiles('C1=CC=CC=C1') % parse_smiles('CC1=CC=CC=C1'))::numeric, 4) AS benzene_vs_toluene;
```
**Output:**
```text
 benzene_identity | benzene_vs_toluene 
------------------+--------------------
                1 |             0.2727
```

### Virtual Screening Workflow
This native integration allows us to build powerful virtual screening queries directly in SQL, avoiding the expensive overhead of transferring millions of rows to Python.

```sql
-- Ultimate SQL Vision: Find all library compounds with > 85% structural similarity to our target scaffold
SELECT compound_name, smiles                                              
FROM library_compounds                                                    
WHERE smiles_to_fingerprint(smiles) % smiles_to_fingerprint('CC(=O)OC1=CC=CC=C1C(=O)O') > 0.85;
```

## Conclusion
With native Spatial Math for 3D proteins, and native Cheminformatics for Small Molecules (SMILES parsing, Morgan Fingerprints, and Tanimoto Similarity), we now have all the primitives needed to perform end-to-end biological simulations and drug discovery inside a SQL query. 

You can now cross-reference billions of ChEMBL small molecules against AlphaFold proteins!
