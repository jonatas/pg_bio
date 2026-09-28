# Tutorial 5: Cheminformatics & Small Molecules

While proteins are the biological machines of a cell, **small molecules** (like pharmaceutical drugs, metabolites, and ligands) are what interact with them. 

The universal language used to represent small chemical structures in a database is called **SMILES** (Simplified Molecular-Input Line-Entry System). 

In this tutorial, we will use our new `pg_bio` native `Molecule` type to parse, store, and query these chemical structures directly inside PostgreSQL, laying the foundation for in-database docking algorithms.

## 1. The `Molecule` Data Type
Instead of just storing chemistry strings in standard `TEXT` columns and hoping they are syntactically valid, `pg_bio` provides the native `Molecule` Postgres type.

When you attempt to cast a string into a Molecule using our `parse_smiles()` Rust function, it uses a blazing-fast pure Rust parser (`purr`) to validate the chemical graph instantly!

```sql
-- Valid Aspirin SMILES
SELECT parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O');
```
**Output:**
```text
                  parse_smiles                  
------------------------------------------------
 {"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","is_valid":true}
```

If we pass it an invalid chemical graph (like an unclosed aromatic ring `C(C`), the native Rust parser correctly flags it:
```sql
-- Invalid Chemical Graph
SELECT parse_smiles('C(C');
```
**Output:**
```text
       parse_smiles       
--------------------------
 {"smiles":"C(C","is_valid":false}
```

## 2. Chemical Substructure Searching
Finding specific chemical scaffolds (like a Benzene ring `C1=CC=CC=C1`) across millions of candidate drugs is the core of virtual screening.

We can use the `smiles_contains()` function to do structural inclusion queries natively against our `Molecule` type.

```sql
SELECT smiles_contains(
    parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O'), -- Aspirin
    'C1=CC=CC=C1' -- Benzene substructure
) as contains_benzene;
```
**Output:**
```text
 contains_benzene
------------------
 t
```

## 3. Morgan Fingerprints & Tanimoto Similarity
While substructure searching is great for finding exact scaffolds, modern drug discovery relies heavily on finding *similar* compounds. To do this, we transform the molecular graph into a **Morgan Fingerprint** (an Extended-Connectivity Fingerprint, or ECFP), which is a numerical representation of the molecule's substructures.

`pg_bio` can natively compute a 1024-dimensional Morgan Fingerprint from a SMILES string, returning it as a vector of floats, which allows us to natively score structural similarity using our custom **Tanimoto Operator (`%`)**.

```sql
-- Generate fingerprints and compute their Tanimoto similarity
SELECT 
    ROUND((smiles_to_fingerprint('C1=CC=CC=C1') % smiles_to_fingerprint('CC1=CC=CC=C1'))::numeric, 4) AS tanimoto_similarity;
```
**Output:**
```text
 tanimoto_similarity 
---------------------
              0.1143
```
*Note: Benzene (`C1=CC=CC=C1`) and Toluene (`CC1=CC=CC=C1`) share some ring substructure invariants, yielding a measurable structural overlap natively calculated in PostgreSQL.*

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
