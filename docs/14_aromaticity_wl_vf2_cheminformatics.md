# Tutorial 14: Hückel $4n+2$ Aromaticity Perception, Weisfeiler-Lehman Hashing, VF2 Subgraph Isomorphism & Salt/Tautomer Normalization

## 1. The Scientist's Problem: The Kekulé vs. Aromatic Trap

In real-world chemical databases (ChEMBL, PubChem, ZINC) and enzyme-substrate discovery pipelines, identical molecules are routinely encoded in incompatible surface notations:
- **Kekulé vs. Aromatic Rings:** Benzene written as `C1=CC=CC=C1` (alternating single/double bonds) vs. `c1ccccc1` (delocalized aromatic $\pi$-system).
- **Salts & Counter-Ions:** PET monomer MHET stored as sodium salt `O=C(OCCO)C1=CC=C(C(=O)[O-])C=C1.[Na+]` vs. neutral acid `O=C(OCCO)c1ccc(C(=O)O)cc1`.
- **Prototropic Tautomers:** Aliphatic enols/imidic acids (`C=CO`, `CC(=N)O`) vs. keto/amide forms (`CC=O`, `CC(=O)N`).

Without a graph normalization pass between the SMILES lexer and the Weisfeiler-Lehman (WL) / Morgan ECFP4 hasher, `C1=CC=CC=C1` and `c1ccccc1` produce two completely different hashes—resulting in **0% recall** on canonical lookups and substructure filters.

---

## 2. The 4-Stage `pg_bio` Molecular Graph Engine

### Stage 1: Salt Stripping & Prototropic Tautomer Normalization
1. **Connected-Component Salt Stripping:** Splits disconnected `.` fragments via BFS, retains the largest heavy-atom component, and neutralizes counter-ion salt charges (`[O-].[Na+]` $\rightarrow$ `OH`, `[NH3+].[Cl-]` $\rightarrow$ `NH2`).
2. **1,3-Prototropic Tautomer Canonicalization:** Normalizes non-ring enols (`C=C-OH` $\rightarrow$ `CH-C=O`) and imidic acids (`N=C-OH` $\rightarrow$ `NH-C=O`) while preserving true aromatic phenols (`Oc1ccccc1`).

### Stage 2: Hückel’s $4n+2$ Aromaticity Perception
For every simple and fused cycle ($3..=10$), `pgbio_atom_pi_electrons` sums each ring atom's $\pi$-electron contribution using pure integer arithmetic:

| Ring Atom State | Chemical Examples | $\pi$-Electrons |
| :--- | :--- | :---: |
| **Has in-ring `=` or `:` bond** (or $sp^2$ `C`/`N`) | `C`, `N`, `c`, `n` in Benzene (`C1=CC=CC=C1`), Pyridine (`c1ccncc1`) | **1** |
| **Single-bonded heteroatom with lone pair** or `[C-]` | `[nH]` in Pyrrole (`C1=CC=CN1`), `O` in Furan (`c1ccoc1`), `S` in Thiophene, `[cH-]` | **2** |
| **Ring Carbon with exocyclic `=O` / `=N` / `=S`** or `[C+]` / `B` | Carbonyl `C(=O)` in 2-Pyridone (`O=C1C=CC=CN1`), Caffeine, Uracil | **0** |
| **Saturated $sp^3$ atom** (`sigma_degree >= 4`) or triple bond `#` | Cyclohexane (`C1CCCCC1`), Cyclopentadiene (`C1=CC=CC1`) | **Abort (`-1`)** |

If all atoms in the cycle are conjugated and `(pi_sum >= 6) && ((pi_sum - 2) & 3) == 0` ($\pi \in \{6, 10, 14, 18\}$), all cycle atoms are promoted to `is_aromatic = true` and all cycle bonds to `BOND_AROMATIC`. Anti-aromatic $4n$ rings (`C1=CC=C1` with $4\pi$, `C1=CC=CC=CC=C1` with $8\pi$) remain non-aromatic.

### Stage 3: Weisfeiler-Lehman (WL) 64-Bit Canonical Hash & 1024-Bit ECFP4 Bitset
- **Radius 0:** Computes vertex-permutation-invariant atom hashes from `(atomic_num, is_aromatic, h_count, formal_charge, degree, in_ring)` using MurmurHash3's `fmix64`.
- **Radius 1..2 (ECFP4):** Sorts neighbor `(bond_type, neighbor_hash)` pairs per atom and folds hashes into both a 64-bit canonical `wl_hash` (`molecule_wl_hash`) and a `16 × u64` (1024-bit) `ecfp_bits` array for hardware `POPCNT` Tanimoto similarity (`Molecule % Molecule`).

### Stage 4: VF2 Subgraph Isomorphism Recheck (`@>`)
The `@>` operator (`molecule_contains_molecule`) and `smiles_contains` execute a backtracking VF2 state machine matching query nodes `Q[i]` to target nodes `T[j]` by `(atomic_num, is_aromatic, degree, in_ring)` and verifying exact bond-type compatibility across all incident edges.

---

## 3. SQL Examples

```sql
-- 1. Canonical equivalence across Kekulé, Aromatic, Salted, and Tautomeric SMILES
SELECT
    molecule_wl_hash('C1=CC=CC=C1') AS kekule_hash,
    molecule_wl_hash('c1ccccc1')    AS aromatic_hash,
    molecules_are_identical('C1=CC=CC=C1', 'c1ccccc1') AS benzene_eq,
    molecules_are_identical('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C') AS caffeine_eq,
    molecules_are_identical('[Na+].O=C([O-])C1=CC=CC=C1', 'O=C(O)c1ccccc1') AS benzoate_salt_eq,
    molecules_are_identical('CC(=N)O', 'CC(=O)N') AS tautomer_eq;

-- 2. VF2 Substructure Search (@>) and Hardware-Popcount Tanimoto Similarity (%)
SELECT
    parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') @> parse_smiles('c1ccccc1') AS aspirin_has_benzene,
    parse_smiles('O=C(OCCO)C1=CC=C(C(=O)[O-])C=C1.[Na+]') % parse_smiles('O=C(OCCO)c1ccc(C(=O)O)cc1') AS mhet_tanimoto;
```

---

## 4. Benchmark Summary (`scripts/benchmark_aromaticity_wl_vf2.py`)

| Benchmark Metric (10,000 Molecules in PostgreSQL 18) | Unnormalized Baseline | `pg_bio` Native Engine | Speedup / Gain |
| :--- | :--- | :--- | :--- |
| **14-System Equivalence Recall** (Kekulé/Salt/Tautomer) | `0 / 14 (0.0%)` | **`14 / 14 (100.0%)`** | **+100.0% Recall** |
| **Swarm Enzyme-Substrate Family Recall** | `0 / 17 (0.0%)` | **`17 / 17 (100.0%)`** | **+100.0% Recall** |
| **B-Tree `wl_hash` Lookup** (10,000 rows) | `1,018.64 ms` (Python pull) | **`0.673 ms`** | **`1,513.7x` faster** |
| **1024-Bit Popcount Tanimoto Scan** (`Molecule % Molecule`) | `1,018.64 ms` (Python pull) | **`434.46 ms`** | **`2.3x` faster** |
| **VF2 Subgraph Isomorphism Scan** (`Molecule @> Molecule`) | `1,018.64 ms` (Python pull) | **`437.85 ms`** | **`2.3x` faster** |
