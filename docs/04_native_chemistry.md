# Tutorial 4: Native Chemistry & Sequence Math

In the previous tutorials, we utilized `pg_bio` for high-performance vector search and Z-Order spatial indexing. 

In this tutorial, we will explore the **Native Chemistry Operations**. By pushing computational chemistry directly into the database using our Rust extensions, we eliminate the need to export data into Pandas or BioPython. This massively reduces I/O bottlenecks and memory overhead.

## 1. On-The-Fly Molecular Weight
Typically, querying for proteins within a specific mass range requires pre-computing the molecular weight and storing it in a separate column. 

With `pg_bio`, the `molecular_weight` function executes blazing-fast native Rust code to calculate the mass (in Daltons) of the amino acid sequence directly during the `SELECT` execution.

### How it works under the hood
The function accurately calculates mass using standard biochemical rules:
1. **Terminals (Water):** It starts with a base weight of **18.015 Daltons**. This accounts for the extra Hydrogen atom at the N-terminus and the Hydroxyl (OH) group at the C-terminus (which together form an intact water molecule).
2. **Per-Residue Mass:** It iterates through the sequence, adding the specific monoisotopic/average mass of each amino acid *residue* (e.g., `A` adds 71.079, `W` adds 186.213).
3. **Error Tolerance:** If the sequence contains unrecognized characters or gaps, they are safely ignored (adding `0.0`), preventing queries from crashing.

Let's find proteins weighing between 15 kDa and 50 kDa:

```sql
SELECT uniprot_id, molecular_weight(sequence) as mass_daltons
FROM proteins
WHERE molecular_weight(sequence) BETWEEN 15000 AND 50000
LIMIT 3;
```

**Output:**
```text
 uniprot_id | mass_daltons
------------+------------------
 A0A0U1RQG5 | 35144.83
 A0A0U1RR37 | 20280.51
 A0A1B0GTI1 | 20621.01
```

## 2. Biological Pattern Matching (PROSITE)
Proteins contain evolutionary conserved motifs (e.g., kinase phosphorylation sites, zinc fingers). The standard format for these is **PROSITE**. 

Instead of writing complex regular expressions in SQL, `pg_bio` provides `prosite_match`. Under the hood, this function compiles the PROSITE pattern into a highly-optimized Rust regex.

Let's search for proteins that contain a sequence matching `[ST]-x(2)-[RK]` (a Serine or Threonine, followed by exactly 2 of any amino acid, followed by an Arginine or Lysine).

```sql
SELECT uniprot_id, substring(sequence from 1 for 20) as seq_start
FROM proteins
WHERE prosite_match(sequence, '[ST]-x(2)-[RK]')
LIMIT 3;
```

**Output:**
```text
 uniprot_id |      seq_start       
------------+----------------------
 A0A0J9YXQ4 | MALAMLRDWCRWMGANAERS
 A0A0U1RQG5 | MSATGDQDLIQEDQEAPVNQ
 A0A0U1RR37 | MNQAFWKTYKSKVLQTLSGE
```

## 3. Spatial Aggregation: Center of Mass
Combining sequence matching with spatial analysis is where `pg_bio` truly shines. 

Suppose we found a target protein (`A6ZM04`) that matched our pattern, and we now want to find the exact geometric center of its atomic structure to set up a bounding box for Molecular Docking.

We can aggregate all coordinates of that protein using Postgres's native `array_agg`, and feed it into our custom `center_of_mass` function:

```sql
SELECT center_of_mass(array_agg(coord))
FROM protein_atoms
WHERE uniprot_id = 'A6ZM04';
```

**Output:**
```text
                      center_of_mass                          
-----------------------------------------------------------
 {"x":-12.783,"y":4.263,"z":-5.510,"name":"COM"}
```

## Conclusion
By embedding chemistry algorithms directly into PostgreSQL via Rust `pgrx`, your database is no longer just a storage engine—it is a live, high-performance computational biology platform!
