# Tutorial 11: Sequence Analysis & Central Dogma

`pg_bio` now natively supports basic sequence analysis and transcription/translation inside PostgreSQL. This eliminates the need to export large sequence tables to Python and Biopython for standard transformations.

## 1. GC Content Calculation
GC content is crucial for primer design, predicting melting temperatures, and assessing genomic regions (like CpG islands).

```sql
SELECT 
    locus_name,
    sequence,
    bio_gc_content(sequence) AS gc_percentage
FROM genomic_regions
WHERE bio_gc_content(sequence) > 0.6;
```

## 2. In-Database DNA Translation
You can instantly translate DNA sequences into Amino Acid sequences using standard codon tables.

```sql
-- Convert a raw DNA sequence into a protein natively in SQL
SELECT bio_translate_dna('ATGCGTTAG'); 
-- Output: "MR*"
```

This becomes incredibly powerful when combined with sequence variant analysis. You can mutate a DNA sequence and instantly predict the resulting protein missense mutation without any external tools:

```sql
WITH mutations AS (
    SELECT 'ATGCGTTAG' AS wild_type, 'ATGCATTAG' AS mutant
)
SELECT 
    wild_type,
    mutant,
    bio_translate_dna(wild_type) AS wt_protein,
    bio_translate_dna(mutant) AS mutant_protein
FROM mutations;
```
