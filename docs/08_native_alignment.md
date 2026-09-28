# Tutorial 8: Ultra-Fast Native Sequence Alignment

One of the largest bottlenecks in molecular biology is pulling massive genomic or proteomic databases over the network into memory to run sequence alignment tools like BLAST or Smith-Waterman.

With `pg_bio`, we leverage the highly optimized `rust-bio` crate natively inside PostgreSQL. This runs SIMD-accelerated Smith-Waterman local alignment at C-like speeds without moving the data.

## 1. Native Smith-Waterman Search

You can score any two DNA/RNA or Protein sequences directly in SQL:

```sql
SELECT 
    target_name,
    smith_waterman_score(sequence, 'ATCGGCTA', 3, -1, -2) as score
FROM viral_genomes
WHERE smith_waterman_score(sequence, 'ATCGGCTA', 3, -1, -2) > 50
ORDER BY score DESC;
```
* **Performance:** Because this runs via `rust-bio` and is compiled with LLVM optimizations, it handles thousands of alignments per second natively against your rows.
* **Flexibility:** Custom match scores, mismatch penalties, and gap penalties are supported.

## 2. Using the Python SDK

We've exposed this inside the Python SDK so you can seamlessly query homologous genes:

```python
from pgbio import PgBioClient

client = PgBioClient("postgresql://localhost:28818/pg_bio")

# Align a query sequence against the database
matches = client.align_sequence("MKTLL", threshold=10, limit=5)
for m in matches:
    print(f"{m['name']}: Score {m['alignment_score']}")
```
