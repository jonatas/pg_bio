# On-the-Fly 3D Structural Folding (ESMFold API)

The dark proteome is massive, and while platforms like AlphaFold have pre-computed structures for millions of proteins, they haven't caught everything. Many extremophile orphans or highly novel synthetic constructs do not exist in the AlphaFold database.

To solve this, `pg_bio` natively embeds Meta's **Evolutionary Scale Modeling (ESMFold)** API directly into PostgreSQL.

## `bio_fold_sequence`

The `bio_fold_sequence(sequence TEXT) RETURNS TEXT` function takes a raw amino acid sequence string, requests an on-the-fly folding simulation via the ESMAtlas API, and streams back the physical atomic coordinates in standard `.pdb` format.

### Example: Predicting an Orphan's Structure

If you have discovered an orphan protein in your `proteins` or `orphan_discoveries` table, you can generate its 3D structure on demand:

```sql
SELECT 
    uniprot_id,
    sequence,
    bio_fold_sequence(sequence) AS pdb_data
FROM orphan_discoveries
WHERE is_structurally_characterized = false
LIMIT 1;
```

This allows you to construct completely autonomous discovery pipelines inside PostgreSQL. You can combine this with `bio_search_uniprot()` to fetch a raw sequence dynamically, fold it, and insert the resulting PDB directly into your structural tables—all without ever leaving SQL or writing a single line of Python.

### Performance Considerations

The ESMFold API computes the structure using a massive Language Model (LLM) designed for biology. 
- For short sequences (e.g., < 100 amino acids), the API typically responds in under a second.
- For extremely large proteins (e.g., > 1000 amino acids), the folding simulation can take several minutes.
- Because `bio_fold_sequence` relies on an external HTTP call via `ureq`, it blocks the current SQL transaction. If you are predicting structures for thousands of sequences, consider using `CROSS JOIN LATERAL` in batches, or wrapping it in a background worker queue to avoid hanging your primary database sessions.
