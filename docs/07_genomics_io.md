# Tutorial 7: Direct Genomics I/O (VCF Parsing)

A massive bottleneck in modern biology is data loading. Genomics pipelines often require you to run external tools (like `bcftools`) to parse Variant Call Format (VCF) files, massage the data into CSVs, and load it into PostgreSQL using `COPY`.

`pg_bio` eliminates this entire pipeline by allowing you to query genomic files **natively from disk** directly in SQL.

## 1. Native VCF Parsing (`parse_vcf`)

Using the highly optimized Rust `noodles` crate, `pg_bio` exposes a Set Returning Function (SRF) `parse_vcf` that acts as a Virtual Table.

```sql
SELECT chrom, pos, id, ref_allele, alt_allele, qual, filter, info 
FROM parse_vcf('/data/patient_001.vcf.gz')
LIMIT 5;
```
* **Streaming Engine**: Because this uses a natively streaming Rust iterator (`VcfIterator`), it requires zero memory buildup. You can query a 50GB compressed VCF file without running out of RAM.
* **Format Support**: Automatically detects and streams both plain `.vcf` and block-gzipped `.vcf.gz` (BGZF) files.

## 2. Advanced Multiomics: Spatial VCF Injection

The true power of native SQL parsing is chaining it with `pg_bio`'s 1D and 3D spatial indexing architectures. You can insert variants directly into `int4range` arrays to find functional overlaps.

```sql
-- Inject patient variants straight from disk into a GiST-indexed table
CREATE TABLE patient_variants (
    chrom TEXT,
    genomic_range int4range,
    id TEXT,
    info TEXT
);

INSERT INTO patient_variants (chrom, genomic_range, id, info)
SELECT 
    chrom, 
    int4range(pos, pos + length(alt_allele)), 
    id, 
    info
FROM parse_vcf('/data/patient_001.vcf.gz')
WHERE info LIKE '%Pathogenic%';
```

You can now instantly cross-reference these variants with known promoters, enhancers, or 3D protein structures using `pg_bio`'s `<->` adjacency operators.
