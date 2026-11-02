# Mass Spectrometry & Top-Down Proteomics in `pg_bio`

Top-Down proteomics focuses on the analysis of intact proteins (proteoforms), which is crucial for identifying Post-Translational Modifications (PTMs) and understanding true biological function. 

`pg_bio` brings extreme performance to Top-Down proteomics by allowing you to generate theoretical fragment masses, predict isotopic envelopes, score spectral matches, and estimate False Discovery Rates (FDR)—all natively within PostgreSQL!

## 1. The `Spectrum` Type & Static Datasets

To support experimental data, `pg_bio` introduces the highly-optimized `Spectrum` data type for storing mass-to-charge (`m/z`) and intensity vectors. 

To benchmark against static experimental datasets, you can directly parse standard MGF (Mascot Generic Format) files:

```sql
SELECT title, pepmass, charge, spectrum 
FROM parse_mgf('/data/static_test_set.mgf')
LIMIT 5;
```

## 2. Theoretical Fragment Generation

Different mass spectrometry collision methods break protein backbones differently. 
- **HCD/CID:** Breaks standard peptide bonds, resulting in **b-ions** (N-terminal) and **y-ions** (C-terminal).
- **ETD:** Uses electron transfer to break proteins much more gently. Crucially, ETD leaves fragile PTMs intact while generating **c-ions** and **z-ions**.

`pg_bio` allows you to generate these theoretical fragments on the fly using exact monoisotopic masses:

```sql
-- Generate HCD fragments
SELECT * FROM generate_theoretical_fragments('MFEGFERRLVD', 2, 'HCD');

-- Generate ETD fragments (better for PTM preservation)
SELECT * FROM generate_theoretical_fragments('MFEGFERRLVD', 2, 'ETD');
```

## 3. High-Precision Spectral Matching (Isotope Fit Score)

Comparing mass spectra isn't just about matching a single mass. Large fragments have complex *isotopic envelopes* (e.g., M, M+1, M+2 peaks primarily driven by natural Carbon-13 abundance). 

Instead of relying on basic mass deconvolution, `pg_bio` implements an **Averagine Model** Isotope Fit Score directly in SQL. It estimates the carbon count from the neutral mass, predicts the isotopic envelope, searches the observed `Spectrum` within a strict parts-per-million (ppm) tolerance, and scores the match using Cosine Similarity (0.0 to 1.0).

```sql
-- Compare theoretical heavy chain masses against an experimental spectrum
SELECT 
    m.title,
    isotope_fit_score(
        m.spectrum, 
        theoretical_mz := 1000.5, 
        charge := 2, 
        ppm_tolerance := 3.2
    ) as score
FROM parse_mgf('/data/experiment.mgf') m
WHERE isotope_fit_score(m.spectrum, 1000.5, 2, 3.2) > 0.85;
```

## 4. In-Database False Discovery Rate (FDR)

When scoring millions of permutations, random noise will inevitably create false positive matches. The "Target-Decoy" strategy is the gold standard for calculating the False Discovery Rate (FDR).

`pg_bio` provides a deterministic pseudo-random generator to create decoys on-the-fly:

```sql
-- Standard reversal
SELECT generate_decoy_sequence('MFEGFERRLVD', 'reverse');
-- Returns: DVLVRREGEFM

-- Deterministic shuffle (safe for indexing)
SELECT generate_decoy_sequence('MFEGFERRLVD', 'shuffle');

-- Completely random fake sequence of the same length
SELECT generate_decoy_sequence('MFEGFERRLVD', 'random');
```

You can use these decoys in a `UNION ALL` with your real targets. By measuring the distribution of high-scoring decoys vs high-scoring targets, you can automatically filter your Top-Down proteomics results to a strict 1% or 5% FDR natively!
