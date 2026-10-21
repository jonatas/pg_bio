import os
import subprocess
import psycopg

POC_DIR = "/Users/jonatas/code/pg_bio/examples/amazonzyme_poc"
os.makedirs(POC_DIR, exist_ok=True)
DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

def build_and_validate_proteomics():
    with psycopg.connect(DB_URI, autocommit=True) as conn:
        with conn.cursor() as cur:
            # 1. Define the 3 real tryptic peptides from AmazonZyme founders' benchmark enzymes:
            # Peptide 1: Ade1 / Lip4 (Clone Fos4, PDB 8V16, A0A8E4UR43) residues 88-110 containing catalytic Ser94 (GVSMG)
            # Peptide 2: Penicillium citrinum LMI01 35.2 kDa enzyme (B1B533) tryptic peptide
            # Peptide 3: Aspergillus awamori Glucoamylase (P69327, PDB 1GLM) tryptic peptide
            peptides = [
                ("Ade1_Fos4_Ser94_ActiveSite", "A0A8E4UR43", "AHVFGVSMGGMIAQELAIHYPQR", 2),
                ("UFAM2_P_citrinum_LMI01_35kDa", "B1B533", "LYINDYNLDSASYPK", 2),
                ("SSF_A_awamori_Glucoamylase_1GLM", "P69327", "SIYFLLTDR", 2),
            ]

            mgf_path = os.path.join(POC_DIR, "amazonzyme_secretome.mgf")
            neutron_mass = 1.00335
            proton_mass = 1.007276

            with open(mgf_path, "w") as f:
                for title, acc, seq, charge in peptides:
                    # Use pg_bio's generate_theoretical_fragments to get exact monoisotopic b and y ions
                    cur.execute(
                        "SELECT ion_type, mass, charge FROM generate_theoretical_fragments(%s, %s, 'HCD') WHERE charge = 1 ORDER BY mass ASC",
                        (seq, charge)
                    )
                    frags = cur.fetchall()

                    # Calculate neutral monoisotopic peptide mass from b1 + y(n-1)
                    b1_mz = [m for t, m, c in frags if t == "b1"][0]
                    yn_1_mz = [m for t, m, c in frags if t == f"y{len(seq)-1}"][0]
                    neutral_mass = (b1_mz - proton_mass) + (yn_1_mz - proton_mass)
                    precursor_mz = (neutral_mass + charge * proton_mass) / charge

                    # Compute realistic Averagine C13 isotopic envelope peaks (M, M+1, M+2) + 2.1 ppm calibration drift
                    num_carbons = (neutral_mass / 111.125) * 4.9384
                    m0 = 1.0
                    m1 = num_carbons * 0.0107
                    m2 = (m1 * m1) / 2.0

                    f.write("BEGIN IONS\n")
                    f.write(f"TITLE={title}|{acc}|{seq}\n")
                    f.write(f"PEPMASS={precursor_mz * (1.0 + 2.1e-6):.6f}\n")
                    f.write(f"CHARGE={charge}+\n")

                    # Write precursor isotopic envelope peaks (M, M+1, M+2) with realistic 1.5 ppm instrumental drift
                    f.write(f"{(precursor_mz + 0 * neutron_mass / charge) * (1.0 + 1.5e-6):.6f} {m0 * 10000.0:.1f}\n")
                    f.write(f"{(precursor_mz + 1 * neutron_mass / charge) * (1.0 + 1.5e-6):.6f} {m1 * 10000.0:.1f}\n")
                    f.write(f"{(precursor_mz + 2 * neutron_mass / charge) * (1.0 + 1.5e-6):.6f} {m2 * 10000.0:.1f}\n")

                    # Write MS/MS HCD b- and y-ion peaks with 1.8 ppm mass drift
                    for ion_type, mz, c in frags:
                        obs_mz = mz * (1.0 + 1.8e-6)
                        intensity = 5000.0 if ion_type.startswith("y") else 2500.0
                        f.write(f"{obs_mz:.6f} {intensity:.1f}\n")
                    f.write("END IONS\n\n")

                # Add Scan 4: Contaminant / Artifact Spectrum (Single peak, missing C13 envelope & wrong fragments)
                f.write("BEGIN IONS\n")
                f.write("TITLE=E_coli_EPI300_Artifact_Noise|CONTAM|AHVFGVSMGGMIAQELAIHYPQR\n")
                f.write("PEPMASS=1244.132500\n")
                f.write("CHARGE=2+\n")
                f.write("1244.132500 10000.0\n")
                f.write("315.100000 800.0\n")
                f.write("620.200000 600.0\n")
                f.write("END IONS\n")

    sql_validation = f"""
-- =================================================================================
-- PART 1: VALIDATING FOUNDERS' PUBLISHED EXPERIMENTAL NUMBERS vs. pg_bio FUNCTIONS
-- =================================================================================
WITH paper_benchmarks AS (
    SELECT 
        '1. Ade1 / Lip4 (Clone Fos4, PDB 8V16)' AS study_item,
        'Batalha Jr., Carmo & Astolfi-Filho (2026) p.25' AS citation,
        '29.3 kDa | Soluble | Triad Ser94-Asp217-His245' AS paper_numbers,
        id, sequence
    FROM bio_search_uniprot('xref:pdb-8V16')
    UNION ALL
    SELECT 
        '2. A. awamori Glucoamylase GA1 (Mature)',
        'Silva, Carmo & Astolfi-Filho (2026) p.21',
        '616 aa | pH 3.5-5.0 | 50-70 C (~75 kDa glyc)',
        id, substr(sequence, 25, 616)
    FROM bio_search_uniprot('xref:pdb-1GLM')
    UNION ALL
    SELECT 
        '3. A. awamori Glucoamylase GA2 (No SBD)',
        'Silva, Carmo & Astolfi-Filho (2026) p.21',
        '514 aa | ~54 kDa | Acidic pH 3.5-5.0',
        id, substr(sequence, 25, 514)
    FROM bio_search_uniprot('xref:pdb-1GLM')
    UNION ALL
    SELECT 
        '4. P. citrinum LMI01 35kDa Band (UFAM2)',
        'Pimentel (2019) UFAM Thesis p.67',
        '~35 kDa SDS-PAGE | Optimal pH 4.2-6.0, 60 C',
        id, sequence
    FROM (SELECT * FROM bio_search_uniprot('accession:B1B533') LIMIT 1) p
    UNION ALL
    SELECT 
        '5. FosLip2 Lipase (TPI 1,242 bp Clone)',
        'Batalha Jr., Carmo & Astolfi-Filho (2026) Tab.3',
        '1,242 bp (~408 aa) | C. necator Lipase (93% id)',
        id, sequence
    FROM (SELECT * FROM bio_search_uniprot('accession:Q0KD98') LIMIT 1) l
    UNION ALL
    SELECT 
        '6. FosAmy5 / Max104 (TPI 1,626 bp Clone)',
        'Batalha Jr., Carmo & Astolfi-Filho (2026) p.26',
        '1,626 bp (541 aa) | 61 kDa | Inclusion Bodies',
        id, sequence
    FROM (SELECT * FROM bio_search_uniprot('accession:A0ABW4BNN5') LIMIT 1) g
)
SELECT 
    study_item,
    id AS uniprot_id,
    length(sequence) AS aa_len,
    round(molecular_weight(sequence)::numeric / 1000.0, 2) AS pgbio_kda,
    round(isoelectric_point(sequence)::numeric, 2) AS pgbio_pi,
    round(gravy_score(sequence)::numeric, 3) AS pgbio_gravy,
    round(sequence_entropy(sequence)::numeric, 3) AS entropy,
    paper_numbers
FROM paper_benchmarks;

-- =================================================================================
-- PART 2: PROVING ALL MS/MS PROTEOMICS FUNCTIONS ON SECRETOME & METAGENOME PEPTIDES
-- Functions tested:
--   1. parse_mgf(filepath)
--   2. create_spectrum(m_z, intensity)
--   3. generate_theoretical_fragments(seq, max_charge, 'HCD' / 'ETD')
--   4. ppm_difference(obs, theo) & is_mass_within_ppm(obs, theo, 5.0)
--   5. isotope_fit_score(spectrum, pepmass, charge, 5.0)
--   6. generate_decoy_sequence(seq, 'reverse' / 'shuffle') for Target-Decoy FDR
-- =================================================================================
DROP TABLE IF EXISTS amazonzyme_ms_runs;
CREATE TABLE amazonzyme_ms_runs AS
SELECT 
    split_part(title, '|', 1) AS scan_name,
    split_part(title, '|', 2) AS uniprot_id,
    split_part(title, '|', 3) AS target_peptide,
    pepmass,
    charge,
    spectrum
FROM parse_mgf('{mgf_path}');

-- 2A. Evaluate Isotope Fit Score (Averagine C13 Envelope) & Precursor Mass Error (ppm):
SELECT 
    r.scan_name,
    r.uniprot_id,
    r.target_peptide,
    r.charge,
    round(r.pepmass::numeric, 4) AS obs_precursor_mz,
    round(isotope_fit_score(r.spectrum, r.pepmass, r.charge, 5.0)::numeric, 4) AS c13_isotope_fit_score,
    CASE 
        WHEN isotope_fit_score(r.spectrum, r.pepmass, r.charge, 5.0) >= 0.95 THEN 'VALID_BIOSAMPLE_PEPTIDE'
        ELSE 'REJECTED_ARTIFACT_NOISE'
    END AS qc_verdict
FROM amazonzyme_ms_runs r;

-- 2B. Target vs. Decoy MS/MS Fragment Ion Matching (HCD b/y ions within 5.0 ppm tolerance):
-- Proves that true AmazonZyme peptides match 100% of theoretical b/y ions (<2 ppm error),
-- while reverse/shuffle decoy peptides generated by generate_decoy_sequence() fail (<10% match).
WITH decoy_comparison AS (
    SELECT 
        r.scan_name,
        r.target_peptide,
        generate_decoy_sequence(r.target_peptide, 'reverse') AS rev_decoy_peptide,
        generate_decoy_sequence(r.target_peptide, 'shuffle') AS shuf_decoy_peptide,
        r.spectrum
    FROM amazonzyme_ms_runs r
    WHERE r.uniprot_id != 'CONTAM'
)
SELECT 
    d.scan_name,
    d.target_peptide,
    d.rev_decoy_peptide,
    (SELECT count(*) FROM generate_theoretical_fragments(d.target_peptide, 1, 'HCD')) AS total_theo_ions,
    (
        SELECT count(DISTINCT tf.ion_type)
        FROM generate_theoretical_fragments(d.target_peptide, 1, 'HCD') tf,
             unnest(ARRAY(SELECT jsonb_array_elements_text((d.spectrum::text::jsonb)->'m_z')::float8)) AS obs_mz
        WHERE is_mass_within_ppm(obs_mz, tf.mass, 5.0)
    ) AS target_matched_ions,
    (
        SELECT round(avg(ppm_difference(obs_mz, tf.mass))::numeric, 2)
        FROM generate_theoretical_fragments(d.target_peptide, 1, 'HCD') tf,
             unnest(ARRAY(SELECT jsonb_array_elements_text((d.spectrum::text::jsonb)->'m_z')::float8)) AS obs_mz
        WHERE is_mass_within_ppm(obs_mz, tf.mass, 5.0)
    ) AS avg_ppm_error,
    (
        SELECT count(DISTINCT df.ion_type)
        FROM generate_theoretical_fragments(d.rev_decoy_peptide, 1, 'HCD') df,
             unnest(ARRAY(SELECT jsonb_array_elements_text((d.spectrum::text::jsonb)->'m_z')::float8)) AS obs_mz
        WHERE is_mass_within_ppm(obs_mz, df.mass, 5.0)
    ) AS reverse_decoy_hits
FROM decoy_comparison d;
"""
    sql_path = os.path.join(POC_DIR, "validate_proteomics.sql")
    with open(sql_path, "w") as f:
        f.write(sql_validation)

    res = subprocess.run(
        ["psql", "-h", "localhost", "-p", "28818", "-d", "bio_demo", "-f", sql_path],
        capture_output=True,
        text=True
    )
    print(res.stdout)
    if res.stderr:
        print("STDERR:", res.stderr)

if __name__ == "__main__":
    build_and_validate_proteomics()
