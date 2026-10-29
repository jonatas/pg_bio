import time
import psycopg

DB_URI = "postgresql://jonatas@localhost:28818/bio_demo"

# 14 Real-World Chemical Equivalence Pairs (Kekulé vs Aromatic, Salted vs Neutral, Enol/Imidic vs Keto/Amide)
VALIDATION_PAIRS = [
    ("Benzene (Kekulé vs Aromatic)", "C1=CC=CC=C1", "c1ccccc1", "c1ccccc1"),
    ("Pyridine (6π Heterocycle)", "C1=CC=NC=C1", "c1ccncc1", "c1ccncc1"),
    ("Pyrrole (Lone-Pair [nH] 6π)", "C1=CC=CN1", "c1cc[nH]c1", "c1cc[nH]c1"),
    ("Furan (Lone-Pair O 6π)", "O1C=CC=C1", "c1ccoc1", "c1ccoc1"),
    ("Thiophene (Lone-Pair S 6π)", "S1C=CC=C1", "c1ccsc1", "c1ccsc1"),
    ("Imidazole (His Catalytic Ring)", "C1=CN=CN1", "c1c[nH]cn1", "c1c[nH]cn1"),
    ("2-Pyridone (Exocyclic C=O 0π)", "O=C1C=CC=CN1", "O=c1cccc[nH]1", "O=c1cccc[nH]1"),
    ("Cyclopentadienyl Anion ([C-] 2π)", "[CH-]1C=CC=C1", "[cH-]1cccc1", "[cH-]1cccc1"),
    ("Naphthalene (10π Fused Kekulé A/B)", "C1=CC2=CC=CC=C2C=C1", "c1ccc2ccccc2c1", "c1ccccc1"),
    ("Caffeine (Fused Purine-2,6-dione)", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "c1ncn(C)c1"),
    ("Aspirin (Kekulé vs Reversed Arom)", "CC(=O)OC1=CC=CC=C1C(=O)O", "O=C(O)c1ccccc1OC(=O)C", "c1ccccc1C(=O)O"),
    ("Sodium Benzoate Salt -> Benzoic Acid", "[Na+].O=C([O-])C1=CC=CC=C1", "O=C(O)c1ccccc1", "c1ccccc1"),
    ("Glycine HCl Salt -> Neutral Glycine", "[NH3+]CC(=O)O.[Cl-]", "NCC(=O)O", "NCC(=O)O"),
    ("Acetimidic Acid -> Acetamide Tautomer", "CC(=N)O", "CC(=O)N", "CC(=O)N"),
]

# Cognate enzyme-substrate & cofactor library mapped directly to Swarm Discovery enzyme families
SWARM_SUBSTRATES = [
    (
        "Cutinase",
        "MHET (Mono-(2-hydroxyethyl)terephthalate)",
        "O=C(OCCO)C1=CC=C(C(=O)[O-])C=C1.[Na+]",
        "O=C(OCCO)c1ccc(C(=O)O)cc1",
        "c1ccc(C(=O)O)cc1",
        "PET Polyester Depolymerization",
    ),
    (
        "Hydrolase",
        "BHET (Bis(2-hydroxyethyl)terephthalate)",
        "OCCOC(=O)C1=CC=C(C(=O)OCCO)C=C1",
        "OCCOC(=O)c1ccc(C(=O)OCCO)cc1",
        "c1ccccc1",
        "Ester Bond Hydrolysis",
    ),
    (
        "Carboxylesterase",
        "p-Nitrophenyl Butyrate",
        "CCCC(=O)OC1=CC=C(C=C1)N(=O)=O",
        "CCCC(=O)Oc1ccc(N(=O)=O)cc1",
        "c1ccccc1",
        "Carboxylester Cleavage",
    ),
    (
        "Lipase",
        "p-Nitrophenyl Acetate",
        "CC(=O)OC1=CC=C(C=C1)N(=O)=O",
        "CC(=O)Oc1ccc(N(=O)=O)cc1",
        "CC(=O)Oc1ccccc1",
        "Acyl-Ester Hydrolysis",
    ),
    (
        "Peroxidase",
        "Guaiacol (Lignin Monolignol Model)",
        "COC1=CC=CC=C1O",
        "COc1ccccc1O",
        "Oc1ccccc1",
        "Radical Phenolic Oxidation",
    ),
    (
        "Manganese peroxidase",
        "Syringol (2,6-Dimethoxyphenol)",
        "COC1=CC=CC(OC)=C1O",
        "COc1cccc(OC)c1O",
        "Oc1ccccc1",
        "Lignin Depolymerization",
    ),
    (
        "Catalase",
        "Salicylic Acid Phenolic Donor",
        "[Na+].O=C([O-])C1=CC=CC=C1O",
        "O=C(O)c1ccccc1O",
        "c1ccccc1O",
        "Peroxide Detoxification",
    ),
    (
        "Dehalogenase",
        "Pentachlorophenolate Sodium Salt",
        "[Na+].[O-]C1=C(Cl)C(Cl)=C(Cl)C(Cl)=C1Cl",
        "Oc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl",
        "Oc1ccccc1Cl",
        "Organohalide Reductive/Hydrolytic Dechlorination",
    ),
    (
        "Carbonic anhydrase",
        "Acetazolamide (1,3,4-Thiadiazole Inhibitor)",
        "CC(=O)NC1=NN=C(S(N)(=O)=O)S1",
        "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
        "c1nncs1",
        "Zn2+ Active-Site Sulfonamide Chelation",
    ),
    (
        "Kinase",
        "Adenine Core (ATP Nucleobase)",
        "NC1=NC=NC2=C1N=CN2",
        "Nc1ncnc2[nH]cnc12",
        "c1ncnc2[nH]cnc12",
        "Phosphoryl Transfer Cofactor Recognition",
    ),
    (
        "Nuclease",
        "2-Pyridone / Uracil-like Nucleobase Analog",
        "O=C1C=CC=CN1",
        "O=c1cccc[nH]1",
        "O=c1cccc[nH]1",
        "Phosphodiester Cleavage Recognition",
    ),
    (
        "Ligase",
        "Nicotinamide (NAD+ Cofactor Headgroup)",
        "NC(=O)C1=CC=CN=C1",
        "NC(=O)c1cccnc1",
        "c1ccncc1",
        "Adenylation / NAD+ Cofactor Binding",
    ),
    (
        "Luciferase",
        "Benzothiazole Luciferin Core",
        "OC1=CC=C2N=CSC2=C1",
        "Oc1ccc2ncsc2c1",
        "c1ccc2ncsc2c1",
        "Bioluminescent Adenylation & Oxidation",
    ),
    (
        "Amidase",
        "Acetanilide / Paracetamol Amide Substrate",
        "CC(=O)NC1=CC=C(O)C=C1",
        "CC(=O)Nc1ccc(O)cc1",
        "CC(=O)Nc1ccccc1",
        "Aromatic Amide C-N Hydrolysis",
    ),
    (
        "Protease",
        "N-Acetyl-L-tyrosine Amide Model",
        "CC(=O)NC(CC1=CC=C(O)C=C1)C(=O)N",
        "CC(=O)NC(Cc1ccc(O)cc1)C(=O)N",
        "c1ccc(O)cc1",
        "Endopeptidase Peptide Bond Cleavage",
    ),
    (
        "Nitrite reductase",
        "Phenazine Electron Mediator",
        "C1=CC=C2N=C3C=CC=CC3=NC2=C1",
        "c1ccc2nc3ccccc3nc2c1",
        "c1ccc2nc3ccccc3nc2c1",
        "Denitrification Electron Transfer",
    ),
    (
        "Methane monooxygenase",
        "Phenol Hydroxylation Product",
        "OC1=CC=CC=C1",
        "Oc1ccccc1",
        "c1ccccc1",
        "Selective C-H Oxidation",
    ),
]


def seed_swarm_substrates_and_validate(conn):
    print("=" * 92)
    print("1. VALIDATION: HUCKEL 4n+2 AROMATICITY, WL HASH, VF2 SUBGRAPH (@>), AND SALT/TAUTOMER")
    print("=" * 92)

    with conn.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS swarm_enzyme_substrates CASCADE;")
        cur.execute("""
            CREATE TABLE swarm_enzyme_substrates (
                id SERIAL PRIMARY KEY,
                family_name TEXT UNIQUE NOT NULL,
                substrate_name TEXT NOT NULL,
                ingested_kekule_or_salt_smiles TEXT NOT NULL,
                canonical_query_smiles TEXT NOT NULL,
                warhead_substructure_smiles TEXT NOT NULL,
                catalytic_role TEXT NOT NULL,
                mol Molecule NOT NULL,
                wl_hash BIGINT NOT NULL
            );
            CREATE INDEX idx_swarm_substrates_wl_hash ON swarm_enzyme_substrates (wl_hash);
            CREATE INDEX idx_swarm_substrates_family ON swarm_enzyme_substrates (family_name);
        """)

        for fam, sub_name, kek_salt, aro_query, warhead, role in SWARM_SUBSTRATES:
            cur.execute(
                """
                INSERT INTO swarm_enzyme_substrates
                    (family_name, substrate_name, ingested_kekule_or_salt_smiles, canonical_query_smiles,
                     warhead_substructure_smiles, catalytic_role, mol, wl_hash)
                VALUES (%s, %s, %s, %s, %s, %s, parse_smiles(%s), molecule_wl_hash(%s));
                """,
                (fam, sub_name, kek_salt, aro_query, warhead, role, kek_salt, kek_salt),
            )
        conn.commit()

        # Validate all 14 chemical equivalence pairs
        print(f"{'Chemical System':<38} | {'Naive Str Eq':<12} | {'WL Hash Match':<13} | {'Tanimoto':<8} | {'VF2 @>':<6}")
        print("-" * 92)
        naive_hits = 0
        pgbio_hits = 0
        for label, s_kek, s_aro, sub_q in VALIDATION_PAIRS:
            cur.execute(
                """
                SELECT
                    (%s = %s) AS naive_eq,
                    (molecule_wl_hash(%s) = molecule_wl_hash(%s)) AS wl_eq,
                    molecules_are_identical(%s, %s) AS ident,
                    (parse_smiles(%s) %% parse_smiles(%s)) AS tanimoto,
                    (parse_smiles(%s) @> parse_smiles(%s)) AS vf2_sub
                """,
                (s_kek, s_aro, s_kek, s_aro, s_kek, s_aro, s_kek, s_aro, s_kek, sub_q),
            )
            naive_eq, wl_eq, ident, tanimoto, vf2_sub = cur.fetchone()
            if naive_eq:
                naive_hits += 1
            if wl_eq and ident and abs(tanimoto - 1.0) < 1e-9 and vf2_sub:
                pgbio_hits += 1
            print(
                f"{label:<38} | {str(naive_eq):<12} | {str(wl_eq and ident):<13} | {tanimoto:<8.4f} | {str(vf2_sub):<6}"
            )

        print("-" * 92)
        print(
            f"Recall Without Normalization: {naive_hits}/{len(VALIDATION_PAIRS)} ({100.0*naive_hits/len(VALIDATION_PAIRS):.1f}%)  vs.  "
            f"pg_bio Normalized Recall: {pgbio_hits}/{len(VALIDATION_PAIRS)} ({100.0*pgbio_hits/len(VALIDATION_PAIRS):.1f}%)\n"
        )
        assert pgbio_hits == len(VALIDATION_PAIRS), "All validation pairs must match!"

        # Verify Swarm Enzyme Substrate table queries (ingested as Kekule/Salt, queried as Aromatic)
        cur.execute("""
            SELECT COUNT(*)
            FROM swarm_enzyme_substrates
            WHERE wl_hash = molecule_wl_hash(canonical_query_smiles)
              AND (mol % parse_smiles(canonical_query_smiles)) = 1.0
              AND (mol @> parse_smiles(warhead_substructure_smiles));
        """)
        matched_swarm = cur.fetchone()[0]
        print(
            f"[Swarm Substrate Table] Verified {matched_swarm}/{len(SWARM_SUBSTRATES)} enzyme families match across "
            f"Kekulé/Salted ingestion -> Aromatic WL hash + VF2 warhead `@>` screening!\n"
        )
        assert matched_swarm == len(SWARM_SUBSTRATES)


def run_10k_benchmark(conn, n_rows: int = 10000):
    print("=" * 92)
    print(f"2. PERFORMANCE BENCHMARK: {n_rows:,} MOLECULES IN POSTGRESQL 18 (COMPUTE-TO-DATA vs PYTHON)")
    print("=" * 92)

    templates = [
        "C1=CC=CC=C1",                           # Kekulé Benzene
        "CC(=O)OC1=CC=CC=C1C(=O)[O-].[Na+]",      # Sodium Aspirin Salt (Kekulé)
        "O=C(OCCO)C1=CC=C(C(=O)[O-])C=C1.[Na+]",  # Sodium MHET (PET monomer, Kekulé)
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",           # Caffeine (Kekulé)
        "CC(=O)NC1=NN=C(S(N)(=O)=O)S1",           # Acetazolamide (Kekulé)
        "NC1=NC=NC2=C1N=CN2",                     # Adenine (Kekulé)
        "C1=CC2=CC=CC=C2C=C1",                    # Naphthalene (Kekulé)
        "COC1=CC=CC=C1O",                         # Guaiacol (Kekulé)
        "CC(=O)NC1=CC=C(O)C=C1",                  # Paracetamol (Kekulé)
        "C1CCCCC1",                               # Cyclohexane (Aliphatic sp3 control)
    ]

    with conn.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS cheminfo_bench_10k CASCADE;")
        cur.execute("""
            CREATE TABLE cheminfo_bench_10k (
                id SERIAL PRIMARY KEY,
                raw_smiles TEXT NOT NULL,
                mol Molecule NOT NULL,
                wl_hash BIGINT NOT NULL
            );
        """)

        t0 = time.perf_counter()
        # Generate 10,000 rows inside PostgreSQL using set-based INSERT + normalization
        cur.execute(
            """
            WITH raw_templates AS (
                SELECT unnest(%s::text[]) AS tpl, generate_series(1, %s) AS idx
            )
            INSERT INTO cheminfo_bench_10k (raw_smiles, mol, wl_hash)
            SELECT
                rt.tpl,
                parse_smiles(rt.tpl),
                molecule_wl_hash(rt.tpl)
            FROM raw_templates rt, generate_series(1, %s) g;
            """,
            (templates, len(templates), n_rows // len(templates)),
        )
        cur.execute("CREATE INDEX idx_bench_10k_wl_hash ON cheminfo_bench_10k (wl_hash);")
        conn.commit()
        ingest_s = time.perf_counter() - t0
        print(f"Ingested & normalized {n_rows:,} molecules (Salt Strip + Hückel 4n+2 + WL + ECFP4): {ingest_s*1000:.1f} ms")

        # Benchmark 1: Indexed WL Hash Lookup (Querying Aromatic `O=C(O)c1ccccc1OC(=O)C` against Kekulé/Salted DB)
        arom_aspirin = "O=C(O)c1ccccc1OC(=O)C"
        t0 = time.perf_counter()
        for _ in range(20):
            cur.execute(
                """
                SELECT COUNT(*) FROM cheminfo_bench_10k
                WHERE wl_hash = molecule_wl_hash(%s);
                """,
                (arom_aspirin,),
            )
            cnt_wl = cur.fetchone()[0]
        wl_ms = ((time.perf_counter() - t0) / 20.0) * 1000.0

        # Benchmark 2: Full-table 1024-bit Hardware Popcount Tanimoto Similarity (`mol % parse_smiles(...)`)
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (SELECT parse_smiles(%s) AS qmol)
            SELECT COUNT(*) FROM cheminfo_bench_10k b, q
            WHERE (b.mol %% q.qmol) >= 0.85;
            """,
            (arom_aspirin,),
        )
        cnt_tan = cur.fetchone()[0]
        tan_ms = (time.perf_counter() - t0) * 1000.0

        # Benchmark 3: In-Database VF2 Subgraph Isomorphism (`mol @> parse_smiles('c1ccccc1C(=O)O')`)
        salicylate_query = "c1ccccc1C(=O)O"
        t0 = time.perf_counter()
        cur.execute(
            """
            WITH q AS (SELECT parse_smiles(%s) AS qmol)
            SELECT COUNT(*) FROM cheminfo_bench_10k b, q
            WHERE b.mol @> q.qmol;
            """,
            (salicylate_query,),
        )
        cnt_vf2 = cur.fetchone()[0]
        vf2_ms = (time.perf_counter() - t0) * 1000.0

        # Benchmark 4: Python Client-Side Pull + Python Graph/Fingerprint Reconstruction
        t0 = time.perf_counter()
        cur.execute("SELECT id, raw_smiles, (mol)::text FROM cheminfo_bench_10k;")
        rows = cur.fetchall()
        import json
        q_wl = None
        py_matches = 0
        for _, raw_s, mol_json in rows:
            m = json.loads(mol_json)
            # Re-evaluate bitset popcount in Python loop
            bits = m.get("ecfp_bits", [])
            if sum(bin(w).count("1") for w in bits) > 0 and m.get("is_valid"):
                py_matches += 1
        py_ms = (time.perf_counter() - t0) * 1000.0

        print("-" * 92)
        print(f"1. B-Tree Canonical WL Hash Lookup (Kekulé/Salt <- Aromatic Query): {wl_ms:7.3f} ms  (Hits: {cnt_wl:,}) -> {py_ms/wl_ms:6.1f}x faster than Python pull")
        print(f"2. Native Rust 1024-bit Popcount Tanimoto Scan (10,000 rows):       {tan_ms:7.3f} ms  (Hits: {cnt_tan:,}) -> {py_ms/tan_ms:6.1f}x faster than Python pull")
        print(f"3. Native Rust VF2 Subgraph Isomorphism `@>` Scan (10,000 rows):    {vf2_ms:7.3f} ms  (Hits: {cnt_vf2:,}) -> {py_ms/vf2_ms:6.1f}x faster than Python pull")
        print(f"4. Python Client-Side Fetch + JSON Graph Deserialization Baseline:  {py_ms:7.3f} ms")
        print("=" * 92)


if __name__ == "__main__":
    conn = psycopg.connect(DB_URI)
    seed_swarm_substrates_and_validate(conn)
    run_10k_benchmark(conn, 10000)
